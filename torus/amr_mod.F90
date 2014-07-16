module amr_mod

  
 
  ! 21 nov
  ! routines for adaptive mesh refinement. nhs
  ! twod stuff added by tjh started 25/08/04

  use amr_utils_mod
  use discwind_class
  USE octal_mod, only: OCTAL, wrapperArray, octalWrapper, subcellCentre, cellVolume, &
       allocateattribute, copyattribute, deallocateattribute, returndphi, subcellCorners
  use utils_mod, only: blackbody, logint, loginterp, locate, solvequaddble, spline, splint, regular_tri_quadint
  use romanova_class, only: romanova
  use gridtype_mod, only:   gridtype, hydrospline
#ifdef SPH
  USE cluster_class, only:  cluster
#endif
  use mpi_global_mod, only: myRankGlobal
  IMPLICIT NONE

  type STREAMTYPE
     integer :: nSamples
     type(VECTOR),pointer :: position(:) => null()
     type(VECTOR), pointer :: direction(:) => null()
     type(VECTOR), pointer :: velocity(:) => null()
     real, pointer :: speed(:) => null()
     real(double),pointer :: distanceAlongStream(:)   => null()
     real(double), pointer :: rho(:)  => null()
     real, pointer :: temperature(:)  => null()
     real(double),pointer :: streamRadius(:)  => null()
  end type STREAMTYPE

  type(STREAMTYPE),save :: globalStream(5000)
  integer,save :: globalnStream

  real(double), parameter :: amr_min_rho = 1.0e-30_db 

  integer :: mass_split, mass_split2, density_split, velocity_split, both_split, maxdensity_split
  integer :: scaleheighta_count, scaleheightb_count, scaleheightc_count
  TYPE(octal), POINTER :: recentOctal

CONTAINS
  SUBROUTINE fill_Velocity_Corners(this,thisOctal, debug)
    ! store the velocity values at the subcell corners of an octal so
    !   that they can be used for interpolation.

    IMPLICIT NONE
    type(discwind_type) :: this
    TYPE(octal), INTENT(INOUT) :: thisOctal
    logical, optional :: debug
    logical :: writedebug
    real(oct)      :: r1, r2, r3
    real(oct)      :: phi1, phi2, phi3
    real(oct)      :: x1, x2, x3
    real(oct)      :: y1, y2, y3
    real(oct)      :: z1, z2, z3
    

    writedebug = .false.
    if (present(debug)) writedebug=debug

    if (thisOctal%oneD) then
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       y1 = 0.d0
       z1 = 0.d0
       thisOctal%cornerVelocity(1) = discwind_velocity(this,vector(x1,y1,z1))
       thisOctal%cornerVelocity(2) = discwind_velocity(this,vector(x2,y1,z1))
       thisOctal%cornerVelocity(3) = discwind_velocity(this,vector(x3,y1,z1))
       goto 666
    endif


    if (thisOctal%threed) then
       if (.not.thisOctal%cylindrical) then ! 3d cartesian case
          ! we first store the values we use to assemble the position vectors
          
          x1 = thisOctal%centre%x - thisOctal%subcellSize
          x2 = thisOctal%centre%x
          x3 = thisOctal%centre%x + thisOctal%subcellSize
          
          y1 = thisOctal%centre%y - thisOctal%subcellSize
          y2 = thisOctal%centre%y
          y3 = thisOctal%centre%y + thisOctal%subcellSize
          
          z1 = thisOctal%centre%z - thisOctal%subcellSize
          z2 = thisOctal%centre%z
          z3 = thisOctal%centre%z + thisOctal%subcellSize
                    
           ! now store the 'base level' values
          
          thisOctal%cornerVelocity(1) = discwind_velocity(this,vector(x1,y1,z1))
          thisOctal%cornerVelocity(2) = discwind_velocity(this,vector(x2,y1,z1))
          thisOctal%cornerVelocity(3) = discwind_velocity(this,vector(x3,y1,z1))
          thisOctal%cornerVelocity(4) = discwind_velocity(this,vector(x1,y2,z1))
          thisOctal%cornerVelocity(5) = discwind_velocity(this,vector(x2,y2,z1))
          thisOctal%cornerVelocity(6) = discwind_velocity(this,vector(x3,y2,z1))
          thisOctal%cornerVelocity(7) = discwind_velocity(this,vector(x1,y3,z1))
          thisOctal%cornerVelocity(8) = discwind_velocity(this,vector(x2,y3,z1))
          thisOctal%cornerVelocity(9) = discwind_velocity(this,vector(x3,y3,z1))
          
          ! middle level
          
          thisOctal%cornerVelocity(10) = discwind_velocity(this,vector(x1,y1,z2))
          thisOctal%cornerVelocity(11) = discwind_velocity(this,vector(x2,y1,z2))
          thisOctal%cornerVelocity(12) = discwind_velocity(this,vector(x3,y1,z2))
          thisOctal%cornerVelocity(13) = discwind_velocity(this,vector(x1,y2,z2))
          thisOctal%cornerVelocity(14) = discwind_velocity(this,vector(x2,y2,z2))
          thisOctal%cornerVelocity(15) = discwind_velocity(this,vector(x3,y2,z2))
          thisOctal%cornerVelocity(16) = discwind_velocity(this,vector(x1,y3,z2))
          thisOctal%cornerVelocity(17) = discwind_velocity(this,vector(x2,y3,z2))
          thisOctal%cornerVelocity(18) = discwind_velocity(this,vector(x3,y3,z2))
          
          ! top level
          
          thisOctal%cornerVelocity(19) = discwind_velocity(this,vector(x1,y1,z3))
          thisOctal%cornerVelocity(20) = discwind_velocity(this,vector(x2,y1,z3))
          thisOctal%cornerVelocity(21) = discwind_velocity(this,vector(x3,y1,z3))
          thisOctal%cornerVelocity(22) = discwind_velocity(this,vector(x1,y2,z3))
          thisOctal%cornerVelocity(23) = discwind_velocity(this,vector(x2,y2,z3))
          thisOctal%cornerVelocity(24) = discwind_velocity(this,vector(x3,y2,z3))
          thisOctal%cornerVelocity(25) = discwind_velocity(this,vector(x1,y3,z3))
          thisOctal%cornerVelocity(26) = discwind_velocity(this,vector(x2,y3,z3))
          thisOctal%cornerVelocity(27) = discwind_velocity(this,vector(x3,y3,z3))

       else ! cylindrical 
          if (thisOctal%splitAzimuthally) then
             z1 = thisOctal%centre%z - thisOctal%subcellSize
             z2 = thisOctal%centre%z
             z3 = thisOctal%centre%z + thisOctal%subcellSize
             phi1 = thisOctal%phiMin
             phi2 = thisOctal%phi 
             phi3 = thisOctal%phiMax
             r1 = thisOctal%r - thisOctal%subcellSize
             r2 = thisOctal%r
             r3 = thisOctal%r + thisOctal%subcellSize

             ! bottom level

             thisOctal%cornerVelocity(1) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z1))
             thisOctal%cornerVelocity(2) = discwind_velocity(this,vector(r2*cos(phi1),r2*sin(phi1),z1))
             thisOctal%cornerVelocity(3) = discwind_velocity(this,vector(r3*cos(phi1),r3*sin(phi1),z1))
             thisOctal%cornerVelocity(4) = discwind_velocity(this,vector(r1*cos(phi2),r1*sin(phi2),z1))
             thisOctal%cornerVelocity(5) = discwind_velocity(this,vector(r2*cos(phi2),r2*sin(phi2),z1))
             thisOctal%cornerVelocity(6) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z1))
             thisOctal%cornerVelocity(7) = discwind_velocity(this,vector(r1*cos(phi3),r1*sin(phi3),z1))
             thisOctal%cornerVelocity(8) = discwind_velocity(this,vector(r2*cos(phi3),r2*sin(phi3),z1))
             thisOctal%cornerVelocity(9) = discwind_velocity(this,vector(r3*cos(phi3),r3*sin(phi3),z1))

             thisOctal%cornerVelocity(10) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z2))
             thisOctal%cornerVelocity(11) = discwind_velocity(this,vector(r2*cos(phi1),r2*sin(phi1),z2))
             thisOctal%cornerVelocity(12) = discwind_velocity(this,vector(r3*cos(phi1),r3*sin(phi1),z2))
             thisOctal%cornerVelocity(13) = discwind_velocity(this,vector(r1*cos(phi2),r1*sin(phi2),z2))
             thisOctal%cornerVelocity(14) = discwind_velocity(this,vector(r2*cos(phi2),r2*sin(phi2),z2))
             thisOctal%cornerVelocity(15) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z2))
             thisOctal%cornerVelocity(16) = discwind_velocity(this,vector(r1*cos(phi3),r1*sin(phi3),z2))
             thisOctal%cornerVelocity(17) = discwind_velocity(this,vector(r2*cos(phi3),r2*sin(phi3),z2))
             thisOctal%cornerVelocity(18) = discwind_velocity(this,vector(r3*cos(phi3),r3*sin(phi3),z2))

             thisOctal%cornerVelocity(19) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z3))
             thisOctal%cornerVelocity(20) = discwind_velocity(this,vector(r2*cos(phi1),r2*sin(phi1),z3))
             thisOctal%cornerVelocity(21) = discwind_velocity(this,vector(r3*cos(phi1),r3*sin(phi1),z3))
             thisOctal%cornerVelocity(22) = discwind_velocity(this,vector(r1*cos(phi2),r1*sin(phi2),z3))
             thisOctal%cornerVelocity(23) = discwind_velocity(this,vector(r2*cos(phi2),r2*sin(phi2),z3))
             thisOctal%cornerVelocity(24) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z3))
             thisOctal%cornerVelocity(25) = discwind_velocity(this,vector(r1*cos(phi3),r1*sin(phi3),z3))
             thisOctal%cornerVelocity(26) = discwind_velocity(this,vector(r2*cos(phi3),r2*sin(phi3),z3))
             thisOctal%cornerVelocity(27) = discwind_velocity(this,vector(r3*cos(phi3),r3*sin(phi3),z3))


          else

             z1 = thisOctal%centre%z - thisOctal%subcellSize
             z2 = thisOctal%centre%z
             z3 = thisOctal%centre%z + thisOctal%subcellSize
             phi1 = thisOctal%phi - thisOctal%dPhi/2.d0
             phi2 = thisOctal%phi + thisOctal%dPhi/2.d0
             r1 = thisOctal%r - thisOctal%subcellSize
             r2 = thisOctal%r
             r3 = thisOctal%r + thisOctal%subcellSize


             ! bottom level

             thisOctal%cornerVelocity(1) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z1))
             thisOctal%cornerVelocity(2) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z1))
             thisOctal%cornerVelocity(3) = discwind_velocity(this,vector(r2*cos(phi1),r2*sin(phi1),z1))
             thisOctal%cornerVelocity(4) = discwind_velocity(this,vector(r2*cos(phi2),r2*sin(phi2),z1))
             thisOctal%cornerVelocity(5) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z1))
             thisOctal%cornerVelocity(6) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z1))

             ! middle level

             thisOctal%cornerVelocity(7) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z2))
             thisOctal%cornerVelocity(8) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z2))
             thisOctal%cornerVelocity(9) = discwind_velocity(this,vector(r2*cos(phi1),r2*sin(phi1),z2))
             thisOctal%cornerVelocity(10) = discwind_velocity(this,vector(r2*cos(phi2),r2*sin(phi2),z2))
             thisOctal%cornerVelocity(11) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z2))
             thisOctal%cornerVelocity(12) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z2))

             ! top level

             thisOctal%cornerVelocity(13) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z3))
             thisOctal%cornerVelocity(14) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z3))
             thisOctal%cornerVelocity(15) = discwind_velocity(this,vector(r2*cos(phi1),r2*sin(phi1),z3))
             thisOctal%cornerVelocity(16) = discwind_velocity(this,vector(r2*cos(phi2),r2*sin(phi2),z3))
             thisOctal%cornerVelocity(17) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z3))
             thisOctal%cornerVelocity(18) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z3))

          endif
       endif
    else       
       
    ! we first store the values we use to assemble the position vectors
       
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       
       z1 = thisOctal%centre%z - thisOctal%subcellSize
       z2 = thisOctal%centre%z
       z3 = thisOctal%centre%z + thisOctal%subcellSize
       
       ! now store the 'base level' values
       
       thisOctal%cornerVelocity(1) = discwind_velocity(this,vector(x1,0.d0,z1))
       thisOctal%cornerVelocity(2) = discwind_velocity(this,vector(x2,0.d0,z1))
       thisOctal%cornerVelocity(3) = discwind_velocity(this,vector(x3,0.d0,z1))
       thisOctal%cornerVelocity(4) = discwind_velocity(this,vector(x1,0.d0,z2))
       thisOctal%cornerVelocity(5) = discwind_velocity(this,vector(x2,0.d0,z2))
       thisOctal%cornerVelocity(6) = discwind_velocity(this,vector(x3,0.d0,z2))
       thisOctal%cornerVelocity(7) = discwind_velocity(this,vector(x1,0.d0,z3))
       thisOctal%cornerVelocity(8) = discwind_velocity(this,vector(x2,0.d0,z3))
       thisOctal%cornerVelocity(9) = discwind_velocity(this,vector(x3,0.d0,z3))
    endif
666 continue

!    if(isnan(thisOctal%cornerVelocity(1)%x)) then
!          write(*,*) "cornervel",thisOctal%cornerVelocity(1)
!          write(*,*) discwind_velocity(this,vector(x1,0.d0,z1),grid)
!          write(*,*) x1,z1
!          write(*,*) x2,z2
!          write(*,*) x3,z3
!       enddo
!    endif
    
  END SUBROUTINE fill_Velocity_Corners

  RECURSIVE SUBROUTINE add_discwind(thisOctal,grid,this, limitscalar)
    ! uses an external function to decide whether to split a subcell of
    !   the current octal. 

    use inputs_mod, only : maxDepthAMR
    IMPLICIT NONE
    TYPE(OCTAL), POINTER :: thisOctal
    TYPE(gridtype), INTENT(INOUT) :: grid   ! need to pass the grid through to the 
    type(discwind_type), intent(in)  :: this     
    ! the value the decideSplit function uses to  decide whether or not to split cell.
    real(double), intent(in) :: limitscalar 
    !
    !
    TYPE(OCTAL), POINTER :: child
    INTEGER              :: i, j, k    ! loop counters
    integer :: iindex
    integer :: isubcell


    DO iSubcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(isubcell)) cycle
       IF ((need_to_split2(thisOctal,isubcell,this).and.(thisOctal%nDepth < maxDepthAMR))) then

          CALL add_new_children_discwind(thisOctal, isubcell, grid, this)

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
       child => thisOctal%child(iIndex)
       CALL add_discwind(child,grid,this,limitscalar)      
   END DO
  END SUBROUTINE add_discwind

  RECURSIVE SUBROUTINE addDisc(thisOctal,grid)
    ! uses an external function to decide whether to split a subcell of
    !   the current octal. 

    use inputs_mod, only : maxDepthAMR, height, rinner, router, betaDisc
    use inputs_mod, only : heightSplitFac
    IMPLICIT NONE
    type(romanova) :: romData
    TYPE(OCTAL), POINTER :: thisOctal
    TYPE(gridtype), INTENT(INOUT) :: grid   ! need to pass the grid through to the 
    TYPE(OCTAL), POINTER :: child
    INTEGER              :: i, j, k    ! loop counters
    integer :: iindex
    integer :: isubcell
    logical :: split
    real(double) :: cellSize, r, thisHeightSplitFac, hr
    type(VECTOR) :: cellCentre

    DO iSubcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(isubcell)) cycle

       split = .false.

       cellSize = thisOctal%subcellSize 
       cellCentre = subcellCentre(thisOctal,isubCell)
       r = sqrt(cellcentre%x**2 + cellcentre%y**2)
       if ((r < 2.*autocm/1.d10).or.(r > 12.*autocm/1.d10)) then

          call addDiscDensity(thisOctal, isubcell)

          thisHeightSplitFac = heightSplitFac
          hr = height * (r / (100.d0*autocm/1.d10))**betadisc
             
          if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > thisheightSplitFac)) split = .true.
          
          if ((abs(cellcentre%z)/hr > 2.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.
             
          if (((r-cellsize/2.d0) < rOuter).and. ((r+cellsize/2.d0) > rOuter)) then
             if ((thisOctal%subcellSize/rOuter > 0.01) .and. (abs(cellCentre%z/hr) < 7.d0)) then
                split = .true.
             endif
          endif
          
          if ((r+cellsize/2.d0) < rInner) split = .false.
          if ((r-cellsize/2.d0) > Router) split = .false.
       endif



       IF (split.and.(thisOctal%nDepth < maxDepthAMR)) then

        CALL addNewChild(thisOctal, iSubcell, grid, adjustGridInfo=.TRUE., &
                         inherit=.false., interp=.false.,splitAzimuthally=.false., romData=romData)

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
       child => thisOctal%child(iIndex)
       CALL addDisc(child,grid)
   END DO
  END SUBROUTINE addDisc



  !
  !  Given a grid object, this routine will add extra grid and assign density (and etc) 
  !  to the grid using the density function defined in this module.
  !  


  subroutine addDiscWind(grid)
    use inputs_mod, only : DW_d, DW_Rmin, DW_Rmax, DW_Tmax, DW_gamma, &
       DW_Mdot, DW_alpha, DW_beta, DW_Rs, DW_f, DW_Twind, limitscalar, ttauriMstar
    real(double) :: DW_Hdisc
    type(GRIDTYPE) :: grid
    type(DISCWIND_type) :: myDiscWind

    call new(mydiscWind, DW_d, DW_Rmin, DW_Rmax, DW_Tmax, DW_gamma, &
       DW_Mdot, DW_alpha, DW_beta, DW_Rs, DW_f, DW_Twind, dble(ttauriMstar)/msol, DW_Hdisc)
    globalDiscWind = myDiscWind
    call add_discwind(grid%octreeRoot, grid, myDiscWind, limitscalar)
    call assignDensitiesDiscwind(grid, grid%octreeRoot, myDiscWind)

  end subroutine addDiscWind



  !
  ! Split using the log-scaled radial grid.
  logical function need_to_split2(thisOctal,subcell, this)
    use utils_mod, only: locate 

    IMPLICIT NONE

    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    type(discwind_type), intent(in) :: this
    
    real(oct)  :: cellSize, d
    TYPE(VECTOR)     :: cellCentre
!    integer, parameter :: nr = 150  ! normal resolution
!    integer, parameter :: nr = 180  ! normal resolution
    integer, parameter :: nr = 40  ! low resolution

    real(double) :: r
    integer :: i
    !
    logical, save :: first_time = .true.
    real(double) , save:: rGrid(nr)
!    real(double) :: Rmax = 1.5d3  ! [10^10cm] = 1 AU
    TYPE(VECTOR)     :: VecInnerEdge
    real(double) :: wi

    need_to_split2 = .false.

    if (first_time) then  ! setup ref. r grid
       do i = 1, nr
          ! this should not depend on the size of the model boundary box.
!          rGrid(i) = log10(this%Rmin)+dble(i-1)/dble(nr-1)*(log10(Rmax)-log10(this%Rmin))
          rGrid(i) = log10(this%Rmin)+dble(i-1)/dble(nr-1)*(log10(10.d0*autocm/1.d10)-log10(this%rMin))
       enddo
       do i = 1, nr
          rGrid(i) = 10.d0**rGrid(i)
       end do
       first_time = .false.
    end if

    cellSize = (thisOctal%subcellSize)*2.0d0
    cellCentre = subcellCentre(thisOctal, subcell)


!    if (.not.in_jet_flow(this, cellCentre) ) then
    if (.not.in_discwind(this, cellCentre%x, cellCentre%y, cellCentre%z, thisOctal%subcellSize) ) then
       need_to_split2 = .false.
    else
       ! get the size and the position of the centre of the current cell
!       r = modulus(cellCentre)  ! used in the paper
       wi = sqrt(cellCentre%x*cellCentre%x + cellCentre%y*cellCentre%y)
       VecInnerEdge = VECTOR(cellCentre%x, cellCentre%y, 0.0d0)* (this%Rmin/wi)
       r = modulus(cellCentre-VecInnerEdge)  ! shift it to the inner edge of the disc
!       r = ABS(wi-this%Rmin)  ! just a cylindical radius
       call locate(rGrid,nr,r,i)
       if (i > (nr-1)) i = nr-1
       d = rGrid(i+1) - rGrid(i)
       if (cellSize > d ) then
          need_to_split2 = .true.
       end if
    endif
  end function need_to_split2




  SUBROUTINE add_new_children_discwind(parent, ichild, grid, this, splitAzimuthally)
    use memory_mod
    use inputs_mod, only : cylindrical!, vturb
    ! adds all eight new children to an octal
    IMPLICIT NONE
    logical, optional :: splitAzimuthally
    type(VECTOR) :: rVec
    TYPE(octal), POINTER :: parent     ! pointer to the parent octal 
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid to routines that
                                          !   calculate the variables stored in
                                          !   the tree.
    TYPE(discwind_type), INTENT(IN)  :: this
    INTEGER       :: subcell           ! loop counter
                                       ! - this isn't very clever. might change it. 
    INTEGER       :: newChildIndex     ! the storage location for the new child
    INTEGER       :: iChild
    integer :: parentSubcell
    integer :: nChildren
    TYPE(OCTAL), POINTER :: thisOctal

    nChildren = parent%nChildren

    parentSubcell = iChild

    ! safety checks of child array
    IF ( ASSOCIATED(parent%child) ) THEN
      IF ( ( nChildren == 0 ) .OR.                  &
           ( nChildren /= SIZE(parent%child) ) ) THEN
        PRINT *, 'Panic: in addNewChild, %child array wrong size'
        PRINT *, 'nChildren:',nChildren,' SIZE %child:', SIZE(parent%child)
        STOP
      END IF
    END IF
    IF ( (.NOT. ASSOCIATED(parent%child)) .AND. (nChildren > 0) ) THEN
      PRINT *, 'Panic: in addNewChild, %child array wrong size'
      PRINT *, 'nChildren:',nChildren,' ASSOCIATED %child:', ASSOCIATED(parent%child)
      STOP
    END IF

    ! check that new child does not already exist
    IF ( parent%hasChild(iChild) .EQV. .TRUE. ) THEN
      PRINT *, 'Panic: in addNewChild, attempted to add a child ',&
               '       that already exists'
      STOP
    ENDIF

!    call interpFromParent(subcellCentre(parent, iChild, parent%subcellSize, &
!         grid, temperature, density, dusttypeFraction)

    CALL growChildArray(parent, nNewChildren=1, grid=grid )

    ! update the bookkeeping
    nChildren = nChildren + 1
    newChildIndex = nChildren

    ! update the parent octal
    parent%nChildren = nChildren
    parent%hasChild(iChild) = .TRUE.
    parent%indexChild(newChildIndex) = iChild

    NULLIFY(parent%child(newChildIndex)%child)

    parent%child(newChildIndex)%nDepth = parent%nDepth + 1
    if (parent%child(newChildIndex)%nDepth  > grid%maxDepth) then
       grid%maxDepth = grid%maxDepth + 1
       CALL setSmallestSubcell(grid)
    endif
    ! set up the new child's variables
    parent%child(newChildIndex)%threeD = parent%threeD
    parent%child(newChildIndex)%twoD = parent%twoD
    parent%child(newChildIndex)%oneD = parent%oneD
    parent%child(newChildIndex)%maxChildren = parent%maxChildren
    parent%child(newChildIndex)%cylindrical = parent%cylindrical

    if (cylindrical) then  
       if (parent%splitAzimuthally) then
          rVec =  subcellCentre(parent,iChild)
          parent%child(newChildIndex)%phi = atan2(rvec%y, rVec%x)
          if (parent%child(newChildIndex)%phi < 0.d0) then
              parent%child(newChildIndex)%phi = parent%child(newChildIndex)%phi + twoPi 
          endif
          parent%child(newChildIndex)%dphi = parent%dphi/2.d0
          parent%child(newChildIndex)%phimin = parent%child(newChildIndex)%phi - parent%dPhi/4.d0
          parent%child(newChildIndex)%phimax = parent%child(newChildIndex)%phi + parent%dPhi/4.d0
       else
          parent%child(newChildIndex)%phi = parent%phi
          parent%child(newChildIndex)%dphi = parent%dphi
          parent%child(newChildIndex)%phimin = parent%phimin
          parent%child(newChildIndex)%phimax = parent%phimax
       endif
       if (parent%child(newChildIndex)%phimin < 1.d-10) parent%child(newChildIndex)%phimin = 0.d0 ! fixed!!!!
       parent%child(newChildIndex)%splitAzimuthally = .false.
       parent%child(newChildIndex)%maxChildren = 4

       if (PRESENT(splitAzimuthally)) then
          if (splitAzimuthally) then
             parent%child(newChildIndex)%splitAzimuthally = .true.
             parent%child(newChildIndex)%maxChildren = 8
             rVec =  subcellCentre(parent,iChild)
             parent%child(newChildIndex)%phi = atan2(rvec%y, rVec%x)
             if (parent%child(newChildIndex)%phi < 0.d0) then
                parent%child(newChildIndex)%phi = parent%child(newChildIndex)%phi + twoPi 
             endif
             if (parent%splitAzimuthally) then
                parent%child(newChildIndex)%dphi = parent%dphi/2.d0
             else
                parent%child(newChildIndex)%dphi = parent%dphi
             endif
          else
             if (parent%splitAzimuthally) then
                rVec =  subcellCentre(parent,iChild)
                parent%child(newChildIndex)%phi = atan2(rvec%y, rVec%x)
                if (parent%child(newChildIndex)%phi < 0.d0) then
                   parent%child(newChildIndex)%phi = parent%child(newChildIndex)%phi + twoPi 
                endif
                parent%child(newChildIndex)%dphi = parent%dphi/2.d0
             else
                parent%child(newChildIndex)%phi = parent%phi
                parent%child(newChildIndex)%dphi = parent%dphi
             endif
             parent%child(newChildIndex)%splitAzimuthally = .false.
             parent%child(newChildIndex)%maxChildren = 4
          endif
       endif
    endif
    

    parent%child(newChildIndex)%inFlow = parent%inFlow
    parent%child(newChildIndex)%parent => parent
    parent%child(newChildIndex)%parentSubcell = iChild
    parent%child(newChildIndex)%subcellSize = parent%subcellSize / 2.0_oc
    parent%child(newChildIndex)%hasChild = .false.
    parent%child(newChildIndex)%nChildren = 0
    parent%child(newChildIndex)%indexChild = -999 ! values are undefined
    parent%child(newChildIndex)%centre = subcellCentre(parent,iChild)
    if (parent%cylindrical) then
       parent%child(newChildIndex)%r = subcellRadius(parent,iChild)
    endif


    parent%child(newChildIndex)%xMin = parent%child(newChildIndex)%centre%x - parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%yMin = parent%child(newChildIndex)%centre%y - parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%zMin = parent%child(newChildIndex)%centre%z - parent%child(newChildIndex)%subcellSize

    parent%child(newChildIndex)%xMax = parent%child(newChildIndex)%centre%x + parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%yMax = parent%child(newChildIndex)%centre%y + parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%zMax = parent%child(newChildIndex)%centre%z + parent%child(newChildIndex)%subcellSize


    thisOctal => parent%child(newChildIndex)
    call allocateOctalAttributes(grid, thisOctal)

    globalMemoryFootprint = globalMemoryFootprint + octalMemory(thisOctal)

    call fill_velocity_corners(this, thisOctal)


    do subcell = 1, thisOctal%maxChildren
       thisOctal%temperature(subcell) = parent%temperature(parentSubcell)
       thisOctal%rho(subcell) = parent%rho(parentSubcell)
       thisOctal%velocity(subcell) = parent%velocity(parentSubcell)

!       rVec = subcellCentre(thisOctal, subcell)
!       x = rVec%x
!       y = rVec%y
!       z = rVec%z
!       if (all_in_discwind(thisOctal, subcell, this)) then
!          thisOctal%inFlow(subcell) = .true.
!          thisOctal%temperature(subcell) = real(this%Twind)
!          thisOctal%rho(subcell) = ave_discwind_density(thisOctal, subcell, this)
!          thisOctal%velocity(subcell) = discwind_velocity(this, vector(x,y,z))
!          if (associated(thisOctal%microturb)) thisOctal%microturb(subcell) = vturb
!       endif
    enddo
  end SUBROUTINE add_new_children_discwind


  SUBROUTINE calcValuesAMR(thisOctal,subcell,grid, inherit, interp, &
       romData, stream)
    ! calculates the variables describing one subcell of an octal.
    ! each geometry that can be used with AMR should be described here, 
    !   otherwise the program will print a warning and exit.
    use inputs_mod, only : geometry
    use luc_cir3d_class, only: calc_cir3d_mass_velocity
    use cmfgen_class, only: cmfgen_mass_velocity
    use jets_mod, only: calcJetsMassVelocity
    use romanova_class, only: calc_romanova_mass_velocity
    use vh1_mod, only: assign_from_vh1, vh1FileRequired
    use gridFromFlash, only: assign_from_flash, flashFileRequired
#ifdef USECFITSIO
    use gridFromFitsFile, only: assign_from_fitsfile
#endif
    use angularImage_utils, only: calcAngImgTest

    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT)    :: thisOctal   ! the octal being changed
    INTEGER, INTENT(IN)           :: subcell   ! the subcell being changed
    TYPE(gridtype), INTENT(INOUT) :: grid      ! the grid
    !
    type(STREAMTYPE), optional :: stream
    LOGICAL, OPTIONAL :: inherit               ! inherit densities, temp, etc of parent
    LOGICAL, OPTIONAL :: interp                ! interpolate densities, temp, etc of parent
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry
    !
    TYPE(octal), POINTER :: parentOctal
!    real(double) :: rhoDouble, r
    INTEGER :: parentSubcell
!    type(VECTOR) :: rVec
    LOGICAL :: inheritProps
    LOGICAL :: interpolate


    inheritProps = .false.
    if (present(inherit)) then
       inheritProps = inherit
    endif

    interpolate = .false.
    if (present(interp)) then
       interpolate = interp
    endif
       parentOctal => thisOctal%parent
       parentSubcell = thisOctal%parentSubcell

    if (inheritProps) then
       if (associated (thisOctal%boundaryCondition)) &
            thisOCtal%boundaryCondition(subcell) = parentOctal%boundaryCondition(parentSubcell)
       if (associated (thisOctal%boundaryCell)) &
            thisOCtal%boundaryCell(subcell) = parentOctal%boundaryCell(parentSubcell)
       if (associated(thisOctal%gamma)) thisOctal%gamma(subcell)  = parentOctal%gamma(parentSubcell)
       if (associated(thisOctal%iEquationOfState)) &
            thisOctal%iEquationOfState(subcell)  = parentOctal%iEquationOfState(parentSubcell)

       if (associated(thisOctal%iAnalyticalVelocity)) &
            thisOctal%iAnalyticalVelocity(subcell)  = parentOctal%iAnalyticalVelocity(parentSubcell)

       thisOctal%rho(subcell) = parentOctal%rho(parentSubcell)
       thisOctal%temperature(subcell) = parentOctal%temperature(parentSubcell)
       if (associated(thisOctal%etaCont)) thisOctal%etaCont(subcell) = parentOctal%etaCont(parentSubcell)
       thisOctal%inFlow(subcell) = parentOctal%inFlow(parentSubcell)
       thisOctal%velocity(subcell) = parentOctal%velocity(parentSubcell)
       if (associated(thisOctal%biasCont3d)) thisOctal%biasCont3D(subcell) = parentOctal%biasCont3D(parentSubcell)
       if (associated(thisOctal%etaLine)) thisOctal%etaLine(subcell) = parentOctal%etaLine(parentSubcell)
       if (associated(thisOctal%chiLine)) thisOctal%chiLine(subcell) = parentOctal%chiLine(parentSubcell)
       if (associated(thisOctal%boundaryCondition)) thisOCtal%boundaryCondition(subcell) &
            = parentOctal%boundaryCondition(parentSubcell)
       if (associated(thisOctal%dustTypeFraction)) then
          thisOctal%dustTypeFraction(subcell,:) = parentOctal%dustTypeFraction(parentSubcell,:)
       endif
       if (associated(thisOctal%origdustTypeFraction)) then
          thisOctal%origdustTypeFraction(subcell,:) = parentOctal%origdustTypeFraction(parentSubcell,:)
       endif
       if (associated(thisOctal%oldFrac)) thisOctal%oldFrac(subcell) = parentOctal%oldFrac(parentSubcell)
       if (associated(thisOctal%ionFrac)) then
          thisOctal%ionFrac(subcell,:) = parentOctal%ionFrac(parentsubcell,:)
       endif
       if (associated(thisOctal%nh)) thisOctal%nh(subcell) = parentOctal%nh(parentsubcell)
       if (associated(thisOctal%ne)) thisOctal%ne(subcell) = parentOctal%ne(parentsubcell)

       if (associated(thisOctal%divV)) thisOctal%divV(subcell) = parentOctal%divV(parentSubcell)
       if (associated(thisOctal%rhou)) thisOctal%rhou(subcell) = parentOctal%rhou(parentSubcell)
       if (associated(thisOctal%rhov)) thisOctal%rhov(subcell) = parentOctal%rhov(parentSubcell)
       if (associated(thisOctal%rhow)) thisOctal%rhow(subcell) = parentOctal%rhow(parentSubcell)
       if (associated(thisOctal%rhoe)) thisOctal%rhoe(subcell) = parentOctal%rhoe(parentSubcell)
       if (associated(thisOctal%rhoeLastTime)) thisOctal%rhoeLastTime(subcell) = &
            parentOctal%rhoeLastTime(parentSubcell)
       if (associated(thisOctal%phi_i)) thisOctal%phi_i(subcell) = parentOctal%phi_i(parentSubcell)
       if (associated(thisOctal%phi_gas)) thisOctal%phi_gas(subcell) = parentOctal%phi_gas(parentSubcell)
       if (associated(thisOctal%phi_stars)) thisOctal%phi_stars(subcell) = parentOctal%phi_stars(parentSubcell)
       if (associated(thisOctal%energy)) thisOctal%energy(subcell) = parentOctal%energy(parentSubcell)

    else if (interpolate) then
       if (associated(thisOctal%boundaryCondition)) &
            thisOCtal%boundaryCondition(subcell) = parentOctal%boundaryCondition(parentSubcell)
       if (associated(thisOctal%etaCont)) Thisoctal%etacont(subcell) = parentOctal%etaCont(parentSubcell)
       thisOctal%inFlow(subcell) = parentOctal%inFlow(parentSubcell)
       thisOctal%velocity(subcell) = parentOctal%velocity(parentSubcell)
       if (associated(thisOctal%biasCont3d)) thisOctal%biasCont3D(subcell) = parentOctal%biasCont3D(parentSubcell)
       if (associated(thisOctal%etaline)) thisOctal%etaLine(subcell) = parentOctal%etaLine(parentSubcell)
       if (associated(thisOctal%chiLine)) thisOctal%chiLine(subcell) = parentOctal%chiLine(parentSubcell)
!       thisOctal%dustTypeFraction(subcell,:) = parentOctal%dustTypeFraction(parentSubcell,:)
!       thisOctal%oldFrac(subcell) = parentOctal%oldFrac(parentSubcell)
       if (associated(thisOctal%ionFrac)) then
          thisOctal%ionFrac(subcell,:) = parentOctal%ionFrac(parentsubcell,:)
       endif
       if (associated(thisOctal%nh)) &
            thisOctal%nh(subcell) = parentOctal%nh(parentsubcell)
       if (associated(thisOctal%ne)) &
            thisOctal%ne(subcell) = parentOctal%ne(parentsubcell)

       parentOctal%hasChild(parentsubcell) = .false.
       call interpFromParent(subcellCentre(thisOctal, subcell), thisOctal%subcellSize, grid, &
            thisOctal%temperature(subcell), thisOctal%rho(subcell), thisOctal%dusttypeFraction(subcell, :), &
            thisOctal%etaLine(subcell))
       parentOctal%hasChild(parentsubcell) = .true.

    else

!    write(*,'(a,a,a)') "trim(grid%geometry):",trim(grid%geometry),"*"
    SELECT CASE (geometry)!trim(grid%geometry))

    CASE("pathtest")
       call calcPathTestDensity(thisOctal,subcell)

    CASE ("ttauri")
      CALL calcTTauriMassVelocity(thisOctal,subcell, grid)
      
    CASE ("jets")
      CALL calcJetsMassVelocity(thisOctal,subcell,grid)

    CASE ("luc_cir3d")
      CALL calc_cir3d_mass_velocity(thisOctal,subcell)

    CASE ("cmfgen")
      CALL cmfgen_mass_velocity(thisOctal,subcell)

    CASE ("romanova")
      CALL calc_romanova_mass_velocity(romData, thisOctal,subcell)

    CASE ("testamr")
       CALL calcTestDensity(thisOctal,subcell,grid)

    CASE("toydisc")
       call calcToyDiscDensity(thisOctal, subcell)
!       call calcToyDiscDensity(thisOctal, subcell)

    CASE("lexington", "lexpdr")
       CALL calcLexington(thisOctal, subcell, grid)
       if (thisOctal%nDepth > 1) then
          thisOctal%ionFrac(subcell,:) = parentOctal%ionFrac(parentsubcell,:)
       endif

    CASE("point")
       CALL calcPointSource(thisOctal, subcell)
       if (thisOctal%nDepth > 1) then
          thisOctal%ionFrac(subcell,:) = parentOctal%ionFrac(parentsubcell,:)
       endif

    CASE("imgTest")
       CALL calcCylinderTest(thisOctal, subcell)
       if (thisOctal%nDepth > 1) then
          thisOctal%ionFrac(subcell,:) = parentOctal%ionFrac(parentsubcell,:)
       endif

    CASE("angImgTest")
       call calcAngImgTest(thisOctal, subcell)

    CASE("starburst")
       CALL calcStarburst(thisOctal, subcell)
       if (thisOctal%nDepth > 1) then
          thisOctal%ionFrac(subcell,:) = parentOctal%ionFrac(parentsubcell,:)
       endif

    CASE("symbiotic")
       CALL calcSymbiotic(thisOctal, subcell)
       if (thisOctal%nDepth > 1) then
          thisOctal%ionFrac(subcell,:) = parentOctal%ionFrac(parentsubcell,:)
       endif

    CASE("proto")
       CALL calcProtoDensity(thisOctal,subcell,grid)

    CASE("kleyring")
       CALL calcKleyRingDensity(thisOctal,subcell)


    CASE("wrshell")
       CALL calcWRShellDensity(thisOctal,subcell,grid)


    CASE("lighthouse")
       CALL calcLighthouseDensity(thisOctal,subcell)


    CASE("slab")
       CALL calcSlabDensity(thisOctal,subcell)

    CASE("hydro1d")
       call calcHydro1DDensity(thisOctal, subcell)

    CASE("bowshock")
       call calcBowShock(thisOctal, subcell)

    CASE("fluxTest")
       call assign_gaussian(thisOctal,subcell)

    CASE("diagSod")
       call calcHydro1DDensity(thisOctal, subcell)

    CASE("kelvin")
       call calcKelvinDensity(thisOctal, subcell)

    CASE("rcTest")
       call calcRCTestDensity(thisOctal, subcell)

    CASE("rtaylor")
       call calcRTaylorDensity(thisOctal, subcell)

    CASE("turbulence")
       call calcRadiativeRoundUpDensity(thisOctal, subcell)

    CASE("bonnor")
       call calcBonnorEbertDensity(thisOctal, subcell)

    CASE("bubble")
       call calcBubbleDensity(thisOctal, subcell)

    CASE("SB_CD_1Da")
       call calcContactDiscontinuityOneDDensity(thisOctal, subcell, v1=.true.)

    CASE("SB_CD_1Db")
       call calcContactDiscontinuityOneDDensity(thisOctal, subcell, v1=.false.)

    CASE("SB_CD_2Da")
       call calcContactDiscontinuityTwoDDensity(thisOctal, subcell, v1=.true.)

    CASE("SB_CD_2Db")
       call calcContactDiscontinuityTwoDDensity(thisOctal, subcell, v1=.false.)

    CASE("SB_instblt")
       call calcPlanarIfrontDensity(thisOctal, subcell)

    CASE("SB_gasmix")
       call calcMixingGasDensity(thisOctal, subcell)

    CASE("SB_WNHII")
       call calcWhalenNormanHIIExpansionDensity(thisOctal, subcell)

    CASE("SB_Dtype")
       call DtypeDensity(thisOctal, subcell)

    CASE("SB_offCen")
       call calcOffCentreExpansionDensity(thisOctal, subcell)
       
    CASE("SB_isoshck")
       call calcIsothermalShockDensity(thisOctal, subcell)

    CASE("SB_coolshk")
       call calcCoolingShockDensity(thisOctal, subcell)

    CASE("SB_runaway")
       call calcSB_runaway(thisOctal, subcell)

    CASE("isosphere")
       call calcIsoSphereDensity(thisOctal, subcell)

    CASE("rv1test")
       call calcrv1TestDensity(thisOctal, subcell)

    CASE("rv2test")
       call calcrv2TestDensity(thisOctal, subcell)

    CASE("rv3test")
       call calcrv3TestDensity(thisOctal, subcell)

    CASE("rv4test")
       call calcrv4TestDensity(thisOctal, subcell)

    CASE("brunt")
       call calcBruntDensity(thisOctal, subcell) 

    CASE("blobtest")
       call calcBlobTestDensity(thisOctal, subcell)

    CASE("radcloud")
       call calcRadialClouds(thisOctal, subcell)

    CASE("unisphere")
       call calcUniformsphere(thisOctal, subcell)

    CASE("unimed")
       call calcUniMed(thisOctal, subcell)


    CASE("sphere")
       call calcSphere(thisOctal, subcell)


    CASE("triangle")
       call calcTriangle(thisOctal, subcell)

    CASE("arbitrary")
       call calcArbitrary(thisOctal, subcell)

    CASE("spiral")
       call calcSpiral(thisOctal, subcell)


    CASE("envelope")
       call calcEnvelope(thisOctal, subcell)

    CASE("interptest")
       call calcinterptest(thisOctal, subcell)

    CASE("gravtest")
       call calcGravtest(thisOctal, subcell)

    CASE("maclaurin")
       call maclaurinSpheroid(thisOctal, subcell)

    CASE("sedov")
       call calcSedovDensity(thisOctal, subcell)

    CASE("protobin")
       call calcProtoBinDensity(thisOctal, subcell)

    CASE("nbody")
       CALL calcnBodyDensity(thisOctal,subcell)

    CASE("bondi")
       CALL calcBondiDensity(thisOctal,subcell)

    CASE("bondihoyle")
       CALL calcBondiHoyleDensity(thisOctal,subcell)

    CASE("krumdisc")
       CALL calcKrumholzDiscDensity(thisOctal,subcell)

    CASE("shu")
       CALL calcShuDensity(thisOctal,subcell)

    CASE("empty")
       CALL calcempty(thisOctal,subcell)


    CASE("gammavel")
       CALL calcGammaVel(thisOctal,subcell,grid)

    CASE ("spiralwind")
       CALL spiralWindSubcell(thisOctal, subcell)

    CASE ("wind")
       CALL WindSubcell(thisOctal, subcell)
       
    CASE("sphfile","cluster","molcluster","theGalaxy","wr104", "dale")

       ! Flag as missing data to be filled in by FinishGrid
       thisoctal%rho(subcell) = -9.9d99
       if (associated (thisoctal%cornervelocity)) thisoctal%cornervelocity = VECTOR(-9.9d99,-9.9d99,-9.9d99)

    CASE ("benchmark")
       CALL benchmarkDisk(thisOctal, subcell)

    CASE ("RHDDisc")
       CALL RHDDisc(thisOctal, subcell)

    CASE ("simpledisc")
       CALL simpledisc(thisOctal, subcell)

    CASE ("parker")
       CALL parkerwind(thisOctal, subcell)

    CASE ("fontdisc")
       CALL fontdisc(thisOctal, subcell)

    CASE ("molebench")
       thisoctal%rho(subcell) = -9.9d99
       CALL molecularBenchmark(thisOctal, subcell)


    CASE ("h2obench1")
       CALL WaterBenchmark1(thisOctal, subcell)

    CASE ("h2obench2")
       CALL WaterBenchmark2(thisOctal, subcell)

    CASE ("agbstar")
       CALL AGBStarBenchmark(thisOctal, subcell)

    CASE ("clumpyagb")
       CALL clumpyagb(thisOctal, subcell)

    CASE ("filament")
       CALL molecularFilamentFill(thisOctal, subcell)

    CASE ("ggtau")
       CALL ggtauFill(thisOctal, subcell)

    CASE ("shakara","aksco","circumbin")
       CALL shakaraDisk(thisOctal, subcell ,grid)

    CASE ("adddisc")
       CALL addDiscDensity(thisOctal, subcell)

    CASE ("iras04158")
       CALL iras04158(thisOctal, subcell)

    CASE ("turbbox")
       CALL turbBox(thisOctal, subcell)

    CASE ("warpeddisc")
       CALL warpedDisk(thisOctal, subcell ,grid)

    CASE("ppdisk")
       CALL calcPPDiskDensity(thisOctal,subcell,grid)

    CASE("benchi")
       CALL calcBenchI(thisOctal,subcell)


    CASE("melvin")
       CALL assign_melvin(thisOctal,subcell,grid)

    CASE("gaussian")
       CALL assign_gaussian(thisOctal,subcell)

    CASE("hii_test")
       CALL assign_hii_test(thisOctal,subcell)

    CASE("radpress")
       CALL assign_radpresstest(thisOctal,subcell)

    CASE("whitney")
       CALL assign_whitney(thisOctal,subcell,grid)

    CASE("planetgap")
       CALL assign_planetgap(thisOctal,subcell,grid)

    CASE("toruslogo")
       CALL assign_toruslogo(thisOctal,subcell)

    CASE("clumpydisc")
       CALL assign_clumpydisc(thisOctal, subcell, grid)

    CASE ("windtest")
      CALL calcWindTestValues(thisOctal,subcell,grid)

   CASE ("fractal")
      thisOctal%rho = 100.d0 * mHydrogen
      thisOctal%temperature = 8000.

   CASE ("runaway")
      call calcRunaway

   CASE("NLR")
      call calcNarrowLineRegion(thisOctal, subcell)

#ifdef USECFITSIO
   CASE("fitsfile")
      call assign_from_fitsfile(thisOctal, subcell)
#endif

    CASE ("magstream")

!      thisOctal%rho(subcell) = 1.d-30
!      thisOctal%temperature(subcell) = 6000.
!
!       rVec = subcellCentre(thisOctal, subcell)
!       if (iGlobalSample /= 0) then
!          r = modulus(rVec - globalStream%position(iGlobalSample))
!          if (r < globalStream%streamradius(iGlobalSample)) then
!             thisOctal%rho(subcell) = globalStream%rho(iGlobalSample)
!             write(*,*) thisOctal%rho(subcell),iglobalsample
!             thisOctal%temperature(subcell) = 7500.
!          endif
!       endif


!      CALL getMagStreamValues2(point=subcellCentre(thisOctal,subcell),&
!                              grid=grid,                            &
!                              rho=rhoDouble,                        &
!                              temperature=thisOctal%temperature(subcell),&
!                              velocity=thisOctal%velocity(subcell),  &
!                              inFlow=thisOctal%inFlow(subcell))

      call getMagStreamValues3(thisOctal, subcell, stream)


!      thisOctal%rho(subcell) = REAL(rhoDouble)
!      IF (subcell == thisOctal%maxChildren) CALL fillVelocityCorners(thisOctal,grid,magStreamVelocity)
!       thisOctal%microturb(subcell) = (20.d5/cSpeed) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

    CASE DEFAULT
      WRITE(*,*) "! Unrecognised grid geometry in calcValuesAMR: ",TRIM(grid%geometry)
      STOP

    END SELECT
 
!    CALL fillGridDummyValues(thisOctal,subcell, grid)
 
    end if
 
    contains

      subroutine calcRunaway
        use inputs_mod, only: sourcePos
        TYPE(vector) :: sourceToCell, thisVel
        real(db) :: vOutflow
        real(db), parameter :: vWind = 1800.0_db
        logical, parameter :: runawayDustFromOutflow=.false.

        if (vh1FileRequired()) then 
           call assign_from_vh1(thisOctal, subcell)
        elseif(flashFileRequired()) then
           call assign_from_flash(thisOctal, subcell)
        else
           call torus_abort("No way to assign density for runaway geometry")
        endif

! Any temperature set in call to assign_from_* is overwritten here
        thisOctal%temperature(subcell) = 10000.
        thisOctal%etaCont(subcell) = 0.
        thisOctal%nh(subcell)      = thisOctal%rho(subcell) / mHydrogen
        thisOctal%ne(subcell)      = thisOctal%nh(subcell)
        thisOctal%nhi(subcell)     = 1.e-8
        thisOctal%nhii(subcell)    = thisOctal%ne(subcell)
        thisOctal%inFlow(subcell)  = .true.
        thisOctal%biasCont3D       = 1.0
        thisOctal%etaLine          = 1.e-30
        if (thisOctal%nDepth > 1) then
           thisOctal%ionFrac(subcell,:) = parentOctal%ionFrac(parentsubcell,:)
        endif

! This section decides where to put dust
! Vector from first source to cell centre
        sourceToCell = subcellcentre(thisOctal, subcell) - sourcePos(1)
        thisVel = thisOctal%velocity(subcell)
! Project velocity onto sourceToCell vector to get the outflow velocity
        if (runawayDustFromOutflow) then
           vOutflow = (thisVel.dot.sourceToCell)/modulus(sourceToCell)
        else
           vOutflow = modulus(thisVel)
        endif
! Convert to km/s
        vOutflow = vOutflow * (cspeed / 1.0e5_db)
! Set up initial dust distribution with no dust in outflow
        if (vOutflow > vWind) then 
           thisOctal%dustTypeFraction(subcell,:) = 0.0
        else
           thisOctal%dustTypeFraction(subcell,:) = 0.01
        endif

      end subroutine calcRunaway

#ifdef USECFITSIO
      subroutine calcPion(thisOctal, subcell)
        use gridFromFitsFile, only : assign_from_fitsfile_interp
        type(OCTAL) :: thisOctal
        integer :: subcell

        call assign_from_fitsfile_interp(thisOctal, subcell)
      end subroutine calcPion
#endif
  END SUBROUTINE calcValuesAMR

  SUBROUTINE initFirstOctal(grid, centre, size, oned, twod, threed, romData ,&
       stream)
    ! creates the first octal of a new grid (the root of the tree).
    ! this should only be used once; use addNewChild for subsequent
    !  additions.
    use inputs_mod, only : cylindrical
    use memory_mod, only : globalMemoryFootprint, octalMemory
    IMPLICIT NONE
    type(OCTAL), pointer :: thisOctal
    TYPE(gridtype), INTENT(INOUT)    :: grid 
    TYPE(vector), INTENT(IN)    :: centre ! coordinates of the grid centre
    REAL, INTENT(IN)                 :: size 
      ! 'size' should be the vertex length of the cube that contains the whole
      !   of the simulation space, *not* the size of a subcell.
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry
    type(streamtype), optional :: stream
    !
    LOGICAL :: oned, twod, threed  ! true if this is a twoD amr grid
    INTEGER :: subcell ! loop counter 
    integer :: i

#ifdef SPH
    character(len=100) :: message
#endif

    ALLOCATE(grid%octreeRoot)

    do i = 1, 8
       grid%octreeRoot%mpiThread(i) = i 
    enddo
    
    ! allocate any variables that need to be 

    grid%octreeRoot%cylindrical = .false.

    if (oneD) then
       grid%octreeRoot%oneD = .true.
       grid%octreeRoot%twoD = .false.
       grid%octreeRoot%threeD = .false.
       grid%octreeRoot%maxChildren = 2
    endif

    if (twoD) then
       grid%octreeRoot%oneD = .false.
       grid%octreeRoot%twoD = .true.
       grid%octreeRoot%threeD = .false.
       grid%octreeRoot%maxChildren = 4
    else if (threed) then
       grid%octreeRoot%oneD = .false.
       grid%octreeRoot%twoD = .false.
       grid%octreeRoot%threeD = .true.
       grid%octreeRoot%maxChildren = 8
       if (cylindrical) then
          grid%octreeRoot%oneD = .false.
          grid%octreeRoot%twoD = .false.
          grid%octreeRoot%threeD = .true.
          grid%octreeRoot%cylindrical = .true.
          grid%octreeRoot%splitAzimuthally = .false.
          grid%octreeRoot%maxChildren = 4
          grid%octreeRoot%phi = pi
          grid%octreeRoot%dPhi = twoPi
          grid%octreeRoot%phiMin = 0.d0
          grid%octreeRoot%phiMax = twoPi
       endif
    endif


    grid%octreeRoot%nDepth = 1
    grid%octreeRoot%nChildren = 0
    grid%octreeRoot%hasChild = .FALSE.
    grid%octreeRoot%subcellSize = dble(size)/2.0_oc
    grid%octreeRoot%centre = centre
    if (cylindrical) then
       grid%octreeRoot%centre%x = -size / 2.d0
       grid%octreeRoot%r        = size / 2.d0
       grid%octreeRoot%centre%y = 0.d0
       grid%octreeRoot%centre%z = 0.d0
    endif
    grid%octreeRoot%xMin = grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize
    grid%octreeRoot%yMin = grid%octreeRoot%centre%y - grid%octreeRoot%subcellSize
    grid%octreeRoot%zMin = grid%octreeRoot%centre%z - grid%octreeRoot%subcellSize
    grid%octreeRoot%xMax = grid%octreeRoot%centre%x + grid%octreeRoot%subcellSize
    grid%octreeRoot%yMax = grid%octreeRoot%centre%y + grid%octreeRoot%subcellSize
    grid%octreeRoot%zMax = grid%octreeRoot%centre%z + grid%octreeRoot%subcellSize
    grid%octreeRoot%indexChild = -999 ! values are undefined
    NULLIFY(grid%octreeRoot%parent)   ! tree root does not have a parent
    NULLIFY(grid%octreeRoot%child)    ! tree root does not yet have children
    ! initialize some values
    grid%octreeRoot%rho = 1.e-30

    thisOctal => grid%octreeRoot
    call allocateOctalAttributes(grid, thisOctal)

    globalMemoryFootprint = globalMemoryFootprint + octalMemory(thisOctal)

    select case (grid%geometry)
#ifdef SPH
       case("sphfile","cluster","molcluster","theGalaxy","wr104", "dale")
          ! Initially we copy the idecies of particles (in SPH data)
          ! to the root node. The indecies will copy over to
          ! to the subcells if the particles are in the subcells.
          ! This will allow us to work with the subsets of gas particle
          ! list hence reduces the computation time when we are 
          ! splitting/constructing the octree data structure. 
          !
          ! Using the routine in grid_mod.f90

          write(message,*) "Copying SPH index to root"
          call writeinfo(message, TRIVIAL)
          call copy_sph_index_to_root(grid)
          write(message,*) "Initializing OctreeRoot"
          call writeinfo(message, TRIVIAL)
          !
          DO subcell = 1, grid%octreeRoot%maxChildren
            ! calculate the values at the centre of each of the subcells
             CALL calcValuesAMR(grid%octreeRoot,subcell,grid)
            ! label the subcells
             grid%octreeRoot%label(subcell) = subcell
          END DO

          write(message,*) "Done."
          call writeinfo(message, TRIVIAL)
#endif
       case("romanova")
          ! just checking..
          if (.not. PRESENT(romData)) then
             print *, "Error:: romData is not present in initFirstOctal!"
             stop
          end if
          DO subcell = 1, grid%octreeRoot%maxChildren
             ! calculate the values at the centre of each of the subcells
             CALL calcValuesAMR(grid%octreeRoot,subcell,grid, romData=romData)
             ! label the subcells
             grid%octreeRoot%label(subcell) = subcell
          END DO

       case("magstream")
          ! just checking..
          if (.not. PRESENT(stream)) then
             print *, "Error:: stream is not present in initFirstOctal!"
             stop
          end if
          DO subcell = 1, grid%octreeRoot%maxChildren
             ! calculate the values at the centre of each of the subcells
             CALL calcValuesAMR(grid%octreeRoot,subcell,grid, stream=stream)
             ! label the subcells
             grid%octreeRoot%label(subcell) = subcell
          END DO
       case DEFAULT

!!!!!!!!!! EDITTED OUT BY TJH
!          DO subcell = 1, grid%octreeRoot%maxChildren
            ! calculate the values at the centre of each of the subcells
!             CALL calcValuesAMR(grid%octreeRoot,subcell,grid)
            ! label the subcells
!             grid%octreeRoot%label(subcell) = subcell
!          END DO
       end select
    


    ! we keep track of the maximum depth of the grid...
    grid%maxDepth = 1
    ! ...and the size of the smallest subcells in the grid.
    ! we will actually store the value which is half the size of the 
    !   smallest subcell because this is more useful for later
    !   calculations.
    grid%halfSmallestSubcell = grid%octreeRoot%subcellSize / &
                                2.0_oc**REAL(grid%maxDepth,kind=oct)

  END SUBROUTINE initFirstOctal


  SUBROUTINE addNewChild(parent, iChild, grid, adjustGridInfo, &
                         inherit, interp, amrHydroInterp, splitAzimuthally, romData, &
                         stream)
    ! adds one new child to an octal

    USE inputs_mod, ONLY : cylindrical, maxMemoryAvailable
!    use mpi_global_mod, only: nThreadsGlobal
    use octal_mod, only: subcellRadius
    use memory_mod, only : octalMemory, globalMemoryFootprint, humanReadableMemory, globalMemoryChecking
#ifdef SPH
    USE cluster_class, only: update_particle_list
    USE sph_data_class, only: isAlive
#endif

    IMPLICIT NONE
    
    TYPE(octal), TARGET, INTENT(INOUT) :: parent ! the parent octal 
    type(octal), pointer :: thisOctal
    INTEGER, INTENT(IN)  :: iChild     ! the label (1-8) of the subcell gaining the child 
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid to routines that
                                          !   calculate the variables stored in
                                          !   the tree.
    type(STREAMTYPE), optional :: stream
    LOGICAL, INTENT(IN) :: adjustGridInfo
    LOGICAL, optional :: splitAzimuthally
      ! whether these variables should be updated: 
      !   grid%nOctals, grid%maxDepth, grid%halfSmallestSubcell    
    LOGICAL, OPTIONAL :: inherit       ! inherit densities, temps, etc from parent
    LOGICAL, OPTIONAL :: interp        ! interpolate densities, temps, etc from parent    
    
    INTEGER       :: subcell           ! loop counter
    INTEGER, SAVE :: counter = 9       ! keeps track of the current subcell label
                                       ! - this isn't very clever. might change it. 
    INTEGER       :: nChildren         ! number of children the parent octal has
    INTEGER       :: newChildIndex     ! the storage location for the new child
    integer :: i
    logical :: inheritProps, interpolate
    logical, optional, intent(in) :: amrHydroInterp
    logical :: doAmrHydroInterp
!   logical, save :: firstTimeMem = .true.
    character(len=80) :: message
    type(VECTOR) :: rVec
    ! array of octals that may be needed for temporarily storing child octals

    ! For "romanova" geometry
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry
   
    if (globalMemoryChecking.and.(globalMemoryFootprint > maxMemoryAvailable)) then
       write(message,'(a)') "Child added while memory exceeded for grid :"//humanReadableMemory(globalMemoryFootprint)
       call writeWarning(message)
    endif

    doAMRhydroInterp = .false.
    if (PRESENT(amrHydroInterp)) doAMRhydroInterp = amrHydroInterp

    inheritProps = .false.
    if (present(inherit)) then
       inheritProps = inherit
    endif

    interpolate = .false.
    if (present(interp)) then
       interpolate = interp
    endif


    ! store the number of children that already exist
    nChildren = parent%nChildren

    ! safety checks of child array
    IF ( ASSOCIATED(parent%child) ) THEN
      IF ( ( nChildren == 0 ) .OR.                  &
           ( nChildren /= SIZE(parent%child) ) ) THEN
        PRINT *, 'Panic: in addNewChild, %child array wrong size'
        PRINT *, 'nChildren:',nChildren,' SIZE %child:', SIZE(parent%child)
        STOP
      END IF
    END IF
    IF ( (.NOT. ASSOCIATED(parent%child)) .AND. (nChildren > 0) ) THEN
      PRINT *, 'Panic: in addNewChild, %child array wrong size'
      PRINT *, 'nChildren:',nChildren,' ASSOCIATED %child:', ASSOCIATED(parent%child)
      STOP
    END IF

    ! check that new child does not already exist
    IF ( parent%hasChild(iChild) .EQV. .TRUE. ) THEN
      PRINT *, 'Panic: in addNewChild, attempted to add a child ',&
               '       that already exists'
      STOP
    ENDIF

!    call interpFromParent(subcellCentre(parent, iChild, parent%subcellSize, &
!         grid, temperature, density, dusttypeFraction)

    CALL growChildArray(parent, nNewChildren=1, grid=grid )

    ! update the bookkeeping
    nChildren = nChildren + 1
    newChildIndex = nChildren

    ! update the parent octal
    parent%nChildren = nChildren
    parent%hasChild(iChild) = .TRUE.
    parent%indexChild(newChildIndex) = iChild

    ! allocate any variables that need to be  
!    IF (.NOT.grid%oneKappa) THEN
!       ! The kappa arrays should be allocated with grid%nopacity instead of grid%nlambda
!       ! because for line calculation, there is only one kappa needed.
!       ! (but grid%nlambda is not 1). If you allocate the arrays with grid%nlambda,
!       ! it will be a huge waste of RAM. ---  (RK) 
!       ALLOCATE(parent%child(newChildIndex)%kappaAbs(8,grid%nopacity))
!       ALLOCATE(parent%child(newChildIndex)%kappaSca(8,grid%nopacity))
!       ! ALLOCATE(parent%child(newChildIndex)%kappaAbs(8,grid%nlambda))
!       ! ALLOCATE(parent%child(newChildIndex)%kappaSca(8,grid%nlambda))
!       parent%child(newChildIndex)%kappaAbs = 1.e-30
!       parent%child(newChildIndex)%kappaSca = 1.e-30
!    ENDIF
    NULLIFY(parent%child(newChildIndex)%child)

    parent%child(newChildIndex)%nDepth = parent%nDepth + 1

! setup mpiThread values

    if ( ((parent%twoD)  .and.((nHydroThreadsGlobal) == 4)) .or. &
         ((parent%threed).and.((nHydroThreadsGlobal) == 8)).or. &
         ((parent%oneD)  .and.((nHydroThreadsGlobal) == 2)) ) then
       parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
    else

       if (parent%child(newChildIndex)%nDepth > 2) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          if (parent%oneD) then
             do i = 1, 2
                parent%child(newChildIndex)%mpiThread(i) = 2 * (parent%mpiThread(iChild) - 1) + i
             enddo
          else if (parent%twoD) then
             do i = 1, 4
                parent%child(newChildIndex)%mpiThread(i) = 4 * (parent%mpiThread(iChild) - 1) + i
             enddo
          else if (parent%Threed) then
             do i = 1, 8
                parent%child(newChildIndex)%mpiThread(i) = 8 * (parent%mpiThread(iChild) - 1) + i
             enddo
          endif
       endif
    endif

    if ((parent%twod).and.nHydrothreadsGlobal == 16) then
       if (parent%child(newChildIndex)%nDepth > 2) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          do i = 1, 4
             parent%child(newChildIndex)%mpiThread(i) = 4 * (parent%mpiThread(iChild) - 1) + i
          enddo
       endif
    endif

    if ((parent%threed).and.nHydroThreadsGlobal == 64) then
       if (parent%child(newChildIndex)%nDepth > 2) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          do i = 1, 8
             parent%child(newChildIndex)%mpiThread(i) = 8 * (parent%mpiThread(iChild) - 1) + i
          enddo
       endif
    endif

    if ((parent%twoD).and.nHydroThreadsGlobal == 64) then
       if (parent%child(newChildIndex)%nDepth > 3) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          do i = 1, 4
             parent%child(newChildIndex)%mpiThread(i) = 4 * (parent%mpiThread(iChild) - 1) + i
          enddo
       endif
    endif

    if ((parent%threed).and.nHydrothreadsGlobal == 512) then
       if (parent%child(newChildIndex)%nDepth > 3) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          do i = 1, 8
             parent%child(newChildIndex)%mpiThread(i) = 8 * (parent%mpiThread(iChild) - 1) + i
          enddo
       endif
    endif

    ! set up the new child's variables
    parent%child(newChildIndex)%threeD = parent%threeD
    parent%child(newChildIndex)%twoD = parent%twoD
    parent%child(newChildIndex)%oneD = parent%oneD
    parent%child(newChildIndex)%maxChildren = parent%maxChildren
    parent%child(newChildIndex)%cylindrical = parent%cylindrical
    
    ! if splitAzimuthally is not present then we assume we are not

    if (cylindrical) then  
       if (parent%splitAzimuthally) then
          rVec =  subcellCentre(parent,iChild)
          parent%child(newChildIndex)%phi = atan2(rvec%y, rVec%x)
          if (parent%child(newChildIndex)%phi < 0.d0) then
              parent%child(newChildIndex)%phi = parent%child(newChildIndex)%phi + twoPi 
          endif
          parent%child(newChildIndex)%dphi = parent%dphi/2.d0
          parent%child(newChildIndex)%phimin = parent%child(newChildIndex)%phi - parent%dPhi/4.d0
          parent%child(newChildIndex)%phimax = parent%child(newChildIndex)%phi + parent%dPhi/4.d0
       else
          parent%child(newChildIndex)%phi = parent%phi
          parent%child(newChildIndex)%dphi = parent%dphi
          parent%child(newChildIndex)%phimin = parent%phimin
          parent%child(newChildIndex)%phimax = parent%phimax
       endif
       if (parent%child(newChildIndex)%phimin < 1.d-10) parent%child(newChildIndex)%phimin = 0.d0 ! fixed!!!!
       parent%child(newChildIndex)%splitAzimuthally = .false.
       parent%child(newChildIndex)%maxChildren = 4

       if (PRESENT(splitAzimuthally)) then
          if (splitAzimuthally) then
             parent%child(newChildIndex)%splitAzimuthally = .true.
             parent%child(newChildIndex)%maxChildren = 8
             rVec =  subcellCentre(parent,iChild)
             parent%child(newChildIndex)%phi = atan2(rvec%y, rVec%x)
             if (parent%child(newChildIndex)%phi < 0.d0) then
                parent%child(newChildIndex)%phi = parent%child(newChildIndex)%phi + twoPi 
             endif
             if (parent%splitAzimuthally) then
                parent%child(newChildIndex)%dphi = parent%dphi/2.d0
             else
                parent%child(newChildIndex)%dphi = parent%dphi
             endif
          else
             if (parent%splitAzimuthally) then
                rVec =  subcellCentre(parent,iChild)
                parent%child(newChildIndex)%phi = atan2(rvec%y, rVec%x)
                if (parent%child(newChildIndex)%phi < 0.d0) then
                   parent%child(newChildIndex)%phi = parent%child(newChildIndex)%phi + twoPi 
                endif
                parent%child(newChildIndex)%dphi = parent%dphi/2.d0
             else
                parent%child(newChildIndex)%phi = parent%phi
                parent%child(newChildIndex)%dphi = parent%dphi
             endif
             parent%child(newChildIndex)%splitAzimuthally = .false.
             parent%child(newChildIndex)%maxChildren = 4
          endif
       endif
    endif

    parent%child(newChildIndex)%inFlow = parent%inFlow
    parent%child(newChildIndex)%parent => parent
    parent%child(newChildIndex)%parentSubcell = iChild
    parent%child(newChildIndex)%subcellSize = parent%subcellSize / 2.0_oc
    parent%child(newChildIndex)%hasChild = .false.
    parent%child(newChildIndex)%nChildren = 0
    parent%child(newChildIndex)%indexChild = -999 ! values are undefined
    parent%child(newChildIndex)%centre = subcellCentre(parent,iChild)
    if (parent%cylindrical) then
       parent%child(newChildIndex)%r = subcellRadius(parent,iChild)
    endif

    parent%child(newChildIndex)%xMin = parent%child(newChildIndex)%centre%x - parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%yMin = parent%child(newChildIndex)%centre%y - parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%zMin = parent%child(newChildIndex)%centre%z - parent%child(newChildIndex)%subcellSize

    parent%child(newChildIndex)%xMax = parent%child(newChildIndex)%centre%x + parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%yMax = parent%child(newChildIndex)%centre%y + parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%zMax = parent%child(newChildIndex)%centre%z + parent%child(newChildIndex)%subcellSize


    thisOctal => parent%child(newChildIndex)
    call allocateOctalAttributes(grid, thisOctal)

    globalMemoryFootprint = globalMemoryFootprint + octalMemory(thisOctal)

#ifdef SPH
    if (isAlive()) then
       ! updates the sph particle list.           
       CALL update_particle_list(parent, iChild, newChildIndex)
    endif
#endif

    ! put some data in the four/eight subcells of the new child


    if (.not.doAMRhydroInterp) then
       DO subcell = 1, parent%child(newChildIndex)%maxChildren
          CALL calcValuesAMR(parent%child(newChildIndex),subcell,grid, &
               inherit=inheritProps, &
               interp=interpolate,  &
               romData=romData, stream=stream)
          parent%child(newChildIndex)%label(subcell) = counter
          counter = counter + 1
       END DO
    else

       thisOctal => parent%child(newChildIndex)
       call AMRHydroInterpFromParent(thisOctal, grid)
       DO subcell = 1, parent%child(newChildIndex)%maxChildren
          parent%child(newChildIndex)%label(subcell) = counter
          counter = counter + 1
       END DO
    endif
    IF ( adjustGridInfo ) THEN
      grid%nOctals = grid%nOctals + 1

      ! check for a new maximum depth 
      IF (parent%child(newChildIndex)%nDepth > grid%maxDepth) THEN
        grid%maxDepth = parent%child(newChildIndex)%nDepth
        CALL setSmallestSubcell(grid)
      END IF
    END IF

    IF (( COUNT(parent%hasChild(:)) /= parent%nChildren ) .OR. &     
        (  SIZE(parent%child(:))    /= parent%nChildren )) THEN
        PRINT *, "Problem in addNewChild"
        PRINT *, "  nchildren: ",parent%nChildren
        PRINT *, "  haschild: ",parent%hasChild
        PRINT *, "  indexchild: ",parent%indexchild
        PRINT *, "  ASSOCIATED(child): ",ASSOCIATED(parent%child)
        PRINT *, "  SIZE(child): ",SIZE(parent%child)
        do;enddo
    END IF

    !CALL checkAMRgrid(grid,checkNoctals=.FALSE.)
  END SUBROUTINE addNewChild
  

  RECURSIVE SUBROUTINE splitGrid(thisOctal,amrLimitScalar,amrLimitScalar2,grid, wvars,&
       setChanged, romData)
    ! uses an external function to decide whether to split a subcell of
    !   the current octal. 
!    use inputs_mod, only : splitOverMPI
    IMPLICIT NONE


    TYPE(OCTAL), intent(inout) :: thisOctal
    
    real(double), INTENT(IN) :: amrLimitScalar, amrLimitScalar2 
      ! 'limitScalar' is the value the decideSplit function uses to
      !   decide whether or not to split cell.
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid through to the 
                                          !   routines that this subroutine calls
    LOGICAL, INTENT(IN) :: wvars !refinegrid on hydro variables
    LOGICAL, INTENT(IN), OPTIONAL :: setChanged
    !
    INTEGER              :: iSubcell, iIndex ! loop counters
    INTEGER              :: i, j, k
    logical :: splitInAzimuth

    !
    ! For "romanova" geometry
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry

!THAW - this line was causing cells to thing they were too big in distanceToCellBoundary at 64 way domain decomposition
!    if ((splitOverMPI).and.(myrankGlobal == 0)) goto 666

    DO iSubcell = 1, thisOctal%maxChildren

      IF (decideSplit(thisOctal,iSubcell,amrLimitScalar,amrLimitScalar2,grid, wvars,&
            splitInAzimuth, &
            romData=romData)) THEN

        CALL addNewChild(thisOctal, iSubcell, grid, adjustGridInfo=.TRUE., &
                         splitAzimuthally=splitInAzimuth, romData=romData)
       
        if (.not.thisOctal%hasChild(isubcell)) then
          write(*,*) "add child failed in splitGrid"
          do
          enddo
       endif
        

       IF (PRESENT(setChanged)) THEN
          IF (setChanged) THEN
             
             DO iIndex = 1, thisOctal%nChildren, 1
                IF (thisOctal%indexChild(iIndex) == iSubcell) &
                     thisOctal%child(iIndex)%changed = .TRUE.
             END DO
             
          END IF
       END IF
        
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
      CALL splitGrid(thisOctal%child(iIndex),amrLimitScalar,amrLimitScalar2,grid, wvars, &
                     setChanged, romData=romData)
      
   END DO
!666  continue
  END SUBROUTINE splitGrid
  
  
  RECURSIVE SUBROUTINE getOctalArray(thisOctal,array,counter,maxOctals) 
    ! returns an array of pointers to all of the subcells in the grid.
    ! NB because fortran cannot create arrays of pointers, the output
    !   array is actually of a derived type which *contains* the 
    !   pointer to an octal.
    ! counter should be set to 0 before this routine is called

    IMPLICIT NONE

    TYPE(octal), POINTER                            :: thisOctal
    TYPE(octalWrapper), DIMENSION(:), INTENT(INOUT) :: array 
    INTEGER, INTENT(INOUT)                          :: counter 
    INTEGER, INTENT(IN), OPTIONAL                   :: maxOctals
  
    INTEGER              :: i
    TYPE(octal), POINTER :: child

    ! if this is the root of the tree, we initialize the counter
    IF (.NOT. ASSOCIATED(thisOctal%parent)) counter = 0
    
    counter = counter + 1 
    array(counter)%content => thisOctal
    !array(counter)%inUse = .TRUE. 
    array(counter)%inUse = .NOT. thisOctal%hasChild 
    IF (PRESENT(maxOctals)) THEN
      IF (counter >= maxOctals) RETURN
    END IF
    
    IF ( thisOctal%nChildren > 0 ) THEN
      DO i = 1, thisOctal%nChildren, 1
        
        ! call this subroutine recursively on each of its children
        child => thisOctal%child(i)
        CALL getOctalArray(child,array,counter)
        
      END DO
    END IF

  END SUBROUTINE getOctalArray


  RECURSIVE SUBROUTINE getOctalArrayLevel(thisOctal,array,counter,ndepth,maxOctals) 
    ! returns an array of pointers to all of the subcells in the grid.
    ! NB because fortran cannot create arrays of pointers, the output
    !   array is actually of a derived type which *contains* the 
    !   pointer to an octal.
    ! counter should be set to 0 before this routine is called

    IMPLICIT NONE

    TYPE(octal), POINTER                            :: thisOctal
    TYPE(octalWrapper), DIMENSION(:), INTENT(INOUT) :: array 
    INTEGER, INTENT(INOUT)                          :: counter 
    INTEGER, INTENT(IN), OPTIONAL                   :: maxOctals
    integer :: nDepth
  
    INTEGER              :: i
    TYPE(octal), POINTER :: child

    ! if this is the root of the tree, we initialize the counter
    IF (.NOT. ASSOCIATED(thisOctal%parent)) counter = 0
    
    if ((thisOctal%nDepth == nDepth).or.((thisOctal%nDepth < nDepth).and.(thisOctal%nChildren==0))) then
       counter = counter + 1 
       array(counter)%content => thisOctal
       array(counter)%inUse = .true.
    endif

    IF (PRESENT(maxOctals)) THEN
      IF (counter >= maxOctals) RETURN
    END IF
    
    IF ( thisOctal%nChildren > 0 .and.(thisOctal%nDepth < nDepth)) THEN
      DO i = 1, thisOctal%nChildren, 1
        
        ! call this subroutine recursively on each of its children
        child => thisOctal%child(i)
        CALL getOctalArrayLevel(child,array,counter,nDepth)
        
      END DO
    END IF

  END SUBROUTINE getOctalArrayLevel

  
  SUBROUTINE sortOctalArray(array,grid)
    ! sorts by radius using HeapSort
    
    IMPLICIT NONE

    TYPE(octalWrapper), DIMENSION(:), INTENT(INOUT) :: array 
    TYPE(gridType), INTENT(IN) :: grid
    
    TYPE(octal), POINTER  :: tempContent
    LOGICAL, DIMENSION(8) :: tempInUse
    
    INTEGER               :: i, n

    n = SIZE(array) 
    
    DO i = n/2, 1, -1 
      CALL sift_down(i,n)
    END DO
    
    DO i = n, 2, -1
      ! swap the array elements
      tempContent      => array(1)%content
      tempInUse        =  array(1)%inUse
      array(1)%content => array(i)%content 
      array(1)%inUse   =  array(i)%inUse 
      array(i)%content => tempContent 
      array(i)%inUse   =  tempInUse 
      CALL sift_down(1,i-1)
    END DO

    CONTAINS

    SUBROUTINE sift_down(l,r)
    
      INTEGER, INTENT(IN) :: l, r
      INTEGER :: j, jold 
      real(oct) :: valueX
      real(oct) :: valueY
      real(oct) :: valueA
      real(oct) :: valueJ
      TYPE(octal), POINTER :: Acontent
      LOGICAL, DIMENSION(8) :: AinUse
      TYPE(vector)    :: starPos
      
      starPos = grid%starPos1
      valueA = modulus(array(l)%content%centre - starPos)
      Acontent => array(l)%content
      AinUse   =  array(l)%inUse
      
      jold = l
      j = l + l 
      DO 
        IF (j > r) EXIT 
        
        IF (j < r) THEN
          valueX = modulus(array(j)%content%centre   - starPos)
          valueY = modulus(array(j+1)%content%centre - starPos)
          IF (valueX < valueY) j = j + 1
        END IF 
        
        valueJ = modulus(array(j)%content%centre - starPos)
        IF (valueA >= valueJ) EXIT

        array(jold)%content => array(j)%content
        array(jold)%inUse   =  array(j)%inUse
        
        jold = j
        j = j + j
      END DO
      array(jold)%content => Acontent
      array(jold)%inUse   =  AinUse

    END SUBROUTINE sift_down 

  END SUBROUTINE sortOctalArray    
  
  RECURSIVE SUBROUTINE finishGrid(thisOctal, grid, romData)
    ! takes the octree grid that has been created using 'splitGrid'
    !   and calculates all the other variables in the model.
    ! this should be called once the structure of the grid is complete.
    
    USE inputs_mod, ONLY : modelwashydro, splitOverMPI !, useHartmannTemp
    USE luc_cir3d_class, ONLY:  calc_cir3d_temperature
    USE cmfgen_class, ONLY:     calc_cmfgen_temperature
    USE jets_mod, ONLY:         calcJetsTemperature
    USE romanova_class, ONLY:   calc_romanova_temperature
#ifdef SPH
    USE cluster_class, ONLY:    assign_grid_values 
#endif
    use h21cm_mod

    IMPLICIT NONE

    TYPE(octal), POINTER   :: thisOctal
    TYPE(gridtype)         :: grid
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry

    TYPE(octal), POINTER   :: child
    logical, save :: firstTIme=.true.
    INTEGER :: subcell, iChild

    if (splitOverMPI.and.(myrankGlobal == 0)) goto 666
    ! all of the work that must be done recursively goes here: 
    DO subcell = 1, thisOctal%maxChildren
   

       SELECT CASE (grid%geometry)

       CASE ("ttauri")
!          IF (.NOT. useHartmannTemp) &
!               CALL calcTTauriTemperature(thisOctal,subcell)
          
       CASE ("jets")
          CALL calcJetsTemperature(thisOctal,subcell, grid)
          
       CASE ("luc_cir3d")
          CALL calc_cir3d_temperature(thisOctal,subcell)
        
       CASE ("cmfgen")
         CALL calc_cmfgen_temperature(thisOctal,subcell)

       CASE ("romanova")
          CALL calc_romanova_temperature(romData, thisOctal,subcell)

#ifdef SPH
! Without corner velocities
       CASE ("sphfile","cluster","wr104")
          call assign_grid_values(thisOctal,subcell)

! With corner velocities
       CASE ("molcluster", "theGalaxy", "dale")
          if( .not. thisoctal%haschild(subcell)) then 

             if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
             if (.not. associated(thisoctal%cornervelocity)) then 
                allocate(thisoctal%cornervelocity(27))
                thisoctal%cornervelocity(:) = VECTOR(-9.9d99,-9.9d99,-9.9d99)
             endif
             if(thisoctal%cornervelocity(14)%x .eq. -9.9d99) then
                if(.not. associated(thisoctal%cornerrho)) Allocate(thisOctal%cornerrho(27))
                recentoctal => thisoctal
                CALL fillDensityCorners(thisOctal, clusterdensity, clustervelocity)
                thisOctal%velocity = thisoctal%cornervelocity(14)
             endif
             call assign_grid_values(thisOctal,subcell)
!             write(*,*) "Assigned grid value ",thisOctal%rho(subcell)
          end if
#endif

!          CASE("molebench")
!             CALL fillDensityCorners(thisOctal,grid, molebenchdensity, molebenchvelocity, thisOctal%threed)

       CASE DEFAULT
          ! Nothing to be done for this geometry so just return. 
          if(.not. (modelWasHydro.or.h21cm)) then
             goto 666
          end if
          
       END SELECT
      
       if(modelwashydro) then
          if(firstTime) then
             print *, "Allocating and populating cell corners"
             firstTime = .false.
          end if
          if( .not. thisoctal%haschild(subcell)) then 
             if (.not. associated(thisoctal%cornervelocity)) then 
                allocate(thisoctal%cornervelocity(27))
                thisoctal%cornervelocity(:) = VECTOR(-9.9d99,-9.9d99,-9.9d99)
             endif
             if(thisoctal%cornervelocity(14)%x .eq. -9.9d99) then
                if(.not. associated(thisoctal%cornerrho)) Allocate(thisOctal%cornerrho(27))
                recentoctal => thisoctal
                CALL fillHydroDensityVelocityCorners(thisOctal, grid)
                thisOctal%velocity = thisoctal%cornervelocity(14)
             endif
          end if
       end if

       if (h21cm) then
          call hi_emop(thisOctal%rho(subcell),    thisOctal%temperature(subcell), &
               thisOctal%etaLine(subcell), thisOctal%chiLine(subcell)      )
       endif

    END DO
    
   
    DO iChild = 1, thisOctal%nChildren, 1
       child => thisOctal%child(iChild)
       CALL finishGrid(child, grid, romData=romData)
    END DO

666 continue
    
  END SUBROUTINE finishGrid
 
  
  SUBROUTINE countVoxels(thisOctal,nOctals,nVoxels)  
    ! count the number of octals in the current section of the grid.
    ! also counts the number of unique volume elements (voxels) i.e.
    !   those subcells that are not subdivided

    IMPLICIT NONE

    TYPE(OCTAL), POINTER  :: thisOctal 
    INTEGER,INTENT(INOUT) :: nOctals   ! number of octals
    INTEGER,INTENT(INOUT) :: nVoxels   ! number of childless subcells
    
    nOctals = 0 
    nVoxels = 0
    CALL countVoxelsPrivate(thisOctal)
    
    CONTAINS
    
      RECURSIVE SUBROUTINE countVoxelsPrivate(thisOctal)
      
        TYPE(OCTAL), POINTER  :: thisOctal 
        TYPE(OCTAL), POINTER  :: child
        INTEGER :: i
        
        nOctals = nOctals + 1
        
        IF ( thisOctal%nChildren > 0 ) THEN
          ! call this subroutine recursively on each of its children
          DO i = 1, thisOctal%nChildren, 1
            child => thisOctal%child(i)
            CALL countVoxelsPrivate(child)
          END DO
        END IF

        ! increment the counter once for each of its childless subcells
        nVoxels = nVoxels + (thisOctal%maxChildren - thisOctal%nchildren)
        
      END SUBROUTINE countVoxelsPrivate

  END SUBROUTINE countVoxels


  RECURSIVE SUBROUTINE fixParentPointers(thisOctal)
      
    TYPE(OCTAL), POINTER  :: thisOctal 
    TYPE(OCTAL), POINTER  :: child
    INTEGER :: i
        
    IF ( thisOctal%nChildren > 0 ) THEN
       ! call this subroutine recursively on each of its children
       DO i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          thisOctal%child(i)%parent => thisOctal
          CALL fixParentPointers(child)
       END DO
    END IF
        
  END SUBROUTINE fixParentPointers


  RECURSIVE SUBROUTINE returnSamples (currentPoint,startPoint,locator,direction, &
             octree,grid,sampleFreq,nSamples,maxSamples,abortRay,lambda,usePops, &
             iLambda,error,margin,distanceLimit,kappaAbs,kappaSca,velocity,      &
             velocityDeriv,chiLine,levelPop,rho,temperature, Ne, inFlow,         &
             etaCont,etaLine)
    ! this uses a recursive ray traversal algorithm to sample the octal
    !   grid at points along the path of the ray. 
    ! no checks are made that the ray lies within the boundaries of the
    !  grid, so this subroutine not be called directly. use  
    !  'startReturnSamples' instead.

    IMPLICIT NONE

    ! note: in the code comments, the terms 'subcell' and 'cell' are
    !  mostly interchangeable.

    TYPE(vector), INTENT(INOUT)    :: currentPoint ! current ray position
    TYPE(vector), INTENT(IN)       :: startPoint   ! initial ray position
    TYPE(vector)                   :: locator, rotloc
                  ! 'locator' is used to indicate a point that lies within the  
                  !   *next* cell of the octree that the ray will interesect.
                  !   initially this will be the same as the currentPoint
      
    TYPE(vector), INTENT(IN)       :: direction    ! ray's direction vector 
    TYPE(octal), INTENT(IN)             :: octree       ! an octal grid
    TYPE(gridtype), INTENT(IN)          :: grid         ! grid structure
    REAL, INTENT(IN)                    :: sampleFreq
!    real(oct), INTENT(IN)    :: sampleFreq
                  ! 'sampleFreq' is the maximum number of samples that will be made
                  !   in any subcell of the octree. 
    INTEGER, INTENT(INOUT)              :: nSamples     ! number of samples made
    INTEGER, INTENT(IN)                 :: maxSamples   ! size of sample arrays 
    LOGICAL, INTENT(INOUT)              :: abortRay     ! flag to signal completion
    REAL, DIMENSION(:), INTENT(INOUT)   :: lambda       ! distance travelled by photon
    INTEGER, INTENT(IN)                 :: iLambda      ! wavelength index
    INTEGER, INTENT(INOUT)              :: error        ! error code
    real(oct)                :: margin
                  ! margin is the size of the region around the edge of a subcell
                  !   where numerical inaccuracies may cause problems.
    real(oct), INTENT(IN)    :: distanceLimit ! max length of ray before aborting
    LOGICAL, INTENT(IN)                 :: usePops      ! whether to use level populations
    
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaAbs     ! continuous absorption opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaSca     ! scattering opacities
    TYPE(vector),DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocity ! sampled velocities
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocityDeriv ! sampled velocity derivatives
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL :: chiLine       ! line opacities
    REAL(double),DIMENSION(:),OPTIONAL               :: rho           ! density at sample points
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL   :: temperature! in [K]
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL   :: Ne         ! electron density
    real(double),DIMENSION(:,:),INTENT(INOUT),OPTIONAL :: levelPop   ! level populations 
    logical, DIMENSION(:),INTENT(INOUT),OPTIONAL       :: inFlow     ! inflow flag
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaCont ! contiuum emissivity
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaLine ! line emissivity



    TYPE(vector)      :: exitPoint      ! where ray leaves current cell
    TYPE(vector)      :: centre         ! centre of current subcell
    real(oct)   :: subcellSize    ! size of current subcell
    real(oct)   :: minWallDistance ! distance to *nearest* wall
    
    real(oct)   :: length         ! distance from start of the ray's path.
    real(oct)   :: sampleLength   ! distance interval between samples 
    real(oct)   :: trialLength    ! trial distance to a possible next sample
    TYPE(vector)      :: trialPoint     ! trial location for a next sample 
    INTEGER                :: trial          ! loop counter for trial points 
    INTEGER                :: subcell        ! current subcell 
    INTEGER                :: subIndex       ! octal's child index for current subcell 
    INTEGER                :: i
    real(oct) :: dL


    ! find which of the subcells the point lies in
    DO 
       if (octree%twod) then
          rotloc = projecttoxz(locator)
       else
          rotloc = locator
       endif
       if (octree%oneD) then
          rotloc  = VECTOR(modulus(locator),0.d0,0.d0)
       endif



      subcell = whichSubcell(octree,rotloc)

      ! the IF statement below is for debugging
      !IF ( looseInOctal(octree,locator) .EQV. .FALSE. ) THEN
      !  PRINT *, "Panic: looseInOctal failed in returnSamples"
      !  STOP
      !ENDIF
 
      ! if the subcell has a child, we will find a pointer to that
      !   child and use recursion to sample it.

      IF (octree%hasChild(subcell)) THEN

        ! find pointer to child
        subindex = -99
        DO i = 1, octree%nChildren, 1
          IF ( octree%indexChild(i) == subcell ) THEN
            subIndex = i  
            EXIT
          ENDIF
        ENDDO
        IF (subindex == -99) THEN
          PRINT *, ' Panic: subindex not found'
          STOP
        ENDIF

        CALL returnSamples(currentPoint,startPoint,locator,direction,        &
                    octree%child(subIndex),grid,sampleFreq,nSamples,         &
                    maxSamples,abortRay,lambda,usePops,iLambda,error,margin, &
                    distanceLimit,kappaAbs=kappaAbs,kappaSca=kappaSca,       &
                    velocity=velocity,velocityDeriv=velocityDeriv,           &
                    chiLine=chiLine,levelPop=levelPop,rho=rho,               &
                    temperature=temperature, Ne=Ne, inFlow=inFlow,           &
                    etaCont=etaCont,etaLine=etaLine)
        
        ! after returning from the recursive subroutine, we may have 
        !   finished tracing the ray's path
        IF ( abortRay .EQV. .TRUE. ) RETURN

        ! otherwise, we check whether the ray is still within the
        !   boundaries of the current subsection of the octree.
        IF (.NOT. inOctal(octree,locator)) RETURN
      
      
      ELSE
        ! we now consider the case of the current subcell being childless

        ! find the subcell centre
        centre = subcellCentre(octree,subcell)
        subcellSize = octree%subcellSize

        ! we find the exit point from the current subcell, 
        !  and also 'locator' - a point that lies in the *next* subcell 
        CALL getExitPoint(currentPoint,direction,locator,abortRay,error, &
                          grid%halfSmallestSubcell,exitPoint,centre,subcellSize,&
                          minWallDistance,margin,grid%octreeRoot%threed)

!       call getExitPoint2(currentPoint, direction ,locator, abortRay, error, &
!       grid%halfSmallestSubcell, exitPoint, minWallDistance, grid, octree)

        IF ( abortRay .EQV. .TRUE. ) RETURN

        ! we now decide where we are going to sample the quantities
        ! we will define an approximate rate that is a fraction of
        !   the current subcell size
        sampleLength = octree%subcellSize / sampleFreq

        ! if there are no previously sampled points, we definitely have to take a
        !   sample here
!        IF ( nSamples == 0 ) then
!          length = modulus(currentPoint - startPoint)
!          CALL takeSample(currentPoint,length,direction,grid,octree,subcell,    &
!                          nSamples,maxSamples,usePops,iLambda,error,lambda,     &
!                          kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,&
!                          velocityDeriv=velocityDeriv,chiLine=chiLine,          &
!                          levelPop=levelPop,rho=rho, Ne=Ne)
!        ENDIF

        
        ! we check whether we should take a sample at currentPoint
        !   (which will usually be the entry point to the subcell).

        length = modulus(currentPoint - startPoint)

! force sampling here        
        !IF ( modulus(currentPoint - (startPoint + (direction *         & 
        !                      REAL(lambda(nSamples),kind=oct)))) &
        !                                       > sampleLength ) THEN
          CALL takeSample(currentPoint,length,direction,grid,octree,subcell,    &
                          nSamples,maxSamples,usePops,iLambda,error,lambda,     &
                          kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,&
                          velocityDeriv=velocityDeriv,chiLine=chiLine,          &
                          levelPop=levelPop,rho=rho,temperature=temperature,    &
                          Ne=Ne, inFlow=inFlow, etaCont=etaCont, etaLine=etaLine)
        !END IF

        ! we add sampleLength to the distance from the last location
        !   that was sampled, and decide whether to take a new sample.
        trial = 1
        DO 
          trialPoint = currentPoint + direction * &
                       ((REAL(lambda(nSamples),kind=oct) + &
                       (REAL(trial,kind=oct) * sampleLength)) - length)
                       
          ! we only want to take a sample if we are still within the subcell
          dL = modulus(trialPoint - currentPoint)
          IF (  dL < minWallDistance ) THEN
          
            trialLength = length + dL
            IF (trialLength > distanceLimit) THEN
              abortRay = .TRUE.
              RETURN
            END IF
              
            IF (trialLength >= 0.0_oc) THEN
              CALL takeSample(trialPoint,trialLength,direction,grid,octree,subcell, &
                         nSamples,maxSamples,usePops,iLambda,error,lambda,     &
                         kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,&
                         velocityDeriv=velocityDeriv,chiLine=chiLine,          &
                         levelPop=levelPop,rho=rho,temperature=temperature,    &
                         Ne=Ne, inFlow=inFlow, etaCont=etaCont, etaLine=etaLine)
            ELSE
              EXIT
            END IF
            
          ELSE 
            EXIT
          END IF
          trial = trial + 1
        END DO    

        ! adjust some variables
        length = modulus(exitPoint - startPoint)
        currentPoint = exitPoint
        
      ! if we have left the boundaries of the simulation space, we are finished
      IF (.NOT. inOctal(octree,locator) .OR. &
               length > distanceLimit ) RETURN 
      
    ENDIF
    
  ENDDO


  END SUBROUTINE returnSamples

  
  SUBROUTINE getExitPoint(currentPoint,direction,locator,abortRay,error,    &
                          halfSmallestSubcell,exitPoint,centre,subcellSize, &
                          minWallDistance,margin,threed)
    ! this subroutine finds the closest subcell wall in the direction the
    !   photon is travelling. 
    use inputs_mod, only : amr2d
    IMPLICIT NONE
          
    TYPE(vector), INTENT(IN)       :: currentPoint ! current ray position
    TYPE(vector), INTENT(IN)       :: direction    ! ray's direction vector 
    TYPE(vector), INTENT(OUT)      :: locator       
                  ! 'locator' is used to indicate a point that lies within the  
                  !   *next* cell of the octree that the ray will interesect.
    LOGICAL, INTENT(INOUT)              :: abortRay     ! flag to signal completion
    INTEGER, INTENT(INOUT)              :: error        ! error code
    real(oct), INTENT(IN)    :: halfSmallestSubcell
    real(oct), INTENT(IN)    :: margin
                  ! margin is the size of the region around the edge of a subcell
                  !   where numerical inaccuracies may cause problems.
    TYPE(vector), INTENT(OUT)      :: exitPoint    ! where ray leaves current cell
    TYPE(vector), INTENT(IN)       :: centre       ! centre of current subcell
    real(oct), INTENT(IN)         :: subcellSize  ! size of current subcell
    real(oct), INTENT(OUT)   :: minWallDistance ! distance to *nearest* wall
    
    
    real(oct)   :: wallDistanceX  ! distance to next x-wall-plane 
    real(oct)   :: wallDistanceY  ! distance to next y-wall-plane
    real(oct)   :: wallDistanceZ  ! distance to next z-wall-plane
                  ! 'wallDistanceX,Y,Z' are the distances between the current 
                  !   position and the intersections with the cell walls (along the 
                  !   'direction' vector)
    LOGICAL                :: found          ! status flag  
    real(oct)   :: wallFromOrigin ! distance of cell wall from the origin
    TYPE(vector)      :: wallNormal     ! normal to plane of cell wall

    INTEGER, parameter :: max_num_err = 10;
    INTEGER, save :: num_err = 0    

    LOGICAL :: threed
    REAL(oct) :: r2, d, cosMu, disttor2, disttoxboundary,disttozboundary
    REAL(oct) :: compZ, currentZ, tval, x1, x2
    REAL(oct) :: r1, theta, mu, disttor1

    TYPE(vector) :: xDir, zDir, rVec
    logical :: ok
    ! Specify the ratio of extra length to give it for "locater" to the 
    ! "halfSmallestSubcell" size.
!    REAL(oct), parameter :: frac =1.e-6_oc  
    REAL(oct), parameter :: frac =1.e-2_oc  


    if (threed) then
       
    ! this code is ugly, but should work OK.   
                  

    ! there are six subcell walls - and each wall is part of a plane    
    !   parallel with an axis. we use the direction of the photon to 
    !   find the three planes that the photon will intersect.


    
    IF ( direction%x > 0.0_oc ) THEN
      walldistanceX = (centre%x - currentPoint%x + subcellSize / 2.0_oc  ) / ABS(direction%x) 
    ELSE IF ( direction%x < 0.0_oc ) THEN
      wallDistanceX = (currentPoint%x - centre%x + subcellSize / 2.0_oc ) / ABS(direction%x) 
    ELSE 
      wallDistanceX = HUGE(wallDistanceX) 
    ENDIF

    IF ( direction%y > 0.0_oc ) THEN
      wallDistanceY = (centre%y - currentPoint%y + subcellSize / 2.0_oc ) / ABS(direction%y) 
    ELSE IF ( direction%y < 0.0_oc ) THEN
      wallDistanceY = (currentPoint%y - centre%y + subcellSize / 2.0_oc ) / ABS(direction%y) 
    ELSE 
      wallDistanceY = HUGE(wallDistanceX) 
   ENDIF
        
    IF ( direction%z > 0.0_oc ) THEN
      wallDistanceZ = (centre%z - currentPoint%z + subcellSize / 2.0_oc ) / ABS(direction%z) 
    ELSE IF ( direction%z < 0.0_oc ) THEN
      wallDistanceZ = (currentPoint%z - centre%z + subcellSize / 2.0_oc ) / ABS(direction%z) 
    ELSE 
      wallDistanceZ = HUGE(wallDistanceX) 
   ENDIF

    minWallDistance = MIN(wallDistanceX,wallDistanceY,wallDistanceZ)

    ! we may have problems if "currentPoint" actually lies slightly outside
    !  the current subcell (probably due to numerical accuracy problems). 

    ! if any of the wallDistances are negative or very small, we will 
    !  abandon this photon.

    IF ( wallDistanceX < margin .OR. &
         wallDistanceY < margin .OR. &
         wallDistanceZ < margin      ) THEN 
      abortRay = .TRUE.
      if (num_err < max_num_err) then
         PRINT *, 'In getExitPoint, aborting ray because distance is less than ''margin'''
         PRINT *, '  wallDistances=',wallDistanceX,wallDistanceY,wallDistanceZ
      elseif (num_err ==  max_num_err) then
         PRINT *, ' '
         PRINT *, 'Surpressing the further message from subroutine getExitPoint....'
         PRINT *, ' '
      else
         continue
      endif
      
      error = -20
      num_err = num_err + 1
      
      RETURN      
    ENDIF
       
       
    IF ( wallDistanceX < wallDistanceY .AND. &
         wallDistanceX < wallDistanceZ ) THEN
      wallNormal =  xHat
      IF ( direction%x > 0.0_oc ) THEN
        wallFromOrigin = centre%x + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          PRINT *, 'wallDistances = ',wallDistanceX,&
          wallDistanceY,wallDistanceZ, margin 
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        ENDIF
        locator = exitpoint + (halfSmallestSubcell*xHat)
      ELSE
        wallFromOrigin = centre%x - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        ENDIF
        locator = exitpoint - (halfSmallestSubcell*xHat)
      ENDIF

    ELSEIF ( wallDistanceY < wallDistanceX .AND. &
             wallDistanceY < wallDistanceZ ) THEN
      wallNormal =  yHat
      IF ( direction%y > 0.0_oc ) THEN
        wallFromOrigin = centre%y + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        ENDIF
        locator = exitpoint + (halfSmallestSubcell*yHat)
      ELSE
        wallFromOrigin = centre%y - (subcellSize / 2.0_oc) 
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        ENDIF
        locator = exitpoint - (halfSmallestSubcell*yHat)
      ENDIF
      
    ELSEIF ( wallDistanceZ < wallDistanceX .AND. &
             wallDistanceZ < wallDistanceY ) THEN
      wallNormal =  zHat
      IF ( direction%z > 0.0_oc ) THEN
       wallFromOrigin = centre%z + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        ENDIF
        locator = exitpoint + (halfSmallestSubcell*zHat)
      ELSE
        wallFromOrigin = centre%z - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        ENDIF
        locator = exitpoint - (halfSmallestSubcell*zHat)
      ENDIF

    ! we now consider the case where the ray is leaving through one
    !   of the vertices of the cell.
    ! we could integrate this into the code above, but it is very
    !   unlikely that it will be executed, so we avoid evaluating a 
    !   few IFs by keeping it separate.
    ELSEIF ( wallDistanceX == wallDistanceY .AND. &
             wallDistanceX /= wallDistanceZ ) THEN
      wallNormal =  xHat
      IF ( direction%x > 0.0_oc ) THEN
        wallFromOrigin = centre%x + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        ENDIF
        locator = exitpoint + (halfSmallestSubcell*xHat)
      ELSE
        wallFromOrigin = centre%x - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        ENDIF
        locator = exitpoint - (halfSmallestSubcell*xHat)
      ENDIF
      IF ( direction%y > 0.0_oc ) THEN             
        locator = locator + (halfSmallestSubcell*yHat)
      ELSE  
        locator = locator - (halfSmallestSubcell*yHat)
      ENDIF  
    ELSEIF ( wallDistanceX == wallDistanceZ .AND. &
             wallDistanceX /= wallDistanceY ) THEN
      wallNormal =  xHat
      IF ( direction%x > 0.0_oc ) THEN
        wallFromOrigin = centre%x + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        ENDIF
        locator = exitpoint + (halfSmallestSubcell*xHat)
      ELSE
        wallFromOrigin = centre%x - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
      IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        ENDIF
        locator = exitpoint - (halfSmallestSubcell*xHat)
      ENDIF
      IF ( direction%z > 0.0_oc ) THEN             
        locator = locator + (halfSmallestSubcell*zHat)
      ELSE  
        locator = locator - (halfSmallestSubcell*zHat)
      ENDIF  

    ELSEIF ( wallDistanceY == wallDistanceZ  .AND. &
             wallDistanceX /= wallDistanceZ ) THEN
      wallNormal =  yHat
      IF ( direction%y > 0.0_oc ) THEN
        wallFromOrigin = centre%y + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        ENDIF
        locator = exitpoint + (halfSmallestSubcell*yHat)
      ELSE
        wallFromOrigin = centre%y - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        ENDIF
        locator = exitpoint - (halfSmallestSubcell*yHat)
      ENDIF
      IF ( direction%z > 0.0_oc ) THEN             
        locator = locator + (halfSmallestSubcell*zHat)
      ELSE  
        locator = locator - (halfSmallestSubcell*zHat)
      ENDIF  
      
    ! we now consider the case where the ray is leaving through one
    !   of the corners of the cell.
         
    ELSE IF ( wallDistanceX == wallDistanceY  .AND. &
             wallDistanceY == wallDistanceZ ) THEN
      wallNormal =  xHat
      IF ( direction%x > 0.0_oc ) THEN
        wallFromOrigin = centre%x + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                   (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        ENDIF
        locator = exitpoint + (halfSmallestSubcell*xHat)
      ELSE
        wallFromOrigin = centre%x - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        ENDIF
        locator = exitpoint - (halfSmallestSubcell*xHat)
     ENDIF
      IF ( direction%y > 0.0_oc ) THEN             
        locator = locator + (halfSmallestSubcell*yHat)
      ELSE  
        locator = locator - (halfSmallestSubcell*yHat)
      ENDIF  
      IF ( direction%z > 0.0_oc ) THEN             
        locator = locator + (halfSmallestSubcell*zHat)
      ELSE  
        locator = locator - (halfSmallestSubcell*zHat)
      ENDIF  

    ! we should now have run out of possibilities
    ELSE
      PRINT *, 'Panic: Fell through direction finder!'
      stop
!      DO ; END DO ; STOP ! sit in loop for debugging purposes
    ENDIF

    else if (amr2d) then

       ! two-d case written by TJH on 27/1/05
       ! this code isn't ugly and it does work !!!!

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

       r1 = centre%x - subcellSize/2.d0
       r2 = centre%x + subcellSize/2.d0
       d = sqrt(currentpoint%x**2+currentpoint%y**2)
       xDir = VECTOR(currentpoint%x, currentpoint%y,0.d0)


       if ((direction%x**2+direction%y**2) /= 0.0) then
          if (modulus(xDir) /= 0.d0) then
             call normalize(xDir)       
             cosmu =((-1.d0)*xDir).dot.direction
             call solveQuadDble(1.d0, -2.d0*d*cosmu, d*d-r2*r2, x1, x2, ok)
             if (.not.ok) then
                write(*,*) "Error:: Quad solver failed in [amr_mod::getExitPoint] (1)"
                write(*,*) " Correctionn forced for now. If problem insistet, you need to fix it."
                write(*,*) "           d = ", d
                write(*,*) "      cosmud = ", cosmu
                write(*,*) "          r2 = ", r2
                write(*,*) "currentPoint = ", currentPoint
                x1 = 1.0d-5; x2=1.0d-5
!                stop
!                do;enddo
             endif
             distTor2 = max(x1,x2)      
         
         
             theta = asin(max(-1.d0,min(1.d0,r1 / d)))
             cosmu = xDir.dot.direction
             mu = acos(max(-1.d0,min(1.d0,cosmu)))
             distTor1 = 1.e30
             if (mu  < theta ) then
                call solveQuadDble(1.d0, -2.d0*d*cosmu, d*d-r1*r1, x1, x2, ok)
                if (.not.ok) then
                   write(*,*) "Quad solver failed in  [amr_mod::getExitPoint] (2)"
                   write(*,*) "           d = ", d
                   write(*,*) "      cosmud = ", cosmu
                   write(*,*) "          r2 = ", r2
                   write(*,*) "currentPoint = ", currentPoint
                   stop
!                   do;enddo
                endif
             endif
             distTor1 = max(x1,x2)
             distToXboundary = min(distTor1, distTor2)
          else
             distToXboundary = 1.e30
          endif
       else
          distToXboundary = 1.e30
       endif

       zDir = VECTOR(0.d0, 0.d0, 1.d0)
       compZ = zDir.dot.direction
       currentZ = currentpoint%z

       if (compZ /= 0.d0 ) then
          if (compZ > 0.d0) then
             distToZboundary = (centre%z + subcellsize/2.d0 - currentZ ) / compZ
          else
             distToZboundary = abs((centre%z - subcellsize/2.d0 - currentZ ) / compZ)
          endif
       else
          disttoZboundary = 1.e30
       endif

   
       tVal = min(distToZboundary, distToXboundary)
       if (tval <= 0.) tval = 1.0e-5
       if (tVal > 1.e29) then
          write(*,*) "Error :: tVal > 1.e29 [amr_mod:getExitPoint]."
          write(*,*) "tVal,compZ, distToZboundary,disttoxboundary = "
          write(*,*) tVal,compZ, distToZboundary,disttoxboundary
          write(*,*) "z = ", currentZ
          stop 
       elseif (tval <= 0.) then
          write(*,*) "Error :: tVal <= 0  [amr_mod:getExitPoint]."
          write(*,*) "tVal, compZ, distToZboundary,disttoxboundary = "
          write(*,*) tVal,compZ, distToZboundary,disttoxboundary
          write(*,*) "z = ", currentZ
          stop
       endif
   
       minWallDistance = tVal
   
       exitPoint = currentPoint + tval * direction
       locator = exitPoint + frac*halfSmallestSubcell * direction

!       write(*,*) tval, exitpoint

    else


       distToR1 = 1.d30
       distToR2 = 1.d30

       rVec = currentPoint
       call normalize(rVec)
       cosmu = ((-1.d0)*direction).dot.rVec
       d = modulus(currentPoint)

       ! distance to outer radius

       r2 = centre%x + subcellSize/2.d0
       call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
       distToR2 = max(x1,x2)
!             write(*,*) "r2",x1,x2,disttor2

       !   inner radius

       r1 = centre%x - subcellSize/2.d0
       theta = asin(max(-1.d0,min(1.d0,r1 / d)))
       cosmu =((-1.d0)*rVec).dot.direction
       mu = acos(max(-1.d0,min(1.d0,cosmu)))
       distTor1 = 1.e30
       if (mu  < theta ) then
          call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
          distTor1 = min(x1,x2)
       endif
!             write(*,*) "r1",x1,x2,disttor1,mu,theta

       tval = min(distTor1, distTor2)

       minWallDistance = tVal
   
       exitPoint = currentPoint + tval * direction
       locator = exitPoint + frac*halfSmallestSubcell * direction

    endif





  END SUBROUTINE getExitPoint 



  SUBROUTINE takeSample(point,length,direction,grid,thisOctal,subcell,nSamples,&
                        maxSamples,usePops,iLambda,error,lambda,kappaAbs,      &
                        kappaSca,velocity,velocityDeriv,chiLine,levelPop,rho,  &
                        temperature,Ne,inFlow,etaCont,etaLine) 
    use inputs_mod, only : molecularPhysics, atomicPhysics
    use jets_mod, only: JetsVelocity, get_jets_parameter, dV_dn_jets    
    use density_mod, only: density

    IMPLICIT NONE
    
    TYPE(vector), INTENT(IN)      :: point     ! place to make sample
    real(oct), INTENT(IN)   :: length    ! 
    TYPE(vector), INTENT(IN)      :: direction ! direction vector
    TYPE(gridtype), INTENT(IN)         :: grid      ! grid structure
    TYPE(octal), TARGET, INTENT(IN)    :: thisOctal ! grid structure
    INTEGER, INTENT(IN)                :: subcell   ! subcell containing 'point'
    INTEGER, INTENT(INOUT)             :: nSamples  ! number of samples so far
    INTEGER, INTENT(IN)                :: maxSamples! size of sample arrays 
    REAL, DIMENSION(:), INTENT(INOUT)  :: lambda    ! distances of samples from startPoint 
    LOGICAL, INTENT(IN)                :: usePops   ! whether to use level populations
    INTEGER, INTENT(IN)                :: iLambda   ! wavelength index
    INTEGER, INTENT(INOUT)             :: error     ! error code
    
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL          :: kappaAbs  ! continuous absorption opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL          :: kappaSca  ! scattering opacities 
    TYPE(vector), DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocity ! sampled velocities
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL          :: velocityDeriv   ! sampled velocity derivatives
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL          :: chiLine   ! line opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL          :: rho       ! density at sample points
    real(double),DIMENSION(:,:),INTENT(INOUT),OPTIONAL :: levelPop  ! level populations
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL   :: Ne        ! electron density
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL           :: temperature ! in [K]
    logical,DIMENSION(:),INTENT(INOUT),OPTIONAL        :: inFlow     ! indicates if the cell is in use.
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaCont ! contiuum emissivity
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaLine ! line emissivity


    TYPE(vector)                       :: directionReal ! direction vector (REAL values)
    TYPE(octal), POINTER               :: localPointer  ! pointer to the current octal
 
    TYPE(vector) :: starPosn
    TYPE(vector) :: pointVec
    REAL :: Vr, r, Rs
    
    starPosn = grid%starPos1
    pointVec = (point - starPosn)

    nSamples = nSamples + 1
    IF (nSamples > maxSamples) THEN
      PRINT *, "Erro:: nSamples > maxSamples in takeSample subroutine"
      PRINT *, "nSamples   = ", nSamples
      PRINT *, "maxSamples = ", maxSamples      
      write(*,*) "lambda(nSamples-10:nSamples) = ", lambda(nSamples-10:nSamples)
      STOP
    END IF

    
    lambda(nSamples) = real(length)
    directionReal = direction

    localPointer => thisOctal

    IF (usePops) THEN
      CALL amrGridValues(grid%octreeRoot,point,startOctal=localPointer,&
                         actualSubcell=subcell,                        &
                         iLambda=iLambda,                              &
                         direction=directionReal,                      &
                         velocity=velocity(nSamples),                  &
                         velocityDeriv=velocityDeriv(nSamples),        &
                         kappaAbs=kappaAbs(nSamples),                  &
                         kappaSca=kappaSca(nSamples),                  &
                         rho=rho(nSamples),                            &
                         N=levelPop(nSamples,:),                       &
                         grid=grid,                                    &
                         temperature=temperature(nSamples),            &
                         Ne=Ne(nSamples),                              &
                         inFlow=inFlow(nSamples)                       &
                         )
    ELSE
       if (molecularPhysics.or.atomicPhysics) then
          CALL amrGridValues(grid%octreeRoot,point,startOctal=localPointer,&
               actualSubcell=subcell,                        &
               iLambda=iLambda,                              &
               direction=directionReal,                      &
               velocity=velocity(nSamples),                  &
               velocityDeriv=velocityDeriv(nSamples),        &
               kappaAbs=kappaAbs(nSamples),                  &
               kappaSca=kappaSca(nSamples),                  &
               rho=rho(nSamples),                            &
               chiLine=chiLine(nSamples),                    &
               grid=grid,                                    &
               temperature=temperature(nSamples),            &
               Ne=Ne(nSamples),                              &
               inFlow=inFlow(nSamples),                      &
               etaCont=etaCont(nSamples),                    &
               etaLine=etaLine(nSamples)                     &                         
               )
       else
          CALL amrGridValues(grid%octreeRoot,point,startOctal=localPointer,&
               actualSubcell=subcell,                        &
               iLambda=iLambda,                              &
               direction=directionReal,                      &
               kappaAbs=kappaAbs(nSamples),                  &
               kappaSca=kappaSca(nSamples),                  &
               rho=rho(nSamples),                            &
               grid=grid,                                    &
               temperature=temperature(nSamples),            &
               Ne=Ne(nSamples),                              &
               inFlow=inFlow(nSamples),                      &
               etaCont=etaCont(nSamples),                    &
               etaLine=etaLine(nSamples)                     &                         
               )
       endif
    END IF

    
    ! special case if the geometry is "jets" (for accuracy)
    If (grid%geometry(1:4) == "jets" ) then
       ! using routines in jets_mod.f90
       rho(nSamples) = Density(pointVec, grid)                 ! in [g/cm^3]
       Vr = real(JetsVelocity(pointVec, grid)/(cSpeed/1.0e5))        ! in [c]
       r = real(modulus(pointVec))
       Rs = get_jets_parameter("Rmin")
       if (r < Rs) r =Rs
       velocity(nSamples) = REAL(Vr, kind=oct)*(pointVec/REAL(r,kind=oct))
       velocityDeriv(nSamples) = dV_dn_jets(VECTOR(pointVec%x, pointVec%y, pointVec%z), &
            VECTOR(direction%x, direction%y, direction%z))    ! [1/sec]
    end If

    ! use variables to silence compiler warnings
    error = error

  END SUBROUTINE takeSample

  
  SUBROUTINE amrGridValues(octalTree,point,startOctal,foundOctal,&
                           foundSubcell,actualSubcell,iLambda,lambda,direction,&
                           velocity,velocityDeriv,temperature,kappaAbs,&
                           kappaSca,rho,chiLine,etaLine,etaCont, &
                           probDistLine,probDistCont,N,Ne,nTot,inflow,grid, &
                           interp, departCoeff,kappaAbsArray,kappaScaArray, dusttypeFraction, rosselandKappa, kappap, &
                           atthistemperature)

    USE inputs_mod, only: hydrodynamics

    ! POINT, direction --> should be in unrotated coordinates for 2D case (not projected onto x-z plane!)
    !

    ! returns one or more physical variables at a given point in the grid.
    ! optional arguments should be specified for all of the variables that
    !   are wanted.
    
    ! if the foundOctal argument is supplied, it is made to point to 
    !   the octal containing 'point'.
    ! if the foundSubcell argument is supplied, it is made to point to 
    !   the subcell containing 'point'.
    ! if the startOctal argument is supplied, the function uses a 
    !   local search for the correct octal starting at that octal.
    ! if actualSubcell and startSubcell are both supplied, these 
    !   locations are assumed to be correct and no search is performed.

    IMPLICIT NONE

    TYPE(octal), POINTER, intent(in)              :: octalTree
    TYPE(vector), INTENT(IN)     :: point
    TYPE(vector)                 :: point2 ! may be projected point
    TYPE(octal), OPTIONAL, POINTER    :: startOctal
    TYPE(octal), OPTIONAL, POINTER    :: foundOctal
    INTEGER, INTENT(OUT), OPTIONAL    :: foundSubcell
    INTEGER, INTENT(IN), OPTIONAL     :: actualSubcell 
    INTEGER, INTENT(IN), OPTIONAL     :: iLambda       ! wavelength index
    TYPE(vector),INTENT(IN),OPTIONAL  :: direction     
    TYPE(gridtype),INTENT(IN),OPTIONAL:: grid          
    LOGICAL, INTENT(IN), OPTIONAL     :: interp        ! use interpolation
                                      !  ^^^^^^ ! not implemented yet???
    REAL,INTENT(IN),OPTIONAL          :: lambda
    TYPE(vector),INTENT(OUT),OPTIONAL :: velocity
    REAL,INTENT(OUT),OPTIONAL         :: velocityDeriv
    REAL,INTENT(OUT),OPTIONAL         :: temperature
    REAL(double),OPTIONAL,intent(out)         :: kappaAbs
    REAL(double),OPTIONAL,intent(out)         :: kappaSca
    REAL(double),OPTIONAL         :: kappaAbsArray(:)
    REAL(double),OPTIONAL         :: kappaScaArray(:)
    REAL(double),OPTIONAL        :: rosselandKappa
    REAL, OPTIONAL        :: kappap
    REAL,INTENT(IN), OPTIONAL         :: atThisTemperature
    REAL(double),INTENT(OUT),OPTIONAL         :: rho
    REAL,INTENT(OUT),OPTIONAL         :: chiLine
    REAL(double),INTENT(OUT),OPTIONAL         :: etaLine
    REAL(double),INTENT(OUT),OPTIONAL         :: etaCont
    real(double), dimension(:), intent(out), optional :: dusttypeFraction
    real(double),INTENT(OUT),OPTIONAL :: probDistLine
    real(double),INTENT(OUT),OPTIONAL :: probDistCont
    real(double),DIMENSION(:),INTENT(OUT),OPTIONAL :: N
    real(double),INTENT(OUT),OPTIONAL :: Ne
    real(double),INTENT(OUT),OPTIONAL :: nTot
    REAL,DIMENSION(:),INTENT(OUT),OPTIONAL     :: departCoeff
    LOGICAL,INTENT(OUT),optional      :: inFlow
    
    TYPE(octal), POINTER              :: resultOctal
    INTEGER                           :: subcell
    LOGICAL                           :: interpolate
    LOGICAL                           :: boundaryProblem = .false.

!    type(vector), save :: lastpoint, lastpoint2

! this for possibility of twoD AMR grid


    point2 = point
!    if(point .eq. lastpoint) then
!       point2 = lastpoint2
!    else
       if (octaltree%twoD) then
          if (.not.hydrodynamics .or. .not. cart2d) then
             point2 = projectToXZ(point)
          else
             point2 = point
          endif
       elseif(octaltree%oneD) then
          point2  = VECTOR(modulus(point),0.d0,0.d0)
       endif
!    endif
!    lastpoint = point
!    lastpoint2 = point2

    IF (PRESENT(interp)) THEN
      interpolate = interp
    ELSE
      interpolate = .FALSE. 
    END IF
    
    ! first we find the correct subcell (if we need to)
    
    IF (PRESENT(startOctal)) THEN
      IF (PRESENT(actualSubcell)) THEN
        subcell = actualSubcell
      ELSE 
    ! called with rotated position if 2d octal
        CALL findSubcellLocal(point2,startOctal,subcell, boundaryProblem)
        IF (PRESENT(foundOctal))   foundOctal   => startOctal
        IF (PRESENT(foundSubcell)) foundSubcell =  subcell
        if (boundaryProblem) then
          boundaryProblem = .false.
          CALL findSubcellTD(point2,octalTree,resultOctal,subcell)
          IF (PRESENT(foundOctal))   foundOctal   => resultOctal
          IF (PRESENT(foundSubcell)) foundSubcell =  subcell
        end if
      END IF
      resultOctal => startOctal
      
    ELSE
      CALL findSubcellTD(point2,octalTree,resultOctal,subcell)
      IF (PRESENT(foundOctal))   foundOctal   => resultOctal
      IF (PRESENT(foundSubcell)) foundSubcell =  subcell

    END IF
 
    IF (interpolate) THEN
       PRINT *, 'Interpolation not implemented!' 
       STOP
    ELSE


!      IF (PRESENT(velocity))         velocity = resultOctal%velocity(subcell)

       ! unrotated poistion should be passed here!

       IF (PRESENT(velocity)) then
          velocity = VECTOR(0.d0, 0.d0, 0.d0) 
          if (associated(resultOctal%cornerVelocity)) then
             velocity = amrGridVelocity(octalTree,point,startOctal=resultOctal,&
                  actualSubcell=subcell) 
          endif
       ENDIF



      IF (PRESENT(temperature))   temperature = resultOctal%temperature(subcell)
      IF (PRESENT(rho))                   rho = resultOctal%rho(subcell)
      IF (PRESENT(chiLine))           chiLine = real(resultOctal%chiLine(subcell))

!      IF (PRESENT(chiLine)) then
!         CALL interpAMR(resultOctal, point, chiLine)
!         write(*,*) point,resultOCTAL%chiLine(subcell), chiline
!      ENDIF

      IF (PRESENT(etaLine))           etaLine = resultOctal%etaLine(subcell)
      IF (PRESENT(etaCont))           etaCont = resultOctal%etaCont(subcell)
      IF (PRESENT(probDistLine)) probDistLine = resultOctal%probDistLine(subcell)
      IF (PRESENT(probDistCont)) probDistCont = resultOctal%probDistCont(subcell)
      IF (PRESENT(Ne).and.associated(resultOctal%ne))                     Ne = resultOctal%Ne(subcell)
      IF (PRESENT(N))                       N = resultOctal%N(subcell,:)
      IF (PRESENT(nTot))                 nTot = resultOctal%nTot(subcell)
      IF (PRESENT(departCoeff))   departCoeff = resultOctal%departCoeff(subcell,:)
      IF (PRESENT(inFlow))           inFlow = resultOctal%inFlow(subcell)     

      IF (PRESENT(dusttypeFraction))           dusttypeFraction = resultOctal%dusttypeFraction(subcell,:)     

      IF (PRESENT(kappaAbsArray)) THEN
         call returnKappa(grid, resultOctal, subcell, kappaAbsArray=kappaAbsArray)
      ENDIF
      IF (PRESENT(kappaScaArray)) THEN
         call returnKappa(grid, resultOctal, subcell, kappaScaArray=kappaScaArray)
      ENDIF

      IF (PRESENT(rosselandKappa)) THEN
         if (PRESENT(atthistemperature)) then
            call returnKappa(grid, resultOctal, subcell, rosselandKappa=rosselandKappa, &
                 atthisTemperature=atthisTemperature)
         else
            call returnKappa(grid, resultOctal, subcell, rosselandKappa=rosselandKappa)
         endif
      ENDIF

      IF (PRESENT(kappap)) THEN
         if (PRESENT(atthistemperature)) then
            call returnKappa(grid, resultOctal, subcell, kappap=kappap, &
                 atthistemperature=atthistemperature)
         else
            call returnKappa(grid, resultOctal, subcell, kappap=kappap)
         endif
      ENDIF

      IF (PRESENT(kappaAbs)) THEN 
        IF (PRESENT(iLambda)) THEN
!           if (.not.grid%oneKappa) then
!              kappaAbs = resultOctal%kappaAbs(subcell,iLambda)
!           else
              call returnKappa(grid, resultOctal, subcell, ilambda,&
              lambda=lambda, kappaAbs=kappaAbs)
!              if (resultOctal%gasOpacity) then
!                 call returnKappaValue(resultOctal%temperature(subcell), lambda, kappaAbs)
!              else
!                 IF (.NOT.PRESENT(lambda)) THEN
!                    kappaAbs = grid%oneKappaAbs(resultOctal%dustType(subcell),iLambda)*resultOctal%rho(subcell)
!                 ELSE
!                    kappaAbs = logint(lambda, grid%lamArray(ilambda), grid%lamArray(ilambda+1), &
!                         grid%oneKappaAbs(resultOctal%dustType(subcell),iLambda)*resultOctal%rho(subcell), &
!                         grid%oneKappaAbs(resultOctal%dustType(subcell),iLambda+1)*resultOctal%rho(subcell))
!                 ENDIF
!           endif
!        endif
        ELSE
          PRINT *, 'In amrGridValues, can''t evaluate ''kappaAbs'' without',&
                   ' ''iLambda''.'
          stop
        END IF
      END IF
      
      IF (PRESENT(kappaSca)) THEN 
        IF (PRESENT(iLambda)) THEN
!           if (.not.grid%oneKappa) then
!              kappaSca = resultOctal%kappaSca(subcell,iLambda)
!           else
              call returnKappa(grid, resultOctal, subcell, &
              ilambda,lambda=lambda, kappaSca=kappaSca)
!              if (resultOctal%gasOpacity) then
!                 call returnKappaValue(resultOctal%temperature(subcell), lambda, kappaSca)
!              else
!                 IF (.NOT.PRESENT(lambda)) THEN
!                    kappaSca = grid%oneKappaSca(resultOctal%dustType(subcell),iLambda)*resultOctal%rho(subcell)
!                 ELSE
!                    kappaSca = logint(lambda, grid%lamArray(ilambda), grid%lamArray(ilambda+1), &
!                         grid%oneKappaSca(resultOctal%dustType(subcell),iLambda)*resultOctal%rho(subcell), &
!                         grid%oneKappaSca(resultOctal%dustType(subcell),iLambda+1)*resultOctal%rho(subcell))
!           ENDIF
!           endif
!        endif
        ELSE
          PRINT *, 'In amrGridValues, can''t evaluate ''kappaSca'' without',&
                   ' ''iLambda''.'
          STOP
       END IF
     END IF

      IF (PRESENT(velocityDeriv)) THEN
        IF (PRESENT(direction) .AND. PRESENT(grid)) THEN
           ! unrotated poistion should be passed here!
          velocityDeriv = amrGridDirectionalDeriv(grid,point,direction,&
                                                  startOctal=resultOctal)
!          velocityDeriv = amrGridDirectionalDeriv(grid,point2,direction,&
!                                                  startOctal=resultOctal)
        ELSE
          PRINT *, 'In amrGridValues, can''t evaluate ''velocityDeriv'' without',&
                   ' ''direction'' and ''grid''.'
          STOP
        END IF
      END IF
    END IF
      

  END SUBROUTINE amrGridValues



  FUNCTION amrGridDirectionalDeriv(grid,position,direction,startOctal,&
                                   foundOctal,foundSubcell) 
    !
    ! POINT --> should be in unrotated coordinates for 2D case (not projected onto x-z plane!)
    !

    ! returns the directional derivative of velocity at a given point
    !   in the grid.
    ! this function can be called with just the first three arguments
    !   and it will start at the root of the octal tree to locate
    !   the correct octal.
    ! if the foundOctal argument is supplied, it is made to point to 
    !   the octal containing 'position'.
    ! if the foundSubcell argument is supplied, it is made to point to 
    !   the subcell containing 'position'.
    ! if the startOctal argument is supplied, the function uses a 
    !   local search for the correct octal starting at that octal.
    !   NOTE that startOctal may be *changed* by this function! 

    ! Note that the results of this function DO have a 1^10 factor in them. 
    !   (they are in cSpeed, but the distance is calculated 1e10cm)

    ! THIS SHOULD WORK ALSO IN TWO-D CASE as long as position and direction
    ! are both are not projected on z-x plane (un-rotated coodinates).

    IMPLICIT NONE

    TYPE(gridtype), INTENT(IN)     :: grid 
    REAL                           :: amrGridDirectionalDeriv
    TYPE(vector), INTENT(IN)  :: position
    TYPE(vector), INTENT(IN)       :: direction
    TYPE(octal), OPTIONAL, POINTER :: startOctal
    TYPE(octal), OPTIONAL, POINTER :: foundOctal
    TYPE(octal),  POINTER         :: thisOctal, currentOctal
    INTEGER, INTENT(OUT), OPTIONAL :: foundSubcell

    TYPE(vector)              :: octalDirection
    TYPE(octal), POINTER           :: firstOctal
    real(oct)           :: dr, dx, dphi
    real(oct)           :: r
    real(oct)           :: phi1, phi2
    TYPE(vector)              :: position1
    TYPE(vector)              :: position2
    INTEGER                        :: subcell


!    if (.not.grid%lineEmission) then
!       amrGridDirectionalDeriv = 1.e30
!       goto 666
!    endif
    
    octalDirection = direction

!! This call to amrGridValues
!    call amrGridValues(grid%octreeRoot,position,foundOctal=thisOctal,&
!                           foundSubcell=subcell)
!! is replaced by
    currentOctal => grid%octreeRoot
    thisOctal => grid%octreeRoot
    CALL findSubcellTD(position, currentOctal, thisOctal, subcell)
!! to avoid a recursive call to the non-recursive amrGridValues.
!! This subroutine is unused hence this change is untested. 
!! DMA 1/12/11

    ! dr is a small increment of distance
    dr = thisOctal%subcellSize * 1.d-1

    ! get a new position a little way back from current position
    position1 = position - (dr * octalDirection)
    
    ! this might be inside core or outside grid - in which case
    !   just use the current position as the first point

    r = modulus(position1)
    IF (.NOT. inOctal(grid%octreeRoot,position1) .OR. (r < grid%rCore)) THEN
      position1 = position
    END IF
      
    ! first line of sight velocity
    phi1 = direction .dot. (AMRgridVelocity(grid%octreeRoot,position1,&
                                            startOctal=startOctal,&
                                            foundOctal=firstOctal))
    
    ! now go forward a bit from current position
    position2 = position + (dr * octalDirection)

    ! check we're still inside grid
    r = modulus(position2)
    IF (.NOT. inOctal(grid%octreeRoot,position2) .OR. (r < grid%rCore)) THEN
      position2 = position
    END IF

    ! the second position l.o.s. velocity
    phi2 = direction .dot. (AMRgridVelocity(grid%octreeRoot,position2,&
                                            startOctal=firstOctal,&
                                            foundOctal=foundOctal,&
                                            foundSubcell=subcell))

    IF (PRESENT(foundSubcell)) foundSubcell = subcell                                        

    dx = modulus(position2 - position1)

    dphi = phi2 - phi1

    ! the line of sight velocity gradient

    IF (ABS(dx) >= 1.e-20) THEN
       amrGridDirectionalDeriv = real(abs(dphi / dx))
    ELSE
       amrGridDirectionalDeriv = 1.e-20
    ENDIF
    IF (.NOT. ABS(amrGridDirectionalDeriv) >= 1.e-10) &
       amrGridDirectionalDeriv = 1.e-10
!666 continue
  END FUNCTION amrGridDirectionalDeriv


  RECURSIVE SUBROUTINE locateContProbAMR(probability,thisOctal,subcell) 
    ! finds the subcell that contains a given value of 'probability'.
    ! each subcell of the tree's octals has a value for continuous emission 
    !   probability which is an upper bound of the subcell's value in the cumulative probability
    !   distribution for the 

    IMPLICIT NONE

    real(double), INTENT(IN) :: probability
    TYPE(octal), POINTER              :: thisOctal
    INTEGER, INTENT(OUT)              :: subcell

    INTEGER              :: i, j 
   
!    write(*,'(9f9.5)') probability,thisOctal%probDistCont(1:8)
    ! we need to treat the first subcell as a special case
    
    IF (probability < thisOctal%probDistCont(1)) THEN
    
      IF (thisOctal%hasChild(1)) THEN
      
        ! find the child
        DO j = 1, thisOctal%nChildren, 1
          IF (thisOctal%indexChild(j) == 1) THEN
            thisOctal => thisOctal%child(j)
            CALL locateContProbAMR(probability,thisOctal,subcell)
            RETURN
          END IF
        END DO
        
      ELSE 
        subcell = 1
        RETURN
        
      END IF
    END IF   
  
   
    DO i = 2, thisOctal%maxChildren, 1
    
      IF (probability > thisOctal%probDistCont(i-1) .AND. &
          probability < thisOctal%probDistCont(i)) THEN
      
        IF (thisOctal%hasChild(i)) THEN
          
          ! find the child
          DO j = 1, thisOctal%nChildren, 1
            IF (thisOctal%indexChild(j) == i) THEN
               thisOctal => thisOctal%child(j)
              CALL locateContProbAMR(probability,thisOctal,subcell)
              RETURN
            END IF
          END DO
          
        ELSE 
          subcell = i
          RETURN
          
        END IF
      END IF
    END DO
      
  END SUBROUTINE locateContProbAMR


  subroutine logSpaceGridCheck(position, dx, splitLog)
    use inputs_mod, only : amrgridsize, amrgridcentrex, npoints, Nmag
    
    integer :: i, nPerDivision!, depth
    real(double) :: DeltaX, thisdx, xmin, xmax, dx
    type(vector) :: position
    logical, intent(out) :: splitLog
    
    splitLog = .false.

    nPerDivision = int(real(npoints/Nmag))
!    print *, "nperdivision", nperdivision
    xmin = amrgridcentrex - (amrgridsize/2.d0)

    do i = Nmag, 0, -1
       thisdx = dble(((amrgridsize/10.d0**i))*(1.d0/nperdivision))
       DeltaX = dble(nPerDivision * thisdx)
       xmax = xmin + DeltaX
!       if((position%x-dx) > xmin .and. position%x < xmax) then
       if((position%x-dx/2.) < xmax) then
          !found where I reside...
 !         print *, "dx ", dx
  !        print *, "thisdx ", dx
          if(dx > thisdx) splitLog = .true.
       endif
       !update Xmin for next iteration
       xmin = xmin + DeltaX
    enddo
!    stop
  end subroutine logSpaceGridCheck


  FUNCTION decideSplit(thisOctal,subcell,amrLimitScalar,amrLimitScalar2,grid, wvars, splitInAzimuth,&
       romData) RESULT(split)
    ! returns true if the current voxel is to be subdivided. 
    ! decision is made by comparing 'amrLimitScalar' to some value
    !   derived from information in the current cell  

    use inputs_mod, only: height, betadisc, rheight, flaringpower, rinner, router, hydrodynamics
    use inputs_mod, only: drInner, drOuter, rStellar, cavangle, erInner, erOuter, rCore, &
         ttauriRouter, sphereRadius, amr1d, amr2d, amr3d
    use inputs_mod, only: warpFracHeight, warpRadius, warpSigma, warpAngle, hOverR
    use inputs_mod, only: solveVerticalHydro, hydroWarp, rsmooth
    use inputs_mod, only: rGap, gapWidth, rStar1, rStar2, mass1, mass2, binarysep, mindepthamr, &
         maxdepthamr, vturbmultiplier, rGapInner, rGapOuter
    use inputs_mod, only: planetgap, heightSplitFac, refineCentre, doVelocitySplit, ttauriRstar
    use inputs_mod, only: DW_rMin, DW_rMax,rSublimation, ttauridisc, ttauriwarp, ttauriRinner, amr2d
    use inputs_mod, only : phiRefine, dPhiRefine, minPhiResolution, SphOnePerCell
    use inputs_mod, only : dorefine, dounrefine, maxcellmass
    use inputs_mod, only : inputnsource, sourcepos, logspacegrid
    use inputs_mod, only : amrtolerance, refineonJeans, rhoThreshold, smallestCellSize, ttauriMagnetosphere, rCavity
    use inputs_mod, only : cavdens, limitscalar, addDisc
    use inputs_mod, only : discWind, planetDisc, sourceMass
    use luc_cir3d_class, only: get_dble_param, cir3d_data
    use cmfgen_class,    only: get_cmfgen_data_array, get_cmfgen_nd, get_cmfgen_Rmin
    use magnetic_mod, only : accretingAreaMahdavi
    use romanova_class, only:  romanova_density
    use mpi_global_mod, only:   myRankGlobal
    use magnetic_mod, only : inflowMahdavi, inflowBlandfordPayne
    use vh1_mod, only: get_density_vh1, vh1FileRequired
    use density_mod, only: density
    use angularImage_utils, only: galaxyInclination, galaxyPositionAngle, intPosX, intPosY, refineQ2Only
    use magnetic_mod, only : safierfits
! Currently commented out. Reinstate if required. 
!    use inputs_mod, only: ttauriwind, smoothinneredge, amrgridsize, amrgridcentrex, amrgridcentrey, amrgridcentrez

#ifdef USECFITSIO
    use gridFromFitsFile, only : checkFitsSplit
#endif
#ifdef SPH
    USE cluster_class, only:   find_n_particle_in_subcell
    use sph_data_class, only:  sphVelocityPresent, get_npart
#endif

    IMPLICIT NONE

    TYPE(octal), target, intent(inout) :: thisOctal
    type(OCTAL), pointer :: pOctal
    type(octal), pointer :: neighbourOctal
    integer :: neighbourSubcell
!    TYPE(octal), POINTER       :: thisOctal
    integer :: iSource
    INTEGER, INTENT(IN)        :: subcell
    LOGICAL, INTENT(INOUT) :: splitInAzimuth
    real(double), INTENT(IN) :: amrLimitScalar, amrLimitScalar2 ! used for split decision
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry
    !
    LOGICAL                    :: split        
    real(double) :: massratio, d1, d2, masstol
    real(oct)  :: cellSize
    TYPE(vector)     :: searchPoint, rVec
    TYPE(vector)     :: cellCentre
    REAL                  :: x, y, z
    REAL(double) :: hr, rd, fac, warpHeight, phi1, phi2, phi
    real(double) :: warpheight1, warpheight2
    real(double) :: warpradius1, warpradius2, height1, height2
    INTEGER               :: i
    real(double)      :: total_mass
    real(double), save :: rgrid(1000)
    real(double)      :: ave_density,  r, dr
    INTEGER               :: nr, nr1, nr2
    real(double)          :: minDensity, maxDensity, hillRadius
    INTEGER               :: nsample = 400
    INTEGER               :: npt_subcell
!    real(double) :: timenow
    real(double),save  :: R_tmp(204)  ! [10^10cm]
    real(double),allocatable, save  :: R_cmfgen(:)  ! [10^10cm]
    real(double),save  :: Rmin_cmfgen  ! [10^10cm]
    real(double) :: rho
    real(double) :: OstrikerRho(2), r0
    logical, save :: first_time=.true.
    logical :: close_to_star
    real(double)      :: thisScale
    real(double) :: h0
    logical :: inflow, insideStar,outSideStar
    logical,save  :: firstTime = .true.
    logical,save  :: firstTimeTTauri = .true.
    logical, intent(in) :: wvars
    real(double) :: lAccretion, thisHeightSplitFac
    real(double), save :: astar
    real(double) :: b, dphi, rhoc
    type(VECTOR) :: centre, dirVec(6), locator
    integer :: nDir
    real(double) :: maxGradient, grad, phiMax
    
    real(double) :: bigJ, rhoJeans, cs
#ifdef SPH
    real(double) :: massPerCell, n_bin_az
    real(double) :: vgradx, vgrady, vgradz, vgrad
    INTEGER      :: nparticle, limit
    real(double) :: chi, zeta0dash, psi, eta, zeta
    real(double) :: dummyDouble
    type(VECTOR) :: minV, maxV
    real(double) :: T, vturb
#endif
    real(double) :: dx, cornerDist(8), d, muval(8), r1, r2, v

    splitInAzimuth = .false.
    split = .false.


    if(wvars) then






       if(thisOctal%hasChild(subcell)) then
          split = .false.
          goto 101
       end if

       if (grid%splitOverMPI) then
          if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
             split = .false.
             goto 101
          endif
       end if


       if(logspacegrid) then
          rvec = subcellCentre(thisOctal, subcell)
          dx = grid%octreeroot%subcellSize*2.d0/(2.d0**thisOctal%nDepth)
          call logspacegridcheck(rvec, dx, split)
          goto 333
       endif


       if(minDepthAMR==maxDepthAMR) then
          split = .false.
          goto 101
       end if





       r = thisOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell
       centre = subcellCentre(thisOctal, subcell)
       if (thisOctal%threed) then
          nDir = 6
          dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
          dirVec(2) = VECTOR( 0.d0,+1.d0,  0.d0)
          dirVec(3) = VECTOR(+1.d0, 0.d0,  0.d0)
          dirVec(4) = VECTOR(-1.d0, 0.d0,  0.d0)
          dirVec(5) = VECTOR( 0.d0,-1.d0,  0.d0)
          dirVec(6) = VECTOR( 0.d0, 0.d0, -1.d0)
       else if (thisOctal%twod) then
          nDir = 4
          dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)
          dirVec(2) = VECTOR(-1.d0,0.d0, 0.d0)
          dirVec(3) = VECTOR( 0.d0, 0.d0,  1.d0)
          dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
       else
          nDir = 2
          dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)
          dirVec(2) = VECTOR(-1.d0, 0.d0, 0.d0)
       endif

       do i = 1, nDir
          maxGradient = 1.d-30
          locator = subcellCentre(thisOctal, subcell) + &
               (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell) * dirVec(i)
          if (inOctal(grid%octreeRoot, locator)) then
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
!             split = .false.

!Initially just checking rho             
             grad = abs((thisOctal%rho(subcell)-neighbourOctal%rho(neighbourSubcell)) / &
                  thisOctal%rho(subcell))
             maxGradient = max(grad, maxGradient)
             if (grad > amrtolerance .and. thisOctal%nDepth < maxDepthAMR) then
                split = .true.
                exit
             endif

             if((neighbourOctal%nDepth - thisOctal%ndepth) > 1 .and. thisOctal%nDepth < maxdepthAMR) then
                split = .true.
                exit
             end if

          end if
       end do


       if (inputnSource>0) then
          centre = subcellCentre(thisOctal, subcell)
          do iSource = 1, inputNsource
             r = modulus(centre - sourcePos(iSource))/smallestCellSize
             if ((r < 12.d0) .and. (thisOctal%nDepth < maxDepthAMR)) then
                split = .true.
                exit
             endif
             if ((r < 24.d0) .and. (thisOctal%nDepth < maxDepthAMR-1)) then
                split = .true.
                exit
             endif
          enddo
       endif

       if (inputnSource > 0) then
          do iSource = 1, inputnSource
             if (inSubcell(thisOctal,subcell, sourcePos(isource)) &
                  .and. (thisOctal%nDepth < maxDepthAMR)) then
                split = .true.
                exit
             endif
          enddo
       endif
             

       if (refineOnJeans) then
          massTol = (1.d0/8.d0)*rhoThreshold*1.d30*smallestCellSize**3
          if (((thisOctal%rho(subcell)*1.d30*thisOctal%subcellSize**3) > massTol) &
               .and.(thisOctal%nDepth < maxDepthAMR)) then 
!             write(*,*) "split ",thisOctal%rho(subcell)*1.d30*thisOctal%subcellSize**3/masstol
             
             split = .true.
          endif
       endif
       
333 continue       


    else

       !
       if(logspacegrid) then
          rvec = subcellCentre(thisOctal, subcell)
          dx = grid%octreeroot%subcellSize*2.d0/(2.d0**thisOctal%nDepth)
          call logspacegridcheck(rvec, dx, split)
          goto 222
       endif


       
       select case(grid%geometry)
          
       case("melvin")
          
          cellSize = thisOctal%subcellSize 
          cellCentre = subcellCentre(thisOctal,subCell)
          if (thisOctal%nDepth < 4) split = .true.
          if ((abs(cellCentre%z) < 2.e6) .and. &
               (cellCentre%x  < 1.e6).and.(cellSize > 1.e4)) split = .true.
          if ((density(cellCentre, grid) > 1.d29).and.(thisOctal%nDepth < 9)) split = .true.
          if ((cellCentre%x < 7.e6).and.(thisOctal%nDepth < 7)) split = .true.
                             
       case("radpress")
          rVec = subcellCentre(thisOctal,subcell)
          if (thisOctal%nDepth < minDepthAMR) split = .true.
          if (rCavity > 0.d0) then
             if (abs(modulus(rvec) - rCavity)/rCavity < 0.1d0) split = .true.
          endif

       case("hii_test")
          rInner = 0.2d0*pctocm/1.d10
          rVec = subcellCentre(thisOctal, subcell)
          if (modulus(rVec) < rInner) then
             split = .true.
          endif
          if ((modulus(rVec)-thisOctal%subcellSize/2.d0*sqrt(3.d0)) < rInner) then
             split = .true.
          endif
          if (thisOctal%nDepth < minDepthAMR) split = .true.
          
       case("shu")
          rVec = subcellCentre(thisOctal, subcell)
          if (thisOctal%nDepth < maxDepthAmr) then
             if (modulus(rVec) < 5.d0*thisOctal%subcellSize) split = .true.
          endif
          v = cellVolume(thisOctal, subcell) * 1.d30
          massTol = 0.1d0*v*rhoThreshold
          if (((thisOctal%rho(subcell)*v) > massTol) &
               .and.(thisOctal%nDepth < maxDepthAMR)) split = .true.

       case("whitney")          
          cellSize = thisOctal%subcellSize * 1.d10
          cellCentre = 1.d10 * subcellCentre(thisOctal,subCell)
          nr1 = 50
          nr2 = 10
          nr = nr1 + nr2
          
          do i = 1, nr1
             rgrid(i) = log10(0.5*drInner)+dble(i)*(log10(drOuter)-log10(0.5*drInner))/dble(nr1)
          end do
          do i = 1, nr2
             rgrid(nr1+i) = log10(drOuter)+dble(i)*(log10(erOuter)-log10(drInner))/dble(nr2)
          end do
          rgrid(1:nr) = 10.d0**rgrid(1:nr)
          r = modulus(cellcentre)
          if (thisOctal%nDepth < 5) split = .true.
          if ((r < rGrid(nr)).and.(r > rGrid(1))) then
             call locate(rGrid, nr, r, i)      
             if (cellsize > (rGrid(i+1)-rGrid(i))) split = .true.
          endif
          r = sqrt(cellcentre%x**2 + cellcentre%y**2)
          if ((r > drInner*0.9).and.(r < drOuter)) then
             hr = 0.01 * rStellar * (r / rStellar)**1.25
             if ((abs(cellcentre%z)/hr < 10.) .and. (cellsize/hr > 1.)) split = .true.
             if ((abs(cellcentre%z)/hr > 5.).and.(abs(cellcentre%z/cellsize) < 1.)) split = .true.
          endif
          dr = tan(cavAngle/2.) * abs(cellCentre%z)
          if ( ((abs(cellCentre%x) - cellsize/2.) < dr).and.(cellSize > dr/4.) .and.(abs(cellCentre%z)>erInner)) then
             split = .true.
          endif
          
          
       case("jets","spiralwind","romanova")
          nsample = 100
          ! the density is only sampled at the centre of the grid
          ! we will search in each subcell to see if any point exceeds the 
          ! threshold density
          
          ! get the size and centre of the current cell
          cellSize = thisOctal%subcellSize 
          cellCentre = subcellCentre(thisOctal,subCell)
          
          ! check the density of random points in the current cell - 
          !   if any of them exceed the critical density, set the flag
          !   indicating the cell is to be split, and return from the function
          
          ave_density = 0.0_db
          minDensity = 1.e30
          maxDensity = -1.e30
          DO i = 1, nsample
             searchPoint = cellCentre
             call randomNumberGenerator(getreal=x)
             call randomNumberGenerator(getreal=z)
             if (thisOctal%threed) then
                
                call randomNumberGenerator(getreal=y)
                
                searchPoint%x = searchPoint%x - (cellSize / 2.0_oc) + cellSize*REAL(x,KIND=oct)
                searchPoint%y = searchPoint%y - (cellSize / 2.0_oc) + cellSize*REAL(y,KIND=oct) 
                searchPoint%z = searchPoint%z - (cellSize / 2.0_oc) + cellSize*REAL(z,KIND=oct) 
             else
                searchPoint%x = searchPoint%x - (cellSize / 2.0_oc) + cellSize*REAL(x,KIND=oct) 
                searchPoint%y = 0.d0
                searchPoint%z = searchPoint%z - (cellSize / 2.0_oc) + cellSize*REAL(z,KIND=oct) 
             endif
             if (grid%geometry=="romanova") then
                rho = romanova_density(romData, searchPoint)
             else
                ! using a generic density function in density_mod.f90
                !           if (grid%geometry == "ttauri") then
                !              searchPoint = rotateY(searchPoint,  dble(dipoleOffset))
                !           endif
                
                rho = density(searchPoint,grid)
             end if
             ave_density  =  rho + ave_density
             minDensity = min(minDensity, rho)
             maxDensity = max(maxDensity, rho)
          END DO
          
          ave_density = ave_density / REAL(nsample,KIND=double)
          
          
          if (thisOctal%threed) then
             total_mass = maxDensity * (cellSize*1.e10_db)**3
          else
             total_mass = maxDensity * pi * ((searchPoint%x+cellsize/2.)**2-(searchPoint%x-cellsize/2.)**2)*cellsize*1.d30
          endif
          
          
          ! weigting toward the smaller radial positions
          if (grid%geometry=="romanova") then
             thisScale = grid%rstar1/modulus(cellCentre)
             total_mass = total_mass * thisScale**3
          elseif (grid%geometry=="ttauri") then
             thisScale = grid%rstar1/modulus(cellCentre)
             total_mass = total_mass * thisScale**4
          end if
          
          
          IF (total_mass > amrLimitScalar) then
             split = .TRUE.
          END IF
                    
          
          ! If 2D and uses large root cell special care must be taken
          if (grid%geometry == "ttauri" .and. .not. thisOctal%threeD) then         
             !         if (cellSize > 40.d0  .and.  &
             !              (subcell == 1 .or. subcell == 3) )  split=.true.
             close_to_star = .false.
             r = modulus(cellcentre)
             if (r < 300.d0) close_to_star =.true.
             if (close_to_star) then
                if (cellSize > 20.d0) split=.true.
             else
                ! Splits the two cells closer to the origin.
                if (r < cellSize*2.0d0)  split = .true.
             end if
             
             !
             ! We use a slightly higher resolution near the edge of the accretion flow
             !         if (ttau_fuzzy_edge) then
             if (in_fuzzy_edge(cellCentre) ) then
                if (total_mass*4.0d0 > amrLimitScalar) then
                   !               if (total_mass*20.0d0 > amrLimitScalar) then  ! for Halpha106
                   !               if (total_mass*40.0d0 > amrLimitScalar) then  ! for Halpha107
                   !           ^^^^^ Note the extra factor here
                   split = .TRUE.
                end if
             end if
             !         end if                                      
          end if          
          
          if (grid%geometry == "romanova") then         
             if (cellSize > 100.d0) split=.true.
          end if
          
          r = sqrt(cellcentre%x**2 + cellCentre%y**2)
          if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 15.)) then
             split = .true.
             splitInAzimuth = .true.
          endif
          
          r = sqrt(cellcentre%x**2 + cellCentre%y**2)
          phi = atan2(cellCentre%y,cellCentre%x)
          height = real(discHeightFunc( phi,hOverR) * r)
          if (((r-cellsize/2.d0) < (ttaurirOuter/1.d10)).and.( (r+cellsize/2.d0) > (ttaurirouter/1.d10))) then
             if (thisOctal%cylindrical.and.(thisOctal%dPhi*radtodeg > 2.)) then
                if (abs(cellCentre%z-height) < 2.d0*cellSize) then
                   split = .true.
                   splitInAzimuth = .true.
                endif
             endif
          endif
          
       case("slab")
          cellCentre = subcellCentre(thisOctal,subcell)
          d = thisOctal%subcellSize/2.d0
          if ((cellCentre%z < -2.d0*pctocm/1.d10).and.(thisOctal%nDepth < 7)) split = .true.
          if  (((cellCentre%z+d) > -2.d0*pctocm/1.d10).and.((cellCentre%z-d) < -2.d0*pctocm/1.d10)) then
             if (thisOctal%nDepth < maxDepthAMR) split = .true.
          endif
       case("lighthouse")
          cellCentre = subcellCentre(thisOctal,subcell)
          d = thisOctal%subcellSize/2.d0
          cornerDist(1) = modulus(VECTOR(cellCentre%x+d, cellCentre%y+d, cellCentre%z+d))
          cornerDist(2) = modulus(VECTOR(cellCentre%x-d, cellCentre%y+d, cellCentre%z+d))
          cornerDist(3) = modulus(VECTOR(cellCentre%x+d, cellCentre%y-d, cellCentre%z+d))
          cornerDist(4) = modulus(VECTOR(cellCentre%x-d, cellCentre%y-d, cellCentre%z+d))
          cornerDist(5) = modulus(VECTOR(cellCentre%x+d, cellCentre%y+d, cellCentre%z-d))
          cornerDist(6) = modulus(VECTOR(cellCentre%x-d, cellCentre%y+d, cellCentre%z-d))
          cornerDist(7) = modulus(VECTOR(cellCentre%x+d, cellCentre%y-d, cellCentre%z-d))
          cornerDist(8) = modulus(VECTOR(cellCentre%x-d, cellCentre%y-d, cellCentre%z-d))

          muval(1) = (cellCentre%z+d)/max(1.d-20,modulus(VECTOR(cellCentre%x+d, cellCentre%y+d, cellCentre%z+d)))
          muval(2) = (cellCentre%z+d)/max(1.d-20,modulus(VECTOR(cellCentre%x-d, cellCentre%y+d, cellCentre%z+d)))
          muval(3) = (cellCentre%z+d)/max(1.d-20,modulus(VECTOR(cellCentre%x+d, cellCentre%y-d, cellCentre%z+d)))
          muval(4) = (cellCentre%z+d)/max(1.d-20,modulus(VECTOR(cellCentre%x-d, cellCentre%y-d, cellCentre%z+d)))
          muval(5) = (cellCentre%z-d)/max(1.d-20,modulus(VECTOR(cellCentre%x+d, cellCentre%y+d, cellCentre%z-d)))
          muval(6) = (cellCentre%z-d)/max(1.d-20,modulus(VECTOR(cellCentre%x-d, cellCentre%y+d, cellCentre%z-d)))
          muval(7) = (cellCentre%z-d)/max(1.d-20,modulus(VECTOR(cellCentre%x+d, cellCentre%y-d, cellCentre%z-d)))
          muval(8) = (cellCentre%z-d)/max(1.d-20,modulus(VECTOR(cellCentre%x-d, cellCentre%y-d, cellCentre%z-d)))
          do i = 1, 8
             muVal(i) = acos(max(-1.d0,min(muval(i),1.d0)))
          enddo
          if (ANY(cornerDist(1:8) > 50.d0*autocm/1.d10).and. &
              ANY(cornerDist(1:8) < 50.d0*autocm/1.d10).and.(thisOctal%nDepth<11).and. &
              ANY(muVal(1:8) > cavAngle)) then
             split = .true.
          endif
          rVec = VECTOR(14960., 29920., -7480.)
          if (inSubcell(thisOctal, subcell, rVec).and.(thisOctal%nDepth<maxDepthAMR)) split = .true.
          if (((thisOctal%rho(subcell)*cellVolume(thisOctal,subcell)*1.d30) > limitScalar).and. &
               (thisOctal%nDepth < maxDepthAMR)) split = .true.

          if (thisOctal%cylindrical) then
             if (returndPhi(thisOctal) > minPhiResolution) then
                split = .true.
                splitinAzimuth = .true.
             endif
          endif
       case("ttauri")
          
          cellCentre = subcellCentre(thisOctal,subcell)
          cellSize = thisOctal%subcellSize
          r0 = modulus(cellCentre)
          
          if (ttauriMagnetosphere) then
             if (inflowMahdavi(cellcentre*1.d10).and.&
                  cellVolume(thisOctal,subcell)*1.d30*density(cellCentre,grid) > maxCellMass) &
                  split=.true.
          endif
          
          !     if (thisOctal%threed) then
          !        cellsize = MAX(cellsize, r0 * thisOctal%dphi)
          !     endif
          
          if (firstTimeTTauri) then
             astar = accretingAreaMahdavi()
             firstTimeTTauri = .false.
          endif
     
          lAccretion = ((astar * ((r0*1.d10)/ttauriRstar)**3) / (twoPi*r0*1.d10))/1.d10
          
          inflow=.false.
          insideStar = .false.
          outsideStar = .false.

          if (ttauriMagnetosphere) then

          do i = 1, 400
             rVec = randomPositionInCell(thisOctal,subcell)
             if (inFLowMahdavi(1.d10*rVec)) then
                inFlow = .true.
             endif
             r = modulus(rVec)
             if (r < (1.05*ttauriRstar/1.d10)) insideStar = .true.
          enddo
          r = sqrt(cellcentre%x**2 + cellCentre%y**2)
          
          fac = (r*1.d10-TTauriRInner)/(TTauriRouter-TTauriRinner)
                    
          fac = 5.d0
          
          if (inFlow) then
             !        write(*,*) "c/l ",cellsize/laccretion

             if (cellsize > lAccretion/fac) split = .true.
!             if (cellSize/(ttauriRstar/1.d10) > 0.05d0*(r0/(TTaurirStar/1.d10))**1.5d0) &
!                  split = .true.

!             if (cellSize/(ttauriRstar/1.d10) > 0.05d0*(r0/(TTaurirStar/1.d10))**1.5d0) &
!                  split = .true.


             !        if (cellSize > 0.0d0*(TTauririnner/1.d10)) split = .true.
!             if (insidestar.and.inflow.and.(thisOctal%dPhi*radtoDeg > 0.5d0)) then
!                split = .true.
!                splitinazimuth = .true.
!             endif
          endif
          
!          if (insidestar) then
!             if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 5.)) then
!                split = .true.
!                splitInAzimuth = .true.
!             endif
!          endif


          
          if (inflow) then
!             if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 7.5d0)) then
             phiMax = 7.5d0 !max(30.d0*(r0/(ttauriRouter/1.d10)),15.d0)
             if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > phiMax)) then
                split = .true.
                splitInAzimuth = .true.
             endif
          endif
       endif

          inSideStar = .true.
          do i = 1, 40
             rVec = randomPositionInCell(thisOctal,subcell)
             r = modulus(rVec)
             if (r > (ttauriRstar/1.d10)) insideStar = .false.
          enddo
          if (inSideStar) split = .false.
          
!          if (ttauriwind) then
!             inflow = .false.
!             do i = 1, 1000
!                rVec = randomPositionInCell(thisOctal,subcell)
!                if (inFLowBlandfordPayne(rVec)) then
!                   inFlow = .true.
!                   exit
!                endif
!             enddo
!        
!             if (inflow) then
!                if ((cellSize/(DW_rMax-DW_rMin)) > 0.05) split = .true.
!                !           if (abs(cellCentre%z) < 0.5*DW_rMax) then
!                !              if ((cellSize/(DW_rMax-DW_rMin)) > 0.02) split = .true.
!                !           endif
!             endif
!          endif
          
          if (ttauriwarp) then
             phi = atan2(cellCentre%y,cellCentre%x)
             r = sqrt(cellcentre%x**2 + cellCentre%y**2)
             height = real(discHeightFunc(phi,hOverR) * r)
             height1 = real(discHeightFunc(phi-thisOctal%dPhi/2.d0,hOverR) * r)
             height2 = real(discHeightFunc(phi+thisOctal%dPhi/2.d0,hOverR) * r)
             if (((r-cellsize/2.d0) < (ttaurirOuter/1.d10)).and.( (r+cellsize/2.d0) > (ttaurirouter/1.d10))) then
                if (thisOctal%cylindrical.and.(thisOctal%dPhi*radtodeg > 2.)) then
                   if (abs(cellCentre%z-height) < cellSize) then
                      split = .true.
                      splitInAzimuth = .true.
                   endif
                   if (abs(cellCentre%z-height1) < cellSize) then
                      split = .true.
                      splitInAzimuth = .true.
                   endif
                   if (abs(cellCentre%z-height2) < cellSize) then
                      split = .true.
                      splitInAzimuth = .true.
                   endif
                endif
             endif
          endif
          
          if (ttauridisc) then
             cellSize = thisOctal%subcellSize 
             cellCentre = subcellCentre(thisOctal,subCell)
             r = sqrt(cellcentre%x**2 + cellcentre%y**2)
             hr = height * (r / (100.d0*autocm/1.d10))**betadisc
             
             if ((r+cellsize/2.d0) > max(ttauririnner/1.d10,dble(rinner))) then
                if (r-cellsize/2.d0 > rSublimation) then
                   if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > 0.2)) split = .true.
                else
                   if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > 1.0)) split = .true.
                endif
                if ((abs(cellcentre%z)/hr > 2.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.
!                if (((r-cellsize/2.d0) < rSublimation).and. ((r+cellsize/2.d0) > rsublimation) .and. &
!                     (thisOctal%nDepth < maxdepthamr) .and. (abs(cellCentre%z/hr) < 5.d0) ) split=.true.
             endif
          endif
       case("lexington", "NLR")
          if (thisOctal%nDepth < mindepthamr) then
             split = .true.
          else
             split = .false.
          endif

       case("lexpdr")
          rVec = subcellCentre(thisoctal, subcell)

          if (thisOctal%nDepth < mindepthamr) then
             split = .true.
          else
             split = .false.
          endif
          
          if(rvec%x*1.d10/pctocm > 4.5 .and. rvec%x*1.d10/pctocm < 5.d0) then
             print *, "got a  max ", rvec%x*1.d10/pctocm
             if(thisOctal%ndepth < maxdepthamr) split = .true.
          endif

       case("toydisc")
          if (thisOctal%nDepth < mindepthamr) then
             split = .true.
          else
             split = .false.
          endif


       case("point")
          if (thisOctal%nDepth < mindepthamr) then
             split = .true.
          else
             split = .false.
          endif                    
          
       case("imgTest")
          if (thisOctal%nDepth < mindepthamr) then
             split = .true.
          else
             split = .false.
          endif          
#ifdef USECFITSIO
       case("fitsfile")
          split = checkFitsSplit(thisOctal,subcell)
#endif
       case("runaway")
          
          if (vh1FileRequired()) then 
             call get_density_vh1(thisOctal, subcell, ave_density, minDensity, maxDensity, npt_subcell)
          
             ! Split on mass per cell 
             total_mass = cellVolume(thisOctal, subcell)  * 1.d30 * ave_density
             if ( total_mass > amrlimitscalar .and. amrlimitscalar > 0.0 ) split = .true.
          
             ! Split on density contrast
             fac = ( maxDensity - minDensity ) / ( maxDensity + minDensity )
             if ( npt_subcell >= 2 .and. fac > amrlimitscalar2 .and. amrlimitscalar2 > 0.0 ) split = .true. 
          endif
          
       case("starburst")
          if (thisOctal%nDepth < mindepthamr) then
             split = .true.
          else
             split = .false.
          endif
          
       case("symbiotic")
          if (thisOctal%nDepth < 5) then
             split = .true.
          else
             split = .false.
          endif
          
       case("gammavel")
          
          if (thisOctal%nDepth < 6) then
             split = .true.
          else
             split = .false.
          endif
          massRatio = mass1/mass2
          
          d1 = binarySep * (1./(massRatio+1.))
          d2 = binarySep - d1
          
          rVec = subcellCentre(thisOctal,subcell)
          
          do i = 1, 100
             rgrid(i) = log10(0.9d0) + (log10(100.d0) - log10(0.9d0))*dble(i-1)/99.d0
          enddo
          rgrid(1:100) = 10.d0**rgrid(1:100)
          
          r = modulus(rVec - VECTOR(0.d0,0.d0,-d1))/rstar1
          if ((r > rgrid(1)).and.(r<rgrid(99))) then
             call locate(rgrid,100,r,i)
             if (thisOctal%subcellsize/rstar1 > (rgrid(i+1)-rgrid(i))) split = .true.
          endif
          
          r = modulus(rVec - VECTOR(0.d0,0.d0,d2))/rstar2
          if ((r > rgrid(1)).and.(r<rgrid(99))) then
             call locate(rgrid,100,r,i)
             if (thisOctal%subcellsize/rstar2 > (rgrid(i+1)-rgrid(i))) split = .true.
          endif
                            


       case ("testamr","proto")
          cellSize = thisOctal%subcellSize 
          cellCentre = subcellCentre(thisOctal,subCell)
          
          nr1 = 8
          nr2 = 100
          rgrid(1) = 0.8
          rgrid(2) = 0.9
          rgrid(3) = 0.999
          rGrid(4) = 1.
          rGrid(5) = 1.001
          rGrid(6) = 1.002
          rGrid(7) = 1.004
          rGrid(8) = 1.008
          rGrid(1:nr1) = log10(rGrid(1:nr1)*grid%rInner)
          nr = nr1 + nr2
          !      do i = 1, nr1
          !         rgrid(i) = log10(grid%rInner)+dble(i-1)*(log10(4.*grid%rInner)-log10(grid%rInner))/dble(nr1-1)
          !      end do
          do i = 1, nr2
             rgrid(nr1+i) = log10(1.01*grid%rInner)+dble(i)*(log10(grid%rOuter)-log10(1.01*grid%rInner))/dble(nr2)
          end do
          rgrid(1:nr) = 10.d0**rgrid(1:nr)
          r = modulus(cellcentre)
          if (thisOctal%nDepth < 4) split = .true.
          if ((r-cellsize/2.) < grid%rOuter) then
             call locate(rGrid, nr, r, i)      
             if (cellsize > (rGrid(i+1)-rGrid(i))) split = .true.
          endif
          
          !      if (thisOctal%nDepth > 3) split = .false.
       case ("wrshell")
          cellSize = thisOctal%subcellSize 
          cellCentre = subcellCentre(thisOctal,subCell)
          
          nr1 = 8
          nr2 = 100
          rgrid(1) = 0.8
          rgrid(2) = 0.9
          rgrid(3) = 0.999
          rGrid(4) = 1.
          rGrid(5) = 1.001
          rGrid(6) = 1.002
          rGrid(7) = 1.004
          rGrid(8) = 1.008
          rGrid(1:nr1) = log10(rGrid(1:nr1)*grid%rInner)
          nr = nr1 + nr2
          do i = 1, nr2
             rgrid(nr1+i) = log10(1.01*grid%rInner)+dble(i)*(log10(grid%rOuter)-log10(1.01*grid%rInner))/dble(nr2)
          end do
          rgrid(1:nr) = 10.d0**rgrid(1:nr)
          r = modulus(cellcentre)-cellsize/2.
          if (thisOctal%nDepth < 4) split = .true.
          call locate(rGrid, nr, r, i)      
          if (cellsize > (rGrid(i+1)-rGrid(i))) split = .true.
          if ( (r+cellsize) < rgrid(1)) split = .false.
      
          
       case("fluxTest")
          rVec = subcellCentre(thisOctal, subcell)
          if (thisOctal%nDepth < minDepthAMR) split = .true.
          if(rVec%x > 0.5 .and. thisOctal%nDepth < maxDepthAMR) split=.true.
          !  if(rVec%x > 0.6 .and. rVec%x < 0.7 .and. thisOctal%nDepth < maxDepthAMR) split=.true.
          
       case("hydro1d")
          if(ghostCell(grid, thisOCtal, subcell) .and. thisOctal%nDepth < maxDepthAMR) split = .true.
          
          if (thisOctal%nDepth < minDepthAMR) split = .true.
          
       case("diagSod")

          if (thisOctal%nDepth < minDepthAMR) split = .true.
          
          if(cornerCell(grid, thisOctal, subcell) .and. &
               thisOctal%nDepth < maxDepthAMR) split = .true.             
          
          if(.not. dorefine .and. .not. dounrefine) then
             if(thisOctal%twoD) then
                if(((rVec%x-0.5)**2 + rvec%z**2) < 0.05 .and. thisOctal%nDepth < maxDepthAMR) split=.true.
                
             else if (thisOctal%threeD) then
                if(((rVec%x-0.5)**2 + rvec%z**2) < 0.05 .and. thisOctal%nDepth < maxDepthAMR) split=.true.
             end if
          end if
          
!         if(dorefine .or. dounrefine) then
!             rVec = subcellCentre(thisOctal, subcell)
!                                       
!             if(thisOctal%twoD) then
!                if ( ((rVec%x+rvec%z) > 0.03).and. (((rVec%x+rvec%z) < 0.07)) .and. &
!                     (thisOctal%nDepth < maxDepthAMR)) split = .true.
!             else if(thisoctal%threeD) then 
!                if ( ((rVec%x+rvec%z) > 0.03).and. (((rVec%x+rvec%z) < 0.07)) .and. &
!                     (thisOctal%nDepth < maxDepthAMR)) split = .true.
!             else
!                print *, "1D diag sod doesn't work"
!                stop
!             end if
          
       case("bonnor", "empty", "unimed", "SB_WNHII", "SB_instblt", "SB_CD_1Da" & 
            ,"SB_CD_2Da" , "SB_CD_2Db", "SB_offCentre", "SB_isoshck", &
            "SB_coolshk", "SB_gasmix", "bubble", "SB_Dtype")

          if (thisOctal%nDepth < minDepthAMR) split = .true.

       case("isosphere")
          if(thisOctal%nDepth < mindepthamr) split = .true.


       case("rv1test", "rv2test", "rv3test", "rv4test")
          if(thisOctal%nDepth < mindepthamr) split = .true.
 
          rVec = subcellCentre(thisOctal, subcell)
          !refine LHS to maxdepth                   
          dx = grid%octreeroot%subcellSize*2.d0/(2.d0**maxDepthAMR)
          !found the xmin 
          if (grid%octreeroot%xmin == thisOctal%xmin .or. & 
               (grid%octreeroot%xmin+dx) == thisOctal%xmin .or. &
               (grid%octreeroot%xmin+(2.d0*dx)) == thisOctal%xmin) then
             if(thisOctal%ndepth < maxdepthamr) split = .true.
          endif

          if(rvec%x < grid%octreeroot%subcellSize/10.d0) then
             if(thisOctal%ndepth < maxdepthamr) split = .true.
          endif
!          if(grid%octreeroot%xmin + 10.d0*dx >= rVec%x) then
!             if(thisOctal%ndepth < maxdepthamr) split = .true.
!          endif
!          
!          
!          if(grid%octreeroot%xmin + 100000.d0*dx < rVec%x .and. grid%octreeroot%xmin + 10.d0*dx >= rVec%x) then
!             if(thisOctal%ndepth < maxdepthamr-1) split = .true.
!          endif
!          
!!          if(grid%octreeroot%xmin + 1000000.d0*dx < rVec%x .and. grid%octreeroot%xmin + 50.d0*dx >= rVec%x) then
 !            if(thisOctal%ndepth < maxdepthamr-2) split = .true.
 !         endif
! 
!          if(rvec%x > (grid%octreeroot%xmin + grid%octreeroot%subcellsize)) then
!             if(thisoctal%ndepth >= mindepthamr) split = .false.
!          endif
         
!          if(grid%octreeroot%xmin + 1000.d0*dx < rVec%x .and. grid%octreeroot%xmin + 10000.d0*dx >= rVec%x) then
!             if(thisOctal%ndepth < maxdepthamr-3) split = .true.
!          endif
!          
!          if(grid%octreeroot%xmin + 10000.d0*dx < rVec%x .and. grid%octreeroot%xmin + 100000.d0*dx >= rVec%x) then
!             if(thisOctal%ndepth < maxdepthamr-4) split = .true.
!          endif
!          
!          if(grid%octreeroot%xmin + 100000.d0*dx < rVec%x .and. grid%octreeroot%xmin + 1000000.d0*dx >= rVec%x) then
!             if(thisOctal%ndepth < maxdepthamr-5) split = .true.
!          endif
          
                         
       case("sphere")
!          if (thisOctal%subcellSize > (rGrid(i+1)-rGrid(i))) then
!             if (thisOctal%nDepth < minDepthAMR) split = .true.
!          endif
          bigJ = 0.25d0
          cs = sqrt(1.d0/(2.33d0*mHydrogen)*kerg*thisOctal%temperature(subcell))
          rhoJeans = max(1.d-30,bigJ**2 * pi * cs**2 / (bigG * (thisOctal%subcellSize*1.d10)**2)) ! krumholz eq 6
          massTol  = 0.011d0*rhoJeans * 1.d30* cellVolume(thisOctal,subcell)
          if (((thisOctal%rho(subcell)*cellVolume(thisOctal,subcell)*1.d30) > massTol) &
               .and.(thisOctal%nDepth < maxDepthAMR)) then
             split = .true.
!             write(*,*) "split on jeans",thisOctal%rho(subcell)*1.d30*thisOCtal%subcellSize**3 / masstol
          endif
       case("spiral")
          call splitSpiral(thisOctal, split, splitInAzimuth)
          if (thisOctal%nDepth < minDepthAMR) split = .true.
          if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 31.)) then
             split = .true.
             splitInAzimuth = .true.
          endif

       case("envelope")
          call calcEnvelope(thisOctal, subcell,checksplit=split)
!          r = modulus(thisOctal%centre)
!          if (firsttime) then
!             nr = 100
!             do i = 1, 100
!                rgrid(i) = log10(rinner) + (dble(i-1)/dble(nr-1))*log10(rOuter/rInner)
!             enddo
!             rgrid(1:nr) = 10.d0**rgrid(1:nr)
!             firsttime = .false.
!          endif
!          if ((r > rInner).and.(r < rOuter)) then
!             call locate(rGrid, 100, r ,i)
!             if (thisOctal%subcellSize > (rgrid(i+1)-rGrid(i))) split = .true.
!          endif
!          if (((r < rInner).and.(r+thisOctal%subcellsize > rInner))) split = .true.
!          if (((r > rInner).and.(r-thisOctal%subcellsize < rInner))) split = .true.

          if (thisOctal%nDepth < minDepthAMR) split = .true.

       case("unisphere")
          if (thisOctal%nDepth < minDepthAMR) split = .true.
          if (thisOctal%nDepth < halfRefined(minDepthAMR, maxDepthAMR)) split = .true.

       case("interptest")
          if (thisOctal%nDepth < minDepthAMR) split = .true.
                    
       case("turbulence")
          if (thisOctal%nDepth < 7) split = .true.

       case("gravtest")
          if (thisOctal%nDepth < minDepthAMR) split = .true.          
          rVec = subcellCentre(thisOctal, subcell)          
          if(modulus(rVec) > (sphereRadius-(sphereRadius*0.15d0)) .and. modulus(rVec) < (sphereRadius+(sphereRadius*0.15d0)) &
               .and. (thisOctal%nDepth < maxDepthAMR)) split = .true.
          
       case("brunt")
          if (thisOctal%nDepth < minDepthAMR) split = .true.
          
       case("radcloud")
          if (thisOctal%nDepth < minDepthAMR) split = .true.
          
       case("blobtest")
          if (thisOctal%nDepth < minDepthAMR) split = .true.
          if(cornerCell(grid, thisOctal, subcell) .and. (thisOctal%nDepth < maxDepthAMR)) split = .true.
          if(edgecell(grid, thisOctal, subcell) .and. (thisOctal%nDepth < maxDepthAMR)) split = .true.

       case("kelvin")
          if (thisOctal%nDepth < minDepthAMR) split = .true.          
          if ( (abs(thisOctal%zMax-0.25d0) < 1.d-5).and.(thisOctal%nDepth < maxDepthAMR)) split = .true.
          if ( (abs(thisOctal%zMin-0.25d0) < 1.d-5).and.(thisOctal%nDepth < maxDepthAMR)) split = .true.
          if ( (abs(thisOctal%zMax+0.25d0) < 1.d-5).and.(thisOctal%nDepth < maxDepthAMR)) split = .true.
          if ( (abs(thisOctal%zMin+0.25d0) < 1.d-5).and.(thisOctal%nDepth < maxDepthAMR)) split = .true.
          if(cornerCell(grid, thisOctal, subcell) .and. (thisOctal%nDepth < maxDepthAMR)) split = .true.
          
       case("rcTest")
          if (thisOctal%nDepth < minDepthAMR) split = .true.          

       case("rtaylor")
          if (thisOctal%nDepth < minDepthAMR) split = .true.
          if ( (abs(thisOctal%zMax) < 1.d-10).and.(thisOctal%nDepth < maxDepthAMR)) split = .true.
          if ( (abs(thisOctal%zMin) < 1.d-10).and.(thisOctal%nDepth < maxDepthAMR)) split = .true.
          
       case("gaussian")
          if (thisOctal%nDepth < minDepthAMR) split = .true.
          rVec = subcellCentre(thisOctal, subcell)
          if ((abs(rVec%x) < 0.25d0).and.(thisOctal%nDepth < maxDepthAMR)) split = .true.
          
       case("sedov")
          rInner = 0.02e0
          rVec = subcellCentre(thisOctal, subcell)
          if (sqrt((rVec%x-0.5d0)**2 + rVec%z**2) < rInner) then
             split = .true.
          endif
          if ((sqrt((rVec%x-0.5d0)**2 + rVec%z**2)-thisOctal%subcellSize/2.d0*sqrt(2.d0)) < rInner) then
             split = .true.
          endif
          
          ! Split is decided using mindepthAMR defined globally
          
       case("protobin")
          
          ! Split is decided using mindepthAMR defined globally
          
          massTol = (1.d0/8.d0)*rhoThreshold*1.d30*smallestCellSize**3
          if ((thisOctal%rho(subcell)*1.d30*thisOctal%subcellSize**3) > massTol) split = .true.

       case("parker")
!          rVec = subcellCentre(thisOctal, subcell)
          if(thisOctal%ndepth < mindepthamr) split = .true.
       case("fontdisc")
 !         rVec = subcellCentre(thisOctal, subcell)
          if(thisOctal%ndepth < mindepthamr) split = .true.

       case("simpledisc")
          if(thisOctal%ndepth < mindepthamr) split = .true.
  
          rVec = subcellCentre(thisOctal, subcell)
          
          if(modulus(rvec) < 5.*autocm/1.d10) then
             if(thisOctal%ndepth < maxdepthamr) split = .true.
          elseif(modulus(rvec) < 10.*autocm/1.d10) then
             if(thisOctal%ndepth < maxdepthamr-1) split = .true.
          elseif(modulus(rvec) < 30.*autocm/1.d10) then
             if(thisOctal%ndepth < maxdepthamr-2) split = .true.
          elseif(modulus(rvec) < 40.*autocm/1.d10) then
             if(thisOctal%ndepth < maxdepthamr-3) split = .true.
          elseif(modulus(rvec) < 80.*autocm/1.d10) then
             if(thisOctal%ndepth < maxdepthamr-4) split = .true.
          endif

       case("RHDDisc")
          rVec = subcellCentre(thisOctal, subcell)

          if(thisOctal%ndepth < mindepthamr) split = .true.


          if(rVec%x < 5.d-10*autocm) then
             if(thisOctal%ndepth < maxdepthamr) split = .true.
          elseif(rVec%x > 5.d-10*autocm  .and. rVec%x < 10.d-10*autocm)  then
             if(thisOctal%ndepth < maxdepthamr-1) split = .true.
          elseif(rVec%x > 10.d-10*autocm  .and. rVec%x < 15.d-10*autocm)  then
             if(thisOctal%ndepth < maxdepthamr-2) split = .true.
          endif

!
!          dx = amrgridsize/2.d0**(thisOctal%ndepth)
!!
!
!          if(rVec%x < (amrgridcentrex - (amrgridsize/2.d0)) + 2.d0*dx) split = .true.
!          if(rVec%x > (amrgridcentrex + (amrgridsize/2.d0)) - 2.d0*dx) split = .true.
!          if(rVec%z < (amrgridcentrez - (amrgridsize/2.d0)) + 2.d0*dx) split = .true.
!          if(rVec%z > (amrgridcentrez + (amrgridsize/2.d0)) - 2.d0*dx) split = .true.!

       case("benchmark")
          
          cellSize = thisOctal%subcellSize 
          cellCentre = subcellCentre(thisOctal,subCell)
          !      if (thisOctal%nDepth < 6) split = .true.
          r = sqrt(cellcentre%x**2 + cellcentre%y**2)
          rd = grid%rOuter / 2. 
          hr = height * (r/rd)**1.125
          
          !      if (.not.thisOctal%cylindrical) then
          if ((abs(cellcentre%z)/hr < 10.) .and. (cellsize/hr > 0.2)) split = .true.
          if ((abs(cellcentre%z)/hr > 5.).and.(abs(cellcentre%z/cellsize) < 0.2)) split = .true.
          !      else
          !         if ((abs(cellcentre%z)/hr < 10.) .and. (cellsize/hr > 0.2)) split = .true.
          !      endif
          
          if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 91.)) then
             splitInAzimuth = .true.
          endif
          if ((r+cellsize/2.d0) < grid%rinner) split = .false.
          if ((r-cellsize/2.d0) > rOuter) split = .false.
          
       case("molebench")
          cellCentre = subcellCentre(thisOctal, subcell)
          
          nr = 50
          if (firsttime) then
             open(31, file="model_1.dat", status="old", form="formatted")
             do i = nr,1,-1
                read(31,*) rgrid(i)
             enddo
             rgrid = rgrid * 1.e-10
             close(31)
             firsttime = .false.
          endif
          rd = modulus(cellCentre)
          call locate(rgrid, nr, rd, i)
          
          if (thisOctal%subcellSize > (rgrid(i+1)-rgrid(i))) split = .true.
          if (rd+0.5d0*thisOctal%subcellSize < rgrid(1)) split = .false.
          if (rd-0.5d0*thisOctal%subcellSize > rgrid(nr)) split = .false.
          
       case("filament")
          
          rhoc = 5.d-19
          r0 = sqrt(2.d0*(10.d0*kerg/(2.3d0*mHydrogen))/(pi*bigG*rhoc))/1.d10

          cellCentre = subcellCentre(thisOctal, subcell)
          rd = modulus(VECTOR(cellCentre%x, 0.d0, 0.d0))
          ! change the parameter
          rd = rd+thisOctal%subcellSize/2.d0
          OstrikerRho(1) = rhoc * (1.d0+(rd/r0)**2)**(-2.d0)
          rd = rd-thisOctal%subcellSize
          OstrikerRho(2) = rhoc * (1.d0+(rd/r0)**2)**(-2.d0)
          if(abs((OstrikerRho(1) - OstrikerRho(2))/rhoc) .ge. 0.02) split = .true.
          
       case("h2obench1")
          cellCentre = subcellCentre(thisOctal, subcell)
          ! mindepthamr replaces value of 8 in line below
          if (thisOctal%nDepth < mindepthamr) split = .true.
          nr = 200
          if (firsttime) then
             open(31, file="grid.dat", status="old", form="formatted")
             do i = 1,nr
                read(31,*) rgrid(i)
             enddo
             rgrid = rgrid * 3.08568025e8
             close(31)
             firsttime = .false.
          endif
          rd = modulus(cellCentre)
          call locate(rgrid, nr, rd, i)
          !      if (thisOctal%subcellSize > (rgrid(i+1)-rgrid(i))) split = .true.
          if (rd+0.5d0*thisOctal%subcellSize < rgrid(1)) split = .false.
          if (rd-0.5d0*thisOctal%subcellSize > rgrid(nr)) split = .false.
          
       case("h2obench2")
          cellCentre = subcellCentre(thisOctal, subcell)
          ! mindepthamr replaces value of 8 in line below
          if (thisOctal%nDepth < mindepthamr) split = .true.
          nr = 200
          if (firsttime) then
             open(31, file="grid.dat", status="old", form="formatted")
             do i = 1,nr
                read(31,*) rgrid(i)
             enddo
             rgrid = rgrid * 3.08568025e8
             close(31)
             firsttime = .false.
          endif
          rd = modulus(cellCentre)
          call locate(rgrid, nr, rd, i)
          !      if (thisOctal%subcellSize > (rgrid(i+1)-rgrid(i))) split = .true.
          if (rd+0.5d0*thisOctal%subcellSize < rgrid(1)) split = .false.
          if (rd-0.5d0*thisOctal%subcellSize > rgrid(nr)) split = .false.
          
       case("agbstar")
          cellCentre = subcellCentre(thisOctal, subcell)
          ! mindepthamr replaces 4
          if (thisOctal%nDepth < mindepthamr) split = .true.
          nr = 100
          if (firsttime) then
             open(31, file="mc_100.dat", status="old", form="formatted")
             do i = nr,1,-1
                read(31,*) rgrid(i)
             enddo
             close(31)
             rgrid = rgrid * 1e-10
             firsttime = .false.
          endif
          
          rd = modulus(cellCentre)
          call locate(rgrid, nr, rd, i)
          if (2. * thisOctal%subcellSize > (rgrid(i+1)-rgrid(i))) split = .true.
          if (rd+0.5d0*thisOctal%subcellSize < rgrid(1)) split = .false.
          if (rd-0.5d0*thisOctal%subcellSize > rgrid(nr)) split = .false.
          
       case("luc_cir3d") 
          if (first_time) then
             open(unit=77, file ="zeus_rgrid.dat", status="old")
             do i = 1, 204
                read(77,*) R_tmp(i)
             end do
             close(77)
             R_tmp(:) = R_tmp(:) * get_dble_param(cir3d_data,"Rs") ! [10^10cm]
             first_time=.false.         
          end if
          cellSize = thisOctal%subcellSize
          cellCentre = subcellCentre(thisOctal,subCell)
          r = modulus(cellcentre)
          nr=204
          call locate(R_tmp, nr, r, i)      
          if (i == 0) i = nr-1
          if (i == nr) i = nr-1
          if (cellsize*amrlimitscalar > (R_tmp(i+1)-R_tmp(i))) then
             split = .true.
          else
             split = .false.
          end if
          
          if (cellSize > 100.0d0)  split=.true.
          
       case("wind") 
          if (first_time) then
             nr = 50
             do i = 1, nr
                r_tmp(i) = log10(grid%rcore) + real(i-1)/real(nr-1) * (log10(2.d0*grid%rcore) - log10(grid%rcore))
             enddo
             nr = 50
             do i = 1, nr
                r_tmp(i+50) = log10(2.*grid%rcore) + &
                     real(i)/real(nr) * (log10(2.*grid%octreeRoot%subcellSize)-log10(2.*grid%rcore))
             enddo
             nr = 100
             r_tmp(1:nr) = 10.d0**r_tmp(1:nr)
             write(*,*) r_tmp(1:nr)/rcore
             first_time=.false.         
          end if
          cellSize = thisOctal%subcellSize
          cellCentre = subcellCentre(thisOctal,subCell)
          nr = 100
          r = modulus(cellcentre)+sqrt(2.)*cellSize/2.
          if (r > grid%rCore) then
             call locate(R_tmp, nr, r, i)      
             if (cellsize > (R_tmp(i+1)-R_tmp(i))) then
                split = .true.
             else
                split = .false.
             end if
          endif
          
       case("cmfgen") 
          nr = get_cmfgen_nd()
          if (first_time) then
             ! retriving the r grid in CMFGEN data.
             ALLOCATE(R_cmfgen(nr))
             call get_cmfgen_data_array("R", R_cmfgen) ! [10^10cm]
             Rmin_cmfgen = get_cmfgen_Rmin()
             first_time=.false.         
          end if
          cellSize = thisOctal%subcellSize
          cellCentre = subcellCentre(thisOctal,subCell)
          r = modulus(cellcentre)
          
          if (r > R_cmfgen(nr) .or. r < Rmin_cmfgen ) then
             split = .false.
          else
             call locate(R_cmfgen, nr, r, i)      
             if (i == 0) i = nr-1
             if (i == nr) i = nr-1
             
             dR = ABS(R_cmfgen(i+1)-R_cmfgen(i))
             thisScale = cellsize
             
             if ( thisScale > dR) then
                split = .true.
             else
                split = .false.
             end if
          end if
          if (cellSize > ABS(R_cmfgen(nr)-Rmin_cmfgen))  split=.true.
          if ((r < Rmin_cmfgen).and.((r+sqrt(2.d0)*cellsize/2.d0)>rMin_cmfgen).and.&
               (cellSize > (r_cmfgen(2)-r_cmfgen(1))) ) split = .true.
          
#ifdef SPH
       case ("sphfile","cluster","molcluster","theGalaxy", "dale")
          pOctal => thisOctal
          if (octalOnThread(pOctal, subcell, myrankGlobal).and.(get_npart() > 0)) then

             ! Switch off velocity splitting if SPH velocities are not available to avoid referencing unallocated pointer
             if (.not. sphVelocityPresent() ) doVelocitySplit =.false.
             
             if ( doVelocitySplit ) then 
                call find_n_particle_in_subcell(nparticle, ave_density, &
                     thisOctal, subcell, rho_min=minDensity, rho_max=maxDensity, n_pt=npt_subcell, v_min=minV, v_max=maxV)
             else
                call find_n_particle_in_subcell(nparticle, ave_density, &
                     thisOctal, subcell, rho_min=minDensity, rho_max=maxDensity, n_pt=npt_subcell)
             end if


             if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > dPhirefine ).and.(nParticle>1)) then
                split = .true.
                splitInAzimuth = .true.
             endif

             
             total_mass = cellVolume(thisOctal, subcell)  * 1.d30
             
             if ( thisOctal%cylindrical ) then
                n_bin_az = nint(twoPi / thisOctal%dPhi)
                massPerCell = ( (twoPi * thisOctal%r * 1.0e10) / n_bin_az ) * ( ave_density ** (1.0/3.0) ) * &
                     ( amrlimitscalar**(2.0/3.0) )
             else if ( amr2d ) then
                cellCentre = subcellCentre(thisOctal,subCell)
                massPerCell = ( (twoPi * cellCentre%x * 1.0e10) ) * ( ave_density ** (1.0/3.0) ) * &
                     ( amrlimitscalar**(2.0/3.0) )
             else
                massPerCell = amrlimitscalar
             end if
             
             ! placeholder for maximum expected smoothing length
             thisOctal%h(subcell) = ((maxdensity * total_mass) / mindensity)**(1.d0/3.d0)
             
             massPerCell = amrlimitscalar
             total_mass = ave_density * total_mass
             if (total_mass > massPerCell) then
                split = .true.
                splitInAzimuth = .true.
                mass_split = mass_split + 1
             endif
             
             ! Split in order to capture density gradients. 
             
             ! Use sign of amrlimitscalar2 to determine which condition to use.
             ! This allows testing of the grid generation method without recompiling the code
             
             if ( amrlimitscalar2 > 0.0 ) then 
                if ( ( (maxDensity-minDensity) / (maxDensity+minDensity) )  > amrlimitscalar2 ) then 
                   if(split) then
                      both_split = both_split + 1
                   else
                      density_split = density_split + 1
                      split = .true.
                   endif
                end if
             elseif (amrlimitscalar2 < 0.0 ) then
                if (  (maxDensity / minDensity) > abs(amrlimitscalar2) ) then 
                   if(split) then
                      both_split = both_split + 1
                   else
                      density_split = density_split + 1
                      split = .true.
                   endif
                end if
             endif
             
             !      if(maxdensity .gt. 1d-13 .and. nparticle .gt. 1) then
             !         split = .true.
             !         maxdensity_split = maxdensity_split + 1
             !      endif
             
! The velocity splitting condition is described in Rundle et al, 2010, MNRAS, 407, 986
! There are hardwired values here so caveat emptor if you use this option
             if(.not. split .and. (nparticle .ge. 2) .and. doVelocitySplit ) then
                if(ave_density .gt. 1d-13) then
                   T = 10.d0 * (ave_density * 1d13)**(0.4d0)
                   ! 5 is fudge factor to make sure condition isn't too strigent ! 28 is mass of CO
                   Vturb = 5.d0 * sqrt(2d-10 * kerg * T / (28.d0 * amu) + 0.3**2) / (cspeed * 1d-5)
                   Vturb = sqrt(5.938e-4 * T + 0.09) / (cspeed * 1d-5)
                else
                   Vturb = 1.03246E-06 ! above calculation with T = 10
                endif
                
                vgradx = maxV%x - minV%x
                vgrady = maxV%y - minV%y
                vgradz = maxV%z - minV%z
                
                vgrad = max(vgradx,max(vgrady,vgradz))
                
                if(vgrad .gt. vturbmultiplier * vturb) then
                   velocity_split = velocity_split + 1
                   split = .true.
                endif
             endif
             
             
             !Jeans mass condition
             !      if(.not. split .and. )then
             !         split = .true.
             !         mass_split2 = mass_split2 + 1
             !      endif
             
             ! Additional refinement at the grid centre used for SPH-Torus discs. 
             if ( refineCentre ) then 
                
                cellSize   = thisOctal%subcellSize
                cellCentre = subcellCentre(thisOctal,subCell)
                
                if ( thisOctal%cylindrical ) then 
                   
                   if ( abs(thisOctal%r) < 1.0e4 .and.  abs(cellCentre%z) < 1.0e4 ) then 
                      if ( cellSize > 1.0e3 ) split = .true.
                   end if
                
                   if ( abs(thisOctal%r) < 4.0e3 .and.  abs(cellCentre%z) < 4.0e3 ) then 
                      if ( cellSize > 4.0e2 ) split = .true.
                   end if
                   
                else
                   
                   if ( abs(cellCentre%x) < 10.0*grid%rCore .and.  abs(cellCentre%y) < 10.0*grid%rCore .and.  abs(cellCentre%z) &
                        < 40.0*grid%rCore ) then 
                      if ( cellSize > grid%rCore ) split = .true. 
                   end if
                   
                   if ( abs(cellCentre%x) < 2.5*grid%rCore .and.  abs(cellCentre%y) < 2.5*grid%rCore .and.  & 
                        abs(cellCentre%z) < 20.0*grid%rCore ) then 
                      if ( cellSize > 0.5*grid%rCore ) split = .true. 
                   end if
                   
                end if
                
             end if

! Split to one SPH particle per cell
             if (nparticle > 1 .and. SphOnePerCell) then
                split = .true. 
! Uncomment these lines if you want to know how many SPH particles are in each cell
!             elseif ( (.not.split).and.(.not.SphOnePerCell) .and. myRankIsZero ) then
!                write(113,*) nparticle
             end if

! Limit refinement to the second quadrant of a Galactic plane survey
             if ( refineQ2Only ) then                 

                ! Find this point on the unmodified grid
                cellCentre  = subcellCentre(thisOctal,subCell)
                cellCentre  = rotateY( cellCentre, -1.0*galaxyInclination*degToRad   )
                cellCentre  = rotateZ( cellCentre, -1.0*galaxyPositionAngle*degToRad )
                
                if ( cellCentre%x < (intPosX - 0.2e12) .or. cellCentre%y < (intPosY - 0.2e12) ) &
                   split = .false.
            
             end if
             
          ! The stellar disc code is retained in case this functionality needs to be used in future 
!!$
!!$      if (include_disc(stellar_cluster)) then
!!$      
!!$         ! If the subcell intersects with the stellar disk.. then do additional check.
!!$         ! This time, the cells will be split when the mass of the cell exceeds a 
!!$         ! critical mass specified by "amrlimitscalar.
!!$
!!$         if (stellar_disc_exists(sphData) .and.  &       
!!$              disc_intersects_subcell(stellar_cluster, sphData, thisOctal, subcell) ) then
!!$         rho_disc_ave = average_disc_density_fast(sphData, thisOctal, &
!!$              subcell, stellar_cluster, scale_length)
!!$
!!$!!         rho_disc_ave = average_disc_density(sphData, thisOctal, &
!!$!!              subcell, stellar_cluster, scale_length)
!!$!            rho_disc_ave = max_disc_density_from_sample(sphData, thisOctal, &
!!$!                 subcell, stellar_cluster, scale_length)
!!$
!!$            total_mass = total_mass + rho_disc_ave * (cellSize*1.e10_db)**3  !  [g]
!!$
!!$            if (cellSize > 1.0e4  .or.  &
!!$                 (cellsize > scale_length) ) then
!!$!                 (cellsize > scale_length .and. rho_disc_ave > 1.0e-20) ) then
!!$               split = .true.  ! units in 10^10cm
!!$            end if
!!$            
!!$         end if
!!$
!!$      end if
          endif

          if (addDisc) then
             rVec = subcellCentre(thisOctal, subcell)
             if (sqrt(rVec%x**2 + rVec%y**2) < 1.5*autocm/1.d10) then
                split = .false.
                splitinazimuth = .false.
             endif
          endif

       case ("wr104")
          ! Splits if the number of particle is more than a critical value (~3).
          limit = nint(amrLimitScalar)
          
          call find_n_particle_in_subcell(nparticle, dummyDouble, thisOctal, subcell)
          
          rVec = subcellCentre(thisOctal,subcell)
          r = sqrt(rVec%x**2 + rVec%y**2)
          
          if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 89.)) then
             splitInAzimuth = .true.
             limit = nint(amrLimitScalar)*1000
          endif
          
          if ((r < 2.e6).and.(thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 31.)) then
             limit = nint(amrLimitScalar)*100
             splitInAzimuth = .true.
          endif
          
          if ((r < 2.e5).and.(thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 15)) then
             splitInAzimuth = .true.
             limit = nint(amrLimitScalar)*1
             
          endif
          
          
          
          !
          if (nparticle > limit ) then 
             split = .TRUE.
          else
             split = .FALSE.
          end if
          cellCentre = subcellCentre(thisOctal,subCell)
                    
          !      write(*,*) nparticle,thisOctal%nDepth,subcell
#endif



       case("ggtau")
          ! used to be 5
          if (thisOctal%ndepth < mindepthamr) split = .true.
          
          cellCentre = subcellCentre(thisOctal,subCell)
          
          r = sqrt(cellcentre%x**2+cellcentre%y**2)
          r = r * 6.68458134e-06! (torus units to 100's of AUs)
          
          cellcentre%z = cellcentre%z * 6.68458134e-04
          cellSize = thisOctal%subcellSize * 6.68458134e-4
          
          H0 = 14.55
          hr = H0 * (r ** (1.25))
          
          if ((abs(cellcentre%z)/hr < 5.) .and. (cellsize/hr > 0.1)) then
             split = .true.
             scaleheighta_count = scaleheighta_count + 1
             goto 1001
          endif
          
          if ((abs(cellcentre%z)/hr > 5.) .and. (abs(cellcentre%z/cellsize) < 5.)) then
             split = .true.
             scaleheightb_count =  scaleheightb_count + 1
             goto 1001
          endif
      
          if ((r .gt. 1.8*0.99).and.(r .lt. 1.01*1.8)) then
             if ((abs(cellcentre%z)/hr < 1.)) then
                if (cellsize > 0.1*1.8) then
                   split = .true.
                   scaleheightc_count = scaleheightc_count + 1
                   goto 1001
                endif
             endif
          endif
      
1001      if ((r+(cellsize/(2.d0*100.))) < 1.8) split = .false.
          if ((r-(cellsize/(2.d0*100.))) > 8.) split = .false.
          
       case("shakara","aksco")
          ! used to be 5
          if (thisOctal%ndepth  < mindepthamr) split = .true.
          cellSize = thisOctal%subcellSize 
          cellCentre = subcellCentre(thisOctal,subCell)
          r = sqrt(cellcentre%x**2 + cellcentre%y**2)
          thisHeightSplitFac = heightSplitFac
          if (r < rSublimation) thisheightSplitFac = 1.
          hr = height * (r / (100.d0*autocm/1.d10))**betadisc
          
          !      if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > 1.)) split = .true.

!          write(*,*) hr, height, r, betadisc
          if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > thisheightSplitFac)) split = .true.
          
          if ((abs(cellcentre%z)/hr > 2.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.
          
!          if (.not.smoothinneredge) then
!             if (((r-cellsize/2.d0) < rSublimation).and. ((r+cellsize/2.d0) > rSublimation) .and. &
!                  (thisOctal%nDepth < maxdepthamr) .and. (abs(cellCentre%z/hr) < 4.d0) .and. &
!                  (.not.thisOctal%cylindrical)) split=.true.
!         endif
          
             if (((r-cellsize/2.d0) < rOuter).and. ((r+cellsize/2.d0) > rOuter)) then
                if ((thisOctal%subcellSize/rOuter > 0.01) .and. (abs(cellCentre%z/hr) < 7.d0)) then
                   if (.not.thisOctal%cylindrical) split = .true.
                endif
             endif
          
          if ((r+cellsize/2.d0) < rInner) split = .false.
          if ((r-cellsize/2.d0) > Router) split = .false.
          
          !      if ((r > grid%rinner).and.(r < 1.01d0*grid%rinner)) then
          !         if ((abs(cellcentre%z)/hr < 1.)) then
          !            if (cellsize > 5.d-1*grid%rinner) split = .true.
          !         endif
          !      endif
          
          if (thisOctal%cylindrical) then
             dphi = returndPhi(thisOctal)
             if ((r > 0.1*rGapInner).and.(r < 5.*rGapOuter)) then
                if (dphi > minPhiResolution) then
                   split = .true.
                   splitinAzimuth = .true.
                endif
             endif
          endif
          
          
          !
          
          if (thisOctal%cylindrical) then
             if ((r > rGapInner*0.95).and.(r < rGapOuter*1.05).and.(abs(cellCentre%z)< 5.d0*hr)) then
                phi = atan2(cellcentre%y, cellcentre%x)
                if (thisOctal%subcellSize > 0.2d0*(rGapOuter-rGapInner)) split = .true.
                if (phi < 0.d0) phi = phi + twopi
                dphi = returndPhi(thisOctal)
                phi1 = phi - dphi
                if (phi1 < 0.d0) phi1 = phi1 + twoPi
                phi2 = phi + dphi
                if (abs(cellCentre%z) < hr) then
                   if (dPhi*radtodeg > dPhiRefine) then
                      if (((phi1*radtodeg > 180.d0-phiRefine/2.d0).and. &
                           (phi1*radtodeg < 180.d0+phiRefine/2.d0)).or. &
                           ((phi2*radtodeg > 180.d0-phiRefine/2.d0).and. &
                           (phi2*radtodeg < 180.d0+phiRefine/2.d0)).or. &
                           ((phi1*radtodeg < 180.d0-phiRefine/2.d0).and. &
                           (phi2*radtodeg > 180.d0+phiRefine/2.d0))) then
                         splitInAzimuth = .true.
                         split = .true.
                      endif
                   endif
                endif
             endif
          endif
          !
          if (((r-cellsize/2.d0) > rOuter*1.1d0).and.(thisOctal%nDepth > 4)) then
             split = .false.
             splitInAzimuth = .false.
          endif
          
          if (hydrodynamics) then
             split = .false.
             if (thisOctal%nDepth < minDepthAMR) split = .true.
          endif
          
          !      if ((r > grid%rinner).and.(r < 1.001d0*grid%rinner)) then
          !         if ((abs(cellcentre%z)/hr < 2.)) then
          !            split = .true.
!         endif
          !      endif

          cellsize = thisOctal%subcellSize
          cellCentre = subcellCentre(thisOctal, subcell)
          if (cavdens > 1.d-30) then
             dr = tan(cavAngle/2.) * abs(cellCentre%z)
             if ( ((abs(cellCentre%x) - cellsize/2.d0) < dr).and.(cellSize > dr/8.d0) .and.(abs(cellCentre%z)>erInner/1.d10)) then
                split = .true.
             endif
          endif


          if (discWind) then
            cellCentre = subcellCentre(thisOctal, subcell)
            cellSize = thisOctal%subcellSize
            zeta0dash = tan(60.d0*degtorad)

            cellCentre = subcellCentre(thisOctal, subcell) + VECTOR(-cellsize/2.d0, 0.d0, -cellSize/2.d0)
            chi = abs(cellcentre%z)/DW_rMin
            call safierFits("C", chi, zeta0dash, 0.d0, zeta, psi, eta)
            r1 = DW_rMin*zeta
            chi = abs(cellcentre%z)/DW_rMax
            call safierFits("C", chi, zeta0dash, 0.d0, zeta, psi, eta)
            r2 = DW_rMax*zeta
            if ((cellCentre%x > r1).and.(cellCentre%x  < r2)) then
               if (thisOctal%subcellSize > (DW_rMax-DW_rMin)/5.d0) split = .true.
            endif

            cellCentre = subcellCentre(thisOctal, subcell) + VECTOR(+cellsize/2.d0, 0.d0, -cellSize/2.d0)
            chi = abs(cellcentre%z)/DW_rMin
            call safierFits("C", chi, zeta0dash, 0.d0, zeta, psi, eta)
            r1 = DW_rMin*zeta
            chi = abs(cellcentre%z)/DW_rMax
            call safierFits("C", chi, zeta0dash, 0.d0, zeta, psi, eta)
            r2 = DW_rMax*zeta
            if ((cellCentre%x > r1).and.(cellCentre%x  < r2)) then
               if (thisOctal%subcellSize > (DW_rMax-DW_rMin)/5.d0) split = .true.
            endif

            cellCentre = subcellCentre(thisOctal, subcell) + VECTOR(-cellsize/2.d0, 0.d0, +cellSize/2.d0)
            chi = abs(cellcentre%z)/DW_rMin
            call safierFits("C", chi, zeta0dash, 0.d0, zeta, psi, eta)
            r1 = DW_rMin*zeta
            chi = abs(cellcentre%z)/DW_rMax
            call safierFits("C", chi, zeta0dash, 0.d0, zeta, psi, eta)
            r2 = DW_rMax*zeta
            if ((cellCentre%x > r1).and.(cellCentre%x  < r2)) then
               if (thisOctal%subcellSize > (DW_rMax-DW_rMin)/5.d0) split = .true.
            endif

            cellCentre = subcellCentre(thisOctal, subcell) + VECTOR(+cellsize/2.d0, 0.d0, +cellSize/2.d0)
            chi = abs(cellcentre%z)/DW_rMin
            call safierFits("C", chi, zeta0dash, 0.d0, zeta, psi, eta)
            r1 = DW_rMin*zeta
            chi = abs(cellcentre%z)/DW_rMax
            call safierFits("C", chi, zeta0dash, 0.d0, zeta, psi, eta)
            r2 = DW_rMax*zeta
            if ((cellCentre%x > r1).and.(cellCentre%x  < r2)) then
               if (thisOctal%subcellSize > (DW_rMax-DW_rMin)/5.d0) split = .true.
            endif
          endif
          if (thisOctal%ndepth  < mindepthamr) split = .true.



          if (planetDisc) then
             if (thisOctal%cylindrical) then


                hillRadius = modulus(sourcePos(2)-sourcePos(1))*sqrt(sourceMass(2)/(3.d0*sourceMass(1)))
                rVec = sourcePos(2) - cellCentre
                r = modulus(rVec)

                rVec = sourcePos(2) + VECTOR(0.d0, hillRadius*0.01d0, 0.d0)
                if (inSubcell(thisOctal, subcell, rVec)) then
                  dphi = returndphi(thisOctal)
                  if (dphi * sqrt(cellCentre%x**2 + cellCentre%y**2) > 0.05d0*hillRadius) then
                     split = .true.
                     splitInAzimuth = .true.
                  endif
                endif
                rVec = sourcePos(2) + VECTOR(0.d0, -hillRadius*0.01d0, 0.d0)
                if (inSubcell(thisOctal, subcell, rVec)) then
                  dphi = returndphi(thisOctal)
                  if (dphi * sqrt(cellCentre%x**2 + cellCentre%y**2) > 0.05d0*hillRadius) then
                     split = .true.
                     splitInAzimuth = .true.
                  endif
                endif


                if (r < hillRadius) then
                   fac = radtodeg*(0.02d0 * hillradius / sqrt(cellCentre%x**2 + cellCentre%y**2))
                   if (dPhi*radtodeg > fac) then
                      splitInAzimuth = .true.
                      split = .true.
                   endif
                  if (thisOctal%subcellsize > 0.02d0*hillRadius) then
                     split = .true.
                  endif
                endif
             endif
          endif


          if (thisOctal%cylindrical) then
             if (split.and.(.not.splitInAzimuth)) then
                phi = atan2(cellcentre%y, cellcentre%x)
                if (phi < 0.d0) phi = phi + twopi
                dphi = returndPhi(thisOctal)
                phi1 = phi - dphi
                if (phi1 < 0.d0) phi1 = phi1 + twoPi
                phi2 = phi + dphi
                if (dPhi*radtodeg > 10.) then
                   if (((phi1*radtodeg > 180.d0-phiRefine/2.d0).and. &
                        (phi1*radtodeg < 180.d0+phiRefine/2.d0)).or. &
                        ((phi2*radtodeg > 180.d0-phiRefine/2.d0).and. &
                        (phi2*radtodeg < 180.d0+phiRefine/2.d0)).or. &
                        ((phi1*radtodeg < 180.d0-phiRefine/2.d0).and. &
                        (phi2*radtodeg > 180.d0+phiRefine/2.d0))) then
                      splitInAzimuth = .true.
                      split = .true.
                   endif
                endif
             endif
          endif
          

       case("circumbin")
          
          if (thisOctal%ndepth  < 5) split = .true.
          cellSize = thisOctal%subcellSize 
          cellCentre = subcellCentre(thisOctal,subCell)
          r = sqrt(cellcentre%x**2 + cellcentre%y**2)
          hr = height * (r / (100.d0*autocm/1.d10))**betadisc
          
          !      if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > 1.)) split = .true.
          if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > 0.2)) split = .true.
          
          if ((abs(cellcentre%z)/hr > 2.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.
          
          if ((r+cellsize/2.d0) < grid%rinner*1.) split = .false.
          if ((r-cellsize/2.d0) > grid%router*1.) split = .false.

          if (((r-cellsize/2.d0) < grid%rinner).and. ((r+cellsize/2.d0) > grid%rInner) .and. &
               (thisOctal%nDepth < maxDepthAmr) .and. (abs(cellCentre%z/hr) < 3.d0) ) split=.true.
          
          !      splitInAzimuth = .false.
          !      if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 181.)) then
          !         splitInAzimuth = .true.
          !         split = .true.
          !      endif
          
          !      if (abs(cellCentre%z) < rinner/2.) then
          !         if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 11.).and.(r < rInner)) then
          !            splitInAzimuth = .true.
          !            split = .true.
          !         endif
          !      endif
          !      if (abs(cellCentre%z) < rinner/5.) then
          !         if ((r < rinner).and.(thisOctal%subcellSize > (0.05*rinner))) split = .true.
          !      endif
          
          if ((r > rOuter*1.1d0).and.(thisOctal%nDepth > 4)) then
             split = .false.
             splitInAzimuth = .false.
          endif
          
       case("turbbox")
          if (thisOctal%nDepth < maxDepthAMR) split = .true.
          if (thisOctal%nDepth < minDepthAMR) split = .true.
          
       case("iras04158")
          
          cellSize = thisOctal%subcellSize 
          cellCentre = subcellCentre(thisOctal,subCell)
          r = sqrt(cellcentre%x**2 + cellcentre%y**2)
          hr = height * (r / (50.d0*autocm*1.d-10))**betadisc
          
          !! This splitting law replaces the commented out code below. It does the same stuff but further out
          !! and more gradually. I hope the splitting in r has become redundant because it should be linked to
          !! the increased resolution in z.
          
          if(abs(cellcentre%z) .lt. 6.d0 * hr) then ! 6 is a reasonable number for this law
             !         if(cellsize .gt. 2.d0**(abs(cellcentre%z)/hr - 3.d0) * hr) split = .true.
             if(8.d0 * cellsize .gt. 2.d0**(abs(cellcentre%z)/hr) * hr) split = .true.
          endif
          
          if((abs(cellcentre%z) .gt. 3.d0 * hr).and.(abs(cellcentre%z) < 4.d0 * cellsize)) then
             split = .true.
          endif
          
          !      if ((abs(cellcentre%z) < 1.d0 * hr) .and. (cellsize > 0.25d0 * hr)) then
          !         split = .true.
          !      elseif ((abs(cellcentre%z) < 2.d0 * hr) .and. (cellsize > 0.5d0 * hr)) then
          !         split = .true.
          !      elseif((r > grid%rInner).and.(r < grid%rInner * 1.02)) then
          !         if ((abs(cellcentre%z) < 3. * hr)) then        
          !            if (cellsize > 0.1d0 * grid%rInner) split = .true.
          !         endif
          !      endif
          
          !      if((r > grid%rInner).and.(r < grid%rInner * 1.1)) then
          !         if ((abs(cellcentre%z) < 1.d0 * hr)) then        
          !            if (cellsize > grid%rInner) split = .true.
          !         endif
          !      endif
          
          if ((r-cellsize*0.5d0) > grid%router*1.005) then
             split = .false.
             goto 555
          elseif ((r+cellsize*0.5d0) < grid%rinner*0.995) then
             split = .false.
             goto 555
          endif
          
555       continue
          
       case("oldshakara")
          
          !      if (thisOctal%ndepth  < 6) split = .true.
          cellSize = thisOctal%subcellSize 
          cellCentre = subcellCentre(thisOctal,subCell)
          r = sqrt(cellcentre%x**2 + cellcentre%y**2)
          hr = height * (r / (100.d0*autocm/1.d10))**betadisc
          
          if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > 0.2)) split = .true.
          
          !      if ((abs(cellcentre%z)/hr < 5.) .and. (cellsize/hr > 0.2)) split = .true.
          
          if ((abs(cellcentre%z)/hr > 5.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.
          
          !      if ((r > grid%rInner).and.(r < grid%rInner * 1.01)) then
          !         if ((abs(cellcentre%z)/hr < 5.)) then
          !            if (cellsize > 1.e-3 * grid%rInner) split = .true.
          !         endif
          !      endif
          
          

          if ((r+cellsize/2.d0) < grid%rinner*0.9) split = .false.
          
          !      if ((r > grid%rinner).and.(r < 1.01d0*grid%rinner)) then
          !         if ((abs(cellcentre%z)/hr < 2.)) then
          !            if (cellsize > 1.d-3*grid%rinner) split = .true.
          !         endif
          !      endif
          
          !      if ((r > grid%rinner).and.(r < 1.001d0*grid%rinner)) then
          !         if ((abs(cellcentre%z)/hr < 2.)) then
          !            split = .true.
          !         endif
          !      endif
          

       case("ppdisk")
          
          cellSize = thisOctal%subcellSize
          cellCentre = subcellCentre(thisOctal,subCell)
          r = sqrt(cellcentre%x**2 + cellcentre%y**2)
          ! The 100 is because Tim's code uses the scale height at 100AU.
          !      hr = height * (r / (100.d0*autocm/1.d10))**1.25
          hr = (auTocm * height * rHeight * (r / (rHeight*autocm/1.e10))**flaringPower)/1.e10
          if ((abs(cellcentre%z)/hr < 5.) .and. (cellsize/hr > 0.5)) split = .true.
          if ((abs(cellcentre%z)/hr > 5.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.
          if ((r+cellsize/2.d0) < rsmooth*autocm/1.e10) split = .false.
          
       case("planetgap")
          
          cellSize = thisOctal%subcellSize
          cellCentre = subcellCentre(thisOctal,subCell)
          r = sqrt(cellcentre%x**2 + cellcentre%y**2)
          hr = height * rCore * (r/rCore)**betaDisc
          if ((abs(cellcentre%z)/hr > 5.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.
          if ((abs(cellcentre%z)/hr < 10.) .and. (cellsize/hr > 0.5)) split = .true.
          ! The next few lines are used in the Varniere models instead of the line
          ! above. The final line tries to prevent splitting in the denser
          ! midplane of the disc.
          !      if (r > 6000.) then
          !        if ((abs(cellcentre%z)/hr < 10.) .and. (cellsize/hr > 0.5)) split = .true.
          !      else
          !        if ((abs(cellcentre%z)/hr < 12.) .and. (cellsize/hr > 0.5)) split = .true.
          !      end if
          !      if (abs(cellcentre%z)/hr < 0.9) split = .false.
          ! Added in to get an extra level of refinement at the outer edge of the
          ! gap. Cells at radii between that of the gap and (gap+gapWidth) must
          ! have a resolution of a quarter of a scale height instead of a half...
          if (planetgap) then
             if ((r > rGap*autocm/1e10) .and. (r < (rGap+1.5*gapWidth)*autocm/1e10) .and.  &
                  (abs(cellcentre%z)/hr < 10.) .and. (cellsize/hr > 0.25)) split = .true.
             ! The alternative version below was used for the models where we want a
             ! high resolution gap (hi-res-gap models).
             !        if ((r > (rGap-gapWidth)*autocm/1e10) .and. (r < (rGap+1.5*gapWidth)*autocm/1e10) .and.  &
             !            (abs(cellcentre%z)/hr < 10.) .and. (cellsize/hr > 0.1)) split = .true.
          end if
          if ((r+cellsize/2.d0) < 0.9*grid%rinner) split = .false.
          
       case("clumpyagb")
          if (thisOctal%nDepth < minDepthAMR) split = .true.
          

       case("clumpydisc")
          
          cellSize = thisOctal%subcellSize 
          cellCentre = subcellCentre(thisOctal,subCell)
          !      if (thisOctal%nDepth < 4) split = .true.
          !      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
          !      phi = atan2(cellcentre%y,cellcentre%x)
          !      warpHeight = 0. !cos(phi) * grid%rInner * sin(30.*degtorad) * sqrt(grid%rinner / r)
          do i = 1, grid%ng
             if (inOctal(thisOctal,grid%gArray(i)%centre).and.(cellSize > grid%gArray(i)%sigma)) then
                split = .true.
             endif
          enddo
          
!      if (r > grid%rInner*0.8) then
!         hr = height * (r / (100.*autocm/1.e10))**1.25
!         fac = cellsize/hr
!         if (abs((cellCentre%z-warpHeight)/hr) < 10.) then
!            if (fac > 2.) split = .true.
!         endif
!      endif

       case("pathtest")
          if (thisOctal%nDepth < minDepthAMR) split = .true.
          if ((modulus(subcellCentre(thisOctal, subcell)) < 0.2).and. &
               (thisOctal%nDepth<maxDepthAMR)) split = .true.
       case("toruslogo")
          !used to be 6
          if (thisOctal%nDepth  < 8) split = .true.
          
          if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 10.)) then
             split = .true.
             splitInAzimuth = .true.
          endif
          cellCentre = subcellCentre(thisOctal, subcell)
          r = sqrt(cellcentre%x**2 + cellcentre%y**2 + cellcentre%z**2)
          if ((r + thisOctal%subcellSize/2.d0)  .lt. 1.5*rsol/1.e10) then
             split = .false.
             splitinAzimuth = .false.
          endif
          if ((r - thisOctal%subcellSize/2.d0)  .gt. 2.*rsol/1.e10) then
             split = .false.
             splitinAzimuth = .false.
          endif
          
          !      if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 90.)) then
          !         split = .true.
          !         splitInAzimuth = .true.
          !      endif
          
          
       case("warpeddisc")
          
          cellSize = thisOctal%subcellSize 
          cellCentre = subcellCentre(thisOctal,subCell)
          r = sqrt(cellcentre%x**2 + cellcentre%y**2)
          phi1 = atan2(cellcentre%y,cellcentre%x)
          phi2 = phi1 - pi
          if (phi1 < 0.d0) phi1 = phi1 + twoPi
          if (phi2 < 0.d0) phi2 = phi2 + twoPi
          b = (1.d0/twoPi)*log(20.d0)
          warpRadius1 = rInner * exp(b * phi1)
          warpRadius2 = rInner * exp(b * phi2)
          
          warpheight1  = sin(0.5d0*phi1+warpAngle) * warpFracHeight * warpradius1 * exp(-0.5d0*((r - warpRadius1)/warpSigma)**2)
          warpheight2  = -sin(0.5d0*phi2+warpAngle) * warpFracHeight * warpradius2 * exp(-0.5d0*((r - warpRadius2)/warpSigma)**2)
          warpheight = warpheight1 + warpheight2
          
          hr = height * (r / rOuter)**betaDisc
          
          
          if ((r+cellsize/2.d0) > rInner*0.8) then
             hr = height * (r / rOuter)**betaDisc
             fac = cellsize/hr
             if (solveVerticalHydro) then
                ! split the grid in a way more similar to planetgap, to give hydro the best chance possible
                if ((abs((cellCentre%z-warpHeight)/hr) < 10.).and.(fac .gt. 0.5)) split = .true.
                if ((abs(cellcentre%z-warpheight)/hr > 5.).and.(abs((cellcentre%z-warpheight)/cellsize) < 2.)) split = .true.
             else
                if (hydroWarp) then
                   fac=hr
                   hr = getHydroWarpScaleHeight(r, grid)
                   !print *, fac, hr
                   if (.not. (hr > 0.)) then
                      hr = fac
                      !               if (r < rOuter*1.1d0) then
                      !                 print *, "Error: scaleheight returned by getHydroWarpScaleHeight = 0."
                      !                 print *, r, hr
                      !                 stop
                      !               end if
                   end if
                   fac = cellsize/hr
                   ! o The first criterion of this if (cell is less than 5 scaleheights from midplane) is fine because
                   !   it overestimates the height above the midplane that needs to be split.
                   ! o The second criterion will not split some cells that should be split (the new scaleheight is smaller,
                   !   so fac will be smaller than it should be and may not be > 1). To compensate, we change 1. to 0.75.
                   ! If we were to do this properly, we would find out what the new scale height of the disc is...
                end if
                if ((abs((cellCentre%z-warpHeight)/hr) < 5.).and.(fac .gt. 1.)) split = .true.
                if ((abs(cellcentre%z-warpheight)/hr > 5.).and.(abs((cellcentre%z-warpheight)/cellsize) < 1.)) split = .true.
                if ((abs(r - warpradius) < warpsigma).and.(abs((cellCentre%z-warpHeight)/hr) < 8.).and.(fac .gt. 0.5)) &
                     split = .true.
             end if
          endif
          
          
          if ((abs(r-rinner) < 0.9*rinner).and.(cellSize > 0.02*rInner).and.(abs(cellCentre%z) < 2.*hr)) then
             split = .true.
          endif

          if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 15.)) then
             split = .true.
             splitInAzimuth = .true.
          endif
          
          if (r > rOuter*1.1d0) then
             split = .false.
             splitInAzimuth = .false.
          endif
          
       case("magstream")
          
          
          IF (thisOctal%nDepth < 5) THEN 
             split = .TRUE.
          endif
!      ELSE
!
!        nsample = 1
!        ! the density is only sampled at the centre of the grid
!        ! we will search in each subcell to see if any point exceeds the 
!        ! threshold density
!
!        ! get the size and centre of the current cell
!        cellSize = thisOctal%subcellSize 
!        cellCentre = subcellCentre(thisOctal,subCell)
!    
!        ! check the density of random points in the current cell - 
!        !   if any of them exceed the critical density, set the flag
!        !   indicating the cell is to be split, and return from the function
!        split = .FALSE.
!        ave_density = 0.0_db
!        DO i = 1, nsample
!          searchPoint = cellCentre
!          searchPoint%x = searchPoint%x - (cellSize / 2.0_oc) + cellSize*REAL(x,KIND=oct) 
!          searchPoint%y = searchPoint%y - (cellSize / 2.0_oc) + cellSize*REAL(y,KIND=oct) 
!          searchPoint%z = searchPoint%z - (cellSize / 2.0_oc) + cellSize*REAL(z,KIND=oct)
!
!          ! using a generic density function in density_mod.f90
!          ave_density  = density(searchPoint,grid,timeNow) + ave_density
!
!        END DO
!
!        ave_density = ave_density / REAL(nsample,KIND=double)
!        total_mass = ave_density * (cellSize*1.e10_db)**3.0_db
!        IF (total_mass > amrLimitScalar) then
!           split = .TRUE.
!        END IF
!
!      END IF


       end select

222 continue
    end if

   if (thisOctal%nDepth == maxDepthAmr) then
      split = .false.
      if (firstTime) then
         if (minDepthAmr == maxDepthAmr) then 
            call writeInfo("Splitting grid to uniform depth")
         else
            call writeWarning("AMR cell depth capped")
         endif
         firstTime = .false.
      endif
   endif

   if (thisOctal%nDepth < minDepthAmr) then

      split = .true.
   endif

   if (grid%splitOverMPI) then
      if ((thisOctal%twoD.and.nHydroThreadsGlobal==4).or. &
           (thisOctal%threed.and.nHydroThreadsGlobal==8).or. &
           (thisOctal%oned  .and.nHydroThreadsGlobal==2) ) then
         if ((thisOctal%nDepth == 1).and.(subcell /= myRankGlobal)) then
            split = .false.
         endif
         if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
            split = .false.
         endif
      else
         if (thisOctal%oneD) then
            if (thisOctal%nDepth == 1) then
               split = .true.
            else if (thisOctal%nDepth > 1) then
               if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
                  split = .false.
               endif
            endif
         endif
         if (thisOctal%twoD) then
            if (thisOctal%nDepth == 1) then
               split = .true.
            else if (thisOctal%nDepth > 1) then
               if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
                  split = .false.
               endif
            endif
         endif

         if (thisOctal%threeD) then
            if (thisOctal%nDepth == 1) then
               split = .true.
            else if (thisOctal%nDepth > 1) then
               if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
                  split = .false.
               endif
            endif
         endif
      endif

      if (thisOctal%twod.and.(nHydroThreadsGlobal==16)) then
         if (thisOctal%nDepth <= 2) then
            split = .true.
         else
            if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
               split = .false.
            endif
         endif
      endif


      if (thisOctal%threed.and.(nHydroThreadsGlobal==64)) then
         if (thisOctal%nDepth <= 2) then
            split = .true.
         else
            if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
               split = .false.
            endif
         endif
      endif

      if (thisOctal%threed.and.(nHydroThreadsGlobal==512)) then
         if (thisOctal%nDepth <= 3) then
            split = .true.
         else
            if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
               split = .false.
            endif
         endif
      endif

      if (thisOctal%twod.and.(nHydroThreadsGlobal==64)) then
         if (thisOctal%nDepth <= 3) then
            split = .true.
         else
            if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
               split = .false.
            endif
         endif
      endif

   endif

101 continue

  END FUNCTION decideSplit

integer FUNCTION halfRefined(minDepthAMR, maxDepthAMR)
  integer :: minDepthAMR, maxDepthAMR

  if(minDepthAMR == maxDepthAMR) then
     halfRefined = minDepthAMR
  else
     halfRefined = int((minDepthAMR + maxDepthAMR)/2)
  end if
  
end FUNCTION halfRefined

logical  FUNCTION cornerCell(grid, thisOctal, subcell)
      type(OCTAL) :: thisOctal
      type(GRIDTYPE) :: grid
      integer :: subcell
      type(VECTOR) :: probe(6), rVec, locator
      integer :: nProbeOutside, nProbes, iProbe
      

      if (thisOctal%oned) then
         nProbes = 2
         probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
         probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
      endif
      if (thisOctal%twod) then
         nProbes = 4
         probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
         probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
         probe(3) = VECTOR(0.d0, 0.d0, 1.d0)
         probe(4) = VECTOR(0.d0, 0.d0, -1.d0)
      endif
      if (thisOctal%threeD) then
         nProbes = 6
         probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
         probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
         probe(3) = VECTOR(0.d0, 1.d0, 0.d0)
         probe(4) = VECTOR(0.d0, -1.d0, 0.d0)
         probe(5) = VECTOR(0.d0, 0.d0, 1.d0)
         probe(6) = VECTOR(0.d0, 0.d0, -1.d0)
      endif

      rVec = subcellCentre(thisOctal, subcell)
      nProbeOutside = 0
      do iProbe = 1, nProbes
         locator = rVec + &
         (thisOctal%subcellsize/2.d0 + 0.01d0*grid%halfSmallestSubcell)*probe(iProbe)                 
         if (.not.inOctal(grid%octreeRoot, locator).or. &
         (locator%x < (grid%octreeRoot%centre%x-grid%octreeRoot%subcellSize))) &                      
         nProbeOutside = nProbeOutside + 1
      enddo
      cornerCell=.false.
      if (thisOctal%twoD.and.(nProbeOutside > 1)) cornerCell = .true.
      if (thisOctal%threeD.and.(nProbeOutside > 2)) cornerCell = .true.

      
  END FUNCTION cornerCell


logical  FUNCTION edgeCell(grid, thisOctal, subcell)
      type(OCTAL) :: thisOctal
      type(GRIDTYPE) :: grid
      integer :: subcell
      type(VECTOR) :: probe(6), rVec, locator
      integer :: nProbeOutside, nProbes, iProbe


      if (thisOctal%oned) then
         nProbes = 2
         probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
         probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
      endif
      if (thisOctal%twod) then
         nProbes = 4
         probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
         probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
         probe(3) = VECTOR(0.d0, 0.d0, 1.d0)
         probe(4) = VECTOR(0.d0, 0.d0, -1.d0)
      endif
      if (thisOctal%threeD) then
         nProbes = 6
         probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
         probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
         probe(3) = VECTOR(0.d0, 1.d0, 0.d0)
         probe(4) = VECTOR(0.d0, -1.d0, 0.d0)
         probe(5) = VECTOR(0.d0, 0.d0, 1.d0)
         probe(6) = VECTOR(0.d0, 0.d0, -1.d0)
      endif

      rVec = subcellCentre(thisOctal, subcell)
      nProbeOutside = 0
      do iProbe = 1, nProbes
         locator = rVec + &
         (thisOctal%subcellsize/2.d0 + 0.01d0*grid%halfSmallestSubcell)*probe(iProbe)                                              
         if (.not.inOctal(grid%octreeRoot, locator).or. &
         (locator%x < (grid%octreeRoot%centre%x-grid%octreeRoot%subcellSize))) &                                                   
          nProbeOutside = nProbeOutside + 1
      enddo
      edgeCell=.false.
      if (thisOctal%twoD.and.(nProbeOutside >= 1)) edgeCell = .true.
      if (thisOctal%threeD.and.(nProbeOutside >= 1)) edgeCell = .true.

  END FUNCTION edgeCell


logical  FUNCTION ghostCell(grid, thisOctal, subcell)
      type(OCTAL) :: thisOctal
      type(GRIDTYPE) :: grid
      integer :: subcell
      type(VECTOR) :: probe(6), rVec, locator
      integer :: nProbeOutside, nProbes, iProbe


      if (thisOctal%oned) then
         nProbes = 2
         probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
         probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
      endif
      if (thisOctal%twod) then
         nProbes = 4
         probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
         probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
         probe(3) = VECTOR(0.d0, 0.d0, 1.d0)
         probe(4) = VECTOR(0.d0, 0.d0, -1.d0)
      endif
      if (thisOctal%threeD) then
         nProbes = 6
         probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
         probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
         probe(3) = VECTOR(0.d0, 1.d0, 0.d0)
         probe(4) = VECTOR(0.d0, -1.d0, 0.d0)
         probe(5) = VECTOR(0.d0, 0.d0, 1.d0)
         probe(6) = VECTOR(0.d0, 0.d0, -1.d0)
      endif

      rVec = subcellCentre(thisOctal, subcell)
      nProbeOutside = 0
      do iProbe = 1, nProbes
         locator = rVec + &
         (thisOctal%subcellsize + 0.01d0*grid%halfSmallestSubcell)*probe(iProbe)
         if (.not.inOctal(grid%octreeRoot, locator).or. &
         (locator%x < (grid%octreeRoot%centre%x-grid%octreeRoot%subcellSize))) &
          nProbeOutside = nProbeOutside + 1
      enddo
      ghostCell=.false.
      if (thisOctal%twoD.and.(nProbeOutside >= 1)) ghostCell = .true.
      if (thisOctal%threeD.and.(nProbeOutside >= 1)) ghostCell = .true.


  END FUNCTION ghostCell
  
  SUBROUTINE fillVelocityCorners(thisOctal,velocityFunc, debug)
    ! store the velocity values at the subcell corners of an octal so
    !   that they can be used for interpolation.

    IMPLICIT NONE
  
    TYPE(octal), INTENT(INOUT) :: thisOctal
    logical, optional :: debug
    logical :: writedebug
    real(oct)      :: r1, r2, r3
    real(oct)      :: phi1, phi2, phi3
    real(oct)      :: x1, x2, x3
    real(oct)      :: y1, y2, y3
    real(oct)      :: z1, z2, z3
    
    INTERFACE 
      TYPE(vector) FUNCTION velocityFunc(point)
        USE vector_mod
        USE gridtype_mod
        TYPE(vector), INTENT(IN) :: point
      END FUNCTION velocityFunc
    END INTERFACE

    writedebug = .false.
    if (present(debug)) writedebug=debug

    if (thisOctal%oneD) then
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       y1 = 0.d0
       z1 = 0.d0
       thisOctal%cornerVelocity(1) = velocityFunc(vector(x1,y1,z1))
       thisOctal%cornerVelocity(2) = velocityFunc(vector(x2,y1,z1))
       thisOctal%cornerVelocity(3) = velocityFunc(vector(x3,y1,z1))
       goto 666
    endif


    if (thisOctal%threed) then
       if (.not.thisOctal%cylindrical) then ! 3d cartesian case
          ! we first store the values we use to assemble the position vectors
          
          x1 = thisOctal%centre%x - thisOctal%subcellSize
          x2 = thisOctal%centre%x
          x3 = thisOctal%centre%x + thisOctal%subcellSize
          
          y1 = thisOctal%centre%y - thisOctal%subcellSize
          y2 = thisOctal%centre%y
          y3 = thisOctal%centre%y + thisOctal%subcellSize
          
          z1 = thisOctal%centre%z - thisOctal%subcellSize
          z2 = thisOctal%centre%z
          z3 = thisOctal%centre%z + thisOctal%subcellSize
                    
           ! now store the 'base level' values
          
          thisOctal%cornerVelocity(1) = velocityFunc(vector(x1,y1,z1))
          thisOctal%cornerVelocity(2) = velocityFunc(vector(x2,y1,z1))
          thisOctal%cornerVelocity(3) = velocityFunc(vector(x3,y1,z1))
          thisOctal%cornerVelocity(4) = velocityFunc(vector(x1,y2,z1))
          thisOctal%cornerVelocity(5) = velocityFunc(vector(x2,y2,z1))
          thisOctal%cornerVelocity(6) = velocityFunc(vector(x3,y2,z1))
          thisOctal%cornerVelocity(7) = velocityFunc(vector(x1,y3,z1))
          thisOctal%cornerVelocity(8) = velocityFunc(vector(x2,y3,z1))
          thisOctal%cornerVelocity(9) = velocityFunc(vector(x3,y3,z1))
          
          ! middle level
          
          thisOctal%cornerVelocity(10) = velocityFunc(vector(x1,y1,z2))
          thisOctal%cornerVelocity(11) = velocityFunc(vector(x2,y1,z2))
          thisOctal%cornerVelocity(12) = velocityFunc(vector(x3,y1,z2))
          thisOctal%cornerVelocity(13) = velocityFunc(vector(x1,y2,z2))
          thisOctal%cornerVelocity(14) = velocityFunc(vector(x2,y2,z2))
          thisOctal%cornerVelocity(15) = velocityFunc(vector(x3,y2,z2))
          thisOctal%cornerVelocity(16) = velocityFunc(vector(x1,y3,z2))
          thisOctal%cornerVelocity(17) = velocityFunc(vector(x2,y3,z2))
          thisOctal%cornerVelocity(18) = velocityFunc(vector(x3,y3,z2))
          
          ! top level
          
          thisOctal%cornerVelocity(19) = velocityFunc(vector(x1,y1,z3))
          thisOctal%cornerVelocity(20) = velocityFunc(vector(x2,y1,z3))
          thisOctal%cornerVelocity(21) = velocityFunc(vector(x3,y1,z3))
          thisOctal%cornerVelocity(22) = velocityFunc(vector(x1,y2,z3))
          thisOctal%cornerVelocity(23) = velocityFunc(vector(x2,y2,z3))
          thisOctal%cornerVelocity(24) = velocityFunc(vector(x3,y2,z3))
          thisOctal%cornerVelocity(25) = velocityFunc(vector(x1,y3,z3))
          thisOctal%cornerVelocity(26) = velocityFunc(vector(x2,y3,z3))
          thisOctal%cornerVelocity(27) = velocityFunc(vector(x3,y3,z3))

       else ! cylindrical 
          if (thisOctal%splitAzimuthally) then
             z1 = thisOctal%centre%z - thisOctal%subcellSize
             z2 = thisOctal%centre%z
             z3 = thisOctal%centre%z + thisOctal%subcellSize
             phi1 = thisOctal%phiMin
             phi2 = thisOctal%phi 
             phi3 = thisOctal%phiMax
             r1 = thisOctal%r - thisOctal%subcellSize
             r2 = thisOctal%r
             r3 = thisOctal%r + thisOctal%subcellSize

             ! bottom level

             thisOctal%cornerVelocity(1) = velocityFunc(vector(r1*cos(phi1),r1*sin(phi1),z1))
             thisOctal%cornerVelocity(2) = velocityFunc(vector(r2*cos(phi1),r2*sin(phi1),z1))
             thisOctal%cornerVelocity(3) = velocityFunc(vector(r3*cos(phi1),r3*sin(phi1),z1))
             thisOctal%cornerVelocity(4) = velocityFunc(vector(r1*cos(phi2),r1*sin(phi2),z1))
             thisOctal%cornerVelocity(5) = velocityFunc(vector(r2*cos(phi2),r2*sin(phi2),z1))
             thisOctal%cornerVelocity(6) = velocityFunc(vector(r3*cos(phi2),r3*sin(phi2),z1))
             thisOctal%cornerVelocity(7) = velocityFunc(vector(r1*cos(phi3),r1*sin(phi3),z1))
             thisOctal%cornerVelocity(8) = velocityFunc(vector(r2*cos(phi3),r2*sin(phi3),z1))
             thisOctal%cornerVelocity(9) = velocityFunc(vector(r3*cos(phi3),r3*sin(phi3),z1))

             thisOctal%cornerVelocity(10) = velocityFunc(vector(r1*cos(phi1),r1*sin(phi1),z2))
             thisOctal%cornerVelocity(11) = velocityFunc(vector(r2*cos(phi1),r2*sin(phi1),z2))
             thisOctal%cornerVelocity(12) = velocityFunc(vector(r3*cos(phi1),r3*sin(phi1),z2))
             thisOctal%cornerVelocity(13) = velocityFunc(vector(r1*cos(phi2),r1*sin(phi2),z2))
             thisOctal%cornerVelocity(14) = velocityFunc(vector(r2*cos(phi2),r2*sin(phi2),z2))
             thisOctal%cornerVelocity(15) = velocityFunc(vector(r3*cos(phi2),r3*sin(phi2),z2))
             thisOctal%cornerVelocity(16) = velocityFunc(vector(r1*cos(phi3),r1*sin(phi3),z2))
             thisOctal%cornerVelocity(17) = velocityFunc(vector(r2*cos(phi3),r2*sin(phi3),z2))
             thisOctal%cornerVelocity(18) = velocityFunc(vector(r3*cos(phi3),r3*sin(phi3),z2))

             thisOctal%cornerVelocity(19) = velocityFunc(vector(r1*cos(phi1),r1*sin(phi1),z3))
             thisOctal%cornerVelocity(20) = velocityFunc(vector(r2*cos(phi1),r2*sin(phi1),z3))
             thisOctal%cornerVelocity(21) = velocityFunc(vector(r3*cos(phi1),r3*sin(phi1),z3))
             thisOctal%cornerVelocity(22) = velocityFunc(vector(r1*cos(phi2),r1*sin(phi2),z3))
             thisOctal%cornerVelocity(23) = velocityFunc(vector(r2*cos(phi2),r2*sin(phi2),z3))
             thisOctal%cornerVelocity(24) = velocityFunc(vector(r3*cos(phi2),r3*sin(phi2),z3))
             thisOctal%cornerVelocity(25) = velocityFunc(vector(r1*cos(phi3),r1*sin(phi3),z3))
             thisOctal%cornerVelocity(26) = velocityFunc(vector(r2*cos(phi3),r2*sin(phi3),z3))
             thisOctal%cornerVelocity(27) = velocityFunc(vector(r3*cos(phi3),r3*sin(phi3),z3))


          else

             z1 = thisOctal%centre%z - thisOctal%subcellSize
             z2 = thisOctal%centre%z
             z3 = thisOctal%centre%z + thisOctal%subcellSize
             phi1 = thisOctal%phi - thisOctal%dPhi/2.d0
             phi2 = thisOctal%phi + thisOctal%dPhi/2.d0
             r1 = thisOctal%r - thisOctal%subcellSize
             r2 = thisOctal%r
             r3 = thisOctal%r + thisOctal%subcellSize


             ! bottom level

             thisOctal%cornerVelocity(1) = velocityFunc(vector(r1*cos(phi1),r1*sin(phi1),z1))
             thisOctal%cornerVelocity(2) = velocityFunc(vector(r1*cos(phi1),r1*sin(phi1),z1))
             thisOctal%cornerVelocity(3) = velocityFunc(vector(r2*cos(phi1),r2*sin(phi1),z1))
             thisOctal%cornerVelocity(4) = velocityFunc(vector(r2*cos(phi2),r2*sin(phi2),z1))
             thisOctal%cornerVelocity(5) = velocityFunc(vector(r3*cos(phi2),r3*sin(phi2),z1))
             thisOctal%cornerVelocity(6) = velocityFunc(vector(r3*cos(phi2),r3*sin(phi2),z1))

             ! middle level

             thisOctal%cornerVelocity(7) = velocityFunc(vector(r1*cos(phi1),r1*sin(phi1),z2))
             thisOctal%cornerVelocity(8) = velocityFunc(vector(r1*cos(phi1),r1*sin(phi1),z2))
             thisOctal%cornerVelocity(9) = velocityFunc(vector(r2*cos(phi1),r2*sin(phi1),z2))
             thisOctal%cornerVelocity(10) = velocityFunc(vector(r2*cos(phi2),r2*sin(phi2),z2))
             thisOctal%cornerVelocity(11) = velocityFunc(vector(r3*cos(phi2),r3*sin(phi2),z2))
             thisOctal%cornerVelocity(12) = velocityFunc(vector(r3*cos(phi2),r3*sin(phi2),z2))

             ! top level

             thisOctal%cornerVelocity(13) = velocityFunc(vector(r1*cos(phi1),r1*sin(phi1),z3))
             thisOctal%cornerVelocity(14) = velocityFunc(vector(r1*cos(phi1),r1*sin(phi1),z3))
             thisOctal%cornerVelocity(15) = velocityFunc(vector(r2*cos(phi1),r2*sin(phi1),z3))
             thisOctal%cornerVelocity(16) = velocityFunc(vector(r2*cos(phi2),r2*sin(phi2),z3))
             thisOctal%cornerVelocity(17) = velocityFunc(vector(r3*cos(phi2),r3*sin(phi2),z3))
             thisOctal%cornerVelocity(18) = velocityFunc(vector(r3*cos(phi2),r3*sin(phi2),z3))

          endif
       endif
    else       
       
    ! we first store the values we use to assemble the position vectors
       
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       
       z1 = thisOctal%centre%z - thisOctal%subcellSize
       z2 = thisOctal%centre%z
       z3 = thisOctal%centre%z + thisOctal%subcellSize
       
       ! now store the 'base level' values
       
       thisOctal%cornerVelocity(1) = velocityFunc(vector(x1,0.d0,z1))
       thisOctal%cornerVelocity(2) = velocityFunc(vector(x2,0.d0,z1))
       thisOctal%cornerVelocity(3) = velocityFunc(vector(x3,0.d0,z1))
       thisOctal%cornerVelocity(4) = velocityFunc(vector(x1,0.d0,z2))
       thisOctal%cornerVelocity(5) = velocityFunc(vector(x2,0.d0,z2))
       thisOctal%cornerVelocity(6) = velocityFunc(vector(x3,0.d0,z2))
       thisOctal%cornerVelocity(7) = velocityFunc(vector(x1,0.d0,z3))
       thisOctal%cornerVelocity(8) = velocityFunc(vector(x2,0.d0,z3))
       thisOctal%cornerVelocity(9) = velocityFunc(vector(x3,0.d0,z3))
    endif
666 continue

!    if(isnan(thisOctal%cornerVelocity(1)%x)) then
!          write(*,*) "cornervel",thisOctal%cornerVelocity(1)
!          write(*,*) velocityFunc(vector(x1,0.d0,z1),grid)
!          write(*,*) x1,z1
!          write(*,*) x2,z2
!          write(*,*) x3,z3
!       enddo
!    endif
    
  END SUBROUTINE fillVelocityCorners

  TYPE(vector) FUNCTION TTauriVelocity(point)
    ! calculates the velocity vector at a given point for a model
    !   of a T Tauri star with magnetospheric accretion
    ! see Hartman, Hewett & Calvet 1994ApJ...426..669H 

    use inputs_mod, only: TTauriRinner, TTauriRouter, TTauriRstar, &
                               TTauriMstar, TTauriDiskHeight, dipoleoffset
                               
    IMPLICIT NONE

    TYPE(vector), INTENT(IN) :: point

    TYPE(vector) :: pointVec
    TYPE(vector)      :: vP
    REAL(double)  :: modVp
    REAL(double)  :: phi
    REAL(double)  :: r, rM, theta, y

    pointVec = point  * 1.e10_oc
    r = modulus( pointVec ) 
    if (r /= 0.d0) then
       theta = ACOS( max(-1.d0,min(1.d0,pointVec%z / dble(r) )))
    else
       theta = pi/2.d0
    endif
    if (theta == 0.d0) theta=1.d-20
    rM  = r / SIN(theta)**2
    y = SIN(theta)**2 

    ! test if the point lies within the star
    IF ( r < TTauriRstar ) THEN
      TTauriVelocity = vector(1.e-15,1.e-15,1.e-15)
      
    ! test if the point lies too close to the disk
    ELSE IF ( ABS(pointVec%z) < TTauriDiskHeight) THEN
      TTauriVelocity = vector(1.e-14,1.e-14,1.e-14)
    ! test if the point lies outside the accretion stream
    ELSE IF ((rM > TTauriRinner) .AND. (rM < TTauriRouter )) THEN
  
      vP = vector(3.0 * SQRT(y) * SQRT(1.0-y) / SQRT(4.0 - (3.0*y)), &
                  0.0, &
                 (2.0 - 3.0 * y) / SQRT(4.0 - 3.0 * y))
      modVp = SQRT((2.0 * bigG * TTauriMstar / TTauriRstar) * &
                     (TTauriRstar/r - TTauriRstar/rM))
      vP = (-1.0 * (modVp/cSpeed)) * vP
      phi = ATAN2(pointVec%y,pointVec%x)
      vP = rotateZ(vP,-phi)
      
      IF (pointVec%z < 0.0) vP%z = -vP%z
      TTauriVelocity = vP

    ELSE
      TTauriVelocity = vector(1.e-12,1.e-12,1.e-12)
    END IF
    TTauriVelocity = rotateY(TTauriVelocity,-dble(dipoleOffset))
  END FUNCTION TTauriVelocity

  SUBROUTINE calcTTauriMassVelocity(thisOctal,subcell,grid) 
    ! calculates some of the variables at a given point for a model
    !   of a T Tauri star with magnetospheric accretion
    ! see Hartman, Hewett & Calvet 1994ApJ...426..669H 

    USE inputs_mod, ONLY : useHartmannTemp, maxHartTemp, TTauriRstar,&
                                TTauriRinner, TTauriRouter, ttau_acc_on, &
                                 vturb, ttauriMagnetosphere
    use magnetic_mod, only : inflowMahdavi, velocityMahdavi
    use density_mod, only:   density, TTauriInFlow

    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(INOUT) :: grid

    TYPE(vector) :: point

    TYPE(vector) :: starPosn
    TYPE(vector) :: pointVec

    REAL :: r, rM, theta, bigR, thetaStar
    REAL :: bigRstar, thetaStarHartmann, bigRstarHartmann, rMnorm
    REAL :: tmp

    starPosn = grid%starPos1

    point = subcellCentre(thisOctal,subcell)
    pointVec = (point - starPosn) * 1.e10_oc

! rotation for warp dipole offset!!!!!!! TJH
!    point = rotateY(point, dble(dipoleOffset))
!    pointvec = rotateY(pointvec, dble(dipoleOffset))



    r = real(modulus(point))
    if (r < (ttauriRstar/1.d10)) thisOctal%inflow(subcell) = .false.

    if (.not.ttauriMagnetosphere) then
    if (.not. ttau_acc_on) then
       ! we don't include the magnetopsherical accretion
       thisOctal%rho(subcell) = 1.e-19
       IF (useHartmannTemp) thisOctal%temperature(subcell) = 6500.0
    else
    
       if (inflowMahdavi(pointvec)) then
          thisOctal%rho(subcell) = 1.d-14
          thisOctal%velocity(subcell) = velocityMahdavi(pointVec/1.d10)
       else
          thisOctal%rho(subcell) = 1.d-25
          thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
       endif
       goto 666

       IF (TTauriInFlow(point,grid)) THEN
          thisOctal%rho(subcell) = Density(point,grid)
          thisOctal%inFlow(subcell) = .TRUE.
          IF (useHartmannTemp) THEN
             ! need to calculate the flow point AS IF the magnetosphere
             !   was the same size as the one in the Hartmann et al paper
             r = real(modulus( pointVec ) )
             theta = real(acos( pointVec%z  / r ))
        
             rM  = r / SIN(theta)**2
        
             ! work out the intersection angle (at the stellar surface) for
             !   the current flow line  
             thetaStar = ASIN(SQRT(TTauriRstar/rM))
             bigRstar = TTauriRstar * SIN(thetaStar) 

             ! normalize the magnetic radius (0=inside, 1=outside)
             rMnorm = (rM-TTauriRinner) / (TTauriRouter-TTauriRinner)
        
             ! convert rM to Hartmann units
             rM = real(2.2*(2.0*rSol) + rMnorm*(0.8*(2.0*rSol)))
        
             bigR = real(SQRT(pointVec%x**2+pointVec%y**2))
             bigR = (bigR-bigRstar) / (TTauriRouter-bigRstar) ! so that we can rescale it
             bigR = min(bigR,TTauriRouter)
             bigR = max(bigR,0.0)

             ! get the equivalent bigR for the Hartmann geometry
             ! first work out the intersection angle (at the stellar surface) for
             !   the current flow line  
             thetaStarHartmann = real(ASIN(SQRT((2.0*rSol)/rM)))
             ! work out bigR for this point
             bigRstarHartmann = real((2.0*rSol) * SIN(thetaStarHartmann) )
             
             bigR = real(bigRstarHartmann + (bigR * (3.0*(2.0*rSol) - bigRstarHartmann)))
             
             r = (rM * bigR**2)**(1./3.) ! get r in the Hartmann setup
             
             ! get theta in the Hartmann setup
             tmp = r/rM
             ! quick fix for out of range argments (R. Kurosawa)
             if (tmp > 1.0) tmp=1.0
             if (tmp < -1.0) tmp = -1.0
             theta = ASIN(SQRT(tmp))
             theta = MIN(theta,REAL(pi/2.0-0.0001))
             theta = MAX(theta,0.0001)
        
             rMnorm = MAX(rMnorm, 0.0001)
             rMnorm = MIN(rMnorm, 0.9999)
             !        thisOctal%temperature(subcell) = hartmannTemp(rMnorm, theta, maxHartTemp)
             thisOctal%temperature(subcell) = hartmannTemp2(rMnorm, theta, maxHartTemp)
          END IF
       ELSE  ! not inflow
          !     thisOctal%rho(subcell) = 1.e-25
          thisOctal%rho(subcell) = 1.e-19
          IF (useHartmannTemp) thisOctal%temperature(subcell) = 6500.0
       END IF
    end if ! ttau_acc_on

  thisOctal%velocity(subcell) = TTauriVelocity(point)
  !thisOctal%velocity(subcell) = TTauriRotation(point,grid)    
endif
  if (associated(thisOctal%microturb)) thisOctal%microturb(subcell) = vturb
  IF ((thisoctal%threed).and.(subcell == 8)) &
       CALL fillVelocityCorners(thisOctal,TTauriVelocity)
  
  IF ((thisoctal%twod).and.(subcell == 4)) &
       CALL fillVelocityCorners(thisOctal,TTauriVelocity)
  
666 continue
  END SUBROUTINE calcTTauriMassVelocity

  !
  ! The function to check if a given point (rvec) is in "fuzzy" edge region
  logical function  in_fuzzy_edge(rvec) 
    use inputs_mod, only: TTauriRinner, TTauriRouter

    IMPLICIT NONE    
    TYPE(vector), intent(in):: rvec
    real(double) :: r, theta, rM,  w
    real(double) :: rM_fuzzy_in, rM_fuzzy_out  ! beginning of the fuzzy edges
    !
    !
    r = modulus(rvec)
    
    if (r/=0.0d0) then
       theta = ACOS( MIN(ABS(rVec%z/r),0.998_oc) )
    else
       theta=0.01
    end if
    
    if (theta == 0.0d0) theta = 0.01
    rM  = r / SIN(theta)**2
    
    
    ! The fuzzy density starts from a 5-th of the thickness (2h) 
    ! below the surface.
    w = get_fuzzy_width() ! [10^10cm] using a fucntion in this module

    rM_fuzzy_in  = TTauriRinner*1.0d-10 + w   ! [10^10cm]
    rM_fuzzy_out = TTauriRouter*1.0d-10 - w   ! [10^10cm]
             
    ! If the point is close to the edge, we make it fuzzy
    if ( rM < rM_fuzzy_in) then
       in_fuzzy_edge = .true.
    elseif (rM > rM_fuzzy_out) then
       in_fuzzy_edge = .true.
    else
       in_fuzzy_edge = .false.
    end if

  end function in_fuzzy_edge

  !
  ! The function to check if a given point (rvec) is in "fuzzy" edge region
  function  get_fuzzy_width()  RESULT(w)
    use inputs_mod, only: TTauriRinner, TTauriRouter

    IMPLICIT NONE    
    real(double) :: w  ! width in  [10^10cm]
    real(double) :: rM_center, h
    !
    ! bias more toward the edges of the accreation stream
    h = 0.5d0*(TTauriRouter - TTauriRinner)*1.0d-10  ! [10^10cm] distance from the centere to edge
    rM_center = 0.5d0*(TTauriRouter + TTauriRinner)*1.0d-10 ! [10^10cm]   
    
    ! The fuzzy density starts from a 5-th of the thickness (2h) 
    ! below the surface.
             w = 0.1d0*h;  ! Halpha081 model
!    w = 0.2d0*h;  ! a fifth for now  ! Halpha067 model
!             w = 0.25d0*h;  ! Halpha083 model
!             w = 0.3d0*h;  ! Halpha080 model


  end function get_fuzzy_width

  subroutine infallEnhancmentAMR(grid, distortionVector, nVec, timeStep, doDistortion, &
                                 particleMass, alreadyDoneInfall)
    ! creates an increase in the accretion rate for a T Tauri geometry.
    ! currently, the setup must be simple (no disc inclination etc.)
                                 
    type(GRIDTYPE), intent(inout) :: grid
    integer, intent(in) :: nVec
    logical, intent(in) :: doDistortion
    type(VECTOR), intent(inout) :: distortionVector(nVec)
    real, intent(in) :: timeStep
    real, intent(in) :: particleMass
    logical, intent(inout) :: alreadyDoneInfall

    real(oct) :: dTime
    type(vector) :: distortionVec(nVec)
    type(vector) :: thisVec
    type(vector) :: thisVel
    type(vector) :: thisVelOc
    logical :: setAllChanged
    
    integer, parameter :: nTimes = 1000
    integer :: i, j
    type(vector) :: starPos

    ! if we are not undoing a previous infall, we will assume that we need to
    ! run stateq on ALL subcells
    setAllChanged = .not. alreadyDoneInfall
    
    do i = 1, nVec
       distortionVec(i) = distortionVector(i)
    end do

    if (alreadyDoneInfall) then 
      ! we first have to undo the previous infall phase's changes
      call infallEnhancePrivate(grid%octreeRoot,distortionVec,undoPrevious=.true.,&
                                setAllChanged=.false.)
    end if
    
    write(*,*) "Time stepping vectors..."
    dTime = timeStep/real(nTimes,kind=db)
    do j = 1, nVec
      do i = 1, nTimes
        thisVec = distortionVec(j) / 1.e10_oc
        call amrGridValues(grid%octreeRoot,thisVec,velocity=thisVel)
        thisVelOc = thisVel
        thisVelOc = cSpeed * thisVel
        distortionVec(j) = distortionVec(j) + (dTime * thisVelOc)
      enddo
    enddo
    write(*,*) "done."
    
    do i = 1, nVec
       distortionVector(i) = distortionVec(i)
    end do
    
    if (doDistortion) then
      write(*,*) "Distorting grid..."    
      starPos = grid%starPos1
      call infallEnhancePrivate(grid%octreeRoot,distortionVec, undoPrevious=.false.,&
                                setAllChanged=(.not.alreadyDoneInfall))
      alreadyDoneInfall = .true.
      write(*,*) "done."
    end if
    
  contains

    recursive subroutine infallEnhancePrivate(thisOctal, distortionVec, undoPrevious,&
                                              setAllChanged)

      type(octal), pointer :: thisOctal
      type(vector), intent(in) :: distortionVec(nVec)
      logical, intent(in) :: undoPrevious
      logical, intent(in) :: setAllChanged
      
      integer, parameter :: nPhi = 360
      real(oct) :: phi
      integer :: containedParticles
      real :: addedDensity
      real :: deltaRho
      integer :: iSubcell, iChild
      type(octal), pointer :: child
      type(vector) :: trueVector
      type(vector) :: thisVec, thisVec2
 
      do iSubcell = 1, thisOctal%maxChildren, 1
        if (thisOctal%hasChild(iSubcell)) then
          
          ! descend to a child
          
          do iChild = 1, thisOctal%nChildren, 1
            if (thisOctal%indexChild(iChild) == iSubcell) exit
          end do 

          child => thisOctal%child(iChild)
          call infallEnhancePrivate(child, distortionVec, undoPrevious,&
                                    setAllChanged)
          
        else

          if (setAllChanged) thisOctal%changed(:) = .true.
          
          ! if we are out of the accretion flow, we can ignore this cell
          if (thisOctal%rho(iSubcell) < 1.e-24) cycle
          
          ! check if any of the 'particles' lie in the current subcell

          containedParticles = 0 
          do i = 1, nVec
            do j = 1, nPhi

              phi = twoPi * real(j-1) / real(nPhi-1)
              thisVec = rotateZ(distortionVec(i), phi)
              thisVec = thisVec * 1.0e-10_oc
              thisVec = thisVec + starPos

              if (inOctal(thisOctal,thisVec)) then
                 if (thisOctal%twod) then
                    thisVec2 = projecttoXZ(thisVec)
                 else
                    thisVec2 = thisVec
                 endif
                if (whichSubcell(thisOctal,thisVec2) == iSubcell) &
                  containedParticles = containedParticles + 1
                  
              else
                trueVector = subcellCentre(thisOctal,iSubcell) - starPos
                
                if (trueVector%z < 0.0_oc) then 

                  ! mirror the particles about the disc, and try again 
                  thisVec%z = starPos%z - (thisVec%z - starPos%z)
                  
                  if (inOctal(thisOctal,thisVec)) then
                 if (thisOctal%twod) then
                    thisVec2 = projecttoXZ(thisVec)
                 else
                    thisVec2 = thisVec
                 endif
                    if (whichSubcell(thisOctal,thisVec2) == iSubcell) &
                      containedParticles = containedParticles + 1
                  end if
                end if

              end if

            end do
          end do

          if (containedParticles > 0) then 
            
            addedDensity = real((real(containedParticles) * particleMass) / &
                            (real(thisOctal%subcellSize,kind=db)*1.e10_db)**3)
            ! calculate the density change in the subcell
            deltaRho = real(addedDensity / thisOctal%rho(iSubcell))

            if (.not. undoPrevious) then
              thisOctal%rho(iSubcell) = thisOctal%rho(iSubcell) + addedDensity
              thisOctal%etaLine(iSubcell) = thisOctal%etaLine(iSubcell) * (1.0+deltaRho)**2
              thisOctal%chiLine(iSubcell) = thisOctal%chiLine(iSubcell) * (1.0+deltaRho)**2
            else  
              ! we are reversing the changes made on a previous run of this
              !   subroutine.
              thisOctal%rho(iSubcell) = thisOctal%rho(iSubcell) - addedDensity
              thisOctal%etaLine(iSubcell) = thisOctal%etaLine(iSubcell) / (1.0+deltaRho)**2
              thisOctal%chiLine(iSubcell) = thisOctal%chiLine(iSubcell) / (1.0+deltaRho)**2
            end if

            thisOctal%changed(iSubcell) = .true.
              
          endif
        end if
      end do
    
    end subroutine infallEnhancePrivate

  end subroutine infallEnhancmentAMR  

  !
  ! This is based on hartmannTemp. The temperature is interpolarted
  ! by using all five field lines instead of just two to avoid 
  ! the jumps in the temperature across the field lines. (RK)
  FUNCTION hartmannTemp2(rM,inTheta,maxHartTemp)
    ! the paper: Hartman, Hewett & Calvet 1994ApJ...426..669H 
    !   quotes data from points along the magnetic field lines marked 
    !   on one of the figures. this function returns the temperature
    !   at any point within the accretion flow, by interpolating
    !   between the values in the published figure 6.

    USE unix_mod, only: unixGetEnv
    USE utils_mod, only: polint

    REAL             :: hartmannTemp2
    REAL, INTENT(IN) :: rM ! normalized (0=inside, 1=outside)
    REAL, INTENT(IN) :: inTheta
    REAL, INTENT(IN) :: maxHartTemp
    
    
    REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: temperatures
    REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: angles
    INTEGER, DIMENSION(5), SAVE             :: nSamples
    
    LOGICAL, SAVE:: alreadyLoaded = .FALSE.
    REAL :: radius
    REAL :: theta
    CHARACTER(LEN=80) :: dataDirectory 
    INTEGER :: errNo, i
    
    INTEGER :: subPos
    INTEGER :: fieldline
    INTEGER :: lowerLine, upperLine
    INTEGER :: iSample
    INTEGER :: maxSamples
    INTEGER :: iStat
    logical, save :: warned_already_01 = .false.
    logical, save :: warned_already_02 = .false.
    real:: T(5), R(5), dT, R_ip
    integer :: indx_T(5)

    dt = 0.; hartmanntemp2 = 0.

    IF (.NOT. alreadyLoaded) THEN
      ! if this is the first time the function has been called, need
      !   to load in the data.
      
      call unixGetenv("TORUS_DATA",dataDirectory,i,errno)
      dataDirectory = trim(dataDirectory)//"/hartmann/"
      OPEN(31,FILE=TRIM(dataDirectory)//"tempProfile1.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF 
      OPEN(32,FILE=TRIM(dataDirectory)//"tempProfile2.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF 
      OPEN(33,FILE=TRIM(dataDirectory)//"tempProfile3.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF  
      OPEN(34,FILE=TRIM(dataDirectory)//"tempProfile4.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF  
      OPEN(35,FILE=TRIM(dataDirectory)//"tempProfile5.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF 

      DO fieldline = 1, 5, 1
        READ(UNIT=(30+fieldline),FMT=*) nSamples(fieldLine)
      END DO

      maxSamples = MAXVAL(nSamples)

      ALLOCATE(temperatures(5,maxSamples))
      ALLOCATE(angles(5,maxSamples))
      ! make sure that there are no uninitialized variables
      temperatures = 0.0
      angles = 0.0
      
      DO fieldline = 1, 5, 1
        DO iSample = 1,  nSamples(fieldLine), 1
          READ(UNIT=(30+fieldline),FMT=*) angles(fieldLine,iSample),     &
                                          temperatures(fieldLine,iSample)
        END DO
        CLOSE(UNIT=(30+fieldline))
      END DO


      ! might need to rescale the temperature distribution to have a new maximum
      temperatures = temperatures * maxHartTemp/MAXVAL(temperatures)
      
      alreadyLoaded = .TRUE.
    END IF

    ! find which fieldlines bracket the point
   
    radius = rM * 0.8 ! because the Hartmann magnetic radii span 0.8 R_star
    radius = radius / 0.2 ! divide by the field line spacing

    IF (radius < 0.0) THEN 
       if (.not. warned_already_01) then
          PRINT *, 'In hartmannTemp, rM value is out of range! (',rM,')'
          PRINT *, '  assuming fieldlines 1-2'
          PRINT *, '  ==> Further warning surpressed.'
          warned_already_01 = .true.
       end if
       lowerLine = 1
       upperLine = 2 
       subPos    = 0
    ELSE IF (radius > 5.0) THEN
       if (.not. warned_already_02) then
          PRINT *, 'In hartmannTemp, rM value is out of range! (',rM,')'
          PRINT *, '  assuming fieldlines 1-2'
          PRINT *, '  ==> Further warning surpressed.'
          warned_already_02 = .true.
       end if
       lowerLine = 4
       upperLine = 5 
    ELSE 
       lowerLine = INT(radius) + 1
       upperLine = lowerLine + 1 
!       subPos = radius - REAL(lowerline)
       subPos = int(radius) - 1
    END IF
      
    IF (inTheta > pi) THEN
      theta = inTheta - real(pi)
    ELSE
      theta = inTheta
    END IF
    
    IF (inTheta > pi/2.0) theta = real(pi - theta)
    
    ! find the index position in each field line 
    do i = 1, 5 
       CALL locate(angles(i,1:nSamples(i)), &
            nSamples(i), ABS(theta), indx_T(i))
       ! Forcing the values to be within the range...
       ! ... quick fix to avoid out of range problem... (R. Kurosawa)
       if (indx_T(i) <=0) indx_T(i)=1
       if (indx_T(i) >=nSamples(i)) indx_T(i)=nSamples(i)-1

       ! get the cooresponding temperature from each lines.       
       T(i) = temperatures(i, indx_T(i))
       R(i) = REAL(i)
    end do
                
    ! interpolate between the fieldLines 
    ! polynomial interpolation 
    ! using a routine in utils_mod.f90
    ! The output here is hartmannTemp2
    R_ip = radius
    if (R_ip<=1.0) R_ip=1.0
    if (R_ip>=5.0) R_ip=5.0
    call polint(R, T, 5, R_ip, hartmannTemp2, dT)

    
  END FUNCTION hartmannTemp2


  SUBROUTINE hartmannLines(fieldline,RAD,phi,grid,point,azVec,polVec,ok)
    ! the paper: Hartman, Hewett & Calvet 1994ApJ...426..669H 
    !   quotes data from points along the magnetic field lines marked 
    !   on one of the figures. this function returns the coordinates 
    !   of points along these lines, for obtaining comparable data. 

    use inputs_mod, only: TTauriRinner, TTauriRouter, TTauriRstar, &
                               TTauriDiskHeight
                               
    INTEGER, INTENT(IN)              :: fieldline
      ! 1 is the inner fieldline, 5 is the outer fieldline
    real(oct), INTENT(IN) :: RAD 
      ! RAD is projection of radius onto z=0 plane (in units of R_star)
    REAL(double), INTENT(IN)                 :: phi    ! azimuth angle (radians)
    TYPE(gridtype), INTENT(IN)       :: grid
    TYPE(vector), INTENT(OUT)   :: point  ! coordinates returned
    TYPE(vector), INTENT(OUT)        :: azVec  ! azimuth vector
    TYPE(vector), INTENT(OUT)        :: polVec ! poloidal vector
    LOGICAL, INTENT(OUT)             :: ok     ! coordinates are valid
    
    TYPE(vector) :: starPosn
    REAL              :: rM, theta, radius
    REAL              :: swap, y

    starPosn = grid%starPos1

    ok = .FALSE.
    
    ! test if the radius is sensible
    IF ( RAD  > (TTauriRouter/TTauriRstar) ) RETURN 

    rM = TTauriRinner + REAL(fieldline-1)*((TTauriRouter - TTauriRinner) / 4.0)
    rM = rM / TTauriRstar
    
    ! test if the radius is greater than this fieldline's max radius
    IF ( RAD  > rM ) RETURN 
    
    radius = real((rM * RAD**2)**(1./3.)  )

    ! test if the point lies inside the star
    IF ( radius < 1.0 ) RETURN
    
    y = radius / rM
    theta = ASIN(SQRT(y))

    ! set up vector in x-z plane
    point = vector(RAD, 0., radius*cos(theta))
    azVec = yHat

    ! test if the point lies too close to the disk
    IF ( ABS(point%z) < (TTauriDiskHeight/TTauriRstar)) RETURN

    ! calculate the poloidal vector
    polVec = vector(3.0 * SQRT(y) * SQRT(1.0-y) / SQRT(4.0 - (3.0*y)), &
                     0.0, &
                    (2.0 - 3.0 * y) / SQRT(4.0 - 3.0 * y))

    polVec = rotateZ(polVec,-phi)

    call randomNumberGenerator(getReal=swap)
    IF ( swap > 0.5 ) THEN 
      point%z  = -1.0_oc * point%z
      polVec%z = -1.0_oc * polVec%z
    END IF

    ! rotate about z-axis
    point  = rotateZ(point,REAL(-phi,kind=oct))
    polVec = rotateZ(polVec,-phi)
    azVec  = rotateZ(azVec,-phi)
    
    CALL normalize(polVec)

    ! convert to grid units (1.e10 cm)
    point = (TTauriRstar * 1.e-10_oc) * point
    
    ! offset to match centre of star
    point = point + starPosn
    
    ok = .TRUE.

  END SUBROUTINE hartmannLines
    
  subroutine initTTauriAMR(grid)

    use constants_mod
    use vector_mod

    use inputs_mod, only: TTauriRinner, TTauriRouter, TTauriRstar, &
                               dipoleOffset, &
                               magStreamFile, thinDiskRin !, TTauriMstar
    USE magField, only: loadMagField

    implicit none
    
    type(GRIDTYPE), intent(inout)      :: grid                
    real :: theta1, theta2
   
    real :: rStar



    rStar  = TTauriRstar / 1.e10
    
    grid%rCore = rStar
    grid%rStar1 = rStar
    grid%rStar2 = 0.
    grid%dipoleOffset = dipoleOffset
!    grid%diskRadius = TTauriRInner * 1.e-10
!    grid%diskRadius = TTauriDiskRin*rStar  ! [10^10cm]
    grid%diskRadius = ThinDiskRin*rStar  ! [10^10cm]
    grid%diskNormal = VECTOR(0.,0.,1.)
    grid%diskNormal = rotateX(grid%diskNormal,dble(grid%dipoleOffSet))
    grid%starPos1 = VECTOR(0.,0.,0.) ! in units of 1.e-10 cm
    grid%starPos2 = VECTOR(9.e9,9.e9,9.e9) ! in units of 1.e-10 cm


    theta1 = asin(sqrt(TTauriRstar/TTauriRouter))
    theta2 = asin(sqrt(TTauriRstar/TTauriRinner))

    !Laccretion = (REAL(bigG,KIND=db)* &
    !              REAL(TTauriMstar,KIND=db)* &
    !              !REAL(TTauriMdot,KIND=db)/ &
    !              REAL(10e-8*Msol*secstoyears,KIND=db)/ &
    !              REAL(TTauriRstar,kind=db))* &
    !        REAL((1.0_db-(2.0_db*TTauriRstar/(TTauriRouter+TTauriRinner))),KIND=db)

    !TaccretionDouble = Laccretion / REAL(((fourPi * TTauriRstar**2)*stefanBoltz* &
    !                                  abs(cos(theta1)-cos(theta2))),kind=db)

!    sAccretion = (fourPi * TTauriRstar**2)*abs(cos(theta1)-cos(theta2))!1.e20
!    tAccretionDouble = max(1.d0, tAccretionDouble)
!    Taccretion = TaccretionDouble**0.25

!    write(*,*) (fourPi*TTauriRstar**2)*tot/lSol," solar luminosities"   
!
!    write(*,*) (fourPi*TTauriRstar**2)*tot/ &
!               (fourPi*TTauriRstar**2*stefanBoltz*4000.**4)

    ! add the accretion luminosity spectrum to the stellar spectrum,
    ! write it out and pass it to the stateq routine.
    

    IF ( grid%geometry == "magstream" ) THEN           
      CALL loadMagField(fileName=magStreamFile,starPosn=grid%starPos1,Rstar=grid%rStar1)           
    END IF

    !open(22,file="star_plus_acc.dat",form="formatted",status="unknown")
    !do i = 1, nNu
    !   fNu(i) = fNu(i) + pi*blackbody(tAccretion, 1.e8*cSpeed/ nuArray(i))* &
    !            ((1.e20*sAccretion)/(fourPi*TTauriRstar**2))
    !   write(22,*) nuArray(i), fNu(i)
    !enddo
    !tot = 0.
    !do i = 1, nnu-1
    !   tot = tot + 0.5*(nuArray(i+1)-nuArray(i))*(fnu(i+1)+fnu(i))
    !enddo
    !write(*,*) "Final star+accretion luminosity (in solar)",(fourPi*TTauriRstar**2)*tot/lSol
    !close(22)
    !newContfile="star_plus_acc.dat"


  end subroutine initTTauriAMR


  subroutine calcTestDensity(thisOctal,subcell,grid)

    use inputs_mod
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    r = real(modulus(rVec))

    thisOctal%rho(subcell) = 1.e-30
    thisOctal%temperature(subcell) = 30.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.

    if ((r > grid%rInner).and.(r < grid%rOuter)) then
       thisOctal%rho(subcell) = rho * (grid%rInner / r)**2 
       thisOctal%temperature(subcell) = 30.
       thisOctal%inFlow(subcell) = .true.
       thisOctal%etaCont(subcell) = 0.
    endif
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

! for test of diffusion zone

!    thisOctal%diffusionApprox(subcell) = .false.
!    if ((r > grid%rinner).and.(r < grid%rinner*2.)) then
!       thisOctal%temperature(subcell) = 1000.
!    endif
!    if ((r > grid%rOuter*0.5).and.(r < grid%rOuter)) then
!       thisOctal%temperature(subcell) = 100.
!    endif
!    if ((r > grid%rinner*2.).and.(r < 0.5*grid%rOuter)) then
!       thisOctal%diffusionApprox(subcell) = .true.
!    endif
    

  end subroutine calcTestDensity

  subroutine calcPathTestDensity(thisOctal,subcell)

    use inputs_mod
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    thisOctal%rho(subcell) = 1.d0

  end subroutine calcPathTestDensity


  subroutine calcToyDiscDensity(thisOctal, subcell)
    use inputs_mod, only :  hOnly, rinner, hydrodynamics
    use inputs_mod, only :  amrgridcentrez
    type(octal) :: thisOctal
    integer :: subcell
    type(vector) :: rvec
    real(double) :: ethermal, gamma
!    type(gridtype) :: grid
    
    rvec = subcellcentre(thisOctal, subcell)
    thisOctal%temperature = 1.d4
!    if(rvec%z > (amrgridsize/2.d0 - 1.d-1*amrgridsize) .and. rVec%z < &
!         (amrgridsize/2.d0 + 1.d-1*amrgridsize) .and. rVec%x > & 
!         amrgridsize/20.d0) then
!       thisOctal%rho(subcell) = 1.d5 * mHydrogen
    if(rvec%z > (amrgridcentrez-2.d5) .and. rvec%z < (amrgridcentrez+2.d5) &
!         .and. rVec%x > amrgridcentrex - 3.5d5) then
         .and. rVec%x > rinner) then
        thisOctal%rho(subcell) = 1.d6 * mHydrogen
    else
       thisOctal%rho(subcell) = 1.d2*mHydrogen
    end if
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-8
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%columnRho(subcell) = 1.d30
    thisOctal%inFlow(subcell) = .true.
    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)

    thisOctal%ionFrac(subcell,1) = 1.e-10
    thisOctal%ionFrac(subcell,2) = 1.
    if(.not. hOnly) then
       thisOctal%ionFrac(subcell,3) = 1.e-10
       thisOctal%ionFrac(subcell,4) = 1.       
    endif
    thisOctal%etaCont(subcell) = 0.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

    if (hydrodynamics) then
       gamma = 1.d0
       ethermal = 1.5d0 * (1.d0/(2.d0*mHydrogen)) * kerg * 10.d0
       thisOctal%rhoe(subcell) = thisOctal%energy(subcell) * thisOctal%rho(subcell)
       thisOctal%pressure_i(subcell) = (gamma-1.d0)* thisOctal%rho(subcell)*ethermal
       thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
!       thisOctal%boundaryCondition(subcell) = 4
    endif


  end subroutine calcToyDiscDensity

!·         I also downloaded the zip file and got an error message that read “The Compressed (zipped) Folder is invalid or corrupted.” I’m not sure if this is because of the file itself or something to do with the way I downloaded it.
!o   Pete, is this something you guys can check?

  subroutine calcLexington(thisOctal,subcell,grid)

    use inputs_mod, only : hydrodynamics
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r
    TYPE(vector) :: rVec
    real(double) :: ethermal, gamma

    ethermal  = 0.;
    gamma = 5.d0/3.d0

    rVec = subcellCentre(thisOctal,subcell)
    r = real(modulus(rVec))

    thisOctal%rho(subcell) = 1.d-30
    thisOctal%temperature(subcell) = 10000.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-8
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%inFlow(subcell) = .true.

    if (r > grid%rinner) then
       thisOctal%inFlow(subcell) = .true.
       thisOctal%rho(subcell) = 100.*mHydrogen


       thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
       thisOctal%ne(subcell) = thisOctal%nh(subcell)
       thisOctal%nhi(subcell) = 1.e-5
       thisOctal%nhii(subcell) = thisOctal%ne(subcell)
       thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)

       thisOctal%ionFrac(subcell,1) = 1.e-10
       thisOctal%ionFrac(subcell,2) = 1.
       thisOctal%ionFrac(subcell,3) = 1.e-10
       thisOctal%ionFrac(subcell,4) = 1.       
       thisOctal%etaCont(subcell) = 0.
       thisOctal%temperature(subcell) = 10000.
    endif
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

    if (hydrodynamics) then
       ethermal = 1.5d0 * (1.d0/(2.d0*mHydrogen)) * kerg * 10.d0
       thisOctal%rhoe(subcell) = thisOctal%energy(subcell) * thisOctal%rho(subcell)
       thisOctal%pressure_i(subcell) = (gamma-1.d0)* thisOctal%rho(subcell)*ethermal
       thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
       thisOctal%boundaryCondition(subcell) = 4
    endif
  end subroutine calcLexington


!Ercolano+2008/Pequignot+2001 
  subroutine calcNarrowLineRegion(thisOctal,subcell)

!    use inputs_mod, only : hydrodynamics
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell   
    TYPE(vector) :: rVec
    real(double) :: ethermal, gamma

    ethermal  = 0.
    gamma = 5.d0/3.d0

    rVec = subcellCentre(thisOctal,subcell)

    thisOctal%temperature(subcell) = 10000.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-8
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%inFlow(subcell) = .true.

    thisOctal%rho(subcell) = 1.d4*mHydrogen
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-5
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    
    thisOctal%ionFrac(subcell,1) = 1.e-10
    thisOctal%ionFrac(subcell,2) = 1.
    thisOctal%ionFrac(subcell,3) = 1.e-10
    thisOctal%ionFrac(subcell,4) = 1.       

    thisOctal%etaCont(subcell) = 0.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine calcNarrowLineRegion

  subroutine calcPointSource(thisOctal, subcell)
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    real :: r
    TYPE(vector) :: rVec

    rVec = subcellCentre(thisOctal,subcell)
    r = real(modulus(rVec))

    thisOctal%rho(subcell) = 1.d-30
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

  end subroutine calcPointSource

  subroutine calcCylinderTest(thisOctal, subcell)
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    real :: r
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    r = real(rVec%y**2 + rVec%z**2)
    r = r**0.5
    thisOctal%rho(subcell) = 1.d-40
    thisOctal%ionFrac(subcell,1) = 1.
    thisOctal%ionFrac(subcell,2) = 1.e-10
    thisOctal%temperature(subcell) = 10.
    if((abs(rVec%x) < 0.3d9) .and. r < 0.3d9) then
       !Inside cylinder
       thisOctal%rho(subcell) = 1000.d0*mHydrogen       
       thisOctal%ionFrac(subcell,1) = 1.e-10
       thisOctal%ionFrac(subcell,2) = 1.
       thisOctal%temperature(subcell) = 10000.
    end if

    thisOctal%etaCont(subcell) = 0.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-8
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%inFlow(subcell) = .true.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

  end subroutine calcCylinderTest


  subroutine calcStarburst(thisOctal,subcell)

    use inputs_mod, only : hydrodynamics
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(vector) :: rVec
    real(double) :: ethermal, gamma = 5.d0/3.d0

    rVec = subcellCentre(thisOctal,subcell)

    thisOctal%rho(subcell) = tiny(thisOctal%rho(subcell))
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-8
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%inFlow(subcell) = .true.
    thisOctal%rho(subcell) = 1.e2*mHydrogen
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-5
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    thisOctal%ionFrac(subcell,1) = 1.e-10
    thisOctal%ionFrac(subcell,2) = 1.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
    thisOctal%dustTypeFraction(subcell,1)=1.d0

    if (hydrodynamics) then
       ethermal = 1.5d0 * (1.d0/(2.d0*mHydrogen)) * kerg * 10.d0
       thisOctal%rhoe(subcell) = thisOctal%energy(subcell) * thisOctal%rho(subcell)
       thisOctal%pressure_i(subcell) = (gamma-1.d0)* thisOctal%rho(subcell)*ethermal
       thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
       thisOctal%boundaryCondition(subcell) = 4
    endif
  end subroutine calcStarburst


  subroutine calcSymbiotic(thisOctal,subcell)

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(vector) :: rVec
    real(double) :: r, v, mdot


    mDot = 1.e-8*mSol/(365.25*24.*3600.)

    rVec = subcellCentre(thisOctal, subcell)

    r = modulus(rVec - VECTOR(-250.d0*rSol/1.d10,0.d0,0.d0))
    v = 10.e5 ! 10 km/s
    thisOctal%rho(subcell) = mDot/(fourPi * r**2 * 1.d20 * v)

    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-8
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%inFlow(subcell) = .true.

    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
    thisOctal%dustTypeFraction(subcell,1)=0.d0

  end subroutine calcSymbiotic
  
  
  subroutine initWindTestAMR(grid)

    use constants_mod
    use vector_mod
    use inputs_mod

    implicit none
    
    type(GRIDTYPE), intent(inout)      :: grid                
   
    real :: rStar
    
    rStar  = real((rSol * 20.) / 1.e10)   
    grid%rCore = rStar
    grid%rStar1 = rStar
    grid%rStar2 = 0.
    grid%rInner = rStar
    grid%rOuter = rStar * 500.
    grid%starPos1 = VECTOR(0.,0.,0.) ! in units of 1.e-10 cm
    grid%starPos2 = VECTOR(9.e9,9.e9,9.e9) ! in units of 1.e-10 cm

  end subroutine initWindTestAMR


  subroutine calcProtoDensity(thisOctal,subcell,grid)

    use inputs_mod
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    r = real(modulus(rVec))

    thisOctal%rho(subcell) = 1.e-30
    thisOctal%temperature(subcell) = 1.e-3
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .false.

    if ((r > grid%rInner).and.(r < grid%rOuter)) then
       thisOctal%rho(subcell) = grid%densityScaleFac * rho * (grid%rInner / r)**rPower 
       thisOctal%temperature(subcell) = 10.
       thisOctal%inFlow(subcell) = .true.
       thisOctal%etaCont(subcell) = 0.
    endif

! put different sort of dust within 4 rinner

    if (r < (grid%rInner*4.)) then
       thisOctal%dustType(subcell) = 1
    else
       thisOctal%dustType(subcell) = 2
    endif

    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine calcProtoDensity

  subroutine calcGammaVel(thisOctal, subcell, grid)
    use inputs_mod
    type(OCTAL) :: thisOctal
    type(GRIDTYPE) :: grid
    integer :: subcell
    integer, parameter :: nsteps = 1000
    real :: ydist(nSteps), xDist(nSteps)
    real :: d1, d2, curlyR, dx, dybydx, ddash1, ddash2
    type(VECTOR) :: direction(nSteps), stagVec, shockdirection
    type(VECTOR) :: rvec,direction2
    real :: stagpoint
    real :: v, r
    integer :: i, j

    thisOctal%inflow(subcell) = .true.
    thisOctal%temperature(subcell) = teff1
    massRatio = mass1/mass2

    d1 = binarySep * (1./(massRatio+1.))
    d2 = binarySep - d1

    momRatio = (mdot1 * vterm1) / (mdot2 * vterm2)

    curlyR = sqrt(momRatio)         ! = d1/d2 (Equ 1.) 
    ! Stevens, Blondin & Pollock 1992
    ! ApJ 386 265

    stagPoint = binarySep * sqrt(momRatio) / (1. + sqrt(momRatio))

    stagVec = VECTOR(0.,0.,-d1) + VECTOR(0., 0., stagPoint)

    shockDirection = VECTOR(0.,0.,0.)


    dx = real((grid%octreeRoot%subcellSize-stagVec%z)/real(nsteps))
    yDist(1) = dx*10.
    xDist(1) = real(stagVec%z )
    direction(1) = VECTOR(0.,1.,0.)
    do j = 2,nsteps
       xDist(j) = real(stagVec%z + real(j-1)*(grid%octreeroot%subcellsize-stagVec%z)/real(nSteps-1))
       ddash1 = sqrt((xDist(j-1) + d1)**2 + yDist(j-1)**2)
       ddash2 = sqrt((xDist(j-1) - d2)**2 + yDist(j-1)**2)

       dybydx = ((curlyR*ddash2**2 + ddash1**2)*yDist(j-1)) / &
            (curlyR *ddash2**2*(xDist(j-1)+d1) + ddash1**2*(xDist(j-1)-d2))

       yDist(j) = yDist(j-1) + dx * dyBydx

       direction(j) = VECTOR(dx,dybydx,0.)
       call normalize(direction(j))
    enddo

    rVec = subcellCentre(thisOctal, subcell)
    
    if (rVec%z > xDist(1)) then
       call locate(xDist, nSteps, real(rVec%z), i)
       if (rvec%x > ydist(i)) then
          direction2 = (rVec - VECTOR(0.d0, 0.d0, -dble(d1)))
          call normalize(direction2)
          r = real(modulus(rVec - VECTOR(0.d0, 0.d0, -dble(d1))))
          if (r > rStar1) then
             v = real(vNought1 + (vterm1-vNought1)*(1.d0 - rstar1/r)**beta1)
             thisOctal%rho(subcell) = mdot1 / (fourPi * (r*1.d10)**2 * v)
             thisOctal%temperature(subcell) = real(0.9d0 * teff1)
             thisOctal%velocity(subcell) = dble(v/cspeed) * direction2
!             thisOctal%atomAbundance(subcell, 1) =  1.d-5 / mHydrogen
             thisOctal%atomAbundance(subcell, :) =  1.d0 / (4.d0*mHydrogen)
!             thisOctal%atomAbundance(subcell, 2) =  1.d0 / (4.d0*mHydrogen)
          endif
       else
          direction2 = (rVec - VECTOR(0.d0, 0.d0, dble(d2)))
          call normalize(direction2)
          r = real(modulus(rVec - VECTOR(0.d0, 0.d0, dble(d2))))
          if (r > rStar2) then
             v = real(vNought2 + (vterm2-vNought2)*(1.d0 - rstar2/r)**beta2)
             thisOctal%rho(subcell) = mdot2 / (fourPi * (r*1.d10)**2 * v)
             thisOctal%temperature(subcell) = real(0.9d0 * teff2)
             thisOctal%velocity(subcell) = dble(v/cspeed) * direction2
!             thisOctal%atomAbundance(subcell, 1) =  0.71d0 / mHydrogen
             thisOctal%atomAbundance(subcell, :) =  0.27d0 / (4.d0*mHydrogen)
!             thisOctal%atomAbundance(subcell, 2) =  0.27d0 / (4.d0*mHydrogen)
          endif
       endif
    else
       direction2 = (rVec - VECTOR(0.d0, 0.d0, -dble(d1)))
       call normalize(direction2)
       r = real(modulus(rVec - VECTOR(0.d0, 0.d0, -dble(d1))))
       if (r > rStar1) then
          v =  real(vNought1 + (vterm1-vNought1)*(1.d0 - rstar1/r)**beta1)
          thisOctal%rho(subcell) = mdot1 / (fourPi * (r*1.d10)**2 * v)
          thisOctal%temperature(subcell) = real(0.9d0 * teff1)
          thisOctal%velocity(subcell) = dble(v/cspeed) * direction2
!          thisOctal%atomAbundance(subcell, 1) =  1.d-5 / mHydrogen
          thisOctal%atomAbundance(subcell, :) =  1.d0 / (4.d0*mHydrogen)
!          thisOctal%atomAbundance(subcell, 2) =  1.d0 / (4.d0*mHydrogen)
       endif
    endif

    thisOctal%microturb(subcell) = sqrt((2.d0*kErg*thisOctal%temperature(subcell))/mhydrogen)/cspeed

  end subroutine calcGammaVel


  subroutine calcWRShellDensity(thisOctal,subcell,grid)
    use density_mod, only: density
    use inputs_mod
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r
    TYPE(vector) :: rVec
    real(double) :: v
    
    rVec = subcellCentre(thisOctal,subcell)
    r = real( modulus(rVec))

    thisOctal%rho(subcell) = 1.e-30
    thisOctal%temperature(subcell) = 30000.
    thisOctal%etaCont(subcell) = 1.e-30
    thisOctal%inFlow(subcell) = .false.
    thisOctal%velocity(subcell) = VECTOR(0.,0.,0.)

    if (((r-thisOctal%subcellSize/2.d0) > grid%rInner)) then !.and.(r < grid%rOuter)) then
       thisOctal%rho(subcell) = density(rVec, grid)
       thisOctal%temperature(subcell) = real(30000.d0 + (100000.d0-30000.d0)*(rcore/r)**4)
       thisOctal%inFlow(subcell) = .true.
       thisOctal%etaCont(subcell) = 0.
       thisOctal%ne(subcell) = thisOctal%rho(subcell)/mHydrogen
       v = 0.001d5+(vterm-0.001d5)*(1.d0 - grid%rinner/r)**beta
       thisOctal%microturb(subcell) = 20.d5/cspeed
       thisOctal%velocity(subcell) = rVec
       thisOctal%inFlow(subcell) = .true.
       call normalize(thisOctal%velocity(subcell))
       thisOctal%velocity(subcell) = thisOctal%velocity(subcell) * v/cSpeed

    endif
    CALL fillVelocityCorners(thisOctal,wrshellVelocity)

    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

!    thisOctal%atomAbundance(subcell, 1) =  1.d0 / (mHydrogen)

    thisOctal%atomAbundance(subcell, 1) =  0.9d0 * 1.d0 / (4.d0*mHydrogen)
    if (natom>1) &
    thisOctal%atomAbundance(subcell, 2) =  0.9d0 * 1.d0 / (4.d0*mHydrogen)
    

  end subroutine calcWRShellDensity

  subroutine calcHydro1DDensity(thisOctal,subcell)

    use inputs_mod
    TYPE(octal) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: xmid, x, z, r , zprime !, gd
    
    rVec = subcellCentre(thisOctal, subcell)
    x = rVec%x
    z = rVec%z

    if (thisOctal%twod) then
       r = modulus(rVec - VECTOR(0.5d0, 0.d0, 0.5d0))
    else
       r = modulus(rVec - VECTOR(0.5d0, 0.d0, 0.0d0))
    endif

    thisOctal%energy(subcell) = 2.5d0
    thisOctal%velocity(subcell) = VECTOR(0., 0., 0.)
    thisOctal%pressure_i(subcell) = 1.d0

    thisOctal%rho(subcell) = 1.d0
    xMid = -0.75d0 - z

    if  (thisOctal%threed) then
       zprime = -1.d0 - (rVec%x+rVec%y)
    else
       zprime = -rVec%x
    endif

    if (rvec%x < 0.5d0) then
       thisOctal%rho(subcell) = 1.d0
       thisOctal%energy(subcell) = 2.5d0
       thisOctal%pressure_i(subcell) = 1.d0
       thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
    else
       thisOctal%rho(subcell) = 0.125d0
       thisOctal%energy(subcell) = 2.d0
       thisOctal%pressure_i(subcell) = 0.1d0
       thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
    endif

    if (thisOctal%oneD) then
       if (rvec%x < 0.5d0) then
          thisOctal%rho(subcell) = 1.d0
          thisOctal%energy(subcell) = 2.5d0
          thisOctal%pressure_i(subcell) = 1.d0
          thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
       else
          thisOctal%rho(subcell) = 0.125d0
          thisOctal%energy(subcell) = 2.d0
          thisOctal%pressure_i(subcell) = 0.1d0
          thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
       endif

    else if(thisOctal%twoD) then
       if((rvec%z+rVec%x) <= 0.05) then
          thisOctal%rho(subcell) = 1.d0
          thisOctal%energy(subcell) = 2.5d0
          thisOctal%pressure_i(subcell) = 1.d0
          thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
       else
          thisOctal%rho(subcell) = 0.1d0
          thisOctal%energy(subcell) = 2.d0
          thisOctal%pressure_i(subcell) = 0.1d0
          thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
       end if
    else if(thisOctal%threeD) then
       if((rvec%x <0.3) .and. (rvec%y>0.2) .and. (rvec%z < -0.2)) then
          thisOctal%rho(subcell) = 1.d0
          thisOctal%energy(subcell) = 2.5d0
          thisOctal%pressure_i(subcell) = 1.d0
          thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
       else
          thisOctal%rho(subcell) = 0.125d0
          thisOctal%energy(subcell) = 2.d0
          thisOctal%pressure_i(subcell) = 0.1d0
          thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
       end if       
    endif
 
    thisOctal%phi_i(subcell) = 0.d0
    thisOctal%boundaryCondition(subcell) = 1
    thisOctal%gamma(subcell) = 7.d0/5.d0
    thisOctal%iEquationOfState(subcell) = 0

    xplusbound = 1
    xminusbound = 1

  end subroutine calcHydro1DDensity

  subroutine calcBowshock(thisOctal,subcell)

    use inputs_mod
    TYPE(octal) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: r, x, z, eKinetic
    
    rVec = subcellCentre(thisOctal, subcell)
    x = rVec%x
    z = rVec%z

    r = sqrt((x-0.75d0)**2 + z**2)

    if (r < 0.0625d0) then
       thisOctal%rho(subcell) = 100.d0
    else
       thisOctal%rho(subcell) = 0.01d0
    endif
    thisOctal%pressure_i(subcell) = 0.0015d0
    thisOctal%phi_i(subcell) = 0.d0
    thisOctal%boundaryCondition(subcell) = 1
    thisOctal%gamma(subcell) = 5.d0/3.d0
    thisOctal%iEquationOfState(subcell) = 0

    thisOctal%velocity(subcell) = VECTOR(1.d0,0.d0,0.d0)
    if ((abs(z) < 0.0625d0).and.(x > 0.75d0)) then
       thisOctal%velocity(subcell) = VECTOR(0.d0,0.d0,0.d0)
    endif
    if (r < 0.0625d0) thisOctal%velocity(subcell) = VECTOR(0.d0,0.d0,0.d0)

    thisOctal%velocity(subcell) = thisOctal%velocity(subcell)/cSpeed
    thisOctal%rhou(subcell) = thisOctal%velocity(subcell)%x*cspeed*thisOctal%rho(subcell)
    thisOctal%rhov(subcell) = thisOctal%velocity(subcell)%y*cspeed*thisOctal%rho(subcell)
    thisOctal%rhow(subcell) = thisOctal%velocity(subcell)%z*cspeed*thisOctal%rho(subcell)

    thisOctal%gamma(subcell) = 5.d0/3.d0
    thisOctal%iEquationOfState(subcell) = 0
    thisOctal%energy(subcell) = thisOctal%pressure_i(subcell) /( (thisOctal%gamma(subcell)-1.d0)*thisOctal%rho(subcell))
    eKinetic = 0.5d0 * &
         (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOCtal%rhow(subcell)**2) &
         /thisOctal%rho(subcell)**2
    thisOctal%energy(subcell) = thisOctal%energy(subcell) + eKinetic
    thisOctal%rhoe(subcell) = thisOctal%energy(subcell) * thisOctal%rho(subcell)
    inflowRho = 0.01d0
    inflowSpeed = 1.d0
    inflowMomentum = inflowSpeed * inflowRho
    inflowPressure = 0.0015d0
    inflowEnergy = (inflowPressure / ((5.d0/3.d0-1.d0)*inflowRho) + 0.5d0*inflowMomentum**2/inflowRho**2)
    inflowRhoe = inflowEnergy * inflowRho



  end subroutine calcBowshock

  subroutine calcLighthouseDensity(thisOctal,subcell)
    use inputs_mod, only : grainFrac, nDustType, cavAngle
    TYPE(octal) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: rhoAmbient, rhoShell, rShell, rMin, rMax, thetaOpening, r, theta, sigmaShell

    rhoAmbient = 2.d-20
    rhoShell = 2.d-14
    rShell = 50.d0*autocm/1.d10
    sigmaShell = 1.d0*autocm/1.d10
    rMin = 1.d0 * auToCm / 1.d10
    rMax = 300.d0 * auToCm / 1.d10
    thetaOpening = cavangle

    rVec = subcellCentre(thisOctal, subcell)
    r = modulus(rVec)

    theta = acos ( rVec%z / r)
    if ((theta >= thetaOpening).and.(r >= rMin).and.(r <= rMax)) then
       thisOctal%rho(subcell) = rhoAmbient + rhoShell * exp(-(r-rShell)**2/(2.d0*sigmaShell**2))
    else if ( (theta < thetaOpening).and.(r >= rMin).and.(r <= rMax)) then
       thisOctal%rho(subcell) = rhoAmbient
    else if ((r < rMin).or.(r > rMax)) then
       thisOctal%rho(subcell) = 1.d-30
    endif
    thisOctal%dustTypeFraction(subcell, 1:nDustType) = grainFrac(1:nDustType)
  end subroutine calcLighthouseDensity

  subroutine calcSlabDensity(thisOctal,subcell)
    use inputs_mod, only : grainFrac, nDustType, tauSlab
    TYPE(octal),target :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    real(double) rhoSlab
    real(double) :: kappa5500
    type(VECTOR) :: rVec

    kappa5500 = 2.9368e4
    rhoSlab = tauSlab / (kappa5500*3.d0*pctocm)
!    Write(*,*) "rho slab " ,rhoSlab


    rVec = subcellCentre(thisOctal, subcell)
    
    thisOctal%inflow(subcell) = .false.
    thisOctal%rho(subcell) = 1.d-30
    if ((rVec%z < -2.d0*pctocm/1.d10).and.(rVec%z > -5.d0*pctocm/1.d10)) then
       thisOctal%rho(subcell) = rhoSlab
    thisOctal%inflow(subcell) = .true.
    endif
    thisOctal%dustTypeFraction(subcell, 1:nDustType) = grainFrac(1:nDustType)


  end subroutine calcSlabDensity


  subroutine calcKelvinDensity(thisOctal,subcell)
    use inputs_mod  
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real :: u1 !, u2
    real(double) :: eKinetic


    rVec = subcellCentre(thisOctal, subcell)
    
    if (abs(rVec%z) > 0.25d0) then
       thisOctal%velocity(subcell) = VECTOR(-0.5d0, 0.d0, 0.d0)
       thisOctal%rho(subcell) = 1.d0
    else
       thisOctal%velocity(subcell) = VECTOR(0.5d0, 0.d0, 0.d0)
       thisOctal%rho(subcell) = 2.d0
    endif

    u1 = 0.d0
    if (abs(rVec%z-0.25d0) < 0.025d0) then
       u1 = real(0.025d0 * sin ( -twoPi*(rvec%x+0.5d0)/(1.d0/6.d0) ))
    else if (abs(rVec%z+0.25d0) < 0.025d0) then
       u1 = real(0.025d0 * sin (  twoPi*(rvec%x+0.5d0)/(1.d0/6.d0) ))
    endif

    thisOctal%pressure_i(subcell) = 2.5d0
    thisOctal%velocity(subcell) = (thisOctal%velocity(subcell) + VECTOR(0.d0, 0.d0, u1))/cSpeed

    thisOctal%rhou(subcell) = thisOctal%velocity(subcell)%x*cspeed*thisOctal%rho(subcell)
    thisOctal%rhov(subcell) = thisOctal%velocity(subcell)%y*cspeed*thisOctal%rho(subcell)
    thisOctal%rhow(subcell) = thisOctal%velocity(subcell)%z*cspeed*thisOctal%rho(subcell)

    thisOctal%gamma(subcell) = 5.d0/3.d0
    thisOctal%iEquationOfState(subcell) = 0
    thisOctal%energy(subcell) = thisOctal%pressure_i(subcell) /( (thisOctal%gamma(subcell)-1.d0) * thisOctal%rho(subcell))
    eKinetic = 0.5d0 * &
         (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOCtal%rhow(subcell)**2) &
         /thisOctal%rho(subcell)**2
    thisOctal%energy(subcell) = thisOctal%energy(subcell) + eKinetic
    thisOctal%rhoe(subcell) = thisOctal%energy(subcell) * thisOctal%rho(subcell)

  end subroutine Calckelvindensity

  subroutine calcRCTestDensity(thisOctal,subcell)
    use inputs_mod  
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: ekinetic


    rVec = subcellCentre(thisOctal, subcell)
    !1,4,6,7
    if(subcell == 1 .or. subcell == 4 .or. subcell == 6 .or. subcell == 7) then
       thisOctal%pressure_i(subcell) = 0.1
    else
       thisOctal%pressure_i(subcell) = 0.2
    end if

!    if (rVec%x < 0.5d0) then
!       thisOctal%pressure_i(subcell) = 0.1
!    else
!       thisOctal%pressure_i(subcell) = 0.2
!    end if
       
    thisOctal%rho(subcell) = 1.d0

    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)/cspeed

    thisOctal%rhou(subcell) = thisOctal%velocity(subcell)%x*thisOctal%rho(subcell)*cspeed
    thisOctal%rhov(subcell) = thisOctal%velocity(subcell)%y*thisOctal%rho(subcell)*cspeed
    thisOctal%rhow(subcell) = thisOctal%velocity(subcell)%z*Thisoctal%rho(subcell)*cspeed

    thisOctal%temperature(subcell) = 1.d0

    thisOctal%gamma(subcell) = 5.d0/3.d0
    thisOctal%iEquationOfState(subcell) = 0
    thisOctal%energy(subcell) = thisOctal%pressure_i(subcell) /( (thisOctal%gamma(subcell)-1.d0) * thisOctal%rho(subcell))
    eKinetic = 0.5d0 * &
         (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOCtal%rhow(subcell)**2) &
         /thisOctal%rho(subcell)**2
    thisOctal%energy(subcell) = thisOctal%energy(subcell) + eKinetic
    thisOctal%rhoe(subcell) = thisOctal%energy(subcell) * thisOctal%rho(subcell)

  end subroutine CalcRCTestDensity

  subroutine calcRTaylorDensity(thisOctal,subcell)
    use inputs_mod
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real :: u1 !, u2
    real(double) :: g, eKinetic
    real(double) :: zpos
    real(double) :: ghostSize ! - perturbation needs to be inside the ghost layers
    !Thaw- ghost stuff, only works for fixed grid
    ghostSize = 2.d0 * thisoctal%subcellsize * griddistancescale
    g = 0.1d0
    zpos = 0.0d0

    rVec = subcellCentre(thisOctal, subcell)
    
    if (rVec%z > zPos) then
       thisOctal%rho(subcell) = 2.d0
    else
       thisOctal%rho(subcell) = 1.d0
    endif
    
    u1 = 0.d0
    if (thisOctal%threed) then
       u1 = real(0.01d0*(1.d0+cos(twoPi*rVec%x))*(1.d0+cos(twoPi*rVec%y))*(1.d0+cos(twoPi*rVec%z))/8.d0)
    else
          if(rVec%x > 0.20 .and. rVec%x < 0.30 ) then
             u1 = -0.055
          end if
          if(rVec%x > 0.70 .and. rVec%x < 0.80 ) then
             u1 = -0.055
          end if
    endif
    
    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, u1)/cSpeed

    thisOctal%rhou(subcell) = thisOctal%velocity(subcell)%x*cspeed*thisOctal%rho(subcell)
    thisOctal%rhov(subcell) = thisOctal%velocity(subcell)%y*cspeed*thisOctal%rho(subcell)
    thisOctal%rhow(subcell) = thisOctal%velocity(subcell)%z*cspeed*thisOctal%rho(subcell)

    thisOctal%pressure_i(subcell) = 2.5d0 - g * thisOctal%rho(subcell) * rVec%z
    thisOctal%phi_i(subcell) = rvec%z * g 

    thisOctal%boundaryCondition(subcell) = 2

    thisOctal%gamma(subcell) = 7.d0/5.d0
    thisOctal%iEquationOfState(subcell) = 0
    thisOctal%energy(subcell) = thisOctal%pressure_i(subcell) /( (thisOctal%gamma(subcell)-1.d0) * thisOctal%rho(subcell))
    eKinetic = 0.5d0 * &
         (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOCtal%rhow(subcell)**2) &
         /thisOctal%rho(subcell)**2
    thisOctal%energy(subcell) = thisOctal%energy(subcell) + eKinetic
    thisOctal%rhoe(subcell) = thisOctal%energy(subcell) * thisOctal%rho(subcell)
    zplusbound = 1
    zminusbound = 1
    xplusbound = 2
    xminusbound = 2
    yplusbound = 1
    yminusbound = 1
  end subroutine calcRTaylorDensity


!Sets up a turbulent medium. To perform radiative round up models:
!1. Run this with all boundaries reflective or periodic
!2. Then run a radiation hydro calculation with periodic ±y, ±z and reflecting -x, freeoutnoin +x
  subroutine calcRadiativeRoundUpDensity(thisOctal,subcell)
!    use inputs_mod, only : xplusbound, xminusbound, yplusbound, yminusbound, zplusbound, zminusbound
    use inputs_mod, only : readTurb 
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    real(double) :: cs, a, b, c, eThermal, rand

    thisOctal%rho(subcell) = 300.d0*mHydrogen
    thisOctal%temperature(subcell) = 10.d0
    ethermal = (1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)    

    cs = sqrt(ethermal)   

    call randomNumberGenerator(getDouble=a)
    call randomNumberGenerator(getDouble=b)
    call randomNumberGenerator(getDouble=c)
    call randomNumberGenerator(getDouble=rand)

    if(a < 0.4 .and. a > 0.1) then
       a = cs/((a*10.d0)**2.)
       if(rand > 0.5) then
          a = - a
       end if
    else
       a = 0.d0
    end if

    if(b < 0.4 .and. b > 0.1) then
       b = cs/((b*10.d0)**2.)
       if(rand > 0.5) then
          b = - b
       end if
       else 
          b = 0.d0
    end if

    if(c < 0.4 .and. c > 0.1) then
       c = cs/((c*10.d0)**2.)
       if(rand > 0.5) then
          c = - c
       end if
      else
         c = 0.d0
    end if
      
    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)

    thisOctal%pressure_i(subcell) = thisOctal%rho(subcell)*ethermal
!The kinetic energy is handled when the grid is read
    if(readTurb) then
       thisOctal%energy(subcell) = ethermal !+ 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    else
       thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    end if
    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
    thisOctal%gamma(subcell) = 1.0
    thisOctal%iEquationOfState(subcell) = 1

    thisOctal%inFlow(subcell) = .true.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-5
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    
    thisOctal%ionFrac(subcell,1) = 1.               !HI
    thisOctal%ionFrac(subcell,2) = 1.e-10           !HII
    if (SIZE(thisOctal%ionFrac,2) > 2) then      
       thisOctal%ionFrac(subcell,3) = 1.            !HeI
       thisOctal%ionFrac(subcell,4) = 1.e-10        !HeII
    endif
    
    thisOctal%etaCont(subcell) = 0.
  end subroutine calcRadiativeRoundUpDensity

!http://www.astrosim.net/code/doku.php?id=home:codetest:hydrotest:wengen:wengen3
!Agertz et al. (2007) - fundamental differences between SPH and grid methods
  subroutine calcBlobTestDensity(thisOctal, subcell)
    use inputs_mod, only : inflowPressure, inflowRho, inflowMomentum, inflowEnergy, inflowSpeed, inflowRhoe
    use inputs_mod, only : egrad, rhograd, rhoegrad, pgrad, momgrad
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: mUnit, vUnit, lUnit, eThermal, tExt, tInt, gamma, gridSize

    gridSize = 1.d15
    tExt = 1.d7
    tInt = 1.d6
    gamma = 5.d0/3.d0

    mUnit = 2.3262d5 * mSol
    lUnit = 1.d3 * pctocm
    vUnit = 1.d5

    rVec = subcellCentre(thisOctal, subcell)

    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)

    if(sqrt(rVec%x**2.d0 + rVec%y**2.d0 + (rVec%z**2.d0)) < (197.d3*pctocm/1.d10)) then
       thisOctal%rho(subcell) = 3.13e-7*(mUnit/(lUnit**3.d0))
       thisOctal%temperature(subcell) = real(tInt)
    else
       thisOctal%rho(subcell) = 3.13e-8*(mUnit/(lUnit**3.d0))
       thisOctal%temperature(subcell) = real(tExt)
    end if

    thisOctal%gamma(subcell) = gamma
    eThermal = (kErg*thisOctal%temperature(subcell))/(gamma - 1.d0)
    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * eThermal
    thisOctal%iEquationOfState(subcell) = 1
    thisOCtal%pressure_i(subcell) = (gamma-1.d0) * eThermal * thisOctal%rho(subcell)

    inflowRho = 3.13e-8*(mUnit/(lUnit**3.d0))
    inflowSpeed = 1.d8
    inflowMomentum = inflowSpeed * inflowRho
    inflowEnergy = ((kErg*tExt)/(gamma-1.d0)) + (0.5d0*(inflowSpeed**2.d0))
    inflowRhoe = inflowEnergy * inflowRho
    inflowPressure = kErg*tExt*thisOctal%rho(subcell)

    momgrad = (inflowMomentum - (inflowMomentum/10.d0))/gridSize
    egrad = (inflowEnergy - (inflowEnergy/10.d0))/gridSize
    rhoegrad = (inflowRhoe - (inflowRhoe/10.d0))/gridSize
    pgrad = (inflowPressure - (inflowPressure/10.d0))/gridSize
    rhograd = (inflowRho - (inflowRho/10.d0))/gridSize
    
  end subroutine calcBlobTestDensity


 subroutine calcBruntDensity(thisOctal, subcell)

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: eThermal, rMod, fac, distance, rho_o
    logical, save :: firstTime = .true.
    integer, parameter :: nr = 1000
    real(double), save :: r(nr), rho(nr)
    integer :: i, j
    integer, parameter :: numClouds=4
    type(VECTOR) :: centre(numClouds)

    rho_o = 1.d3*mHydrogen
    
    rVec = subcellCentre(thisOctal, subcell)
    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)

    thisOctal%rho(subcell) = rho_o * exp(-((rVec%y+1.54d10)/3.08d10))

    centre(1) = VECTOR(0.d0, 0.d0, 0.d0)
!    centre(2) = VECTOR(-30.d0, 25.d0, 0.d0)
!    centre(3) = VECTOR(15.d0, 25.d0, 0.d0)
!    centre(4) = VECTOR(30.d0, 25.d0, 0.d0)
!    centre(5) = VECTOR(75.d0, 25.d0, 0.d0)
     
   do j = 1, numClouds
       centre(j) = (centre(j) * pctocm)/1.d10
    end do

    if (firstTime) then
       firstTime = .false.
       r = 0.d0; rho = 0.d0

       call bonnorEbertRun(10.d0, 1.d0, 1000.d0*1.d0*mhydrogen,  nr, r, rho) 
       r = r / 1.d10
       if (myrankGlobal==1) then
          do i =1 , nr
             write(55, *) r(i)*1.d10/autocm, rho(i)
          enddo
       endif
    endif

    do j = 1, numClouds
       rVec = subcellCentre(thisOctal, subcell)   
       distance = (rVec%x - centre(j)%x)**2 + (rVec%y - centre(j)%y)**2 + (rVec%z - centre(j)%z)**2 
       if((distance**0.5) < r(nr)) then
          rMod = sqrt((rVec%x - centre(j)%x)**2 + (rVec%y - centre(j)%y)**2 + (rVec%z - centre(j)%z)**2)
          call locate(r, nr, rMod, i)
          fac = (rMod-r(i))/(r(i+1)-r(i))
          thisOctal%rho(subcell) = thisOctal%rho(subcell) + (rho(i) + fac*(rho(i+1)-rho(i)))
          thisOctal%velocity(subcell) = VECTOR(0.d0, -3.d5, 0.d0) !-3kms^-1
       endif
    end do
    
    thisOctal%temperature = 10.d0
    ethermal = (1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%pressure_i(subcell) = thisOctal%rho(subcell)*ethermal
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
    thisOctal%phi_i(subcell) = -bigG * 6.d0 * mSol / (modulus(rVec)*1.d10) 
    thisOctal%gamma(subcell) = 1.0
    thisOctal%iEquationOfState(subcell) = 1

!    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
!    thisOctal%ne(subcell) = thisOctal%nh(subcell)
!    thisOctal%nhi(subcell) = 1.e-5
!    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
!    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)

!    thisOctal%ionFrac(subcell,1) = 1.               !HI 
!    Thisoctal%ionfrac(subcell,2) = 1.e-10           !HII
!    if (SIZE(thisOctal%ionFrac,2) > 2) then
!       thisOctal%ionFrac(subcell,3) = 1.            !HeI
!       thisOctal%ionFrac(subcell,4) = 1.e-10        !HeII
!    endif
!    thisOctal%etaCont(subcell) = 0.

    thisOctal%inFlow(subcell) = .true.

  end subroutine calcBruntDensity


  subroutine calcPlanarIfrontDensity(thisOctal, subcell)
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    real(double) :: eThermal

!Parameters changed to agree with Gritschender et al. 2009 iVINE paper.
    thisOctal%rho(subcell) = 7.3d-23

    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
    thisOctal%temperature(subcell) = 10.d0
    !Thaw - will probably want to change this to use returnMu
    ethermal = (1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%pressure_i(subcell) = thisOctal%rho(subcell)*ethermal
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
!    thisOctal%phi_i(subcell) = -bigG * 6.d0 * mSol / (modulus(rVec)*1.d10)
    thisOctal%gamma(subcell) = 1.0
    thisOctal%iEquationOfState(subcell) = 1
     
    thisOctal%inFlow(subcell) = .true.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-5
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    
    thisOctal%ionFrac(subcell,1) = 1.               !HI
    thisOctal%ionFrac(subcell,2) = 1.e-10           !HII
    if (SIZE(thisOctal%ionFrac,2) > 2) then      
       thisOctal%ionFrac(subcell,3) = 1.            !HeI
       thisOctal%ionFrac(subcell,4) = 1.e-10        !HeII
       
    endif
    thisOctal%etaCont(subcell) = 0.
  end subroutine calcPlanarIfrontDensity

  subroutine calcContactDiscontinuityOneDDensity(thisOctal, subcell, v1)
    use inputs_mod, only : CD_version, amrgridcentrex
    TYPE(octal) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    logical :: v1
!    logical, save :: firstTime = .true.
    real(double) :: ethermal, thisT

    rVec = subcellCentre(thisOctal, subcell)

!    thisOctal%pressure_i(subcell) = 10.d0
!
    if(CD_version == 1) then
       thisOctal%velocity(subcell) = VECTOR(0., 0., 0.)
       thisOctal%velocity(subcell)%x = thisOctal%velocity(subcell)%x/cSpeed
    else if (CD_version == 2) then
       thisOctal%velocity(subcell) = VECTOR(0.5, 0., 0.)
       thisOctal%velocity(subcell)%x = thisOctal%velocity(subcell)%x/cSpeed
    else if (CD_version == 3) then
       thisOctal%velocity(subcell) = VECTOR(2., 0., 0.)
       thisOctal%velocity(subcell)%x = thisOctal%velocity(subcell)%x/cSpeed
    else if (CD_version == 4) then
       thisOctal%velocity(subcell) = VECTOR(20., 0., 0.)
       thisOctal%velocity(subcell)%x = thisOctal%velocity(subcell)%x/cSpeed
    end if

    if(rVec%x > amrgridcentrex) then
       if(v1) then
          thisOctal%rho(subcell) = 10.d0
          thisOctal%pressure_i(subcell) = 10.d0
          thisT = 10.d0*mHydrogen/&
               (thisOctal%rho(subcell)*kerg)
       else
          thisOctal%rho(subcell) = 1000.d0
          thisOctal%pressure_i(subcell) = 1000.d0
          thisT = 1000.d0*mHydrogen/&
               (thisOctal%rho(subcell)*kerg)
       end if
    else
       if(v1) then
          thisOctal%rho(subcell) = 1.d0
          thisOctal%pressure_i(subcell) = 10.d0
          thisT = 10.d0*mHydrogen/&
               (thisOctal%rho(subcell)*kerg)
       else
          thisOctal%rho(subcell) = 1.d0
          thisOctal%pressure_i(subcell) = 1000.d0
          thisT = 1000.d0*mHydrogen/&
               (thisOctal%rho(subcell)*kerg)          
       end if

    endif

    thisOctal%pressure_i(subcell) = (thisOctal%rho(subcell)* &
         kerg*thisT)/(mHydrogen)
    
    ethermal = thisOctal%pressure_i(subcell)/thisOctal%rho(subcell)
    
    
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%rhoe(subcell) = thisOctal%energy(subcell)*thisOctal%rho(subcell)

    thisOctal%phi_i(subcell) = 0.d0

    thisOctal%iEquationOfState(subcell) = 0


  end subroutine calcContactDiscontinuityOneDDensity

  subroutine calcSB_Runaway(thisOctal, subcell)
    use inputs_mod, only : photoionPhysics
    TYPE(octal) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double), parameter :: rSphere = 0.99d0 * pctocm / 1.d10
    real(double), parameter :: rhoAmbient = 2.338d-24
    real(double) :: eThermal
    rVec = subcellCentre(thisOctal, subcell)
    if (modulus(rVec) > rSphere) then
       thisOctal%rho(subcell) = rhoAmbient
       thisOctal%temperature(subcell) = 2.e4
       if (photoionPhysics) then
          thisOctal%ionFrac(subcell,1) = 1.d-2             !HI
          thisOctal%ionFrac(subcell,2) = 1.d0           !HII
       endif
    else
       thisOctal%rho(subcell) = rhoAmbient * 2000.d0
       thisOctal%temperature(subcell) = 7.5
       if (photoionPhysics) then
          thisOctal%ionFrac(subcell,1) = 1.d0           !HI
          thisOctal%ionFrac(subcell,2) = 1.d-2               !HII
       endif
    endif

    ethermal = (1.d0/(2.d0*mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%energy(subcell) = thisOctal%rho(subcell) * ethermal
    thisOctal%rhoe(subcell) = ethermal
    thisOctal%rhou(subcell) = 0.d0
    thisOctal%rhov(subcell) = 0.d0
    thisOctal%rhow(subcell) = 0.d0
    thisOctal%gamma(subcell) = 1.0
    thisOctal%iEquationOfState(subcell) = 1
    thisOctal%inFlow(subcell) = .true.
    if (photoionPhysics) then
       thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
       thisOctal%ne(subcell) = thisOctal%nh(subcell)*thisOctal%ionFrac(subcell,2)
       thisOctal%nhi(subcell) = 1.e-5
       thisOctal%nhii(subcell) = thisOctal%ne(subcell)
       thisOctal%nHeI(subcell) = 0.d0
       thisOctal%dustTypeFraction(subcell,:) = 1.d-20
    endif
    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)

  end subroutine calcSB_Runaway

  subroutine calcContactDiscontinuityTwoDDensity(thisOctal, subcell, v1)
    use inputs_mod, only : CD_version
    TYPE(octal) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    logical :: v1
    real(double) :: x1, y1, x2, y2, x3, y3, x4, y4
    real(double) :: x, y, ycheck14, ycheck34, ycheck23, ycheck21
    real(double) :: m14, m34, m23, m21
    logical :: ok
    
    rVec = subcellCentre(thisOctal, subcell)

    thisOctal%energy(subcell) = 1.d0
    thisOctal%pressure_i(subcell) = 10.d0

    if(CD_version == 1) then
       thisOctal%velocity(subcell) = VECTOR(0., 0., 0.)
       thisOctal%velocity(subcell)%x = thisOctal%velocity(subcell)%x/cSpeed
    else if (CD_version == 2) then
       thisOctal%velocity(subcell) = VECTOR(0.5, 0., 0.)
       thisOctal%velocity(subcell)%x = thisOctal%velocity(subcell)%x/cSpeed
    else if (CD_version == 3) then
       thisOctal%velocity(subcell) = VECTOR(2., 0., 0.)
       thisOctal%velocity(subcell)%x = thisOctal%velocity(subcell)%x/cSpeed
    else if (CD_version == 4) then
       thisOctal%velocity(subcell) = VECTOR(20., 0., 0.)
       thisOctal%velocity(subcell)%x = thisOctal%velocity(subcell)%x/cSpeed
    end if

    x = rVec%x
    y = rVec%z

!         x3
!
!    x2       x4
!
!         x1

    !rotate the square
    !xmax
    x4 = 0.5*cos(1.0) + 0.5*sin(1.d0)
    !ymax
    y3 = 0.5*cos(1.0) + 0.5*sin(1.0)
    !xmin
    x2 = -0.5*cos(1.0) - 0.5*sin(1.0)
    !ymin
    y1 = -0.5*cos(1.0) - 0.5*sin(1.0)
    !others
    x3 = -0.5*cos(1.0) + 0.5*sin(1.0)
    x1 = 0.5*cos(1.0) - 0.5*sin(1.0)
    y2 = -0.5*cos(1.0) + 0.5*sin(1.0)
    y4 = 0.5*cos(1.0) - 0.5*sin(1.0)


    !get the new actual limit coordinates
    x4 = 1.d0 + x4
    x2 = 1.d0 + x2
    y3 = 1.d0 + y3
    y1 = 1.d0 + y1
    x3 = 1.d0 + x3
    y4 = 1.d0 + y4
    x1 = 1.d0 + x1
    y2 = 1.d0 + y2


!get the gradients
    m21 = (y1-y2)/(x1-x2)
    m23 = (y3-y2)/(x3-x2)
    m34 = (y4-y3)/(x4-x3)
    m14 = (y4-y1)/(x4-x1)


    
    ok = .false.
    if(x <= x4 .and. x >= x3) then
       ycheck14 = y1 + m14*(x-x1) 
       ycheck34 = y3 + m34*(x-x3) 

       if(y > ycheck14 .and. y < ycheck34) then
          ok = .true.
       end if
    else if (x < x3 .and. x > x1) then
       ycheck14 = y1 + m14*(x-x1)
       ycheck23 = y2 + m23*(x-x2)
       if(y > ycheck14 .and. y < ycheck23) then
          ok = .true.
       end if
    else if (x < x1 .and. x > x2) then
       ycheck21 = y2 + m21*(x-x2)
       ycheck23 = y2 + m23*(x-x2)
       if(y > ycheck21 .and. y < ycheck23) then
          ok = .true.
       end if

    else 
       ok = .false.
    end if


    if (ok) then
       if(v1) then
          thisOctal%rho(subcell) = 1.d-2

          thisOctal%rhoe(subcell) = thisOctal%rho(subcell)* &
               thisOctal%velocity(subcell)%x
       else
          thisOctal%rho(subcell) = 1.d0

          thisOctal%rhoe(subcell) = thisOctal%rho(subcell)* &
               thisOctal%velocity(subcell)%x

       end if
    else
       if(v1) then
          thisOctal%rho(subcell) = 1.d-3

          thisOctal%rhoe(subcell) = thisOctal%rho(subcell)* &
               thisOctal%velocity(subcell)%x*(10.d0)
       else
          thisOctal%rho(subcell) = 1.d-3

          thisOctal%rhoe(subcell) = thisOctal%rho(subcell)* &
               thisOctal%velocity(subcell)%x*(1000.d0)

       end if
    endif
    thisOctal%energy(subcell) = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)

    thisOctal%phi_i(subcell) = 0.d0
!    thisOctal%boundaryCondition(subcell) = 1
    thisOctal%gamma(subcell) = 1.0001
    thisOctal%iEquationOfState(subcell) = 0

  end subroutine calcContactDiscontinuityTwoDDensity


  subroutine calcBonnorEbertDensity(thisOctal,subcell)
    use inputs_mod, only : xplusbound, xminusbound, yplusbound, yminusbound, zplusbound, zminusbound
    use inputs_mod, only : pdrcalc
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: eThermal, rMod, fac, centre
    logical, save :: firstTime = .true.
    integer, parameter :: nr = 1000
    real(double), save :: r(nr), rho(nr)
    integer :: i

!Parameters changed to agree with Gritschender et al. 2009 iVINE paper.

    if (firstTime) then
       firstTime = .false.
       r = 0.d0; rho = 0.d0
       centre = 0.d0
       call bonnorEbertRun(10.d0, 1.d0, 1000.d0*1.d0*mhydrogen,  nr, r, rho)

       r = r / 1.d10
       if (myrankGlobal==1) then
          do i =1 , nr
             write(55, *) r(i)*1.d10/autocm, rho(i)
          enddo
       endif
    endif

    rVec = subcellCentre(thisOctal, subcell)
    rMod = modulus(rVec)
    if (rMod < r(nr)) then
       call locate(r, nr, rMod, i)
       fac = (rMod-r(i))/(r(i+1)-r(i))
       thisOctal%rho(subcell) = rho(i) + fac*(rho(i+1)-rho(i))
       thisOctal%temperature(subcell) = 10.d0
    else
       thisOctal%rho(subcell) = rho(nr)
       thisOctal%temperature(subcell) = 10.d0
    endif
    
    if (pdrcalc) then
       thisOctal%uvvector(subcell)%x = 100.*Draine/1.d10
       thisOctal%uvvector(subcell)%y = 0.d0
       thisOctal%uvvector(subcell)%z = 0.d0
       !             thisOctal%uvvectorPlus(subcell)%x = 0.d0
       thisOctal%uvvectorPlus(subcell)%x = 100.*Draine/1.d10
       thisOctal%uvvectorPlus(subcell)%y = 0.d0
       thisOctal%uvvectorPlus(subcell)%z = 0.d0
       thisOctal%uvvectorMinus(subcell)%x = 0.d0
       !             thisOctal%uvvectorMinus(subcell)%x = 0.d0
       thisOctal%uvvectorMinus(subcell)%y = 0.d0
       thisOctal%uvvectorMinus(subcell)%z = 0.d0
    endif


!THAW - temporary uniform density to check propagation stability
    thisOctal%rho(subcell) = rho(nr)
    thisOctal%temperature = 10.d0
    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
    !Thaw - will probably want to change this to use returnMu
    ethermal = (1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%pressure_i(subcell) = thisOctal%rho(subcell)*ethermal
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
    thisOctal%phi_i(subcell) = -bigG * 6.d0 * mSol / (modulus(rVec)*1.d10)
    thisOctal%gamma(subcell) = 1.0
    thisOctal%iEquationOfState(subcell) = 1
     
    thisOctal%inFlow(subcell) = .true.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-5
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    
    thisOctal%ionFrac(subcell,1) = 1.               !HI
    thisOctal%ionFrac(subcell,2) = 1.e-10           !HII
    if (SIZE(thisOctal%ionFrac,2) > 2) then      
       thisOctal%ionFrac(subcell,3) = 1.            !HeI
       thisOctal%ionFrac(subcell,4) = 1.e-10        !HeII
       
    endif
    thisOctal%etaCont(subcell) = 0.

    zplusbound = 2
    zminusbound = 2
    xplusbound = 4
    xminusbound = 4
    yplusbound = 2
    yminusbound = 2
  end subroutine calcBonnorEbertDensity


  subroutine calcBubbleDensity(thisOctal,subcell)
    use inputs_mod, only : amrgridcentrez, slabwidth
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: rMod!, fac, centre
!    logical, save :: firstTime = .true.
!    integer, parameter :: nr = 1000
!    real(double), save :: r(nr), rho(nr)
!    integer :: i

    rVec = subcellCentre(thisOctal, subcell)
    rMod = modulus(rVec)
!    if (rMod < (1.5d-10*pctocm)) then
!       thisOctal%rho(subcell) = 10.d0*mhydrogen
!       thisOctal%temperature(subcell) = 10.d0
!    else
!       thisOctal%rho(subcell) = 500.d0*mhydrogen
!       thisOctal%temperature(subcell) = 10.d0
!    endif

    if (rMod < (1.5d-10*pctocm)) then
       print *, "GOT ONE"
       thisOctal%rho(subcell) = 10.d0*mhydrogen
       thisOctal%temperature(subcell) = 10000.d0
    else
       if(abs(rVec%z - amrgridcentrez) < slabwidth) then
          thisOctal%rho(subcell) = 500.d0*mhydrogen
          thisOctal%temperature(subcell) = 10000.d0
       else
          thisOctal%rho(subcell) = 10.d-10*mhydrogen
          thisOctal%temperature(subcell) = 10000.d0
       end if
    endif

!THAW - temporary uniform density to check propagation stability
!    thisOctal%rho(subcell) = rho(nr)

!    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
    !Thaw - will probably want to change this to use returnMu
!    ethermal = (1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)
!    thisOctal%pressure_i(subcell) = thisOctal%rho(subcell)*ethermal
!    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
!    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
!    thisOctal%phi_i(subcell) = -bigG * 6.d0 * mSol / (modulus(rVec)*1.d10)
!    thisOctal%gamma(subcell) = 1.0
!    thisOctal%iEquationOfState(subcell) = 1
     
    thisOctal%inFlow(subcell) = .true.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-5
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    
    thisOctal%ionFrac(subcell,1) = 1.               !HI
    thisOctal%ionFrac(subcell,2) = 1.e-10           !HII
    if (SIZE(thisOctal%ionFrac,2) > 2) then      
       thisOctal%ionFrac(subcell,3) = 1.            !HeI
       thisOctal%ionFrac(subcell,4) = 1.e-10        !HeII
       
    endif
    thisOctal%etaCont(subcell) = 0.
  end subroutine calcBubbleDensity


  subroutine calcOffCentreExpansionDensity(thisOctal, subcell)
    type(octal) :: thisOctal
    integer :: subcell
    type(vector) :: rVec, cen
    real(double) :: R, ethermal

    cen = VECTOR(1.54d9 + (0.4d0*pcToCm/1.d10), 0.d0, 1.54d9)
    rVec = subcellCentre(thisOctal, subcell)

    rVec = rVec - cen

    R = modulus(rVec)

    if (R < (1.d-10*pcToCm)) then
       !inside the sphere
       thisOctal%rho(subcell) = 4.85d-21
       thisOctal%temperature(subcell) = 100.d0
    else
!       thisOctal%rho(subcell) = 1.6733000000000001e-26
       thisOCtal%rho(subcell) = 1.7d-24
       thisOctal%temperature(subcell) = 100.d0
    end if

    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)

    ethermal = (1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%pressure_i(subcell) = thisOctal%rho(subcell)*ethermal
    
    thisOctal%energy(subcell) = thisOctal%pressure_i(subcell)/thisOctal%rho(subcell) &
         + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)

    thisOctal%gamma(subcell) = 1.0
    thisOctal%iEquationOfState(subcell) = 1
     
    thisOctal%inFlow(subcell) = .true.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-5
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
!    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    
    thisOctal%ionFrac(subcell,1) = 1.               !HI
    thisOctal%ionFrac(subcell,2) = 1.e-10           !HII
    if (SIZE(thisOctal%ionFrac,2) > 2) then      
       thisOctal%ionFrac(subcell,3) = 1.            !HeI
       thisOctal%ionFrac(subcell,4) = 1.e-10        !HeII
       
    endif
    thisOctal%etaCont(subcell) = 0.


  end subroutine calcOffCentreExpansionDensity


  subroutine calcCoolingShockDensity(thisOctal, subcell)
    use inputs_mod, only : amrgridcentrex, inflowTemp
    use inputs_mod, only : inflowrho, inflowspeed, inflowmomentum
    use inputs_mod, only : inflowpressure, inflowenergy, inflowrhoe

    TYPE(octal) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: temperatureUnit
!    real(double) :: eKinetic
    
    rVec = subcellCentre(thisOctal, subcell)

    thisOctal%rho(subcell) = 1.d0
    thisOctal%pressure_i(subcell) = 1.d0
!    thisOctal%temperature(subcell) = 1.d0

!    thisOctal%rhoe(subcell) = 3.d0/2.d0

    thisOctal%gamma(subcell) = 5.d0/3.d0

    if (rvec%x < amrgridcentrex) then
       thisOctal%velocity(subcell) = VECTOR(32., 0., 0.)
    else
       thisOctal%velocity(subcell) = VECTOR(-32., 0., 0.)
    end if
    thisOctal%velocity(subcell)%x = thisOctal%velocity(subcell)%x/cSpeed
    thisOctal%rhou(subcell) = thisOctal%rho(subcell) * &
         (cSpeed * (thisOctal%velocity(subcell)%x))

!    ekinetic = thisOctal%rhou(subcell)**2/(thisOctal%rho(subcell)**2)

    thisOctal%rhoe(subcell) = (3.d0/2.d0) !+ ekinetic)

    thisOctal%energy(subcell) = thisOctal%rhoe(subcell)/thisOctal%rho(subcell) 
    
    temperatureUnit = (2.33d0*mHydrogen/(kerg))!*(3.d0/2.d0)*((5.d0/3.d0)-1.d0)))
    thisOctal%temperature(subcell) = real((thisOctal%gamma(subcell) - 1.d0) * &
         thisOctal%energy(subcell)*temperatureUnit)

    thisOctal%phi_i(subcell) = 0.d0
    thisOctal%iEquationOfState(subcell) = 1

    inflowRho = 1.d0
    inflowSpeed = thisOctal%velocity(subcell)%x
    inflowMomentum = inflowSpeed * inflowRho * cspeed
    inflowPressure = 1.d0
    inflowEnergy = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
    inflowRhoe = inflowEnergy * inflowRho
    inflowTemp = thisOctal%temperature(subcell)

  end subroutine calcCoolingShockDensity


  subroutine calcIsothermalShockDensity(thisOctal, subcell)
    use inputs_mod, only : CD_version, amrgridcentrex
    use inputs_mod, only : inflowrho, inflowspeed, inflowmomentum
    use inputs_mod, only : inflowpressure, inflowenergy, inflowrhoe, inflowTemp
    TYPE(octal) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: cs

    rVec = subcellCentre(thisOctal, subcell)

    thisOctal%rho(subcell) = 1.d0
    thisOctal%pressure_i(subcell) = 1.d0
!    thisOctal%temperature(subcell) = 1.d0

    thisOctal%rhoe(subcell) = 3.d0/2.d0

    if (rvec%x < amrgridcentrex) then
       if(CD_version == 1) then
          thisOctal%velocity(subcell) = VECTOR(4., 0., 0.)
       else if (CD_version == 2) then
          thisOctal%velocity(subcell) = VECTOR(8., 0., 0.)
       else if (CD_version == 3) then
          thisOctal%velocity(subcell) = VECTOR(16., 0., 0.)
       else if (CD_version == 4) then
          thisOctal%velocity(subcell) = VECTOR(32., 0., 0.)
       else if (CD_version == 5) then
          thisOctal%velocity(subcell) = VECTOR(64., 0., 0.)
       end if       
    else
       if(CD_version == 1) then
          thisOctal%velocity(subcell) = VECTOR(-4., 0., 0.)
       else if (CD_version == 2) then
          thisOctal%velocity(subcell) = VECTOR(-8., 0., 0.)
       else if (CD_version == 3) then
          thisOctal%velocity(subcell) = VECTOR(-16., 0., 0.)
       else if (CD_version == 4) then
          thisOctal%velocity(subcell) = VECTOR(-32., 0., 0.)
       else if (CD_version == 5) then
          thisOctal%velocity(subcell) = VECTOR(-64., 0., 0.)
       end if
    endif
    thisOctal%gamma(subcell) = 5.d0/3.d0



    thisOctal%temperature(subcell) = real(2.33d0*mHydrogen/(kerg))
!    thisOctal%temperature(subcell) = real(2.33*mHydrogen/(thisOctal%gamma(subcell) * kerg))
!    thisOctal%temperature(subcell) = real((3.d0/2.d0)*(2.33*mHydrogen)/(kerg*(5.d0/3.d0 - 1.d0)))
!    thisOctal%temperature(subcell) = real(thisOctal%subcellSize/(1024.d0*kerg))
    cs =  thisOctal%rho(subcell)
    cs =  cs/((2.33d0*mHydrogen))
    cs =  cs*kerg*thisOctal%temperature(subcell)
!    cs = sqrt(thisOctal%gamma(subcell) * cs/thisOctal%rho(subcell))
   cs = sqrt(cs/thisOctal%rho(subcell))
!    cs = 1.d0
    thisOctal%velocity(subcell)%x = thisOctal%velocity(subcell)%x/cSpeed
    thisOCtal%velocity(subcell)%x = thisOctal%velocity(subcell)%x*cs

    thisOctal%energy(subcell) = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
    thisOctal%phi_i(subcell) = 0.d0
!    thisOctal%boundaryCondition(subcell) = 1

    thisOctal%iEquationOfState(subcell) = 1

    inflowRho = 1.d0
    inflowSpeed = thisOctal%velocity(subcell)%x
    inflowMomentum = inflowSpeed * inflowRho * cspeed
    inflowPressure = 1.d0
    inflowEnergy = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
    inflowRhoe = inflowEnergy * inflowRho
    inflowTemp = thisOctal%temperature(subcell)

  end subroutine calcIsothermalShockDensity

!Pascal Tremblins starbench test
  subroutine calcMixingGasDensity(thisOctal, subcell)

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
!    type(VECTOR) :: rVec, cen
    real(double) :: eThermal!, R

!    cen = VECTOR(1.d0*pcToCm/1.d10, 1.d0*pcToCm/1.d10, 1.d0*pcToCm/1.d10)
!    rVec = subcellCentre(thisOctal, subcell)
!
!    rVec = rVec - cen
!
!    R = modulus(rVec)
!
!    if (R < (0.25d-10*pcToCm)) then
!       thisOctal%rho(subcell) = 1.d4*1.4d0*mHydrogen
!       thisOctal%temperature(subcell) = real(13.34d0)
!       thisOctal%ionFrac(subcell,1) = 1.               !HI
!       thisOctal%ionFrac(subcell,2) = 1.e-10           !HII!
!       if (SIZE(thisOctal%ionFrac,2) > 2) then      !
!          thisOctal%ionFrac(subcell,3) = 1.         !   !HeI
!          thisOctal%ionFrac(subcell,4) = 1.e-10        !HeII          
!       endif      
!    else!
!       t!hisOctal%rho(subcell) = 6.67d0*1.4d0*mHydrogen
!       thisOctal%temperature(subcell) = real(10000.d0)

!       endif
!    endif

    thisOctal%rho(subcell) = 5.d0*mhydrogen
    thisOctal%temperature(subcell) = real(10000.d0)
    thisOctal%ionFrac(subcell,1) = 1.e-10               !HI!
    thisOctal%ionFrac(subcell,2) = 1.           !HII
    if (SIZE(thisOctal%ionFrac,2) > 2) then      
       thisOctal%ionFrac(subcell,3) = 1.e-10            !HeI
       thisOctal%ionFrac(subcell,4) = 1.        !HeII          
    end if

    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)

    ethermal = (1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%pressure_i(subcell) = thisOctal%rho(subcell)*ethermal
   
    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * ethermal
    thisOctal%energy(subcell) = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)

    thisOctal%gamma(subcell) = 1.d0
    thisOctal%iEquationOfState(subcell) = 1
     
    thisOctal%inFlow(subcell) = .true.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-10
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
!    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    

    thisOctal%etaCont(subcell) = 0.
  end subroutine calcMixingGasDensity

  subroutine DtypeDensity(thisOctal, subcell)

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: eThermal, rMod
    real(double), parameter :: rInner = 0.2*(pcTocm/1.d10)

    rVec = subcellCentre(thisOctal, subcell)
    rMod = modulus(rVec)

    thisOctal%rho(subcell) = 5.21d-21
    thisOctal%temperature(subcell) = 1.d4

    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)

    ethermal = (1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%pressure_i(subcell) = (thisOctal%rho(subcell)/(mHydrogen))*kerg*thisOctal%temperature(subcell)

!    ethermal = 1.5d0*(1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
    thisOctal%phi_i(subcell) = 0.d0

    thisOctal%gamma(subcell) = 1.0
    thisOctal%iEquationOfState(subcell) = 1
     
    thisOctal%inFlow(subcell) = .true.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-5
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
!    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    
    thisOctal%ionFrac(subcell,1) = 1.               !HI
    thisOctal%ionFrac(subcell,2) = 1.e-10           !HII
    if (SIZE(thisOctal%ionFrac,2) > 2) then      
       thisOctal%ionFrac(subcell,3) = 1.            !HeI
       thisOctal%ionFrac(subcell,4) = 1.e-10        !HeII
       
    endif
    thisOctal%etaCont(subcell) = 0.

  end subroutine DtypeDensity

!for StarBench code comparison workshop
  subroutine calcWhalenNormanHIIExpansionDensity(thisOctal, subcell)

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: eThermal, rMod, rand
    real(double), parameter :: rInner = 0.2*(pcTocm/1.d10)

    rVec = subcellCentre(thisOctal, subcell)
    rMod = modulus(rVec)
    if (rMod < rInner) then
!       thisOctal%rho(subcell) = 1.d4*mHydrogen
       thisOctal%rho(subcell) = 2.338d-20
       thisOctal%temperature(subcell) = 100.d0
    else
       thisOctal%rho(subcell) = (2.338d-20) * (rMod/rInner)**(-2.0)
       thisOctal%temperature(subcell) = 100.d0
    endif

    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)

    ethermal = (1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%pressure_i(subcell) = thisOctal%rho(subcell)*ethermal
    
    call randomNumberGenerator(getDouble=rand)

    rand = (2.d0*(rand - 0.5d0))/10.d0
    thisOctal%pressure_i(subcell) = thisOctal%pressure_i(subcell)* (1. + rand)

    thisOctal%energy(subcell) = thisOctal%pressure_i(subcell)/thisOctal%rho(subcell) &
         + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)

    thisOctal%gamma(subcell) = 1.0
    thisOctal%iEquationOfState(subcell) = 1
     
    thisOctal%inFlow(subcell) = .true.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-5
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
!    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    
    thisOctal%ionFrac(subcell,1) = 1.               !HI
    thisOctal%ionFrac(subcell,2) = 1.e-10           !HII
    if (SIZE(thisOctal%ionFrac,2) > 2) then      
       thisOctal%ionFrac(subcell,3) = 1.            !HeI
       thisOctal%ionFrac(subcell,4) = 1.e-10        !HeII
       
    endif
    thisOctal%etaCont(subcell) = 0.


  end subroutine calcWhalenNormanHIIExpansionDensity

  subroutine calcIsoSphereDensity(thisOctal,subcell)
    use inputs_mod, only : amrgridcentrex, amrgridcentrey, amrgridcentrez
    use inputs_mod, only : pdrcalc, photoionequilibrium
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: rMod


    rVec = subcellCentre(thisOctal, subcell)
    rVec%x = rVec%x - amrgridcentrex
    rVec%y = rVec%y - amrgridcentrey
    rVec%z = rVec%z - amrgridcentrez
    rMod = modulus(rVec)
    if ((rMod*1.d10) < 3.5d0*pcToCm) then
!    if ((rMod*1.d10) < 20.d0*pcToCm) then
       thisOctal%rho(subcell) = 1.d3*mHydrogen
       thisOctal%temperature(subcell) = 10.d0
       thisOctal%ionFrac(subcell,1) = 1.d0               !HI
       thisOctal%ionFrac(subcell,2) = 1.d-10          !HII
       if (SIZE(thisOctal%ionFrac,2) > 2) then      
          thisOctal%ionFrac(subcell,3) = 1.d0            !HeI
          thisOctal%ionFrac(subcell,4) = 1.d-10        !HeII          
       endif       

!       thisOctal%ionFrac(subcell,1) = 1.d-10               !HI
!       thisOctal%ionFrac(subcell,2) = 1.d0          !HII
!       if (SIZE(thisOctal%ionFrac,2) > 2) then      
!          thisOctal%ionFrac(subcell,3) = 1.d-10            !HeI
!          thisOctal%ionFrac(subcell,4) = 1.d0        !HeII          
!       endif       

       if(pdrcalc .and. .not. photoionequilibrium) then
          thisOctal%uvvector(subcell)%x = 0.d0
          thisOctal%uvvector(subcell)%y = 0.d0
          thisOctal%uvvector(subcell)%z = 0.d0
          thisOctal%uvvectorPlus(subcell)%x = 0.d0
          thisOctal%uvvectorPlus(subcell)%y = 0.d0
          thisOctal%uvvectorPlus(subcell)%z = 0.d0
          thisOctal%uvvectorMinus(subcell)%x = 0.d0
          thisOctal%uvvectorMinus(subcell)%y = 0.d0
          thisOctal%uvvectorMinus(subcell)%z = 0.d0
       end if
       
    else
       thisOctal%rho(subcell) = 10.d0*mHydrogen
       thisOctal%temperature(subcell) = 1.d4
       thisOctal%ionFrac(subcell,1) = 1.e-10               !HI
       thisOctal%ionFrac(subcell,2) = 1.           !HII
       if (SIZE(thisOctal%ionFrac,2) > 2) then      
          thisOctal%ionFrac(subcell,3) = 1.e-10            !HeI
          thisOctal%ionFrac(subcell,4) = 1.        !HeII       
       endif
       
       if(pdrcalc .and. .not. photoionequilibrium) then
          if(rVec%x < (amrgridcentrex)) then! .and. &
 !              rvec%z < amrgridcentrez - 2.5d0*pctocm/1.d10)then
             thisOctal%uvvector(subcell)%x = 10.*Draine*pi/1.d10
             thisOctal%uvvector(subcell)%y = 0.d0
             thisOctal%uvvector(subcell)%z = 0.d0
!             thisOctal%uvvectorPlus(subcell)%x = 0.d0
             thisOctal%uvvectorPlus(subcell)%x = 10.*Draine*pi/1.d10
             thisOctal%uvvectorPlus(subcell)%y = 0.d0
             thisOctal%uvvectorPlus(subcell)%z = 0.d0
             thisOctal%uvvectorMinus(subcell)%x = 0.d0
!             thisOctal%uvvectorMinus(subcell)%x = 0.d0
             thisOctal%uvvectorMinus(subcell)%y = 0.d0
             thisOctal%uvvectorMinus(subcell)%z = 0.d0
!1.d0*Draine/1.d10

             
!          else if (rVec%x > (amrgridcentrex)+ 2.5d0*pcToCm/1.d10 .and. &
!               rvec%z < amrgridcentrez - 2.5d0*pctocm/1.d10)then 
!             thisOctal%uvvector(subcell)%x = -1.d0*Draine/1.d10
!!             thisOctal%uvvector(subcell)%y = 0.d0
 !            thisOctal%uvvector(subcell)%z = 1.d0*Draine/1.d10

          end if
          

       endif
    end if
!    if(pdrcalc) then
!       thisOctal%converged(subcell) = .false.
!       thisOctal%biChop
!    endif


    thisOctal%inFlow(subcell) = .true.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = thisOctal%nh(subcell) * thisOctal%ionfrac(subcell, 1)
    thisOctal%nhii(subcell) = thisOctal%nh(subcell) * thisOctal%ionfrac(subcell, 2)
!    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)

    thisOctal%etaCont(subcell) = 0.

  end subroutine calcIsoSphereDensity


  subroutine calcrv1TestDensity(thisOctal,subcell)
  !  use inputs_mod, only : amrgridcentrex, amrgridcentrey, amrgridcentrez
  !  use inputs_mod, only : pdrcalc, photoionequilibrium
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
!    type(VECTOR) :: rVec
!    real(double) :: rMod

    thisOctal%rho(subcell) = 1.d3*mHydrogen
    thisOctal%temperature(subcell) = 10.d0
    thisOctal%ionFrac(subcell,1) = 1.d0               !HI
    thisOctal%ionFrac(subcell,2) = 1.d-10          !HII
    
!    if(pdrcalc .and. .not. photoionequilibrium) then
    thisOctal%uvvector(subcell)%x = 10.*Draine/1.d10
    thisOctal%uvvector(subcell)%y = 0.d0
    thisOctal%uvvector(subcell)%z = 0.d0
    thisOctal%uvvectorPlus(subcell)%x = 10.*Draine/1.d10
    thisOctal%uvvectorPlus(subcell)%y = 0.d0
    thisOctal%uvvectorPlus(subcell)%z = 0.d0
    thisOctal%uvvectorMinus(subcell)%x = 0.d0
    thisOctal%uvvectorMinus(subcell)%y = 0.d0
    thisOctal%uvvectorMinus(subcell)%z = 0.d0
    !    endif
    
    thisOctal%inFlow(subcell) = .true.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = thisOctal%nh(subcell) * thisOctal%ionfrac(subcell, 1)
    thisOctal%nhii(subcell) = thisOctal%nh(subcell) * thisOctal%ionfrac(subcell, 2)
    !    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    
    thisOctal%etaCont(subcell) = 0.
 
  end subroutine calcrv1TestDensity


  subroutine calcrv2TestDensity(thisOctal,subcell)
 !   use inputs_mod, only : amrgridcentrex, amrgridcentrey, amrgridcentrez
 !!   use inputs_mod, only : pdrcalc, photoionequilibrium
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
!    type(VECTOR) :: rVec
!    real(double) :: rMod



    thisOctal%rho(subcell) = 1.d3*mHydrogen
    thisOctal%temperature(subcell) = 10.d0
    thisOctal%ionFrac(subcell,1) = 1.d0               !HI
    thisOctal%ionFrac(subcell,2) = 1.d-10          !HII
    
!    if(pdrcalc .and. .not. photoionequilibrium) then
    thisOctal%uvvector(subcell)%x = (1.d5)*Draine/1.d10
    thisOctal%uvvector(subcell)%y = 0.d0
    thisOctal%uvvector(subcell)%z = 0.d0
    thisOctal%uvvectorPlus(subcell)%x = (1.d5)*Draine/1.d10
    thisOctal%uvvectorPlus(subcell)%y = 0.d0
    thisOctal%uvvectorPlus(subcell)%z = 0.d0
    thisOctal%uvvectorminus(subcell)%x = (1.d5)*Draine/1.d10
    thisOctal%uvvectorminus(subcell)%y = 0.d0
    thisOctal%uvvectorminus(subcell)%z = 0.d0
    !    endif
    
    thisOctal%inFlow(subcell) = .true.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = thisOctal%nh(subcell) * thisOctal%ionfrac(subcell, 1)
    thisOctal%nhii(subcell) = thisOctal%nh(subcell) * thisOctal%ionfrac(subcell, 2)
    !    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    
    thisOctal%etaCont(subcell) = 0.
 
  end subroutine calcrv2TestDensity


  subroutine calcrv3TestDensity(thisOctal,subcell)
!    use inputs_mod, only : amrgridcentrex, amrgridcentrey, amrgridcentrez
!    use inputs_mod, only : pdrcalc, photoionequilibrium
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
!    type(VECTOR) :: rVec
!    real(double) :: rMod

    thisOctal%rho(subcell) = (1.*(10.d0**(5.5)))*mHydrogen
    thisOctal%temperature(subcell) = 10.d0
    thisOctal%ionFrac(subcell,1) = 1.d0               !HI
    thisOctal%ionFrac(subcell,2) = 1.d-10          !HII
    
!    if(pdrcalc .and. .not. photoionequilibrium) then
    thisOctal%uvvector(subcell)%x = 10.*Draine/1.d10
    thisOctal%uvvector(subcell)%y = 0.d0
    thisOctal%uvvector(subcell)%z = 0.d0
    thisOctal%uvvectorPlus(subcell)%x = (10.)*Draine/1.d10
    thisOctal%uvvectorPlus(subcell)%y = 0.d0
    thisOctal%uvvectorPlus(subcell)%z = 0.d0
    thisOctal%uvvectorminus(subcell)%x = (10.)*Draine/1.d10
    thisOctal%uvvectorminus(subcell)%y = 0.d0
    thisOctal%uvvectorminus(subcell)%z = 0.d0



    !    endif
    
    thisOctal%inFlow(subcell) = .true.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = thisOctal%nh(subcell) * thisOctal%ionfrac(subcell, 1)
    thisOctal%nhii(subcell) = thisOctal%nh(subcell) * thisOctal%ionfrac(subcell, 2)
    !    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    
    thisOctal%etaCont(subcell) = 0.
 
  end subroutine calcrv3TestDensity


  subroutine calcrv4TestDensity(thisOctal,subcell)
!    use inputs_mod, only : amrgridcentrex, amrgridcentrey, amrgridcentrez
!    use inputs_mod, only : pdrcalc, photoionequilibrium
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
!    type(VECTOR) :: rVec
!    real(double) :: rMod

!    thisOctal%rho(subcell) = 1.d5.5*mHydrogen
    thisOctal%rho(subcell) = (1.*(10.d0**(5.5)))*mHydrogen
    thisOctal%temperature(subcell) = 10.d0
    thisOctal%ionFrac(subcell,1) = 1.d0               !HI
    thisOctal%ionFrac(subcell,2) = 1.d-10          !HII
    
!    if(pdrcalc .and. .not. photoionequilibrium) then
    thisOctal%uvvector(subcell)%x = (1.d5)*Draine/1.d10
    thisOctal%uvvector(subcell)%y = 0.d0
    thisOctal%uvvector(subcell)%z = 0.d0
    thisOctal%uvvectorPlus(subcell)%x = (1.d5)*Draine/1.d10
    thisOctal%uvvectorPlus(subcell)%y = 0.d0
    thisOctal%uvvectorPlus(subcell)%z = 0.d0
    thisOctal%uvvectorminus(subcell)%x = (1.d5)*Draine/1.d10
    thisOctal%uvvectorminus(subcell)%y = 0.d0
    thisOctal%uvvectorminus(subcell)%z = 0.d0
    !    endif
    
    thisOctal%inFlow(subcell) = .true.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = thisOctal%nh(subcell) * thisOctal%ionfrac(subcell, 1)
    thisOctal%nhii(subcell) = thisOctal%nh(subcell) * thisOctal%ionfrac(subcell, 2)
    !    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    
    thisOctal%etaCont(subcell) = 0.
 
  end subroutine calcrv4TestDensity



  subroutine calcRadialClouds(thisOctal, subcell)

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: eThermal, rMod, fac, distance
    logical, save :: firstTime = .true.
    integer, parameter :: nr = 1000
    real(double), save :: r(nr), rho(nr)
    integer :: i, j
    integer, parameter :: numClouds=4
    type(VECTOR) :: centre(numClouds)
!    logical :: inSphere = .false.

!THaw - initially trying 30pc box centered at 0,0,0 and extending ± 15pc
!Star at 0,0,0
!BES's at 2.5, 5, 7.5, 10, 12.5pc

!    centre(1) = VECTOR(2.5d0, 0.d0, 0.d0)
!    centre(2) = VECTOR(0.d0, 5.d0, 0.d0)
!    centre(3) = VECTOR(0.d0, 0.d0, 7.5d0)
!    centre(4) = VECTOR(-10.d0, 0.d0, 0.d0)
!    centre(5) = VECTOR(0.d0, 0.d0, -12.5d0)  

!    centre(1) = VECTOR(7.5d8, 0.d0, 2.d9)
    centre(1) = VECTOR(1.7d9, 0.d0, 1.d9)
    centre(2) = VECTOR(2.5d9, 0.d0, 0.d9)
    centre(3) = VECTOR(1.d9, 0.d0, -1.2d9)
!    centre(4) = VECTOR(1.d9, 0.d0, 1.25d9)
!    centre(1) = VECTOR(7.5d8, 0.d0, 5.d9)

!    do j = 1, numClouds
!       centre(j) = (centre(j) * pctocm)/1.d10
 !   end do

    if (firstTime) then
       firstTime = .false.
       r = 0.d0; rho = 0.d0

       call bonnorEbertRun(10.d0, 1.d0, 1000.d0*1.d0*mhydrogen,  nr, r, rho)

       r = r / 1.d10
       if (myrankGlobal==1) then
          do i =1 , nr
             write(55, *) r(i)*1.d10/autocm, rho(i)
          enddo
       endif
    endif

      thisOctal%rho(subcell) = rho(nr)

      do j = 1, numClouds
         rVec = subcellCentre(thisOctal, subcell)

         distance = (rVec%x - centre(j)%x)**2 + (rVec%y - centre(j)%y)**2 + (rVec%z - centre(j)%z)**2
         if((distance**0.5) < r(nr)) then
            rMod = sqrt((rVec%x - centre(j)%x)**2 + (rVec%y - centre(j)%y)**2 + (rVec%z - centre(j)%z)**2)
            call locate(r, nr, rMod, i)
            fac = (rMod-r(i))/(r(i+1)-r(i))
            thisOctal%rho(subcell) = thisOctal%rho(subcell) + (rho(i) + fac*(rho(i+1)-rho(i)))
         endif
      end do

    thisOctal%temperature = 10.d0
    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
    !Thaw - will probably want to change this to use returnMu
    ethermal = (1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%pressure_i(subcell) = thisOctal%rho(subcell)*ethermal
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
    thisOctal%phi_i(subcell) = -bigG * 6.d0 * mSol / (modulus(rVec)*1.d10)
    thisOctal%gamma(subcell) = 1.0
    thisOctal%iEquationOfState(subcell) = 1
!THAW - temporary - for grav test
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-5
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)

    thisOctal%ionFrac(subcell,1) = 1.               !HI
    thisOctal%ionFrac(subcell,2) = 1.e-10           !HII
    if (SIZE(thisOctal%ionFrac,2) > 2) then
       thisOctal%ionFrac(subcell,3) = 1.            !HeI
       thisOctal%ionFrac(subcell,4) = 1.e-10        !HeII
    endif
    thisOctal%etaCont(subcell) = 0.


    thisOctal%inFlow(subcell) = .true.

  end subroutine calcRadialClouds

  subroutine bonnorEbertRun(t, mu, rho0,  nr, r, rho)
    use inputs_mod, only : zetacutoff
    use constants_mod
    implicit none
    real(double) :: t, rho0
    integer :: nr
    real(double) :: r(:), rho(:)
    real(double), allocatable :: zetaB(:), phi(:)
    real(double) ::  r0
    real(double) :: mu, soundSpeed, mass, zinner, eThermal, eGrav, dv
    integer :: i
    real(double) :: dr, drhodr, d2rhodr2
    character(len=80) :: message

    zinner = 0.001d0
    allocate(zetaB(1:nr), phi(1:nr))
    zetaB = 0.d0
    phi = 0.d0

    do i = 1, nr
      zetaB(i) = log10(zinner) + (log10(zetacutoff)-log10(zinner))*dble(i-1)/dble(nr-1)
   enddo
   zetaB(1:nr) = 10.d0**zetaB(1:nr)
   zetaB(1) = 0.d0
   soundSpeed = sqrt((kErg*t) /(mu*mHydrogen))

!Thaw - Gritschneder is @ 1.6pc
!  r0 = 0.89d0 * soundSpeed / sqrt(bigG * rho0)
   r0 = 1.6d0*pctocm


   do i = 1, nr
      r(i) = zetaB(i) * r0
   enddo
   
   drhodr  = 0.d0
   rho(1) = rho0
   d2rhodr2 = (-fourPi*bigG*rho(1)**2 * (mu*mHydrogen) / (kerg * t)) 
 
   do i = 2, nr
!   do i = 1, nr
      dr = r(i) - r(i-1)
      rho(i) = rho(i-1) + drhodr * dr
      d2rhodr2 = (-fourPi*bigG*rho(i)**2 * (mu*mHydrogen) / (kerg * t)) - (2.d0/r(i))*drhodr + (1.d0/rho(i))*drhodr**2
      drhodr = drhodr + d2rhodr2 * dr
   enddo

   mass = 0.d0
   eThermal = 0.d0
   eGrav = 0.d0
   do i = 2, nr
!   do i = 1, nr
      dv = fourPi*r(i)**2*(r(i)-r(i-1))
      mass = mass + rho(i)*dv
      eGrav = eGrav + bigG*dv*rho(i)*mass/r(i)
      eThermal = eThermal + (dv*rho(i)/(mu*mHydrogen))*kerg*t
   enddo
   write(message,'(a,f5.2)') "Outer radius of Bonnor-Ebert sphere is (in pc): ",r0/pctocm
   call writeInfo(message, TRIVIAL)
   write(message,'(a,f7.2)') "Mass contained in Bonnor-Ebert sphere is: ",mass/msol
   call writeInfo(message, TRIVIAL)
   write(message,'(a,1pe12.3)') "Gravitational p.e.  contained in Bonnor-Ebert sphere is: ",eGrav
   call writeInfo(message, TRIVIAL)
   write(message,'(a,1pe12.3)') "Thermal energy contained in Bonnor-Ebert sphere is: ",eThermal
   call writeInfo(message, TRIVIAL)
   write(message,'(a,f6.3)') "Ratio of thermal enery/grav energy: ",eThermal/eGrav
   call writeInfo(message, TRIVIAL)
 end subroutine bonnorEbertRun

  subroutine calcUniformSphere(thisOctal,subcell)

    use inputs_mod, only : sphereRadius, sphereMass, spherePosition, sphereVelocity, hydrodynamics
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: eThermal, rMod,  rhoSphere

    rVec = subcellCentre(thisOctal, subcell)
    rMod = modulus(rVec-spherePosition)
    rhoSphere = sphereMass / ((fourPi/3.d0) * sphereRadius**3 * 1.d30)

    if (rMod < sphereRadius) then
       thisOctal%rho(subcell) = rhoSphere
       thisOctal%temperature(subcell) = 10.d0
    else
       thisOctal%rho(subcell) = 1.d-2 * rhoSphere
       thisOctal%temperature(subcell) = 100.d0
    endif
    thisOctal%velocity(subcell) = sphereVelocity
    if (hydrodynamics) then
       thisOctal%iequationOfState(subcell) = 3 ! n=1 polytrope
       ethermal = 1.5d0*(1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)
       thisOctal%energy(subcell) = eThermal
       thisOctal%gamma(subcell) = 2.d0
    endif
  end subroutine calcUniformSphere


  subroutine calcUnimed(thisOctal,subcell)

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell

    thisOctal%rho(subcell) = 100.d0*mHydrogen
    thisOctal%temperature(subcell) = 10.d0

    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)

    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

    thisOctal%inFlow(subcell) = .true.

    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%nh2(subcell) = thisOctal%rho(subcell) / (2.d0*mHydrogen)
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-8
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)

  end subroutine calcUnimed



  subroutine calcSphere(thisOctal,subcell)

    use inputs_mod, only : sphereRadius, sphereMass, spherePosition
    use inputs_mod, only : beta, omega, hydrodynamics, rhoThreshold, cylindricalHydro
!    use inputs_mod, only : smallestCellSize
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec, vVec
    real(double) :: eThermal, rMod,  rhoSphere, rDash, fac

    rVec = subcellCentre(thisOctal, subcell)
    rMod = modulus(rVec-spherePosition)
    rhoSphere = sphereMass / ((fourPi/3.d0) * sphereRadius**3 * 1.d30)

    rVec = VECTOR(rVec%x, rVec%y, 0.d0)
    rDash = modulus(rVec)
    call normalize(rVec)

    vVec = rVec .cross. VECTOR(0.d0, 0.d0, 1.d0)
    call normalize(vVec)
    

    rhoSphere = sphereMass * (3.d0+beta) / (fourPi * sphereRadius**3 * 1.d30)

    if (hydrodynamics) thisOctal%phi_gas(subcell) = -bigG *sphereMass / (rMod * 1.d10)
    if (rMod < sphereRadius) then
       thisOctal%rho(subcell) = min(rhoThreshold,rhoSphere * (rMod/sphereRadius)**beta)
       thisOctal%temperature(subcell) = 20.d0 ! max(20.0,50.0*real(sqrt(smallestCellSize/rMod)))
       thisOctal%velocity(subcell) = ((rDash * 1.d10)*omega/cSpeed)*vVec
       if (cylindricalHydro) then
          if (hydrodynamics) thisOctal%rhov(subcell) = omega *  (rDash*1.d10) *(rDash*1.d10)*thisOctal%rho(subcell)
       endif
    else
       fac = 1.d-2
       if (rhosphere > 1.d-8) fac = 1.d-20
       thisOctal%rho(subcell) = fac * rhoSphere
       thisOctal%temperature(subcell) = 20.d0
       thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
    endif
    if (hydrodynamics) then
       thisOctal%iequationOfState(subcell) = 1 ! isothermal
       ethermal = 1.5d0*(1.d0/(2.33d0*mHydrogen))*kerg*thisOctal%temperature(subcell)
       thisOctal%energy(subcell) = eThermal
       thisOctal%gamma(subcell) = 2.d0
    endif

  end subroutine calcSphere

  subroutine addDiscDensity(thisOctal,subcell)
    use inputs_mod, only : height, alphaDisc, betaDisc, rho0, rinner, rOuter
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    real(double) :: r, h, fac
    type(VECTOR) :: rVec

    rVec = subcellCentre(thisOctal, subcell)
    r = sqrt(rVec%x**2 + rVec%y**2)
    h = height * (r / (100.d0*autocm/1.d10))**betaDisc
    
    
    fac = -0.5d0 * (rVec%z/h)**2
    fac = max(-50.d0,fac)
    thisOctal%rho(subcell) = 1.d-30
    if ((r > rInner).and.(r < router)) then
       thisOctal%rho(subcell) = max(1.d-30,dble(rho0) * (dble(rInner/r))**dble(alphaDisc) * exp(fac))
    endif

  end subroutine addDiscDensity

  subroutine calcEnvelope(thisOctal,subcell,checkSplit)

    use inputs_mod, only : rInner, rOuter, n2max, amr1d, amr2d, beta
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    real(double) :: fac, zHeight
    type(VECTOR) :: rVec
    integer :: i
    real(double), save :: rArray(250),zArray(250), rho(250), rMod, rho0
    real(double) :: junk,rDash,zDash
    integer, save :: nPoints = 250
    logical, save :: firstTime = .true.
    logical, optional :: checkSplit

    if (present(checkSplit)) checkSplit = .false.

    if (firstTime) then
       open(33, file="envelope.txt", status="old", form="formatted")
       do i = 1, nPoints
          read(33,*) rArray(i), zArray(i), junk, rho(i)
       enddo
       close(33)
       rArray = rArray * autocm/1.d10
       zArray = zArray * autocm/1.d10
       rho = rho * 2.33d0 * mHydrogen
       firstTime = .false.
    endif

    rVec = subcellCentre(thisOctal, subcell)
    rMod = modulus(rVec)
    zDash = abs(rVec%z)
    rDash = sqrt(rVec%x**2 + rVec%y**2)
    thisOctal%rho(subcell) = 1.d-30
    thisOctal%temperature(subcell) = 10.
!    rho0 =   massEnvelope &
!         / (((8.d0/3.d0) * pi * (rInner*1.d10)**1.5d0)*((rOuter*1.d10)**1.5d0 - (rInner*1.d10)**1.5d0))
!
    if (amr1d) then
       rho0 = n2max * 2.33d0 * mHydrogen
       if ((rMod > rInner).and.(rMod < rOuter)) then
          thisOctal%rho(subcell) = rho0 * (rmod/rInner)**beta
       endif
    endif
    if (amr2d) then
       if ((rDash > rArray(1)).and.(rDash < rArray(nPoints))) then
          call locate(rArray, nPoints, rDash, i)
          fac = (rDash - rArray(i))/(rArray(i+1)-rArray(i))
          zHeight = zArray(i) + (zArray(i+1)-zArray(i))*fac
          if (((zDash-thisOctal%subcellSize/2.d0) < zHeight).and.(rDash < rArray(npoints))) then
             if (rDash > rArray(1)) thisOctal%rho(subcell) = rho(i)+fac*(rho(i+1)-rho(i))
             if (present(checkSplit)) then
                if ((rArray(i+1)-rArray(i)) < 5.d0*thisOctal%subcellSize) then
                   checkSplit = .true.
                endif
             endif
          endif
       endif
    endif
  end subroutine calcEnvelope
    


  
  subroutine calcSpiral(thisOctal,subcell)
    use inputs_mod, only : vterm, period
    real(double) :: rMod, spiralFac, theta
    real(double) :: rSpiral, densityFac
    integer :: iSpiral
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec, spiralVec
!    period = 245.d0 * 24.d0 * 3600.d0

    rVec = subcellCentre(thisOctal, subcell)
    rMod = sqrt(rVec%x**2 + rVec%y**2)

    spiralFac = dble(vterm * period) / 1.d10


    theta = atan2(rVec%y, rVec%x)
    if (theta < 0.d0) theta = theta + twoPi

    thisOctal%inflow(subcell) = .true.
    thisOctal%rho(subcell) = 1.d-20
    do iSpiral = 0, 100
       rSpiral = (dble(iSpiral) + theta/twopi) * spiralFac
       densityFac = (spiralFac/rMod)**3
       spiralVec = rSpiral * VECTOR(cos(theta),sin(theta),0.d0)
       if (modulus(rVec-spiralVec) < 0.1d0*rSpiral) then
          thisOctal%rho(subcell) = densityFac
       endif
    enddo


  end subroutine calcSpiral

  subroutine splitSpiral(thisOctal, split, splitInAzimuth)
    use inputs_mod, only : vterm, period
    real(double) :: rMod, spiralFac, theta
    real(double) :: rSpiral, densityFac
    integer :: iSpiral
    logical :: split, splitInAzimuth
    TYPE(octal), INTENT(IN) :: thisOctal
    type(VECTOR) :: rVec, spiralVec


    split = .false.
    splitInAzimuth = .false.
    rVec = thisOctal%centre
    rMod = sqrt(rVec%x**2 + rVec%y**2)

    spiralFac = dble(vterm * period) / 1.d10


    theta = atan2(rVec%y, rVec%x)
    if (theta < 0.d0) theta = theta + twoPi


    do iSpiral = 0, 100
       rSpiral = (dble(iSpiral) + theta/twopi) * spiralFac
       densityFac = (spiralFac/rMod)**3
       spiralVec = rSpiral * VECTOR(cos(theta),sin(theta),0.d0)
       if ((modulus(rVec-spiralVec) < max(2.d0*thisOctal%subcellSize,rMod*thisOctal%dPhi)).or. &
           (modulus(rVec-spiralVec) < 0.2d0*rSpiral)) then
          if (thisOctal%subcellSize > 0.02d0 * rSpiral) split=.true.
          if (iSpiral < 10) then
             if (thisOctal%dPhi*radtodeg > 11.d0) then
                splitInAzimuth = .true.
                split = .true.
             endif
             if (iSpiral < 2) then
                if (thisOctal%dPhi*radtodeg > 2.1d0) then
                   splitInAzimuth = .true.
                   split = .true.
                endif
             endif
          endif
       endif

    enddo

    
  end subroutine splitSpiral


  subroutine calcinterptest(thisOctal,subcell)

!    use inputs_mod,only : amrgridsize
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec,vVec,cen
    real(double) :: eThermal, rMod
!    cen = VECTOR(amrgridSize/4.d0,amrgridSize/4.d0,amrgridSize/4.d0)
    cen = VECTOR(0.d0, 0.d0, 0.d0)
    rVec = subcellCentre(thisOctal, subcell)-cen
    rmod = modulus(rvec)
    if (rMod < 3.1d6) then
       thisOctal%rho(subcell) = 1.e-12 * (1.d0-(modulus(rVec)/3.1e6))**2
       thisOctal%temperature(subcell) = 10.d0
    else
       thisOctal%rho(subcell) = 1.d-16
       thisOctal%temperature(subcell) = 100.d0
    endif
    vVec = rvec
    call normalize(vVec)
    vVec = vVec .cross. VECTOR(0.d0, 0.d0, 1.d0)
    call normalize(vVec)
    vVec = ((1.d0/cspeed)*(rMod/3.1d6))*vVec
    thisOctal%velocity(subcell) = vVec
    thisOctal%iequationOfState(subcell) = 3 ! n=1 polytrope
    ethermal = 1.5d0*(1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%energy(subcell) = eThermal
    thisOctal%gamma(subcell) = 2.d0

  end subroutine calcinterptest

  subroutine calcGravtest(thisOctal,subcell)

    real(double) :: sphereRadius1, sphereMass1
    real(double) :: sphereRadius2, sphereMass2
    type(VECTOR) ::  spherePosition1,  spherePosition2

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: eThermal, rMod1,  rhoSphere1, rMod2,  rhoSphere2

    spherePosition1 = VECTOR(0.d0, 0.d0, 0.d0)
    sphereMass1 = 1.d0 * mSol
    sphereRadius1 = 0.5d0*pctocm/1.d10

    spherePosition2 = VECTOR(0.d0, 0.d0, 5.d8)
    sphereMass2 = 1.d0 * mSol
    sphereRadius2 =  0.5d0*pctocm/1.d10

    rVec = subcellCentre(thisOctal, subcell)
    rMod1 = modulus(rVec-spherePosition1)
    rhoSphere1 = sphereMass1 / ((fourPi/3.d0) * sphereRadius1**3 * 1.d30)
    rMod2 = modulus(rVec-spherePosition2)
    rhoSphere2 = sphereMass2 / ((fourPi/3.d0) * sphereRadius2**3 * 1.d30)


    thisOctal%rho(subcell) = 1.d-30
    thisOctal%temperature(subcell) = 100.d0

    if (rMod1 < sphereRadius1) then
       thisOctal%rho(subcell) = rhoSphere1
       thisOctal%temperature(subcell) = 10.d0
    endif

    if (rMod2 < sphereRadius2) then
       thisOctal%rho(subcell) = rhoSphere2
       thisOctal%temperature(subcell) = 10.d0
    endif
    thisOctal%velocity(subcell) = VECTOR(0.d0,0.d0,0.d0)
    thisOctal%iequationOfState(subcell) = 3 ! n=1 polytrope
    ethermal = 1.5d0*(1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%energy(subcell) = eThermal
    thisOctal%gamma(subcell) = 2.d0
    thisOctal%pressure_i(subcell) = kerg * thisOctal%temperature(subcell) * thisOctal%rho(subcell)/(2.33d0*mHydrogen)
    thisOctal%rhoe(subcell) = sqrt(thisOctal%pressure_i(subcell))*thisOctal%rho(subcell)

  end subroutine calcGravtest

  subroutine maclaurinSpheroid(thisOctal,subcell)

!    use inputs_mod, only : sphereRadius, sphereMass, spherePosition, sphereVelocity
    use utils_mod
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: z, ethermal
    real(double), parameter :: a1 = 1.d0, a2 = 1.d0, a3 = 0.5d0
    rVec = subcellCentre(thisOctal, subcell)


    if (sqrt(rVec%x**2 + rVec%y**2) < a1) then

       z = sqrt(a3**2 * (1.d0 - (rVec%x**2 + rVec%y**2)/a1**2))
       if (abs(rVec%z) <= z) then
          thisOctal%rho(subcell) = 1.d0
       else
          thisOctal%rho(subcell) = tiny(1.d0)
       endif
    endif
    thisOctal%temperature(subcell) = 100.d0
    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
    thisOctal%iequationOfState(subcell) = 3 ! n=1 polytrope
    ethermal = 1.5d0*(1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%energy(subcell) = eThermal
    thisOctal%gamma(subcell) = 2.d0
    if (.not.associated(thisOctal%biasLine3d)) allocate(thisOctal%biasLine3d(1:thisOctal%maxChildren))
    thisOctal%biasLine3d(subcell) = maclaurinPhi(rVec%x*1.d10, rVec%y*1.d10, rVec%z*1.d10, 1.d10, 1.d10, 0.5d10)


  end subroutine maclaurinSpheroid


  subroutine calcSedovDensity(thisOctal,subcell)

    use inputs_mod
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: gamma, ethermal
    logical :: blast

    rinner = 0.01e0

    gamma = 7.d0/5.d0
    thisOctal%gamma(subcell) = gamma
    thisOctal%iEquationOfState(subcell) = 0

    rVec = subcellCentre(thisOctal, subcell)-VECTOR(0.5d0,0.d0,0.d0)
    blast = .false.
    if (modulus(rVec) < rInner) then
       blast = .true.
    endif
    thisOctal%velocity(subcell) = VECTOR(0., 0., 0.)
    thisOctal%rho(subcell) = 1.d0

    if (blast) then
       ethermal = 1.d0
       thisOctal%energy(subcell) = eThermal/(pi*rInner**2)
    else
       eThermal = 1.d-5
       thisOctal%energy(subcell) = eThermal
    endif

    thisOctal%pressure_i(subcell) = (gamma-1.d0)*thisOctal%rho(subcell)*thisOctal%energy(subcell)
    thisOctal%rhoe(subcell) = thisOctal%energy(subcell) * thisOctal%rho(subcell)
    zplusbound = 1
    zminusbound = 1
    xplusbound = 1
    xminusbound = 1
  end subroutine calcSedovDensity

  subroutine calcProtoBinDensity(thisOctal,subcell)

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: omega,r, v, phi, ethermal
    real(double) :: inertia, beta, rCloud, mCloud, eGrav


!    mCloud = 1.d0 * msol
!    rCloud = 5.d16

!    mCloud = 100.d0 * msol
!    rCloud = 0.1d0 * pctocm

    mCloud = mSol
    rCloud = 3.2d16

    inertia = (2.d0/5.d0)*mCloud*rCloud**2
    eGrav = 3.d0/5.d0 * bigG * mCloud**2 / rCloud
    beta = 0.16d0
    thisOctal%velocity(subcell) = VECTOR(0., 0., 0.)
    rVec = subcellCentre(thisOctal,subcell)
    omega = sqrt(2.d0 * beta * eGrav / inertia)
    omega = 1.6d-12
    thisOctal%phi_i(subcell) = -bigG * mCloud / (modulus(rVec)*1.d10)
    thisOctal%velocity(subcell) = vector(0., 0., 0.)

    if (thisOctal%twoD) then
       phi = atan2(rVec%z, rVec%x)
    else
       phi = atan2(rVec%y, rVec%x)
    endif

    if (thisOctal%twoD) then
       r = modulus(rVec)*1.d10
       v = r * omega
       thisOctal%velocity(subcell) = (1./cspeed)*VECTOR(v*cos(phi+pi/2.d0), 0., v*sin(phi+pi/2.d0))
    else
       r = sqrt(rVec%x**2 + rVec%y**2)*1.d10
       v = r * omega
       thisOctal%velocity(subcell) = (1./cspeed)*VECTOR(v*cos(phi+pi/2.d0),  v*sin(phi+pi/2.d0), 0.0)
    endif

    if (modulus(rVec) < (rCloud/1.d10)) then
        thisOctal%temperature(subcell) = 10.d0
       thisOctal%rho(subcell) = mCloud / (4.d0/3.d0*pi*rCloud**3)
       thisOctal%rho(subcell) = thisOctal%rho(subcell) * (1.d0 + 0.5d0 * cos(2.d0 * phi)) !m=2 dens perturbation
       ethermal = 1.5d0 * (1.d0/(2.d0*mHydrogen)) * kerg * 10.d0
    else
       thisOctal%temperature(subcell) = 100.d0
       thisOctal%rho(subcell) = 1.d-1 * mCloud / (4.d0/3.d0*pi*rCloud**3)
       thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
       ethermal = 1.5d0 * (1.d0/(2.d0*mHydrogen)) * kerg * 10.d0
    endif

    thisOctal%temperature(subcell) = 20.d0
    thisOctal%iequationOfState(subcell) = 2 ! bonnell
    ethermal = 1.5d0*(1.d0/(2.33d0*mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%energy(subcell) = eThermal
    thisOctal%gamma(subcell) = 2.d0

!    eThermal = kerg * thisOctal%temperature(subcell)/(2.33d0*mHydrogen)
!    thisOctal%pressure_i(subcell) = kerg * thisOctal%temperature(subcell) * thisOctal%rho(subcell)/(2.33d0*mHydrogen)
!    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
!    thisOctal%boundaryCondition(subcell) = 4


  end subroutine calcProtoBinDensity

  subroutine calcKleyRingDensity(thisOctal,subcell)
    use inputs_mod, only : smallestCellsize
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec, vvec,zaxis
    real(double) :: r, v, ethermal, eKinetic, fac
    real(double) :: theta, r0
    real(double) :: rhoambient, tAmbient
    real(double) :: rhoRing, Tring

    rhoRing = 1.d-15
    rhoAmbient = 1.d-22
    tRing = 1.e-7
    tAmbient = 1.e-2
    zAxis = VECTOR(0.d0, 0.d0, 1.d0)
    rVec = subcellCentre(thisOctal, subcell)
    r = sqrt(rVec%x**2 + rVec%y**2)
    r0 = 1.25d4
    theta = atan2(sqrt(rVec%x**2+rVec%y**2),rVec%z)
    v = sqrt(bigG * msol /(r*1.d10))*sin(theta)
    vVec = rVec .cross. zAxis
    call normalize(vVec)
    vVec = (v/cSpeed) * vVec

    thisOctal%rho(subcell) = rhoambient
    thisOctal%temperature(subcell) = real(Tambient)
    thisOctal%velocity(subcell) = vVec
    if ((abs(rVec%z)) < 0.1d0*Smallestcellsize) then
       fac = exp(-((rVec%x-r0)/(0.1d0*r0))**2)
       thisOctal%rho(subcell) = rhoAmbient + (rhoRing-rhoAmbient)* fac
       thisOctal%temperature(subcell) = real(tAmbient * rhoAmbient/thisOctal%rho(subcell))
       thisOctal%rhoV(subcell) = thisOctal%rho(subcell) * v * (r * 1.d10)
    else
       thisOctal%rho(subcell) = rhoambient
       thisOctal%temperature(subcell) = real(tambient)
       thisOctal%rhoV(subcell) = 0.
    endif

    eKinetic = 0.5d0 * (cspeed*modulus(thisOctal%velocity(subcell)))**2

    eThermal = 1.5d0*kerg * thisOctal%temperature(subcell)/(2.33d0*mHydrogen)
!    thisOctal%gamma(subcell) = 5.d0/3.d0
!    eThermal = (kErg*thisOctal%temperature(subcell))/(thisOctal%gamma(subcell) - 1.d0)
    thisOctal%energy(subcell) = ethermal
    thisOctal%rhoe(subcell) = thisOctal%energy(Subcell) * thisOctal%rho(subcell)
!    thisOCtal%pressure_i(subcell) = (thisOctal%gamma(subcell)-1.d0) * eThermal * thisOctal%rho(subcell)
    thisOctal%iEquationOfState(subcell) = 1


  end subroutine calcKleyRingDensity
    
  subroutine calcnbodyDensity(thisOctal,subcell)

    use inputs_mod, only : inflowPressure, inflowRho, inflowMomentum, inflowEnergy, inflowSpeed, inflowRhoe
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: gamma, ethermal, soundSpeed
    real(double) :: rho0, r0, n

    gamma = 7.d0/5.d0
    rho0 = 1.0d0
    r0 = 0.4d0
    soundSpeed = 0.01d0
    n = 2.d0
    ethermal = 0.1d0
    rVec = subcellCentre(thisOctal, subcell)


    thisOctal%temperature(subcell) = 10.d0
    thisOCtal%rho(subcell) = 1.d-25
    thisOctal%pressure_i(subcell) = (thisOctal%rho(subcell)/(2.33d0*mHydrogen))*kerg*thisOctal%temperature(subcell)

    soundSpeed = sqrt(thisOctal%pressure_i(subcell)/thisOctal%rho(subcell))
    inflowSpeed = 3.d0*soundSpeed
    thisOctal%velocity(subcell) = VECTOR(inflowSpeed/cSpeed, 0.d0, 0.d0)

    eThermal = kerg * thisOctal%temperature(subcell)/(2.33d0*mHydrogen)
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%iEquationOfState(subcell) = 1

    inflowPressure = thisOctal%pressure_i(subcell)
    inflowRho = 1.d-25
    inflowMomentum = inflowRho * inflowSpeed
    inflowEnergy = thisOctal%energy(subcell)
    inflowRhoE = inflowEnergy * inflowRho

  end subroutine calcnbodyDensity

  subroutine calcBondiHoyleDensity(thisOctal,subcell)

    use inputs_mod, only : inflowPressure, inflowRho, inflowMomentum, inflowEnergy, inflowSpeed, inflowRhoe, &
         inflowTemp, amr3d
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: gamma, ethermal, soundSpeed
    real(double) :: rho0, r0, n, mdot, lambda, r
    logical, save :: firstTime = .true.

    gamma = 7.d0/5.d0
    rho0 = 1.0d0
    r0 = 0.4d0
    soundSpeed = 0.01d0
    n = 2.d0
    ethermal = 0.1d0
    rVec = subcellCentre(thisOctal, subcell)
    inflowTemp = 10.d0
    thisOctal%temperature(subcell) = 10.d0
    thisOCtal%rho(subcell) = 1.d-25
    thisOctal%pressure_i(subcell) = (thisOctal%rho(subcell)/(2.33d0*mHydrogen))*kerg*thisOctal%temperature(subcell)

    soundSpeed = sqrt(thisOctal%pressure_i(subcell)/thisOctal%rho(subcell))
    inflowSpeed = 3.d0*soundSpeed
    if (amr3d) then
       thisOctal%velocity(subcell) = VECTOR(inflowSpeed/cSpeed, 0.d0, 0.d0)
    else
       thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, inflowSpeed/cSpeed)
    endif
    eThermal = kerg * thisOctal%temperature(subcell)/(2.33d0*mHydrogen)
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%iEquationOfState(subcell) = 1

    inflowPressure = thisOctal%pressure_i(subcell)
    inflowRho = 1.d-25
    inflowMomentum = inflowRho * inflowSpeed
    inflowEnergy = thisOctal%energy(subcell)
    inflowRhoE = inflowEnergy * inflowRho
    if (firstTime) then
       if (writeoutput) then
          mdot = fourpi * bigG**2 * msol**2 * Inflowrho
          lambda = 1.12d0
          mdot = mdot * ((lambda**2 * soundSpeed**2 + inflowSpeed**2)/(soundspeed**2 + inflowSpeed**2)**4)**0.5d0
          write(*,*) "Expected mass accretion rate is ",(mDot/msol)*365.25d0*24.d0*3600.d0, " solar masses/year"
          r = bigG * msol / (soundSpeed**2)
          write(*,*) "Bondi radius is ",r
       endif
       firstTime = .false.
    endif

  end subroutine calcBondiHoyleDensity

  subroutine calcKrumholzDiscDensity(thisOctal,subcell)

!    use inputs_mod, only : inflowPressure, inflowRho, inflowMomentum, inflowEnergy, inflowSpeed, inflowRhoe
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: ethermal, soundSpeed, sinTheta
    real(double) :: sigma, r, v
    type(VECTOR) :: zAxis = VECTOR(0.d0, 0.d0, 1.d0), vVec

    rVec = subcellCentre(thisOctal, subcell)
    sinTheta = sqrt(1.d0-(abs(rVec%z)/modulus(rVec))**2)
    r = sqrt(rVec%x**2 + rVec%y**2)
    v = sqrt(bigG * mSol /(r*1.d10))
    vVec = rVec .cross. zAxis
    call normalize(vVec)
    vVec = (v/cSpeed) * vVec

    thisOctal%temperature(subcell) = 10.d0
    sigma = 0.1d0 * 2d5/min(r,2d5)

!    sigma = 1000.d0 * 2d5/min(r,2d5)

    if (((abs(rVec%z - thisOctal%subcellsize/2.d0) < thisOctal%subcellSize/10.d0).or. & 
         (abs(rVec%z + thisOctal%subcellsize/2.d0) < thisOctal%subcellSize/10.d0)).and.(r<2d5)) then
       thisOctal%rho(subcell) = sigma / (thisOctal%subcellSize*1.d10)
    else
       thisOctal%rho(subcell) = (sigma / (thisOctal%subcellSize*1.d10))/100.d0
    endif

    thisOctal%pressure_i(subcell) = 1.d-4

    thisOctal%temperature(subcell) =  real(thisOctal%pressure_i(subcell) / ((thisOctal%rho(subcell)/(2.33d0*mHydrogen))*kerg ))

    soundSpeed = sqrt(thisOctal%pressure_i(subcell)/thisOctal%rho(subcell))
    thisOctal%velocity(subcell) = vVec
    thisOctal%rhoV(subcell) = thisOctal%rho(subcell) * v * (r * 1.d10) * sinTheta

    eThermal = kerg * thisOctal%temperature(subcell)/(2.33d0*mHydrogen)
    thisOctal%energy(subcell) = ethermal + 0.5d0*(modulus(thisOctal%velocity(subcell)))**2
    thisOctal%iEquationOfState(subcell) = 1


  end subroutine calcKrumholzDiscDensity

  subroutine calcBondiDensity(thisOctal,subcell)
    use inputs_mod, only : gridDistanceScale, bondiCentre, sourcemass, smallestCellSize
    use utils_mod, only : bondi
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec, vVec
    real(double) :: ethermal, soundSpeed
    real(double) :: x, y, z, r, rBondi, rhoInfty, v
    real(double), parameter :: lambda = 1.12d0
    logical, save :: firstTime = .true.
    rVec = subcellCentre(thisOctal, subcell)-bondiCentre
    vVec = (-1.d0)*rVec
    call normalize(vVec)
    r = modulus(rVec)*gridDistanceScale
    soundSpeed = sqrt(kerg*10.d0/(2.33d0 * mHydrogen))
    rBondi = bigG*sourcemass(1)/ soundSpeed**2
    x = (r / rBondi)
    rhoInfty = 1.d-25
    call bondi(x, y, z)
    v = y * soundSpeed

    if (firstTime) then
       firstTime = .false.
       if (writeoutput) then
          write(*,*) "Bondi radius (cm): ",rBondi
          write(*,*) "Bondi radius/cell size: ",rBondi/(smallestCellSize*1.d10)
          write(*,*) "Mass accretion rate ",fourPi*1.d-25*rBondi**2*soundSpeed/msol*(365.25*24.*3600.)
          write(*,*) "theoretical ",fourPi*1.d-25*rBondi**2 * 1.12d0*soundSpeed/msol*(365.25*24.*3600.)
       endif
    endif
    thisOctal%temperature(subcell) = 10.d0
    thisOCtal%rho(subcell) = rhoinfty * z
    thisOctal%pressure_i(subcell) = (thisOctal%rho(subcell)/(2.33d0*mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%velocity(subcell) = (v/cSpeed) * Vvec

    eThermal = kerg * thisOctal%temperature(subcell)/(2.33d0*mHydrogen)
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%iEquationOfState(subcell) = 1

  end subroutine calcBondiDensity

  subroutine calcShuDensity(thisOctal,subcell)
    use inputs_mod, only : gridDistanceScale, mu
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec, vVec
!    real(double) :: x, r,  v
    logical, save :: firstTime = .true.
!    real(double) :: xArray(21) = (/ 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, &
!         0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0/)
!    real(double) :: alphaArray(21) = (/ 1.d30, 71.5d0, 27.8d0, 16.4d0, 11.5d0, 8.76d0, 7.09d0, 5.95d0, &
!         5.14d0, 4.52d0, 4.04d0, 3.66d0, 3.35d0, 3.08d0, 2.86d0, 2.67d0, 2.50d0, 2.35d0, 2.22d0, 2.10d0, 2.00d0 /)
!    real(double) :: vArray(21) = (/ -1.d30, -5.44d0, -3.47d0, -2.58d0, -2.05d0, -1.68d0, -1.40d0, -1.18d0, -1.01d0, &
!         -0.861d0, -0.735d0, -0.625d0, -0.528d0, -0.442d0, -0.363d0, -0.291d0, -0.225d0, -0.163d0, -0.106d0, -0.051d0, 0.000d0 /)
!    real(double) :: mArray(21) = (/ 0.975d0, 0.981d0, 0.993d0, 1.01d0, 1.03d0, 1.05d0, 1.08d0, 1.12d0, 1.16d0, 1.20d0, 1.25d0, &
!         1.30d0, 1.36d0, 1.42d0, 1.49d0, 1.56d0, 1.64d0, 1.72d0, 1.81d0, 1.9d0, 2.00d0 /)
!    real(double) :: a, t, r0, fac, alpha1, m, u, p, rho, temp, ethermal, mval,rhoArray(21), velArray(21)
    real(double) :: bigA, a, r0, r, rho, u, temp, p, ethermal
!    integer :: i, j

    bigA = 29.3d0
    a = sqrt(kerg*10.d0/(mu * mHydrogen))
    r0 = 5e16
    temp = 10.
    if (firstTime) then
       if (writeoutput) write(*,*) "Radius of sphere is ",r0
       if (writeoutput) write(*,*) "Mass accretion rate should be ",(133.d0 * a**3 / bigG)/msol / secstoYears
       if (writeoutput) write(*,*) "Sound speed ",a/1.d5
       if (writeoutput) write(*,*) "rho(r0) ",a**2 * bigA / (fourPi * bigG *r0**2)
       if (writeoutput) write(*,*) "mu ",mu
       firstTime = .false.
    endif
    rVec = subcellCentre(thisOctal, subcell)
    vVec = rVec
    call normalize(vVec)
    r = modulus(rVec)*gridDistanceScale

    rho = a**2 * bigA / (fourPi*bigG * max(r,1.d-30)**2)
    p = rho/(mu*mHydrogen)*kerg*temp
    u = 0.d0
    temp = 10.d0
    if (r > r0) then 
       rho = 0.01d0 * a**2 * bigA / (fourPi*bigG * r0**2)
       temp = 10.d0
    endif



    thisOctal%temperature(subcell) = real(temp)
    thisOCtal%rho(subcell) = rho
    thisOctal%pressure_i(subcell) =  p
    thisOctal%velocity(subcell) = (u/cSpeed) * Vvec



    eThermal = kerg * thisOctal%temperature(subcell)/(mu*mHydrogen)
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%iEquationOfState(subcell) = 1

  end subroutine calcShuDensity

  subroutine setupInflowParameters()
    use inputs_mod, only : inflowPressure, inflowRho, inflowMomentum, inflowEnergy, inflowSpeed, inflowRhoe, inflowTemp
    real(double) :: soundSpeed, pressure, temperature
    temperature = inflowTemp
    pressure = (inflowRho/(2.33d0*mHydrogen))*kerg*temperature
    soundSpeed = sqrt(pressure / inflowRho)
    inflowSpeed = 3.d0*soundSpeed

    inflowPressure = pressure
    inflowMomentum = inflowRho * inflowSpeed
    inflowEnergy = kerg * temperature/(2.33d0*mHydrogen) + 0.5d0*inflowSpeed**2
    inflowRhoE = inflowEnergy * inflowRho
  end subroutine setupInflowParameters

  subroutine calcEmpty(thisOctal,subcell)
    use inputs_mod, only : centralMass
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    thisOCtal%rho(subcell) = 100.d0*mHydrogen
    thisOctal%velocity = VECTOR(0.d0, 0.d0, 0.d0)
    rVec = VECTOR(0.d0, 0.d0, 0.d0)
    if (inSubcell(thisOctal, subcell, rVec)) then
       thisOctal%rho(subcell) = max(thisOctal%rho(subcell), centralMass/(thisOctal%subcellSize**3 * 1.d30)/8.d0)
    endif
  end subroutine calcEmpty
    
  subroutine benchmarkDisk(thisOctal,subcell)

    use inputs_mod, ONLY : rInner, rOuter, height, rho
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    real :: r, hr, rd
    real(double), parameter :: min_rho = 1.0d-35 ! minimum density
    TYPE(vector) :: rVec

    real :: rInnerGap, rOuterGap
    logical :: gap

    gap = .false.

    rInnerGap = real(2. * auToCm / 1.e10)
    rOuterGap = real(3. * auToCm / 1.e10)
    
    rVec = subcellCentre(thisOctal,subcell)
    r = real(modulus(rVec))

    thisOctal%rho(subcell) = min_rho
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.
    rd = rOuter / 2.
    r = real(sqrt(rVec%x**2 + rVec%y**2))
!    if (gap.and.((r < rInnerGap).or.(r > rOuterGap))) then
       if ((r > rInner).and.(r < rOuter)) then
          hr = height * (r/rd)**1.125
! Calculate density and check the exponential won't underflow
          if ( rVec%z/hr < 20.0 ) THEN
             thisOctal%rho(subcell) = rho * ((r / rd)**(-1.))*exp(-pi/4.*(rVec%z/hr)**2)
          endif
          thisOctal%rho(subcell) = max(thisOctal%rho(subcell), min_rho)
          thisOctal%temperature(subcell) = 100.
          thisOctal%inFlow(subcell) = .true.
          thisOctal%etaCont(subcell) = 0.
       endif
!    endif
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine benchmarkDisk

  subroutine parkerWind(thisOctal,subcell)
    use inputs_mod, only : amrgridsize, maxdepthamr
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(vector) :: rvec
    real(double) :: dx, ethermal
    
    rvec = subcellcentre(thisoctal, subcell)
    dx = amrgridsize/2.**(maxdepthamr)

    if(rvec%x < dx .and. abs(rvec%z) < dx) then
       thisOctal%rho(subcell) = 1.d4*mhydrogen
    else
       thisOctal%rho(subcell) = 1.d0*mhydrogen
    endif
    thisOctal%temperature(subcell) = 10.d0

    ethermal = (1.d0/(1.d0*mHydrogen)) * kerg * thisOctal%temperature(subcell)
    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * ethermal
    thisOctal%pressure_i(subcell) = (thisOctal%rho(subcell)/(0.5d0*mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%iequationOfState(subcell) = 1
    thisOctal%velocity(subcell) = vector(0.d0, 0.d0, 0.d0)

  end subroutine parkerWind


  subroutine simpleDisc(thisOctal,subcell)
    use inputs_mod, ONLY : rho !, rInner, rOuter, height
    use inputs_mod, ONLY : extmass !, stellarMass
    use inputs_mod, ONLY : sourcepos, sourceteff, sourceradius, sourcemass
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(vector) :: rvec
    
!    real(double) :: h, midRho, scaleHeight, r, disttostar!, cs, vkep
    real(double) ::midRho, scaleHeight, r, disttostar!, cs, vkep
    real(double) :: innerrad
    rvec = subcellcentre(thisoctal, subcell)

    r = real(modulus(rVec))
    disttostar = abs(rVec%x-sourcepos(1)%x)
    disttostar = disttostar
    if(disttostar /= 0.d0) then
       thisOctal%dust_t(subcell) = sourceTeff(1)*((sourceRadius(1)*rsol)/(2.d0*&
            abs(rVec%x-sourcepos(1)%x)*1.d10))**0.5d0
    else
       thisOctal%dust_t(subcell) = sourceTeff(1)*((sourceRadius(1)*rsol)/(2.d0*&
            abs((rVec%x+1.d-10)-sourcepos(1)%x)*1.d10))**0.5d0
    endif
    thisOctal%temperature(subcell) = real(thisOCtal%dust_t(subcell))
    thisOctal%rho(subcell) = extmass
    innerrad = 10.d0

    !    print *, "x ", rvec%x, rinner*autocm/1.e10, router*autocm/1.e10
    
    if(abs(rvec%x) > innerrad*autocm/1.e10 .and. abs(rvec%x) < 300.*autocm/1.e10) then
       !calc mid-plane density
       midRho = rho*(abs(rvec%x)/(innerrad*autocm/1.e10))**(-9.d0/4.d0)

       scaleHeight =  sqrt((kerg*thisOctal%temperature(subcell)*abs(rvec%x*1.d10)**3)/(mhydrogen*bigG*sourcemass(1)))       
       print *, "mid rho ", midrho, rvec%x, innerrad*autocm/1.e10
       print *, "H", scaleheight, thisOctal%temperature(subcell)
       thisOctal%rho(subcell) = midRho*exp(-(rvec%z)**2/(2.*(scaleHeight/1.d10)**2))
    endif
    thisOctal%rho(subcell) = max(thisOctal%rho(subcell), extmass)

    thisOctal%iequationofstate = 1


  end subroutine simpleDisc

  subroutine FontDisc(thisOctal,subcell)
    use inputs_mod, only : amrgridsize, maxdepthamr, amrgridcentrez
    use inputs_mod, only : sourcemass
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(vector) :: rvec
    real(double) :: dx, rg, alpha, cs
    real(double) :: upperbound
    
    alpha = 3.d0/2.d0

    
    
    rvec = subcellcentre(thisoctal, subcell)
    dx = amrgridsize/2.**(maxdepthamr)

    cs = sqrt(kerg*10.0/mhydrogen)

    rg = bigG*sourceMass(1)/cs**2

    upperbound = amrgridcentrez - 0.5*amrgridsize + dx
    
    
    if(abs(rvec%z) < dx) then
       thisOctal%rho(subcell) = 1.d8*mhydrogen*(rvec%x/rg)**(-alpha)
    else
       thisOctal%rho(subcell) = 1.d4*mhydrogen
    endif

    thisOctal%temperature(subcell) = 10.d0
    thisOctal%iequationOfState(subcell) = 1
    thisOctal%velocity(subcell) = vector(0.d0, 0.d0, 0.d0)

  end subroutine FontDisc


  subroutine RHDDisc(thisOctal,subcell)


    use inputs_mod, ONLY : rInner, rOuter, rho, hydrodynamics !, height
    use inputs_mod, ONLY : photoionPhysics,  extmass!, stellarMass
    use inputs_mod, ONLY : sourcepos, sourceteff, sourceradius, sourcemass

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    real(double) :: r, hr, rd
!    real(double), parameter :: min_rho = 1.0d-24 ! minimum density
!    real(double), parameter :: min_rho = 1.0d-21 ! minimum density
!    real(double), parameter :: min_rho = 1.0d-22 ! minimum density
    real(double) :: min_rho!, vphi, vkep

    real(double) :: ethermal, gamma, disttostar, rinnergap, routergap

    TYPE(vector) :: rVec

!    real :: rInnerGap, rOuterGap
    logical :: gap

    min_rho = extmass
    gap = .false.
    gamma = 1.d0

    gamma = 5.d0/3.d0
!    rInnerGap = 2. * auToCm / 1.e10

!    rOuterGap = 3. * auToCm / 1.e10
   
    rVec = subcellCentre(thisOctal,subcell)
    r = real(modulus(rVec))
    disttostar = abs(rVec%x-sourcepos(1)%x)
    disttostar = disttostar
    thisOctal%dust_t(subcell) = sourceTeff(1)*((sourceRadius(1)*rsol)/(2.d0*&
         abs(rVec%x-sourcepos(1)%x)*1.d10))**0.5d0

    thisOctal%rho(subcell) = min_rho

!    thisOctal%temperature(subcell) = 1.d4
    thisOctal%temperature(subcell) = real(thisOCtal%dust_t(subcell))
    thisOctal%temperature(subcell) = 500.0

    if(photoionPhysics) then
       thisOctal%etaCont(subcell) = 0.
       thisOctal%inFlow(subcell) = .true.
    end if
    rd = rOuter / 2.
    r = real(sqrt(rVec%x**2))! + rVec%y**2))
!    if (gap.and.((r < rInnerGap).or.(r > rOuterGap))) then
!    vphi= 0.d0
!    rinnerGap = rinner * autocm*1.e10
    rinnerGap = 4. * autocm/1.e10
!       print *, "r ", r
!       print *, "rinner ", rinner
    rInner = 3.*real(autocm)/1.e10
    routerGap  = 300.d0*autocm/1.e10
!    rinnerer = 
    if ((r > rInnerGap) .and.(r < rOuterGap)) then
!       hr = height * (r/rd)**1.125
       !       hr = height * (r/rd)**1.125
       hr = dble(((kerg*dble(thisOctal%temperature(subcell))/mhydrogen)**0.5)*((bigG*sourcemass(1)/(disttostar*1.d10)**3)**(-0.5)))
       print *, "hr is ", hr, rVec%z*1.d10, ((rVec%z*1.d10)/(2.d0*hr))**2
       ! Calculate density and check the exponential won't underflow

!       if ( rVec%z/hr < 20.0 ) THEN
       hr = dble(((kerg*dble(thisOctal%temperature(subcell))/mhydrogen)**0.5)*((bigG*sourcemass(1)/&
            (disttostar*1.d10)**3)**(-0.5)))
          thisOctal%rho(subcell) = rho*exp(-((rVec%z*1.d10)/(2.d0*hr))**2)
          thisOctal%rho(subcell) = max(thisOctal%rho(subcell), min_rho)

!          thisOctal%rho(subcell) = rho * ((r / rd)**(-1.))*exp(-pi/4.*(rVec%z/hr)**2)
!       endif

       thisOctal%rho(subcell) = max(thisOctal%rho(subcell), min_rho)
    endif
    

    !       if(thisOctal%rho(subcell) > min_rho) then
    !          vkep = sqrt(biG*stellarMass(1)/r)
    !          vphi = Vkep*(1.d0-eta*(rvec%z/r)**2)
    !       endif
       !       thisOctal%temperature(subcell) = 100.
    ! endif
    !    endif
    !    thisOctal%nh = thisOctal%rho(subcell)/
    
    thisOctal%velocity = VECTOR(0.,0.,0.)
    if(photoionphysics) then
       thisOctal%biasCont3D = 1.
       thisOctal%etaLine = 1.e-30
    end if
    if (hydrodynamics) then
!       ethermal = 1.5d0 * (1.d0/(2.d0*mHydrogen)) * kerg * 10.d0
       thisOctal%phi_gas(Subcell) = 0.d0
       thisOctal%phi_i(Subcell) = 0.d0
       ethermal = (1.d0/(1.d0*mHydrogen)) * kerg * thisOctal%temperature(subcell)
       thisOctal%rhoe(subcell) = thisOctal%energy(subcell) * thisOctal%rho(subcell)
       thisOctal%pressure_i(subcell) = thisOctal%rho(subcell)*ethermal
       thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
!       thisOctal%boundaryCondition(subcell) = 4
    endif

    thisOctal%iEquationOfState(subcell) = 1

  end subroutine RHDDisc

  subroutine molecularBenchmark(thisOctal,subcell)

    use inputs_mod, only : molAbundance !, amr2d

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    logical, save :: firsttime = .true.
    integer, parameter :: nr = 50
    real(double),save :: r(nr), nh2(nr), junk,t(nr), v(nr) , mu(nr)
    real(double) :: mu1, r1, t1, t2
    real(double) :: v1 !, vDopp
    integer :: i

    type(VECTOR) :: vel

    if (firsttime) then
       open(31, file="model_1.dat", status="old", form="formatted") ! Model 2 in the Hogerheijde 2000 paper. 
       do i = nr,1,-1                                             
          read(31,*) r(i), nh2(i), junk,t(i), v(i) , mu(i)
       enddo
       r = r * 1.e-10
       close(31)
       firsttime = .false.
    endif

    r1 = modulus(subcellCentre(thisOctal,subcell))
    thisOctal%temperature(subcell) = real(tcbr)
    thisOctal%microTurb(subcell) = 1d-8 !0.159e5/cspeed
    thisOctal%velocity(subcell) = VECTOR(-1d4,-1d4,-1d4)
    thisOctal%molAbundance(subcell) = molAbundance

    if(r1 > r(nr) .or. r1 < r(1)) then 
       thisOctal%nh2(subcell) = 1.e-20
       thisOctal%rho(subcell) = 1.e-20 * 2. * mhydrogen
    else

       call locate(r, nr, r1, i)
       t2 = (r1 - r(i))/(r(i+1)-r(i)) ! linear but know its a power law so use better interpolation

       t1 = log(r1/r(i))/log(r(i+1)/r(i))

       thisOctal%nh2(subcell) = exp((1.d0 - t1) * log(nh2(i))  +  t1 * log(nh2(i+1)))
       
       thisOctal%rho(subcell) = thisOctal%nh2(subcell)*2.d0*mhydrogen

!       thisOctal%temperature(subcell) = t(i) + t1 * (t(i+1)-t(i))
!       thisOctal%nh2(subcell) = nh2(i) + t1 * (nh2(i+1)-nh2(i))
!       thisOctal%rho(subcell) = (nh2(i) + t1 * (nh2(i+1)-nh2(i)))*2.*mhydrogen

       thisOctal%temperature(subcell) = real((t(i) + t2 * (t(i+1)-t(i))))
!       thisOctal%temperature(subcell) = 15. ! debug 30/8/08
!       thisOctal%nh2(subcell) = nh2(i) + t1 * (nh2(i+1)-nh2(i))
!       thisOctal%rho(subcell) = (nh2(i) + t1 * (nh2(i+1)-nh2(i)))*2.*mhydrogen


!       if(amr2d) then

!          thisOctal%nh2(subcell) = thisOctal%nh2(subcell) * costheta**2
!          thisOctal%rho(subcell) = thisOctal%rho(subcell) * costheta**2

!       endif

       v1 = (v(i) + t2 * (v(i+1)-v(i)))*1.d5
       vel = subcellCentre(thisOctal, subcell)
       call normalize(vel)
       thisOctal%velocity(subcell) = (v1 * vel)/cspeed
!       thisOctal%velocity(subcell) = vector(0.d0,0.d0,0.d0) ! debug 17/1/08
       mu1 = mu(i) + t2 * (mu(i+1)-mu(i))

      thisOctal%microturb(subcell) = max(1d-8,sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) / (29.0 * amu)) + mu1**2) &
                                     / (cspeed * 1e-5)) ! mu is subsonic turbulence
   endif
  end subroutine molecularBenchmark

  subroutine calcTriangle(thisOctal,subcell)

    use inputs_mod, only : molAbundance, tKinetic, vturb, n2max, amrGridSize

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    real(double) :: x, bigL
    
    bigL = amrGridSize


    x = modulus(subcellCentre(thisOctal,subcell))
    thisOctal%temperature(subcell) = real(tkinetic,si)
    thisOctal%microTurb(subcell) = vturb*1.d5/cspeed
    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
    thisOctal%molAbundance(subcell) = molAbundance

    if (n2max < 0.d0) then
       thisOctal%nh2(subcell) = abs(n2max)
    else

       if (x < bigL/2.d0) then
          thisOctal%nh2(subcell) = 2.d0 * n2max * x / bigL
       else
          thisOctal%nh2(subcell) = 2.d0 * n2max * (bigL - x)/bigL
       endif
    endif
    thisOctal%rho(subcell) = thisOctal%nh2(Subcell) * 2.d0 * mHydrogen

  end subroutine calcTriangle


  subroutine calcArbitrary(thisOctal,subcell)
    use unix_mod, only: unixGetenv
    use inputs_mod, only : molAbundance, tKinetic, vturb, amrGridSize
    use inputs_mod, only : rhofile, nrholines

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    real(double) :: x, bigL
    integer :: ier, i
    real(double) :: xpos(nrholines), rho(nrholines)
    real(double) :: dx, grad
    logical :: found
    character(len=120) :: dataDirectory, file, thisrhofile

!not efficient, but minimal code modification
!    filename = "./", rhofile
!    rhofile = trim(rhofile)
    call unixGetenv("WORKDIR", dataDirectory, i)
    file = (trim(dataDirectory)//"/"//rhofile)
!    print *, "ATTEMPTING TO OPEN FILE"
    thisrhofile = trim(file)
!   print *, thisrhofile

    open(20, file=thisrhofile, status="old",  form="formatted", position="rewind", iostat=ier)
!    open(unit=20, iostat=ier, file=rhofile, form="unformatted", status="old")
    print *, "iostat = ", ier
    print *, "nrholines", nrholines
    do i = 1, nrholines
!       print *, "i ", i 
       read(20, *, iostat=ier) xpos(i), rho(i)
       if (ier /= 0) then
          call torus_abort("read failure in arbitrary geometry")
       end if
!       print *, "x ", xpos(i)
!       print *, "rho ", rho(i)
!       print *, "nlines ", nrholines
!       print *, "rhofile ", rhofile
    end do

    close(20)


    bigL = amrGridSize

    x = modulus(subcellCentre(thisOctal,subcell))
    found = .false.
!inefficient but itll do for now
    i = 1
    do while (i <  nrholines .and. .not. found)
!       print *, "x ", x
!       print *, "thisx ", xpos(i)
       

       if (x >= xpos(i) .and. x <= xpos(i+1)) then
          dx = x - xpos(i)
          grad = (rho(i+1) - rho(i))/(xpos(i+1) - xpos(i))
          thisOctal%nh2(subcell) = rho(i) + (grad*dx)
          found = .true.
!          exit
       end if  

       if(x < xpos(1)) then
          dx = xpos(1) - x
          grad = (rho(2) - rho(1))/(xpos(2) - xpos(1))
          thisOctal%nh2(subcell) = rho(1) - (grad*dx)
          found = .true.
       end if

       if( x > xpos(nrholines)) then
          dx = x - xpos(nrholines)
          grad = (rho(nrholines) - rho(nrholines - 1))/(xpos(nrholines) - xpos(nrholines - 1))
          thisOctal%nh2(subcell) = rho(nrholines) + (grad*dx)
          found = .true.
       end if
       i = i + 1
    end do

    if(.not. found) then
       print *, "x", x
       print *, "xfilemin", xpos(1)
       print *, "xfilemax", xpos(nrholines)
       call torus_abort("Error setting up grid from inputrho file, tell tom h off")
    end if

    thisOctal%temperature(subcell) = real(tkinetic,si)
    thisOctal%microTurb(subcell) = vturb*1.d5/cspeed
    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
    thisOctal%molAbundance(subcell) = molAbundance

!    if (n2max < 0.d0) then
!       thisOctal%nh2(subcell) = abs(n2max)
!    else

!       if (x < bigL/2.d0) then
!          thisOctal%nh2(subcell) = 2.d0 * n2max * x / bigL
!       else
!          thisOctal%nh2(subcell) = 2.d0 * n2max * (bigL - x)/bigL
!       endif
!    endif
    thisOctal%rho(subcell) = thisOctal%nh2(Subcell) * 2.d0 * mHydrogen

  end subroutine calcArbitrary
  
  subroutine WaterBenchmark1(thisOctal, subcell)

    use inputs_mod, only : molAbundance !, amr2d

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    logical, save :: firsttime = .true.
    integer, parameter :: nr = 200
    real,save :: r(nr)
    real :: r1
    integer :: i

    if (firsttime) then
       open(31, file="grid.dat", status="old", form="formatted") ! Model 2 in the Hogerheijde 2000 paper. 
       do i = 1,nr                                             
          read(31,*) r(i)
       enddo
       r = r * 3.08568025e8 ! (pc -> 10^10cm)
       close(31)
       firsttime = .false.
    endif

    r1 = real(modulus(subcellCentre(thisOctal,subcell)))
    thisOctal%temperature(subcell) = 1e-20
    if(r1 > r(nr) .or. r1 < r(1)) then 
       thisOctal%nh2(subcell) = 1.e-20
       thisOctal%rho(subcell) = 1.e-20 * 2. * mhydrogen
    endif

    thisOctal%microTurb(subcell) = 1d-20
    thisOctal%molAbundance(subcell) = molabundance

    if ((r1 > r(1)).and.(r1 < r(nr))) then

       thisOctal%nh2(subcell) = 1e4
       thisOctal%rho(subcell) = thisOctal%nh2(subcell)*2.*mhydrogen
       thisOctal%temperature(subcell) = 40.
       thisOctal%velocity(subcell) = vector(0d-20,0d-20,0d-20)
       thisOctal%microturb(subcell) = max(1d-9,sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) &
                                      / (18.0 * amu))) / (cspeed * 1e-5)) ! mu is subsonic turbulence
    endif
  end subroutine WaterBenchmark1

  subroutine WaterBenchmark2(thisOctal,subcell)

    use inputs_mod, only : molAbundance !, amr2d

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    logical, save :: firsttime = .true.
    integer, parameter :: nr = 200
    real,save :: r(nr)
    real :: r1
    integer :: i

    if (firsttime) then
       open(31, file="grid.dat", status="old", form="formatted") ! Model 2 in the Hogerheijde 2000 paper. 
       do i = 1,nr                                             
          read(31,*) r(i)
       enddo
       r = r * 3.08568025e8 ! (pc -> 10^10cm)
       close(31)
       firsttime = .false.
    endif

    r1 = real(modulus(subcellCentre(thisOctal,subcell)))

    if(r1 > r(nr) .or. r1 < r(1)) then 
       thisOctal%nh2(subcell) = 1.e-20
       thisOctal%rho(subcell) = 1.e-20 * 2. * mhydrogen
    endif

    thisOctal%temperature(subcell) = 40.
    thisOctal%microTurb(subcell) = 1d-10
    thisOctal%molAbundance(subcell) = 1e-20
    thisOctal%velocity(subcell) = vector(1d-20,1d-20,1d-20)

    if ((r1 > r(1)).and.(r1 < r(nr))) then
       thisOctal%molAbundance(subcell) = molabundance
       thisOctal%nh2(subcell) = 1e4
       thisOctal%rho(subcell) = thisOctal%nh2(subcell)*2.*mhydrogen
       thisOctal%velocity(subcell) = WaterBenchmarkvelocity(subcellCentre(thisOctal,subcell))
       thisOctal%microturb(subcell) = max(1d-9,sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) &
                                      / (18.0 * amu))) / (cspeed * 1e-5)) ! mu is subsonic turbulence
    endif

   CALL fillVelocityCorners(thisOctal,WaterBenchmarkVelocity)
  end subroutine WaterBenchmark2

  subroutine AGBStarBenchmark(thisOctal,subcell)

    use inputs_mod, only : molAbundance

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    logical, save :: firsttime = .true.
    integer, parameter :: nr = 100
    real,save :: r(nr), nh2(nr), td(nr), tg(nr), v(nr), mu(nr)
    real :: r1!, mu1 , t1, t2
    real(double) :: v1 !, vDopp
    integer :: i

    type(VECTOR) :: vel

    if (firsttime) then
       open(31, file="mc_100.dat", status="old", form="formatted")
       do i = nr,1,-1                                             
          read(31,*) r(i), nh2(i), tg(i), td(i), v(i), mu(i)
       enddo
       r = r * 1e-10
       mu = mu / 2.35482 ! FWHM -> sigma
       close(31)
       firsttime = .false.
    endif

    r1 = real(modulus(subcellCentre(thisOctal,subcell)))

    if(r1 > r(nr)) then
       thisOctal%nh2(subcell) = 1.e-20
       thisOctal%rho(subcell) = 1.e-20 * 2. * mhydrogen
       thisOctal%temperaturegas(subcell) = real(tcbr)
       thisOctal%temperaturedust(subcell) = real(tcbr)
       thisOctal%microturb(subcell) = 1e-8
    elseif(r1 < r(1)) then
       thisOctal%nh2(subcell) = nh2(1)
       thisOctal%rho(subcell) = thisOctal%nh2(subcell) * 2. * mhydrogen
       thisOctal%temperaturegas(subcell) = 3000.
       thisOctal%temperaturedust(subcell) = 3000.
       thisOctal%microturb(subcell) = sqrt((2.d-10 * kerg * thisOctal%temperaturegas(subcell) / (18.0 * amu)) + mu(1)**2) &
                                     / (cspeed * 1e-5)
    endif

    thisOctal%molAbundance(subcell) = molabundance
    thisOctal%velocity(subcell) = vector(1d-20,1d-20,1d-20)

    if ((r1 > r(1)).and.(r1 < r(nr))) then

       thisOctal%nh2(subcell) = 2.4652426912e19*r1**(-2.5)
       thisOctal%rho(subcell) = thisOctal%nh2(subcell)*2.*mhydrogen

       thisOctal%temperaturegas(subcell) = 4.93117864e8 / sqrt(r1*r1*r1) ! 5e8
       thisOctal%temperaturedust(subcell) = 9.001104e6 / r1 ! 9e6 

       if(r1 .gt. (6e0 * 1.49598e3) ) then ! if r gt than 6 AU out this is outward velocity
          v1 = 8.0673550611e-3 * r1**(0.65) * 1d5 ! km -> cm
          thisOctal%microturb(subcell) = sqrt((2.d-10 * kerg * thisOctal%temperaturegas(subcell) / (18.0 * amu)) + 1.d0**2) &
                                         / (cspeed * 1e-5) ! mu is subsonic turbulence
       else
          v1 = 1d-20
          thisOctal%microturb(subcell) = sqrt((2.d-10 * kerg * thisOctal%temperaturegas(subcell) / (18.0 * amu)) + 3.d0**2) &
                                         / (cspeed * 1e-5) ! mu is subsonic turbulence
       endif

       vel = subcellCentre(thisOctal, subcell)
       call normalize(vel)
       thisOctal%velocity(subcell) = (v1 * vel)/cspeed
    endif

    thisOctal%temperature(subcell) = thisOctal%temperaturegas(subcell)

   CALL fillVelocityCorners(thisOctal,AGBStarVelocity)
  end subroutine AGBStarBenchmark

  subroutine clumpyAgb(thisOctal,subcell)

    use inputs_mod, only : rinner, router, vterm, mdot

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rvec
    real(double) :: r
    rVec = subcellCentre(thisOctal, subcell)
    r = modulus(rVec)
    thisOctal%rho(subcell) = 1.d-30
       thisOctal%temperature(subcell) = 10.
    if ((r > rinner).and.(r < router)) then
       thisOctal%rho(subcell) = mdot / (fourpi* (r*1.d10)**2 * vterm)
    endif


  end subroutine ClumpyAgb

  subroutine turbBox(thisOctal, subcell)

    type(octal) :: thisOctal
    integer :: subcell

    thisOctal%rho(subcell) = 100.d0 * mHydrogen
    thisOctal%velocity(subcell) = VECTOR(0.d0,0.d0,0.d0)
    thisOctal%temperature(subcell) = 10.d0
    thisOctal%iEquationOfState(subcell) = 2

  end subroutine turbBox

  subroutine iras04158(thisOctal, subcell)

    use inputs_mod, only : router, rinner
    use density_mod, only: iras04158Disc
    
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    real :: r
    TYPE(vector) :: rVec
    character(len=10) :: out

    rvec = subcellCentre(thisOctal,subcell)
    r = real(sqrt(rVec%x**2 + rVec%y**2))

    thisOctal%temperature(subcell) = real(tcbr )
    thisOctal%nh2(subcell) = 1d-20
    thisOctal%rho(subcell) = thisOctal%nh2(subcell) * 2.*mhydrogen
    thisOctal%molabundance(subcell) = 1e-20
    thisOctal%temperaturedust(subcell) = real(tcbr)
    thisOctal%temperaturegas(subcell) = real(tcbr)
    thisOctal%microturb(subcell) = 1d-10
    thisOctal%Velocity(subcell) = VECTOR(1d-20,1d-20,1d-20)

    if ((r .lt. rOuter) .and. (r .gt. rinner)) then
      
       out ='nh2'
       thisOctal%rho(subcell) = iras04158Disc(rvec)
       thisOctal%nh2(subcell) = thisOctal%rho(subcell) / (2.*mhydrogen)
       
       out ='abundance'
       thisOctal%molabundance(subcell) = real(readparameterfrom2dmap(rvec,out,.true.))

       out='td'
!       thisOctal%temperaturedust(subcell) = readparameterfrom2dmap(rvec,out,.true.)
       thisOctal%temperaturedust(subcell) = 10.
       thisOctal%temperaturegas(subcell) = thisOctal%temperaturedust(subcell)
       thisOctal%temperature(subcell) = thisOctal%temperaturedust(subcell)
       
       thisOctal%microturb(subcell) = max(1d-8,sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) / (28.0 * amu)) + 0.3**2) &
                                      / (cspeed * 1e-5)) ! mu is 0.3km/s subsonic turbulence
       thisOctal%velocity(subcell) = keplerianVelocity(rvec)
    else
       thisOctal%velocity = VECTOR(1d-20,1d-20,1d-20)
    endif
!    write(*,*) thisOctal%temperature(subcell), thisOctal%rho(subcell), modulus(thisOctal%velocity(subcell))
   CALL fillVelocityCorners(thisOctal,keplerianVelocity)
 end subroutine iras04158
 
 function readparameterfrom2dmap(point,output,dologint) result(out)
   
   TYPE(vector), INTENT(IN) :: point
   logical, save :: firsttime = .true.
   integer, parameter :: nr = 140, nz = 60
   real,save :: rgrid(nr), zgrid(nr,nz), nh2(nr,nz), td(nr,nz), abundance(nr,nz)
   real :: junk, zarray(nz,2)
   integer :: r1, r2, z11, z12, z21, z22
   integer :: ri, zi
   real :: rmax, rmin, z1max, z1min, z2min, z2max, r, z, u, v
   real :: zminint, zmaxint, gradmin, gradmax, d
   real(double) :: out
   character(len=10) :: output

   logical, optional :: dologint
   logical :: logint

   if(present(dologint)) then
      logint = dologint
   else
      logint = .false.
   endif

   if (firsttime) then
      open(31, file="molecularmap.dat", status="old", form="formatted") ! Model 2 in the Hogerheijde 2000 paper. 
      
      do ri = 1, nr
         read(31,*) rgrid(ri), zgrid(ri,1), td(ri,1), nh2(ri,1), abundance(ri,1)
         do zi = 2,nz
            read(31,*) junk, zgrid(ri,zi), td(ri,zi), nh2(ri,zi), abundance(ri,zi)
         enddo
      enddo
             
      rgrid = rgrid * 1e-10
      zgrid = zgrid * 1e-10
      close(31)
      firsttime = .false.
   endif

    r = real(sqrt(point%x**2+point%y**2))
    z = real(abs(point%z))

    if(r .gt. rgrid(1) .and. r .lt. rgrid(nr)) then
       call locate(rgrid,size(rgrid),r,r1)
    elseif(r .lt. rgrid(1)) then
       r1 = 1
    else
       r1 = nr
    endif
       
    if( (r .gt. rgrid(r1))) then
       r2 = min(r1+1,nr)
    else
       r2 = r1
       r1 = max(1,r2-1)
    endif

    zarray(:,1) = zgrid(r1,:)
    zarray(:,2) = zgrid(r2,:)

    if((z .gt. zarray(1,1) .and. z .lt. zarray(nz,1))) then
       call locate(zarray(:,1),size(zarray(:,1)),z,z11)
    elseif(z .lt. zarray(1,1)) then
       z11 = 1
    else
       z11 = nz
    endif

    if(z .gt. zarray(z11,1)) then
       z12 = min(z11 + 1,nz)
    else
       z12 = z11
       z11 = max(1,z12 - 1)
    endif

    if((z .gt. zarray(1,2)) .and. (z .lt. zarray(nz,2))) then
       call locate(zarray(:,2),size(zarray(:,2)),z,z21)
    elseif(z .lt. zarray(1,2)) then
       z21 = 1
    else
       z21 = nz
    endif

    if( (z .gt. zarray(z21,2))) then
       z22 = min(z21 + 1,nz)
    else
       z22 = z21
       z21 = max(1,z22 - 1)
    endif

     rmax = rgrid(r2)
     rmin = rgrid(r1)
     z1min = zarray(z11,1)
     z1max = zarray(z12,1)
     z2min = zarray(z21,2)
     z2max = zarray(z22,2)

     if(logint) then
        r = log(r)
        rmax= log(rmax) 
        rmin= log(rmin)
        z = log(z)
        z1max= log(z1max) 
        z2max= log(z2max) 
        z1min= log(z1min)
        z2min= log(z2min)
     endif

     if(r1 .ne. r2) then
!        d = sqrt((rmax-rmin)**2+(z2min-z1min)**2)
        u = (r - rmin) / (rmax - rmin)
     else
        u = 0
     endif
     
!     if((zarray(z12,1) .ne. zarray(z11,1)) .and. (zarray(z22,2) .ne. zarray(z21,2)) .and. (r1 .ne. r2)) then
     if((z12 .ne. z11) .and. (z22 .ne. z21) .and. (r1 .ne. r2)) then
        if((z2min-z)*(z2max-z) .gt. 0.d0 .or. (z1min-z)*(z1max-z) .gt. 0.d0) stop
        gradmin = (z2min - z1min) / (rmax - rmin)
        gradmax = (z2max - z1max) / (rmax - rmin)
        zminint = z1min + (r - rmin) * gradmin
        zmaxint = z1max + (r - rmin) * gradmax
        d = zmaxint - zminint
        v = (z - zminint) / d

     elseif((z1max .ne. z1min) .and. (z2max .ne. z2min) .and. (r1 .eq. r2)) then
        zminint = z1min
        zmaxint = z1max
        d = zmaxint - zminint
        v = (z - zminint) / d

     elseif(z .gt. z1max .or. z .lt. z2min) then ! possible for v to be greater than 1 so limit it
        if (rmax /= rMin) then
           gradmin = (z2min - z1min) / (rmax - rmin)
           gradmax = (z2max - z1max) / (rmax - rmin)
        else
           gradmin = 0.
           gradmax = 0.
        endif
        zminint = z1min + (r - rmin) * gradmin
        zmaxint = z1max + (r - rmin) * gradmax
        d = zmaxint - zminint
        if (d /= 0.) then
           v = max(min((z - zminint) / d,2.e0),-1.e0)
        else
           v = 0.d0
        endif

     else

        v = 0

     endif

     if(output .eq. 'abundance') then
        if(z .gt. 1.1 * z1max .and. z .gt. 1.1 * z2max) then
           out = 1e-13
        elseif( .not. logint) then
           out = (1-v) * ((1-u) * abundance(r1,z11) + u * abundance(r2,z21)) +  &
                     v * ((1-u) * abundance(r1,z12) + u * abundance(r2,z22))
        else
           out = 10.d0**((1-v) * ((1-u) * log10(abundance(r1,z11)) + u * log10(abundance(r2,z21))) +  &
                             v * ((1-u) * log10(abundance(r1,z12)) + u * log10(abundance(r2,z22))))
        endif
!        out = v
     
       
     elseif(output .eq. 'nh2') then
        if(z .gt. 1.1 * z1max .and. z .gt. 1.1 * z2max) then
           out = 1e-5
        elseif( .not. logint) then
           out = (1-v) * ((1-u) * nh2(r1,z11) + u * nh2(r2,z21)) + &
                     v * ((1-u) * nh2(r1,z12) + u * nh2(r2,z22))        
        else
           out = 10.d0**((1-v) * ((1-u) * log10(nh2(r1,z11)) + u * log10(nh2(r2,z21))) +  &
                             v * ((1-u) * log10(nh2(r1,z12)) + u * log10(nh2(r2,z22))))
        endif

     elseif(output .eq. 'td') then
        if(z .gt. 1.1 * z1max .and. z .gt. 1.1 * z2max) then
           out = tcbr
        elseif( .not. logint) then
           out = 10.d0**((1-v) * ((1-u) * log10(td(r1,z11)) + u * log10(td(r2,z21))) +  &
                             v * ((1-u) * log10(td(r1,z12)) + u * log10(td(r2,z22))))
        else
           out = (1-v) * ((1-u) * td(r1,z11) + u * td(r2,z21)) + &
                     v * ((1-u) * td(r1,z12) + u * td(r2,z22))        
        endif

     else
        write(*,*) "Undefined parameter"
        out = 1d-20
     endif
end function readparameterfrom2dmap

  subroutine molecularFilamentFill(thisOctal,subcell)

    use inputs_mod, only : molAbundance, molecularphysics, photoionPhysics, hydrodynamics

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    real(double) :: r0, r1, rhoc
    TYPE(VECTOR) :: cellCentre


    rhoc = 5.d-19
    r0 = sqrt(2.d0*(10.d0*kerg/(2.3d0*mHydrogen))/(pi*bigG*rhoc))/1.d10

    cellCentre = subcellCentre(thisOctal, subcell)
    r1 = modulus(VECTOR(cellCentre%x, 0.d0, 0.d0))
    thisOctal%temperature(subcell) = 10.d0
    thisOctal%rho(subcell) = rhoc * (1.d0 + (r1/r0)**2)**(-2.d0)
    if (molecularPhysics) then
       thisOctal%molAbundance(subcell) = molAbundance
       thisOctal%microTurb(subcell) = 1.e5/cspeed
       thisOctal%nh2(subcell) = thisOctal%rho(subcell)/(2.*mhydrogen)
       thisOctal%cornerVelocity =  VECTOR(0.d0,0.d0,0.d0)
    endif

    thisOctal%velocity(subcell) = VECTOR(0.d0,0.d0,0.d0)

    thisOctal%inFlow(subcell) = .true.

    if (photoionPhysics) then

       
       thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
       thisOctal%ne(subcell) = thisOctal%nh(subcell)
       thisOctal%nhi(subcell) = 1.e-5
       thisOctal%nhii(subcell) = thisOctal%ne(subcell)
       thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
       
       thisOctal%ionFrac(subcell,1) = 1.               !HI
       thisOctal%ionFrac(subcell,2) = 1.e-10           !HII
       if (SIZE(thisOctal%ionFrac,2) > 2) then      
          thisOctal%ionFrac(subcell,3) = 1.            !HeI
          thisOctal%ionFrac(subcell,4) = 1.e-10        !HeII
          
       endif
    endif
    if (hydrodynamics) then
       thisOctal%gamma(subcell) = 1.0
       thisOctal%iEquationOfState(subcell) = 1
    endif


  end subroutine molecularFilamentFill

  subroutine ggtauFill(thisOctal,subcell)

    use inputs_mod, only : molAbundance 

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    real(double) :: r,t0,nh20,H0,H,z, rmin, rmax
    TYPE(VECTOR) :: cellCentre

    CellCentre = subcellCentre(thisOctal, subcell)
!    r = modulus(cellCentre) ! spherical
    r = sqrt(cellcentre%x**2+cellcentre%y**2) ! cylindrical
    rmin = sqrt((cellcentre%x-thisoctal%subcellsize)**2+(cellcentre%y-thisoctal%subcellsize)**2)
    rmax = sqrt((cellcentre%x+thisoctal%subcellsize)**2+(cellcentre%y+thisoctal%subcellsize)**2)
    r = r * 6.68458134e-06! (torus units to 100's of AUs)
    rmin = rmin * 6.68458134e-06! (torus units to 100's of AUs)
    rmax = rmax * 6.68458134e-06! (torus units to 100's of AUs)

    z = 6.68458134e-04 * Cellcentre%z !in AU

    t0 = 30.
    nh20 = 6.3e9
    H0 = 14.55

    H = H0 * r * sqrt(sqrt(r))  !H0*r^(5/4) r in 100s of AU and H in AU

    thisOctal%molabundance = molAbundance

    if(rmin .gt. 1.8 .and. rmax .lt. 8.) then
       thisOctal%nh2(subcell) = nh20 * (r **(-2.75)) * exp(-((z/H)**2))!(in cm-3) 
       thisOctal%rho(subcell) = thisOctal%nh2(subcell) * 2. * mHydrogen

       thisOctal%temperature(subcell) = real(t0 / sqrt(r) )
    else
       thisOctal%nh2(subcell) = 1.d-60
       thisOctal%rho(subcell) = 1.d-60 * 2. * mhydrogen
       
       thisOctal%temperature(subcell) = real(tcbr)
    endif

    CALL fillVelocityCorners(thisOctal,ggtauVelocity)

  end subroutine ggtauFill


  TYPE(vector)  function wrshellVelocity(point)
    use inputs_mod, only : vterm, beta, rinner
    type(vector), intent(in) :: point
    type(vector) :: rvec
    real(double) :: v, r
    rVec = point

    wrShellVelocity = VECTOR(0.d0, 0.d0, 0.d0)
    r = modulus(rVec)
    If ((r > rInner)) then !.and.(r < rOuter)) then
       v = 0.001d5+(vterm-0.001d5)*(1.d0 - rinner/r)**beta
       call normalize(rvec)
       wrshellvelocity = rvec * (v/cSpeed)
    endif

  end function wrshellVelocity

  TYPE(vector) FUNCTION moleBenchVelocity(point)

    type(VECTOR), intent(in) :: point
    logical, save :: firsttime = .true.
    integer, parameter :: nr = 50
    integer :: i
    real(double), save :: r(nr), nh2(nr), junk,t(nr), v(nr) , mu(nr), OneOverrdiff(nr)
    real(double) :: v1, t1
    real(double) :: r1, pos
    real(double), save :: nstartOverlogdiff, noverlogdiff

    if (firsttime) then
       open(31, file="model_1.dat", status="old", form="formatted")
       do i = nr,1,-1
          read(31,*) r(i), nh2(i), junk,t(i), v(i) , mu(i)
       enddo
       r = r * 1d-10
       close(31)
       firsttime = .false.

       OneOverrdiff(nr) = 1d30
       do i = 1, nr-1
          OneOverrdiff(i) = 1.d0 / (r(i+1)-r(i))
       enddo

       nOverlogdiff = dble(nr-1) / (log(r(nr)) - log(r(1)))
       nstartOverlogdiff = log(r(1)) * dble(nr-1) / (log(r(nr)) - log(r(1)))
    endif

    moleBenchVelocity = VECTOR(0d4,0d4,0d4)
      
    r1 = modulus(point)
    
    if ((r1 > r(1)).and.(r1 < r(nr))) then

!       call locate(r, nr, r1, i)
       
!       t1 = (r1 - r(i)) * OneOverrdiff(i)

       pos = (nOverlogdiff * log(r1) - nstartOverlogdiff)
       
       t1 = fraction(pos)
       i = int(pos) + 1
       v1 = (v(i) + t1 * (v(i+1)-v(i)))*1.d5
       moleBenchVelocity = (v1 * OneOvercSpeed / r1) * point
    endif

  end FUNCTION moleBenchVelocity

  TYPE(vector) FUNCTION WaterBenchmarkVelocity(point)

    type(VECTOR), intent(in) :: point
    real :: v1, r1
    type(VECTOR) :: vel

    r1 = real(modulus(point) * 3.24077649e-9) ! 10^10cm -> pc
    WaterBenchmarkVelocity = VECTOR(1d-20,1d-20,1d-20)

    if ((r1 > 0.001) .and. (r1 < 0.1)) then
       v1 = 100. * r1
       vel = point
       call normalize(vel)
       WaterBenchmarkVelocity = (v1*1e5/cSpeed) * vel
    endif

  end FUNCTION WaterBenchmarkVelocity

  TYPE(vector) FUNCTION AGBStarVelocity(point)

    type(VECTOR), intent(in) :: point
    logical, save :: firsttime = .true.
    integer, parameter :: nr = 100
    integer :: i
    real,save :: r(nr)
    real :: v1, r1
    type(VECTOR) :: vel

    if (firsttime) then
       open(31, file="mc_100.dat", status="old", form="formatted")
       do i = nr,1,-1                                             
          read(31,*) r(i)
       enddo
       r = r * 1e-10
       close(31)
       firsttime = .false.
    endif

    r1 = real(modulus(point))
    AGBStarVelocity = VECTOR(1d-20,1d-20,1d-20)

    if ((r1 > r(1)).and.(r1 < r(nr))) then
       v1 = 8.0673550611E-03 * r1**(0.65)
       vel = point
       call normalize(vel)
       AGBStarVelocity = (v1/cSpeed) * vel
    endif

  end FUNCTION AGBStarVelocity

  recursive subroutine hydroVelocityconvert(thisOctal) 
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call hydroVelocityConvert(child)
                exit
             end if
          end do
       else 
          thisOctal%velocity(subcell)%x = (thisOctal%rhou(subcell)/thisOctal%rho(subcell))/cSpeed
          thisOctal%velocity(subcell)%y = (thisOctal%rhov(subcell)/thisOctal%rho(subcell))/cSpeed
          thisOctal%velocity(subcell)%z = (thisOctal%rhow(subcell)/thisOctal%rho(subcell))/cSpeed
       endif
    enddo
 end subroutine hydroVelocityconvert


  TYPE(vector)  function keplerianVelocity(point)
    use inputs_mod, only :  mcore
    type(vector), intent(in) :: point
    type(vector) :: rvec
    real(double) :: v, r
    rVec = point

    keplerianVelocity = VECTOR(1.d-20, 1.d-20, 1.d-20)

    r = sqrt(point%x**2+point%y**2)

!    write(*,*) "r", r
!    write(*,*) "rvec", rvec, modulus(rvec)

!    if ((r > rInner)) then !.and.(r < rOuter)) then
       v = sqrt(2.d0*6.672d-8*mcore/(r*1d10)) ! G in cgs and M in g (from Msun)
!    write(*,*) "v", v
 
       keplerianvelocity = VECTOR(0.d0,0.d0,1.d0) .cross. rvec
      call normalize(keplerianvelocity) 
      keplerianvelocity = (v/cSpeed) * keplerianvelocity        

!    endif
  end function keplerianVelocity


  TYPE(vector)  function ggtauVelocity(point)

    type(vector), intent(in) :: point
    type(vector) :: rvec
    real(double) :: v, r,v0

    rvec = point

    r = sqrt(point%x**2+point%y**2)
    r = r * 6.68458134d-06! (torus units to 100's of AUs)

    ggtauvelocity = VECTOR(1.d-20,1.d-20,1.d-20)
    v0 = 3.3d5

    if(r .gt. 1.6d0 .and. r .lt.  10.d0) then
       v = v0 / sqrt(r)
       call normalize(rvec)
       ggtauvelocity = VECTOR(0.d0,0.d0,1.d0) .cross. rvec
       call normalize(ggtauvelocity)
       ggtauvelocity = (v/cSpeed) * ggtauvelocity        
    else
       ggtauvelocity = VECTOR(v0/cspeed,v0/cspeed,v0/cspeed)
    endif    

  end function ggtauVelocity
 
#ifdef SPH
 TYPE(vector) FUNCTION ClusterVelocity(point)

   use sph_data_class, only : clusterparameter

    type(VECTOR), intent(in) :: point
    integer :: subcell

    call findSubcellLocal(point, recentOctal,subcell)

    clustervelocity = Clusterparameter(point, recentoctal, subcell, theparam = 1) ! use switch for storing velocity
  end FUNCTION ClusterVelocity

  real(double) FUNCTION ClusterDensity(point)
    
    use sph_data_class, only : clusterparameter

    type(VECTOR) :: out
    type(VECTOR), intent(in) :: point
    integer :: subcell

    call findSubcellLocal(point, recentOctal,subcell)
    
    out = Clusterparameter(point, recentoctal, subcell, theparam = 2) ! use switch for storing velocity
    clusterdensity = out%x
  end FUNCTION ClusterDensity
#endif

  subroutine assign_melvin(thisOctal,subcell,grid)

    use density_mod, only: melvinDensity
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    thisOctal%rho(subcell) = melvinDensity(rVec, grid)
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

    thisOctal%ionFrac(subcell,1) = 1.e-10
    thisOctal%ionFrac(subcell,2) = 1.
    thisOctal%ionFrac(subcell,3) = 1.e-10
    thisOctal%ionFrac(subcell,4) = 1.       
    thisOctal%etaCont(subcell) = 0.

!thap this if statement is to try and empty the core                                                                      
    If (thisOctal%rho(subcell) > 3.00d30) then !In cavity and in bad range                                                
       thisOctal%rho(subcell) = 1.d-10 * mhydrogen
       thisOctal%dustTypeFraction(subcell, :) = 1.d-20
       thisOctal%temperature(subcell) = 1.d4
    endif

    if (thisOctal%rho(subcell) > 1.e29)  then ! in cavity
       thisOctal%rho(subcell) = 100.d0 * mHydrogen
       thisOctal%dustTypeFraction(subcell, :) = 1.d-20 ! no dust in cavity
       thisOctal%temperature(subcell) = 1.d4 ! 10,000K in cavity
    endif

    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)

  end subroutine assign_melvin

!Thaw 2D gaussian distribution
  subroutine assign_gaussian(thisOctal,subcell)
    use inputs_mod, only : xplusbound, xminusbound, zplusbound, zminusbound
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    type(VECTOR) :: rVec
    real(double) :: x, gd, z

    x = 0.d0
    z = 0.d0
    xplusbound = 2
    xminusbound = 2
    zplusbound = 1
    zminusbound = 1
    rVec = subcellCentre(thisOctal, subcell)
    x = rVec%x
    z = rVec%z
    gd = 0.1d0
    if(thisOctal%twoD) then
       thisOctal%rho(subcell) =1.d0 + 0.3d0 * exp(-(((x-0.5d0)**2 + (z-0.5d0)**2))/gd**2)
    else if (thisOctal%oneD) then
       thisOctal%rho(subcell) =1.d0 + 0.3d0 * exp(-(((x-0.5d0)**2)/gd**2))
    end if

    thisOctal%energy(subcell) = 1.d0
    thisOctal%velocity(subcell) = VECTOR(0., 0., 0.)
    thisOctal%pressure_i(subcell) = 1.d0
    thisOctal%phi_i(subcell) = 0.d0
    thisOctal%gamma(subcell) = 7.d0/5.d0
    thisOctal%iEquationOfState(subcell) = 0

  end subroutine assign_gaussian

  subroutine assign_hii_test(thisOctal,subcell)
    use inputs_mod, only : xplusbound, xminusbound, yplusbound, yminusbound, zplusbound, zminusbound
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(vector) :: rVec
    real(double) :: eThermal!, numDensity
   
    rVec = subcellCentre(thisOctal,subcell)
    thisOctal%rho(subcell) = (5.d2)*mHydrogen
    thisOctal%temperature(subcell) = 1.d4
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)


    thisOctal%ionFrac(subcell,1) = 1.
    thisOctal%ionFrac(subcell,2) = 1.e-10
    
    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
    
    ethermal = 1.5d0*(1.d0/(0.5d0*mHydrogen))*kerg*thisOctal%temperature(subcell)
!    thisOctal%gamma(subcell) = 5.d0/3.d0
    thisOctal%gamma(subcell) = 1.d0


    thisOctal%pressure_i(subcell) = (thisOctal%rho(subcell)/(0.5d0*mHydrogen))*kerg*thisOctal%temperature(subcell)
   
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
    thisOctal%phi_i(subcell) = 0.d0
    thisOctal%iEquationOfState(subcell) = 1

    zplusbound = 1
    zminusbound = 1
    xplusbound = 1
    xminusbound = 1
    yplusbound = 1
    yminusbound = 1
  end subroutine assign_hii_test

  subroutine assign_radpresstest(thisOctal,subcell)
    use inputs_mod, only : rCavity
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(vector) :: rVec
    real(double) :: eThermal!, numDensity
   
    rVec = subcellCentre(thisOctal,subcell)
    thisOctal%rho(subcell) = (1.d2)*mHydrogen
    if (modulus(rVec) < rCavity) thisOctal%rho(subcell) = 1.d-30
    thisOctal%temperature(subcell) = 10.d0
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)


    thisOctal%ionFrac(subcell,1) = 1.
    thisOctal%ionFrac(subcell,2) = 1.e-10
!    thisOctal%ionFrac(subcell,3) = 1.e-10
!    thisOctal%ionFrac(subcell,4) = 1.       
!    thisOctal%etaCont(subcell) = 0.
    
    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
    
    ethermal = 1.5d0*(1.d0/(mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%gamma(subcell) = 5.d0/3.d0


    thisOctal%pressure_i(subcell) = (thisOctal%rho(subcell)/(mHydrogen))*kerg*thisOctal%temperature(subcell)
   
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
    thisOctal%phi_i(subcell) = 0.d0
    thisOctal%iEquationOfState(subcell) = 1

  end subroutine assign_radpresstest

  subroutine assign_whitney(thisOctal,subcell,grid)

    use density_mod, only: whitneyDensity
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    thisOctal%rho(subcell) = whitneyDensity(rVec, grid)
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine assign_whitney

  subroutine assign_planetgap(thisOctal,subcell,grid)
    use density_mod, only: planetgapDensity

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    thisOctal%rho(subcell) = planetgapDensity(rVec, grid)
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine assign_planetgap

  subroutine assign_toruslogo(thisOctal,subcell)

    use density_mod, only: toruslogoDensity
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    thisOctal%rho(subcell) = toruslogoDensity(rVec)
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine assign_toruslogo

  subroutine calcBenchI(thisOctal,subcell)
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(vector) :: rVec
    real(double) :: r1, r2, rho1, rho2, mass1, mass2, rhoAmbient, rsub
    type(VECTOR) :: pos1, pos2, sourcepos

    rsub = autocm/1.d10
    r1 = 1000.d0*autocm/1.d10
    r2 = 500.d0*autocm/1.d10

    sourcePos = VECTOR(-5.d6, 0.d0, 2.5d6)
    pos1 = VECTOR(2000.d0*autocm/1.d10, 0.d0, 0.d0)
    pos2 = VECTOR(-2000.d0*autocm/1.d10, 0.d0, 0.d0)

    mass1 = 10.d0 * msol
    mass2 = 5.d0 * msol
    rhoambient = 1.d-30


    rVec = subcellCentre(thisOctal,subcell)

    rho1 = mass1 / (fourPi/3.d0 * (r1*1.d10)**3)
    rho2 = mass2 / (fourPi/3.d0 * (r2*1.d10)**3)

    thisOctal%rho(subcell) = rhoAmbient
    thisOctal%dustTypeFraction(subcell,1) = 0.01d0
    if (modulus(rVec-pos1) < r1) then
       thisOctal%rho(subcell) = rho1
    endif
    if (modulus(rVec-pos2) < r2) then
       thisOctal%rho(subcell) = rho2
    endif

    if (modulus(rVec-sourcepos) < rSub) then
       thisOctal%dustTypeFraction(subcell,1) = 1.d-30
    endif

    thisOctal%temperature(subcell) = 10.
  end subroutine calcBenchI


  subroutine shakaraDisk(thisOctal,subcell,grid)
    use density_mod, only: density, shakaraSunyaevDisc
    use inputs_mod, only : rOuter, betaDisc !, rInner, erInner, erOuter, alphaDisc
    use inputs_mod, only : curvedInnerEdge, nDustType, grainFrac, gridDistanceScale, rInner
    use inputs_mod, only : height, hydrodynamics, dustPhysics, mCore, molecular, photoionization
    use inputs_mod, only : rSublimation

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real(double) :: r, h, rd, ethermal, rhoFid, thisRSub,z,fac, rho, sinTheta,v
    TYPE(vector) :: rVec
    
    type(VECTOR),save :: velocitysum
    logical,save :: firsttime = .true.

    if(firsttime) then
       velocitysum = VECTOR(1d-20,1d-20,1d-20)
       firsttime = .false.
    endif

    rVec = subcellCentre(thisOctal,subcell)
    r = real(modulus(rVec))



    thisOctal%inflow(subcell) = .true.
    thisOctal%temperature(subcell) = 100. 
    rd = rOuter / 2.

    if (associated(thisOctal%dustTypeFraction)) thisOctal%dustTypeFraction(subcell,:) = 1.d-20

    thisOctal%rho(subcell) = 1.d-30
    if (associated(thisOctal%nh)) &
         thisOctal%nh(subcell) =  thisOctal%rho(subcell) / mHydrogen
    if (associated(thisOctal%ne)) &
         thisOctal%ne(subcell) = 1.d-5 !thisOctal%nh(subcell)
    if (photoionization) then
       thisOctal%nhi(subcell) = 1.e-5
       thisOctal%nhii(subcell) = thisOctal%ne(subcell)
       thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    endif
    if (associated(thisOctal%ionFrac)) then
       thisOctal%ionFrac(subcell,1) = 1.
       thisOctal%ionFrac(subcell,2) = 1.e-5
       thisOctal%ionFrac(subcell,3) = 1.
       thisOctal%ionFrac(subcell,4) = 1.e-5
       thisOctal%etaCont(subcell) = 0.
    endif

    if (photoionization) then
       thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
       thisOctal%ne(subcell) = 1.d-5!thisOctal%nh(subcell)
       thisOctal%nhi(subcell) = 1.e-5
       thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    endif

    r = real(sqrt(rVec%x**2 + rVec%y**2))

!    if (r < rOuter) then
       thisOctal%rho(subcell) = density(rVec, grid)
! tinkered from 10K - I figured the cooler bits will gently drop but a 
! large no. of cells are close to this temp.
       thisOctal%temperature(subcell) = 10.
       thisOctal%inflow(subcell) = .true.

       if (photoionization) then
          thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
          thisOctal%ne(subcell) = 1.e-5!thisOctal%nh(subcell)
          thisOctal%nhi(subcell) = 1.e-5
          thisOctal%nhii(subcell) = thisOctal%ne(subcell)
          thisOctal%biasCont3D = 1.
          thisOctal%etaLine = 1.e-30
       endif

       h = real(height * (r / (100.d0*autocm/1.d10))**betaDisc)
    
!    endif

    if(molecular) then
 !      if(modulus(rvec) .lt. 1000.) then
          thisOctal%velocity(subcell) = keplerianVelocity(rvec)
          CALL fillVelocityCorners(thisOctal, keplerianVelocity)
 !      else
 !         thisOctal%velocity(subcell) = vector(0.,0.,0.)
 !      endif
    else
       thisOctal%velocity(Subcell) = VECTOR(0.,0.,0.)
    endif

    if (molecular) then
       thisOctal%microturb(subcell) = sqrt((2.d0*kErg*thisOctal%temperature(subcell))/(29.0 * amu))/cspeed
       thisOctal%nh2(subcell) = thisOctal%rho(subcell)/(2.*mhydrogen)
    endif

    if (hydrodynamics) then
       r = sqrt(rVec%x**2 + rVec%y**2)
!       if (r > rOuter*0.99d0) then
!          thisOctal%rho(subcell) = thisOctal%rho(subcell) * exp(-(r-rOuter*0.99d0)/(0.01d0*rOuter))
!       endif
       sinTheta = sqrt(1.d0-(abs(rVec%z)/modulus(rVec))**2)
       v = sqrt(bigG * mSol /(r*1.d10))


       thisOctal%rho(subcell) = max(thisOctal%rho(subcell), 1.e-20_db)
       thisOctal%temperature(subcell) = min(100.,100. * real((r/Rinner)**(-4.d0/3.d0)))
       thisOctal%velocity(subcell) = keplerianVelocity(rvec)
       thisOctal%iEquationOfState(subcell) = 1
       thisOctal%phi_i(subcell) = -bigG * mCore / (modulus(rVec)*1.d10)
       thisOctal%gamma(subcell) = 7.d0/5.d0
       ethermal = real(1.5d0 * (1.d0/(2.d0*mHydrogen)) * kerg * thisOCtal%temperature(subcell))
       thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
       thisOctal%rhoe(subcell) = thisOctal%energy(subcell) * thisOctal%rho(subcell)
       thisOctal%rhov(subcell) = modulus(keplerianVelocity(rVec))*thisOctal%rho(subcell) * (r *gridDistanceScale)*cSpeed * sinTheta
    endif

    if (dustPhysics) then
       rVec = VECTOR(rSublimation*1.001d0, 0.d0, 0.d0)
       rhoFid = shakaraSunyaevDisc(rVec, grid)
       
       thisOctal%DustTypeFraction(subcell,:) = 1.d-10
       rVec = subcellCentre(thisOctal, subcell)
       rho = shakaraSunyaevDisc(rVec, grid)
       thisRsub = 1.01d0 * rSublimation * max(1.d0,(1.d0/(rho/rhoFid)**0.0195)**2)
       r = sqrt(rVec%x**2+rVec%y**2)
       z = rVec%z
       thisOctal%dustTypeFraction(subcell,1:nDustType) = grainFrac(1:nDustType)
       !          write(*,*) r/1496.,thisRsub/1496.d0
       if (modulus(rVec) < rSublimation) then
          thisOctal%dustTypeFraction(subcell,1:nDustType) = 1.d-25
       endif
       
       if (curvedInnerEdge.and.(r < thisRsub).and.(modulus(rVec) < 2.d0*rsublimation)) then
          fac = (thisRsub-r)/(0.002d0*rSublimation)
          thisOctal%dustTypeFraction(subcell,1:nDustType) = max(1.d-20,grainFrac(1:nDustType)*exp(-fac))
       endif
    endif


  end subroutine shakaraDisk


  subroutine fillgridSafier(grid)
    use inputs_mod, only : DW_rMin, DW_rMax, grainfrac, ndusttype, sourceMass, DW_mdot, rSublimation 
    use magnetic_mod, only : safierFits
    type(OCTAL), pointer :: thisOctal
    type(VECTOR) :: rVec
    integer :: subcell
    type(GRIDTYPE) :: grid


    integer :: i, j
    real(double) :: r0, chi, eta, zeta0dash, zeta, psi, rho1, thisRho
    logical :: doDust

    rho1 = 3.5d-15*(dw_mdot/5.d-8)*(sourceMass(1)/msol)**(-0.5d0) * (log(DW_rmax/DW_rmin)/2.5d0)**(-1.d0)
          

    zeta0Dash = tan(60.d0*degToRad)
    thisOctal => grid%octreeRoot
    do i = 1, 1000
       r0 = DW_rMin + (DW_rMax - DW_rMin)*dble(i-1)/999.d0
       do j = 1, 10000
          chi = 1.d-6+10.d0 * dble(j-1)/9999.d0

          call safierFits("C", chi, zeta0dash, 0.d0, zeta, psi, eta)
          
          rVec = VECTOR(r0 * zeta, 0.d0, r0*chi)
          if (inOctal(grid%octreeRoot, rVec)) then
             call findSubcellLocal(rVec, thisOctal,subcell)

             thisRho = eta * rho1 * (r0 / (autocm/1.d10))**(-1.5d0)
             doDust = .false.
             if (thisRho > thisOctal%rho(subcell)) then
                thisOctal%rho(subcell) = max(thisOctal%rho(subcell), thisRho)
                doDust = .true.
             endif

             rVec%z = -rVec%z
             call findSubcellLocal(rVec, thisOctal,subcell)

             if (thisRho > thisOctal%rho(subcell)) then
                thisOctal%rho(subcell) = max(thisOctal%rho(subcell), thisRho)
             endif

             if (doDust) then
                rVec = subcellCentre(thisOctal, subcell)
                r0 = rSublimation
                chi = abs(rVec%z)/r0
                call safierFits("C", chi, zeta0dash, 0.d0, zeta, psi, eta)
                r0 = r0*zeta
                if (rVec%x > r0) then
                   thisOctal%dustTypeFraction(subcell,1:nDustType) = grainFrac(1:nDustType)
                else
                   thisOctal%dustTypeFraction(subcell,1:nDustType) = 1.d-20
                   write(*,*) "seeting to zero"
                endif
             endif
          endif
       enddo
    enddo
  end subroutine fillgridSafier

  subroutine warpedDisk(thisOctal,subcell,grid)
    use density_mod, only: density
    use inputs_mod
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r, h, rd
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    r = real(modulus(rVec))
    thisOctal%inflow(subcell) = .true.
    thisOctal%rho(subcell) = 1.e-30
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    rd = rOuter / 2.
    r = real(sqrt(rVec%x**2 + rVec%y**2))
    if ((r > rInner).and.(r < rOuter)) then
       if (hydroWarp) then
         thisOctal%rho(subcell) = hydroWarpDensity(rVec, grid)
         thisOctal%temperature(subcell) = real(hydroWarpTemperature(rVec, grid))
       else
         thisOctal%rho(subcell) = density(rVec, grid)
         thisOctal%temperature(subcell) = 20.
       end if
       thisOctal%etaCont(subcell) = 0.
       thisOctal%inflow(subcell) = .true.
       h = real(height * (r / (100.d0*autocm/1.d10))**betaDisc)
    endif


    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine warpedDisk

  ! chris (26/05/04)
  subroutine calcPPDiskDensity(thisOctal, subcell, grid)
    use density_mod, only: density
    use inputs_mod
    type(octal), intent(inout) :: thisOctal
    integer, intent(in) :: subcell
    type(gridtype), intent(in) :: grid

    type(vector) :: rVec

    rVec = subcellCentre(thisOctal,subcell)

    thisOctal%rho(subcell) = density(rVec,grid)
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 0.
    thisOctal%chiLine = 0.

  end subroutine calcPPDiskDensity

  subroutine assign_clumpydisc(thisOctal,subcell,grid)
    use density_mod, only: density
    use inputs_mod
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector) :: rVec
    real :: r
    
    rVec = subcellCentre(thisOctal,subcell)
    r = real(modulus(rVec))

    thisOctal%rho(subcell) = 1.e-33
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .false.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

    thisOctal%rho(subcell) = density(rVec, grid)
   if (thisOctal%rho(subcell) > 1.e-33) thisOctal%inflow(subcell) = .true.
  end subroutine assign_clumpydisc

  !
  ! Recursively deletes the sph_particle list
  ! Call after the grid has been set up
  !
  recursive subroutine delete_particle_lists(thisoctal)
    implicit none
    
    type(octal), pointer :: thisoctal
    type(octal), pointer :: pChild  
    integer              :: i             ! loop counters
    
    IF (ASSOCIATED(thisOctal%gas_particle_list)) then
       DEALLOCATE(thisOctal%gas_particle_list)
       nullify(thisOctal%gas_particle_list)
    endif

    if ( thisOctal%nChildren > 0) then
       do i = 1, thisOctal%nChildren
          pChild => thisOctal%child(i)
          call delete_particle_lists(pChild)
       end do
    end if
    
  end subroutine delete_particle_lists
  
  
#ifdef SPH
  !
  ! Assign the indecies to root node in the grid
  !
  
  subroutine copy_sph_index_to_root(grid)
    use sph_data_class, only: get_npart
    implicit none
    type(gridtype), intent(inout) :: grid    
    !
    integer :: i, n

    ! extracting the number of gas particles in sph
    n = get_npart()

    ! initialize the list
    ALLOCATE(grid%octreeRoot%gas_particle_list(n))

    ! Actually the initial indecies are just 1, 2, 3 ..n
    do i = 1, n
       ! assigin the values to the root
       grid%octreeRoot%gas_particle_list(i) = i
    end do
       
  end subroutine copy_sph_index_to_root
#endif
    
  recursive subroutine scaleDensityAMR(thisOctal, scaleFac)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  real(double) :: scaleFac
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call scaleDensityAMR(child, scaleFac)
                exit
             end if
          end do
       else
          thisOctal%rho(subcell) = max(1.d-29,thisOctal%rho(subcell) * scaleFac)
       endif
    enddo
  end subroutine scaleDensityAMR

  recursive subroutine findTotalMass(thisOctal, totalMass, totalMassTrap, totalMassMol, minRho, maxRho)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double), intent(inout) :: totalMass
  real(double),optional, intent(inout) :: totalMassTrap, totalMassMol
  real(double),optional, intent(inout) :: minRho, maxRho
  real(double) :: dv, rhoMol
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                if (present(totalMassMol)) then
                   call findtotalMass(child, totalMass, totalmasstrap, totalMassMol=totalMassMol, minRho=minRho, maxRho=maxRho)
                else
                   call findtotalMass(child, totalMass, totalmasstrap, minRho, maxRho)
                endif
                exit
             end if
          end do
       else

          dv = cellVolume(thisOctal, subcell)
          totalMass = totalMass + (1.d30)*thisOctal%rho(subcell) * dv
          if(present(totalmasstrap) .and. associated(thisOctal%cornerRho)) then 
             totalMassTrap = totalMassTrap + (1.d30) * averagerhofromoctal(thisoctal,subcell) * dv
          end if

          if (present(totalMassMol).and.associated(thisOctal%molabundance)) then
             rhoMol =  thisOctal%molabundance(subcell) * thisOctal%nh2(subcell)*( 28.0_db * mhydrogen )
             totalMassMol = totalMassMol + (1.d30) * dv * rhoMol
          endif

          if (PRESENT(minRho)) then
             if (thisOctal%rho(subcell) > 1.d-40) then
                minRho = min(dble(thisOctal%rho(subcell)), minRho)
             endif
          endif
          if (PRESENT(maxRho)) maxRho = max(dble(thisOctal%rho(subcell)), maxRho)
       endif
    enddo
  end subroutine findTotalMass


  ! find the sum of temperature
  recursive subroutine sum_temp_mass(thisOctal, totalTemperature, &
       totalMassTemperature, TotalMass, ncell)
    implicit none
    type(octal), pointer   :: thisOctal
    real(double), intent(inout) :: totalTemperature      ! simple total mass
    real(double), intent(inout) :: totalMassTemperature
    real(double), intent(inout) :: totalMass
    integer, intent(inout) :: ncell ! number of cells used for this calculation
    !
    type(octal), pointer  :: child
    real(double) :: M, T
    integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call sum_temp_mass(child, totalTemperature, totalMassTemperature, &
                     totalMass, ncell)
                exit
             end if
          end do
       else
          if (thisOctal%inFlow(subcell)) then
             M = (1.d30)*thisOctal%rho(subcell) * thisOctal%subcellSize**3
             T = thisOctal%temperature(subcell)
             totalTemperature = totalTemperature + T
             totalMassTemperature = totalMassTemperature  + M*T
             totalMass = totalMass  + M
             ncell = ncell + 1
          endif
       endif
    enddo
  end subroutine sum_temp_mass


  subroutine spiralWindSubcell(thisOctal, subcell)
    use ostar_mod, only: returnSpiralFactor
    use inputs_mod
    real(double) :: v, r, rhoOut
    type(VECTOR) :: rVec, rHat
    type(VECTOR) :: octVec
    integer :: subcell
    type(OCTAL) :: thisOctal
    character(len=80) :: opacityFilename
    real, allocatable :: rgrid(:), sigma(:), t(:), vgrid(:), eta(:), chi(:), etal(:), chil(:), escat(:)
    integer :: i, j

    rVec = subcellCentre(thisOctal,subcell)
    opacityFilename = "opacity.dat"

    open(20,file=opacityFilename, status="old", form="formatted")
    read(20,*) nr

    allocate(rgrid(1:nr))
    allocate(sigma(1:nr))
    allocate(t(1:nr))
    allocate(vgrid(1:nr))
    allocate(eta(1:nr))
    allocate(etal(1:nr))
    allocate(chi(1:nr))
    allocate(chil(1:nr))
    allocate(escat(1:nr))

!    write(*,*) "Reading opacities from: ",trim(opacityFilename)
    do j = 1, nr
       i = nr - j + 1
       read(20,*) rgrid(i), t(i), sigma(i), vgrid(i), eta(i), chi(i), escat(i), &
            etal(i), chil(i)
    enddo
    close(20)

    r = modulus(rVec)

    thisOctal%inFlow(subcell) = .false.
    thisOctal%velocity(subcell) = VECTOR(0., 0., 0.)
    if (r > rCore) then
       r = modulus(rVec)
       v = v0 + (vTerm - v0) * (1. - rCore/r)**beta
       rhoOut = mdot / (fourPi * (r*1.e10)**2 * v)
       rHat = rVec / r
       octVec = rVec
       thisOctal%rho(subcell) = rhoOut* returnSpiralFactor(octVec, 10.*rCore/real(twoPi))
       thisOctal%temperature(subcell) = Teff
       thisOctal%velocity(subcell)  = (v/cSpeed) * rHat
       thisOctal%inFlow(subcell) = .true.
       thisOctal%etaLine(subcell) =  logInterp(etal, nr, rgrid, real(r))
       thisOctal%etaCont(subcell) =  1.e-30
       thisOctal%chiLine(subcell) =  logInterp(chil, nr, rgrid, real(r))
       thisOctal%kappaSca(subcell,1) = logInterp(etal, nr, rgrid,  real(r))/1.e10
       thisOctal%kappaAbs(subcell,1) = 1.e-30
    endif

    CALL fillVelocityCorners(thisOctal,ostarVelocity)

  end subroutine spiralWindSubcell

  subroutine windSubcell(thisOctal, subcell)
    use inputs_mod, only : v0, vterm, beta, twind, mdot, rcore
    real(double) :: v, r, rhoOut
    type(VECTOR) :: rHat, rvec
    type(VECTOR) :: octVec
    integer :: subcell
    type(OCTAL) :: thisOctal

    rVec = subcellCentre(thisOctal, subcell)
    r = modulus(rVec)

    thisOctal%inFlow(subcell) = .false.
    thisOctal%velocity(subcell) = VECTOR(0., 0., 0.)
    thisOctal%temperature(subcell) = Twind
    thisOctal%rho(subcell) = 1.d-30
    if (r > rCore) then
       r = modulus(rVec)
       v = v0 + (vTerm - v0) * (1. - rCore/r)**beta
       rhoOut = mdot / (fourPi * (r*1.e10)**2 * v)
       rHat = rVec / r
       octVec = rVec
       thisOctal%rho(subcell) = rhoOut
       thisOctal%velocity(subcell)  = (v/cSpeed) * rHat
       thisOctal%inFlow(subcell) = .true.
    endif

    CALL fillVelocityCorners(thisOctal,ostarVelocity)

  end subroutine WindSubcell


  TYPE(vector) FUNCTION ostarVelocity(point)
    use inputs_mod, only : v0, vterm, beta, rcore
    real(double) :: r1
    type(VECTOR), intent(in) :: point
    type(VECTOR) rHat
    real :: v

    r1 = modulus(point)

    rhat = point
    call normalize(rHat)
    if (r1 > rCore) then
       v = real(v0 + (vTerm - v0) * (1. - rCore/r1)**beta)
       ostarVelocity  = (dble(v/cSpeed)) * rHat
    else
       ostarVelocity  = VECTOR(0., 0., 0.)
    endif
  end FUNCTION ostarVelocity



  SUBROUTINE checkAMRgrid(grid,checkNoctals)
    ! does some checking that the cells in an AMR grid are
    !   set up and linked to correctly.
    use inputs_mod, only : amr1d, spherical, hydrodynamics
    TYPE(gridType), INTENT(IN) :: grid
    LOGICAL, INTENT(IN) :: checkNoctals ! whether to confirm grid%nOctals
    
    TYPE(OCTAL), POINTER :: rootOctal
    INTEGER :: nOctals

    nOctals = 1
    rootOctal => grid%octreeRoot
    
    CALL checkAMRgridPrivate(grid=grid,              &
                             thisOctal=rootOctal,    &
                             thisDepth=1,            &
                             thisParent=NULL(),      &
                             thisParentSubcell=-999, &
                             nOctals=nOctals)

    IF ( checkNoctals .AND. nOctals /= grid%nOctals ) THEN
      PRINT *, "In checkAMRgrid, nOctals mismatch:"
      PRINT *, "  nOctals = ", nOctals
      PRINT *, "  grid%nOctals = ", grid%nOctals
      PRINT *, "Ignoring..."
    END IF

  CONTAINS
  
    RECURSIVE SUBROUTINE checkAMRgridPrivate(grid,thisOctal, thisDepth ,&
                                             thisParent,thisParentSubcell,  &
                                             nOctals)
      TYPE(gridType), INTENT(IN) :: grid
      TYPE(OCTAL), INTENT(IN), TARGET :: thisOctal
      INTEGER, INTENT(IN) :: thisDepth
      TYPE(OCTAL), POINTER :: thisParent
      INTEGER, INTENT(IN) :: thisParentSubcell
      INTEGER, INTENT(INOUT) :: nOctals

      INTEGER :: iSubcell, iIndex
      TYPE(OCTAL), POINTER :: thisOctalPointer
      REAL(oct) :: sizeRatio

      nOctals = nOctals + 1
      
      IF ( thisOctal%nDepth /= thisDepth ) THEN
        PRINT *, "Error: In checkAMRgridPrivate, depth mismatch"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
      END IF
      
      ! check parent (except at root of tree)
      IF ( thisDepth /= 1 ) THEN
        IF ( .NOT. ASSOCIATED(thisOctal%parent,thisParent) ) THEN
          PRINT *, "Error: In checkAMRgridPrivate, mismatch with parent argument"
          CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
        END IF
      END IF
       
      IF ( thisDepth /= 1 ) THEN
         IF ( thisOctal%parentSubcell /= thisParentSubcell ) THEN
            PRINT *, "Error: In checkAMRgridPrivate, parentSubcell mismatch"
            CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
         END IF
      END IF
      
      IF ( thisOctal%nDepth > grid%maxDepth ) THEN
        PRINT *, "Error: In checkAMRgridPrivate, maxDepth exceeded"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
      END IF
      
      IF ( thisOctal%subcellSize < grid%halfSmallestSubcell ) THEN
        PRINT *, "Error: In checkAMRgridPrivate, halfSmallestSubcell exceeded"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals) 
      END IF
       
      IF ( thisOctal%nChildren > thisOctal%maxChildren ) THEN
        PRINT *, "Error: In checkAMRgridPrivate, maxChildren exceeded"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals) 
      END IF

      IF ( ASSOCIATED(thisOctal%child) .AND. (thisOctal%nChildren == 0) ) THEN
        PRINT *, "Error: In checkAMRgridPrivate, thisOctal%child shouldn't be associated"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
      END IF
      
      IF ( .NOT. ASSOCIATED(thisOctal%child) .AND. (thisOctal%nChildren > 0) ) THEN
        PRINT *, "Error: In checkAMRgridPrivate, thisOctal%child should be associated"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
      END IF

      ! check that invalid children are not set:
      IF ( ANY( thisOctal%hasChild(thisOctal%maxChildren+1:SIZE(thisOctal%hasChild)) )) THEN
        PRINT *, "Error: In checkAMRgridPrivate, invalid hasChild variables set"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals) 
      END IF
      
      ! check that the number of children agree 
      IF ( COUNT(thisOctal%hasChild) /= thisOctal%nChildren ) THEN
        PRINT *, "Error: In checkAMRgridPrivate, nChildren does not match hasChild"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals) 
      END IF

      IF ( thisOctal%nChildren > 0 ) THEN

        ! see if %child is sized correctly
        IF ( SIZE(thisOctal%child) /= thisOctal%nChildren ) THEN
          PRINT *, "Error: In checkAMRgridPrivate, %child has wrong size"
          CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals) 
        END IF

        DO iSubcell = 1, thisOctal%maxChildren

          IF ( .NOT. thisOctal%hasChild(iSubcell) ) CYCLE
          
          ! find the correct index of %child for this child
          DO iIndex = 1, thisOctal%nChildren
            IF ( thisOctal%indexChild(iIndex) == iSubcell ) EXIT
            
            IF ( iIndex == thisOctal%nChildren ) THEN
              ! shouldn't get here
              PRINT *, "Error: In checkAMRgridPrivate, indexChild not found"
              PRINT *, "       for iSubcell = ",iSubcell
              CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals) 
            END IF
          END DO

          ! now we know the correct %child variable, let's do some tests

          IF ( thisOctal%child(iIndex)%parentSubcell /= iSubcell ) THEN
            PRINT *, "Error: In checkAMRgridPrivate, child's parentSubcell doesn't match:"
            PRINT *, "       thisOctal%child(iIndex)%parentSubcell = ",thisOctal%child(iIndex)%parentSubcell
            print *,         thisOctal%child(1:thisOctal%nChildren)%parentSubcell
            PRINT *, "       iIndex = ", iIndex
            PRINT *, "       iSubcell = ", iSubcell
            CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
          END IF

          ! see if the child's coordinates really lie in the parent subcell
          IF ( .NOT. inSubcell(thisOctal,iSubcell,point=thisOctal%child(iIndex)%centre) ) THEN
             if (.not.(hydrodynamics.and.amr1d.and.spherical)) then
            PRINT *, "Error: In checkAMRgridPrivate, child isn't in parentSubcell"
            PRINT *, "       thisOctal%centre = ",thisOctal%centre 
            PRINT *, "       iSubcell = ", iSubcell
            PRINT *, "       subcellCentre(thisOctal,iSubcell) = ",subcellCentre(thisOctal,iSubcell)
            PRINT *, "       thisOctal%child(iIndex)%centre = ",thisOctal%child(iIndex)%centre 
            PRINT *, "       iIndex = ", iIndex
            CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
            endif
          END IF

          ! see if the child is of the correct size
          sizeRatio = thisOctal%subcellSize / thisOctal%child(iIndex)%subcellSize
          sizeRatio = sizeRatio / 2.0_oc
          IF ( ABS(sizeRatio-1.0_oc) > 0.1 ) THEN
            PRINT *, "Error: In checkAMRgridPrivate, size of child wrong:"
            PRINT *, "       thisOctal%subcellSize = ", thisOctal%subcellSize
            PRINT *, "       thisOctal%child(iIndex)%subcellSize = ", thisOctal%child(iIndex)%subcellSize
          END IF

          ! check the child's parent pointer
          IF ( .NOT. ASSOCIATED(thisOctal%child(iIndex)%parent,thisOctal) .AND. (thisOctal%nDepth > 1)) THEN
            PRINT *, "Error: In checkAMRgridPrivate, child has wrong %parent"
            CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
          END IF

          ! call recursively on this child
          thisOctalPointer => thisOctal
          CALL checkAMRgridPrivate(grid,                              &
                                   thisOctal=thisOctal%child(iIndex), &
                                   thisDepth=thisDepth+1,             &
                                   thisParent=thisOctalPointer,       &
                                   thisParentSubcell=iSubcell,        &
                                   nOctals=nOctals)
                                   
        END DO                           
                                   
      END IF ! ( thisOctal%nChildren > 0 )
      
    END SUBROUTINE checkAMRgridPrivate
    
    SUBROUTINE printErrorPrivate(grid,thisOctal,thisDepth,nOctals)  
        
      TYPE(gridType), INTENT(IN) :: grid
      TYPE(OCTAL), INTENT(IN) :: thisOctal
      INTEGER, INTENT(IN) :: thisDepth
      INTEGER, INTENT(IN) :: nOctals

      PRINT *, "  thisOctal%nDepth = ", thisOctal%nDepth
      PRINT *, "  thisDepth = ", thisDepth
      PRINT *, "  grid%maxDepth = ", grid%maxDepth
      PRINT *, "  grid%halfSmallestSubcell = ", grid%halfSmallestSubcell
      PRINT *, "  thisOctal%subcellSize = ", thisOctal%subcellSize
      PRINT *, "  thisOctal%nChildren = ", thisOctal%nChildren
      PRINT *, "  thisOctal%maxChildren = ", thisOctal%maxChildren
      PRINT *, "  thisOctal%twoD = ", thisOctal%twoD
      PRINT *, "  thisOctal%threeD = ", thisOctal%threeD
      PRINT *, "  thisOctal%hasChild = ", thisOctal%hasChild
      PRINT *, " ASSOCIATED(thisOctal%child) = ", ASSOCIATED(thisOctal%child) 
      IF ( ASSOCIATED(thisOctal%child) ) THEN
        PRINT *, "  SIZE(thisOctal%child) = ", SIZE(thisOctal%child)
      END IF
      PRINT *, "  grid%nOctals = ", grid%nOctals
      PRINT *, "  nOctals (counted so far) = ", nOctals
      !STOP
      PRINT *, "Entering infinite loop..."
      DO ; END DO
      
    END SUBROUTINE printErrorPrivate       
    
  END SUBROUTINE checkAMRgrid

  recursive subroutine unrefineThinCells(thisOctal, grid, ilambda, converged)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: ilambda
    real(double) :: kappaAbs, kappaSca, tau
    integer :: subcell, i
    logical :: unrefine, converged
    kappaAbs =0.d0; kappaSca = 0.d0
    unrefine = .true.

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unrefineThinCells(child, grid, ilambda, converged)
                exit
             end if
          end do
       else
          if (.not.ASSOCIATED(thisOctal%dustTypeFraction)) then
             write(*,*) "unalloc dusttypefraction!!"
          endif
          call returnKappa(grid, thisOctal, subcell, ilambda, kappaAbs=kappaAbs,kappaSca=kappaSca)
          tau = thisOctal%subcellSize*(kappaAbs+kappaSca)
          if (tau > 1.e-10) unrefine = .false.
       endif
    enddo

    if ((thisOctal%nChildren == 0).and.unrefine.and.converged) then
       call deleteChild(thisOctal%parent, thisOctal%parentSubcell, adjustParent = .true., &
            grid = grid, adjustGridInfo = .true.)
       converged = .false.
    endif

  end subroutine unrefineThinCells

  SUBROUTINE shrinkChildArray(parent, childrenToDelete, adjustParent )
    ! removes children from an octal.
    ! you probably don't want to call this directly - use the 'deleteChild' wrapper instead.
    ! NB this subroutine doesn't update grid%nOctals etc.

    IMPLICIT NONE
    TYPE(octal), TARGET, INTENT(INOUT) :: parent ! the parent octal 
    LOGICAL, INTENT(IN), DIMENSION(parent%maxChildren) :: childrenToDelete
      ! mask defining which children to get rid of 
      ! NB childrenToDelete does not map to the index number in the 
      !   %child array, it is the "real" number of the children (the 
      !   number of the subcell that was refined when the child was
      !   created).
    LOGICAL, INTENT(IN) :: adjustParent 
      ! whether the physical parameters stored in the parent's subcells
      !   should be filled with data derived from the children being deleted.
    
    TYPE(wrapperArray) :: tempChildStorage  
      ! holder for remaining children, while we shrink the %child array 
                                       
    LOGICAL, DIMENSION(parent%maxChildren) :: checkMask 
      ! used for testing the validity of the 'childrenToDelete' input
      
    INTEGER, DIMENSION(SIZE(parent%indexChild)) :: temporaryIndexChild
      ! used for assembling a valid indexChild array, to be used after
      !   the children have been deleted
      
    LOGICAL, DIMENSION(parent%maxChildren) :: deleteMask
      ! which of the elements of %child will be deleted
        
    INTEGER :: iChild ! loop counter
    INTEGER :: nChildren ! number of children the parent octal has
    INTEGER :: error
    INTEGER :: nChildrenToDelete
    TYPE(octal), POINTER :: thisChild ! convenient alias to current child 
    LOGICAL :: deleteALLchildren ! are we getting rid of ALL the children
    INTEGER :: nChildrenStay ! how many children will be left when done
    INTEGER :: insertLocation ! the next location to use in tempChildStorage
    
    NULLIFY(thisChild)
    temporaryIndexChild = -999
    
    ! setup some useful accounting variables
    nChildren = parent%nChildren
    nChildrenToDelete = COUNT(childrenToDelete)
    nChildrenStay = ( nChildren - nChildrenToDelete )
    deleteALLchildren = ( nChildrenStay == 0 )
    
    ! some safety checks
    error = 0
    
    IF ( nChildrenStay < 0 ) error = -1

    ! the following lines check that all the children to be deleted have 
    !   their %hasChild flag set.
    checkMask = childrenToDelete .AND. parent%hasChild(1:SIZE(childrenToDelete))
    checkMask = checkMask .NEQV. childrenToDelete ! (exclusive OR operation)
    IF ( ANY(checkMask) ) error = -2

    IF (error /= 0) THEN
      PRINT *, "In shrinkChildArray, attempting to delete a "
      PRINT *, "child that doesn't exist."
      PRINT *, error, childrenToDelete
      write(*,*) "nchildren ",parent%nChildren
      write(*,*) "haschild ",parent%hasChild(1:SIZE(childrenToDelete))
      write(*,*) "mask ",checkmask
      write(*,*) " and ",childrenToDelete(1:SIZE(childrenToDelete)) .AND. parent%hasChild(1:SIZE(childrenToDelete))
      write(*,*) "xor ",childrentodelete.neqv.(childrenToDelete .AND. parent%hasChild(1:SIZE(childrenToDelete)))
      write(*,*) nChildrentodelete,nchildrenstay,deleteallchildren

      STOP
    END IF
 
    ! we transform 'childrenToDelete' which is of SIZE(1:maxChildren) into
    !   an array of SIZE(1:nChildren), so that it matches the %child
    !   array
    deleteMask(:) = .FALSE.
    FORALL ( iChild = 1:parent%nChildren ) &
      deleteMask(iChild) = childrenToDelete(parent%indexChild(iChild))
    
    IF ( .NOT. deleteALLchildren ) THEN
      ! we need to allocate some temporary storage for the children that 
      !   are not going to be deleted so that we can reposition them
      !   in the %child array.
      ALLOCATE(tempChildStorage%wrappers(nChildrenStay), STAT=error)
      IF ( error /= 0 ) THEN
        PRINT *, 'Panic: allocation failed in shrinkChildArray. (A)'
        STOP
      END IF
      
      FORALL ( iChild = 1:nChildrenStay ) &
        tempChildStorage%wrappers(iChild)%inUse = .FALSE.
        
      insertlocation = 1
        
    END IF  
      
    DO iChild = 1, parent%nChildren
      thisChild => parent%child(iChild)

      IF ( deleteMask(iChild) ) THEN
        ! we want to delete this octal from the %child array

        CALL deleteOctal(thisChild, deleteChildren=.TRUE.,     &
                         adjustParent=adjustParent, adjustMem=.true. )

     ELSE 
        ! we do not want to delete this child. 
        ! instead we move it to the temporary storage array.

        ALLOCATE(tempChildStorage%wrappers(insertlocation)%content, STAT=error)
        IF ( error /= 0 ) THEN
          PRINT *, 'Panic: allocation failed in shrinkChildArray. (B)'
          STOP
        END IF           
        tempChildStorage%wrappers(insertlocation)%inUse = .TRUE.
               
        CALL deleteOctreeBranch(thisOctal=thisChild,                           &
               onlyChildren=.FALSE.,                                           &
               deletedBranch=tempChildStorage%wrappers(insertlocation)%content,&
               adjustParent=.FALSE.)
               
        temporaryIndexChild(insertLocation) = parent%indexChild(iChild)
        
        insertlocation = insertlocation + 1

      END IF
      
    END DO
   
    ! all the unwanted children have now been deleted.
    ! we now get rid of the old %child array.
    DEALLOCATE(parent%child)
    NULLIFY(parent%child) 
    
    ! if there are any remaining children, we need to create a new
    !   %child array, and copy them back there.
    IF ( .NOT. deleteALLchildren ) THEN    
      
      ALLOCATE(parent%child(nChildrenStay), STAT=error)
      IF ( error /= 0 ) THEN
        PRINT *, 'Panic: allocation failed in shrinkChildArray. (C)'
        STOP
      END IF

      DO iChild = 1, nChildrenStay, 1
        
        CALL insertOctreeBranch(parent%child(iChild),               &
               branch=tempChildStorage%wrappers(iChild)%content,    &
               onlyChildren=.FALSE.)                              
               
        parent%child(iChild)%parent => parent
               
        DEALLOCATE(tempChildStorage%wrappers(iChild)%content)
        NULLIFY(tempChildStorage%wrappers(iChild)%content)
        tempChildStorage%wrappers(iChild)%inUse = .FALSE.

      END DO
      
      ! can now clean up the temporary storage
      DEALLOCATE(tempChildStorage%wrappers)
      NULLIFY(tempChildStorage%wrappers)

    END IF

    ! some bookkeeping
    parent%nChildren = nChildrenStay
    parent%indexChild(:) = temporaryIndexChild(:)

    parent%hasChild( 1:SIZE(childrenToDelete) ) =  &
         parent%hasChild(1:SIZE(childrenToDelete)) .NEQV. childrenToDelete

  END SUBROUTINE shrinkChildArray
  
  SUBROUTINE deleteOctreeBranch(thisOctal,onlyChildren,deletedBranch,&
                                adjustParent, grid, adjustGridInfo)
    ! recursively deletes an octal's contents (and any children it has). 
    ! if 'deletedBranch' pointer is supplied, the branch is copied there as it
    !   is deleted. the %parent variable is NOT set here - you must do it 
    !   elsewhere.
    ! if onlyChildren is set, the octal's variables are not changed, only
    !   its children.
    ! *** WARNINGS:                                                   ***
    ! *** This does not delete the top octal itself -                 ***
    ! ***   you must do this yourself after you call this subroutine! ***

    ! HOWTO use this subroutine:
    ! 1. to "unrefine" part of the grid: call without a 'deletedBranch'
    !    and the children of 'thisOctal' will all be removed.
    ! 2. ...

    IMPLICIT NONE

    TYPE(octal), TARGET, INTENT(INOUT) :: thisOctal ! top of branch to be deleted
    LOGICAL, INTENT(IN)            :: onlyChildren ! only delete this octal's *children*
    TYPE(octal), INTENT(INOUT), TARGET, OPTIONAL :: deletedBranch ! optional copy of deleted branch
    LOGICAL, INTENT(IN) :: adjustParent 
      ! whether the physical parameters stored in the parent's subcells
      !   should be filled with data derived from the children being deleted.
    TYPE(gridtype), INTENT(INOUT), OPTIONAL :: grid 
    LOGICAL, INTENT(IN), OPTIONAL :: adjustGridInfo
      ! whether these variables should be updated: 
      !   grid%nOctals, grid%maxDepth, grid%halfSmallestSubcell    
      
    INTEGER :: iChild ! loop counter
    
    IF ( PRESENT(deletedBranch) ) THEN

       IF ( PRESENT(adjustGridInfo) ) THEN
         IF ( adjustGridInfo ) THEN
           PRINT *, "Sorry, adjustGridInfo is not supported with"
           PRINT *, "  a deletedBranch in deleteOctreeBranch"
           STOP
         END IF
       END IF

       IF (.NOT. onlyChildren) THEN
          ! need to copy variables from thisOctal to deletedBranch
          CALL copyOctalComponents(source=thisOctal, &
                                   dest=deletedBranch)
          NULLIFY(deletedBranch%parent)
          NULLIFY(deletedBranch%child)
       END IF
       
       ! now need to copy any children from thisOctal to deletedBranch
       IF ( thisOctal%nChildren > 0 ) THEN 
          deletedBranch%child => thisOctal%child
          
          DO iChild = 1, thisOctal%nChildren, 1
             IF ( onlyChildren ) THEN
                NULLIFY(deletedBranch%child(iChild)%parent)
             ELSE
                deletedBranch%child(iChild)%parent => deletedBranch
             END IF
          END DO
          
       END IF ! octal has children
       
       deletedBranch%nChildren = thisOctal%nChildren
       deletedBranch%indexChild = thisOctal%indexChild
       deletedBranch%hasChild = thisOctal%hasChild
       deletedBranch%maxChildren = thisOctal%maxChildren
       
       CALL deleteOctal(thisOctal, deleteChildren=.FALSE.,    &
                        adjustParent=adjustParent, grid=grid, & 
                        adjustGridInfo=.FALSE. )
       
    ELSE ! no deletedBranch
       
       ! can delete everything
       CALL deleteOctal(thisOctal, deleteChildren=.TRUE.,     &
                        adjustParent=adjustParent, grid=grid, & 
                        adjustGridInfo=adjustGridInfo )
       
    END IF ! PRESENT(deletedBranch)
    
    thisOctal%nChildren = 0
    thisOctal%indexChild(:) = -999
    thisOctal%hasChild(:) = .FALSE.
    
  END SUBROUTINE deleteOctreeBranch
  
  SUBROUTINE deleteOctal(thisOctal, deleteChildren,          &
                         adjustParent, grid, adjustGridInfo, adjustMem, &
                         newMaxDepth )
    ! deallocates the variables in an octal.
    ! optionally deallocates the children of the octal.
    use memory_mod, only : octalMemory, globalMemoryFootprint
    TYPE(octal), INTENT(INOUT) :: thisOctal
    LOGICAL, INTENT(IN) :: deleteChildren 
    LOGICAL, INTENT(IN) :: adjustParent 
      ! whether the physical parameters stored in the parent's subcells
      !   should be filled with data derived from the children being deleted.
    TYPE(gridtype), INTENT(INOUT), OPTIONAL :: grid 
    LOGICAL, INTENT(IN), OPTIONAL :: adjustGridInfo
    LOGICAL, INTENT(IN), OPTIONAL :: adjustMem
      ! whether these variables should be updated: 
      !   grid%nOctals, grid%maxDepth, grid%halfSmallestSubcell
    LOGICAL, INTENT(OUT), OPTIONAL :: newMaxDepth ! true if grid depth has changed
      
    LOGICAL :: doAdjustGridInfo
    LOGICAL :: doAdjustMem

    INTEGER :: maxDeletionDepth
      ! used for tracking the depth in the tree that has been altered   
    INTEGER, PARAMETER :: hugeInt = HUGE(hugeInt)

    doAdjustGridInfo = .FALSE.
    doAdjustMem = .false.
    
    IF ( PRESENT(adjustGridInfo) ) THEN
      IF ( adjustGridInfo .AND. (.NOT. (PRESENT(grid))) ) THEN
        PRINT *, "Panic: you must supply a grid if you"
        PRINT *, "       want to 'adjustGridInfo' in deleteOctal"
        STOP
      END IF
      doAdjustGridInfo = adjustGridInfo
    END IF

    if ( present(adjustMem)) doAdjustMem = .true.
    
    maxDeletionDepth = -99

    CALL deleteOctalPrivate(thisOctal, deleteChildren,        &
                            adjustParent,                     &
                            adjustGridInfo=doAdjustGridInfo,  &
                            adjustMem = doAdjustMem,          &
                            grid=grid, maxDeletionDepth=maxDeletionDepth )

    IF ( doAdjustGridInfo ) THEN
      ! we should see whether the maximum depth of the grid has shrunk

      IF ( .NOT. ( maxDeletionDepth < grid%maxDepth ) ) THEN
        CALL updateMaxDepth(grid,searchLimit=hugeInt, &
                            changeMade=newMaxDepth)
        IF ( newMaxDepth ) CALL setSmallestSubcell(grid)
      END IF
    END IF

  CONTAINS

    RECURSIVE SUBROUTINE deleteOctalPrivate(thisOctal, deleteChildren,      &
                                            adjustParent, adjustGridInfo, adjustMem,  &
                                            grid, maxDeletionDepth )
    
      TYPE(octal), INTENT(INOUT), target :: thisOctal
      type(octal), pointer :: pOctal
      LOGICAL, INTENT(IN) :: deleteChildren 
      LOGICAL, INTENT(IN) :: adjustParent 
      LOGICAL, INTENT(IN) :: adjustGridInfo
      LOGICAL, INTENT(IN) :: adjustMem
      TYPE(gridtype), INTENT(INOUT), OPTIONAL :: grid 
      INTEGER, INTENT(INOUT), OPTIONAL :: maxDeletionDepth

      INTEGER :: iChild
      INTEGER :: error
      error = 0

      maxDeletionDepth = MAX( thisOctal%nDepth, maxDeletionDepth )

      IF (deleteChildren) THEN

        DO iChild = 1, thisOctal%nChildren, 1

          CALL deleteOctalPrivate(thisOctal%child(iChild),         &
                           deleteChildren=deleteChildren,          &
                           adjustParent=adjustParent, grid=grid,   &
                           adjustGridInfo=adjustGridInfo,          &
                           adjustMem=adjustMem,                    &
                           maxDeletionDepth=maxDeletionDepth )
        END DO
        IF (ASSOCIATED(thisOctal%child)) DEALLOCATE(thisOctal%child)
        IF ( error /= 0 ) CALL deallocationError(error,location=1) 
        NULLIFY(thisOctal%child)

      END IF ! (deleteChildren)

      IF ( adjustParent ) CALL updateParentFromChild(childOctal=thisOctal)

      IF ( adjustGridInfo ) grid%nOctals = grid%nOctals - 1

      if (adjustMem .or. adjustGridinfo) then
         pOctal => thisOctal
         globalMemoryFootprint = globalMemoryFootprint - octalMemory(pOctal)
      endif

      call deallocateOctalDynamicAttributes(thisOctal)

    END SUBROUTINE deleteOctalPrivate

    SUBROUTINE deallocationError(error,location)
      INTEGER, INTENT(IN) :: error, location
      PRINT *, "DEALLOCATE error in deleteOctal"
      PRINT *, "Error:", error, " Location:", location
      STOP
    END SUBROUTINE deallocationError
    
  END SUBROUTINE deleteOctal


    SUBROUTINE growChildArray(parent, nNewChildren, grid)
    ! adds storage space for new children to an octal.
    ! 
    ! some of the bookkeeping must be done by the routine that calls
    !   this one (e.g. making sure the new children are indexed).

    IMPLICIT NONE
    
    TYPE(octal), TARGET, INTENT(INOUT) :: parent ! the parent octal 
    INTEGER, INTENT(IN) :: nNewChildren ! number of children to add
    TYPE(GRIDTYPE), INTENT(INOUT) :: grid
    
    TYPE(wrapperArray):: tempChildStorage 
      ! holder for existing children, while we shuffle them around to 
      !   make room for new ones.
                                       
    INTEGER       :: iChild            ! loop counter
    INTEGER       :: nChildren         ! number of children the parent octal has
    INTEGER       :: error
   
    ! store the number of children that already exist
    nChildren = parent%nChildren

    IF ( nChildren == 0 ) THEN
      ! if there are no existing children, we can just allocate
      ! the 'child' array with the new size
      ALLOCATE(parent%child(nNewChildren), STAT=error)
      IF ( error /= 0 ) THEN
        PRINT *, 'Panic: allocation failed in growChildArray.(A)'
        STOP
      END IF

    ELSE ! there are existing children
     
      ! check that there is not a full quota of children
      IF ( (nChildren + nNewChildren) > parent%maxChildren ) THEN
        PRINT *, 'Panic: in growChildArray, attempted to have too ',&
                 '       many children: ',(nChildren + nNewChildren) 
        STOP
      ENDIF
      
      ! if there are existing children, we must enlarge the allocated   
      ! array. we need to use temporary octals[1] and copy the  
      ! existing children into them[2]; then increase the 'child' array size 
      ! [3]; then copy the children back in[4]. 

      ! [1]
      ALLOCATE(tempChildStorage%wrappers(nChildren), STAT=error)
      IF ( error /= 0 ) THEN
        PRINT *, 'Panic: allocation failed in growChildArray. (B)'
        STOP
      END IF     

      ! [2]
      DO iChild = 1, nChildren, 1
        
        ALLOCATE(tempChildStorage%wrappers(iChild)%content, STAT=error)
        IF ( error /= 0 ) THEN
          PRINT *, 'Panic: allocation failed in growChildArray. (C)'
          STOP
        END IF           
        tempChildStorage%wrappers(iChild)%inUse = .TRUE.
               
        CALL deleteOctreeBranch(parent%child(iChild),                   &
               onlyChildren=.FALSE.,                                    &
               deletedBranch=tempChildStorage%wrappers(iChild)%content, &
               adjustParent=.FALSE., grid=grid, adjustGridInfo=.FALSE. )
      END DO

      ! [3]
      IF ( ASSOCIATED(parent%child) ) THEN
        DEALLOCATE(parent%child)
        NULLIFY(parent%child)
        ALLOCATE(parent%child( nChildren + nNewChildren ), STAT=error)
        IF ( error /= 0 ) THEN
          PRINT *, 'Panic: allocation failed in growChildArray. (D)'
          STOP
        END IF
      ELSE
        PRINT *, "Error in growChildArray:"
        PRINT *, "  parent%child not associated."
        STOP
      END IF

      ! [4]
      DO iChild = 1, nChildren, 1
        
        CALL insertOctreeBranch(parent%child(iChild),               &
               branch=tempChildStorage%wrappers(iChild)%content,    &
               onlyChildren=.FALSE.)

        parent%child(iChild)%parent => parent       
               
        if (associated(tempChildStorage%wrappers(iChild)%content)) then
           DEALLOCATE(tempChildStorage%wrappers(iChild)%content, STAT=error)
           if (error /= 0) then
              write(*,*) "error",error
              write(*,*) ichild
              write(*,*) tempChildStorage%wrappers(ichild)%content%ndepth
              parent%child(iChild)%rho = 1.d30
!              do;enddo
           endif
        endif

        NULLIFY(tempChildStorage%wrappers(iChild)%content)
        tempChildStorage%wrappers(iChild)%inUse = .FALSE.

      END DO
      
      ! can now clean up the temporary storage
      DEALLOCATE(tempChildStorage%wrappers)
      NULLIFY(tempChildStorage%wrappers)
      
    END IF ! ( nChildren == 0 )

  END SUBROUTINE growChildArray

  SUBROUTINE insertOctreeBranch(thisOctal,branch,onlyChildren)
    ! adds one octree into another.
    ! use with care, it's not sensible to insert a tree anywhere - this 
    !   subroutine is meant to be used to replace a tree after it
    !   was temporarily removed from the same location.
    ! the %parent variable of thisOctal is not set - you must do this 
    !   elsewhere
    ! if onlyChildren is set, the octal's variables are not changed, only
    !   its children.
    ! the contents of 'branch' are deleted. you should delete 'branch'
    !   itself after calling this subroutine.
    
    ! *** WARNINGS:                                                   ***
    ! *** This does not delete the top octal of the branch itself -   ***
    ! ***   you must do this yourself after you call this subroutine! ***
    ! *** If the octree being changed is part of a 'grid', you'll     ***
    ! ***   have to fix grid%nOctals etc.                             ***

    IMPLICIT NONE

    TYPE(octal), TARGET, INTENT(INOUT) :: thisOctal ! octal where branch is to
!                                                    !   be inserted
    TYPE(octal), INTENT(INOUT) :: branch ! branch being inserted
    LOGICAL, INTENT(IN) :: onlyChildren ! only insert the *children* on the
                                        !   branch, leaving the other variables
                                        !   of thisOctal unchanged
    INTEGER :: iChild ! loop counter
    
    IF (ASSOCIATED(thisOctal%child)) THEN
      WRITE(*,*) "Error in insertOctreeBranch, attempt to overwrite existing children"
      STOP
    END IF

    IF (.NOT. onlyChildren) CALL copyOctalComponents(source=branch,dest=thisOctal)

    IF (ASSOCIATED(branch%child)) THEN
      thisOctal%child => branch%child

      DO iChild = 1, SIZE(thisOctal%child), 1
        thisOctal%child(iChild)%parent => thisOctal
      END DO
      thisOctal%nChildren = branch%nChildren
      thisOctal%indexChild = branch%indexChild
      thisOctal%hasChild = branch%hasChild
      thisOctal%maxChildren = branch%maxChildren
      branch%nChildren = 0 
      branch%hasChild(:) = .FALSE. 
    END IF

    CALL deleteOctal(branch, deleteChildren=.FALSE., adjustParent=.FALSE.)
        
  END SUBROUTINE insertOctreeBranch
  
  SUBROUTINE copyOctalComponents(source,dest)
    ! copy the components within an octal variable to a new octal variable
    !
    ! WARNING: this does not change the parent and child variables - you
    !   must update those yourself elsewhere.
 
    TYPE(octal), INTENT(IN) :: source
    TYPE(octal), INTENT(INOUT) :: dest

    ! first make sure that 'dest' is empty
    IF ( ASSOCIATED(dest%child) ) THEN
      PRINT *, "Problem in copyOctalComponents:"
      PRINT *, "destination seems to have some children"
      do;enddo
    END IF
    CALL deleteOctal(dest, deleteChildren=.FALSE., adjustParent=.FALSE. )

    dest%nDepth =  source%nDepth
    dest%mpiThread =  source%mpiThread
    dest%nChildren = source%nChildren
    dest%indexChild = source%indexChild
    dest%threeD =   source%threeD 
    dest%twod = source%twoD  
    dest%oneD = source%oneD   

    dest%cylindrical = source%cylindrical
    dest%splitAzimuthally = source%splitAzimuthally
    dest%maxChildren =  source%maxChildren
    dest%hasChild = source%hasChild
    dest%centre =  source%centre
    dest%rho =  source%rho
    dest%phi =  source%phi
    dest%dphi = source%dphi

    dest%phimin = source%phimin
    dest%phimax = source%phimax

    dest%r =   source%r
    dest%velocity = source%velocity
    dest%temperature = source%temperature
    dest%inFlow = source%inFlow
    dest%label =   source%label
    dest%parentSubcell =   source%parentSubcell
    dest%subcellSize =  source%subcellSize
    dest%gasOpacity =  source%gasOpacity 

    dest%xMax = source%xMax
    dest%yMax = source%yMax
    dest%zMax = source%zMax


    dest%xMin = source%xMin
    dest%yMin = source%yMin
    dest%zMin = source%zMin

    dest%inStar = source%inStar

    call copyAttribute(dest%nCrossings, source%nCrossings)
    call copyAttribute(dest%chiLine, source%chiLine)
    call copyAttribute(dest%etaLine, source%etaLine)
    call copyAttribute(dest%etaCont, source%etaCont)
    call copyAttribute(dest%biasCont3d, source%biasCont3d)
    call copyAttribute(dest%biasLine3d, source%biasLine3d)
    call copyAttribute(dest%distanceGrid, source%distanceGrid)
    call copyAttribute(dest%cornerVelocity, source%cornerVelocity)
    call copyAttribute(dest%cornerRho, source%cornerRho)

    call copyAttribute(dest%probDistLine, source%probDistLine)
    call copyAttribute(dest%probDistCont,  source%probDistCont)
    call copyAttribute(dest%Ne,  source%Ne)
    call copyAttribute(dest%Ntot,  source%Ntot)
    call copyAttribute(dest%NH,  source%NH)
    call copyAttribute(dest%NHI,  source%NHI)
    call copyAttribute(dest%boundaryCondition,  source%boundaryCondition)
    call copyAttribute(dest%boundaryCell,  source%boundaryCell)
    call copyAttribute(dest%boundaryPartner,  source%boundaryPartner)
    call copyAttribute(dest%GravboundaryPartner,  source%GravboundaryPartner)
    call copyAttribute(dest%NHII,  source%NHII)
    call copyAttribute(dest%NHeI,  source%NHeI)
    call copyAttribute(dest%NHeII,  source%NHeII)
    call copyAttribute(dest%Hheating,  source%Hheating)
    call copyAttribute(dest%tDust,  source%tDust)
    call copyAttribute(dest%Heheating,  source%Heheating)
    call copyAttribute(dest%changed,  source%changed)

    call copyAttribute(dest%fixedTemperature,  source%fixedTemperature)


    call copyAttribute(dest%scatteredIntensity, source%scatteredIntensity)
    call copyAttribute(dest%meanIntensity, source%meanIntensity)
    call copyAttribute(dest%temperaturedust, source%temperatureDust)
    call copyAttribute(dest%temperaturegas, source%temperaturegas)
    call copyAttribute(dest%dustType, source%dustType)
    call copyAttribute(dest%oldFrac, source%oldFrac)
    call copyAttribute(dest%diffusionApprox, source%diffusionApprox)
    call copyAttribute(dest%diffusionCoeff, source%diffusionCoeff)
    call copyAttribute(dest%undersampled, source%undersampled)
    call copyAttribute(dest%nDiffusion, source%nDiffusion)
    call copyAttribute(dest%nDirectPhotons, source%nDirectPhotons)
    call copyAttribute(dest%oldtemperature, source%oldtemperature)
    call copyAttribute(dest%eDens, source%eDens)
    call copyAttribute(dest%oldeDens, source%oldeDens)
    call copyAttribute(dest%kappaRoss, source%kappaRoss)
 
    call copyAttribute(dest%nh2, source%nh2)
    call copyAttribute(dest%microturb, source%microturb)
    call copyAttribute(dest%molmicroturb, source%molmicroturb)

    call copyAttribute(dest%rho_i_minus_1, source%rho_i_minus_1)
    call copyAttribute(dest%rho_i_plus_1, source%rho_i_plus_1)

    call copyAttribute(dest%rhorv_i_minus_1, source%rhorv_i_minus_1)
    call copyAttribute(dest%rhorv_i_plus_1, source%rhorv_i_plus_1)


    call copyAttribute(dest%x_i_minus_1, source%x_i_minus_1)
    call copyAttribute(dest%x_i_minus_2, source%x_i_minus_2)

    call copyAttribute(dest%x_i, source%x_i)
    call copyAttribute(dest%x_i_plus_1, source%x_i_plus_1)
    call copyAttribute(dest%x_i_minus_2, source%x_i_minus_2)

    call copyAttribute(dest%q_i_minus_2, source%q_i_minus_2)
    call copyAttribute(dest%q_i_minus_1, source%q_i_minus_1)
    call copyAttribute(dest%q_i, source%q_i)
    call copyAttribute(dest%fViscosity, source%fViscosity)
    call copyAttribute(dest%q_i_plus_1, source%q_i_plus_1)

    call copyAttribute(dest%flux_i_minus_1, source%flux_i_minus_1)
    call copyAttribute(dest%flux_i, source%flux_i)
    call copyAttribute(dest%flux_i_plus_1, source%flux_i_plus_1)

    call copyAttribute(dest%pressure_i_minus_1, source%pressure_i_minus_1)
    call copyAttribute(dest%pressure_i, source%pressure_i)
    call copyAttribute(dest%phi_i, source%phi_i)
    call copyAttribute(dest%phi_gas, source%phi_gas)
    call copyAttribute(dest%phi_stars, source%phi_stars)
    call copyAttribute(dest%phi_i_plus_1, source%phi_i_plus_1)
    call copyAttribute(dest%phi_i_minus_1, source%phi_i_minus_1)
    call copyAttribute(dest%pressure_i_plus_1, source%pressure_i_plus_1)




    call copyAttribute(dest%q_amr_i_minus_1, source%q_amr_i_minus_1)
    call copyAttribute(dest%q_amr_i_plus_1, source%q_amr_i_plus_1)
    call copyAttribute(dest%u_interface, source%u_interface)
    call copyAttribute(dest%u_amr_interface, source%u_amr_interface)
    call copyAttribute(dest%u_amr_interface_i_plus_1, source%u_amr_interface_i_plus_1)
    call copyAttribute(dest%flux_amr_i, source%flux_amr_i)
    call copyAttribute(dest%flux_amr_i_plus_1, source%flux_amr_i_plus_1)
    call copyAttribute(dest%flux_amr_i_minus_1, source%flux_amr_i_minus_1)
    call copyAttribute(dest%phiLimit_amr, source%phiLimit_amr)


    call copyAttribute(dest%u_interface, source%u_interface)
    call copyAttribute(dest%u_i, source%u_i)
    call copyAttribute(dest%u_i_plus_1, source%u_i_plus_1)
    call copyAttribute(dest%u_i_minus_1, source%u_i_minus_1)
    call copyAttribute(dest%rLimit, source%rLimit)
    call copyAttribute(dest%phiLimit, source%phiLimit)
    call copyAttribute(dest%ghostCell, source%ghostCell)
    call copyAttribute(dest%corner, source%corner)
    call copyAttribute(dest%edgeCell, source%edgeCell)
    call copyAttribute(dest%feederCell, source%feederCell)
    call copyAttribute(dest%energy, source%energy)
    call copyAttribute(dest%divV, source%divV)
    call copyAttribute(dest%rhou, source%rhou)
    call copyAttribute(dest%rhov, source%rhov)
    call copyAttribute(dest%rhow, source%rhow)
    call copyAttribute(dest%rhoe, source%rhoe)
    call copyAttribute(dest%rhoeLastTime, source%rhoeLastTime)
    call copyAttribute(dest%qViscosity, source%qViscosity)
    call copyAttribute(dest%refinedLastTime, source%refinedLastTime)

    call copyAttribute(dest%photoIonCoeff, source%photoIonCoeff)
    call copyAttribute(dest%sourceContribution, source%sourceContribution)
    call copyAttribute(dest%diffuseContribution, source%diffuseContribution)
    call copyAttribute(dest%normSourceContribution, source%normSourceContribution)
    call copyAttribute(dest%molecularLevel, source%molecularLevel)

    call copyAttribute(dest%molAbundance, source%molAbundance)

    call copyAttribute(dest%atomAbundance, source%atomAbundance)

    call copyAttribute(dest%atomLevel, source%atomLevel)

    call copyAttribute(dest%newatomLevel, source%newatomLevel)
    call copyAttribute(dest%dust_T, source%dust_T)
    call copyAttribute(dest%ionFrac, source%ionFrac)
    call copyAttribute(dest%columnRho, source%columnRho)

#ifdef PDR
    call copyAttribute(dest%AV, source%AV)
    call copyAttribute(dest%thisColRho, source%thisColRho)
    call copyAttribute(dest%relch, source%relch)

    call copyAttribute(dest%abundance, source%abundance)
    call copyAttribute(dest%radsurface, source%radsurface)    
    call copyAttribute(dest%Tlast, source%Tlast)    
    call copyAttribute(dest%Tmin, source%Tmin)    
    call copyAttribute(dest%Tmax, source%Tmax)    
    call copyAttribute(dest%Tlow, source%Tlow)    
    call copyAttribute(dest%Thigh, source%Thigh)    
    call copyAttribute(dest%Tminarray, source%Tminarray)    
    call copyAttribute(dest%Tmaxarray, source%Tmaxarray)    
    call copyAttribute(dest%UV, source%UV)
    call copyAttribute(dest%converged, source%converged)
    call copyAttribute(dest%level_converged, source%level_converged)
    call copyAttribute(dest%biChop, source%biChop)
    call copyAttribute(dest%expanded, source%expanded)
    call copyAttribute(dest%lastChange, source%lastChange)
    call copyAttribute(dest%tPrev, source%tPrev)


    call copyAttribute(dest%coolingRate, source%coolingRate)
    call copyAttribute(dest%heatingRate, source%heatingRate)

    call copyAttribute(dest%ciiLine, source%ciiLine)
    call copyAttribute(dest%ciiTransition, source%ciiTransition)
    call copyAttribute(dest%ciLine, source%ciLine)
    call copyAttribute(dest%ciTransition, source%ciTransition)
    call copyAttribute(dest%oiLine, source%oiLine)
    call copyAttribute(dest%oiTransition, source%oiTransition)
    call copyAttribute(dest%c12oLine, source%c12oLine)
    call copyAttribute(dest%c12oTransition, source%c12oTransition)


    call copyAttribute(dest%cii_pop, source%cii_pop)
    call copyAttribute(dest%ci_pop, source%ci_pop)
    call copyAttribute(dest%oi_pop, source%oi_pop)
    call copyAttribute(dest%c12o_pop, source%c12o_pop)
#endif
    call copyAttribute(dest%gamma, source%gamma)
    call copyAttribute(dest%iEquationOfState, source%iEquationOfState)

    call copyAttribute(dest%iAnalyticalVelocity, source%iAnalyticalVelocity)

    call copyAttribute(dest%uDens, source%udens)
    call copyAttribute(dest%aDot, source%aDot)
    call copyAttribute(dest%distancegridaDot, source%distancegridaDot)
    call copyAttribute(dest%distanceGridPhotonFromGas, source%distanceGridPhotonFromGas)
    call copyAttribute(dest%distanceGridPhotonFromSource, source%distanceGridPhotonFromSource)
    call copyAttribute(dest%photonEnergyDensityFromGas, source%photonEnergyDensityFromGas)
    call copyAttribute(dest%photonEnergyDensityFromSource, source%photonEnergyDensityFromSource)

    call copyAttribute(dest%photonEnergyDensity, source%photonEnergyDensity)
    call copyAttribute(dest%oldphotonEnergyDensity, source%oldphotonEnergyDensity)

    call copyAttribute(dest%radiationMomentum, source%radiationMomentum)

    call copyAttribute(dest%kappaTimesFlux, source%kappaTimesFlux)
    call copyAttribute(dest%UVvector, source%UVvector)
    call copyAttribute(dest%UVvectorplus, source%UVvectorplus)
    call copyAttribute(dest%UVvectorminus, source%UVvectorminus)


    IF (ASSOCIATED(source%mpiboundaryStorage)) THEN                   
      ALLOCATE(dest%mpiboundaryStorage( SIZE(source%mpiboundaryStorage,1),       &
                              SIZE(source%mpiboundaryStorage,2), &
                              SIZE(source%mpiboundaryStorage,3)))
      dest%mpiboundaryStorage = source%mpiboundaryStorage
    END IF  

    IF (ASSOCIATED(source%mpiCornerStorage)) THEN                   
      ALLOCATE(dest%mpiCornerStorage( SIZE(source%mpiCornerStorage,1),       &
                              SIZE(source%mpiCornerStorage,2), &
                              SIZE(source%mpiCornerStorage,3)))
      dest%mpiCornerStorage = source%mpiCornerStorage
    END IF  

    IF (ASSOCIATED(source%kappaAbs)) THEN                   
      ALLOCATE(dest%kappaAbs( SIZE(source%kappaAbs,1),       &
                              SIZE(source%kappaAbs,2)))
      dest%kappaAbs = source%kappaAbs
    END IF  

    IF (ASSOCIATED(source%kappaSca)) THEN
      ALLOCATE(dest%kappaSca( SIZE(source%kappaSca,1),       &
                              SIZE(source%kappaSca,2)))
      dest%kappaSca = source%kappaSca
    END IF  

    IF (ASSOCIATED(source%N)) THEN 
      ALLOCATE(dest%N( SIZE(source%N,1),              &
                       SIZE(source%N,2)))
      dest%N = source%N
    END IF  

    IF (ASSOCIATED(source%departCoeff)) THEN
      ALLOCATE(dest%departCoeff( SIZE(source%departCoeff,1),    &
                                 SIZE(source%departCoeff,2)))
      dest%departCoeff = source%departCoeff
    END IF  

    IF (ASSOCIATED(source%gas_particle_list)) THEN
      ALLOCATE(dest%gas_particle_list(SIZE(source%gas_particle_list)))
      dest%gas_particle_list = source%gas_particle_list
    END IF  

    IF (ASSOCIATED(source%dustTypeFraction)) THEN                   
      ALLOCATE(dest%dustTypeFraction( SIZE(source%dustTypeFraction,1), &
                                      SIZE(source%dustTypeFraction,2)))
      dest%dustTypeFraction = source%dustTypeFraction
    END IF  

    IF (ASSOCIATED(source%origdustTypeFraction)) THEN                   
      ALLOCATE(dest%origdustTypeFraction( SIZE(source%origdustTypeFraction,1), &
                                      SIZE(source%origdustTypeFraction,2)))
      dest%origdustTypeFraction = source%origdustTypeFraction
    END IF  

  END SUBROUTINE copyOctalComponents

    
  SUBROUTINE deleteChild(parent, childToDelete, adjustParent, &
                         grid, adjustGridInfo)
    ! removes one child from an octal.
    ! can be used to "unrefine" part of the grid.
    ! this is a simpler wrapper for the shrinkChildArray subroutine

    IMPLICIT NONE
    TYPE(octal), INTENT(INOUT) :: parent ! the parent octal 
    INTEGER :: childToDelete ! number of the child to delete
      ! NB childToDelete does not map to the index number in the 
      !   %child array, it is the "real" number of the child (the 
      !   number of the subcell that was refined when the child was
      !   created).
    LOGICAL, INTENT(IN) :: adjustParent 
      ! whether the physical parameters stored in the parent's subcells
      !   should be filled with data derived from the children being deleted.
    TYPE(gridtype), INTENT(INOUT) :: grid 
    LOGICAL, INTENT(IN) :: adjustGridInfo
      ! whether these variables should be updated: 
      !   grid%nOctals, grid%maxDepth, grid%halfSmallestSubcell

    LOGICAL, DIMENSION(parent%maxChildren) :: deletionMask
    INTEGER, PARAMETER :: hugeInt = HUGE(hugeInt)
    LOGICAL :: newMaxDepth ! true if grid depth has changed
    
    IF ( (childToDelete > parent%maxChildren) .OR.   &
         (childToDelete < 0)                      ) THEN
      PRINT *, "Invalid child number passed to deleteChild: ", childToDelete 
      STOP
    END IF
    
    deletionMask(:) = .FALSE.
    deletionMask(childToDelete) = .TRUE.
    
    CALL shrinkChildArray(parent=parent,                 &
                          childrenToDelete=deletionMask, &
                          adjustParent=adjustParent)

    IF ( adjustGridInfo ) THEN
      grid%nOctals = grid%nOctals - 1

      IF ( ( parent%nChildren == 0 ) .AND.              &
           ( grid%maxDepth == (parent%nDepth+1) ) ) THEN
        ! we should see whether the maximum depth of the grid has shrunk
        CALL updateMaxDepth(grid,searchLimit=hugeInt, &
                            changeMade=newMaxDepth)
        IF ( newMaxDepth ) CALL setSmallestSubcell(grid)
      END IF
    END IF  

    ! only for debugging - comment out later:
!    CALL checkAMRgrid(grid,checkNoctals=.FALSE.)
        
      
  END SUBROUTINE deleteChild

  SUBROUTINE updateParentFromChild(childOctal)
    ! uses the 4 or 8 subcells in an octal to compute the physical
    !   variables in the appropriate subcell of that octal's parent.
    TYPE(OCTAL), INTENT(INOUT) :: childOctal 
    
    TYPE(OCTAL), POINTER :: parentOctal
    INTEGER :: parentSubcell !, iSubcell
    INTEGER :: nVals, i
    REAL(double) :: nValsREAL

!    REAL(double) :: dv, oldMass, newMass, oldEnergy, newEnergy
 !   REAL(double) :: factor

    IF ( childOctal%nDepth == 1 ) THEN
      ! we're at the root of the tree
      PRINT *, "Warning: Attempted to update parent, but already at tree root."
      PRINT *, "         Ignoring this and continuing..."
      RETURN
    END IF
    
    parentOctal => childOctal%parent
    parentSubcell = childOctal%parentSubcell
    
    if (childOctal%oneD) then
       nVals = 2
    else if (childOctal%twoD) then
       nVals = 4
    else if (childOctal%threeD) then
       nVals = 8
    endif

    nValsREAL = REAL( nVals, KIND=double )

    if (associated(childOctal%boundaryCondition)) &
         parentOctal%boundaryCondition(parentSubcell) = childOctal%boundaryCondition(1)
    if (associated(childOctal%gamma)) &
         parentOctal%gamma(parentSubcell) = childOctal%gamma(1)
    if (associated(childOctal%iEquationOfState)) &
         parentOctal%iEquationofState(parentSubcell) = childOctal%iEquationofState(1)
    
    parentOctal%rho(parentSubcell) =                    &
    SUM(childOctal%rho(1:nVals)) / nValsREAL

    if (associated(parentOctal%nh)) parentOctal%nh(parentSubcell) =                    &
    SUM(childOctal%nh(1:nVals)) / nValsREAL

    if (associated(parentOctal%ne)) parentOctal%ne(parentSubcell) =                    &
    SUM(childOctal%ne(1:nVals)) / nValsREAL

    if (associated(parentOctal%rhoe)) parentOctal%rhoe(parentSubcell) =  &
    SUM(childOctal%rhoe(1:nVals)) / nValsREAL

    if (associated(parentOctal%rhou))     parentOctal%rhou(parentSubcell) = &
    SUM(childOctal%rhou(1:nVals)) / nValsREAL
 
    if (associated(parentOctal%rhov))     parentOctal%rhov(parentSubcell) = &
    SUM(childOctal%rhov(1:nVals)) / nValsREAL

    if (associated(parentOctal%rhow))     parentOctal%rhow(parentSubcell) = &
    SUM(childOctal%rhow(1:nVals)) / nValsREAL
   
    parentOctal%temperature(parentSubcell) = &
     real(SUM(childOctal%temperature(1:nVals)) / nValsREAL)

    if (associated(parentOctal%ionFrac)) then
       do i = 1, SIZE(parentOctal%ionFrac,2)
          parentOctal%ionFrac(parentSubcell,i) =            &
               SUM(childOctal%ionFrac(1:nVals,i)) / nValsREAL
       enddo
    endif
    
    if (associated(parentOctal%dustTypeFraction)) then
       do i = 1, SIZE(childOCtal%dustTypeFraction,2)
          parentOctal%dustTypeFraction(parentSubcell,i) =     &
               SUM(childOctal%dustTypeFraction(1:nVals,i)) / nValsREAL
       enddo
    endif

    if (associated(childOctal%mpiBoundaryStorage)) then
       if (.not.associated(parentOctal%mpiBoundaryStorage)) then
          allocate(parentOctal%mpiBoundaryStorage( &
               size(childOctal%mpiBoundaryStorage,1), &
               size(childOctal%mpiBoundaryStorage,2), &
               size(childOctal%mpiBoundaryStorage,3)))
          parentOctal%mpiBoundaryStorage = childOctal%mpiBoundaryStorage
       endif
    endif

  END SUBROUTINE updateParentFromChild
  
  SUBROUTINE updateMaxDepth(grid,searchLimit,changeMade)
    ! if octals have been deleted from the grid, we may have to 
    !   check whether the grid%maxDepth variable should be updated
  
    TYPE(gridtype), INTENT(INOUT)  :: grid
    INTEGER, INTENT(IN), OPTIONAL :: searchLimit 
      ! the depth at which to abandon the search.
      ! e.g. if octals have been deleted from levels 2 and 3, pass 
      !   "3" in here, and the search can be abandoned if an octal
      !   of depth 4 is encountered anywhere
    LOGICAL, INTENT(OUT), OPTIONAL :: changeMade
      ! true if grid%maxDepth is actually changed

    TYPE(octal), POINTER  :: thisOctal
    INTEGER :: maxDepthFound
    INTEGER :: oldMaxDepth
    LOGICAL :: searchLimitReached

    oldMaxDepth = grid%maxDepth
    thisOctal => grid%octreeRoot
    maxDepthFound = 0
    searchLimitReached = .FALSE.
    
    CALL updateMaxDepthPrivate(thisOctal,maxDepthFound,       &
                               searchLimitReached,searchLimit)
    IF ( searchLimitReached .OR. ( maxDepthFound == oldMaxDepth ) ) THEN
      IF ( PRESENT(changeMade) ) changeMade = .FALSE.
    ELSE
      IF ( PRESENT(changeMade) ) changeMade = .TRUE.
      grid%maxDepth = maxDepthFound
    END IF
    
    CONTAINS
    
      RECURSIVE SUBROUTINE updateMaxDepthPrivate(thisOctal,maxDepthFound,&
                                          searchLimitReached,searchLimit)
      
        TYPE(OCTAL), INTENT(IN) :: thisOctal 
        INTEGER, INTENT(INOUT) :: maxDepthFound
        LOGICAL, INTENT(INOUT) :: searchLimitReached
        INTEGER, INTENT(IN), OPTIONAL :: searchLimit 
        INTEGER :: iChild
        
        maxDepthFound = MAX( maxDepthFound, thisOctal%nDepth )
        
        IF ( thisOctal%nChildren > 0 ) THEN

          IF ( PRESENT(searchLimit) ) THEN
            IF ( (thisOctal%nDepth+1) > searchLimit ) THEN
              searchLimitReached = .TRUE.
              RETURN
            END IF
          END IF
          
          ! call this subroutine recursively on each of its children
          DO iChild = 1, thisOctal%nChildren, 1
            CALL updateMaxDepthPrivate(thisOctal%child(iChild),maxDepthFound,&
                                          searchLimitReached,searchLimit)
                                          
            IF ( searchLimitReached ) RETURN
          END DO
        END IF

      END SUBROUTINE updateMaxDepthPrivate

  END SUBROUTINE updateMaxDepth

  
  SUBROUTINE setSmallestSubcell(grid)
    ! calculates and stores the grid%halfSmallestSubcell value.
    ! grid%maxDepth must already be set correctly.
  
    TYPE(gridtype), INTENT(INOUT) :: grid 
    
    ! we actually store the value which is half the size of the 
    !   smallest subcell because this is more useful for later
    !   calculations.
    grid%halfSmallestSubcell = grid%octreeRoot%subcellSize / &
                                2.0_oc**REAL(grid%maxDepth,kind=oct)

  END SUBROUTINE setSmallestSubcell

  RECURSIVE SUBROUTINE setAllUnchanged(thisOctal)
    ! goes through an octree and sets all the %changed variables
    !   to .FALSE.
  
    TYPE(octal), INTENT(INOUT) :: thisOctal 
    INTEGER :: iChild
    
    thisOctal%changed(:) = .FALSE.
    
    DO iChild = 1, thisOctal%nChildren
      CALL setAllUnchanged(thisOctal%child(iChild))
    END DO

  END SUBROUTINE setAllUnchanged


  SUBROUTINE calcWindTestValues(thisOctal,subcell,grid) 
    ! calculates some of the variables for a spherical wind flow

    USE constants_mod

    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(INOUT) :: grid

    TYPE(vector) :: point

    TYPE(vector) :: starPosn
    TYPE(vector) :: pointVec
    TYPE(vector) :: pointVecNorm
    TYPE(Vector)      :: pointVecNormSingle

    REAL :: r, rStar
    REAL :: velocity, rho
    REAL, PARAMETER :: v0         = 100.e5
    REAL, PARAMETER :: vTerminal  = 2000.e5
    REAL, PARAMETER :: tEff       = 30.e3
    REAL, PARAMETER :: mDot       = real(1.e-7 * mSol * secsToYears)
!    REAL, PARAMETER :: mDot       = 1.e-2 * mSol * secsToYears 
    REAL, PARAMETER :: vBeta      = 1.0
    

    rStar = grid%rStar1
    starPosn = grid%starPos1

    point = subcellCentre(thisOctal,subcell)
    pointVec = (point - starPosn) 
    
    r = real(modulus( pointVec ) )
    pointVecNorm = pointVec 
    CALL normalize(pointVecNorm)

    IF (r > grid%rInner .AND. r < grid%rOuter) THEN

       ! calculate the velocity
       velocity = v0 + (vTerminal - v0) * (1. - rStar / r)**vBeta

       ! calculate the density
       rho = real(mDot / ( fourPi * (r*1.e10)**2 * velocity))

       ! store the data 
       thisOctal%inFlow(subcell) = .TRUE.
       thisOctal%temperature(subcell) = 0.8 * tEff
       pointVecNormSingle = pointVecNorm
       thisOctal%velocity(subcell) = (velocity / cSpeed) * pointVecNormSingle
       thisOctal%rho(subcell) = rho

    ELSE
       thisOctal%inFlow(subcell) = .FALSE.
       thisOctal%rho(subcell) = 1.e-20
       thisOctal%temperature(subcell) = 0.8 * tEff
       thisOctal%velocity(subcell) = vector(1.e-25,1.e-25,1.e-25)
    END IF
    
    IF (subcell == 8) CALL fillVelocityCorners(thisOctal,windTestVelocity)
      
  END SUBROUTINE calcWindTestValues

  
  TYPE(vector) FUNCTION windTestVelocity(point)

    use inputs_mod, only : rinner, router
    TYPE(vector), INTENT(IN) :: point
    
    TYPE(vector) :: starPosn
    TYPE(vector) :: pointVec
    TYPE(vector) :: pointVecNorm
    TYPE(Vector)      :: pointVecNormSingle

    REAL :: r, rStar
    REAL :: velocity
    REAL, PARAMETER :: v0         = 100.e5
    REAL, PARAMETER :: vTerminal  = 2000.e5
    REAL, PARAMETER :: vBeta      = 1.0
    
    rStar = real(rsol)
    starPosn = VECTOR(0.d0, 0.d0,0.d0)

    pointVec = (point - starPosn) 
    
    r = real(modulus( pointVec ) )
    pointVecNorm = pointVec 
    CALL normalize(pointVecNorm)
    pointVecNormSingle = pointVecNorm
    
    IF (r > rInner .AND. r < rOuter) THEN
       velocity = v0 + (vTerminal - v0) * (1. - rStar / r)**vBeta
       windTestVelocity = (velocity / cSpeed) * pointVecNormSingle
    ELSE
       windTestVelocity = vector(1.e-25,1.e-25,1.e-25)
    END IF
    
  END FUNCTION windTestVelocity

  

  RECURSIVE SUBROUTINE getIntersectedOctals(thisOctal,listHead,grid,nOctals,onlyChanged)  
    ! returns a linked list of all the octals which are intersected by a
    !   plane. we assume that the z-axis lies in the plane and the plane
    !   is at 45deg from the x and y axes. we only need to consider half
    !   the plane (one quadrant of the simulation space).
    !   
    ! if onlyChanged is True, we ignore cells without the Changed flag

    USE vector_mod
    USE octal_mod, only: octalListElement
    
    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT), TARGET :: thisOctal
    TYPE(octalListElement), POINTER :: listHead
    TYPE(gridType), INTENT(IN)      :: grid
    INTEGER,INTENT(INOUT)           :: nOctals   ! number of octals
    LOGICAL, INTENT(IN), OPTIONAL   :: onlyChanged
    
    LOGICAL               :: octalAdded
    TYPE(vector)     :: centrePoint
    TYPE(vector)     :: lineOrigin
    TYPE(vector)     :: starPos
    TYPE(vector)     :: nearestPlanePoint
    TYPE(octal), POINTER  :: child
    INTEGER               :: i, iSubcell, subIndex
    TYPE(vector)     :: planeVec
    TYPE(octalListElement), POINTER :: newElement
    
    octalAdded = .FALSE.
    starPos = grid%starPos1
    lineOrigin = starPos
    planeVec = vector(1.0_oc/SQRT(2.0_oc),1.0_oc/SQRT(2.0_oc),0.0_oc)

    DO iSubcell = 1, thisOctal%maxChildren, 1 
    
      centrePoint = subcellCentre(thisOctal,iSubcell) - starPos
      lineOrigin%z = centrePoint%z
      centrePoint%z = 0.0_oc
      
      IF (centrePoint%x < 0.0 .OR. centrePoint%y < 0.0) THEN
        ! we are in the wrong quadrant.
        ! this makes some assumptions about the numerical accuracy of the
        !   grid, but it should be mostly OK.
        CYCLE 
      END IF

      nearestPlanePoint = lineOrigin + (planeVec * (centrePoint .dot. planeVec))
    
      IF (inSubcell(thisOctal,iSubcell,nearestPlanePoint)) THEN
                            
        IF (thisOctal%hasChild(iSubcell)) THEN
           ! call the recursive subroutine on its child
           
          subIndex = -99
          DO i = 1, thisOctal%nChildren, 1
            IF ( thisOctal%indexChild(i) == iSubcell ) THEN
              subIndex = i  
              EXIT
            ENDIF
          ENDDO
          IF (subIndex == -99) THEN
            PRINT *, ' Panic: subindex not found'
            STOP
          ENDIF

          child => thisOctal%child(subIndex)
           
          CALL getIntersectedOctals(child,listHead,grid,nOctals) 
                            
        ELSE ! childless
          
          IF (PRESENT(onlyChanged)) THEN
            IF (onlyChanged .AND. thisOctal%changed(iSubcell)) CYCLE
          END IF
          
          IF (.NOT. octalAdded) THEN
            nOctals = nOctals + 1
            ALLOCATE(newElement)
            newElement%content => thisOctal
            newElement%next    => listHead
            newElement%inUse   = .FALSE.
            newElement%inUse(iSubcell) = .TRUE.
            listHead => newElement
            octalAdded = .TRUE.
          ELSE
            newElement%inUse(iSubcell) = .TRUE.
          END IF

        END IF
      END IF
    END DO
        
  END SUBROUTINE getIntersectedOctals

  !
  ! Based on getIntersectedOctals, but for all the octals in the first 
  ! octant (x>0, y>0 and z>0 domain).
  ! 
  ! Returns a linked list of all the octals which are in x>0, y>0 and z>0 
  ! domain. 
  RECURSIVE SUBROUTINE getOctalsInFirstOctant(thisOctal,listHead,grid,nOctals,onlyChanged)  
    !   
    ! if onlyChanged is True, we ignore cells without the Changed flag

    USE vector_mod
    USE octal_mod, only: octalListElement

    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT), TARGET :: thisOctal
    TYPE(octalListElement), POINTER :: listHead
    TYPE(gridType), INTENT(IN)      :: grid
    INTEGER,INTENT(INOUT)           :: nOctals   ! number of octals
    LOGICAL, INTENT(IN), OPTIONAL   :: onlyChanged
    
    LOGICAL               :: octalAdded
    TYPE(vector)     :: centrePoint
    TYPE(vector)     :: starPos
    TYPE(octal), POINTER  :: child
    INTEGER               :: i, iSubcell, subIndex
    TYPE(octalListElement), POINTER :: newElement
    
    octalAdded = .FALSE.
    starPos = grid%starPos1

    MAINLOOP: DO iSubcell = 1, thisOctal%maxChildren, 1 
    
      centrePoint = subcellCentre(thisOctal,iSubcell) - starPos
      
      IF (centrePoint%x < 0.0 .OR. centrePoint%y < 0.0 .OR. centrePoint%z < 0.0) THEN
         ! we are in the wrong octant.
         ! this makes some assumptions about the numerical accuracy of the
         !   grid, but it should be mostly OK.
         CYCLE  MAINLOOP

      ELSE
         IF (thisOctal%hasChild(iSubcell)) THEN
            ! call the recursive subroutine on its child           
            subIndex = -99
            DO i = 1, thisOctal%nChildren, 1
               IF ( thisOctal%indexChild(i) == iSubcell ) THEN
                  subIndex = i  
                  EXIT
               ENDIF
            ENDDO
            IF (subIndex == -99) THEN
               PRINT *, ' Panic: subindex not found'
               STOP
            ENDIF

            child => thisOctal%child(subIndex)
           
            CALL getOctalsInFirstOctant(child,listHead,grid,nOctals) 
                            
         ELSE ! childless
          
            IF (PRESENT(onlyChanged)) THEN
               IF (onlyChanged .AND. thisOctal%changed(iSubcell)) CYCLE
            END IF
            
            IF (.NOT. octalAdded) THEN
               nOctals = nOctals + 1
               ALLOCATE(newElement)
               newElement%content => thisOctal
               newElement%next    => listHead
               newElement%inUse   = .FALSE.
               newElement%inUse(iSubcell) = .TRUE.
               listHead => newElement
               octalAdded = .TRUE.
            ELSE
               newElement%inUse(iSubcell) = .TRUE.
            END IF
            
         END IF

      END IF

   END DO MAINLOOP
        
 END SUBROUTINE getOctalsInFirstOctant
  

  SUBROUTINE moveOctalListToArray(listHead,octalArray)
    ! copies all the pointers (to octals) from a linked list into an array.
    ! assumes that the SIZE of octalArray is the number of list elements.
    USE octal_mod, only: octalListElement

    IMPLICIT NONE 
    
    TYPE(octalListElement), POINTER :: listHead
    TYPE(octalWrapper),DIMENSION(:) :: octalArray

    TYPE(octalListElement), POINTER :: oldHead
    INTEGER :: element
    
    DO element = 1, SIZE(octalArray), 1

      octalArray(element)%content => listHead%content
      octalArray(element)%inUse   =  listHead%inUse
      
      oldHead => listHead
      listHead => listHead%next
      DEALLOCATE(oldHead)

    END DO
    NULLIFY(listHead)
  
  END SUBROUTINE moveOctalListToArray


  SUBROUTINE amrUpdateGrid(amrLimitScalar,amrLimitScalar2,grid)
    ! checks whether each octal has changed significantly since the last phase
    !   of the simulation, and splits/deletes it as necessary.

    IMPLICIT NONE

    real(double), INTENT(IN) :: amrLimitScalar, amrLimitScalar2 
      ! 'limitScalar' is the value the decideSplit function uses to
      !   decide whether or not to split cell.
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid through to the 
                                          !   routines that this subroutine calls
    
    TYPE(OCTAL), POINTER :: thisOctal

    IF ((.NOT. grid%adaptive) .OR. (.NOT.grid%geometry(1:6)=="ttauri")) THEN
      print *, 'amrUpdateGrid probably doesn''t work with this setup'
      stop
    END IF
  
    ! first we check whether we are going to delete any child octals
   
    !print *, 'Deleting cells no longer needed...'
    !thisOctal => grid%octreeRoot
    !CALL amrUpdateGridDelete(thisOctal)
   
    ! then we check if we are adding any new octals
    
    print *, 'Adding new cells to grid...'
    thisOctal => grid%octreeRoot
    CALL amrUpdateGridAdd(thisOctal)

    ! finally, check whether any subcells have changed significantly
    print *, 'Flagging updated cells...'
    thisOctal => grid%octreeRoot
    CALL amrUpdateGridChanged(thisOctal)
    
  CONTAINS

    RECURSIVE SUBROUTINE amrUpdateGridDelete(thisOctal)
      ! if a child should no longer exist, we delete it.

      TYPE(octal), TARGET, INTENT(INOUT) :: thisOctal
      TYPE(octal), POINTER :: thisChild
      logical :: splitInAzimuth
      INTEGER :: iChild, iSubcell

      DO iChild = 1, thisOctal%nChildren, 1
        
        thisChild => thisOctal%child(iChild)
        
        DO iSubcell = 1, thisChild%maxChildren, 1
          IF (thisOctal%indexChild(iSubcell) == iChild) EXIT
        END DO 
        
        IF (.NOT. decideSplit(thisOctal,iSubcell,amrLimitScalar,amrLimitScalar2,grid,.false.,splitInAzimuth)) THEN
          PRINT *, 'Deleting unneeded children'
          CALL deleteChild(parent=thisOctal, childToDelete=iSubcell, &
                           adjustParent=.FALSE., grid=grid,          &
                           adjustGridInfo=.TRUE.)
          ! setting adjustParent and adjustGridInfo to be true might
          !   not be sensible here?

          thisOctal%changed(iSubcell) = .TRUE.
        ELSE
          CALL amrUpdateGridDelete(thisChild)
        END IF

      END DO  
    
    END SUBROUTINE amrUpdateGridDelete
    
    RECURSIVE SUBROUTINE amrUpdateGridAdd(thisOctal)
      ! subdivide any octals that now exceed the threshold

      TYPE(octal), TARGET, INTENT(INOUT) :: thisOctal
      TYPE(octal), POINTER :: thisChild
      INTEGER :: iSubcell, j
      logical :: splitInAzimuth

      DO iSubcell = 1, thisOctal%maxChildren, 1
        IF (thisOctal%hasChild(iSubcell)) THEN 
          ! find the child
          DO j = 1, thisOctal%nChildren, 1
            IF (thisOctal%indexChild(j) == iSubcell) THEN
              thisChild => thisOctal%child(j)
              CALL amrUpdateGridAdd(thisChild)
              EXIT
            END IF
          END DO
        ELSE 
          IF (decideSplit(thisOctal,iSubcell,amrLimitScalar,amrLimitScalar2,grid,.false.,splitinAzimuth)) THEN
            PRINT *, 'Adding new child'
            CALL addNewChild(thisOctal,iSubcell,grid,adjustGridInfo=.TRUE.,splitAzimuthally=splitinAzimuth)
            thisOctal%child(iSubcell)%changed = .TRUE.
            grid%nOctals = grid%nOctals + 1
          END IF
        END IF 
          
      END DO 
        
    END SUBROUTINE amrUpdateGridAdd

    RECURSIVE SUBROUTINE amrUpdateGridChanged(thisOctal)
      use density_mod, only: density
      ! update the octals in the grid
      ! flag any octals that have changed "significantly" since the last phase

      ! NB This is only for T Tauri models. It should be made more elegant and general!
      
      TYPE(octal), TARGET, INTENT(INOUT) :: thisOctal
      TYPE(octal), POINTER :: thisChild
      INTEGER :: iSubcell
      REAL :: newDensity
      INTEGER :: iChild

      DO iSubcell = 1, thisOctal%maxChildren, 1

!        newDensity = TTauriDensity(subcellCentre(thisOctal,iSubcell),grid)
         newDensity = real(Density(subcellCentre(thisOctal,iSubcell),grid))
          IF ( ABS((newDensity/(MAX(thisOctal%rho(iSubcell),1.d-25))-1.0)) > 0.1 ) &
            thisOctal%changed(iSubcell) = .TRUE.
         
        CALL calcTTauriMassVelocity(thisOctal,iSubcell,grid)

      END DO

      DO iChild = 1, thisOctal%nChildren, 1
        thisChild => thisOctal%child(iChild)
        CALL amrUpdateGridChanged(thisChild)
      END DO 
      
    END SUBROUTINE amrUpdateGridChanged

  END SUBROUTINE amrUpdateGrid

  subroutine returnKappa(grid, thisOctal, subcell, ilambda, lambda, kappaSca, kappaAbs, kappaAbsArray, kappaScaArray, &
       rosselandKappa, kappap, atthistemperature, kappaAbsDust, kappaAbsGas, kappaScaDust, kappaScaGas, debug, reset_kappa)
    use inputs_mod, only: nDustType, mie, includeGasOpacity, lineEmission, dustPhysics
    use atom_mod, only: bnu
    use gas_opacity_mod, only: returnGasKappaValue
#ifdef PHOTOION
    use inputs_mod, only: photoionization, hOnly
    use phfit_mod, only : phfit2
#endif
    implicit none
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer, optional :: ilambda
    real, optional :: lambda
    real(double), intent(out), optional :: kappaSca, kappaAbs
    real(double), optional, intent(out) :: kappaAbsArray(:), kappaScaArray(:)
    real(double), optional, intent(out) :: rosselandKappa
    real(double), optional, intent(out) :: kappaAbsDust, kappaScaDust, kappaAbsGas, kappaScaGas
    logical, optional :: debug
    real, optional, intent(out) :: kappap
    real, optional :: atthistemperature
    logical, optional, intent(in) :: reset_kappa
    real :: temperature
    real :: frac
    real :: tlambda
!    real, parameter :: sublimationTemp = 1500., subRange = 100.
    real(double) :: tArray(1000)
    real(double) :: freq, dfreq, norm !,  bnutot
    integer :: i,j,m,itemp
    real :: fac
    real(double), allocatable, save :: tgasArray(:)
    real(double), allocatable, save :: oneKappaAbsT(:,:)
    real(double), allocatable, save :: oneKappaScaT(:,:)

    logical,save :: firsttime = .true. 
    integer(double),save :: nlambda
    real(double) :: tgas
    
#ifdef PHOTOION
    real(double) :: kappaH, kappaHe
    real :: e, h0, he0
#endif
#ifdef _OPENMP
!$OMP THREADPRIVATE (firstTime, nlambda, tgasArray, oneKappaAbsT, oneKAppaScaT)
#endif


    if ( present(reset_kappa) ) then 
       if ( reset_kappa ) then 
          if ( allocated(tgasArray)    ) deallocate ( tgasArray    )
          if ( allocated(oneKappaAbsT) ) deallocate ( oneKappaAbsT )
          if ( allocated(oneKappaScaT) ) deallocate ( oneKappaScaT )
          firsttime = .true.
          return
       end if
    end if

    !! Commented out by DAR 14/10/08 as I don't think this is doing anything here. Moved to
    !! a more relevant place - rosselandkappa and others
!    temperature = thisOctal%temperature(subcell)


!    if (temperature < sublimationTemp) frac = 1.
!    if (temperature > (sublimationTemp+subRange)) frac = 0.
!    
!    if ((temperature > sublimationTemp).and.(temperature < (sublimationTemp+subRange))) then
!       frac = 1.-(temperature-sublimationTemp)/subRange
!    endif

!    frac = max(1.e-20,frac)

    frac = 1.

    if (allocated(oneKappaAbsT)) then
       if  ((SIZE(oneKappaAbsT,1) /= SIZE(grid%oneKappaAbs,2)).or. &
            (SIZE(oneKappaAbsT,2) /= SIZE(grid%oneKappaAbs,1))) then
          deallocate(oneKappaAbsT)
          deallocate(oneKappaScaT)
          if (allocated(tGasArray)) deallocate(tGasArray)
          firstTime = .true.
       endif
    endif

    if(firsttime.and.mie) then

       allocate(oneKappaAbsT(grid%nlambda, ndusttype))
       allocate(oneKappaScaT(grid%nlambda, ndusttype))
       if (allocated(tGasArray)) then
          if (size(tGasArray) /= grid%nLambda) then
             deallocate(tGasArray)
             allocate(tGasArray(grid%nLambda))
          else
             allocate(tGasArray(grid%nLambda))
          endif
       endif
       do i = 1, ndusttype
          do m = 1, grid%nLambda
             oneKappaAbsT(m, i) = grid%oneKappaAbs(i,m)
             oneKappaScaT(m, i) = grid%oneKappaSca(i,m)
          enddo
       enddo

       nlambda = grid%nlambda
       firsttime = .false.
    endif

    if (PRESENT(kappaAbsArray)) then

      if(grid%onekappa) then

         if (nDustType .eq. 1) then
            kappaAbsArray(1:nLambda) = thisOctal%rho(subcell) * oneKappaAbsT(1:nlambda,1) &
                 * thisOctal%dustTypeFraction(subcell,1)
         else
            kappaAbsArray(1:grid%nLambda) = 0.            
            do i = 1, nDustType
               kappaAbsArray(1:nLambda) = kappaAbsArray(1:nLambda) + thisOctal%dustTypeFraction(subcell, i) * & 
                    oneKappaAbsT(1:nLambda,i)*thisOctal%rho(subcell) * frac
            enddo
         endif
      else
         kappaAbsArray(1:nlambda) = thisOctal%kappaAbs(subcell,:)
      endif
       if (includeGasOpacity) then
          call returnGasKappaValue(grid,thisOctal%temperature(subcell), thisOctal%rho(subcell),  kappaAbsArray=tarray)
          kappaAbsArray(1:grid%nLambda) = kappaAbsArray(1:grid%nLambda) + tarray(1:grid%nLambda)*thisOctal%rho(subcell)
       endif
!       write(*,*) nDustType,thisOctal%dusttypeFraction(subcell,1), grid%oneKappaAbs(1,1:grid%nLambda)

    endif
       
    if (PRESENT(kappaScaArray)) then
       if(grid%onekappa) then
 
          if (nDustType .eq. 1) then
             kappaScaArray(1:nLambda) = thisOctal%rho(subcell) * oneKappaScaT(1:nlambda,1) * &
                  thisOctal%dustTypeFraction(subcell,1)
          else
             kappaScaArray(1:grid%nLambda) = 0.   
             do i = 1, nDustType
                kappaScaArray(1:nLambda) = kappaScaArray(1:nLambda) + thisOctal%dustTypeFraction(subcell, i) * & 
                     oneKappaScaT(1:nLambda,i)*thisOctal%rho(subcell) * frac
             enddo
          endif
       else
          kappaScaArray(1:nlambda) = thisOctal%kappaSca(subcell,:)
       endif
       if (includeGasOpacity) then
          call returnGasKappaValue(grid, temperature, thisOctal%rho(subcell),  kappaScaArray=tarray)
          kappaScaArray(1:grid%nLambda) = kappaScaArray(1:grid%nLambda) + tarray(1:grid%nLambda) !*thisOctal%rho(subcell)
       endif
    endif


    if (PRESENT(kappaSca)) then
       kappaSca = 0.
       if (dustPhysics) then
          IF (.NOT.PRESENT(lambda)) THEN
             if (grid%oneKappa) then
                kappaSca = 0.
                
                if(ndusttype .eq. 1) then
                   kappaSca = OneKappaScaT(iLambda,1) * thisOctal%rho(subcell)  * thisOctal%dustTypeFraction(subcell, 1)
                   if (present(debug)) write(*,*) "kappasca1 ",kappaSca, oneKappaScaT(ilambda,1), & 
                        thisOctal%rho(subcell), thisOctal%dusttypeFraction(subcell,1)
                else
                   do i = 1, nDustType
                      kappaSca = kappaSca + thisOctal%dustTypeFraction(subcell, i) * grid%oneKappaSca(i,iLambda) * &
                           thisOctal%rho(subcell)
                      if (present(debug)) write(*,*) "kappasca2 ",kappaSca, grid%oneKappaSca(i,ilambda), &
                           thisOctal%rho(subcell), thisOctal%dusttypeFraction(subcell,i)
                      
                   enddo
                endif
             else 
                ! For line computation (without dust).
                ! Needs modification for a model which include gas and dust here.
                ! This is a temporarily solution
                if (lineEmission) then
                   kappaSca = thisOctal%kappaSca(subcell,1)
                else
                   kappaSca = thisOctal%kappaSca(subcell,iLambda)
                endif
             end if
          else
             kappaSca = 0.
             itemp = ilambda
             if (ilambda == grid%nLambda) itemp = itemp - 1
             do i = 1, nDustType
                if (grid%nLambda == 1) then
                   kappaSca = kappaSca +  thisOctal%dustTypeFraction(subcell, i) * &
                        grid%oneKappaSca(i,1) * thisOctal%rho(subcell)
                else
                   kappaSca = kappaSca + thisOctal%dustTypeFraction(subcell, i) * &
                        logint(dble(lambda), dble(grid%lamArray(itemp)), dble(grid%lamArray(itemp+1)), &
                        grid%oneKappaSca(i,itemp)*thisOctal%rho(subcell), &
                        grid%oneKappaSca(i,itemp+1)*thisOctal%rho(subcell))
                endif
             enddo
          endif
          kappaSca = kappaSca * frac
          if (present(debug)) write(*,*) "kappasca xxx ", kappasca
       endif
    endif
    if (PRESENT(kappaScaDust)) kappaScaDust = kappaSca

    if (PRESENT(kappaAbs)) then
       kappaAbs = 0.
       
       IF (.NOT.PRESENT(lambda)) THEN
          if (grid%oneKappa) then
             kappaAbs = 0.
             if(ndustType .eq. 1) then
                kappaAbs = oneKappaAbsT(iLambda,1)*thisOctal%rho(subcell) * thisOctal%dustTypeFraction(subcell, 1)
             else
                do i = 1, nDustType
                   kappaAbs = kappaAbs + thisOctal%dustTypeFraction(subcell, i) * grid%oneKappaAbs(i,iLambda)*thisOctal%rho(subcell)
                enddo
             endif
          else
             ! For line computation (without dust).
             ! Needs modification for a model which include gas and dust here.
             ! This is a temporarily solution
             if (lineEmission) then
                kappaAbs = thisOctal%kappaAbs(subcell,1)
             else
                kappaAbs = thisOctal%kappaAbs(subcell,iLambda)
             endif
          end if
       else
          kappaAbs = 0.
          itemp = ilambda
          if (ilambda == grid%nLambda) itemp = itemp - 1
          if (dustPhysics) then
             if (ndusttype .eq. 1) then
                if (grid%nLambda == 1) then
                   kappaAbs = thisOctal%dustTypeFraction(subcell, 1) *  &
                        grid%oneKappaAbs(1,1) * thisOctal%rho(subcell)
                else
                   kappaAbs = thisOctal%dustTypeFraction(subcell, 1) *  &
                        logint(dble(lambda), dble(grid%lamArray(itemp)), dble(grid%lamArray(itemp+1)), &
                        oneKappaAbsT(itemp,1)*thisoctal%rho(subcell), &
                        oneKappaAbsT(itemp+1,1)*thisoctal%rho(subcell))
                endif
                
             else
                itemp = ilambda
                if (ilambda == grid%nLambda) itemp = itemp - 1          
                do i = 1, nDustType
                   !write(*,*) grid%oneKappaAbs(i,iLambda),grid%oneKappaAbs(i,iLambda+1),thisOctal%rho(subcell), &
                   !     dble(grid%lamArray(ilambda)), dble(grid%lamArray(ilambda+1)), ilambda, dble(lambda)
                   !write(*,*) logint(dble(lambda), dble(grid%lamArray(ilambda)), dble(grid%lamArray(ilambda+1)), &
                   !     grid%oneKappaAbs(i,iLambda)*thisOctal%rho(subcell), &
                   !     grid%oneKappaAbs(i,iLambda+1)*thisOctal%rho(subcell))
                   !write(*,*) i, subcell
                   !write(*,*) nDusttype,kappaAbs 
                   !write(*,*) thisOctal%dustTypeFraction(subcell, i)
                   if (grid%nLambda == 1) then
                      kappaAbs = kappaAbs + thisOctal%dustTypeFraction(subcell, i) * &
                           grid%oneKappaAbs(i,1)*thisOctal%rho(subcell)
                      
                   else
                      kappaAbs = kappaAbs + thisOctal%dustTypeFraction(subcell, i) * &
                           logint(dble(lambda), dble(grid%lamArray(itemp)), dble(grid%lamArray(itemp+1)), &
                           grid%oneKappaAbs(i,itemp)*thisOctal%rho(subcell), &
                           grid%oneKappaAbs(i,itemp)*thisOctal%rho(subcell))
                   endif
                enddo
             endif
          endif
       endif
!       write(*,*) nDustType,thisOctal%dusttypeFraction(subcell,1), grid%oneKappaAbs(1,1:grid%nLambda)
       
       kappaAbs = kappaAbs * frac

    endif
   if (PRESENT(kappaAbsDust)) kappaAbsDust = kappaAbs


   if (includeGasOpacity.and.(.not.PRESENT(rosselandkappa))) then

      if (PRESENT(kappaAbs)) then
         temperature = thisOctal%temperature(subcell)
         
         if (PRESENT(atthistemperature)) then
            temperature = atthistemperature
         endif
         if (.not.PRESENT(lambda)) then
            tlambda = grid%lamArray(iLambda)
         else
            tlambda = lambda
         endif
         call returnGasKappaValue(grid, temperature, thisOctal%rho(subcell), tlambda, ilambda = ilambda, kappaAbs=tGas)
         kappaAbs = kappaAbs + tGas*thisOctal%rho(subcell)
      endif
      
      if (PRESENT(kappaSca)) then
         if (.not.PRESENT(lambda)) then
            tlambda = grid%lamArray(iLambda)
         else
            tlambda = lambda
         endif
         call returnGasKappaValue(grid, temperature, thisOctal%rho(subcell), tlambda, kappaSca=tGas)
         kappaSca = kappaSca + tGas *thisOctal%rho(subcell)
      endif

      if (PRESENT(kappaScaArray)) then
         call returnGasKappaValue(grid, temperature, thisOctal%rho(subcell), tlambda, kappaScaArray=tGasArray)
         kappaScaArray = kappaScaArray + tGasArray *thisOctal%rho(subcell)
      endif
   endif

   if (PRESENT(rosselandKappa)) then
      temperature = thisOctal%temperature(subcell)

      if (PRESENT(atthistemperature)) then
         temperature = atthistemperature
      endif
      if (temperature < grid%tempRossArray(grid%nTempRossArray)) then
         call locate(grid%tempRossArray, grid%nTempRossArray, temperature, m)
         fac = (temperature - grid%tempRossArray(m))/(grid%tempRossArray(m+1)-grid%tempRossArray(m))
      else
         m = grid%nTempRossArray-1
         fac = 1.
      endif
      rosselandKappa = 0.
      
      do i = 1, nDustType
         rosselandKappa = rosselandKappa + (grid%kappaRossArray(i, m) + &
              fac*(grid%kappaRossArray(i,m+1)-grid%kappaRossArray(i,m))) * &
              max(1.d-30,thisOctal%dustTypeFraction(subcell, i))
      enddo

   endif
   

   if (PRESENT(kappap)) then
      temperature = thisOctal%temperature(subcell)

      if (PRESENT(atthistemperature)) then
         temperature = atthistemperature
      endif
      kappaP = 0.d0
      norm = 0.d0
      if (includeGasOpacity) then
         call returnGasKappaValue(grid,real(temperature), thisOctal%rho(subcell),  kappaAbsArray=tarray)
      endif
      do i = 2, grid%nLambda
         freq = cSpeed / (grid%lamArray(i)*1.e-8)
         dfreq = cSpeed / (grid%lamArray(i-1)*1.e-8) - cSpeed / (grid%lamArray(i)*1.e-8)
         do j = 1, nDustType
            kappaP = kappaP + real(thisOctal%dustTypeFraction(subcell, j) * dble(grid%oneKappaAbs(j,i)) * &
                 thisOctal%rho(subcell) *&
                 dble(bnu(dble(freq),dble(temperature)))  * dfreq)

            if (includeGasOpacity) then
               kappaP = kappaP + real(tarray(i)*thisOctal%rho(subcell) * dble(bnu(dble(freq),dble(temperature)))  * dfreq)
            endif

         enddo
         norm = norm + dble(bnu(dble(freq),dble(temperature)))  * dfreq
      enddo
      if (norm /= 0.d0) then
         kappaP = real((kappaP / norm) /1.d10)
      else
         kappaP = tiny(kappap)
      endif
   endif
   
#ifdef PHOTOION
   if (photoionization) then

      if (PRESENT(kappaAbs)) then
         if (present(lambda)) then
            e = real((hCgs * (cSpeed / (lambda * 1.e-8))) * ergtoev)
         else
            e = real((hCgs * (cSpeed / (grid%lamArray(iLambda) * 1.e-8))) * ergtoev)
            !            write(*,*) "! using rough grid"
         endif
         call phfit2(1, 1, 1 , e , h0)
         call phfit2(2, 2, 1 , e , he0)
         kappaH =  thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,1) * h0
         if (.not. hOnly) then
            kappaHe = thisOctal%nh(subcell)*grid%ion(3)%abundance*thisOctal%ionFrac(subcell,3) * he0
         else
            kappaHe = 0.d0
         end if
         kappaAbs = kappaAbs + (kappaH + kappaHe)
      endif
      if (PRESENT(kappaAbsGas)) kappaAbsGas = (kappaH + kappaHe)
      if (PRESENT(kappaSca)) then
         kappaSca = kappaSca + thisOctal%ne(subcell) * sigmaE * 1.e10
         if (present(debug)) write(*,*) "kappasca3 ",kappasca, thisOctal%ne(subcell) * sigmaE * 1.e10, thisOctal%ne(subcell),&
              thisOctal%rho(subcell)/mHydrogen,thisOctal%nh(subcell)
      endif
      if (PRESENT(kappaScaGas)) kappaScaGas = thisOctal%ne(subcell) * sigmaE * 1.e10 
   endif
#else
   if (PRESENT(kappaAbsGas)) kappaAbsGas = 0.0
   if (PRESENT(kappaScaGas)) kappaScaGas = 0.0
#endif

  end subroutine returnKappa
   
  subroutine interpFromParent(centre, cellSize, grid, temperature, density, dusttypeFraction, thisEtaLine)
    type(GRIDTYPE) :: grid
    real(double) :: cellSize
    real :: temperature
    real(double) :: density, thisEtaLine
    real(double) :: dusttypeFraction(:)
    real(double), allocatable :: tdusttype(:,:)
    real(double), allocatable :: rho(:), etaline(:)
    real, allocatable :: temp(:)
    type(VECTOR) :: centre, octVec
    real(double) :: r
    integer :: j, i

    if (grid%octreeRoot%oneD) then
       write(*,*) "interp from parent not implemented for one-d case"
       stop
    endif

    r = cellSize/2.d0 + 0.1d0 * grid%halfSmallestSubcell
    if (grid%octreeRoot%twod) then
       j = 4
       allocate(tDustType(1:j,1:SIZE(dusttypeFraction)))
       allocate(rho(1:j), temp(1:j), etaline(1:j))

       octVec = centre + r * VECTOR(1.d0, 0.d0, 0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(1), rho=rho(1), dusttypeFraction=tdusttype(1,:), &
            etaline=etaLine(1))
       
       octVec = centre + r * VECTOR(-1.d0,0.d0,0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(2), rho=rho(2), dusttypeFraction=tdusttype(2,:), &
            etaline=etaline(2))
       
       octVec = centre + r * VECTOR(0.d0,0.d0,-1.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(3), rho=rho(3), dusttypeFraction=tdusttype(3,:), &
            etaline=etaLine(3))
       
       octVec = centre + r * VECTOR(0.d0,0.d0,+1.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(4), rho=rho(4), dusttypeFraction=tdusttype(4,:), &
            etaline=etaline(4))
    else
       j = 6
       allocate(tDustType(1:j,1:SIZE(dusttypeFraction)))
       allocate(rho(1:j), temp(1:j),etaline(1:j))
       octVec = centre + r * VECTOR(1.d0, 0.d0, 0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(1), rho=rho(1), dusttypeFraction=tdusttype(1,:), &
            etaline=etaline(1))
       
       octVec = centre + r * VECTOR(-1.d0,0.d0,0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(2), rho=rho(2), dusttypeFraction=tdusttype(2,:), &
            etaline=etaline(2))
       
       octVec = centre + r * VECTOR(0.d0,0.d0,-1.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(3), rho=rho(3), dusttypeFraction=tdusttype(3,:), &
            etaline=etaline(3))
       
       octVec = centre + r * VECTOR(0.d0,0.d0,+1.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(4), rho=rho(4), dusttypeFraction=tdusttype(4,:), &
            etaline=etaline(4))

       octVec = centre + r * VECTOR(0.d0,+1.d0,0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(5), rho=rho(5), dusttypeFraction=tdusttype(5,:), &
            etaline=etaline(5))

       octVec = centre + r * VECTOR(0.d0,-1.d0,0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(6), rho=rho(6), dusttypeFraction=tdusttype(6,:), &
            etaline=etaline(6))

    endif
  

    density = sum(rho) / dble(j)

    thisetaline = sum(etaline) / dble(j)

    temperature = sum(temp) / real(j)

    do i = 1, SIZE(dustTypefraction)
       dusttypeFraction(i) = SUM(tdusttype(1:j,i)) / dble(j)
    enddo
  end subroutine interpFromParent

  recursive subroutine myTauSmooth(thisOctal, grid, ilambda, converged, inheritProps, interpProps, photosphereSplit)
    use inputs_mod, only : tauSmoothMin, tauSmoothMax, erOuter, router, maxDepthAmr, rinner, rGapInner
    use inputs_mod, only : maxMemoryAvailable
    use memory_mod, only : humanreadablememory, globalMemoryFootprint
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal
    logical, optional :: inheritProps, interpProps, photoSphereSplit
    integer :: subcell, i, ilambda
    logical :: converged
    real(double) :: kabs, ksca, r, fac,  kabsDust, kScaDust
    type(VECTOR) :: dirVec(6), centre, octVec, aHat, rVec
    real :: thisTau, neighbourTau
    integer :: neighbourSubcell, j, nDir
    logical :: split, outofMemory
    logical, save :: firsttime = .true., firstTimeMem = .true.
    character(len=80) :: message

!$OMP THREADPRIVATE (firstTime, firstTimeMem)

    kabs = 0.d0; ksca = 0.d0

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call myTauSmooth(child, grid, ilambda, converged, inheritProps, interpProps, photosphereSplit)
                exit
             end if
          end do
       else

          split = .true.

          call returnKappa(grid, thisOctal, subcell, ilambda=ilambda,&
               kappaSca=ksca, kappaAbs=kabs, kappaScaDust=kscaDust, kappaAbsDust=kabsDust)

          thisTau  = real(thisOctal%subcellSize * (kscaDust + kabsDust))

          r = thisOctal%subcellSize/2.d0 + grid%halfSmallestSubcell * 0.1d0
          centre = subcellCentre(thisOctal, subcell)
          if (thisOctal%cylindrical) then 
             aHat = (0.1*grid%halfSmallestsubcell)*randomUnitVector()
             centre = centre + aHat
          endif
          if (thisOctal%threed) then
             if (.not.thisOctal%cylindrical) then
                nDir = 6
                dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
                dirVec(2) = VECTOR( 0.d0,+1.d0,  0.d0)
                dirVec(3) = VECTOR(+1.d0, 0.d0,  0.d0)
                dirVec(4) = VECTOR(-1.d0, 0.d0,  0.d0)
                dirVec(5) = VECTOR( 0.d0,-1.d0,  0.d0)
                dirVec(6) = VECTOR( 0.d0, 0.d0, -1.d0)
             else
                nDir = 4
                aHat = centre
                call normalize(aHat)
                dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
                dirVec(2) = aHat
                dirVec(3) = (-1.d0)*aHat
                dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
             endif
          else if (thisOctal%twod) then
             nDir = 4
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)  
             dirVec(2) = VECTOR(-1.d0,0.d0, 0.d0)
             dirVec(3) = VECTOR( 0.d0, 0.d0,  1.d0)
             dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
          else
             nDir = 2
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(-1.d0, 0.d0, 0.d0)
          endif

          do j = 1, nDir
             octVec = centre + r * dirvec(j)
             if (inOctal(grid%octreeRoot, octVec)) then
                neighbourOctal => thisOctal
                neighbourSubcell = subcell
                call findSubcellLocal(octVec, neighbourOctal, neighbourSubcell)
                call returnKappa(grid, neighbourOctal, neighboursubcell, ilambda=ilambda,&
                     kappaSca=ksca, kappaAbs=kabs, kappaScaDust=kscaDust, kappaAbsDust=kabsDust)

                neighbourTau = real(neighbourOctal%subcellSize * (kscaDust + kabsDust))

                if ((grid%geometry.eq."whitney").and.&
                     (modulus(subcellCentre(thisOctal,subcell)) > 0.9*erouter/1.e10)) split = .false.

                rVec = subcellCentre(thisOctal,subcell)

                if ((grid%geometry.eq."warpeddisc").and.&
                     (sqrt(rVec%x**2+rVec%y**2)) > 0.9*router) split = .false.


                if ((grid%geometry.eq."shakara").and.&
                     ((sqrt(rVec%x**2+rVec%y**2)+thisOctal%subcellsize/2.d0) < rinner)) split = .false.

                if ((grid%geometry.eq."shakara").and.&
                     (sqrt(rVec%x**2+rVec%y**2) > 0.9*rGapInner)) split = .false.

                if ((grid%geometry.eq."shakara").and.&
                     (sqrt(rVec%x**2+rVec%y**2) > 0.9*rOuter)) split = .false.



                if (thisOctal%nDepth == maxDepthamr) then
                   split = .false.
                   if (firstTime) then
                      write(message,'(a,i3)') "AMR cell depth capped at: ",maxDepthamr
                      call writeWarning(message)
                      firstTime = .false.
                   endif
                endif

                outofmemory = .false.
                if (globalMemoryFootprint > maxMemoryAvailable) then
                   split = .false.
                   outofmemory = .true.
                   if (firstTimeMem) then
                      write(message,'(a)') "Maxmimum memory exceeded for grid :"//humanReadableMemory(globalMemoryFootprint)
                      call writeWarning(message)
                      firstTimeMem = .false.
                   endif
                endif


                if ((min(thisTau, neighbourTau) < tauSmoothMin).and.(max(thisTau, neighbourTau) > tauSmoothMax).and.split) then
                   if (thisTau > neighbourTau) then
                      call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                           inherit=inheritProps, interp=interpProps)
                      converged = .false.
                      return
                   endif
                endif

!                call returnKappa(grid, thisOctal, subcell, ilambda, rosselandKappa = kAbs)
!                tauross = thisOctal%subcellSize*kAbs*thisOctal%rho(subcell)*1.d10
!                if ((tauRoss > 500.d0).and.split) then
!!                   if (myrankGlobal==1) write(*,*) "Splitting with tauRoss ",tauross
!                   call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
!                        inherit=inheritProps, interp=interpProps)
!                   converged = .false.
!                   return
!                endif


                fac = abs(thisOctal%temperature(subcell)-neighbourOctal%temperature(neighbourSubcell)) / &
                     thisOctal%temperature(subcell)

                if (split.and.(thisTau > neighbourTau).and.(fac > 0.5d0).and. &
                     (thisOctal%temperature(subcell)>100.d0).and.(thisTau > 1.d-4)) then
                   call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                        inherit=inheritProps, interp=interpProps)
                   converged = .false.
                   return
                endif

                if (PRESENT(photosphereSplit)) then
                   if (photosphereSplit.and.(.not.outofmemory)) then
                      if ((thisOctal%etaLine(subcell) /= 0.d0).and.(neighbourOctal%etaLine(neighbourSubcell)/=0.d0)) then
                         if ((j==3).or.(j==4)) then
                            fac = abs(neighbourOctal%etaLine(neighbourSubcell)-thisOctal%etaLine(subcell))
                            if (split.and.(fac > 0.2d0).and.(thisOctal%etaLine(subcell) > 0.1d0).and. & 
                                 (thisOctal%etaLine(subcell) < 1.d0).and. &
                                 (thisOctal%etaLine(subcell)>neighbourOctal%etaLine(neighbourSubcell))) then
!                               if (myrankGlobal == 1) write(*,*) &
!                               " tau split ", fac, " eta ",thisOctal%etaline(subcell), "depth ",thisOctal%ndepth
                               call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                                    inherit=inheritProps, interp=interpProps)
                               converged = .false.
                               return
                            endif
                         endif
                      endif
                   endif
                endif
             endif
          enddo
       endif
    end do

  end subroutine myTauSmooth

  subroutine myScaleSmooth(factor, grid, converged, &
       inheritProps, interpProps)
    use memory_mod, only : globalMemoryFootprint
    use inputs_mod, only : maxMemoryAvailable
    type(gridtype) :: grid
    real :: factor
    integer :: nTagged
    logical, optional :: inheritProps, interpProps
    !
    logical :: converged

    if ( factor <= 0 ) then 
       call writewarning ("Splitting factor in myScaleSmooth is <= 0")
       return
    end if
    call zeroChiLineLocal(grid%octreeRoot)
    nTagged = 0
    call tagScaleSmooth(nTagged, factor, grid%octreeRoot, grid,  converged, &
         inheritProps, interpProps)
    call splitTagged(grid%octreeRoot, grid, inheritProps, interpProps)
    if (globalMemoryFootprint > maxMemoryAvailable) converged = .true.

  end subroutine myScaleSmooth

  recursive subroutine zeroChiLineLocal(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call zeroChiLineLocal(child)
                exit
             end if
          end do
       else

          call allocateAttribute(thisOctal%chiLine,thisOctal%maxChildren)
          thisOctal%chiLine = 0.d0

       endif
    enddo
  end subroutine zeroChiLineLocal

  recursive subroutine zeroDensity(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call zeroDensity(child)
                exit
             end if
          end do
       else

          thisOctal%rho(subcell) = 1.e-25
          thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)

       endif
    enddo
  end subroutine zeroDensity

  recursive subroutine zeroEtaCont(thisOctal)

    implicit none
    type(octal), pointer              :: thisOctal
    integer :: subcell
    type(octal), pointer  :: child 

    if (thisOctal%nChildren > 0) then
       ! call this subroutine recursively on each of its children
       do subcell = 1, thisOctal%nChildren, 1 
          child => thisOctal%child(subcell)
          call zeroEtaCont(child)
       end do 
    end if
    thisOctal%biasCont3D = 1.d0
    thisOctal%etaCont = 0.d0
  end subroutine zeroEtaCont

  recursive subroutine assignDensitiesBlandfordPayne(grid,thisOctal)
    use inputs_mod, only :  vturb, dw_temperature
    use magnetic_mod, only : inflowBlandfordPayne, velocityBlandfordPayne, rhoBlandfordPayne
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(VECTOR) :: cellCentre
    integer :: subcell, i
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call assignDensitiesBlandfordPayne(grid, child)
                exit
             end if
          end do
       else
          thisOctal%inflow(subcell) = .false.
!          cellCentre = subcellCentre(thisOctal, subcell)
!          call subcellCorners(thisOctal, subcell, corner)
!          if (inflowBlandfordPayne(corner(1:thisoctal%maxChildren))&
!               .and.inflowBlandfordPayne(cellCentre)) then

             thisOctal%inflow(subcell) = .true.
             thisOctal%rho(subcell) = rhoBlandfordPayne(cellCentre)
             thisOctal%velocity(subcell) = velocityBlandfordPayne(cellCentre)
             thisOctal%temperature(subcell) = real(DW_temperature)
             thisOctal%fixedTemperature(subcell) = .true.

!             IF ((thisoctal%threed).and.(subcell == 8)) &
!                  CALL fillVelocityCorners(thisOctal,velocityBlandfordPayne)
          
!             IF ((thisoctal%twod).and.(subcell == 4)) &
!                  CALL fillVelocityCorners(thisOctal,VelocityBlandfordPayne)

          endif

          if (associated(thisOctal%microturb)) thisOctal%microturb(subcell) = vturb



    enddo
  end subroutine assignDensitiesBlandfordPayne

  recursive subroutine assignDensitiesAlphaDisc(grid, thisOctal)
    use magnetic_mod, only : rhoAlphaDisc, velocityAlphaDisc
    use inputs_mod, only : rSublimation, alphaDiscTemp, grainFrac, rinner, router
    type(GRIDTYPE) :: grid
    real(double) :: thisRho, r, fac
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(VECTOR) :: cellCentre
    integer :: subcell, i
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call assignDensitiesAlphaDisc(grid, child)
                exit
             end if
          end do
       else
          cellCentre = subcellCentre(thisOctal, subcell)
          thisRho  = rhoAlphaDisc(grid, cellCentre)
          r = sqrt(cellCentre%x**2 + cellCentre%y**2)
          thisOctal%dustTypeFraction(subcell,:) = 0.d0
          thisOctal%fixedTemperature(subcell) = .true.
          if ((r > rSublimation).and.(thisRho >= thisOctal%rho(subcell))) then
             thisOctal%dustTypeFraction(subcell,1) = grainFrac(1)
             thisOctal%fixedTemperature(subcell) = .false.
             thisOctal%inflow(subcell) = .true.
             thisOctal%temperature(subcell) = 10.
             if (r < rSublimation*1.05d0) then
                fac = ((rSublimation*1.05d0 - r)/(0.005d0*rSublimation))
                thisOctal%dustTypeFraction(subcell,1) = grainFrac(1)*exp(-fac)
             endif
          endif
          if ((thisRho > thisOctal%rho(subcell)).and.(r < rSublimation)) then
             thisOctal%temperature(subcell) = alphaDiscTemp
          endif
          if ( (r > Rinner).and.(r < rOuter).and.(thisRho > thisOctal%rho(subcell))) then
             thisOctal%velocity(subcell) = velocityAlphaDisc(cellcentre)
             CALL fillVelocityCorners(thisOctal,velocityAlphaDisc)
             thisOctal%iAnalyticalVelocity(subcell) = 2
          endif

          thisOCtal%rho(subcell) = max(thisRho, thisOctal%rho(subcell))
       endif
    enddo
  end subroutine assignDensitiesAlphaDisc


  recursive subroutine assignDensitiesMahdavi(grid, thisOctal, astar, mdot, minrCubedRhoSquared)
    use inputs_mod, only :  vturb, isothermTemp, ttauriRstar
    use inputs_mod, only : TTauriDiskHeight
    use magnetic_mod, only : inflowMahdavi, velocityMahdavi
    real(double) :: astar, mdot, thisR, thisRho, thisV, minRcubedRhoSquared
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(VECTOR) :: cellCentre, corner(8)
    integer :: subcell, i, j
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child =>thisOctal%child(i)
                call assignDensitiesMahdavi(grid, child, astar, mdot, minRCubedRhoSquared)
                exit
             end if
          end do
       else
          cellCentre = subcellCentre(thisOctal, subcell)
          thisR = modulus(cellCentre)*1.d10
          call subcellCorners(thisOctal, subcell, corner)
          do j = 1, 8
             corner(j) = 1.d10 * corner(j)
          enddo
          thisOctal%inflow(subcell) = .false.
          if (inflowMahdavi(corner(1:8)).and.inflowMahdavi(1.d10*cellCentre)) then
!          if (inflowMahdavi(1.d10*cellCentre)) then

!          do j = 1, 8
!             write(*,*) j, " modulus ",modulus(corner(j) - (1.d10*cellCentre)) / &
!                  (thisOctal%subcellSize*1.d10),inflowMahdavi(corner(j)), &
!                  returndPhi(thisOctal)*thisR/thisOctal%subcellSize
!          enddo
             
             thisOctal%inFlow(subcell) = .true.
             thisOctal%velocity(subcell) = velocityMahdavi(cellCentre)
             thisOctal%fixedTemperature(subcell) = .true.
             thisV = modulus(thisOctal%velocity(subcell))*cSpeed
             if (thisV /= 0.d0) then
                thisRho =  mdot /(aStar * thisV)  * (ttauriRstar/thisR)**3 
                thisOctal%inflow(subcell) = .true.
                thisOctal%rho(Subcell) = thisRho
                thisOctal%temperature(subcell) = isothermTemp
                if ((thisR**3*thisRho**2) < minRCubedRhoSquared) then
                   minRCubedRhoSquared = thisR**3*thisRho**2
                endif

               if (abs(cellCentre%z) < TTauriDiskHeight/1.d10) then
                   thisOctal%inflow(subcell) = .false.
                   thisrho = 1.d-30
                endif
 
                if (associated(thisOctal%microturb)) thisOctal%microturb(subcell) = vturb
             
                CALL fillVelocityCorners(thisOctal,velocityMahdavi)
 
               if ((subcell==1).and.(thisOCtal%inflow(subcell))) then
                   if ((modulus(thisOctal%cornerVelocity(1)) + &
                       modulus(thisOctal%cornerVelocity(2)) + &
                       modulus(thisOctal%cornerVelocity(3)) + &
                       modulus(thisOctal%cornerVelocity(4)) + &
                       modulus(thisOctal%cornerVelocity(7)) + &
                       modulus(thisOctal%cornerVelocity(8)) + &
                       modulus(thisOctal%cornerVelocity(9)) + &
                       modulus(thisOctal%cornerVelocity(10)) + &
                       modulus(thisOctal%cornerVelocity(2))) < 1.d-20) write(*,*) "corner bug"
                endif               
             endif
          endif
       endif
    enddo
  end subroutine assignDensitiesMahdavi

  recursive subroutine assignTemperaturesMahdavi(grid, thisOctal, astar, mdot, minRcubedRhoSquared)
    use inputs_mod, only : maxHartTemp, isotherm, isothermtemp
    use magnetic_mod, only : inflowMahdavi
    real(double) :: astar, mdot, thisR,  minrCubedRhoSquared
    real(double) :: requiredMaxHeating, thisHeating, localCooling
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(VECTOR) :: cellCentre
    real(double) :: fac
    integer :: subcell, i, j
    real(double), parameter :: logT(8) = (/3.70d0, 3.8d0, 3.9d0, 4.0d0, 4.2d0, 4.6d0, 4.9d0, 5.4d0/)
    real(double), parameter :: gamma(8) = (/-28.3d0, -26.0d0, -24.5d0, -23.6d0, -22.6d0, -21.8d0, &
         -21.2d0, -21.1999d0/)


    if (.not.isoTherm) then
       call locate(logT, 8, log10(dble(maxHartTemp)), j)
       requiredmaxHeating = gamma(j) + (gamma(j+1)-gamma(j))*(log10(dble(maxHartTemp))-logT(j))/(logT(j+1)-logT(j))
       requiredMaxHeating = 10.d0**requiredMaxHeating /mhydrogen**2
       fac = requiredMaxHeating * minRcubedrhoSquared
    endif

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then 
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call assignTemperaturesMahdavi(grid, child, astar, mdot, minRCubedRhoSquared)
                exit
             end if
          end do
       else
          cellCentre = subcellCentre(thisOctal, subcell)
          thisR = modulus(cellCentre)*1.d10
          if (thisOctal%inflow(subcell).and.inflowMahdavi(1.d10*cellcentre)) then
             if (isoTherm) then
                thisOctal%temperature(subcell) = isothermTemp
             else
                thisHeating = fac / thisR**3
                localCooling = log10(thisHeating / (thisOctal%rho(subcell)/mHydrogen)**2)
                call locate(gamma, 8, localCooling, j)
                thisOctal%temperature(subcell) = real(logT(j) + (logT(j+1)-logT(j))*(localCooling - gamma(j))/(gamma(j+1)-gamma(j)))
                thisoctal%temperature(subcell) = max(3.30, min(5.0, thisOctal%temperature(subcell)))
                thisOctal%temperature(subcell) = real(10.d0**thisOctal%temperature(subcell))
             endif
          endif
       endif
    enddo
  end subroutine assignTemperaturesMahdavi
  
  recursive subroutine splitTagged(thisOctal, grid, inheritProps, interpProps)
    use memory_mod, only : globalMemoryFootprint, humanReadableMemory
    use inputs_mod, only : maxMemoryAvailable
    type(GRIDTYPE) :: grid
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  logical, save :: firstTimeMem
  logical :: outOfMemory
  integer :: subcell, i
    logical, optional :: inheritProps, interpProps
    character(len=80) :: message
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call splitTagged(child, grid, inheritProps, interpProps)
                exit
             end if
          end do
       else
          outofmemory = .false.
          if (globalMemoryFootprint > maxMemoryAvailable) then
             outofmemory = .true.
             if (firstTimeMem) then
                write(message,'(a)') "Maxmimum memory exceeded for grid :"//humanReadableMemory(globalMemoryFootprint)
                call writeWarning(message)
                write(message,'(a)') "Not splitting further."
                call writeWarning(message)
                firstTimeMem = .false.
             endif
          endif
          if ((thisOctal%chiline(subcell) > 0.d0).and.(.not.outOfMemory)) then
             thisOctal%chiLine(subcell) = 0.d0
             call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                  inherit=inheritProps, interp=interpProps)
             return
          endif
       endif
    enddo
  end subroutine splitTagged


  recursive subroutine tagScaleSmooth(ntagged, factor, thisOctal, grid,  converged, &
       inheritProps, interpProps)
    type(gridtype) :: grid
    integer :: nTagged
    real :: factor
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    logical, optional :: inheritProps, interpProps
    !
    integer :: subcell, i
    logical :: converged
    real(double) :: r
    type(VECTOR) :: dirVec(6), centre, octVec, aHat
    integer :: neighbourSubcell, j, nDir

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call tagScaleSmooth(nTagged, factor, child, grid, converged, inheritProps, interpProps)
                exit
             end if
          end do
       else

          r = thisOctal%subcellSize/2.d0 + 0.1d0*grid%halfSmallestSubcell
          centre = subcellCentre(thisOctal, subcell) 
          if (thisOctal%threed) then
             if (.not.thisOctal%cylindrical) then
                nDir = 6
                dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
                dirVec(2) = VECTOR( 0.d0,+1.d0,  0.d0)
                dirVec(3) = VECTOR(+1.d0, 0.d0,  0.d0)
                dirVec(4) = VECTOR(-1.d0, 0.d0,  0.d0)
                dirVec(5) = VECTOR( 0.d0,-1.d0,  0.d0)
                dirVec(6) = VECTOR( 0.d0, 0.d0, -1.d0)
             else
                nDir = 4
                aHat = VECTOR(centre%x+0.3*grid%halfSmallestSubcell, centre%y+0.321*grid%halfSmallestSubcell, 0.d0)
                call normalize(aHat)
                dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
                dirVec(2) = aHat
                dirVec(3) = (-1.d0)*aHat
                dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
             endif
          else if (thisOctal%twod) then
             nDir = 4
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)  
             dirVec(2) = VECTOR(-1.d0,0.d0, 0.d0)
             dirVec(3) = VECTOR( 0.d0, 0.d0,  1.d0)
             dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
          else
             nDir = 2
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(-1.d0, 0.d0, 0.d0)
          endif
          do j = 1, nDir
             octVec = centre+VECTOR(0.01*grid%halfSmallestSubcell, 0.01d0*grid%halfSmallestSubcell, &
                  +0.001d0*grid%halfSmallestSubcell) + r * dirvec(j)
             if (inOctal(grid%octreeRoot, octVec)) then

                startOctal => thisOctal
                call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                     foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)

                if ((thisOctal%subcellSize/neighbourOctal%subcellSize) > factor) then
                   thisOctal%chiLine(subcell) = 1.d0
                   converged = .false.
                   nTagged = nTagged + 1
                endif

                if ((neighbourOctal%subcellSize/thisOctal%subcellSize) > factor) then
                   neighbourOctal%chiLine(neighboursubcell) = 1.d0
                   converged = .false.
                   nTagged = nTagged + 1
                endif

             endif
          enddo
       endif
    end do

  end subroutine tagScaleSmooth


  !-----------------------------------------------------------------------------------
  ! Base on startReturnSample of NHS. Optimized for a solving formal solution.
  !    --- (R. Kurosawa)
  ! Returns the amr grid values along a ray in the direction (direction) originated
  ! from "startPoint".
  ! Main differences from startReturnSample is:
  !   1. If a ray encontour "thin disc", starting point will be shifted to the point
  !      of insertersection. 
  !------------------------------------------------------------------------------------
  SUBROUTINE amr_values_along_ray (startPoint,direction,grid,          &
             sampleFreq,nSamples,maxSamples,thin_disc_on, opaqueCore,hitCore, fromDisc, &
             usePops,iLambda,error,lambda,kappaAbs,kappaSca,velocity,&
             velocityDeriv,chiLine,levelPop,rho, temperature, Ne, inflow, &
             etaCont, etaLine)
    ! samples the grid at points along the path.
    ! this should be called by the program, instead of calling 
    !   returnSamples directly, because this checks the start and finish
    !   points are within the grid bounds. returnSamples *assumes* this 
    !   criterion is met. also, the path of the photon is checked for 
    !   intersection with the stellar surface(s), disc etc. 
 
    IMPLICIT NONE

    TYPE(vector), INTENT(INOUT)   :: startPoint ! photon start point
    TYPE(vector), INTENT(IN)      :: direction  ! photon direction 
    TYPE(gridtype), INTENT(IN)         :: grid       ! the entire grid structure
    REAL, INTENT(IN)                   :: sampleFreq ! the maximum number of
!    real(oct), INTENT(IN)   :: sampleFreq ! the maximum number of 
                       ! samples that will be made in any subcell of the octree. 
                       
    INTEGER, INTENT(OUT)               :: nSamples   ! number of samples made
    INTEGER, INTENT(IN)                :: maxSamples ! size of sample arrays 
    logical, intent(in)                :: thin_disc_on   ! T to include thin disc
    LOGICAL, INTENT(IN)                :: opaqueCore ! true if the core is opaque
    LOGICAL, INTENT(OUT)               :: hitCore    ! true if the core is opaque
    LOGICAL, INTENT(OUT)               :: fromDisc   ! starting point the ray was adjusted to the disc
    LOGICAL, INTENT(IN)                :: usePops    ! whether to use level populations
    INTEGER, INTENT(IN)                :: iLambda    ! wavelength index
    INTEGER, INTENT(INOUT)             :: error      ! error code
    REAL, DIMENSION(:), INTENT(INOUT)  :: lambda     ! path distances of sample locations

    
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaAbs   ! continuous absorption opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaSca   ! scattering opacities
    TYPE(vector),DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocity ! sampled velocities
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: velocityDeriv ! sampled velocity derivatives
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: chiLine    ! line opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: rho        ! density at sample points
    real(double),DIMENSION(:,:),INTENT(INOUT),OPTIONAL:: levelPop ! level populations
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: temperature! temperature [K] at sample points
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: Ne ! electron density
    logical,DIMENSION(:),INTENT(INOUT),OPTIONAL  ::inFlow   ! flag to tell if the point is inflow
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaCont ! contiuum emissivity
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaLine ! line emissivity

    TYPE(vector)       :: locator 
                       ! 'locator' is used to indicate a point that lies within the  
                       !   *next* cell of the octree that the ray will interesect.
                       !   initially this will be the same as the startPoint
      
    TYPE(vector)       :: currentPoint ! current position of ray 
    LOGICAL                 :: abortRay     ! flag to signal completion of ray trace
    TYPE(octal)             :: octree       ! the octree structure within 'grid'
    TYPE(vector)       :: directionNormalized
    real(oct)    :: distanceLimit ! max length of ray before aborting
    ! margin is the size of the region around the edge of a subcell
    !   where numerical inaccuracies may cause problems.
    real(oct)    :: margin
 
    ! variables for testing special cases (stellar intersections etc.)
    TYPE(vector)       :: starPosition       ! position vector of stellar centre
    TYPE(vector)       :: diskNormal         ! disk normal vector
    real(oct)    :: rStar              ! stellar radius
    real(oct)    :: endLength          ! max path length of photon
    TYPE(vector)       :: endPoint           ! where photon leaves grid
    LOGICAL                 :: absorbPhoton       ! photon will be absorbed
    real(oct)    :: distanceFromOrigin ! closest distance of plane from origin
    TYPE(vector)       :: diskIntersection   ! point of photon intersection with disk
    real(oct)    :: starIntersectionDistance1 ! distance to first  intersection
    real(oct)    :: starIntersectionDistance2 ! distance to second intersection
    LOGICAL                 :: starIntersectionFound1 ! 1st star intersection takes place      
    LOGICAL                 :: starIntersectionFound2 ! 2nd star intersection takes place      
    LOGICAL                 :: intersectionFound  ! true when intersection takes place      
    real(oct)    :: intersectionRadius ! disk radius when photon intersects
    real(oct)    :: diskDistance       ! distance to disk intersection
    real(oct)    :: distanceThroughStar! distance of chord through star
    TYPE(vector)       :: dummyStartPoint    ! modified start point 
!    real(oct), PARAMETER :: fudgefactor = 1.00001 ! overestimates stellar size
    real(oct), PARAMETER :: fudgefactor = 1.000001 ! overestimates stellar size
    
    
    ! we will abort tracking a photon just before it reaches the edge of the
    !   simulation space. This is the fraction of the total distance to use:
    real(oct), PARAMETER :: distanceFraction = 0.999_oc 

    ! used internally
    TYPE(vector)       ::startPointNew       ! 
!    TYPE(vector)       ::entryPoint       ! 
 
    
    ! we will abort tracking any rays which are too close to 
    !   to a cell wall. The distance to use is defined by:
    margin = 6.0_oc * REAL(grid%maxDepth,kind=oct) * EPSILON(1.0_oc)
    ! some more experimentation is required to find out the best value for
    !   margin.
    
    
    ! set up some variables
    octree = grid%octreeRoot  
    abortRay = .FALSE.
    hitCore = .FALSE.
    fromDisc = .false.
    directionNormalized = direction
    CALL normalize(directionNormalized)
    distanceLimit = HUGE(distanceLimit)
    locator = startPoint
    nSamples = 0
    startPointnew = startPoint

    currentPoint = startPointNew
    
    IF (.NOT. inOctal(octree,startpointnew)) THEN
       
      PRINT *, 'Attempting to find path between point(s) outwith the grid.'
      PRINT *, ' in [amr_mod::amr_values_along_ray].'
      PRINT *, ' ==> StartPoint = (', startPoint, ')'
      error = -30

      !
      print *, "Set the size of the box should be more than the distance between "
      print *, "the center of the star and the edge of the density field."

      stop

    ENDIF
   
   
    ! geometry-specific tests should go here
    IF (grid%geometry(1:6) == "ttauri" .OR. grid%geometry(1:4) == "jets" .or. &
         grid%geometry(1:9) == "luc_cir3d" .or. grid%geometry(1:6) == "cmfgen" .or. &
         grid%geometry(1:8) == "romanova") THEN
      
       ! need to test for both star and disc intersections 
      
       ! we will find out when and where the photon leaves the simulation space 
       CALL getExitPoint(currentPoint,directionNormalized,locator,abortRay,error,&
                         grid%halfSmallestSubcell,endPoint,octree%centre,        &
                         (octree%subcellSize*2.0_oc),endLength,margin,grid%octreeRoot%threed) 

       endLength = endLength * distanceFraction
       
       ! need to reset some of the variables
       abortRay = .FALSE.
       locator = startpointnew
       
       ! for the moment, we will just assume a 2-d disc.
       absorbPhoton = .FALSE.
       diskNormal = grid%diskNormal
       starPosition = grid%starPos1
                                  
       ! find the (geometrycally thin) disk intersection point
       if (grid%geometry == "luc_cir3d" .or. grid%geometry == "cmfgen" ) then
          intersectionFound = .false.
       else
          if (thin_disc_on) then
             distanceFromOrigin = modulus(grid%starPos1)
             diskIntersection = intersectionLinePlane(startpointnew, directionNormalized,&
                  diskNormal, distanceFromOrigin, intersectionFound)
          else
             intersectionFound = .false.
          end if
       end if

       
       IF (intersectionFound) THEN
       
         ! we need to check whether the photon passes through the disk's
         !   central hole, or the outside of the outer radius of accretion disc.
         intersectionRadius =  modulus(diskIntersection - starPosition)
         IF (intersectionRadius > grid%diskRadius .and.  &
              intersectionRadius < 1.5d5) THEN  ! assuming 100 AU disc size

           ! we need to check whether the intersection occurs within our
           !   simulation space.
           diskDistance = modulus(diskIntersection-startpointnew)
           IF (diskDistance < endLength) THEN
              ! Ajdust the current poistion to the insersecting point
              abortRay = .FALSE.
              fromDisc = .true.
              currentPoint = diskIntersection + (fudgefactor*direction)
              startpointnew = currentPoint
              endLength = fudgefactor*(endLength - diskDistance)
              locator = currentPoint                  

           else
              ! no need to do any calculation for this ray.
              abortRay = .TRUE.

           END IF
         END IF
       END IF
                  
       ! now we check for intersections with the star

       rStar = grid%rStar1 
       CALL intersectionLineSphere(startpointnew,directionNormalized,endLength,starPosition, &
                                   rStar,starIntersectionFound1,starIntersectionFound2,   &
                                   starIntersectionDistance1,starIntersectionDistance2)
       ! by passing a line segment to intersectionLineSphere, we ensure that we
       !   do not find intersections with the star that take place after the 
       !   photon has been absorbed by the disk.
       IF (starIntersectionFound1) THEN
       
         endLength = starIntersectionDistance1
         IF (opaqueCore) absorbPhoton = .TRUE.
         
       END IF

       ! we trace the photon path until it encounters the disk, the star,
       !   or the edge of the simulation space.
       error = 0
       CALL returnSamples(currentPoint,startpointnew,locator,directionNormalized,&
                     octree,grid,sampleFreq,nSamples,maxSamples,abortRay,     &
                     lambda,usePops,iLambda,error,margin,endLength,           &
                     kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,   &
                     velocityDeriv=velocityDeriv,chiLine=chiLine,             &
                     levelPop=levelPop,rho=rho, temperature=temperature,      &
                     Ne=Ne, inFlow=inFlow, etaCont=etaCont, etaLine=etaLine)
      
       IF (error < 0)  RETURN

       IF (.NOT. opaqueCore .AND. starIntersectionFound2) THEN 
         ! the photon passes through the star and continues on the other side.

         ! we have to adjust the arguments to returnSamples to make the 
         !   output arrays correct.
         distanceThroughStar = fudgeFactor * &
                      (starIntersectionDistance2 - starIntersectionDistance1)
         currentPoint = currentPoint + distanceThroughStar * directionNormalized
         locator = currentPoint
         dummystartpoint = startpointnew + distanceThroughStar * directionNormalized
         distanceLimit = endLength - distanceThroughStar
           
        CALL returnSamples(currentPoint,dummyStartPoint,locator,             &
                    directionNormalized,octree,grid,sampleFreq,nSamples,     &
                    maxSamples,abortRay,lambda,usePops,iLambda,error,margin, &
                    distanceLimit,kappaAbs=kappaAbs,kappaSca=kappaSca,       &
                    velocity=velocity,velocityDeriv=velocityDeriv,           &
                    chiLine=chiLine,levelPop=levelPop,rho=rho,               &
                    temperature=temperature, Ne=Ne, inFlow=inFlow,           &
                    etaCont=etaCont, etaLine=etaLine)

       END IF

       ! if the photon ends up in the star, we make sure it absorbed.
       IF (absorbPhoton) THEN
       
         nSamples = nSamples + 1
         IF (nSamples > maxSamples) THEN
           PRINT *, "Error:: nSamples > maxSamples in [amr_mod::amr_values_along_ray] "
           PRINT *, "        nSamples   = ", nSamples
           PRINT *, "        maxSamples = ", maxSamples
           STOP
         END IF
         
         lambda(nSamples) = lambda(nSamples-1)
         hitCore = .TRUE.
         IF (PRESENT(kappaSca))      kappaSca(nSamples) = 0.
         IF (PRESENT(kappaAbs))      kappaAbs(1:nSamples) = 1.e20
         IF (PRESENT(velocity))      velocity(nSamples) = vector(0.,0.,0.)
         IF (PRESENT(velocityDeriv)) velocityDeriv(nSamples) = 1.0
         IF (PRESENT(chiLine))       chiLine(1:nSamples) = 1.e20
         IF (PRESENT(levelPop))      levelPop(nSamples,:) = 0.0
         IF (PRESENT(Ne))            Ne(nSamples) = 0.0
         IF (PRESENT(rho))           rho(nSamples) = 0.0
         IF (PRESENT(temperature))   temperature(nSamples) = 0.0
         IF (PRESENT(inFlow))        inFlow(nSamples) = .true.
         IF (PRESENT(etaCont))       etaCont(nSamples) = 0.0
         IF (PRESENT(etaLine))       etaLine(nSamples) = 0.0
       END IF
       
    ELSE
            
       CALL getExitPoint(currentPoint,directionNormalized,locator,abortRay,error,&
                         grid%halfSmallestSubcell,endPoint,octree%centre,        &
                         (octree%subcellSize*2.0_oc),endLength,margin,grid%octreeRoot%threed) 
       distanceLimit = endLength * distanceFraction
       
       ! need to reset some of the variables
       abortRay = .FALSE.
       locator = startpointnew
       
       CALL returnSamples(currentPoint,startpointnew,locator,directionNormalized,&
                   octree,grid,sampleFreq,nSamples,maxSamples,abortRay,lambda,&
                   usePops,iLambda,error,margin,distanceLimit,                &
                   kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,     &
                   velocityDeriv=velocityDeriv,chiLine=chiLine,               &
                   levelPop=levelPop,rho=rho,temperature=temperature,         &
                   Ne=Ne,inFlow=inFlow,etaCont=etaCont,etaLine=etaLine)
    END IF

    ! The new starting point will be returned to a parent routine.
    startpoint = startpointnew
      
  end SUBROUTINE amr_values_along_ray


  !
  ! What is the difference between this and the original startReturnSamples routine??????
  ! 
  SUBROUTINE startReturnSamples2(startPoint,direction,grid,          &
             nSamples,maxSamples,thin_disc_on, opaqueCore,hitCore,      &
             usePops,iLambda,error,lambda,nSource,source,kappaAbs,kappaSca,velocity,&
             velocityDeriv,chiLine,levelPop,rho, temperature, Ne, inflow, etaCont, etaLine, startOctal, startSubcell)
    ! samples the grid at points along the path.
    ! this should be called by the program, instead of calling 
    !   returnSamples directly, because this checks the start and finish
    !   points are within the grid bounds. returnSamples *assumes* this 
    !   criterion is met. also, the path of the photon is checked for 
    !   intersection with the stellar surface(s), disc etc. 
 
    USE magField, only: innerDiskData_Radii, innerDiskData_Phi, maxSizeMagFieldGrid
    USE source_mod, only: SOURCETYPE, distanceToSource

    IMPLICIT NONE

    type(octal), pointer, optional :: startOctal
    integer, optional :: startSubcell
    TYPE(vector), INTENT(IN)      :: startPoint ! photon start point
    TYPE(vector), INTENT(IN)      :: direction  ! photon direction 
    TYPE(gridtype), INTENT(IN)         :: grid       ! the entire grid structure
!    real(oct), INTENT(IN)   :: sampleFreq ! the maximum number of 
                       ! samples that will be made in any subcell of the octree. 
                       
    INTEGER, INTENT(OUT)             :: nSamples   ! number of samples made
    INTEGER, INTENT(IN)                :: maxSamples ! size of sample arrays 
    logical, intent(in)                :: thin_disc_on   ! T to include thin disc
    LOGICAL, INTENT(IN)                :: opaqueCore ! true if the core is opaque
    LOGICAL, INTENT(OUT)               :: hitCore    ! true if the core is opaque
    LOGICAL, INTENT(IN)                :: usePops    ! whether to use level populations
    INTEGER, INTENT(IN)                :: iLambda    ! wavelength index
    INTEGER, INTENT(INOUT)             :: error      ! error code
    REAL, DIMENSION(:), INTENT(INOUT)  :: lambda     ! path distances of sample locations
    integer, intent(in), optional      :: nSource
    type(SOURCETYPE),intent(in), optional :: source(:)
    
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaAbs   ! continuous absorption opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaSca   ! scattering opacities
    TYPE(vector),DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocity ! sampled velocities
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: velocityDeriv ! sampled velocity derivatives
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: chiLine    ! line opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: rho        ! density at sample points
    real(double),DIMENSION(:,:),INTENT(INOUT),OPTIONAL:: levelPop ! level populations
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: temperature! temperature [K] at sample points
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: Ne ! electron density
    logical,DIMENSION(:),INTENT(INOUT),OPTIONAL  ::inFlow   ! flag to tell if the point is inflow
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaCont ! contiuum emissivity
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaLine ! line emissivity

    TYPE(vector)       :: locator 
                       ! 'locator' is used to indicate a point that lies within the  
                       !   *next* cell of the octree that the ray will interesect.
                       !   initially this will be the same as the startPoint
      
    LOGICAL                 :: abortRay     ! flag to signal completion of ray trace
    TYPE(vector)       :: directionNormalized
 
    ! variables for testing special cases (stellar intersections etc.)
    TYPE(vector)       :: starPosition       ! position vector of stellar centre
    TYPE(vector)       :: diskNormal         ! disk normal vector
    real(oct)    :: rStar              ! stellar radius
    real(oct)    :: endLength          ! max path length of photon
    LOGICAL                 :: absorbPhoton       ! photon will be absorbed
    real(oct)    :: distanceFromOrigin ! closest distance of plane from origin
    TYPE(vector)       :: diskIntersection   ! point of photon intersection with disk
    real(oct)    :: starIntersectionDistance1 ! distance to first  intersection
    real(oct)    :: starIntersectionDistance2 ! distance to second intersection
    LOGICAL                 :: starIntersectionFound1 ! 1st star intersection takes place      
    LOGICAL                 :: starIntersectionFound2 ! 2nd star intersection takes place      
    LOGICAL                 :: intersectionFound  ! true when intersection takes place      
    real(oct)    :: intersectionRadius ! disk radius when photon intersects
    real(oct)    :: diskDistance       ! distance to disk intersection
    
    type(vector) :: currentPosition
    type(octal), pointer :: sOctal, thisOctal
    integer :: subcell
    real(double) :: length, distToNextCell
    ! we will abort tracking a photon just before it reaches the edge of the
    !   simulation space. This is the fraction of the total distance to use:
    real(double) :: fudgeFac = 0.1d0
    integer :: nSource_local
    TYPE(vector)       :: vectorToIntersectionXYZ ! intersection with disk (cartesian)
    TYPE(spVector)    :: vectorToIntersectionSP ! intersection with disk (sph. polar)
    REAL(KIND=oct)    :: intersectionPhi ! disk azimuth angle where photon intersects
    INTEGER                 :: lowerBound ! index returned from 'hunt' with lower bound of phi bin
    INTEGER                 :: upperBound ! usually lowerBound+1
    REAL(KIND=oct)    :: localDiskRadius ! disk inner radius at current azimuth (1.e10cm)
    REAL(KIND=oct)    :: lowerDistance ! angular distance to lower phi bin
    REAL(KIND=oct)    :: upperDistance ! angular distance to upper phi bin
    LOGICAL           :: boundaryProblem = .false.

    currentPosition = startPoint
    locator = startPoint
    nSamples = 0
    length = 0.d0
    endLength = 1.d30
    directionNormalized = direction
    CALL normalize(directionNormalized)

    absorbPhoton = .false.
    hitcore = .false.

    ! Since not all gemoetry uses a "source" object, you need this 
    ! statement (RK).
    if (PRESENT(nSource)) then
       nSource_local = nSource
    else
       nSource_local = 0
    end if

    if (nSource_local == 0) then
       ! first we check for intersections with the star
       starPosition = grid%starPos1
       rStar = grid%rStar1 
       CALL intersectionLineSphere(startPoint,direction,1.d30,starPosition, &
            rStar,starIntersectionFound1,starIntersectionFound2,   &
            starIntersectionDistance1,starIntersectionDistance2)
       ! by passing a line segment to intersectionLineSphere, we ensure that we
       !   do not find intersections with the star that take place after the 
       !   photon has been absorbed by the disk.
       IF (starIntersectionFound1) THEN
          endLength = starIntersectionDistance1
          IF (opaqueCore) absorbPhoton = .TRUE.
       END IF
       IF (.NOT. opaqueCore .AND. starIntersectionFound2) then
          absorbPhoton = .false.
          endLength = 1.d30
       endif
    else
       call distanceToSource(source, nSource_local, startPoint, direction, absorbPhoton, endLength)
    endif


    !
    ! geometry-specific tests should go here
    !
    IF (grid%geometry(1:6) == "ttauri" .OR. grid%geometry(1:4) == "jets" .or. &
         grid%geometry(1:8) == "romanova") THEN
      
       ! need to reset some of the variables
       abortRay = .FALSE.
       locator = startPoint
       
       ! for the moment, we will just assume a 2-d disc.
       absorbPhoton = .FALSE.
       diskNormal = grid%diskNormal
       starPosition = grid%starPos1
                                  
       ! find the (geometrycally thin) disk intersection point
       if (thin_disc_on) then
          distanceFromOrigin = modulus(grid%starPos1)
          diskIntersection = intersectionLinePlane(startPoint, directionNormalized,&
               diskNormal, distanceFromOrigin, intersectionFound)
       else
          intersectionFound = .false.
       end if
       
       IF (intersectionFound) THEN
       
         ! we need to check whether the photon passes through the disk's
         !   central hole, or the outside of the outer radius of accretion disc.
         intersectionRadius =  modulus(diskIntersection - starPosition)

         IF (.NOT. grid%geometry == "magstream") THEN
           localDiskRadius = grid%diskRadius
         ELSE ! grid%geometry == "magstream"
           vectorToIntersectionXYZ = diskIntersection - starPosition
           ! convert to spherical polars
            call XYZtoSPvector(vectorToIntersectionSP, vectorToIntersectionXYZ)
           intersectionPhi = vectorToIntersectionSP%phi
           
           ! need to find the appropriate disk truncation radius at this azimuth
           
           IF (.NOT. ALLOCATED(innerDiskData_Phi) .OR. &
               .NOT. ALLOCATED(innerDiskData_Radii)) THEN
             PRINT *, "Panic in startReturnSamples: innerDiskData not ALLOCATED"  
             STOP
           END IF

           CALL locate(innerDiskData_Phi,SIZE(innerDiskData_Phi),intersectionPhi,lowerBound)
           
           ! need to check for special cases - close to 0 / 2pi
           IF (lowerBound == 0) THEN
             upperBound = 1
             lowerBound = SIZE(innerDiskData_Phi)
             ! find out whether phi is closer to the lower or upper value
             lowerDistance = intersectionPhi - (innerDiskData_Phi(lowerBound) - twoPi)
             upperDistance = innerDiskData_Phi(upperBound) - intersectionPhi
           ELSE IF (lowerBound == SIZE(innerDiskData_Phi)) THEN
             upperBound = 1
             lowerBound = SIZE(innerDiskData_Phi)
             lowerDistance = intersectionPhi - innerDiskData_Phi(lowerBound)
             upperDistance = (innerDiskData_Phi(upperBound) + twoPi) - intersectionPhi
           ELSE
             ! not a special case
             upperBound = lowerBound + 1
             lowerDistance = intersectionPhi - innerDiskData_Phi(lowerBound)
             upperDistance = innerDiskData_Phi(upperBound) - intersectionPhi
           END IF

           ! if there is an angular distance of more than 5 degrees to the 
           !   nearest sample point, we will assume we are in a region
           !   with no accretion streams, and adopt the maximum disk radius
           IF ( MIN(lowerDistance,upperDistance) > (5.0_oct * degToRad) ) THEN
             localDiskRadius = maxSizeMagFieldGrid
           ! else, we have a close sample
           ELSE IF (lowerDistance <= upperDistance) THEN
             localDiskRadius = innerDiskData_Radii(lowerBound)
           ELSE
             localDiskRadius = innerDiskData_Radii(upperBound)
           END IF
           
         END IF

         IF (intersectionRadius > localDiskRadius) THEN  ! assuming 100 AU disc size
!              intersectionRadius < grid%octreeRoot%subcellsize*2.0) THEN
           absorbPhoton = .TRUE.

           ! we need to check whether the intersection occurs within our
           !   simulation space.
           diskDistance = modulus(diskIntersection-startPoint)
           ! now compare this distance with the endLength computed
           ! earlier. 
           endLength = min(endLength, diskDistance)
         END IF
       END IF
       
    end IF

    if (present(startOctal)) then
       thisOctal => startOctal
       subcell = startSubcell
    else
       CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)
    endif

    do while (inOctal(grid%octreeRoot, currentPosition).and.(length < endLength))

       call findSubcellLocal(currentPosition,thisOctal,subcell,boundaryProblem)
       if (boundaryProblem) then
         boundaryProblem = .false.
         call findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)
       end if

       sOctal => thisOctal
       call distanceToCellBoundary(grid, currentPosition, direction, DisttoNextCell, sOctal)

       if (nSamples < maxSamples) then
          call takeSample(currentPosition,length,direction,grid,thisOctal,subcell,nSamples,&
               maxSamples,usePops,iLambda,error,lambda,&
               kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity, velocityDeriv=velocityDeriv, &
               chiLine=chiLine,levelPop=levelPop,rho=rho,  &
               temperature=temperature,Ne=Ne,inFlow=inFlow, etaCont=etaCont, etaLine=etaLine) 
       else
          call writeWarning("Reached maxSamples limit in ray trace")
          exit
       endif
  
       currentPosition = currentPosition + (distToNextCell+fudgeFac*grid%halfSmallestSubcell)*direction
       length = length + distToNextCell+fudgeFac*grid%halfSmallestSubcell
    end do


    ! if the photon ends up in the disk or the star, we make sure it absorbed.
    IF (absorbPhoton) THEN
       
       nSamples = nSamples + 1
       IF (nSamples > maxSamples) THEN
          PRINT *, "Error:: nSamples > maxSamples in startReturnSamples subroutine"
          PRINT *, "        nSamples   = ", nSamples
          PRINT *, "        maxSamples = ", maxSamples
          STOP
       END IF
       
       lambda(nSamples) = real(endLength)
       hitCore = .TRUE.
       IF (PRESENT(kappaSca))      kappaSca(nSamples) = 0.
       IF (PRESENT(kappaAbs))      kappaAbs(nSamples) = 1.e20
       IF (PRESENT(velocity))      velocity(nSamples) = vector(0.,0.,0.)
       IF (PRESENT(velocityDeriv)) velocityDeriv(nSamples) = 1.0
       IF (PRESENT(chiLine))       chiLine(nSamples) = 1.e20
       IF (PRESENT(levelPop))      levelPop(nSamples,:) = 0.0
       IF (PRESENT(Ne))            Ne(nSamples) = 0.0
       IF (PRESENT(rho))           rho(nSamples) = 0.0
       IF (PRESENT(temperature))   temperature(nSamples) = 0.0
       IF (PRESENT(inFlow))        inFlow(nSamples) = .true.
       IF (PRESENT(etaLine))       etaLine(nSamples) = 0.0d0
       IF (PRESENT(etaCont))       etaCont(nSamples) = 0.0d0

    END IF

  end SUBROUTINE startReturnSamples2

  recursive subroutine getxValues(thisOctal, nx, xAxis)

    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    type(VECTOR) :: rVec
    integer :: nx, subcell, i
    real(double) :: xAxis(:)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call getxValues(child, nx, xAxis)
                exit
             end if
          end do
       else

          rVec = subcellCentre(thisOctal, subcell)
          nx = nx + 1
          xAxis(nx) = rVec%x
       end if
    end do

  end subroutine getxValues

  !
  ! hydroWarp stuff
  !
  function hydroWarpTemperature(point, grid) result (tempOut)
    use inputs_mod, only : warpFracHeight, warpRadius, warpSigma, warpAngle
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector), INTENT(IN) :: point
    real(double) :: tempOut
    TYPE(vector) :: rVec
    real(double) :: r, phi, warpHeight
    real :: thisTemp

    tempOut = 20.

    rVec = point
    r = sqrt(rVec%x**2 + rVec%y**2)
    phi = atan2(rVec%y,rVec%x)
    warpheight  = cos(phi+warpAngle) * warpFracHeight * warpradius * exp(-0.5d0*((r - warpRadius)/warpSigma)**2)
    rVec%z = rVec%z - warpheight

    if (inOctal(grid%hydroGrid%octreeRoot, rVec)) then
      call amrGridValues(grid%hydroGrid%octreeRoot, rVec, temperature=thisTemp)
      tempOut = thisTemp
    end if
  end function hydroWarpTemperature

  function hydroWarpDensity(point, grid) result (rhoOut)
    use inputs_mod, only : warpFracHeight, warpRadius, warpSigma, warpAngle
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector), INTENT(IN) :: point
    real(double) :: rhoOut
    TYPE(vector) :: rVec
    real(double) :: r, phi, warpHeight

    integer :: nx, thisX
    type(hydroSpline) :: thisSpline

    rhoOut = 1.d-30

    rVec = point
    r = sqrt(rVec%x**2 + rVec%y**2)
    phi = atan2(rVec%y,rVec%x)
    warpheight  = cos(phi+warpAngle) * warpFracHeight * warpradius * exp(-0.5d0*((r - warpRadius)/warpSigma)**2)
    rVec%z = rVec%z - warpheight

    ! Determine which x-value most closely matches our r-value
    nx = size(grid%hydroGrid%hydroSplines)
    call locate_hydroSpline(grid%hydroGrid,nx,r,thisX)
    if (thisX < 1) then
      thisX = 1
    else if (thisX < nx) then
      if (abs(grid%hydroGrid%hydroSplines(thisX+1)%x-r) .lt. abs(r-grid%hydroGrid%hydroSplines(thisX)%x)) thisX = thisX + 1
    end if

    ! Interpolate over the appropriate spline to get the log(rho) at the point
    thisSpline = grid%hydroGrid%hydroSplines(thisX)
    call splint(thisSpline%z, thisSpline%rho, thisSpline%rhoDD, thisSpline%nz, rVec%z, rhoOut)

    ! The spline function stores and interpolates on the log of the density, so we
    ! must unlog the densities to give our result.
    rhoOut = 10**rhoOut
  end function hydroWarpDensity

  function getHydroWarpScaleHeight(r, grid) result (scaleheight)
    type(gridtype), intent(in) :: grid
    real(double), intent(in) :: r
    real(double) :: scaleheight

    integer :: nx, thisX

    ! Determine which x-value most closely matches our r-value
    nx = size(grid%hydroGrid%hydroSplines)
    call locate_hydroSpline(grid%hydroGrid,nx,r,thisX)
    if (thisX < 1) then
      thisX = 1
    else if (thisX < nx) then
      if (abs(grid%hydroGrid%hydroSplines(thisX+1)%x-r) .lt. abs(r-grid%hydroGrid%hydroSplines(thisX)%x)) thisX = thisX + 1
    end if

    scaleheight = grid%hydroGrid%hydroSplines(thisX)%scaleheight
  end function getHydroWarpScaleHeight


    SUBROUTINE LOCATE_hydroSpline(grid,N,X,J)
    type(gridtype), intent(in) :: grid
    integer, intent(in)              :: n
    real(double), intent(in) :: x
    integer,intent(out)              :: j
    integer :: jl, ju,jm
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((grid%hydroSplines(N)%x.GT.grid%hydroSplines(1)%x).EQV.(X.GE.grid%hydroSplines(JM)%x))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF

      ! Will force to be between 1 and the array size -1.
      if(x <= grid%hydroSplines(1)%x)then
        j=1
      else if(x>=grid%hydroSplines(n)%x)then
        j=n-1
      else
        j=jl
      end if

      return

    END SUBROUTINE LOCATE_hydroSpline

  subroutine hydroWarpScaleHeight(spline)
    type(hydroSpline) :: spline
    real(double) :: max_rho, max_rho_z
    integer :: max_rho_i
    real(double) :: target_rho

    integer :: i

    integer :: sv1_i, sv2_i
    real(double) :: sv1, sv2, step1, step2
    real(double) :: sh1, sh2

    call hydroWarpSplinePeak(spline, max_rho, max_rho_z, max_rho_i)
    target_rho = max_rho - log10(2.718281828)

! improve this conditional
    if (max_rho .lt. log10(1.d-29)) then
      ! We have a straight line
      spline%scaleHeight = 0.
      return
    endif

    ! assume the spline is fairly stable around one scaleheight from the peak
    ! find the cell before which we pass through the target density
    sv1_i = 0
    sv2_i = 0
    do i=1, max_rho_i
      if (spline%rho(i) .gt. target_rho) then
        if (i == 1) then
          print *, "Error finding scaleheight: sv1_i = 1!"
          stop
        end if
        sv1_i = i-1
        sv1 = spline%z(sv1_i)
        step1 = abs(spline%z(sv1_i+1)-sv1)/2.
        exit
      end if
    end do
    if (sv1_i == 0) then
      sv1_i = max_rho_i
      sv1 = spline%z(sv1_i)
      step1 = abs(max_rho_z-sv1)/2.
    end if

    do i=spline%nz, max_rho_i+1, -1
      if (spline%rho(i) .gt. target_rho) then
        if (i == spline%nz) then
          print *, "Error finding scaleheight: sv2_i = nz!"
          stop
        end if
        sv2_i = i+1
        sv2 = spline%z(sv2_i)
        step2 = -1.*abs(sv2-spline%z(sv2_i-1))/2.
        exit
      end if
    end do
    if (sv2_i == 0) then
      sv2_i = max_rho_i+1
      sv2 = spline%z(sv2_i)
      step2 = -1.*abs(sv2-max_rho_z)/2.
    end if

    sh1 = hydroWarpSplinty(spline, target_rho, sv1, step1)
    sh2 = hydroWarpSplinty(spline, target_rho, sv2, step2)

    spline%scaleHeight = abs(sh2-sh1)/2.

    ! perform a quick sanity check
    if (abs(sh1+sh2-2.*max_rho_z)/spline%scaleHeight .gt. 0.05) then
      print *, "Warning: Scaleheights either side of midplane differ by"
      print *, "         more than 5%:", sh1, sh2
    end if
  end subroutine hydroWarpScaleHeight

  real(double) function hydroWarpSplinty(spline, target, startval, step)
    implicit none

    type(hydroSpline) :: spline
    real(double) :: target    ! target density value
    real(double) :: startval  ! initial value of z to try
    real(double) :: step      ! initial stepping in value of z

    ! initial values
    real(double) :: tol = 0.0001   ! +/- tol * target
    real(double) :: reduxFac = 0.5 ! step is reduced by this factor when homing in
    integer :: maxIter = 100       ! maximum number of iterations before we give up
    integer :: i                   ! iteration counter
   
    real(double) :: xNew, yNew, xOld, yOld

    xNew = startval
    call splint(spline%z, spline%rho, spline%rhoDD, spline%nz, xNew, yNew)
    i = 0
    do
      i = i + 1
      if (i > maxIter) then
         write (*,*) "Numerical solver: exceeded maximum iterations allowed (", maxIter, ")"
         stop
      end if

      xOld = xNew
      yOld = yNew
      xNew = xNew + step
      call splint(spline%z, spline%rho, spline%rhoDD, spline%nz, xNew, yNew)
      if (abs(yNew-target) < abs(tol * target)) then
         exit
      else if (yNew < target) then
         if (yNew < yOld) then
            if (yOld < target) then
               step = -1. * step
               xNew = xOld
               yNew = yOld
            else ! (yOld > target)
               step = -1. * (reduxFac * step)
            end if
         end if
      else ! (yNew > target)
         if (yNew > yOld) then
            if (yOld < target) then
               step = -1. * (reduxFac * step)
            else ! (yOld > target)
               step = -1. * step
               xNew = xOld
               yNew = yOld
            end if
         end if
      end if
    end do

    hydroWarpSplinty = xNew
  end function hydroWarpSplinty

  ! returns, the max density, z value of max density, and index of element
  ! immediately preceeding the max density location
  subroutine hydroWarpSplinePeak(spline, rho_max, rho_max_z, rho_max_i)
    type(hydroSpline), intent(in) :: spline
    real(double), intent(out) :: rho_max, rho_max_z
    integer, intent(out) :: rho_max_i

    integer :: i

    real(double) :: step, z_test, rho_test
    integer :: nsteps
    nsteps = 100

    rho_test = 0.d0
    ! desnity is log(10) density, so we need a very low value here
    rho_max = -350.
    rho_max_i = 0

    ! find cell with the max density
    do i=1, spline%nz
      if (spline%rho(i) .gt. rho_max) then
        rho_max = spline%rho(i)
        rho_max_i = i
      end if
    end do
    rho_max_z = spline%z(rho_max_i)
    
    ! the peak in the spline should be somewhere between this max cell value
    ! and the two neighbouring cells
    if (rho_max_i == 1) then
      step = (spline%z(2) - spline%z(1)) / (nsteps-1)
      z_test = spline%z(1)
    else if (rho_max_i == spline%nz) then
      step = (spline%z(rho_max_i) - spline%z(rho_max_i-1)) / (nsteps-1)
      z_test = spline%z(rho_max_i-1)
    else
      step = (spline%z(rho_max_i+1) - spline%z(rho_max_i-1)) / (nsteps-1)
      z_test = spline%z(rho_max_i-1)
    end if

    do i=1, nsteps
        call splint(spline%z, spline%rho, spline%rhoDD, spline%nz, z_test, rho_test)
        if (rho_test .gt. rho_max) then
          rho_max = rho_test
          rho_max_z = z_test
        end if
        z_test = z_test + step
    end do

    ! set the index for the density peak to the index of the element
    ! immediately preceeding the peak value
    if (rho_max_z < spline%z(rho_max_i)) rho_max_i = rho_max_i - 1
  end subroutine hydroWarpSplinePeak
  !
  ! end hydroWarp stuff
  !

  subroutine allocateStream(thisStream, nSamples)
    type(STREAMTYPE) :: thisStream
    integer :: nSamples
!    call freeStream(thisStream)
    allocate(thisStream%position(1:nSamples))
    allocate(thisStream%direction(1:nSamples))
    allocate(thisStream%velocity(1:nSamples))
    allocate(thisStream%rho(1:nSamples))
    allocate(thisStream%temperature(1:nSamples))
    allocate(thisStream%speed(1:nSamples))
    allocate(thisStream%distanceAlongStream(1:nSamples))
    allocate(thisStream%streamRadius(1:nSamples))

  end subroutine allocateStream

  subroutine freeStream(thisStream)
    type(STREAMTYPE) :: thisStream
    if (associated(thisStream%position)) deallocate(thisStream%position)
    if (associated(thisStream%velocity)) deallocate(thisStream%velocity)
    if (associated(thisStream%rho)) deallocate(thisStream%rho)
    if (associated(thisStream%temperature)) deallocate(thisStream%temperature)

    if (associated(thisStream%direction)) deallocate(thisStream%direction)
    if (associated(thisStream%speed)) deallocate(thisStream%speed)
    if (associated(thisStream%streamRadius)) deallocate(thisStream%streamRadius)
    if (associated(thisStream%distanceAlongStream)) deallocate(thisStream%distanceAlongStream)

  end subroutine freeStream

  subroutine createOctalStream(thisOctal,  thisStream, octalStream)
    type(OCTAL),pointer :: thisOctal
    type(STREAMTYPE) :: thisStream, octalStream
    integer :: i

    call allocateStream(octalStream, thisStream%nSamples)

    octalStream%nSamples = 0
    do i = 1, thisStream%nSamples
       if (inOctal(thisOctal, thisStream%position(i))) then
          octalStream%nSamples = octalStream%nSamples + 1
          octalStream%rho(octalStream%nSamples) = thisStream%rho(i)
          octalStream%position(octalStream%nSamples) = thisStream%position(i)
          octalStream%velocity(octalStream%nSamples) = thisStream%velocity(i)
          octalStream%temperature(octalStream%nSamples) = thisStream%temperature(i)
       endif
    end do
    write(*,*) "Octal stream created with ",octalStream%nSamples, " samples"
  end subroutine createOctalStream

  subroutine createSubcellStream(thisOctal,  subcell, thisStream, octalStream)
    integer :: subcell
    type(OCTAL),pointer :: thisOctal
    type(STREAMTYPE) :: thisStream, octalStream
    integer :: i
    call allocateStream(octalStream, thisStream%nSamples)
    octalStream%nSamples = 0 
    do i = 1, thisStream%nSamples
       if (inSubcell(thisOctal, subcell, thisStream%position(i))) then
          octalStream%nSamples = octalStream%nSamples + 1
          octalStream%rho(octalStream%nSamples) = thisStream%rho(i)
          octalStream%position(octalStream%nSamples) = thisStream%position(i)
          octalStream%velocity(octalStream%nSamples) = thisStream%velocity(i)
          octalStream%temperature(octalStream%nSamples) = thisStream%temperature(i)
       endif
    end do
  end subroutine createSubcellStream

  subroutine tauAlongPath(ilambda, grid, rVec, direction, tau, tauMax, ross, startOctal, startSubcell, nTau, &
       xArray, tauArray, distanceToEdge, subRadius, stopatGap, stopatdistance)
    use inputs_mod, only : rGap, rGapInner, rGapOuter
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec, direction, currentPosition, beforeVec, afterVec
    real(double), optional,intent(out) :: xArray(:), tauArray(:)
    real(double), optional :: stopAtDistance
    integer, optional, intent(out) :: nTau
    logical, optional :: stopAtGap
    integer :: iLambda
    real(double), optional :: subRadius
    real(double), intent(out) :: tau
    real(double), optional, intent(out) :: distanceToEdge
    real(double) :: distToNextCell
    real(double), optional :: tauMax
    type(OCTAL), pointer :: thisOctal, sOctal
    type(OCTAL), pointer, optional :: startOctal
    integer, optional :: startSubcell
    integer :: nArray
    real(double) :: fudgeFac = 1.d-1
    real(double) :: kappaSca, kappaAbs, kappaExt
    real(double) :: r, rStart
    logical :: outwards, inwards
    integer :: subcell
    logical, optional :: ross
    logical :: planetGap

! Set subradius to missing data in case tau > tauMax condition is never triggered
    if (present(subradius)) then 
       subradius=-1.0d30
    end if

    kappaAbs = 0.d0; kappaSca = 0.d0
    tau = 0.d0
    currentPosition = rVec
    if (PRESENT(nTau)) then
       xArray = 0.d0
       tauArray = 0.d0
       ntau = 1
       nArray = size(tauArray)
    endif
    planetGap  = .false.
    if (grid%geometry == "planetgap") planetgap = .true.
    rStart = modulus(rVec)
    if ((rVec .dot. direction) > 0.d0) then
       outwards = .true.
       inwards = .false.
    else
       inwards = .true.
       outwards=.false.
    endif

    if (PRESENT(startOctal)) then
       thisOctal => startOctal
       subcell = startSubcell
!       if (.not.inOctal(thisOctal, currentPosition)) then
!          write(*,*) "bug in startreturnsamples2"
!          stop
!       endif
    else
       CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)
    endif

    if (PRESENT(distanceToEdge)) distanceToEdge = 0.d0
    do while (inOctal(grid%octreeRoot, currentPosition))

       call findSubcellLocal(currentPosition, thisOctal,subcell)
       if (.not.PRESENT(ross)) then
          call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, kappaSca=kappaSca, kappaAbs=kappaAbs)
          kappaExt = kappaAbs + kappaSca
       else
          call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, rosselandKappa=kappaExt)
          kappaExt = kappaExt * thisOctal%rho(subcell) * 1.d10
       endif
       sOctal => thisOctal
       call distanceToCellBoundary(grid, currentPosition, direction, DisttoNextCell, sOctal)
       
       beforeVec = VECTOR(currentPosition%x, currentPosition%y,0.d0)
       currentPosition = currentPosition + (distToNextCell+fudgeFac*grid%halfSmallestSubcell)*direction
       afterVec = VECTOR(currentPosition%x, currentPosition%y,0.d0)
       
       if (thisOctal%twod.and.(direction%x < 0.d0).and.(currentPosition%x < 0.d0)) exit

       if ((thisOctal%cylindrical).and.((beforeVec.dot.afterVec).lt.0.d0)) then
          exit
       endif
          

       if (planetGap) then
          if ((direction%x < 0.d0).and.(rVec%x > rGap*autocm/1.d10) &
               .and.(currentPosition%x < rGap*autocm/1.d10)) exit
          if ((direction%x > 0.d0).and.(rVec%x < rGap*autocm/1.d10) &
               .and.(currentPosition%x > rGap*autocm/1.d10)) exit
       endif

       if (PRESENT(stopAtGap)) then
          r = modulus(currentPosition)
          if (inwards.and.(rStart > rGapOuter).and.(r < rGapOuter)) exit
          if (outwards.and.(rStart < rGapInner).and.(r > rGapInner)) exit
       endif

       if (PRESENT(stopAtDistance)) then
          r = modulus(currentPosition)
          if (r > stopAtDistance) exit
       endif


       tau = tau + distToNextCell*kappaExt
       if (PRESENT(nTau)) then
          nTau = nTau + 1
          if (nTau > nArray) then
             call writeFatal("Tau array size exceeded")
             stop
          endif
          xArray(nTau) = xArray(nTau-1) + distToNextCell
          tauArray(nTau) = tau
       endif
       if (PRESENT(distanceToEdge)) distanceToEdge = distanceToEdge + distToNextCell
       if (PRESENT(tauMax)) then
          if (tau > tauMax) then
             if (PRESENT(subRadius)) then
                subRadius = modulus(subcellCentre(thisOctal,subcell))
             endif
             exit
          endif
          
       endif
    end do
  end subroutine tauAlongPath

  subroutine tauAlongPathFast(ilambda, grid, rVec, direction, tau, tauMax, ross, startOctal, startSubcell, nTau, &
       xArray, tauArray, distanceToEdge, debug)
    use source_mod, only : globalSourceArray, globalNsource, distanceToSource
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec, direction, currentPosition
    logical, optional :: debug
    real, optional,intent(out) :: xArray(:), tauArray(:)
    integer, optional, intent(out) :: nTau
    integer :: iLambda
    real(double), intent(out) :: tau
    real(double), optional, intent(out) :: distanceToEdge
    real(double) :: distToNextCell
    real(double), optional :: tauMax
    type(OCTAL), pointer :: thisOctal, sOctal
    type(OCTAL), pointer, optional :: startOctal
    integer, optional :: startSubcell
    real(double), parameter :: fudgeFac = 1.d-1
    real(double) :: kappaSca, kappaAbs, kappaExt, distance
    integer :: subcell
    logical, optional :: ross
    logical :: hitGrid, hitSource
    integer :: sourceNumber
    kappaAbs = 0.d0; kappaSca = 0.d0
    tau = 0.d0
    currentPosition = rVec


    if (PRESENT(nTau)) then
       xArray(1) = 0.d0
       tauArray(1) = 0.d0
       ntau = 1
    endif

    if (PRESENT(debug)) write(*,*) "inoctal ",inOctal(grid%octreeRoot, currentPosition)
    if (.not.inOctal(grid%octreeRoot, currentPosition)) then
       distToNextCell = distanceToGridFromOutside(grid, currentPosition, direction, hitGrid) 
       if (hitGrid) then
          currentPosition = currentPosition + (distToNextCell+fudgeFac*grid%halfSmallestSubcell)*direction
          ntau = ntau + 1
          tauArray(ntau) = tauArray(nTau-1)
          xArray(nTau) = real(xArray(nTau-1) + distToNextCell)
       else
          goto 666
       endif
    endif




    if (PRESENT(startOctal)) then
       thisOctal => startOctal
       subcell = startSubcell
    else
       CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)
    endif

    if (PRESENT(distanceToEdge)) distanceToEdge = 0.d0
    do while (inOctal(grid%octreeRoot, currentPosition))

       call findSubcellLocal(currentPosition, thisOctal,subcell)
       if (.not.PRESENT(ross)) then
          call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, kappaSca=kappaSca, kappaAbs=kappaAbs)
          kappaExt = kappaAbs + kappaSca
       else
          call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, rosselandKappa=kappaExt)
          kappaExt = kappaExt * thisOctal%rho(subcell) * 1.d10
       endif
       sOctal => thisOctal
       call distanceToCellBoundary(grid, currentPosition, direction, DisttoNextCell, sOctal)
       
       currentPosition = currentPosition + (distToNextCell+fudgeFac*grid%halfSmallestSubcell)*direction
       

       tau = tau + distToNextCell*kappaExt
       if (PRESENT(nTau)) then
          nTau = nTau + 1
          xArray(nTau) = real(xArray(nTau-1) + distToNextCell)
          tauArray(nTau) = real(tau)
          if (nTau == SIZE(xArray)) then
             call writeWarning("tau array size reached. aborting tau run")
             goto 666
          endif
       endif
       if (PRESENT(distanceToEdge)) distanceToEdge = distanceToEdge + distToNextCell
       if (PRESENT(tauMax)) then
          if (tau > tauMax) exit
       endif
    end do

666 continue
    if (present(debug)) then
       write(*,*) "current pos ",currentposition,modulus(currentposition)
       write(*,*) "direction ", direction
       write(*,*) "source pos ",globalSourceArray(globalnSource)%position
    endif

    call distanceToSource(globalSourceArray, globalnSource, currentPosition, direction, hitSource, distance, sourcenumber)
       if (present(debug)) write(*,*) "hit source ",hitsource
    if (hitSource) then
       nTau = nTau + 1
       tauArray(nTau) = 1.e20
       xArray(nTau) = xArray(ntau-1) + real(distance)
       tau = 1.d20
    endif
  end subroutine tauAlongPathFast

  subroutine tauAlongPath2(ilambda, grid, rVec, direction, tau, tauMax, ross, startOctal, startSubcell)
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec, direction, currentPosition, beforeVec, afterVec
    integer :: iLambda
    real(double), intent(out) :: tau
    real(double) :: distToNextCell
    real(double), optional :: tauMax
    type(OCTAL), pointer :: thisOctal, sOctal
    type(OCTAL), pointer, optional :: startOctal
    integer, optional :: startSubcell
    real(double) :: fudgeFac = 1.d-1
    real(double) :: kappaSca, kappaAbs, kappaExt
    integer :: subcell
    logical, optional :: ross
    kappaAbs = 0.d0; kappaSca = 0.d0
    tau = 0.d0
    currentPosition = rVec
    

    if (PRESENT(startOctal)) then
       thisOctal => startOctal
       subcell = startSubcell
    else
       CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)
    endif

    do while (inOctal(grid%octreeRoot, currentPosition))

       call findSubcellLocal(currentPosition, thisOctal,subcell)
       if (.not.PRESENT(ross)) then
          call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, kappaSca=kappaSca, kappaAbs=kappaAbs)
          kappaExt = kappaAbs + kappaSca
       else
          call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, rosselandKappa=kappaExt)
          kappaExt = kappaExt * thisOctal%rho(subcell) * 1.d10
       endif
       sOctal => thisOctal
       call distanceToCellBoundary(grid, currentPosition, direction, DisttoNextCell, sOctal)
       
       beforeVec = VECTOR(currentPosition%x, currentPosition%y,0.d0)
       currentPosition = currentPosition + (distToNextCell+fudgeFac*grid%halfSmallestSubcell)*direction
       afterVec = VECTOR(currentPosition%x, currentPosition%y,0.d0)
       

       tau = tau + distToNextCell*kappaExt
       if (PRESENT(tauMax)) then
          if (tau > tauMax) exit
       endif
    end do
  end subroutine tauAlongPath2

  subroutine getMagStreamValues3(thisOctal, subcell, stream)
    type(OCTAL) :: thisOctal
    integer :: subcell
    type(STREAMTYPE) :: stream
    real(double) :: rho
    real  :: temperature
    type(VECTOR) :: vel
    integer :: n, i

    temperature = 0.
    rho = 0.d0
    vel = VECTOR(0.d0, 0.d0, 0.d0)
    
    n = 0
    do i = 1, stream%nSamples
       if (inSubcell(thisOctal, subcell, stream%position(i))) then
          n = n + 1
          rho = rho + stream%rho(i)
          temperature = temperature + stream%temperature(i)
          vel = vel + stream%velocity(i)
       endif
    enddo
    if (n ==0 ) then
       thisOctal%inFlow(subcell) = .false.
       thisOctal%rho(subcell) = 1.d-25
       thisOctal%temperature(subcell) = 6000.
       thisOctal%velocity(subcell) = VECTOR(0.,0.,0.)
    else
       thisOctal%inFlow(subcell) = .true.
       thisOctal%rho(subcell) = rho / dble(n)
       thisOctal%temperature(subcell) = temperature/real(n)
       thisOctal%velocity(subcell) = vel / dble(n)
       if (associated(thisOctal%microturb)) thisOctal%microturb(subcell) = 50.d5/cspeed!!!!!!!!!!!!!!!!!!!!

       if (subcell == thisOctal%maxchildren) then
          call fillVelocityCorners(thisOctal, magstreamvelocity)
!          write(*,*) " "
!          write(*,*) thisOctal%velocity(subcell)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(14)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(15)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(17)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(18)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(23)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(24)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(26)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(27)*real(cspeed/1.e5)
       endif
    endif

  end subroutine getMagStreamValues3

  type(VECTOR)  function magStreamVelocity(point) result(vel)
    type(VECTOR),intent(in) :: point
    real(double) :: minDist, dist, dist1, dist2, t
    integer :: i, j, isample, istream
    minDist = 1.d30

    do iStream = 1, globalnStream
       do iSample = 2, globalStream(iStream)%nSamples
          dist = distancePointLineSegment(globalStream(iStream)%position(iSample-1), &
               globalStream(iStream)%position(iSample), point)
          if (dist < minDist) then
             minDist = dist
             i = iStream
             j = iSample-1
          endif
       enddo
    enddo
    dist1 = modulus(point - globalStream(i)%position(j))
    dist2 = modulus(point - globalStream(i)%position(j+1))
    t = (minDist - dist1)/ (dist2 - dist1)
    vel = globalStream(i)%velocity(j) + t * (globalStream(i)%velocity(j+1) - globalStream(i)%velocity(j))
!    write(*,*) modulus(vel)*cspeed/1.d5
  end function magStreamVelocity


  subroutine findNearestSample(thisStream, position, iSample, t, outsideStream)
    type(STREAMTYPE) :: thisStream
    type(VECTOR) :: rVec, position
    real(double) :: minDist, dist
    integer :: iSample
    real(double) :: t
    integer :: i
    logical :: outsideStream

    outsideStream = .true.
    minDist = 1.d30
    iSample = 0

    rVec = thisStream%position(thisStream%nSamples) - position
    if ((rVec .dot. thisStream%direction(thisStream%nSamples)) < 0.d0) then
       outsideStream = .true.
       iSample = thisStream%nSamples
       t = 0.d0
       goto 666
    endif


    if (modulus(position) > modulus(thisStream%position(thisStream%nSamples))) then
       iSample = thisStream%nSamples
       t = 0.d0
       goto 666
    endif
    
    if (modulus(position) < modulus(thisStream%position(1))) then
       iSample = 1
       t = 0.d0
       goto 666
    endif

    do i = 1, thisStream%nSamples
       rVec = thisStream%position(i) - position
       dist = modulus(rVec)
       if (dist < thisStream%streamRadius(i)) outsideStream = .false.
       if ( (rVec .dot. thisStream%direction(i)) >= 0.d0 ) then
          if (dist < minDist) then
             minDist = dist
             iSample = i
          endif
       endif

    enddo
    if ((iSample == 0).and.(.not.outsideStream)) then
       write(*,*) "error in find nearest sample"
       write(*,*) "position",position
       write(*,*) modulus(position)/modulus(thisStream%position(1))
       iSample = 1
    endif
    if (outSideStream) iSample = 1
    rVec = thisStream%position(iSample) - position
    if (modulus(thisStream%position(iSample+1)-thisStream%position(iSample)) /= 0.d0) then
       t = (rVec.dot.thisStream%direction(iSample)) / &
            modulus(thisStream%position(iSample+1)-thisStream%position(iSample))
    else
       t = 0.d0
    endif
666 continue

  end  subroutine findNearestSample

    subroutine genericAccretionSurface(surface, lineFreq,coreContFlux,fAccretion,totalLum)

    USE surface_mod, only: createProbs, sumSurface, SURFACETYPE
    use inputs_mod, only : tHotSpot, mDotParameter1, tTauriRstar
    use magnetic_mod, only : accretingAreaMahdavi, velocityMahdavi, inflowMahdavi
    type(SURFACETYPE) :: surface
    type(VECTOR) :: rVec
    real(double) :: v, area, T, flux, power, totalArea, accretingArea, mdot, totalMdot
    integer :: i
    real(double) :: totalLum
    REAL(double), INTENT(IN) :: coreContFlux
    REAL, INTENT(IN) :: lineFreq
    REAL, INTENT(OUT) :: fAccretion ! erg s^-1 Hz^-1
    real(double) :: astar, thisR, thisMdot, thisRho

    if (Writeoutput) write(*,*) "calculating generic accretion surface ",surface%nElements
    astar = accretingAreaMahdavi()
    thismdot = mDotparameter1*mSol/(365.25d0*24.d0*3600.d0)

    if (writeoutput.and.(Thotspot > 0.)) write(*,*) "Setting hot spot temperature to: ",thotspot
    accretingArea = 0.d0
    totalArea = 0.d0
    totallum = 0.d0
    totalmdot = 0.d0

    do i = 1, surface%nElements
       rVec = (modulus(surface%Element(i)%position-surface%centre)*1.001d0) &
            *surface%element(i)%norm
       surface%element(i)%hot = .false.
       area = (surface%element(i)%area*1.d20)
       totalArea = totalArea + area

       if (inflowMahdavi(rVec*1.d10)) then


          thisR = modulus(rVec)*1.d10

          v = modulus(velocityMahdavi(rVec))*cSpeed
          
          thisRho = 0.d0
          if (v /= 0.d0) then
             thisRho =  thismdot /(aStar * v)  * (ttauriRstar/thisR)**3 
          endif
          
          mdot = thisRho * v * area



          totalMdot = totalMdot + mdot
          power = 0.5d0 * mdot * v**2
          if (area /= 0.d0) then
             flux = power / area
          else
             flux = 0.d0
          endif
          totalLum = totalLum + power
          
          
          T = max((flux/stefanBoltz)**0.25d0,3.d0)


          surface%element(i)%hot = .true.
          allocate(surface%element(i)%hotFlux(surface%nNuHotFlux))

          if (Thotspot > 0.) t = thotspot
          
          surface%element(i)%hotFlux(:) = &
               pi*blackbody(REAL(T), 1.e8*REAL(cSpeed/surface%nuArray(:)))
          surface%element(i)%temperature = real(T)
          accretingArea = accretingArea + area
       end if
    enddo


    if (writeoutput) then
       write(*,*) "Spot fraction is: ",100.d0*accretingArea/totalArea, "%"
       write(*,'(a,1pe12.3,a)') "Mass accretion rate is: ", &
            (totalMdot/mSol)*(365.25d0*24.d0*3600.d0), " solar masses/year"
       
       if (accretingArea > 0.d0) then
          t = (totalLum/(accretingArea*stefanBoltz))**0.25d0
          write(*,*) "Approx accretion temperature is ",t, " kelvin"
       endif
    endif

    CALL createProbs(surface,lineFreq,coreContFlux,fAccretion)
    CALL sumSurface(surface)

  end subroutine genericAccretionSurface

    subroutine hotSpotSurface(surface, lineFreq,coreContFlux,fAccretion,totalLum)

    USE surface_mod, only: createProbs, sumSurface, SURFACETYPE
    use magnetic_mod, only : accretingAreaMahdavi, velocityMahdavi, inflowMahdavi
    type(SURFACETYPE) :: surface
    type(VECTOR) :: rVec
    real(double) :: area, T,  totalArea, accretingArea,  totalMdot
    integer :: i
    real(double) :: totalLum
    REAL(double), INTENT(IN) :: coreContFlux
    REAL, INTENT(IN) :: lineFreq
    REAL, INTENT(OUT) :: fAccretion ! erg s^-1 Hz^-1
    real(double) :: thetaSpot, phiSpot, rSpot, ang

    if (Writeoutput) write(*,*) "calculating generic hotspot surface ",surface%nElements

    accretingArea = 0.d0
    totalArea = 0.d0
    totallum = 0.d0
    totalmdot = 0.d0

    thetaSpot = 65.d0 * degtorad
    phiSpot = 0.d0
    T = 10000.d0
    rSpot = 20.d0 * degToRad

    rVec = VECTOR(cos(phiSpot)*sin(thetaSpot),sin(phiSpot)*sin(thetaSpot),cos(thetaSpot))
    do i = 1, SIZE(surface%element)
       area = (surface%element(i)%area*1.d20)
       totalArea = totalArea + area

       ang = acos(rVec.dot.surface%element(i)%norm)
       if (ang < rSpot) then
          surface%element(i)%hot = .true.

          allocate(surface%element(i)%hotFlux(surface%nNuHotFlux))

          
          surface%element(i)%hotFlux(:) = &
               pi*blackbody(REAL(T), 1.e8*REAL(cSpeed/surface%nuArray(:)))
          surface%element(i)%temperature = real(T)
          accretingArea = accretingArea + area

       else
          surface%element(i)%hot = .false.
       endif
    enddo


    if (writeoutput) then
       write(*,*) "Spot fraction is: ",100.d0*accretingArea/totalArea, "%"
    endif

    CALL createProbs(surface,lineFreq,coreContFlux,fAccretion)
    CALL sumSurface(surface)

  end subroutine hotSpotSurface


  subroutine allocateOctalAttributes(grid, thisOctal)
    use inputs_mod, only : mie,  nDustType, molecular, TminGlobal, &
         photoionization, hydrodynamics, timeDependentRT, nAtom, &
         lineEmission, atomicPhysics, photoionPhysics, dustPhysics, molecularPhysics, cmf!, storeScattered
    use inputs_mod, only : grainFrac, pdrcalc, xraycalc, useionparam
    use gridtype_mod, only: statEqMaxLevels
    use h21cm_mod, only: h21cm
#ifdef PDR
    use inputs_mod, only :  hlevel
#endif
    type(OCTAL), pointer :: thisOctal
    type(GRIDTYPE) :: grid
!    integer, parameter :: nTheta = 10 , nphi = 10
    integer :: i
#ifdef PDR
    integer :: nrays, nside
#endif
    thisOctal%rho = amr_min_rho
    thisOctal%gasOpacity = .false.
    thisOctal%temperature = TMinGlobal


    

    if (atomicPhysics.or.molecularPhysics.or.h21cm) then
       call allocateAttribute(thisOctal%iAnalyticalVelocity,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%cornerVelocity, 27)
    endif


    if ( h21cm ) then 
       call allocateAttribute(thisOctal%etaLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%chiLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%NH2, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%molabundance, thisOctal%maxChildren)
       return
    end if

    if (mie.or.dustPhysics) then
       call allocateAttribute(thisOctal%iAnalyticalVelocity,thisOctal%maxChildren)

       call allocateAttribute(thisOctal%oldFrac, thisOctal%maxChildren)
       thisOctal%oldFrac = 1.e-30
       call allocateAttribute(thisOctal%fixedTemperature, thisOctal%maxChildren)
       thisOctal%fixedTemperature = .false.
       call allocateAttribute(thisOctal%dustType, thisOctal%maxChildren)
       thisOctal%dustType = 1
       ALLOCATE(thisOctal%dusttypefraction(thisOctal%maxchildren,  nDustType))
       thisOctal%dustTypeFraction = 0.d0
       ALLOCATE(thisOctal%origdusttypefraction(thisOctal%maxchildren,  nDustType))
       thisOctal%origdustTypeFraction = 0.d0

       do i = 1, thisOctal%maxChildren
          thisOctal%dustTypeFraction(i,1:nDustType) = grainFrac(1:nDustType)
          thisOctal%origdustTypeFraction(i,1:nDustType) = grainFrac(1:nDustType)
       enddo
       thisOctal%inflow = .true.

!       if (storescattered) allocate(thisOctal%scatteredIntensity(thisOctal%maxChildren, ntheta, nPhi))

       call allocateAttribute(thisOctal%meanIntensity, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%diffusionApprox, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%changed, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nDiffusion, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%eDens, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%diffusionCoeff, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%oldeDens, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nDirectPhotons, thisOctal%maxChildren)
       
       call allocateAttribute(thisOctal%corner,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%boundaryPartner,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%boundaryCondition,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%edgeCell,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%underSampled, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%oldTemperature, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%kappaRoss, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%distanceGrid, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nCrossings, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nTot, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%etaCont, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%etaLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%chiLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%biasCont3D, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%biasLine3D, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%probDistLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%probDistCont, thisOctal%maxChildren)


    endif



    if (atomicPhysics) then
       call allocateAttribute(thisOctal%iAnalyticalVelocity,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%microturb, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%atomAbundance, thisOctal%maxChildren, nAtom)
       call allocateAttribute(thisOctal%biasCont3D, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%biasLine3D, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%etaLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%etaCont, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%chiLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%changed, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nTot, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%ne, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%kappaAbs, thisOctal%maxChildren,1)
       call allocateAttribute(thisOctal%kappaSca, thisOctal%maxChildren,1)
       call allocateAttribute(thisOctal%fixedTemperature, thisOctal%maxChildren)
       thisOctal%fixedTemperature = .false.


       if (.not.cmf) then
          call allocateAttribute(thisOctal%changed, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%nTot, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%N, thisOctal%maxChildren,stateqMaxLevels)
       endif
    endif

    if (molecular) then
       call allocateAttribute(thisOctal%iAnalyticalVelocity,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%molAbundance, thisOctal%maxChildren)
       thisOctal%molAbundance(:) = 1.e-30
       call allocateAttribute(thisOctal%temperatureGas, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%temperatureDust, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nh2, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%microturb, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%molmicroturb, thisOctal%maxChildren)
    endif


    if (photoionization.or.photoionPhysics) then
       call allocateAttribute(thisOctal%oldFrac, thisOctal%maxChildren)
       thisOctal%oldFrac = 1.e-30
       call allocateAttribute(thisOctal%dustType, thisOctal%maxChildren)
       thisOctal%dustType = 1
       ALLOCATE(thisOctal%dusttypefraction(thisOctal%maxchildren,  nDustType))
       thisOctal%dustTypeFraction = 0.d0
       do i = 1, thisOctal%maxChildren
          thisOctal%dustTypeFraction(i,1:nDustType) = grainFrac(1:nDustType)
       enddo
       ALLOCATE(thisOctal%origdusttypefraction(thisOctal%maxchildren,  nDustType))
       thisOctal%origdustTypeFraction = 0.d0
       do i = 1, thisOctal%maxChildren
          thisOctal%origdustTypeFraction(i,1:nDustType) = grainFrac(1:nDustType)
       enddo

       call allocateAttribute(thisOctal%diffusionApprox, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%changed, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nDiffusion, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%eDens, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%diffusionCoeff, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%oldeDens, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nDirectPhotons, thisOctal%maxChildren)
       
       call allocateAttribute(thisOctal%underSampled, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%oldTemperature, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%kappaRoss, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%distanceGrid, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nCrossings, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nTot, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%chiLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%biasLine3D, thisOctal%maxChildren)

       call allocateAttribute(thisOctal%corner,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%boundaryPartner,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%boundaryCondition,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%edgeCell,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%etaCont, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nh, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%ne, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nhi, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nhei, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nhii, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%biasCont3D, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%etaLine, thisOctal%maxChildren)

       call allocateAttribute(thisOctal%HHeating, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%tDust, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%HeHeating, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%radiationMomentum,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%kappaTimesFlux, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%UVvector, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%UVvectorPlus, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%UVvectorminus, thisOctal%maxChildren)


       allocate(thisOctal%ionFrac(1:thisOctal%maxchildren, 1:grid%nIon))
       allocate(thisOctal%photoionCoeff(1:thisOctal%maxchildren, 1:grid%nIon))
       allocate(thisOctal%sourceContribution(1:thisOctal%maxchildren, 1:grid%nIon))
       allocate(thisOctal%diffuseContribution(1:thisOctal%maxchildren, 1:grid%nIon))
       allocate(thisOctal%normSourceContribution(1:thisOctal%maxchildren, 1:grid%nIon))

       call allocateAttribute(thisOctal%UVvector,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%UVvectorPlus,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%UVvectorMinus,thisOctal%maxChildren)

    endif
    call allocateAttribute(thisOctal%dust_T, thisOctal%maxChildren)
    if(xraycalc .or. pdrcalc .or. useionparam) then
       call allocateAttribute(thisOctal%columnRho, thisOctal%maxChildren)
    end if


#ifdef PDR
    if(pdrcalc) then
       nside=2**hlevel
       nrays = 12*nside**2
       if(thisOctal%oned) nrays = 1 
       call allocateAttribute(thisOctal%UV, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%converged, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%level_converged, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%biChop, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%expanded, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%lastChange, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%tPrev, thisOctal%maxChildren)

       call allocateAttribute(thisOctal%Tlast, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%Tmin, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%Tmax, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%Tlow, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%Thigh, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%Tminarray, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%Tmaxarray, thisOctal%maxChildren)
!       call allocateAttribute(thisOctal%radsurface, thisOctal%maxChildren)

       allocate (thisOctal%coolingRate(1:thisOctal%maxChildren,1:4))
!       call allocateAttribute(thisOctal%heatingRate, thisOctal%maxChildren)
       allocate(thisOctal%heatingRate(1:thisOctal%maxchildren, 1:12))
       allocate(thisOctal%cii_pop(1:thisOctal%maxchildren, 1:5))
       allocate(thisOctal%ci_pop(1:thisOctal%maxchildren, 1:5))
       allocate(thisOctal%oi_pop(1:thisOctal%maxchildren, 1:5))
       allocate(thisOctal%c12o_pop(1:thisOctal%maxchildren, 1:41))
       allocate(thisOctal%relch(1:thisOctal%maxchildren, 1:4, 1:41))

       allocate(thisOctal%radsurface(1:thisOctal%maxchildren, 1:nrays))
       allocate(thisOctal%AV(1:thisOctal%maxchildren, 1:nrays))
       allocate(thisOctal%abundance(1:thisOctal%maxchildren, 1:33))
       allocate(thisOctal%thisColRho(1:thisOctal%maxchildren, 1:nrays, 1:33))

       allocate(thisOctal%ciiLine(1:thisOctal%maxchildren, 1:5, 1:5 ))
       allocate(thisOctal%ciiTransition(1:thisOctal%maxchildren, 1:5, 1:5))!, 1:41))
       allocate(thisOctal%ciLine(1:thisOctal%maxchildren, 1:5, 1:5 ))
       allocate(thisOctal%ciTransition(1:thisOctal%maxchildren, 1:5, 1:5))!, 1:41))
       allocate(thisOctal%oiLine(1:thisOctal%maxchildren, 1:5, 1:5 ))
       allocate(thisOctal%oiTransition(1:thisOctal%maxchildren, 1:5, 1:5))!, 1:41))
       allocate(thisOctal%c12oLine(1:thisOctal%maxchildren, 1:41, 1:41 ))
       allocate(thisOctal%c12oTransition(1:thisOctal%maxchildren, 1:41, 1:41))!, 1:41))


    end if
#endif
    if (lineEmission) then
       call allocateAttribute(thisOctal%probDistLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%probDistCont, thisOctal%maxChildren)
    endif


    if (associated(thisOctal%ionFrac)) thisOctal%ionFrac = 1.e-30
    
    if (associated(thisOctal%photoIonCoeff)) then
       thisOctal%photoIonCoeff = 0.d0
       thisOctal%sourceContribution = 0.d0
       thisOctal%diffuseContribution = 0.d0
       thisOctal%normSourceContribution = 0.d0
    endif

    if (timeDependentRT) then
       call allocateAttribute(thisOctal%uDens, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%aDot, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%distancegridaDot, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%distanceGridPhotonFromGas, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%distanceGridPhotonFromSource, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%photonEnergyDensityFromGas, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%photonEnergyDensityFromSource, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%photonEnergyDensity, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%oldphotonEnergyDensity, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%probDistCont, thisOctal%maxChildren)
    endif
    if (hydrodynamics) then

       call allocateAttribute(thisOctal%q_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%q_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%q_i_minus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%q_i_minus_2,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%q_amr_i_minus_1,thisOctal%maxchildren,2)
       call allocateAttribute(thisOctal%q_amr_i_plus_1,thisOctal%maxchildren,2)
       call allocateAttribute(thisOctal%u_interface,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%u_amr_interface,thisOctal%maxchildren,2)
       call allocateAttribute(thisOctal%u_amr_interface_i_plus_1,thisOctal%maxchildren,2)
       call allocateAttribute(thisOctal%flux_amr_i,thisOctal%maxchildren,2)
       call allocateAttribute(thisOctal%flux_amr_i_plus_1,thisOctal%maxchildren,2)
       call allocateAttribute(thisOctal%flux_amr_i_minus_1,thisOctal%maxchildren,2)
       call allocateAttribute(thisOctal%phiLimit_amr,thisOctal%maxchildren, 2)

       

       call allocateAttribute(thisOctal%fviscosity,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%x_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%x_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%x_i_minus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%x_i_minus_2,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%u_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%u_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%u_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%flux_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%flux_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%flux_i_minus_1,thisOctal%maxchildren)


       call allocateAttribute(thisOctal%phiLimit,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%qViscosity,thisOctal%maxchildren,3,3)

       call allocateAttribute(thisOctal%ghostCell,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%corner,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%feederCell,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%edgeCell,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%refinedLastTime,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%pressure_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%pressure_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%pressure_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%rho_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%divV,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%rhou,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%rhov,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%rhow,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%rhoe,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%rhoeLastTime,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%energy,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%phi_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%phi_gas,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%phi_stars,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%phi_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%phi_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%rho_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%rho_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%rhorv_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%rhorv_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%boundaryCondition,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%boundaryCell,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%boundaryPartner,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%gravboundaryPartner,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%changed,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%rLimit,thisOctal%maxChildren)

       call allocateAttribute(thisOctal%iEquationOfState,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%gamma,thisOctal%maxChildren)


       call allocateAttribute(thisOctal%radiationMomentum,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%kappaTimesFlux,thisOctal%maxChildren)


    endif
  end  subroutine allocateOctalAttributes


  subroutine deallocateOctalDynamicAttributes(thisOctal)
    type(OCTAL) :: thisOctal

    call deallocateAttribute(thisOctal%iEquationOfState)
    call deallocateAttribute(thisOctal%iAnalyticalVelocity)
    call deallocateAttribute(thisOctal%gamma)
    call deallocateAttribute(thisOctal%neighbourOctal)
    call deallocateAttribute(thisOctal%neighbourSubcell)
    call deallocateAttribute(thisOctal%cornerVelocity)
    call deallocateAttribute(thisOctal%cornerrho)
    call deallocateAttribute(thisOctal%qViscosity)
    call deallocateAttribute(thisOctal%diffusionApprox)
    call deallocateAttribute(thisOctal%fixedTemperature)
    call deallocateAttribute(thisOctal%nDiffusion)
    call deallocateAttribute(thisOctal%eDens)
    call deallocateAttribute(thisOctal%diffusionCoeff)
    call deallocateAttribute(thisOctal%oldeDens)
    call deallocateAttribute(thisOctal%nDirectPhotons)
    call deallocateAttribute(thisOctal%undersampled)
    call deallocateAttribute(thisOctal%oldTemperature)
    call deallocateAttribute(thisOctal%kappaRoss)
    call deallocateAttribute(thisOctal%distanceGrid)
    call deallocateAttribute(thisOctal%scatteredIntensity)
    call deallocateAttribute(thisOctal%meanIntensity)
    call deallocateAttribute(thisOctal%nCrossings)
    call deallocateAttribute(thisOctal%nTot)
    call deallocateAttribute(thisOctal%oldFrac)
    call deallocateAttribute(thisOctal%dusttype)
    call deallocateAttribute(thisOctal%dustTypeFraction)
    call deallocateAttribute(thisOctal%origdustTypeFraction)
    call deallocateAttribute(thisOctal%departCoeff)
    call deallocateAttribute(thisOctal%kappaAbs)
    call deallocateAttribute(thisOctal%kappaSca)
    call deallocateAttribute(thisOctal%chiLine)
    call deallocateAttribute(thisOctal%etaLine)
    call deallocateAttribute(thisOctal%etaCont)
    call deallocateAttribute(thisOctal%biasLine3D)
    call deallocateAttribute(thisOctal%biasCont3D)
    call deallocateAttribute(thisOctal%probDistLine)
    call deallocateAttribute(thisOctal%probDistCont)
    call deallocateAttribute(thisOctal%N)
    call deallocateAttribute(thisOctal%Ne)
    call deallocateAttribute(thisOctal%NH)
    call deallocateAttribute(thisOctal%molecularLevel)
    call deallocateAttribute(thisOctal%molcellparam)
    call deallocateAttribute(thisOctal%newmolecularLevel)
    call deallocateAttribute(thisOctal%oldmolecularLevel)
    call deallocateAttribute(thisOctal%oldestmolecularLevel)
    call deallocateAttribute(thisOctal%Adot)
    call deallocateAttribute(thisOctal%uDens)
    call deallocateAttribute(thisOctal%distanceGridAdot)
    call deallocateAttribute(thisOctal%distanceGridPhotonFromSource)
    call deallocateAttribute(thisOctal%distanceGridPhotonFromGas)
    call deallocateAttribute(thisOctal%photonEnergyDensityFromSource)
    call deallocateAttribute(thisOctal%photonEnergyDensityFromGas)
    call deallocateAttribute(thisOctal%photonEnergyDensity)
    call deallocateAttribute(thisOctal%oldphotonEnergyDensity)
    call deallocateAttribute(thisOctal%temperaturedust)
    call deallocateAttribute(thisOctal%temperaturegas)
    call deallocateAttribute(thisOctal%NH2)
    call deallocateAttribute(thisOctal%microturb)
    call deallocateAttribute(thisOctal%molmicroturb)
    call deallocateAttribute(thisOctal%atomLevel)
    call deallocateAttribute(thisOctal%atomAbundance)
    call deallocateAttribute(thisOctal%newatomLevel)
    call deallocateAttribute(thisOctal%jnu)
    call deallocateAttribute(thisOctal%jnuCont)
    call deallocateAttribute(thisOctal%jnuLine)
    call deallocateAttribute(thisOctal%tau)
    call deallocateAttribute(thisOctal%bnu)
    call deallocateAttribute(thisOctal%molAbundance)
    call deallocateAttribute(thisOctal%convergence)
    call deallocateAttribute(thisOctal%levelconvergence)
    call deallocateAttribute(thisOctal%nsplit)
    call deallocateAttribute(thisOctal%NHI)
    call deallocateAttribute(thisOctal%NHII)
    call deallocateAttribute(thisOctal%NHeI)
    call deallocateAttribute(thisOctal%NHeII)
    call deallocateAttribute(thisOctal%Hheating)
    call deallocateAttribute(thisOctal%tDust)
    call deallocateAttribute(thisOctal%Heheating)
    call deallocateAttribute(thisOctal%ionFrac)
#ifdef PDR
    call deallocateAttribute(thisOctal%ciiLine)
    call deallocateAttribute(thisOctal%relch)
    call deallocateAttribute(thisOctal%ciiTransition)
    call deallocateAttribute(thisOctal%ciLine)
    call deallocateAttribute(thisOctal%ciTransition)
    call deallocateAttribute(thisOctal%oiLine)
    call deallocateAttribute(thisOctal%oiTransition)
    call deallocateAttribute(thisOctal%c12oLine)
    call deallocateAttribute(thisOctal%c12oTransition)

    call deallocateAttribute(thisOctal%coolingRate)
    call deallocateAttribute(thisOctal%heatingRate)
    call deallocateAttribute(thisOctal%thisColRho)
    call deallocateAttribute(thisOctal%abundance)
    call deallocateAttribute(thisOctal%UV)
    call deallocateAttribute(thisOctal%converged)
    call deallocateAttribute(thisOctal%level_converged)
    call deallocateAttribute(thisOctal%biChop)
    call deallocateAttribute(thisOctal%expanded)
    call deallocateAttribute(thisOctal%lastChange)
    call deallocateAttribute(thisOctal%tPreV)

    call deallocateAttribute(thisOctal%radsurface)
    call deallocateAttribute(thisOctal%Tlast)
    call deallocateAttribute(thisOctal%Tmin)
    call deallocateAttribute(thisOctal%Tlow)
    call deallocateAttribute(thisOctal%Thigh)
    call deallocateAttribute(thisOctal%Tmax)
    call deallocateAttribute(thisOctal%Tminarray)
    call deallocateAttribute(thisOctal%Tmaxarray)
    call deallocateAttribute(thisOctal%cii_pop)
    call deallocateAttribute(thisOctal%ci_pop)
    call deallocateAttribute(thisOctal%oi_pop)
    call deallocateAttribute(thisOctal%c12o_pop)
#endif
    call deallocateAttribute(thisOctal%dust_T)
    call deallocateAttribute(thisOctal%photoIonCoeff)
    call deallocateAttribute(thisOctal%sourceContribution)
    call deallocateAttribute(thisOctal%diffuseContribution)
    call deallocateAttribute(thisOctal%normSourceContribution)
    call deallocateAttribute(thisOctal%gas_particle_list)
    call deallocateAttribute(thisOctal%changed)
    call deallocateAttribute(thisOctal%mpiBoundaryStorage)
    call deallocateAttribute(thisOctal%mpiCornerStorage)
    call deallocateAttribute(thisOctal%q_i)
    call deallocateAttribute(thisOctal%fViscosity)
    call deallocateAttribute(thisOctal%q_i_plus_1)
    call deallocateAttribute(thisOctal%q_i_minus_1)
    call deallocateAttribute(thisOctal%q_i_minus_2)
    call deallocateAttribute(thisOctal%x_i)
    call deallocateAttribute(thisOctal%x_i_plus_1)
    call deallocateAttribute(thisOctal%x_i_minus_1)
    call deallocateAttribute(thisOctal%x_i_minus_2)
    call deallocateAttribute(thisOctal%u_interface)
    call deallocateAttribute(thisOctal%u_amr_interface)
    call deallocateAttribute(thisOctal%u_amr_interface_i_plus_1)
    call deallocateAttribute(thisOctal%u_i)
    call deallocateAttribute(thisOctal%u_i_plus_1)
    call deallocateAttribute(thisOctal%u_i_minus_1)
    call deallocateAttribute(thisOctal%flux_i)
    call deallocateAttribute(thisOctal%flux_i_plus_1)
    call deallocateAttribute(thisOctal%flux_i_minus_1)
    call deallocateAttribute(thisOctal%phiLimit)

    call deallocateAttribute(thisOctal%q_amr_i_minus_1)
    call deallocateAttribute(thisOctal%q_amr_i_plus_1)
    call deallocateAttribute(thisOctal%u_interface)
    call deallocateAttribute(thisOctal%u_amr_interface)
    call deallocateAttribute(thisOctal%flux_amr_i)
    call deallocateAttribute(thisOctal%flux_amr_i_plus_1)
    call deallocateAttribute(thisOctal%flux_amr_i_minus_1)
    call deallocateAttribute(thisOctal%phiLimit_amr)


    call deallocateAttribute(thisOctal%rLimit)
    call deallocateAttribute(thisOctal%ghostCell)
    call deallocateAttribute(thisOctal%columnRho)
    call deallocateAttribute(thisOctal%feederCell)
    call deallocateAttribute(thisOctal%corner)
    call deallocateAttribute(thisOctal%edgeCell)
    call deallocateAttribute(thisOctal%refinedLastTime)
    call deallocateAttribute(thisOctal%divV)
    call deallocateAttribute(thisOctal%rhou)
    call deallocateAttribute(thisOctal%rhov)
    call deallocateAttribute(thisOctal%rhow)
    call deallocateAttribute(thisOctal%rhoE)
    call deallocateAttribute(thisOctal%rhoeLastTime)
    call deallocateAttribute(thisOctal%energy)
    call deallocateAttribute(thisOctal%pressure_i)
    call deallocateAttribute(thisOctal%pressure_i_plus_1)
    call deallocateAttribute(thisOctal%pressure_i_minus_1)
    call deallocateAttribute(thisOctal%rho_i_minus_1)
    call deallocateAttribute(thisOctal%tempStorage)
    call deallocateAttribute(thisOctal%boundaryPartner)
    call deallocateAttribute(thisOctal%gravboundaryPartner)
    call deallocateAttribute(thisOctal%radiationMomentum)
    call deallocateAttribute(thisOctal%kappaTimesFlux)

    call deallocateAttribute(thisOctal%UVvector)
    call deallocateAttribute(thisOctal%UVvectorPlus)
    call deallocateAttribute(thisOctal%UVvectorMinus)
    call deallocateAttribute(thisOctal%phi_i)
    call deallocateAttribute(thisOctal%phi_i_plus_1)
    call deallocateAttribute(thisOctal%phi_i_minus_1)
    call deallocateAttribute(thisOctal%phi_stars)
    call deallocateAttribute(thisOctal%phi_gas)
    call deallocateAttribute(thisOctal%rho_i_minus_1)
    call deallocateAttribute(thisOctal%rho_i_plus_1)
    call deallocateAttribute(thisOctal%rhorv_i_minus_1)
    call deallocateAttribute(thisOctal%rhorv_i_plus_1)
    call deallocateAttribute(thisOctal%boundaryCondition)
    call deallocateAttribute(thisOctal%boundaryCell)

  end subroutine deallocateOctalDynamicAttributes


  recursive subroutine setupNeighbourPointers(grid, thisOctal)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer :: neighbourOctal
    integer :: neighbourSubcell
    type(octal), pointer  :: child 
    type(VECTOR) :: centre, rvec
    real(double) :: d
    integer :: subcell, i
    type(VECTOR), parameter :: zAxis = VECTOR(0.d0, 0.d0, 1.d0)
    type(VECTOR), parameter :: yAxis = VECTOR(0.d0, 1.d0, 0.d0)
    type(VECTOR), parameter :: xAxis = VECTOR(1.d0, 0.d0, 0.d0)
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call  setupNeighbourPointers(grid, child)
                exit
             end if
          end do
       else

          if (.not.associated(thisOctal%neighbourOctal)) then
             allocate(thisOctal%neighbourOctal(thisOctal%maxChildren, 6, 4))
          endif
          if (.not.associated(thisOctal%neighbourSubcell)) then
             allocate(thisOctal%neighbourSubcell(thisOctal%maxChildren, 6, 4))
             thisOctal%neighbourSubcell = -666
          endif

          centre = subcellCentre(thisOctal, subcell)
          d = thisOctal%subcellSize/2.d0

          ! top face (plus z)

          rVec = centre + zAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) - yAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 1, 1)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 1, 1) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 1, 1)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 1, 1) = -1
          endif

          if (inOctal(grid%octreeRoot, rVec)) then
             rVec = centre + zAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) - yAxis*(d/2.d0)
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 1, 2)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 1, 2) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 1, 2)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 1, 2) = -1
          endif


          rVec = centre + zAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) + yAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 1, 3)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 1, 3) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 1, 3)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 1, 3) = -1
          endif

          rVec = centre + zAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) + yAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 1, 4)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 1, 4) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 1, 4)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 1, 4) = -1
          endif

          ! bottom face (minus z)

          rVec = centre - zAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) - yAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 2, 1)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 2, 1) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 2, 1)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 2, 1) = -1
          endif

          rVec = centre - zAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) - yAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 2, 2)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 2, 2) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 2, 2)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 2, 2) = -1
          endif

          rVec = centre - zAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) + yAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 2, 3)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 2, 3) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 2, 3)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 2, 3) = -1
          endif

          rVec = centre - zAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) + yAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 2, 4)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 2, 4) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 2, 4)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 2, 4) = -1
          endif


          ! left face (minus x)

          rVec = centre - xAxis * (d+grid%halfSmallestSubcell*00.1d0) - yAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 3, 1)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 3, 1) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 3, 1)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 3, 1) = -1
          endif

          rVec = centre - xAxis * (d+grid%halfSmallestSubcell*00.1d0) + yAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 3, 2)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 3, 2) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 3, 2)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 3, 2) = -1
          endif

          rVec = centre - xAxis * (d+grid%halfSmallestSubcell*00.1d0) - yAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 3, 3)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 3, 3) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 3, 3)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 3, 3) = -1
          endif

          rVec = centre - xAxis * (d+grid%halfSmallestSubcell*00.1d0) + yAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 3, 4)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 3, 4) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 3, 4)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 3, 4) = -1
          endif

          ! right face (plus x)

          rVec = centre + xAxis * (d+grid%halfSmallestSubcell*00.1d0) - yAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 4, 1)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 4, 1) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 4, 1)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 4, 1) = -1
          endif

          rVec = centre + xAxis * (d+grid%halfSmallestSubcell*00.1d0) + yAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 4, 2)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 4, 2) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 4, 2)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 4, 2) = -1
          endif

          rVec = centre + xAxis * (d+grid%halfSmallestSubcell*00.1d0) - yAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 4, 3)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 4, 3) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 4, 3)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 4, 3) = -1
          endif

          rVec = centre + xAxis * (d+grid%halfSmallestSubcell*00.1d0) + yAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 4, 4)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 4, 4) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 4, 4)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 4, 4) = -1
          endif


          ! front face (plus y)

          rVec = centre + yAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 5, 1)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 5, 1) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 5, 1)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 5, 1) = -1
          endif

          rVec = centre + yAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 5, 2)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 5, 2) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 5, 2)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 5, 2) = -1
          endif

          rVec = centre + yAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 5, 3)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 5, 3) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 5, 3)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 5, 3) = -1
          endif

          rVec = centre + yAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 5, 4)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 5, 4) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 5, 4)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 5, 4) = -1
          endif

          ! back face (minus y)

          rVec = centre - yAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 6, 1)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 6, 1) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 6, 1)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 6, 1) = -1
          endif

          rVec = centre - yAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 6, 2)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 6, 2) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 6, 2)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 6, 2) = -1
          endif

          rVec = centre - yAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 6, 3)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 6, 3) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 6, 3)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 6, 3) = -1
          endif

          rVec = centre - yAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 6, 4)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 6, 4) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 6, 4)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 6, 4) = -1
          endif
!          write(*,*) "subcell ",subcell
!          write(*,*) "top ",thisOctal%neighbourSubcell(subcell, 1, 1:4)
!          write(*,*) "bot ",thisOctal%neighbourSubcell(subcell, 2, 1:4)
!          write(*,*) "left ",thisOctal%neighbourSubcell(subcell, 3, 1:4)
!          write(*,*) "right ",thisOctal%neighbourSubcell(subcell, 4, 1:4)
!          write(*,*) "front ",thisOctal%neighbourSubcell(subcell, 5, 1:4)
!          write(*,*) "back ",thisOctal%neighbourSubcell(subcell, 6, 1:4)
       endif
    enddo
  end subroutine setupNeighbourPointers

  subroutine howmanysplits()
    implicit none
    character(len = 100) :: message    

    if(mass_split .ne. 0) then
       write(message, *) "There were ", mass_split, " splits by mass"
       call writeinfo(message, TRIVIAL)
    endif

    if(mass_split2 .ne. 0) then
       write(message, *) "There were ", mass_split2, " splits by mass2"
       call writeinfo(message, TRIVIAL)
    endif
    
    if(density_split .ne. 0) then
       write(message, *) "There were ", density_split, " splits by density"
       call writeinfo(message, TRIVIAL)
    endif

    if(both_split .ne. 0) then
       write(message, *) "There were ", both_split, " splits by both"
       call writeinfo(message, TRIVIAL)
    endif

    if(maxdensity_split .ne. 0) then
       write(message, *) "There were ", maxdensity_split, " splits by maximum density"
       call writeinfo(message, TRIVIAL)
    endif

    if(velocity_split .ne. 0) then
       write(message, *) "There were ", velocity_split, " splits by velocity"
       call writeinfo(message, TRIVIAL)
    endif

    if(scaleheighta_count .ne. 0) then
       write(message, *) "There were ", scaleheighta_count, " splits by scaleheight condition A"
       call writeinfo(message, TRIVIAL)
    endif

    if(scaleheightb_count .ne. 0) then
       write(message, *) "There were ", scaleheightb_count, " splits by scaleheight condition B"
       call writeinfo(message, TRIVIAL)
    endif

    if(scaleheightc_count .ne. 0) then
       write(message, *) "There were ", scaleheightc_count, " splits by scaleheight condition C"
       call writeinfo(message, TRIVIAL)
    endif

  end subroutine howmanysplits

  function returnScatteredIntensity(position, thisOctal, subcell, uHat) result(intensity)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    type(VECTOR) :: uHat, thisVec, position
    real(double) :: intensity, thisTheta, thisPhi, ang
    integer :: nTheta, nPhi
    integer :: iTheta, iPhi
    
    nTheta = SIZE(thisOctal%scatteredIntensity,2)
    nPhi = SIZE(thisOctal%scatteredIntensity,3)

    ang = atan2(position%y, position%x)

    thisVec = uHat
    
    if (thisOctal%twoD) then
       thisVec = rotateZ(uHat, -ang)
    endif

    thisTheta = acos(thisvec%z)
    thisPhi = atan2(thisVec%y,thisVec%x)

    if (thisPhi < 0.d0) thisPhi = thisPhi + twoPi 
    iTheta = nint((thisTheta / pi) * dble(nTheta-1))+1
    iphi = nint((thisPhi / twoPi) * dble(nPhi-1))+1

    intensity = thisOctal%scatteredIntensity(subcell,iTheta, iPhi) 
!    intensity = SUM(thisOctal%scatteredIntensity(subcell,:,:))/100.d0
  end function returnScatteredIntensity

  subroutine nnint(t1,t2,t3,weights)

    real(double) :: t1, t2, t3
!    real(double), save :: t1old, t2old, t3old
    real(double) :: weights(27)
!    real(double), save :: oldweights(27)
      
!    if(t1 .eq. t1old) then
!       if(t2 .eq. t2old) then
!          if(t3 .eq. t3old) then
!             weights = oldweights
!             return
!          endif
!       endif
!    endif

    weights(:) = 0.d0

    if(t3 .le. 0.25) then
       if(t2 .le. 0.25) then
          if(t1 .le. 0.25) then
             weights(1) = 1.d0
          elseif(t1 .ge. 0.75) then
             weights(3) = 1.d0
          else
             weights(2) = 1.d0
          endif
       elseif(t2 .ge. 0.75) then
          if(t1 .le. 0.25) then
             weights(7) = 1.d0
          elseif(t1 .ge. 0.75) then
             weights(9) = 1.d0
          else
             weights(8) = 1.d0
          endif
       else
          if(t1 .le. 0.25) then
             weights(4) = 1.d0
          elseif(t1 .ge. 0.75) then
             weights(6) = 1.d0
          else
             weights(5) = 1.d0
          endif
       endif
    elseif(t3 .ge. 0.75) then
       if(t2 .le. 0.25) then
          if(t1 .le. 0.25) then
             weights(19) = 1.d0
          elseif(t1 .ge. 0.75) then
             weights(21) = 1.d0
          else
             weights(20) = 1.d0
          endif
       elseif(t2 .ge. 0.75) then
          if(t1 .le. 0.25) then
             weights(25) = 1.d0
          elseif(t1 .ge. 0.75) then
             weights(27) = 1.d0
          else
             weights(26) = 1.d0
          endif
       else
          if(t1 .le. 0.25) then
             weights(22) = 1.d0
          elseif(t1 .ge. 0.75) then
             weights(24) = 1.d0
          else
             weights(23) = 1.d0
          endif
       endif
    else
       if(t2 .le. 0.25) then
          if(t1 .le. 0.25) then
             weights(10) = 1.d0
          elseif(t1 .ge. 0.75) then
             weights(12) = 1.d0
          else
             weights(11) = 1.d0
          endif
       elseif(t2 .ge. 0.75) then
          if(t1 .le. 0.25) then
             weights(16) = 1.d0
          elseif(t1 .ge. 0.75) then
             weights(18) = 1.d0
          else
             weights(17) = 1.d0
          endif
       else
          if(t1 .le. 0.25) then
             weights(13) = 1.d0
          elseif(t1 .ge. 0.75) then
             weights(15) = 1.d0
          else
             weights(14) = 1.d0
          endif
       endif
    endif

  end subroutine nnint

  real(double) function amrgriddensity(position, grid, linearinterp) result(rho)

    type(GRIDTYPE) :: grid
    type(OCTAL),save, pointer :: resultoctal
    integer :: subcell
    type(VECTOR) :: position, centre, point_local
    real(double) :: t1, t2, t3
    real(double) :: inc, fac
    real(double) :: weights(27)
!    real(double) :: rhoa(1)

    logical, save :: firsttime = .true.
    logical, optional :: linearinterp
    logical :: linear
    integer :: i
!$OMP THREADPRIVATE (firsttime, resultoctal)

    weights = 0.d0
    rho = 0.d0

    if(present(linearinterp)) then
       linear = linearinterp
    else
       linear = .true.
    endif

    if(firsttime) then
       CALL findSubcellTD(position,grid%octreeroot, resultOctal,subcell)
       firsttime = .false.
    endif

 !   if(.not. inoctal(resultoctal,position)) then
!       CALL findSubcellTD(position,grid%octreeroot, resultOctal,subcell)
       CALL findSubcellLocal(position,resultOctal,subcell)
!    endif

    if (resultoctal%threeD) then
       point_local = position
    elseif (resultoctal%twoD) then
       point_local = projectToXZ(position)
    else !oneD
       point_local = VECTOR(modulus(position), 0.d0, 0.d0)
    end if

    inc = resultOctal%subcellSize
    centre = resultoctal%centre
    fac = 0.5d0 / resultOctal%subcellsize

    t1 = fac * (point_local%x - (centre%x - inc))
    t2 = fac * (point_local%y - (centre%y - inc))
    t3 = fac * (point_local%z - (centre%z - inc))

    if(resultoctal%oneD) then
       centre = VECTOR(modulus(resultoctal%centre), 0.d0, 0.d0)
       t1 = fac * (point_local%x - (centre%x - inc))

       rho = (2.d0 * t1**2 - 3.d0 * t1 + 1.d0) * resultOctal%cornerrho(1) + &
            (-4.d0 * t1 * (t1 - 1.d0))         * resultOctal%cornerrho(2) + &
            (t1 * (2.d0 * t1 - 1.d0))          * resultOctal%cornerrho(3)

       return
    endif

    if(linear) then

       inc = 0.5 * resultOctal%subcellSize
       centre = subcellCentre(resultOctal,subcell)
       fac = 1. / resultOctal%subcellsize
      
       t1 = MAX(0.0_oc, fac * (point_local%x - (centre%x - inc)))
       t2 = MAX(0.0_oc, fac * (point_local%y - (centre%y - inc)))
       t3 = MAX(0.0_oc, fac * (point_local%z - (centre%z - inc)))

       SELECT CASE(subcell)
               
       CASE(1)
          rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 1) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 2) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 4) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 5) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(10) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(11) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(13) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(14)
               
            CASE(2)
               rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 2) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 3) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 5) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 6) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(11) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(12) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(14) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(15)
               
            CASE(3)
               rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 4) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 5) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 7) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 8) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(13) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(14) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(16) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(17)
               
            CASE(4)
               rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 5) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 6) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 8) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 9) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(14) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(15) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(17) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(18)
               
            CASE(5)
               rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(10) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(11) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(13) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(14) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(19) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(20) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(22) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(23)
               
            CASE(6)
               rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(11) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(12) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(14) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(15) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(20) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(21) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(23) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(24)
            
            CASE(7)
               rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(13) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(14) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(16) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(17) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(22) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(23) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(25) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(26)
               
            CASE(8)
               rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(14) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(15) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(17) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(18) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(23) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(24) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(26) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(27)
            CASE DEFAULT
               PRINT *, 'Invalid subcell in amrGridDensity'
               
            end SELECT
            
         else
!            call regular_tri_quadint(t1,t2,t3,weights)
!            if(rho .lt. 0.d0) then
!               rhoa = resultoctal%cornerrho(maxloc(weights))
!               rho = rhoa(1)
!            endif

            call nnint(t1,t2,t3,weights)
!            rho = sum(weights(:) * resultoctal%cornerrho(:))
            do i = 1, 27
               rho = rho + weights(i) * resultoctal%cornerrho(i)
            enddo
         endif
    
       end function amrgriddensity

  SUBROUTINE fillDensityCorners(thisOctal,densityFunc, velocityfunc)
    ! store the density values at the subcell corners of an octal so
    !   that they can be used for interpolation.

    IMPLICIT NONE
  
    TYPE(octal), INTENT(INOUT) :: thisOctal
    real(oct)      :: r1, r2, r3
    real(oct)      :: phi1, phi2, phi3
    real(oct)      :: x1, x2, x3
    real(oct)      :: y1, y2, y3
    real(oct)      :: z1, z2, z3
    
    INTERFACE 
      real(double) FUNCTION densityFunc(point)
        USE vector_mod
        USE gridtype_mod
        TYPE(vector), INTENT(IN) :: point
      END FUNCTION densityFunc
    END INTERFACE

    INTERFACE
      type(VECTOR) FUNCTION velocityFunc(point)
        USE vector_mod
        USE gridtype_mod
        TYPE(vector), INTENT(IN) :: point
      END FUNCTION velocityFunc
    END INTERFACE

    if (thisOctal%oneD) then
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       y1 = 0.d0
       z1 = 0.d0
       thisOctal%cornerrho(1) = densityFunc(vector(x1,y1,z1))
       thisOctal%cornerrho(2) = densityFunc(vector(x2,y1,z1))
       thisOctal%cornerrho(3) = densityFunc(vector(x3,y1,z1))
       goto 667
    endif

    if (thisOctal%threed) then
       if (.not.thisOctal%cylindrical) then ! 3d cartesian case
          ! we first store the values we use to assemble the position vectors
          
          x1 = thisOctal%centre%x - thisOctal%subcellSize
          x2 = thisOctal%centre%x
          x3 = thisOctal%centre%x + thisOctal%subcellSize
          
          y1 = thisOctal%centre%y - thisOctal%subcellSize
          y2 = thisOctal%centre%y
          y3 = thisOctal%centre%y + thisOctal%subcellSize
          
          z1 = thisOctal%centre%z - thisOctal%subcellSize
          z2 = thisOctal%centre%z
          z3 = thisOctal%centre%z + thisOctal%subcellSize

                    
          ! now store the 'base level' values
                              
          thisOctal%cornerrho(1)      =  densityFunc(vector(x1,y1,z1))
          thisOctal%cornerVelocity(1) = velocityFunc(vector(x1,y1,z1))
          thisOctal%cornerrho(2)      =  densityFunc(vector(x2,y1,z1))
          thisOctal%cornerVelocity(2) = velocityFunc(vector(x2,y1,z1))
          thisOctal%cornerrho(3)      =  densityFunc(vector(x3,y1,z1))
          thisOctal%cornerVelocity(3) = velocityFunc(vector(x3,y1,z1))
          thisOctal%cornerrho(4)      =  densityFunc(vector(x1,y2,z1))
          thisOctal%cornerVelocity(4) = velocityFunc(vector(x1,y2,z1))
          thisOctal%cornerrho(5)      =  densityFunc(vector(x2,y2,z1))
          thisOctal%cornerVelocity(5) = velocityFunc(vector(x2,y2,z1))
          thisOctal%cornerrho(6)      =  densityFunc(vector(x3,y2,z1))
          thisOctal%cornerVelocity(6) = velocityFunc(vector(x3,y2,z1))
          thisOctal%cornerrho(7)      =  densityFunc(vector(x1,y3,z1))
          thisOctal%cornerVelocity(7) = velocityFunc(vector(x1,y3,z1))
          thisOctal%cornerrho(8)      =  densityFunc(vector(x2,y3,z1))
          thisOctal%cornerVelocity(8) = velocityFunc(vector(x2,y3,z1))
          thisOctal%cornerrho(9)      =  densityFunc(vector(x3,y3,z1))
          thisOctal%cornerVelocity(9) = velocityFunc(vector(x3,y3,z1))

          ! middle level
          
          thisOctal%cornerrho(10)      =  densityFunc(vector(x1,y1,z2))
          thisOctal%cornerVelocity(10) = velocityFunc(vector(x1,y1,z2))
          thisOctal%cornerrho(11)      =  densityFunc(vector(x2,y1,z2))
          thisOctal%cornerVelocity(11) = velocityFunc(vector(x2,y1,z2))
          thisOctal%cornerrho(12)      =  densityFunc(vector(x3,y1,z2))
          thisOctal%cornerVelocity(12) = velocityFunc(vector(x3,y1,z2))
          thisOctal%cornerrho(13)      =  densityFunc(vector(x1,y2,z2))
          thisOctal%cornerVelocity(13) = velocityFunc(vector(x1,y2,z2))
          thisOctal%cornerrho(14)      =  densityFunc(vector(x2,y2,z2))
          thisOctal%cornerVelocity(14) = velocityFunc(vector(x2,y2,z2))

          thisOctal%cornerrho(15)      =  densityFunc(vector(x3,y2,z2))
          thisOctal%cornerVelocity(15) = velocityFunc(vector(x3,y2,z2))
          thisOctal%cornerrho(16)      =  densityFunc(vector(x1,y3,z2))
          thisOctal%cornerVelocity(16) = velocityFunc(vector(x1,y3,z2))
          thisOctal%cornerrho(17)      =  densityFunc(vector(x2,y3,z2))
          thisOctal%cornerVelocity(17) = velocityFunc(vector(x2,y3,z2))
          thisOctal%cornerrho(18)      =  densityFunc(vector(x3,y3,z2))
          thisOctal%cornerVelocity(18) = velocityFunc(vector(x3,y3,z2))
          ! top level         

          thisOctal%cornerrho(19)      =  densityFunc(vector(x1,y1,z3))
          thisOctal%cornerVelocity(19) = velocityFunc(vector(x1,y1,z3))
          thisOctal%cornerrho(20)      =  densityFunc(vector(x2,y1,z3))
          thisOctal%cornerVelocity(20) = velocityFunc(vector(x2,y1,z3))
          thisOctal%cornerrho(21)      =  densityFunc(vector(x3,y1,z3))
          thisOctal%cornerVelocity(21) = velocityFunc(vector(x3,y1,z3))
          thisOctal%cornerrho(22)      =  densityFunc(vector(x1,y2,z3))
          thisOctal%cornerVelocity(22) = velocityFunc(vector(x1,y2,z3))
          thisOctal%cornerrho(23)      =  densityFunc(vector(x2,y2,z3))
          thisOctal%cornerVelocity(23) = velocityFunc(vector(x2,y2,z3))
          thisOctal%cornerrho(24)      =  densityFunc(vector(x3,y2,z3))
          thisOctal%cornerVelocity(24) = velocityFunc(vector(x3,y2,z3))
          thisOctal%cornerrho(25)      =  densityFunc(vector(x1,y3,z3))
          thisOctal%cornerVelocity(25) = velocityFunc(vector(x1,y3,z3))
          thisOctal%cornerrho(26)      =  densityFunc(vector(x2,y3,z3))
          thisOctal%cornerVelocity(26) = velocityFunc(vector(x2,y3,z3))
          thisOctal%cornerrho(27)      =  densityFunc(vector(x3,y3,z3))
          thisOctal%cornerVelocity(27) = velocityFunc(vector(x3,y3,z3))

       else ! cylindrical 
          if (thisOctal%splitAzimuthally) then
             z1 = thisOctal%centre%z - thisOctal%subcellSize
             z2 = thisOctal%centre%z
             z3 = thisOctal%centre%z + thisOctal%subcellSize
             phi1 = thisOctal%phi - thisOctal%dPhi/2.d0
             phi2 = thisOctal%phi 
             phi3 = thisOctal%phi + thisOctal%dPhi/2.d0
             r1 = thisOctal%r - thisOctal%subcellSize
             r2 = thisOctal%r
             r3 = thisOctal%r + thisOctal%subcellSize

             ! bottom level

             thisOctal%cornerrho(1) = densityFunc(vector(r1*cos(phi1),r1*sin(phi1),z1))
             thisOctal%cornerrho(2) = densityFunc(vector(r1*cos(phi2),r1*sin(phi2),z1))
             thisOctal%cornerrho(3) = densityFunc(vector(r1*cos(phi3),r1*sin(phi3),z1))
             thisOctal%cornerrho(4) = densityFunc(vector(r2*cos(phi1),r2*sin(phi1),z1))
             thisOctal%cornerrho(5) = densityFunc(vector(r2*cos(phi2),r2*sin(phi2),z1))
             thisOctal%cornerrho(6) = densityFunc(vector(r2*cos(phi3),r2*sin(phi3),z1))
             thisOctal%cornerrho(7) = densityFunc(vector(r3*cos(phi1),r3*sin(phi1),z1))
             thisOctal%cornerrho(8) = densityFunc(vector(r3*cos(phi2),r3*sin(phi2),z1))
             thisOctal%cornerrho(9) = densityFunc(vector(r3*cos(phi3),r3*sin(phi3),z1))

             ! middle level

             thisOctal%cornerrho(10) = densityFunc(vector(r1*cos(phi1),r1*sin(phi1),z2))
             thisOctal%cornerrho(11) = densityFunc(vector(r1*cos(phi2),r1*sin(phi2),z2))
             thisOctal%cornerrho(12) = densityFunc(vector(r1*cos(phi3),r1*sin(phi3),z2))
             thisOctal%cornerrho(13) = densityFunc(vector(r2*cos(phi1),r2*sin(phi1),z2))
             thisOctal%cornerrho(14) = densityFunc(vector(r2*cos(phi2),r2*sin(phi2),z2))
             thisOctal%cornerrho(15) = densityFunc(vector(r2*cos(phi3),r2*sin(phi3),z2))
             thisOctal%cornerrho(16) = densityFunc(vector(r3*cos(phi1),r3*sin(phi1),z2))
             thisOctal%cornerrho(17) = densityFunc(vector(r3*cos(phi2),r3*sin(phi2),z2))
             thisOctal%cornerrho(18) = densityFunc(vector(r3*cos(phi3),r3*sin(phi3),z2))

             ! top level

             thisOctal%cornerrho(19) = densityFunc(vector(r1*cos(phi1),r1*sin(phi1),z3))
             thisOctal%cornerrho(20) = densityFunc(vector(r1*cos(phi2),r1*sin(phi2),z3))
             thisOctal%cornerrho(21) = densityFunc(vector(r1*cos(phi3),r1*sin(phi3),z3))
             thisOctal%cornerrho(22) = densityFunc(vector(r2*cos(phi1),r2*sin(phi1),z3))
             thisOctal%cornerrho(23) = densityFunc(vector(r2*cos(phi2),r2*sin(phi2),z3))
             thisOctal%cornerrho(24) = densityFunc(vector(r2*cos(phi3),r2*sin(phi3),z3))
             thisOctal%cornerrho(25) = densityFunc(vector(r3*cos(phi1),r3*sin(phi1),z3))
             thisOctal%cornerrho(26) = densityFunc(vector(r3*cos(phi2),r3*sin(phi2),z3))
             thisOctal%cornerrho(27) = densityFunc(vector(r3*cos(phi3),r3*sin(phi3),z3))

          else

             z1 = thisOctal%centre%z - thisOctal%subcellSize
             z2 = thisOctal%centre%z
             z3 = thisOctal%centre%z + thisOctal%subcellSize
             phi1 = thisOctal%phi - thisOctal%dPhi/2.d0
             phi2 = thisOctal%phi + thisOctal%dPhi/2.d0
             r1 = thisOctal%r - thisOctal%subcellSize
             r2 = thisOctal%r
             r3 = thisOctal%r + thisOctal%subcellSize


             ! bottom level

             thisOctal%cornerrho(1) = densityFunc(vector(r1*cos(phi1),r1*sin(phi1),z1))
             thisOctal%cornerrho(2) = densityFunc(vector(r1*cos(phi2),r1*sin(phi2),z1))
             thisOctal%cornerrho(3) = densityFunc(vector(r2*cos(phi1),r2*sin(phi1),z1))
             thisOctal%cornerrho(4) = densityFunc(vector(r2*cos(phi2),r2*sin(phi2),z1))
             thisOctal%cornerrho(5) = densityFunc(vector(r3*cos(phi1),r3*sin(phi1),z1))
             thisOctal%cornerrho(6) = densityFunc(vector(r3*cos(phi2),r3*sin(phi2),z1))

             ! middle level

             thisOctal%cornerrho(7) = densityFunc(vector(r1*cos(phi1),r1*sin(phi1),z2))
             thisOctal%cornerrho(8) = densityFunc(vector(r1*cos(phi2),r1*sin(phi2),z2))
             thisOctal%cornerrho(9) = densityFunc(vector(r2*cos(phi1),r2*sin(phi1),z2))
             thisOctal%cornerrho(10) = densityFunc(vector(r2*cos(phi2),r2*sin(phi2),z2))
             thisOctal%cornerrho(11) = densityFunc(vector(r3*cos(phi1),r3*sin(phi1),z2))
             thisOctal%cornerrho(12) = densityFunc(vector(r3*cos(phi2),r3*sin(phi2),z2))

             ! top level

             thisOctal%cornerrho(13) = densityFunc(vector(r1*cos(phi1),r1*sin(phi1),z3))
             thisOctal%cornerrho(14) = densityFunc(vector(r1*cos(phi2),r1*sin(phi2),z3))
             thisOctal%cornerrho(15) = densityFunc(vector(r2*cos(phi1),r2*sin(phi1),z3))
             thisOctal%cornerrho(16) = densityFunc(vector(r2*cos(phi2),r2*sin(phi2),z3))
             thisOctal%cornerrho(17) = densityFunc(vector(r3*cos(phi1),r3*sin(phi1),z3))
             thisOctal%cornerrho(18) = densityFunc(vector(r3*cos(phi2),r3*sin(phi2),z3))

          endif
       endif
    else
       
       
    ! we first store the values we use to assemble the position vectors
       
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       
       z1 = thisOctal%centre%z - thisOctal%subcellSize
       z2 = thisOctal%centre%z
       z3 = thisOctal%centre%z + thisOctal%subcellSize
       
       ! now store the 'base level' values
       
       thisOctal%cornerrho(1) = densityFunc(vector(x1,0.d0,z1))
       thisOctal%cornerrho(2) = densityFunc(vector(x2,0.d0,z1))
       thisOctal%cornerrho(3) = densityFunc(vector(x3,0.d0,z1))
       thisOctal%cornerrho(4) = densityFunc(vector(x1,0.d0,z2))
       thisOctal%cornerrho(5) = densityFunc(vector(x2,0.d0,z2))
       thisOctal%cornerrho(6) = densityFunc(vector(x3,0.d0,z2))
       thisOctal%cornerrho(7) = densityFunc(vector(x1,0.d0,z3))
       thisOctal%cornerrho(8) = densityFunc(vector(x2,0.d0,z3))
       thisOctal%cornerrho(9) = densityFunc(vector(x3,0.d0,z3))
    endif
667 continue
    
  END SUBROUTINE fillDensityCorners
           
  real(double) function molebenchDensity(position) result(rho)

    type(VECTOR),intent(in) :: position
    logical, save :: firsttime = .true.
    integer, parameter :: nr = 50
    real(double),save :: r(nr), nh2(nr), junk
    real(double) :: r1, t1, nh2out
    integer :: i

    if (firsttime) then
       open(31, file="model_1.dat", status="old", form="formatted") ! Model 2 in the Hogerheijde 2000 paper. 
       do i = nr,1,-1                                             
          read(31,*) r(i), nh2(i), junk
       enddo
       r = r * 1.d-10
       close(31)
       firsttime = .false.
    endif

    r1 = modulus(position)

    if(r1 > r(nr) .or. r1 < r(1)) then 
       nh2out = 1.d-20
       rho = 1.d-20 * 2. * mhydrogen
!       nh2out = nh2(1)
!       rho = nh2(1) * 2.d0 * mhydrogen
    endif

    if ((r1 > r(1)).and.(r1 < r(nr))) then
       call locate(r, nr, r1, i)
!       t2 = (r1 - r(i))/(r(i+1)-r(i)) ! linear but know its a power law so use better interpolation
       
       t1 = log(r1/r(i))/log(r(i+1)/r(i))
       
       nh2out = exp((1.d0 - t1) * log(nh2(i))  +  t1 * log(nh2(i+1)))
       rho = nh2out * 2.d0 * mhydrogen
       
    endif
  end function molebenchDensity

  real(double) function ggtauDensity(position) result(rho)
    
    type(vector) :: position
    real(double) :: r,nh20,H0,H,z,nh2

    r = sqrt(position%x**2+position%y**2) ! cylindrical
    r = r * 6.68458134e-06! (torus units to 100's of AUs)

    z = 6.68458134e-04 * position%z !in AU

    nh20 = 6.3e9
    H0 = 14.55

    H = H0 * r * sqrt(sqrt(r))  !H0*r^(5/4) r in 100s of AU and H in AU

    if(r .gt. 1.8 .and. r .lt. 8.) then
       nh2 = nh20 * (r **(-2.75)) * exp(-((z/H)**2))!(in cm-3) 
       rho = nh2  * 2. * mHydrogen
    else
       nh2 = 1.d-60
       rho = 1.d-60 * 2. * mhydrogen
    endif

  end function ggtauDensity
  
  real(double) function averagerhofromoctal(thisoctal,subcell) result(rho)
    
    type(OCTAL) :: thisoctal
    integer :: subcell

       SELECT CASE(subcell)
               
       CASE(1)
          rho =( &
               thisoctal%cornerrho( 1) + &
               thisoctal%cornerrho( 2) + &
               thisoctal%cornerrho( 4) + &
               thisoctal%cornerrho( 5) + &
               thisoctal%cornerrho(10) + &
               thisoctal%cornerrho(11) + &
               thisoctal%cornerrho(13) + &
               thisoctal%cornerrho(14)) * 1.25d-1
          
       CASE(2)
          rho =( &
               thisoctal%cornerrho( 2) + &
               thisoctal%cornerrho( 3) + &
               thisoctal%cornerrho( 5) + &
               thisoctal%cornerrho( 6) + &
               thisoctal%cornerrho(11) + &
               thisoctal%cornerrho(12) + &
               thisoctal%cornerrho(14) + &
               thisoctal%cornerrho(15)) * 1.25d-1
          
       CASE(3)
          rho =( &
               thisoctal%cornerrho( 4) + &
               thisoctal%cornerrho( 5) + &
               thisoctal%cornerrho( 7) + &
               thisoctal%cornerrho( 8) + &
               thisoctal%cornerrho(13) + &
               thisoctal%cornerrho(14) + &
               thisoctal%cornerrho(16) + &
               thisoctal%cornerrho(17)) * 1.25d-1
          
       CASE(4)
          rho =( &
               thisoctal%cornerrho( 5) + &
               thisoctal%cornerrho( 6) + &
               thisoctal%cornerrho( 8) + &
               thisoctal%cornerrho( 9) + &
               thisoctal%cornerrho(14) + &
               thisoctal%cornerrho(15) + &
               thisoctal%cornerrho(17) + &
               thisoctal%cornerrho(18)) * 1.25d-1
          
       CASE(5)
          rho =( &
               thisoctal%cornerrho(10) + &
               thisoctal%cornerrho(11) + &
               thisoctal%cornerrho(13) + &
               thisoctal%cornerrho(14) + &
               thisoctal%cornerrho(19) + &
               thisoctal%cornerrho(20) + &
               thisoctal%cornerrho(22) + &
               thisoctal%cornerrho(23)) * 1.25d-1
          
       CASE(6)
          rho =( &
               thisoctal%cornerrho(11) + &
               thisoctal%cornerrho(12) + &
               thisoctal%cornerrho(14) + &
               thisoctal%cornerrho(15) + &
               thisoctal%cornerrho(20) + &
               thisoctal%cornerrho(21) + &
               thisoctal%cornerrho(23) + &
               thisoctal%cornerrho(24)) * 1.25d-1
          
       CASE(7)
          rho =( &
               thisoctal%cornerrho(13) + &
               thisoctal%cornerrho(14) + &
               thisoctal%cornerrho(16) + &
               thisoctal%cornerrho(17) + &
               thisoctal%cornerrho(22) + &
               thisoctal%cornerrho(23) + &
               thisoctal%cornerrho(25) + &
               thisoctal%cornerrho(26)) * 1.25d-1
          
       CASE(8)
          rho =( &
               thisoctal%cornerrho(14) + &
               thisoctal%cornerrho(15) + &
               thisoctal%cornerrho(17) + &
               thisoctal%cornerrho(18) + &
               thisoctal%cornerrho(23) + &
               thisoctal%cornerrho(24) + &
               thisoctal%cornerrho(26) + &
               thisoctal%cornerrho(27)) * 1.25d-1
          
       end select

  end function averagerhofromoctal

  subroutine AMRHydroInterpFromParent(thisOctal, grid)
    type(OCTAL), pointer :: thisOctal
    type(GRIDTYPE) :: grid
    integer :: parentSubcell
    type(OCTAL), pointer :: parentOctal, probeOctal1, probeOctal2
    integer :: probeSubcell1, probeSubcell2
    type(VECTOR) :: parentCentre, probe1, probe2, rVec
    real(double) :: x, x1, x2
    real(double) :: y1rho, y1rhou, y1rhoe, y1energy
    real(double) :: y2rho, y2rhou, y2rhoe, y2energy

    real(double) :: r, u, v, y1, y2, y
    real(double) :: v11(4), v12(4), v21(4), v22(4)
    integer :: iSubcell

    parentOctal => thisOctal%parent
    parentSubcell = thisOctal%parentSubcell
    parentCentre = subcellCentre(parentOctal, parentSubcell)

    thisOctal%refinedLastTime = .true.
    thisOctal%boundaryCondition = parentOctal%boundaryCondition(parentSubcell)
    thisOctal%boundaryCell = parentOctal%boundaryCell(parentSubcell)
    thisOctal%gamma = parentOctal%gamma(parentSubcell)
    thisOctal%iEquationOfState = parentOctal%iEquationofState(parentSubcell)

    thisOctal%iAnalyticalVelocity = parentOctal%iAnalyticalVelocity(parentSubcell)
    thisOctal%velocity = parentOctal%velocity(parentSubcell)
    thisOctal%dustTypeFraction = parentOctal%dustTypeFraction
    if (thisOctal%oneD) then

       thisOctal%rhov = 0.d0
       thisOctal%rhow = 0.d0

       probe1 = VECTOR(parentOctal%xMin - 0.01d0*grid%halfSmallestSubcell, 0.d0, 0.d0)
       if (inOctal(grid%octreeRoot, probe1)) then
          probeOctal1 => thisOctal
          call findSubcellLocal(probe1, probeOctal1, probeSubcell1)
          rVec = subcellCentre(probeOctal1, probeSubcell1)
          x1 = rVec%x
          y1rho = probeOctal1%rho(probeSubcell1)
          y1rhou = probeOctal1%rhou(probeSubcell1)
          y1rhoe = probeOctal1%rhoe(probeSubcell1)
          y1energy = probeOctal1%energy(probeSubcell1)
       else
          rVec = subcellCentre(thisOctal,1)
          x1 = rVec%x
          y1rho = parentOctal%rho(parentSubcell)
          y1rhou = parentOctal%rhou(parentSubcell)
          y1rhoe = parentOctal%rhoe(parentSubcell)
          y1energy = parentOctal%energy(parentSubcell)
       endif
       probe2 = VECTOR(parentOctal%xMax + 0.01d0*grid%halfSmallestSubcell, 0.d0, 0.d0)
       if (inOctal(grid%octreeRoot, probe2)) then
          probeOctal2 => thisOctal
          call findSubcellLocal(probe2, probeOctal2, probeSubcell2)
          rVec = subcellCentre(probeOctal2, probeSubcell2)
          x2 = rVec%x
          y2rho = probeOctal2%rho(probeSubcell2)
          y2rhou = probeOctal2%rhou(probeSubcell2)
          y2rhoe = probeOctal2%rhoe(probeSubcell2)
          y2energy = probeOctal2%energy(probeSubcell2)
       else
          rVec = subcellCentre(thisOctal,2)
          x2 = rVec%x
          y2rho = parentOctal%rho(parentSubcell)
          y2rhou = parentOctal%rhou(parentSubcell)
          y2rhoe = parentOctal%rhoe(parentSubcell)
          y2energy = parentOctal%energy(parentSubcell)
       endif
       
       do iSubcell = 1, thisOctal%maxChildren
          rVec = subcellCentre(thisOctal, iSubcell)
          x = rvec%x
          thisOctal%rho(iSubcell)  = y1rho +  (y2rho -  y1rho)  * (x - x1)/(x2 - x1)
          thisOctal%rhou(iSubcell) = (y1rhou + (y2rhou - y1rhou) * (x - x1)/(x2 - x1))
          thisOctal%rhoe(iSubcell) = (y1rhoe + (y2rhoe - y1rhoe) * (x - x1)/(x2 - x1))
          thisOctal%energy(iSubcell) = y1energy + (y2energy - y1energy) * (x - x1)/(x2 - x1)
       enddo

    else if (thisOctal%twoD) then
       r = parentOctal%subcellSize+grid%halfSmallestSubcell * 1.d-3
       x1 = parentOctal%xmin
       x2 = parentOctal%xMax
       y1 = parentOctal%zMin
       y2 = parentOctal%zMax
       probeoctal1 => thisOctal

       v11(1) = parentOctal%rho(parentSubcell)
       v11(2) = parentOctal%rhoe(parentSubcell)
       v11(3) = parentOctal%rhou(parentSubcell)
       v11(4) = parentOctal%rhow(parentSubcell)
       probe1 = parentOctal%centre + r * vector(-1.d0, 0.d0, -1.d0)
       if (inOctal(grid%octreeRoot, probe1)) then
          call findSubcellLocal(probe1, probeOctal1, probeSubcell1)
!          if (octalOnThread(probeOctal1, probeSubcell1, myRankGlobal)) then
             v11(1) = probeOctal1%rho(probesubcell1)
             v11(2) = probeOctal1%rhoe(probesubcell1)
             v11(3) = probeOctal1%rhou(probesubcell1)
             v11(4) = probeOctal1%rhow(probesubcell1)
!          endif
       endif


       v12(1) = parentOctal%rho(parentSubcell)
       v12(2) = parentOctal%rhoe(parentSubcell)
       v12(3) = parentOctal%rhou(parentSubcell)
       v12(4) = parentOctal%rhow(parentSubcell)
       probe1 = parentOctal%centre + r * vector(-1.d0, 0.d0, 1.d0)
       if (inOctal(grid%octreeRoot, probe1)) then
          call findSubcellLocal(probe1, probeOctal1, probeSubcell1)
!          if (octalOnThread(probeOctal1, probeSubcell1, myRankGlobal)) then
             v12(1) = probeOctal1%rho(probesubcell1)
             v12(2) = probeOctal1%rhoe(probesubcell1)
             v12(3) = probeOctal1%rhou(probesubcell1)
             v12(4) = probeOctal1%rhow(probesubcell1)
!          endif
       endif

       v22(1) = parentOctal%rho(parentSubcell)
       v22(2) = parentOctal%rhoe(parentSubcell)
       v22(3) = parentOctal%rhou(parentSubcell)
       v22(4) = parentOctal%rhow(parentSubcell)
       probe1 = parentOctal%centre + r * vector(1.d0, 0.d0, 1.d0)
       if (inOctal(grid%octreeRoot, probe1)) then
          call findSubcellLocal(probe1, probeOctal1, probeSubcell1)
!          if (octalOnThread(probeOctal1, probeSubcell1, myRankGlobal)) then
             v22(1) = probeOctal1%rho(probesubcell1)
             v22(2) = probeOctal1%rhoe(probesubcell1)
             v22(3) = probeOctal1%rhou(probesubcell1)
             v22(4) = probeOctal1%rhow(probesubcell1)
!          endif
       endif

       v21(1) = parentOctal%rho(parentSubcell)
       v21(2) = parentOctal%rhoe(parentSubcell)
       v21(3) = parentOctal%rhou(parentSubcell)
       v21(4) = parentOctal%rhow(parentSubcell)
       probe1 = parentOctal%centre + r * vector(1.d0, 0.d0, -1.d0)
       if (inOctal(grid%octreeRoot, probe1)) then
          call findSubcellLocal(probe1, probeOctal1, probeSubcell1)
!          if (octalOnThread(probeOctal1, probeSubcell1, myRankGlobal)) then
             v21(1) = probeOctal1%rho(probesubcell1)
             v21(2) = probeOctal1%rhoe(probesubcell1)
             v21(3) = probeOctal1%rhou(probesubcell1)
             v21(4) = probeOctal1%rhow(probesubcell1)
!          endif
       endif

       do iSubcell = 1, thisOctal%maxChildren
          rVec = subcellCentre(thisOctal, iSubcell)
          x = rvec%x
          y = rvec%z

          u = (x - x1)/(x2 - x1)
          v = (y - y1)/(y2 - y1)

          thisOctal%rho(iSubcell)  = v11(1) * (1.d0-u) * (1.d0-v) + &
                                     v12(1) * (1.d0-u) * (     v) + &
                                     v22(1) * (     u) * (     v) + &
                                     v21(1) * (     u) * (1.d0-v)
          write(*,*) "u ", u, " v ", v, "rho corners ", v11(1), v12(1), v22(1), v21(1), "rho ",thisOctal%rho(iSubcell)

          thisOctal%rhoe(iSubcell)  =v11(2) * (1.d0-u) * (1.d0-v) + &
                                     v12(2) * (1.d0-u) * (     v) + &
                                     v22(2) * (     u) * (     v) + &
                                     v21(2) * (     u) * (1.d0-v)

          thisOctal%rhou(iSubcell)  =v11(3) * (1.d0-u) * (1.d0-v) + &
                                     v12(3) * (1.d0-u) * (     v) + &
                                     v22(3) * (     u) * (     v) + &
                                     v21(3) * (     u) * (1.d0-v)

          thisOctal%rhow(iSubcell)  =v11(4) * (1.d0-u) * (1.d0-v) + &
                                     v12(4) * (1.d0-u) * (     v) + &
                                     v22(4) * (     u) * (     v) + &
                                     v21(4) * (     u) * (1.d0-v)
       enddo



    else
       call writeFatal("AMRhydrointerpfromparent not implemented (apart from 1-d case)")
    endif
  end subroutine AMRHydroInterpFromParent

  function getPhiValue(thisOctal, subcell, geometry) result(phi)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    character(len=*) :: geometry
    real(double) :: phi
    type(VECTOR) :: rVec

    rVec = subcellCentre(thisOctal,subcell)
    select case (geometry)
       case("rtaylor")
          phi = rvec%z * 0.1d0
       case DEFAULT
          call writeWarning("Unrecognized geometry: "//trim(geometry)//" in get phi value. zero returned.")
          phi = 0.d0
    end select
  end function getPhiValue

  subroutine tauAlongPathCMF(grid, rVec, direction, tau, tauMax, startOctal, startSubcell)
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec, direction, currentPosition, beforeVec, afterVec
    real(double), intent(out) :: tau
    real(double) :: distToNextCell
    real(double), optional :: tauMax
    type(OCTAL), pointer :: thisOctal, sOctal
    type(OCTAL), pointer, optional :: startOctal
    integer, optional :: startSubcell
    real(double) :: fudgeFac = 1.d-1
    real(double) :: kappaSca, kappaAbs, kappaExt
    integer :: subcell
    kappaAbs = 0.d0; kappaSca = 0.d0
    tau = 0.d0
    currentPosition = rVec

    if (PRESENT(startOctal)) then
       thisOctal => startOctal
       subcell = startSubcell
    else
       CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)
    endif

    do while (inOctal(grid%octreeRoot, currentPosition))

       call findSubcellLocal(currentPosition, thisOctal,subcell)
       kappaExt = thisOctal%kappaSca(subcell,1) + thisOctal%kappaAbs(subcell,1)
!       kappaExt = thisOctal%kappaAbs(subcell,1)
!        kappaExt = sqrt(thisOctal%kappaSca(subcell,1) * thisOctal%kappaAbs(subcell,1))
       sOctal => thisOctal
       call distanceToCellBoundary(grid, currentPosition, direction, DisttoNextCell, sOctal)
       
       beforeVec = VECTOR(currentPosition%x, currentPosition%y,0.d0)
       currentPosition = currentPosition + (distToNextCell+fudgeFac*grid%halfSmallestSubcell)*direction
       afterVec = VECTOR(currentPosition%x, currentPosition%y,0.d0)
       

       tau = tau + distToNextCell*kappaExt
       if (PRESENT(tauMax)) then
          if (tau > tauMax) exit
       endif
    end do
  end subroutine tauAlongPathCMF


  recursive subroutine addWarpedDisc(thisOctal)
    use inputs_mod, only : ttauriROuter, hOverR, ttauriRinner
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: r, phi, height, cellsize
    type(VECTOR) :: rVec
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call addWarpedDisc(child)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal, subcell)
          r = sqrt(rVec%x**2 + rVec%y**2)
          phi = atan2(rVec%y, rVec%x)
          cellsize = thisOctal%subcellSize
          if (thisOctal%cylindrical) then
             height = discHeightFunc(phi, hOverR) * r
          else
             height = 0.d0
          endif
          if ((rvec%z < 0.d0).and.((rvec%z+1.0001d0*cellSize/2.d0)>0.d0)) then
             if (r > ttaurirInner/1.d10) then
                thisOctal%rho(subcell) = 1.d0
                thisOctal%inflow(subcell) = .false.
             endif
          endif
          if (((r-cellsize/2.d0) < (ttaurirOuter/1.d10)).and.( (r+cellsize/2.d0) > (ttaurirouter/1.d10))) then
             if (rVec%z < height) then
                thisOctal%rho(subcell) = 1.d0
                thisOctal%inflow(subcell) = .false.
             endif
          endif
       endif
    enddo
  end subroutine addWarpedDisc
  
          
  function discHeightFunc(phi, hOverR) result (height)

    real(double) :: phi, hOverR, height, phiShift
    
    phiShift = phi 
    if (phiShift > twoPi) phiShift = phiShift - twoPi
    height = hOverR * abs(cos(phiShift/2.d0))
  end function discHeightFunc
    
  function octalOnThread(thisOctal, subcell, myRank) result(check)
    use inputs_mod, only : splitOverMPI
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer :: myRank
    logical :: check
    integer :: nFirstLevel

    check = .true.

    if (splitOverMPI) then

    if (thisOctal%twoD) then
       if (nHydroThreadsGlobal == 4) then
          if (thisOctal%mpiThread(subcell) == myRank) then
             check = .true.
          else
             check = .false.
          endif
       endif
       if (nHydroThreadsGlobal == 16) then
          nFirstLevel = (myRank-1) / 4 + 1
          if (thisOctal%nDepth == 1) then
             if (thisOctal%mpiThread(subcell) == nFirstLevel) then
                check = .true.
             else
                check = .false.
             endif
          else
             if (thisOctal%mpiThread(subcell) == myRank) then
                check = .true.
             else
                check = .false.
             endif
          endif
       endif

       if (nHydroThreadsGlobal == 64) then
          nFirstLevel = (myRank-1) / 16 + 1
          if (thisOctal%nDepth == 1) then
             if (thisOctal%mpiThread(subcell) == nFirstLevel) then
                check = .true.
             else
                check = .false.
             endif
          else if (thisOctal%nDepth == 2) then
             nFirstLevel = (myRank-1) / 4  + 1
             if (thisOctal%mpiThread(subcell) == nFirstLevel) then
                check = .true.
             else
                check = .false.
             endif
          else
             if (thisOctal%mpiThread(subcell) == myRank) then
                check = .true.
             else
                check = .false.
             endif
          endif
       endif
    endif

    if (thisOctal%threeD) then
       if (nHydroThreadsGlobal == 8) then

          if (thisOctal%mpiThread(subcell) == myRank) then
             check = .true.
          else
             check = .false.
          endif
       endif

       if (nHydroThreadsGlobal == 64) then
          nFirstLevel = (myRank-1) / 8 + 1
          if (thisOctal%nDepth == 1) then
             if (thisOctal%mpiThread(subcell) == nFirstLevel) then
                check = .true.
             else
                check = .false.
             endif
          else
             if (thisOctal%mpiThread(subcell) == myRank) then
                check = .true.
             else
                check = .false.
             endif
          endif
!          write(*,*) "thread ", myrankGlobal, " depth ",thisOctal%ndepth, " mpithread ", thisOctal%mpiThread(subcell), check
       endif

       if (nHydroThreadsGlobal == 512) then
          nFirstLevel = (myRank-1) / 64 + 1
          if (thisOctal%nDepth == 1) then
             if (thisOctal%mpiThread(subcell) == nFirstLevel) then
                check = .true.
             else
                check = .false.
             endif
          else if (thisOctal%nDepth == 2) then
             nFirstLevel = (myRank-1) / 8  + 1
             if (thisOctal%mpiThread(subcell) == nFirstLevel) then
                check = .true.
             else
                check = .false.
             endif
          else
             if (thisOctal%mpiThread(subcell) == myRank) then
                check = .true.
             else
                check = .false.
             endif
          endif
!          write(*,*) "thread ", myrankGlobal, " depth ",thisOctal%ndepth, " mpithread ", thisOctal%mpiThread(subcell), check
       endif



    endif

    if (thisOctal%oned) then
       if (nHydroThreadsGlobal == 2) then
          if (thisOctal%mpiThread(subcell) == myRank) then
             check = .true.
          else
             check = .false.
          endif
       endif
       if (nHydroThreadsGlobal == 4) then
          nFirstLevel = (myRank-1) / 2 + 1
          if (thisOctal%nDepth == 1) then
             if (thisOctal%mpiThread(subcell) == nFirstLevel) then
                check = .true.
             else
                check = .false.
             endif
          else
             if (thisOctal%mpiThread(subcell) == myRank) then
                check = .true.
             else
                check = .false.
             endif
          endif
       endif
    endif
 endif
  end function octalOnThread

  !
  !
  !

END MODULE amr_mod
