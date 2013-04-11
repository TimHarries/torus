module amr_utils_mod
 

  use vector_mod
  use messages_mod
  USE constants_mod
  USE octal_mod
  use gridtype_mod, only: gridtype
  use mpi_global_mod, only: myRankGlobal
  use parallel_mod, only : torus_abort
  use utils_mod, only : solveQuadDble
  public

  contains


    real(double) function gridArea(grid)
      use inputs_mod, only : cylindrical, spherical
      type(GRIDTYPE) :: grid

      if (grid%octreeRoot%oneD.and.spherical) then
         gridArea = fourPi*(grid%octreeRoot%subcellSize*2.d0)*1.d20
      endif

      if (grid%octreeRoot%oneD.and.(.not.spherical)) then
         gridArea = (grid%octreeRoot%subcellSize*2.d0)**2*1.d20
      endif
      if (grid%octreeRoot%twoD) then
         gridArea = twoPi * (2.d0*grid%octreeRoot%subcellSize)**2 + &
              2.d0*pi*(2.d0*grid%octreeRoot%subcellSize)**2
         gridArea = gridArea * 1.d20
      endif
      if (grid%octreeRoot%threed) then
         if (cylindrical) then
            gridArea = twoPi * (2.d0*grid%octreeRoot%subcellSize)**2 + &
                 2.d0*pi*(2.d0*grid%octreeRoot%subcellSize)**2
            gridArea = gridArea * 1.d20
         else
            gridArea = 6.d0 * (2.d0*grid%octreeRoot%subcellSize)**2 * 1.d20
         endif
      endif
    end function gridArea



  function distanceToGridFromOutside(grid, posVec, direction, hitGrid) result (tval)
    use inputs_mod, only : suppressWarnings, spherical
    type(GRIDTYPE) :: grid
    type(VECTOR) :: subcen, direction, posVec, point, hitVec, rdirection, xhat
    type(OCTAL), pointer :: thisOctal
    real(double) :: tval
    real(double) :: distTor1, theta, mu
    real(double) :: distToRboundary, compz,currentZ
    real(double) :: distToZboundary
    type(VECTOR) ::  zHat, rhat
    real(double) ::  compx, gridRadius
    real(oct) :: t(6),denom(6), r, r1, r2, d, cosmu,x1,x2
    integer :: i,j
    logical :: ok, thisOk(6)
    logical, optional :: hitGrid
    logical :: debug
    
    real(double) :: subcellsize, r1overd
    type(VECTOR) :: normdiff, test

    if (PRESENT(hitGrid)) hitGrid = .true.
   tval = HUGE(tval)

   point = posVec

   subcen = grid%OctreeRoot%centre

   thisOctal => grid%octreeRoot

   if (thisOctal%oneD.and.spherical) then
   
      r1 = thisOctal%subcellSize*2.d0
      d = modulus(point)
      rHat = (-1.d0)*posVec
      call normalize(rhat)
      theta = asin(max(-1.d0,min(1.d0,r1 / d)))
      cosmu = rHat.dot.direction
      mu = acos(max(-1.d0,min(1.d0,cosmu)))
      distTor1 = 1.e30
      if (mu  < theta ) then
         call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
         if (.not.ok) then
            if (PRESENT(hitGrid)) then
               hitGrid = .false.
            else
               write(*,*) "Quad solver failed in distanceToGridFromOutside 1d"
               x1 = thisoctal%subcellSize
               x2 = 0.d0
            endif
         endif
         tval = min(x1,x2)
      endif
      goto 666
   endif

   if (thisOctal%oneD.and.(.not.spherical)) then
   
      r1 = subcen%x - thisOctal%subcellSize
      r2 = subcen%x + thisOctal%subcellSize

      if (posVec%x > r2) then
         if (direction%x > 0.d0) then
            hitgrid = .false.
            tval = 1.d30
         else
            tval = (posVec%x - r2)/abs(direction%x)
         endif
         goto 666
      endif

      if (posVec%x < r1) then
         if (direction%x < 0.d0) then
            hitgrid = .false.
            tval = 1.d30
         else
            tval = (r1 - posVec%x)/abs(direction%x)
         endif
         goto 666
      endif
      write(*,*) "error in 1d spherical case"


   endif

    if (thisOctal%twod.and.(cart2d)) then

          ! cube

         ok = .true.

         subcellsize = thisOctal%subcellSize

         if(direction%x .ne. 0.d0) then            
            denom(1) = 1.d0 / direction%x
         else
            denom(1) = 0.d0
         endif
         denom(4) = -denom(1)

!         if(direction%y .ne. 0.d0) then            
!            denom(2) = 1.d0 / direction%y
!         else
!            denom(2) = 0.d0
!         endif
!         denom(5) = -denom(2)

         if(direction%z .ne. 0.d0) then            
            denom(3) = 1.d0 / direction%z
         else
            denom(3) = 0.d0
         endif
         denom(6) = -denom(3)

         normdiff = subcen - posvec
         t = 0.d0
         t(1) =  (normdiff%x + subcellsize) * denom(1)
 !        t(2) =  (normdiff%y + subcellsize) * denom(2)
         t(3) =  (normdiff%z + subcellsize) * denom(3)
         t(4) =  (normdiff%x - subcellsize) * denom(1)
!         t(5) =  (normdiff%y - subcellsize) * denom(2)
         t(6) =  (normdiff%z - subcellsize) * denom(3)

         thisOk = .true.

         do i = 1, 6
            if (denom(i) .ge. 0.d0) thisOK(i) = .false.
            if (t(i) < 0.) thisOk(i) = .false.
         enddo
         
         j = 0
         do i = 1, 6
            if (thisOk(i)) j=j+1
         enddo
         
         if (j == 0) ok = .false.
         
         if (.not.ok) then
            if (PRESENT(hitGrid)) then
               hitGrid = .false.
            else
               write(*,*) "Error: j=0 (no intersection???) in lucy_mod::distancetoGridFromOutside. 2d"
               write(*,*) direction%x,direction%y,direction%z
               write(*,*) t(1:6)
            endif
         endif
         
         tval = maxval(t, mask=thisOk)
         
!         write(*,*) t(1:6),thisOK(1:6)

         if (tval == 0.) then
            if (PRESENT(hitGrid)) then
               hitGrid = .false.
            else
               write(*,*) " tval=0 (no intersection???) in lucy_mod::distancetoGridFromOutside. "
               write(*,*) posVec
               write(*,*) direction%x,direction%y,direction%z
               write(*,*) t(1:6)
               stop
            endif
         endif

         test = posVec + (tval+1.d-3*grid%halfSmallestSubcell) * direction
         if (.not.inOctal(grid%octreeRoot, test)) then
            if (PRESENT(hitGrid)) then
               hitGrid = .false.
            endif
         endif
         goto 666
      end if


    if (thisOctal%threed.and.(.not.thisOctal%cylindrical)) then

          ! cube

         ok = .true.

         subcellsize = thisOctal%subcellSize

         if(direction%x .ne. 0.d0) then            
            denom(1) = 1.d0 / direction%x
         else
            denom(1) = 0.d0
         endif
         denom(4) = -denom(1)

         if(direction%y .ne. 0.d0) then            
            denom(2) = 1.d0 / direction%y
         else
            denom(2) = 0.d0
         endif
         denom(5) = -denom(2)

         if(direction%z .ne. 0.d0) then            
            denom(3) = 1.d0 / direction%z
         else
            denom(3) = 0.d0
         endif
         denom(6) = -denom(3)

         normdiff = subcen - posvec

         t(1) =  (normdiff%x + subcellsize) * denom(1)
         t(2) =  (normdiff%y + subcellsize) * denom(2)
         t(3) =  (normdiff%z + subcellsize) * denom(3)
         t(4) =  (normdiff%x - subcellsize) * denom(1)
         t(5) =  (normdiff%y - subcellsize) * denom(2)
         t(6) =  (normdiff%z - subcellsize) * denom(3)

         thisOk = .true.

         do i = 1, 6
            if (denom(i) .ge. 0.d0) thisOK(i) = .false.
            if (t(i) < 0.) thisOk(i) = .false.
         enddo
         
         j = 0
         do i = 1, 6
            if (thisOk(i)) j=j+1
         enddo
         
         if (j == 0) ok = .false.
         
         if (.not.ok) then
            if (PRESENT(hitGrid)) then
               hitGrid = .false.
            else
               write(*,*) "Error: j=0 (no intersection???) in lucy_mod::distancetoGridFromOutside. "
               write(*,*) direction%x,direction%y,direction%z
               write(*,*) t(1:6)
            endif
         endif
         
         tval = maxval(t, mask=thisOk)
         
!         write(*,*) t(1:6),thisOK(1:6)

         if (tval == 0.) then
            if (PRESENT(hitGrid)) then
               hitGrid = .false.
            else
               write(*,*) " tval=0 (no intersection???) in lucy_mod::distancetoGridFromOutside. "
               write(*,*) posVec
               write(*,*) direction%x,direction%y,direction%z
               write(*,*) t(1:6)
               stop
            endif
         endif

         test = posVec + (tval+1.d-3*grid%halfSmallestSubcell) * direction
         if (.not.inOctal(grid%octreeRoot, test)) then
            if (PRESENT(hitGrid)) then
               hitGrid = .false.
            endif
         endif
         

      else


! now look at the cylindrical case

         ! first do the inside and outside curved surfaces

         if (thisOctal%twoD) then
            r = sqrt(subcen%x**2 + subcen%y**2)
            r1 = r + thisOctal%subcellSize
            gridRadius = r1
         else
            gridRadius = 2.d0*thisOctal%subcellSize
            r1 = gridRadius
         endif


         d = sqrt(point%x**2+point%y**2)
         xHat = (-1.d0)*VECTOR(point%x, point%y,0.d0)
         call normalize(xHat)
         rDirection = VECTOR(direction%x, direction%y, 0.d0)
         compX = modulus(rDirection)
         call normalize(rDirection)
               
         r1overd = 1.d30
         if (d /= 0.d0) r1Overd = r1/d
         theta = asin(max(-1.d0,min(1.d0,r1overd)))
         cosmu = xHat.dot.rdirection
         mu = acos(max(-1.d0,min(1.d0,cosmu)))
         distTor1 = 1.e30
         if (compx /= 0.d0) then
            if (mu  < theta ) then
               call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
               if (.not.ok) then
                  if (PRESENT(hitGrid)) then
                     hitGrid = .false.
                  else
                     write(*,*) "Quad solver failed in distanceToGridFromOutside"
                     x1 = thisoctal%subcellSize
                     x2 = 0.d0
                  endif
               endif
               distTor1 = min(x1,x2)/CompX
               hitVec = posVec + distToR1 * direction
               if (abs(hitVec%z) > thisOctal%subcellSize) distToR1 = 1.d30
            endif
         endif
         
         distToRboundary = distTor1
         if (distToRboundary < 0.d0) then
            distToRboundary = 1.e30
         endif

         ! now do the upper and lower (z axis) surfaces

         zHat = VECTOR(0.d0, 0.d0, 1.d0)
         compZ = zHat.dot.direction
         currentZ = point%z

         debug = .false.
         if(debug) then
            write(*,*) "subcen",subcen
            write(*,*) "thisoctal%subcellsize",thisoctal%subcellsize
            write(*,*) "point", point
            write(*,*) "posvec", posvec
            write(*,*) "direction",direction
         endif

         if (compZ /= 0.d0 ) then
            if (compZ > 0.d0) then
               distToZboundary= (subcen%z - thisOctal%subcellsize - currentZ ) / compZ                              
               hitVec = posvec + disttoZboundary * direction
               if (sqrt(hitVec%x**2 + hitVec%y**2) > gridRadius) distToZboundary = 1.d30
            else
               distToZboundary = abs((subcen%z + thisOctal%subcellsize - currentZ ) / compZ)
               hitVec = posvec + disttoZboundary * direction
               if (sqrt(hitVec%x**2 + hitVec%y**2) > gridRadius) distToZboundary = 1.d30
            endif
         else
            disttoZboundary = 1.e30
         endif
      
         tVal = min(distToZboundary, distToRboundary)
         if(.not. suppressWarnings) then

!            if (tVal > 1.e29) then
!               write(*,*) "Cylindrical"
!               write(*,*) "posVec",posvec
!               write(*,*) "direction",direction
!               write(*,*) tVal,compX,compZ, distToZboundary,disttorboundary
!               write(*,*) "subcen",subcen
!               write(*,*) "z",currentZ
!               write(*,*) "x1,x2",x1,x2
!            endif

            if (tval < 0.) then
               if (PRESENT(hitGrid)) then
                  hitGrid = .false.
               else
                  write(*,*) "Cylindrical"
                  write(*,*) tVal,distToZboundary,disttorboundary
                  write(*,*) "subcen",subcen
                  write(*,*) "z",currentZ
               endif
            endif
         endif
         test = posVec + (tval+1.d-3*grid%halfSmallestSubcell) * direction
         if (.not.inOctal(grid%octreeRoot, test)) then
            if (PRESENT(hitGrid)) then
               hitGrid = .false.
            endif
         endif

      endif
      
666   continue
    end function distanceToGridFromOutside

  SUBROUTINE findSubcellTD(point,currentOctal,resultOctal,subcell)
  ! finds the octal (and that octal's subcell) containing a point.
  !   only searches in downwards direction (TD = top-down) , so
  !   probably best to start from root of tree

    use inputs_mod, only : hydrodynamics, cylindricalHydro, spherical

    IMPLICIT NONE
    TYPE(vector), INTENT(IN) :: point
    type(vector) :: point_local
    TYPE(octal), POINTER :: currentOctal
    TYPE(octal), POINTER :: resultOctal
    INTEGER, INTENT(OUT) :: subcell

    if (currentoctal%threeD) then
       point_local = point
    elseif (currentoctal%twoD) then
       if(cart2d) then
          point_local = point
       elseif (.not.cylindricalHydro) then
          point_local = projectToXZ(point)
       else
          point_local = point
       endif
    else !oneD
       if (spherical) then
          point_local = VECTOR(modulus(point), 0.d0, 0.d0)
       else
          point_local = point
       endif
    end if

    CALL findSubcellTDprivate(point_local,currentOctal,resultOctal,subcell)


  contains


  RECURSIVE SUBROUTINE findSubcellTDPrivate(point,currentOctal,resultOctal,subcell)
  ! finds the octal (and that octal's subcell) containing a point.
  !   only searches in downwards direction (TD = top-down) , so
  !   probably best to start from root of tree

    IMPLICIT NONE

    TYPE(vector), INTENT(IN) :: point
    TYPE(octal), POINTER :: currentOctal
    TYPE(octal), POINTER :: resultOctal
    INTEGER, INTENT(OUT) :: subcell
    TYPE(octal), POINTER :: child

    INTEGER :: i

    resultOctal => currentOctal
    subcell = whichSubcell(currentOctal,point)

    IF ( currentOctal%hasChild(subcell) ) THEN
      ! search the index to see where it is stored
      DO i = 1, currentOctal%maxChildren, 1
        IF ( currentOctal%indexChild(i) == subcell ) THEN
          child => currentOctal%child(i)
          CALL findSubcellTDprivate(point,child,resultOctal,subcell)
          EXIT
        END IF
      END DO
    END IF

  END SUBROUTINE findSubcellTDPrivate

  END SUBROUTINE findSubcellTD

  RECURSIVE SUBROUTINE findSubcellTDLevel(point,currentOctal,resultOctal,subcell,nDepth)
  ! finds the octal (and that octal's subcell) containing a point.
  !   only searches in downwards direction (TD = top-down) , so
  !   probably best to start from root of tree

    IMPLICIT NONE

    TYPE(vector), INTENT(IN) :: point
    TYPE(octal), POINTER :: currentOctal
    TYPE(octal), POINTER :: resultOctal
    INTEGER, INTENT(OUT) :: subcell
    integer :: nDepth
    TYPE(octal), POINTER :: child

    INTEGER :: i

    resultOctal => currentOctal
    subcell = whichSubcell(currentOctal,point)

    IF ( currentOctal%hasChild(subcell).and.(currentOctal%nDepth < nDepth) ) THEN
      ! search the index to see where it is stored
      DO i = 1, currentOctal%maxChildren, 1
        IF ( currentOctal%indexChild(i) == subcell ) THEN
          child => currentOctal%child(i)
          CALL findSubcellTDLevel(point,child,resultOctal,subcell,ndepth)
          EXIT
        END IF
      END DO
    END IF

  END SUBROUTINE findSubcellTDLevel


  SUBROUTINE findSubcellLocal(point,thisOctal,subcell,  prob)
    ! finds the octal (and that octal's subcell) containing a point.
    !   starts searching from the current octal, and goes up and down the
    !   tree as needed to find the correct octal.
    use inputs_mod, only : hydrodynamics, suppresswarnings, cylindricalHydro, spherical

    IMPLICIT NONE
    TYPE(vector), INTENT(IN) :: point
    TYPE(vector) :: point_local
    TYPE(octal),POINTER    :: thisOctal
    INTEGER, INTENT(OUT)   :: subcell
    LOGICAL, INTENT(OUT),optional   :: prob
    
    LOGICAL                :: haveDescended    ! see comments below
    LOGICAL                :: boundaryProblem  ! see comments below
    
                             
    haveDescended = .FALSE.   ! if the 'point' lies very close to an 
    boundaryProblem = .FALSE. !   boundary, the program may go into 
                              !   a loop going up and down the tree.
                              ! we will keep track of the progress of
                              !   the search using these flags.
                             
    if (thisOctal%twoD) then
       if(cart2d) then
          point_local = point
       else if (.not.cylindricalHydro) then
          point_local = projectToXZ(point)
       else
          point_local = point
       endif


!       if (.not.cylindricalHydro) then
!          point_local = projectToXZ(point)
!       else
!          point_local = point
!       endif
    else
       point_local = point
    endif
    if (thisOctal%oneD) then
       if (spherical) then
          point_local = VECTOR(modulus(point), 0.d0, 0.d0)
       else
          point_local = point
       endif
    endif


    CALL findSubcellLocalPrivate(point_local,thisOctal,subcell,&
                                 haveDescended,boundaryProblem)

    if (present(prob)) then
      prob = boundaryProblem
    else
      if (boundaryProblem .and. .not. suppresswarnings) then
        call torus_abort("Torus aborting due to panic in findSubcellLocal")
     else
        CALL findSubcellTD(point_local,thisoctal,thisOctal,subcell)
     endif
    endif
                                 
  CONTAINS

    RECURSIVE SUBROUTINE findSubcellLocalPrivate(point,thisOctal,subcell,&
                                                 haveDescended,boundaryProblem)
      use inputs_mod, only : suppressWarnings
      TYPE(vector), INTENT(IN) :: point
      TYPE(octal),POINTER    :: thisOctal
      real(double) :: phi, phimax, phimin,r
      INTEGER, INTENT(OUT)   :: subcell
      LOGICAL, INTENT(INOUT) :: haveDescended
      LOGICAL, INTENT(INOUT) :: boundaryProblem
!      type(vector) :: rVec
      INTEGER :: i
      
      IF ( inOctal(thisOctal,point,alreadyRotated=.true.) ) THEN

!         write(*,*) "this point is inOctal, setting havedescended to true"


        haveDescended = .TRUE. ! record that we have gone down the tree.
      
        ! if the point lies within the current octal, we identify the
        !   subcell
        subcell = whichSubcell(thisOctal,point)
!        write(*,*) "this point is in subcell ", subcell

        ! if a problem has been detected, this is where we complete the search
        IF (boundaryProblem) RETURN 
      
        ! if the subcell has a child, we look in the child for the point
        IF ( thisOctal%hasChild(subcell) ) THEN
!           write(*,*) "this subcell has a child"
                
          ! search the index to see where it is stored
          DO i = 1, thisOctal%nChildren, 1
            IF ( thisOctal%indexChild(i) == subcell ) THEN
                    
              thisOctal => thisOctal%child(i)

!              write(*,*) "calling findsubcelllocalprivate recursively"
              CALL findSubcellLocalPrivate(point,thisOctal,subcell,haveDescended,boundaryProblem)
              RETURN
              
            END IF
         END DO
          
      ELSE 
          RETURN
       END IF

      ELSE
        ! if the point is outside the current octal, we look in its
        !   parent octal
!         write(*,*) "this point is not in octal "
!         write(*,*) "point ",point
!         write(*,*) "xmin/max, zmin/max ",thisOctal%xmin, thisOctal%xmax, thisOctal%zmin,thisOctal%zmax
        ! first check that we are not outside the grid
        IF ( thisOctal%nDepth == 1 ) THEN
           if(.not. suppresswarnings) then
           write(*,*) "octal", thisOctal%ndepth
           write(*,*) "inoctal min", thisOctal%xMin,thisOctal%yMin,thisOctal%zMin
           write(*,*) "inoctal max", thisOctal%xMax,thisOctal%yMax,thisOctal%zMax

          PRINT *, 'Panic: In findSubcellLocalPrivate, point is outside the grid'
          write(*,*) point
          write(*,*) sqrt(point%x**2+point%y**2)
          write(*,*) atan2(point%y,point%x)*radtodeg
          write(*,*) " "
          write(*,*) "cylindrical hydro ",cylindricalhydro
          write(*,*) thisOctal%centre
          write(*,*) thisOctal%subcellSize
!          write(*,*) thisOctal%phi*radtodeg,thisOctal%dphi*radtodeg
          write(*,*) sqrt(thisOctal%centre%x**2+thisOctal%centre%y**2)
       endif
          if(.not. suppresswarnings) then
             r = -2.d0
             r = sqrt(r)
                STOP
          endif
          boundaryProblem = .TRUE.
          RETURN
       END IF
     
        ! if we have previously gone down the tree, and are now going back up, there
        !   must be a problem.
       
        IF (haveDescended) then
           boundaryProblem = .TRUE.
           phi =  atan2(point%y,point%x)
           phimin = thisOctal%phi - thisOctal%dphi/2.d0
           phimax = thisOctal%phi + thisOctal%dphi/2.d0
           PRINT *, 'Panic: In findSubcellLocalPrivate, have descended and are now going back up'
           write(*,*) "rank ",myrankglobal, thisOctal%mpithread(1:8)
           write(*,*) "split az ",thisOctal%splitAzimuthally
           write(*,*) point
           write(*,*) atan2(point%y,point%x)*radtodeg
           write(*,*) sqrt(point%x**2 + point%y**2)
           write(*,*) " "
           write(*,*) thisOctal%nDepth
           write(*,*) thisOctal%centre
           write(*,*) thisOctal%subcellSize
           write(*,*) thisOctal%phi*radtodeg,thisOctal%dphi*radtodeg
           write(*,*) sqrt(thisOctal%centre%x**2+thisOctal%centre%y**2)
           write(*,*) atan2(thisOctal%centre%y,thisOctal%centre%x)*radtodeg
           write(*,*) " x min/max, z min max ",thisOctal%xMin, thisOctal%xMax, thisOctal%zMin, thisOctal%zMax
           write(*,*) " r min/max ",thisOctal%r-thisOctal%subcellsize,thisOctal%r+thisOctal%subcellsize

           write(*,*) "x > xMin ",point%x > thisOctal%xMin
           write(*,*) "x < xMax ",point%x < thisOctal%xMax
           write(*,*) "z > zMin ",point%z > thisOctal%zMin
           write(*,*) "x < zMax ",point%z < thisOctal%zMax
           write(*,*) "parent x min/max, z min max ",thisOctal%parent%xMin, thisOctal%parent%xMax, thisOctal%parent%zMin, &
                thisOctal%parent%zMax
           write(*,*) "cen ",thisOctal%centre
           write(*,*) "size ",thisOctal%subcellsize
           write(*,*) "inoctal ",inoctal(thisOctal,point), thisOctal%phimin*radtodeg, &
                thisOctal%phimax*radtodeg,phi < phimin, phi > phimax, &
                phi < thisOctal%phimax, phi > thisOctal%phimin
!           rVec = subcellCentre(thisOctal,subcell)
!           write(*,*) rVec%x+thisOctal%subcellSize/2.
!           write(*,*) rVec%x-thisOctal%subcellSize/2.
!           write(*,*) rVec%y+thisOctal%subcellSize/2.
!           write(*,*) rVec%y-thisOctal%subcellSize/2.
!           write(*,*) rVec%z+thisOctal%subcellSize/2.
!           write(*,*) rVec%z-thisOctal%subcellSize/2.
           do ; enddo
!           STOP
           return
        endif
        
        IF ( thisOctal%nDepth /= 1 ) THEN
!           write(*,*) "ascending to octal above"
           thisOctal => thisOctal%parent
        ENDIF

        CALL findSubcellLocalPrivate(point,thisOctal,subcell,haveDescended,boundaryProblem)
       
     END IF
    
    END SUBROUTINE findSubcellLocalPrivate

  END SUBROUTINE findSubcellLocal

  SUBROUTINE findSubcellLocalLevel(point,thisOctal,subcell, nDepth, prob)
    ! finds the octal (and that octal's subcell) containing a point.
    !   starts searching from the current octal, and goes up and down the
    !   tree as needed to find the correct octal.
    use inputs_mod, only : hydrodynamics, cylindricalHydro, spherical

    IMPLICIT NONE
    integer :: nDepth
    TYPE(vector), INTENT(IN) :: point
    TYPE(vector) :: point_local
    TYPE(octal),POINTER    :: thisOctal
    INTEGER, INTENT(OUT)   :: subcell
    LOGICAL, INTENT(OUT),optional   :: prob
    
    LOGICAL                :: haveDescended    ! see comments below
    LOGICAL                :: boundaryProblem  ! see comments below
    
                             
    haveDescended = .FALSE.   ! if the 'point' lies very close to an 
    boundaryProblem = .FALSE. !   boundary, the program may go into 
                              !   a loop going up and down the tree.
                              ! we will keep track of the progress of
                              !   the search using these flags.
                             
    if (thisOctal%twoD) then
       if(cart2d) then
          point_local = point
          elseif (.not.cylindricalHydro) then
          point_local = projectToXZ(point)
       else
          point_local = point
       endif

!
!       if (.not.cylindricalhydro) then
!          point_local = projectToXZ(point)!
!       else!
!          point_local = point
!       endif
    else
       point_local = point
    endif
    if (thisOctal%oneD) then
       if (spherical) then
          point_local = VECTOR(modulus(point), 0.d0, 0.d0)
       else
          point_local = point
       endif
    endif


    CALL findSubcellLocalPrivateLevel(point_local,thisOctal,subcell,&
                                 haveDescended,boundaryProblem, nDepth)
    if (present(prob)) then
      prob = boundaryProblem
    else
      if (boundaryProblem) then
        stop 1
      endif
    endif
                                 
  CONTAINS

    RECURSIVE SUBROUTINE findSubcellLocalPrivateLevel(point,thisOctal,subcell,&
                                                 haveDescended,boundaryProblem, nDepth)
      TYPE(vector), INTENT(IN) :: point
      TYPE(octal),POINTER    :: thisOctal
      INTEGER, INTENT(OUT)   :: subcell
      integer :: nDepth
      LOGICAL, INTENT(INOUT) :: haveDescended
      LOGICAL, INTENT(INOUT) :: boundaryProblem
!      type(vector) :: rVec
      INTEGER :: i
      
      IF ( inOctal(thisOctal,point,alreadyRotated=.true.) ) THEN

        haveDescended = .TRUE. ! record that we have gone down the tree.
      
        ! if the point lies within the current octal, we identify the
        !   subcell
        subcell = whichSubcell(thisOctal,point)

        ! if a problem has been detected, this is where we complete the search
        IF (boundaryProblem) RETURN 
      
        ! if the subcell has a child, we look in the child for the point
        IF ( thisOctal%hasChild(subcell).and.(thisOctal%ndepth < nDepth) ) THEN
                
          ! search the index to see where it is stored
          DO i = 1, thisOctal%nChildren, 1
            IF ( thisOctal%indexChild(i) == subcell ) THEN
                    
              thisOctal => thisOctal%child(i)


              CALL findSubcellLocalPrivateLevel(point,thisOctal,subcell,haveDescended,boundaryProblem,nDepth)
              RETURN
              
            END IF
          END DO
          
        ELSE 
          RETURN
        END IF

      ELSE
        ! if the point is outside the current octal, we look in its
        !   parent octal

        ! first check that we are not outside the grid
        IF ( thisOctal%nDepth == 1 ) THEN
          PRINT *, 'Panic: In findSubcellLocal, point is outside the grid'
          write(*,*) point
          write(*,*) sqrt(point%x**2+point%y**2)
          write(*,*) atan2(point%y,point%x)*radtodeg
          write(*,*) " "
          write(*,*) thisOctal%centre
          write(*,*) thisOctal%subcellSize
          write(*,*) thisOctal%phi*radtodeg,thisOctal%dphi*radtodeg
          write(*,*) sqrt(thisOctal%centre%x**2+thisOctal%centre%y**2)
           DO ; END DO
          STOP
          boundaryProblem = .TRUE.
          RETURN
        END IF
     
        ! if we have previously gone down the tree, and are now going back up, there
        !   must be a problem.
        IF (haveDescended) then
           boundaryProblem = .TRUE.
           PRINT *, 'Panic: In findSubcellLocal, have descended and are now going back up'
           write(*,*) point
           write(*,*) atan2(point%y,point%x)*radtodeg
           write(*,*) sqrt(point%x**2 + point%y**2)
           write(*,*) " "
           write(*,*) thisOctal%nDepth
           write(*,*) thisOctal%centre
           write(*,*) thisOctal%subcellSize
           write(*,*) thisOctal%phi*radtodeg,thisOctal%dphi*radtodeg
           write(*,*) sqrt(thisOctal%centre%x**2+thisOctal%centre%y**2)
           
!           rVec = subcellCentre(thisOctal,subcell)
!           write(*,*) rVec%x+thisOctal%subcellSize/2.
!           write(*,*) rVec%x-thisOctal%subcellSize/2.
!           write(*,*) rVec%y+thisOctal%subcellSize/2.
!           write(*,*) rVec%y-thisOctal%subcellSize/2.
!           write(*,*) rVec%z+thisOctal%subcellSize/2.
!           write(*,*) rVec%z-thisOctal%subcellSize/2.
!           do ; enddo
!           STOP
           return
        endif
        
        IF ( thisOctal%nDepth /= 1 ) THEN
           thisOctal => thisOctal%parent
        ENDIF

        CALL findSubcellLocalPrivateLevel(point,thisOctal,subcell,haveDescended,boundaryProblem,nDepth)
       
      END IF    
    
    END SUBROUTINE findSubcellLocalPrivateLevel

  END SUBROUTINE findSubcellLocalLevel

  FUNCTION whichSubcell(thisOctal,point) RESULT (subcell)
    ! POINT --> MUST be pre-rotated for 2d case!!!!!!!!!!!!!!
    !
    ! returns the identification number (1-8) of the subcell of the 
    ! current octal which contains a given point
    ! NB this does NOT check that the point lies within the bounds of the octal!



    IMPLICIT NONE

    TYPE(octal), INTENT(IN)       :: thisOctal
    TYPE(vector), INTENT(IN) :: point
    INTEGER                       :: subcell
    real(double) :: r, phi

    if (thisOctal%oneD) then
       if (point%x <= thisOctal%centre%x) then
          subcell = 1
       else
          subcell = 2
       endif
       goto 666
    endif


    if (thisOctal%threed) then ! threed case 

       if (.not.thisOctal%cylindrical) then ! cartesian case

          IF ( point%x < thisOctal%centre%x ) THEN
             IF ( point%y < thisOctal%centre%y ) THEN
                IF ( point%z < thisOctal%centre%z ) THEN
                   subcell = 1
                ELSE 
                   subcell = 5
                   
                ENDIF
             ELSE 
                IF (point%z < thisOctal%centre%z) THEN
                   subcell = 3
                ELSE 
                   subcell = 7
                ENDIF
             END IF
          ELSE
             IF (point%y < thisOctal%centre%y) THEN
                IF (point%z < thisOctal%centre%z) THEN
                   subcell = 2
                ELSE 
                   subcell = 6
                ENDIF
             ELSE 
                IF (point%z < thisOctal%centre%z) THEN
                   subcell = 4
                ELSE 
                   subcell = 8
                ENDIF
             END IF
          ENDIF ! cartesian case

       else ! cylindrical case

          r = sqrt(point%x**2+point%y**2)
          phi = atan2(point%y, point%x)
          if (phi < 0.d0) phi = phi + twoPi

          if (thisOctal%splitAzimuthally) then ! azimuthal split case

             if (phi <= thisOctal%phi) then ! small y
                IF ( r <= thisOctal%r) THEN ! small x
                   IF ( point%z <= thisOctal%centre%z ) THEN
                      subcell = 1    ! small x, small y, small z
                   ELSE 
                      subcell = 5    ! small x, smally, big z
                   ENDIF
                ELSE
                   IF (point%z <= thisOctal%centre%z) THEN 
                      subcell = 2  ! big x, small y, small z
                   ELSE 
                      subcell = 6  ! big x, small y, big z
                   ENDIF
                END IF
             else
                IF ( r <= thisOctal%r ) THEN
                   IF ( point%z <= thisOctal%centre%z ) THEN
                      subcell = 3  ! small x, big y, small z
                   ELSE 
                      subcell = 7  ! small x, big y, big z
                   ENDIF
                ELSE
                   IF (point%z <= thisOctal%centre%z) THEN
                      subcell = 4  ! big x, big y, small z
                   ELSE 
                      subcell = 8 ! big x, big y, big z
                   ENDIF
                END IF
             endif
          else

             IF ( r <= thisOctal%r ) THEN
                IF ( point%z <= thisOctal%centre%z ) THEN
                   subcell = 1
                ELSE 
                   subcell = 3
                ENDIF
                ELSE
                   IF (point%z <= thisOctal%centre%z) THEN
                      subcell = 2
                   ELSE 
                      subcell = 4
                   ENDIF
                END IF
             endif ! azi case
          endif ! cylindrical

    else ! twoD case
       
       IF ( point%x <= thisOctal%centre%x ) THEN
          IF ( point%z <= thisOctal%centre%z ) THEN
             subcell = 1
          ELSE 
             subcell = 3
          ENDIF
       ELSE
          IF (point%z <= thisOctal%centre%z) THEN
             subcell = 2
          ELSE 
             subcell = 4
          ENDIF
       END IF
    endif

666 continue

  END FUNCTION whichSubcell    


  FUNCTION inOctal(thisOctal,point,alreadyRotated) 
    ! true if the point lies within the boundaries of the current octal
  
    use inputs_mod, only : hydrodynamics, cylindricalHydro, photoionPhysics, spherical
    use vector_mod, only : projectToXZ
    IMPLICIT NONE
    LOGICAL                       :: inOctal
    TYPE(octal), INTENT(IN)       :: thisOctal
    TYPE(vector), INTENT(IN) :: point
    TYPE(vector)             :: octVec2D
    real(double)                  :: r, phi, dphi, eps
    logical, intent(in), optional :: alreadyRotated
    logical :: doRotate

    doRotate = .true.
    if (PRESENT(alreadyRotated)) doRotate = .not.alreadyRotated

    if (thisOctal%threeD) then
       
       if (.not.thisOctal%cylindrical) then
          IF (point%x < thisOctal%xMin) THEN ; inOctal = .FALSE. ; goto 101
          ELSEIF (point%x > thisOctal%xMax) THEN ; inOctal = .FALSE.; goto 101
          ELSEIF (point%y < thisOctal%yMin) THEN ; inOctal = .FALSE.; goto 101
          ELSEIF (point%y > thisOctal%yMax) THEN ; inOctal = .FALSE.; goto 101
          ELSEIF (point%z < thisOctal%zMin) THEN ; inOctal = .FALSE.; goto 101
          ELSEIF (point%z > thisOctal%zMax) THEN ; inOctal = .FALSE.; goto 101
          ELSE  
             inOctal = .TRUE.
          ENDIF
101       continue
       else
          phi = atan2(point%y,point%x)
          if (phi < 0.d0) phi = phi + twoPi
          dphi = abs(phi - thisOctal%phi)
          r = sqrt(point%x**2 + point%y**2)
          eps =  0.d0 ! sqrt(epsilon(thisOctal%subcellsize))
          IF     (r < thisOctal%r - thisOctal%subcellSize - eps) THEN ; inOctal = .FALSE. 
          ELSEIF (r > thisOctal%r + thisOctal%subcellSize + eps) THEN ; inOctal = .FALSE.
!          ELSEIF (dphi > thisOctal%dphi/2.d0) THEN ; inOctal = .FALSE.
          ELSEIF (phi > thisOctal%phimax) THEN ; inOctal = .FALSE.
          ELSEIF (phi < thisOctal%phimin) THEN ; inOctal = .FALSE.
          ELSEIF (point%z < thisOctal%zMin) THEN ; inOctal = .FALSE.
          ELSEIF (point%z > thisOctal%zMax) THEN ; inOctal = .FALSE.
          ELSE  
             inOctal = .TRUE.
          ENDIF
       endif
    else ! twoD case
       if(cart2d) then
          octVec2d = point
       else if (.not.cylindricalHydro) then
          if (doRotate) then
             octvec2d = projectToXZ(point)
          else
             octVec2d = point
          end if
       else
          octvec2d = point
       endif


!       if (.not.cylindricalHydro) then
!          if (doRotate) then
!             octVec2D = projectToXZ(point)
!          else
!             octVec2D = point   
!          endif
!       else
!          octVec2D = point
!       endif
       if(hydrodynamics .and. (.not. cylindricalHydro).and.(.not.photoionPhysics)) then
          octVec2D = point
       end if
       IF (octVec2D%x < thisOctal%xMin) THEN ; inOctal = .FALSE. 
       ELSEIF (octVec2D%x > thisOctal%xMax) THEN ; inOctal = .FALSE.
       ELSEIF (octVec2D%z < thisOctal%zMin) THEN ; inOctal = .FALSE.
       ELSEIF (octVec2D%z > thisOctal%zMax) THEN ; inOctal = .FALSE.
       ELSE  
          inOctal = .TRUE.
       ENDIF
    endif

    if (thisOctal%oneD) then
       if (spherical) then
          r = modulus(point)
          if ( r < thisOctal%centre%x  - thisOctal%subcellSize) then ; inoctal = .false.
          else if (r > thisOctal%centre%x + thisOctal%subcellSize) then; inOctal = .false.
          else
             inOctal = .true.
          endif
       else
          if ( point%x < thisOctal%centre%x  - thisOctal%subcellSize) then ; inoctal = .false.
          else if (point%x > thisOctal%centre%x + thisOctal%subcellSize) then; inOctal = .false.
          else
             inOctal = .true.
          endif
       endif
       goto 666
    endif


!    if(thisOctal%threeD) then
!       if(isnan(point%x) .or. isnan(point%y) .or. isnan(point%z)) then
!          inOctal = .false.
!       end if
!    else if (thisOCtal%twoD) then
!       if(isnan(point%x)  .or. isnan(point%z)) then
!          inOctal = .false.
!       end if
!    else
!       if(isnan(point%x)) then
!          inOctal = .false.
!       end if
!    end if


666 continue
  END FUNCTION inOctal

  FUNCTION inSubcell(thisOctal,thisSubcell,point) 
    ! true if the point lies within the boundaries of the current octal
  
    IMPLICIT NONE
 
    LOGICAL                       :: inSubcell
    TYPE(octal), INTENT(IN)       :: thisOctal
    INTEGER, INTENT(IN)           :: thisSubcell
    TYPE(vector), INTENT(IN) :: point

    IF (inOctal(thisOctal,point)) THEN
      inSubcell = whichSubcell(thisOctal,point) == thisSubcell
    ELSE
      inSubcell = .FALSE.
    END IF
  
  END FUNCTION inSubcell

  FUNCTION looseInOctal(thisOctal,point) 
    ! true if the point lies 'loosely' in the current octal
    ! ( a 10% margin of error is allowed )
    ! this is useful only for testing purposes!
  
    IMPLICIT NONE
 
    LOGICAL                       :: looseInOctal
    TYPE(octal), INTENT(IN)       :: thisOctal
    TYPE(vector), INTENT(IN) :: point

    if (thisOctal%threeD) then
       IF ((point%x <= thisOctal%centre%x - 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%x >= thisOctal%centre%x + 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%y <= thisOctal%centre%y - 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%y >= thisOctal%centre%y + 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%z <= thisOctal%centre%z - 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%z >= thisOctal%centre%z + 1.1_oc * thisOctal%subcellSize )) THEN
          looseInOctal = .FALSE.
       ELSE  
          looseInOctal = .TRUE.
       ENDIF
    else
       IF ((point%x <= thisOctal%centre%x - 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%x >= thisOctal%centre%x + 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%z <= thisOctal%centre%z - 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%z >= thisOctal%centre%z + 1.1_oc * thisOctal%subcellSize )) THEN
          looseInOctal = .FALSE.
       ELSE  
          looseInOctal = .TRUE.
       ENDIF
    endif
  END FUNCTION looseInOctal


  subroutine distanceToCellBoundary(grid, posVec, direction, tVal, sOctal, sSubcell)
    use inputs_mod, only : spherical
    use octal_mod, only: returndPhi

    implicit none
    type(GRIDTYPE), intent(in)    :: grid
    type(VECTOR), intent(in) :: posVec
    type(VECTOR), intent(in) :: direction
    type(OCTAL), pointer, optional :: sOctal
    integer, optional :: sSubcell
    real(oct), intent(out) :: tval
    !
    type(VECTOR) :: rDirection
    type(OCTAL),pointer :: thisOctal
    real(double) :: distTor1, distTor2, theta, mu
    real(double) :: distToRboundary, compz,currentZ
    real(double) :: phi, distToZboundary, ang1, ang2
    type(VECTOR) :: subcen, point, xHat, rVec, rplane, rnorm, xVec
    integer :: subcell
    real(double) :: distToSide1, distToSide2, distToSide
    real(double) ::  compx,disttoxBoundary, halfCellSize, d2, fac
    real(oct) :: t(6),denom(6), r, r1, r2, d, cosmu,x1,x2, halfSubCellsize
    real(double) :: a, b, c 
    logical :: ok, thisOk(6)
!    integer :: jarray(6)

    type(VECTOR) :: normdiff

    point = posVec

    if (PRESENT(sOctal)) then
       if (PRESENT(sSubcell)) then
          subcell = sSubcell
          thisOctal => sOctal
       else
!          call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid, startOctal=sOctal)
          thisOctal => sOctal
          call findSubcellLocal(point, thisOctal, subcell)
       endif
    else
!       call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
       call findSubcellTD(point, grid%octreeRoot, thisOctal, subcell)
    endif
    subcen =  subcellCentre(thisOctal,subcell)

    if (thisOctal%oneD.and.spherical) then

       distToR1 = 1.d30
       distToR2 = 1.d30

       rVec = posVec
       call normalize(rVec)
       cosmu = ((-1.d0)*direction).dot.rVec
       d = modulus(posVec)

       ! distance to outer radius

       r2 = subcen%x + thisOctal%subcellSize/2.d0
       call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
       distToR2 = max(x1,x2)
!             write(*,*) "r2",x1,x2,disttor2

       !   inner radius

       r1 = subcen%x - thisOctal%subcellSize/2.d0
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
       goto 666
    endif

    if (thisOctal%oneD.and.(.not.spherical)) then
       if (direction%x > 0.d0) then
          tVal = (subcen%x + thisOctal%subcellSize/2.d0 - posVec%x)/direction%x
       else
          tVal = (posVec%x - (subcen%x - thisOctal%subcellSize/2.d0))/abs(direction%x)
       endif
       goto 666
    endif


       if (thisOctal%twoD .and. cart2d) then
          ok = .true.

          halfSubCellsize = thisOctal%subcellsize * 0.5d0
          denom = 0.d0

          if(direction%x .ne. 0.d0) then            
             denom(1) = 1.d0 / direction%x
          else
             denom(1) = 0.d0
          endif
          denom(4) = -denom(1)

!          if(direction%y .ne. 0.d0) then            
!             denom(2) = 1.d0 / direction%y
!          else
!             denom(2) = 0.d0
!          endif
!          denom(5) = -denom(2)
          
          if(direction%z .ne. 0.d0) then            
             denom(3) = 1.d0 / direction%z
          else
             denom(3) = 0.d0
          endif
          denom(6) = -denom(3)

          normdiff = subcen - posvec
          t = 0.d0
          thisOK = .false.
          ok = .false.
          t(1) =  (normdiff%x + halfsubcellsize) * denom(1)
          if (t(1) > 0.d0) then
             thisOK(1) = .true.
             ok = .true.
          endif
          t(2) =  (normdiff%y + halfsubcellsize) * denom(2)
          if (t(2) > 0.d0) then
             thisOK(2) = .true.
             ok = .true.
         endif
          t(3) =  (normdiff%z + halfsubcellsize) * denom(3)
          if (t(3) > 0.d0) then
             thisOK(3) = .true.
             ok = .true.
          endif
          t(4) =  (normdiff%x - halfsubcellsize) * denom(1)
          if (t(4) > 0.d0) then
             thisOK(4) = .true.
             ok = .true.
          endif
          t(5) =  (normdiff%y - halfsubcellsize) * denom(2)
          if (t(5) > 0.d0) then
             thisOK(5) = .true.
             ok = .true.
          endif
          t(6) =  (normdiff%z - halfsubcellsize) * denom(3)
          if (t(6) > 0.d0) then
             thisOK(6) = .true.
             ok = .true.
          endif

!          where(t > 0.d0)
!             jarray = 1
!             thisOk = .true.
!          elsewhere
!             jarray = 0
!             thisOk = .false.
!          end where
!          
!          j = sum(jarray)

!          if (j .eq. 0) ok = .false.
 
          if (.not.ok) then
             write(*,*) "Error: j=0 (no intersection???) in amr_mod::distanceToCellBoundary.2d "
             write(*,*) direction%x,direction%y,direction%z
             write(*,*) t(1:6)
             write(*,*) "denom: ", denom(1:6)
             write(*,*) "subcen", subcen
             write(*,*) "posvec", posvec
             write(*,*) "t", t
             write(*,*) "normdiff%x: ", normdiff%x
             write(*,*) "normdiff%y: ", normdiff%y
             write(*,*) "normdiff%z: ", normdiff%z
             write(*,*) "halfsubcellsize ", halfsubcellsize
             call torus_abort
          endif
          
          tval = minval(t, mask=thisOk)
          
! Commented out by Dave Acreman, October 2008
! tval == 0 is handled at the end of this subroutine    
!          if (tval == 0.) then
!             write(*,*) posVec
!             write(*,*) direction%x,direction%y,direction%z
!             write(*,*) t(1:6)
!             call torus_abort("tval==0 in distanceToCellBoundary")
!          endif

          !if (tval > sqrt(3.)*thisOctal%subcellsize) then
             !     write(*,*) "tval too big",tval/(sqrt(3.)*thisOctal%subcellSize)
             !     write(*,*) "direction",direction
             !     write(*,*) t(1:6)
             !     write(*,*) denom(1:6)
          !endif
          goto 666
       end if


    if (thisOctal%threed) then

       if (.not.thisOctal%cylindrical) then
          ok = .true.

          halfSubCellsize = thisOctal%subcellsize * 0.5d0

          if(direction%x .ne. 0.d0) then            
             denom(1) = 1.d0 / direction%x
          else
             denom(1) = 0.d0
          endif
          denom(4) = -denom(1)

          if(direction%y .ne. 0.d0) then            
             denom(2) = 1.d0 / direction%y
          else
             denom(2) = 0.d0
          endif
          denom(5) = -denom(2)
          
          if(direction%z .ne. 0.d0) then            
             denom(3) = 1.d0 / direction%z
          else
             denom(3) = 0.d0
          endif
          denom(6) = -denom(3)

          normdiff = subcen - posvec
          
          thisOK = .false.
          ok = .false.
          t(1) =  (normdiff%x + halfsubcellsize) * denom(1)
          if (t(1) > 0.d0) then
             thisOK(1) = .true.
             ok = .true.
          endif
          t(2) =  (normdiff%y + halfsubcellsize) * denom(2)
          if (t(2) > 0.d0) then
             thisOK(2) = .true.
             ok = .true.
          endif
          t(3) =  (normdiff%z + halfsubcellsize) * denom(3)
          if (t(3) > 0.d0) then
             thisOK(3) = .true.
             ok = .true.
          endif
          t(4) =  (normdiff%x - halfsubcellsize) * denom(1)
          if (t(4) > 0.d0) then
             thisOK(4) = .true.
             ok = .true.
          endif
          t(5) =  (normdiff%y - halfsubcellsize) * denom(2)
          if (t(5) > 0.d0) then
             thisOK(5) = .true.
             ok = .true.
          endif
          t(6) =  (normdiff%z - halfsubcellsize) * denom(3)
          if (t(6) > 0.d0) then
             thisOK(6) = .true.
             ok = .true.
          endif

!          where(t > 0.d0)
!             jarray = 1
!             thisOk = .true.
!          elsewhere
!             jarray = 0
!             thisOk = .false.
!          end where
!          
!          j = sum(jarray)

!          if (j .eq. 0) ok = .false.
 
          if (.not.ok) then
             write(*,*) "Error: j=0 (no intersection???) in amr_mod::distanceToCellBoundary. "
             write(*,*) direction%x,direction%y,direction%z
             write(*,*) t(1:6)
             write(*,*) "denom: ", denom(1:6)
             write(*,*) "subcen", subcen
             write(*,*) "posvec", posvec
             write(*,*) "t", t
             write(*,*) "normdiff%x: ", normdiff%x
             write(*,*) "normdiff%y: ", normdiff%y
             write(*,*) "normdiff%z: ", normdiff%z
             write(*,*) "halfsubcellsize ", halfsubcellsize
             call torus_abort
          endif
          
          tval = minval(t, mask=thisOk)
          
! Commented out by Dave Acreman, October 2008
! tval == 0 is handled at the end of this subroutine    
!          if (tval == 0.) then
!             write(*,*) posVec
!             write(*,*) direction%x,direction%y,direction%z
!             write(*,*) t(1:6)
!             call torus_abort("tval==0 in distanceToCellBoundary")
!          endif

          !if (tval > sqrt(3.)*thisOctal%subcellsize) then
             !     write(*,*) "tval too big",tval/(sqrt(3.)*thisOctal%subcellSize)
             !     write(*,*) "direction",direction
             !     write(*,*) t(1:6)
             !     write(*,*) denom(1:6)
          !endif

       else

          ! now look at the cylindrical case

          halfCellSize = thisOctal%subcellSize/2.d0
          rVec = subcellCentre(thisOctal,subcell)
          r = sqrt(rVec%x**2 + rVec%y**2)
          r1 = r - halfCellSize
          r2 = r + halfCellSize

          distToR1 = 1.d30
          distToR2 = 1.d30
          d = sqrt(point%x**2+point%y**2)
          xHat = VECTOR(point%x, point%y,0.d0)
          if (modulus(xhat)/=0.d0)    call normalize(xHat)
          rDirection = VECTOR(direction%x, direction%y,0.d0)
          compX = modulus(rDirection)
          if (modulus(rDirection) /= 0.d0) call normalize(rDirection)
          if (compX /= 0.d0) then
             cosmu =((-1.d0)*xHat).dot.rdirection
             call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
             if (.not.ok) then
                write(*,*) "Quad solver failed in distanceToCellBoundary I",d,cosmu,r2
                write(*,*) "xhat",xhat
                write(*,*) "dir",direction
                write(*,*) "point",point
                do;enddo
                endif
                if ((x1.lt.0.d0).and.(x2.lt.0.d0)) then
                   write(*,*) "x1, x2 ",x1,x2
                   write(*,*) "rdirection ",rdirection
                   write(*,*) "xhat ",xhat
                   write(*,*) "compx ",compx
                   write(*,*) "cosmu ",cosmu
                endif
                distTor2 = max(x1,x2)/compX

                if ((d .ne. 0.).and.(r1 > 0.1d0*grid%halfSmallestSubcell)) then
                   theta = asin(max(-1.d0,min(1.d0,r1 / d)))
                   cosmu = ((-1.d0)*xHat).dot.rdirection
                   mu = acos(max(-1.d0,min(1.d0,cosmu)))
                   distTor1 = 1.e30
                   if (mu  < theta ) then
                      call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
                      if (.not.ok) then
                         write(*,*) "Quad solver failed in distanceToCellBoundary II",d,cosmu,r1,x1,x2
                         write(*,*) "coeff b",-2.d0*d*cosmu, "coeff c", d**2-r2**2
                         write(*,*) "direction ",direction
                         write(*,*) "mu ",mu*radtodeg, "theta ",theta*radtodeg
                         x1 = thisoctal%subcellSize/2.d0
                         x2 = 0.d0
                      endif
                      distTor1 = min(x1,x2)/compX
                   endif
                else
                   distTor1 = 1.d30
                end if
             endif
             distToRboundary = min(distTor1, distTor2)

             ! now do the upper and lower (z axis) surfaces

             compZ = zHat.dot.direction
             currentZ = point%z

             if (compZ /= 0.d0 ) then
                if (compZ > 0.d0) then
                   distToZboundary = (subcen%z + halfCellSize - currentZ ) / compZ
                else
                   distToZboundary = abs((subcen%z - halfCellSize - currentZ ) / compZ)
                endif
             else
                disttoZboundary = 1.e30
             endif

             ! ok now we have to tackle the two angled sides...
         thisOk = .true.
         
! Commented out by Dave Acreman, October 2008
! tval == 0 is handled at the end of this subroutine         
!         if (tval == 0.) then
!            write(*,*) posVec
!            write(*,*) direction%x,direction%y,direction%z
!            write(*,*) t(1:6)
!            call torus_abort("tval==0 in distanceToCellBoundary")
!         endif
         
         !if (tval > sqrt(3.)*thisOctal%subcellsize) then
            !     write(*,*) "tval too big",tval/(sqrt(3.)*thisOctal%subcellSize)
            !     write(*,*) "direction",direction
            !     write(*,*) t(1:6)
            !     write(*,*) denom(1:6)
         !endif

             rVec = subcellCentre(thisOctal, subcell)
             phi = atan2(rVec%y, rVec%x)
             if (phi < 0.d0) phi = phi + twoPi

             ang1 = phi - returndPhi(thisOctal)
             rPlane = VECTOR(cos(ang1),sin(ang1),0.d0)
             rnorm = rplane .cross. VECTOR(0.d0, 0.d0, 1.d0)
             call normalize(rnorm)
             distToSide1 = 1.d30
             if ((rnorm .dot. direction) /= 0.d0) then
                distToSide1 = (rnorm.dot.(rPlane-posVec))/(rnorm.dot.direction)
                if (distToSide1 < 0.d0) distToSide1 = 1.d30
             endif

             ang2 = phi + returndPhi(thisOctal)
             rPlane = VECTOR(cos(ang2),sin(ang2),0.d0)
             rnorm = rplane .cross.  VECTOR(0.d0, 0.d0, 1.d0)
             call normalize(rnorm)
             distToSide2 = 1.d30
             if ((rnorm .dot. direction) /= 0.d0) then
                distToSide2 = (rnorm.dot.(rPlane-posVec))/(rnorm.dot.direction)
                if (distToSide2 < 0.d0) distToSide2 = 1.d30
             endif

             distToSide = min(distToSide1, distToSide2)


             tVal = min(distToZboundary, distToRboundary, distToSide)
             if (tVal > 1.e29) then
                write(*,*) "Cylindrical ",tval
                write(*,*) tVal,compX,compZ, distToZboundary,disttorboundary, disttoside
                write(*,*) "subcen",subcen
                write(*,*) "z", currentZ
             endif
             if (tval < 0.) then
                write(*,*) "Cylindrical ",tval
                write(*,*) tVal,distToZboundary,disttorboundary, disttoside
                write(*,*) "subcen",subcen
                write(*,*) "z", currentZ
                write(*,*) "disttor1, disttor2 ",disttor1,disttor2
             endif

          endif

       else ! two-d grid case below

          halfCellSize = thisOctal%subcellSize/2.d0
          r1 = max(0.d0, subcen%x - halfCellSize)
          r2 = subcen%x + halfCellSize

          distToR1 = 1.d30
          distToR2 = 1.d30
          d2 = point%x**2+point%y**2
          d = sqrt(d2)
          xVec = VECTOR(point%x, point%y,0.d0)
          xHat = xVec
          call normalize(xHat)
          rDirection = VECTOR(direction%x, direction%y,0.d0)
          compX = modulus(rDirection)
          if (compX /= 0.d0) call normalize(rDirection)

          if (compX /= 0.d0) then
             cosmu =((-1.d0)*xHat).dot.rdirection
             call solveQuadDble(1.d0, -2.d0*d*cosmu, d2-r2**2, x1, x2, ok)
             if (.not.ok) then
                write(*,*) "Quad solver failed in distanceToCellBoundary I",d,cosmu,r2
                write(*,*) "xhat",xhat
                write(*,*) "dir",direction
                write(*,*) "point",point
                do
                enddo
                x1 = thisoctal%subcellSize/2.d0
                x2 = 0.d0
             endif
             distTor2 = max(x1,x2)/compX

             distTor1 = 1.e30

             if (cosmu > 0.d0) then
                a = 1.d0
                b = -2.d0*d*cosmu
                c = d2-r1**2
                fac = b*b-4.d0*a*c
                if (fac > 0.d0) then
                   call solveQuadDble(a, b, c, x1, x2, ok)
                   !               if(ok) then
                   !                  write(*,*) "All good",d,cosmu,r1,x1,x2
                   !                  write(*,*) "coeff b",-2.d0*d*cosmu, "coeff c", d**2-r2**2
                   !               endif
                   
                   if (.not.ok) then
                      write(*,*) "mu", mu, "theta", theta
                      write(*,*) "Quad solver failed in distanceToCellBoundary IIb",d,cosmu,r1,x1,x2
                      write(*,*) "coeff b",-2.d0*d*cosmu, "coeff c", d**2-r1**2
                      write(*,*) "xhat",xhat
                      write(*,*) "dir",direction
                      write(*,*) "point",point
                      
                      x1 = thisoctal%subcellSize/2.d0
                      x2 = 0.d0
                   endif
                   distTor1 = min(x1,x2)/compX
                endif
             endif
          endif
          distToXboundary = min(distTor1, distTor2)


          compZ = zHat.dot.direction
          currentZ = point%z

          if (compZ /= 0.d0 ) then
             if (compZ > 0.d0) then
                distToZboundary = (subcen%z + halfCellSize - currentZ ) / compZ
             else
                distToZboundary = abs((subcen%z - halfCellSize - currentZ ) / compZ)
             endif
          else
             disttoZboundary = 1.e30
          endif

          tVal = min(distToZboundary, distToXboundary)
          if (tVal > 1.e29) then
             write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
             write(*,*) "subcen",subcen
             write(*,*) "z", currentZ
             write(*,*) "TVAL", tval
             write(*,*) "direction", direction
             call torus_abort
          endif

!          if (tval < 0.) then
!             write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
!             write(*,*) "subcen",subcen
             !         write(*,*) "x,z",currentX,currentZ
!          endif
      
   endif

666    continue

!       tVal = max(tVal, 0.001d0*grid%halfSmallestSubcell) ! avoid sticking on a cell boundary


     end subroutine distanceToCellBoundary

  subroutine distanceToNearestWall(posVec, tVal, sOctal, sSubcell)
    use inputs_mod, only : amr1d, amr2d, amr3d, cylindrical
    implicit none
    type(VECTOR), intent(in) :: posVec
    type(OCTAL), pointer :: sOctal
    integer :: sSubcell
    real(double), intent(out) :: tval
    real(double) :: d, p, r
    type(VECTOR) :: cen
    
    if (amr3d.and.cylindrical) then
       write(*,*) "distanceToNearestWall not implemented for this geometry"
       stop
    endif
    d = sOctal%subcellSize/2.d0
    cen = subcellCentre(sOctal,ssubcell)

    tVal = 1.d30
    if (amr3d) then
       tVal = min((cen%x + d) - posVec%x, tVal)
       tVal = min((cen%y + d) - posVec%y, tVal)
       tVal = min((cen%z + d) - posVec%z, tVal)
       
       tVal = min(posVec%x - (cen%x - d), tVal)
       tVal = min(posVec%y - (cen%y - d), tVal)
       tVal = min(posVec%z - (cen%z - d), tVal)
    else if (amr2d) then
       r = sqrt(posVec%x**2 + posVec%y**2)
       tVal = min((cen%x + d) - r, tVal)
       tVal = min((cen%z + d) - posVec%z, tVal)
       
       tVal = min(r - (cen%x - d), tVal)
       tVal = min(posVec%z - (cen%z - d), tVal)
    else if (amr1d) then
       p = modulus(posVec)
       tVal = min((cen%x + d) - p, tval)
       tVal = min(p - (cen%x - d), tval)
    endif
  end subroutine distanceToNearestWall

  subroutine distanceToGridEdge(grid, posVec, direction, tVal)

   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(VECTOR), intent(in) :: posVec
   type(VECTOR), intent(in) :: direction
   real(oct), intent(out) :: tval
   !
   real(double) :: distTor1, distTor2, theta, mu
   real(double) :: distToRboundary, compz,currentZ
   real(double) :: distToZboundary !, ang1, ang2 , phi
   type(VECTOR) :: subcen, point, xHat, zHat !, rVec
   real(double) :: distToSide  !, distToSide1, distToSide2
   real(double) :: disttoxBoundary, subcellsize, halfsubcellsize
   real(oct) :: t(6),denom(6), r, r1, r2, d, cosmu,x1,x2
   integer :: i,j
   logical :: ok, thisOk(6)

   type(VECTOR) :: normdiff

   point = posVec

   subcen =  grid%octreeRoot%centre
   subcellsize = grid%octreeRoot%subcellSize * 2.d0
   halfsubcellsize = grid%octreeRoot%subcellSize

   if (grid%octreeRoot%threed) then

      if (.not.grid%octreeRoot%cylindrical) then
         ok = .true.
         
         if(direction%x .ne. 0.d0) then            
            denom(1) = 1.d0 / direction%x
         else
            denom(1) = 0.d0
         endif
         
         if(direction%y .ne. 0.d0) then            
            denom(2) = 1.d0 / direction%y
         else
            denom(2) = 0.d0
         endif
         
         if(direction%z .ne. 0.d0) then            
            denom(3) = 1.d0 / direction%z
         else
            denom(3) = 0.d0
         endif
                     
         normdiff = subcen - posvec

         t(1) =  (normdiff%x + halfsubcellsize) * denom(1)
         t(2) =  (normdiff%y + halfsubcellsize) * denom(2)
         t(3) =  (normdiff%z + halfsubcellsize) * denom(3)
         t(4) =  (normdiff%x - halfsubcellsize) * denom(1)
         t(5) =  (normdiff%y - halfsubcellsize) * denom(2)
         t(6) =  (normdiff%z - halfsubcellsize) * denom(3)

         thisOk = .true.
         
         do i = 1, 6
            
            if (denom(i) .eq. 0.0d0) then
               thisOk(i) = .false.
            endif

            if (t(i) < 0.) thisOk(i) = .false.
         enddo
                  
         j = 0
         do i = 1, 6
            if (thisOk(i)) j=j+1
         enddo
         
         if (j == 0) ok = .false.
         
         if (.not.ok) then
            write(*,*) "Error: j=0 (no intersection???) in amr_mod::distanceToGridEdge. "
            write(*,*) direction%x,direction%y,direction%z
            write(*,*) t(1:6)
            stop
         endif
         
         tval = minval(t, mask=thisOk)

         if (tval == 0.) then
            write(*,*) posVec
            write(*,*) direction%x,direction%y,direction%z
            write(*,*) t(1:6)
            stop
         endif
         
         if (tval > sqrt(3.)*subcellsize) then
            !     write(*,*) "tval too big",tval/(sqrt(3.)*thisOctal%subcellSize)
            !     write(*,*) "direction",direction
            !     write(*,*) t(1:6)
            !     write(*,*) denom(1:6)
         endif

      else

! now look at the cylindrical case

         ! first do the inside and outside curved surfaces
         r = sqrt(subcen%x**2 + subcen%y**2)
         r1 = r -subcellSize/2.d0
         r2 = r +subcellSize/2.d0
         d = sqrt(point%x**2+point%y**2)
         xHat = VECTOR(point%x, point%y,0.d0)
         call normalize(xHat)
      
         cosmu =((-1.d0)*xHat).dot.direction
         call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
         if (.not.ok) then
            write(*,*) "Quad solver failed in distanceToGridEdge"
            x1 = subcellSize/2.d0
            x2 = 0.d0
         endif
         distTor2 = max(x1,x2)
         
         theta = asin(max(-1.d0,min(1.d0,r1 / d)))
         cosmu = xHat.dot.direction
         mu = acos(max(-1.d0,min(1.d0,cosmu)))
         distTor1 = 1.e30
         if (mu  < theta ) then
            call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
            if (.not.ok) then
               write(*,*) "Quad solver failed in distanceToGridEdge"
               x1 = subcellSize/2.d0
               x2 = 0.d0
            endif
            distTor1 = max(x1,x2)
         endif
      
         distToRboundary = min(distTor1, distTor2)

         ! now do the upper and lower (z axis) surfaces
      
         zHat = VECTOR(0.d0, 0.d0, 1.d0)
         compZ = zHat.dot.direction
         currentZ = point%z
      
         if (compZ /= 0.d0 ) then
            if (compZ > 0.d0) then
               distToZboundary = (subcen%z + subcellsize/2.d0 - currentZ ) / compZ
            else
               distToZboundary = abs((subcen%z - subcellsize/2.d0 - currentZ ) / compZ)
            endif
         else
            disttoZboundary = 1.e30
         endif
      
        
!         ! ok now we have to tackle the two angled sides...
!
!         ! find posvec to surface centre
!
!         phi = atan2(posVec%y,posVec%x)
!         if (phi < 0.d0) phi = phi + twoPi
!
!         rVec = VECTOR(r, 0.d0, 0.d0)
!         if (grid%octreeRoot%splitAzimuthally) then
!            if (phi < grid%octreeRoot%phi) then
!               ang1 = grid%octreeRoot%phi - grid%octreeRoot%dPhi/2.d0
!               ang2 = grid%octreeRoot%phi
!            else
!               ang1 = grid%octreeRoot%phi
!               ang2 = grid%octreeRoot%phi + grid%octreeRoot%dPhi/2.d0
!            endif
!         else
!            ang1 = grid%octreeRoot%phi - grid%octreeRoot%dPhi/2.d0
!            ang2 = grid%octreeRoot%phi + grid%octreeRoot%dPhi/2.d0
!         endif
!
!         rVec = VECTOR(r, 0.d0, 0.d0)
!         rVec = rotateZ(rVec, -ang1)
!         thisnorm = rVec .cross. zHat
!         call normalize(thisnorm)
!         if ((thisnorm.dot.direction) /= 0.d0) then
!            distToSide1 = (thisnorm.dot.(rVec-posVec))/(thisnorm.dot.direction)
!            if (distToSide1 < 0.d0) distToSide1 = 1.d30
!         endif
!
!         rVec = VECTOR(r, 0.d0, 0.d0)
!         rVec = rotateZ(rVec, -ang2)
!         thisnorm = rVec .cross. zHat
!         call normalize(thisnorm)
!         if ((thisnorm.dot.direction) /= 0.d0) then
!            distToSide2 = (thisnorm.dot.(rVec-posVec))/(thisnorm.dot.direction)
!            if (distToSide2 < 0.d0) distToSide2 = 1.d30
!         endif
!
!         distToSide = min(distToSide1, distToside2)
         distToSide = 1.d30

         tVal = min(distToZboundary, distToRboundary, distToSide)

         write(*,*) disttoside,disttoZboundary,disttoRboundary
         if (tVal > 1.e29) then
            write(*,*) "Cylindrical"
            write(*,*) tVal,compZ, distToZboundary,disttorboundary, disttoside
            write(*,*) "subcen",subcen
            write(*,*) "z", currentZ
         endif
         if (tval < 0.) then
            write(*,*) "Cylindrical"
            write(*,*) tVal,distToZboundary,disttorboundary, disttoside
            write(*,*) "subcen",subcen
            write(*,*) "z", currentZ
         endif

      endif

   else ! two-d grid case below

      r1 = subcen%x - subcellSize/2.d0
      r2 = subcen%x + subcellSize/2.d0
      d = sqrt(point%x**2+point%y**2)
      xHat = VECTOR(point%x, point%y,0.d0)
      call normalize(xHat)
      
      cosmu =((-1.d0)*xHat).dot.direction
      call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
      if (.not.ok) then
         write(*,*) "Quad solver failed in distanceToGridEdge"
         x1 = subcellSize/2.d0
         x2 = 0.d0
      endif
      distTor2 = max(x1,x2)
      
      theta = asin(max(-1.d0,min(1.d0,r1 / d)))
      cosmu = xHat.dot.direction
      mu = acos(max(-1.d0,min(1.d0,cosmu)))
      distTor1 = 1.e30
      if (mu  < theta ) then
         call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
         if (.not.ok) then
            write(*,*) "Quad solver failed in distanceToGridEdge"
            x1 = subcellSize/2.d0
            x2 = 0.d0
         endif
         distTor1 = max(x1,x2)
      endif
      
      distToXboundary = min(distTor1, distTor2)
      
      
      zHat = VECTOR(0.d0, 0.d0, 1.d0)
      compZ = zHat.dot.direction
      currentZ = point%z
      
      if (compZ /= 0.d0 ) then
         if (compZ > 0.d0) then
            distToZboundary = (subcen%z + subcellsize/2.d0 - currentZ ) / compZ
         else
            distToZboundary = abs((subcen%z - subcellsize/2.d0 - currentZ ) / compZ)
         endif
      else
         disttoZboundary = 1.e30
      endif
      
      tVal = min(distToZboundary, distToXboundary)
      if (tVal > 1.e29) then
         write(*,*) tVal,compZ, distToZboundary,disttoxboundary
         write(*,*) "subcen",subcen
         write(*,*) "z",currentZ
      endif
      if (tval < 0.) then
         write(*,*) tVal,compZ, distToZboundary,disttoxboundary
         write(*,*) "subcen",subcen
         write(*,*) "z", currentZ
      endif
      
   endif

   tVal = max(tVal, 1.d-4*grid%halfSmallestSubcell) ! avoid sticking on a cell boundary

 end subroutine distanceToGridEdge


  type(VECTOR) function randomPositionInCell(thisOctal, subcell)
    use inputs_mod, only : spherical
    use octal_mod, only: returndPhi

    type(OCTAL) :: thisOctal
    integer :: subcell
    type(VECTOR) :: octalCentre
    real(double) :: r1, r2, r3, r
    real(double) :: xOctal, yOctal, zOctal
    real(double) :: ang, ang1, ang2, phi


    octalCentre = subcellCentre(thisOctal,subcell)
    
!!! we will just choose a random point within the subcell.
!!! this *should* be done in a better way.


    if (thisOctal%oneD) then
       call randomNumberGenerator(getDouble=r1)
       r1 = r1 - 0.5d0
       r1 = r1 * 0.9999
       r = r1 * thisOctal%subcellSize + octalCentre%x
       if (spherical) then
          randomPositionInCell = r * randomUnitVector()
       else
          randomPositionInCell = VECTOR(r, 0.d0, 0.d0)
       endif
       goto 666
    endif


    if (thisOctal%threed) then

       if (.not.thisOctal%cylindrical) then

          call randomNumberGenerator(getDouble=r1)
          r1 = r1 - 0.5  ! shift value mean value to zero
          r1 = r1 * 0.9999 ! to avoid any numerical accuracy problems
          xOctal = r1 * thisOctal%subcellSize + octalCentre%x
       
          call randomNumberGenerator(getDouble=r2)
          r2 = r2 - 0.5                                  
          r2 = r2 * 0.9999                                          
          yOctal = r2 * thisOctal%subcellSize + octalCentre%y
          
          
          call randomNumberGenerator(getDouble=r3)
          r3 = r3 - 0.5                                  
          r3 = r3 * 0.9999                                          
          zOctal = r3 * thisOctal%subcellSize + octalCentre%z
          
          randomPositionInCell = VECTOR(xOctal,yOctal,zOctal)
          
       else


          call randomNumberGenerator(getDouble=r1)
          r1 = r1 - 0.5  ! shift value mean value to zero
          r1 = r1 * 0.9999 ! to avoid any numerical accuracy problems
          xOctal = r1 * thisOctal%subcellSize + sqrt(octalCentre%x**2+octalCentre%y**2)


          call randomNumberGenerator(getDouble=r3)
          r3 = r3 - 0.5                                  
          r3 = r3 * 0.9999                                          
          zOctal = r3 * thisOctal%subcellSize + octalCentre%z

          randomPositionInCell = VECTOR(xOctal,0.,zOctal)
          
          call randomNumberGenerator(getDouble=r2)
          phi = atan2(octalCentre%y, octalCentre%x)
          if (phi < 0.d0) phi = phi + twoPi

          ang1 = phi - returndPhi(thisOctal)
          ang2 = phi + returndPhi(thisOctal)
          
!          if (thisOctal%splitAzimuthally) then
!             ang1 = thisOctal%phi - thisOctal%dphi/4.
!             ang2 = thisOctal%phi + thisOctal%dphi/4.
!          Else
!             ang1 = thisOctal%phi - thisOctal%dphi/2.
!             ang2 = thisOctal%phi + thisOctal%dphi/2.
!          endif
          ang = ang1 + r2 * (ang2 - ang1)
          randomPositionInCell = rotateZ(randomPositionInCell, -ang)

       endif
    else
       
       call randomNumberGenerator(getDouble=r1)
       r1 = r1 - 0.5  ! shift value mean value to zero
       r1 = r1 * 0.9999 ! to avoid any numerical accuracy problems
       xOctal = r1 * thisOctal%subcellSize + octalCentre%x
       
          
       call randomNumberGenerator(getDouble=r3)
       r3 = r3 - 0.5                                  
       r3 = r3 * 0.9999                                          
       zOctal = r3 * thisOctal%subcellSize + octalCentre%z
       
       randomPositionInCell = VECTOR(xOctal,0.d0,zOctal)

       if (thisOctal%twod) then

          call randomNumberGenerator(getDouble=ang)
          ang = ang * twoPi
          if(cart2d) then
             randomPositionInCell = rotateZ(randomPositionInCell, ang)
          end if
       endif
    endif
666 continue
  end function randomPositionInCell


subroutine returnVelocityVector2(grid, position, velocity)
!For use in molecular_mod.

integer :: subcell
type(octal), pointer :: thisOctal
type(gridtype) :: grid
type(vector) :: position 
type(vector) :: velocity
real(double) :: rho


if(inOctal(grid%octreeRoot, position)) then
   thisOctal => grid%octreeRoot
   call findsubcelllocal(position, thisoctal,subcell)

   rho = thisOctal%rho(subcell)
   velocity = VECTOR(thisOctal%rhou(subcell)/(rho*cspeed), &
        thisOctal%rhov(subcell)/(rho*cspeed), &
        thisOctal%rhow(subcell)/(rho*cspeed))
else
end if

end subroutine returnVelocityVector2


SUBROUTINE fillHydroDensityVelocityCorners(thisOctal, grid)

  type(gridtype) :: grid
  TYPE(octal), pointer :: thisOctal
  TYPE(octal), pointer :: probeOctal
  type(vector), allocatable :: rVecArray(:), probeArray(:)
  type(vector) :: position
  integer ::  j
!  real(double) :: radius, x, y, z, xmin, zmin, dx, dz
  integer, parameter :: maxpts = 8
  integer :: probeSubcell
  real(double) :: x1, x2, x3, y1, y2, y3, z1, z2, z3
  real(double) :: rhoPoints(maxpts)
  real(double) :: rhouPoints(maxpts)
  real(double) :: rhovPoints(maxpts)
  real(double) :: rhowPoints(maxpts)
  type(vector), allocatable :: sourcePoints(:)
  integer :: nPoints, f, k
  logical :: replica 

  if (thisOctal%threed) then
     
     allocate(rVecArray(27))
     allocate(probeArray(8))
     allocate(sourcePoints(8))
     
     probeArray(1) = VECTOR(1.d0, 1.d0, 1.d0)
     probeArray(2) = VECTOR(1.d0, 1.d0, -1.d0)
     probeArray(3) = VECTOR(1.d0, -1.d0, 1.d0)
     probeArray(4) = VECTOR(-1.d0, 1.d0, 1.d0)
     probeArray(5) = VECTOR(1.d0, -1.d0, -1.d0)
     probeArray(6) = VECTOR(-1.d0, -1.d0, 1.d0)
     probeArray(7) = VECTOR(-1.d0, 1.d0, -1.d0)
     probeArray(8) = VECTOR(-1.d0, -1.d0, -1.d0)
     
    
     x1 = thisOctal%centre%x - thisOctal%subcellSize
     x2 = thisOctal%centre%x
     x3 = thisOctal%centre%x + thisOctal%subcellSize
     
     y1 = thisOctal%centre%y - thisOctal%subcellSize
     y2 = thisOctal%centre%y
     y3 = thisOctal%centre%y + thisOctal%subcellSize
     
     z1 = thisOctal%centre%z - thisOctal%subcellSize
     z2 = thisOctal%centre%z
     z3 = thisOctal%centre%z + thisOctal%subcellSize
    
     rVecArray(1) = VECTOR(x1,y1,z1)
     rVecArray(2) = VECTOR(x2,y1,z1)
     rVecArray(3) = VECTOR(x3,y1,z1)
     rVecArray(4) = VECTOR(x1,y2,z1)
     rVecArray(5) = VECTOR(x2,y2,z1)
     rVecArray(6) = VECTOR(x3,y2,z1)
     rVecArray(7) = VECTOR(x1,y3,z1)
     rVecArray(8) = VECTOR(x2,y3,z1)
     rVecArray(9) = VECTOR(x3,y3,z1)
     rVecArray(10) = VECTOR(x1,y1,z2)
     rVecArray(11) = VECTOR(x2,y1,z2)
     rVecArray(12) = VECTOR(x3,y1,z2)
     rVecArray(13) = VECTOR(x1,y2,z2)
     rVecArray(14) = VECTOR(x2,y2,z2)
     rVecArray(15) = VECTOR(x3,y2,z2)
     rVecArray(16) = VECTOR(x1,y3,z2)
     rVecArray(17) = VECTOR(x2,y3,z2)
     rVecArray(18) = VECTOR(x3,y3,z2)
     rVecArray(19) = VECTOR(x1,y1,z3)
     rVecArray(20) = VECTOR(x2,y1,z3)
     rVecArray(21) = VECTOR(x3,y1,z3)
     rVecArray(22) = VECTOR(x1,y2,z3)
     rVecArray(23) = VECTOR(x2,y2,z3)
     rVecArray(24) = VECTOR(x3,y2,z3)
     rVecArray(25) = VECTOR(x1,y3,z3)
     rVecArray(26) = VECTOR(x2,y3,z3)
     rVecArray(27) = VECTOR(x3,y3,z3)

     do j = 1, 27
        nPoints = 0
        rhoPoints = 0.d0
        rhouPoints = 0.d0
        rhovPoints = 0.d0
        rhowPoints = 0.d0
        do k = 1, 8
           replica = .false.
           position = rVecArray(j) + rmult(0.01d0*grid%halfsmallestsubcell,probeArray(k))
           probeOctal => thisOctal
           if(inOctal(grid%octreeRoot, position)) then
              call findSubcellLocal(position, probeOctal, probeSubcell)
              position = subcellCentre(probeOctal, probeSubcell)
              do f = 1, k
                 if(vectorEquivalence(position,sourcePoints(f))) then
                    replica = .true.
                 end if
              enddo
              if(.not. replica) then
                 sourcePoints(k) = position
                 rhoPoints(k) = probeOctal%rho(probeSubcell)
                 rhouPoints(k) = probeOctal%rhou(probeSubcell)
                 rhovPoints(k) = probeOctal%rhov(probeSubcell)
                 rhowPoints(k) = probeOctal%rhow(probeSubcell)                            
                 nPoints = nPoints + 1
              end if
!           else
!              rhoPoints(k) = 0.d0
!              rhouPoints(k) = 0.d0
!              rhovPoints(k) = 0.d0
!              rhowPoints(k) = 0.d0
!              nPoints = nPoints + 1
           end if

        enddo
        
        if(nPoints /= 0) then
           thisOctal%cornerRho(j) = SUM(rhoPoints(1:nPoints))/dble(nPoints)
!           print *, "thisOctal%cornerRho(j)", thisOctal%cornerRho(j)
           thisOctal%cornerVelocity(j)%x = SUM(rhouPoints(1:nPoints))/(dble(nPoints)*thisOctal%cornerRho(j)*cspeed)
           thisOctal%cornerVelocity(j)%y = SUM(rhovPoints(1:nPoints))/(dble(nPoints)*thisOctal%cornerRho(j)*cspeed)
           thisOctal%cornerVelocity(j)%z = SUM(rhowPoints(1:nPoints))/(dble(nPoints)*thisOctal%cornerRho(j)*cspeed)
        else
           thisOctal%cornerRho(j) = 0.d0
           thisOctal%cornerVelocity(j)%x = 0.d0
           thisOctal%cornerVelocity(j)%y = 0.d0
           thisOctal%cornerVelocity(j)%z = 0.d0

!           print *, "zero cells surrounding points ", position
!           call torus_abort("aborting...")
        end if

!        print *, "rhoPoints", rhoPoints
!        print *, "vel ", thisOctal%cornervelocity(j)
!        print *, "dble(nPoints)", dble(nPoints)
!        print *, "thisOctal%cornerRho(j)", thisOctal%cornerRho(j)
     end do
!     stop
  else 
     call torus_abort("Corner velocities only available in 3D calculations at present")

  end if
     
END SUBROUTINE fillHydroDensityVelocityCorners

end module amr_utils_mod
