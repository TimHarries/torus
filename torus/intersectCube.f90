
  subroutine intersectCube(grid, posVec, direction, intersection, ok)
   use vector_mod
   use gridtype_mod
   use grid_mod
   implicit none
   type(GRIDTYPE) :: grid
   type(VECTOR) :: posVec, direction, intersection, norm(6)
   real :: d(6),t(6),tval,denom
   real :: tol = 1.00001
   integer :: i,j
   logical :: ok, thisOk(6)

   ok = .true.

   norm(1) = VECTOR(1., 0., 0.)
   norm(2) = VECTOR(0., 1., 0.)
   norm(3) = VECTOR(0., 0., 1.)
   norm(4) = VECTOR(-1., 0., 0.)
   norm(5) = VECTOR(0., -1., 0.)
   norm(6) = VECTOR(0., 0., -1.)


   d(1) = -grid%xAxis(grid%nx)
   d(2) = -grid%yAxis(grid%ny)
   d(3) = -grid%zAxis(grid%nz)
   d(4) = grid%xAxis(1)
   d(5) = grid%yAxis(1)
   d(6) = grid%zAxis(1)


   thisOk = .true.

   do i = 1, 6

    denom = norm(i) .dot. direction
    if (denom /= 0.) then
    t(i) = &
     -( (norm(i) .dot. posVec) + d(i))/denom
      else
     thisOk(i) = .false.
     t(i) = 0.
   endif
      
    if (t(i) < 0.) thisOk(i) = .false.
    if (denom > 0.) thisOK(i) = .false.
   enddo


  do i = 1, 6
    if (thisOk(i)) then
       intersection = posVec + t(i) * direction
       if ( (intersection%x < tol*grid%xAxis(1)) .or. &
            (intersection%y < tol*grid%yAxis(1)) .or. &
            (intersection%z < tol*grid%zAxis(1)) .or. &
            (intersection%x > tol*grid%xAxis(grid%nx)) .or. &
            (intersection%y > tol*grid%yAxis(grid%ny)) .or. &
            (intersection%z > tol*grid%zAxis(grid%nz)) ) then
           thisOk(i) = .false.
       endif
    endif
  enddo

  tval = minval(t, mask=thisOk)

  
  j = 0
  do i = 1, 6
    if (thisOk(i)) j=j+1
  enddo
   
  if (j == 0) ok = .false.


  intersection = posVec + tval * direction

!  write(*,*) tval, intersection,ok
!  write(*,*) t

  intersection%x = min(max(intersection%x, grid%xAxis(1)),grid%xAxis(grid%nx))
  intersection%y = min(max(intersection%y, grid%yAxis(1)),grid%yAxis(grid%ny))
  intersection%z = min(max(intersection%z, grid%zAxis(1)),grid%zAxis(grid%nz))


  end subroutine intersectCube 


