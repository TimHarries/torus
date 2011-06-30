!THaw - 29/06/2011

program checkFlux

  implicit none
  
  character(len=*), parameter :: filename = "image_flux.dat"
  integer :: ierr, i
  real :: diff
  real :: tolerance = 1.0e-8
  real :: flux(2)
  open(1, file=filename, status="old", iostat=ierr) 
  
  if(ierr /= 0) then
     print *, "Error opening "//filename
     print *, "Test failed"
     stop
  end if

  print *, "Tolerance ", tolerance

  do i = 1, 2
     read(1, *, iostat=ierr) flux(i)
     if(ierr /= 0) then
        print *, "Error reading from "//filename
        stop
     end if
  end do
     
  diff = abs(flux(1)-flux(2))/flux(2)

  print *, "Fractional difference ", diff
  if(diff < tolerance) then
     print *, "TORUS: Test successful"
  else
     print *, "TORUS: Test failed"
  end if

end program
