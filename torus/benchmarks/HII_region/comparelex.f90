program compareLex
  implicit none
  integer :: i, j
  real :: r, newT, oldT, newFrac(14), oldFrac(14), val
  logical :: failed

  real, parameter :: tolerance = 0.2
  open(20,file="lexington.dat", form="formatted",status="old")
  open(21,file="lexington_benchmark.dat", form="formatted",status="old")
  failed = .false.
  val = 0.
  do i = 1, 500
     read(20,*) r, newT, newfrac(1:14)
     read(21,*) r, oldT, oldfrac(1:14)
     if ((newT) < 1.) cycle
     if (abs((newT-oldT)/oldT) > tolerance) then
        failed = .true.
        val = max(val, abs((newT-oldT)/oldT))
     endif
     do j = 1, 1
        if ((10.d0**abs(newFrac(j)-oldFrac(j))-1.d0) > tolerance) then
           failed = .true.
           val = max(val,10.d0**abs(newFrac(j)-oldFrac(j))-1.d0)
           write(*,*) newFrac(j),oldFrac(j),val
        endif
     enddo
  enddo
  close(20)
  close(21)
  if (.not.failed) then
     write(*,*) "TORUS: Test successful"
  else
     write(*,*) "TORUS: Test failed ",100.*val, "%"
  endif
end program compareLex


