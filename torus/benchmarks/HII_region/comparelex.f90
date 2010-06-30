program compareLex
  implicit none
  integer :: i, j
  real :: r, newT, oldT, newFrac(14), oldFrac(14), val
  logical :: failed

  integer :: nTempGood, nIonGood   ! number of points actually used for comparisons

  integer, parameter :: nlines=500     ! number of lines to read from files
  integer, parameter :: nion_check = 1 ! number of species to check 
  real, parameter :: tolerance = 0.33  ! Maximum tolerance
  character(len=*), parameter :: torus_file = "lexington.dat"
  character(len=*), parameter :: ref_file   = "lexington_benchmark.dat"

  open(20,file=torus_file, form="formatted",status="old")
  open(21,file=ref_file, form="formatted",status="old")

  failed = .false.
  val = 0.
  nTempGood = 0 
  nIonGood  = 0 

  write(*,*) "Reading ", nlines, " lines"
  write(*,*) "Torus file: ", torus_file
  write(*,*) "Reference file: ", ref_file
  write(*,*) "Tolerance= ", tolerance

lines:  do i = 1, nlines
     read(20,*) r, newT, newfrac(1:14)
     read(21,*) r, oldT, oldfrac(1:14)
     if ((newT) < 1.) cycle
     if (oldT /= 0.) then
        nTempGood = nTempGood + 1 
        if (abs((newT-oldT)/oldT) > tolerance) then
           failed = .true.
           val = max(val, abs((newT-oldT)/oldT))
           write(*,*) "Fail on temperature:", i, newT, oldT, val
        endif
     endif
     do j = 1, nion_check
        nIonGood = nIonGood + 1 
        if ((10.d0**abs(newFrac(j)-oldFrac(j))-1.d0) > tolerance) then
           failed = .true.
           val = max(val,10.e0**abs(newFrac(j)-oldFrac(j))-1.e0)
           write(*,*) "Fail on fraction: ", i, j, newFrac(j), oldFrac(j), val
        endif
     enddo
  enddo lines

  close(20)
  close(21)
  if (.not.failed) then
     write(*,*) "Compared ", nTempGood, " temperatures"
     write(*,*) "Compared ", nIonGood, "ionization fractions"
     write(*,*) "TORUS: Test successful"
  else
     write(*,*) "TORUS: Test failed ",100.*val, "%"
  endif
end program compareLex


