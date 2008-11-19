program compare_molbench

  implicit none 

  character(len=*), parameter :: model_file="results.dat"
  character(len=*), parameter :: bench_file="moltest.dat"

! Maximum allowable fractional difference 
  real, parameter :: max_diff = 0.1

! Total number of J columns
  integer, parameter :: ncols=8

! No. of columns to check
  integer, parameter :: ncheck=3

  real :: model_R, model_J(ncols)
  real :: bench_R, bench_J(ncols)
  real :: diff(ncols)

  integer :: nlines, status

  open (unit=60, file=bench_file, status='old')
  open (unit=61, file=model_file, status='old')

  nlines=0

  do

     read(60, *, iostat=status) bench_R, bench_J(:)
     if (status /= 0 ) then
        write(*,*) "Reached end of ", bench_file
        write(*,*) "Test passed"
        exit
     end if

     read(61, *, iostat=status) model_R, model_J(:)
     if (status /= 0 ) then
        write(*,*) "Reached end of ", model_file, "before end of ", bench_file
        write(*,*) "Test failed"
        exit
     end if

     diff(:) = abs(model_J(:) - bench_J(:)) / bench_J(:) 

     if ( any(diff(1:ncheck) > max_diff) .and. model_R .lt. 4e17) then 
        write(*,*) "Difference of more than ", max_diff, "found."
        write(*,*) diff(1:ncheck)
        write(*,*) "Test failed"
        exit
     end if

     nlines = nlines + 1

  end do

  close(60)
  close(61)
  
  write(*,*) "Read ", nlines, "lines"

end program compare_molbench
