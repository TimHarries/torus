program compare_molbench

  implicit none 

  character(len=*), parameter :: model_file="results.dat"
  character(len=*), parameter :: bench_file="moltest.dat"

! Maximum allowable fractional difference 
  real :: max_diff,diffmax

! Total number of J columns
  integer, parameter :: ncols=8

! No. of columns to check
  integer, parameter :: ncheck=6

  real :: model_R, model_J(ncols), model_Rarray(100)
  real :: bench_R, bench_J(ncols)
  real :: diff(ncols,200)

  integer :: diffmaxloc(2)
  integer :: nlines, status

  max_diff = 0.01 * sqrt(89.) ! sqrt(Nvoxels)
  diff = -999.
  diffmax = -1.

  open (unit=60, file=bench_file, status='old')
  open (unit=61, file=model_file, status='old')

  nlines=0

  do
     nlines = nlines + 1

     read(60, *, iostat=status) bench_R, bench_J(:)
     if (status /= 0 ) then
        write(*,*) "Reached end of ", bench_file
        exit
     end if

     read(61, *, iostat=status) model_R, model_J(:)
     model_Rarray(nlines) = model_R
     if (status /= 0 ) then
        write(*,*) "Reached end of ", model_file, "before end of ", bench_file
        write(*,*) "TORUS: Test failed"
        status = -1
        exit
     end if

     diff(:,nlines) = abs(model_J(:) - bench_J(:)) / bench_J(:) 
     diffmax = maxval(diff(:,:))
     diffmaxloc = maxloc(diff(:,:))

     write(*,'(6(tr2,es12.5))') diff(1:ncheck,nlines)

  end do

  write(*,*) "Maximum difference = ",diffmax
  write(*,*) "Radius = ", model_Rarray(diffmaxloc(2))
  write(*,*) "Level ", diffmaxloc(1) - 1

  if(status .ne. -1) then
     if ( any(diff(1:ncheck,:) > max_diff) .and. model_R .lt. 4e17) then 
       write(*,*) "Difference of more than ", max_diff, "found."
       write(*,*) "TORUS: Test failed"
        status = 0
     else
       write(*,*) "Difference of more than ", max_diff, "not found."
       write(*,*) "TORUS: Test successful"
        status = 1
     endif
  endif
  
  write(*,*) "Read ", nlines-1, "lines"

  close(60)
  close(61)

end program compare_molbench
