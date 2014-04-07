program compareLex
  implicit none
  integer :: i, j
  integer, parameter :: nIonsOld=14
  integer, parameter :: nIonsNew=11
  real :: r, newT, oldT, newFrac(nIonsNew), oldFrac(nIonsOld), val
  logical :: failed

  integer :: nTempGood, nIonGood   ! number of points actually used for comparisons

  integer, parameter :: nlines=500     ! number of lines to read from files
  integer, parameter :: nion_check = 1 ! number of species to check 
  real, parameter :: tolerance = 0.10  ! Maximum tolerance
  real :: fac1, fac2, avcheck, tav, junk
  integer :: nav, nt
  character(len=*), parameter :: torus_file = "lexington.dat"
  character(len=*), parameter :: ref_file   = "lexington_benchmark.dat"

  integer :: status

  open(20,file=torus_file, form="formatted",status="old", iostat=status)
  if ( status /= 0 ) then 
     write(*,*) "Error opening "//torus_file
     stop
  end if

  open(21,file=ref_file, form="formatted",status="old", iostat=status)
  if ( status /= 0 ) then 
     write(*,*) "Error opening "//ref_file
     stop
  end if

  failed = .false.
  val = 0.
  nTempGood = 0 
  nIonGood  = 0 

  write(*,*) "Reading ", nlines, " lines"
  write(*,*) "Torus file: ", torus_file
  write(*,*) "Reference file: ", ref_file
  write(*,*) "Tolerance= ", tolerance

! Check the number of species we are comparing
  if ( nion_check > min(nIonsOld,nIonsNew) ) then
     write(*,*) "Error: trying to compare more species than are available"
     STOP
  else
     write(*,'(a,i2,a)') "Comparing ", nion_check, " ion(s)"
  end if

  avcheck = 0.
  nav = 0
  tav = 0.
  nt = 0
lines:  do i = 1, nlines
     read(20,*,iostat=status) r, newT, newfrac(1:nIonsNew)
     if ( status /= 0 ) then 
        write(*,*) "Error reading from "//torus_file
        stop
     end if
     read(21,*) r, oldT, oldfrac(1:nIonsOld)
     if ( status /= 0 ) then 
        write(*,*) "Error reading from "//ref_file
        stop
     end if
     if ((newT) < 1.) cycle
     if (oldT /= 0.) then
        nTempGood = nTempGood + 1 
        tav = tav + abs((newT-oldT)/oldT)
        nt = nt + 1
     endif
     do j = 1, nion_check
        nIonGood = nIonGood + 1 
        fac1 = 10.d0**newFrac(j)
        fac2 = 10.d0**oldFrac(j)
        avcheck = avcheck + abs((fac1-fac2)/fac2)
        nav = nav + 1
     enddo
  enddo lines
  avcheck = avcheck / real(nav)
  tav = tav / real(nt)
  write(*,*) "Average percentage difference in ions is ",100.*avcheck
  write(*,*) "Average percentage difference in temp is ",100.*tav
  if (avcheck > tolerance) then
     failed =.true.
     write(*,*) "Failed on ions"
  endif
  if (tav > tolerance) then
     failed =.true.
     write(*,*) "Failed on temperatures"
  endif


  close(20)
  close(21)
  if (.not.failed) then
     write(*,*) "Compared ", nTempGood, " temperatures"
     write(*,*) "Compared ", nIonGood, "ionization fractions"
     write(*,*) "TORUS: Test successful"
  else
     write(*,*) "TORUS: Test failed"
  endif
end program compareLex


