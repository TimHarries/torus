program check

  implicit none

  logical :: fail

  fail=.false.

  write(*,*) "Checking angular image results"

! Check that all the fits files have been generated successfully
  call checkFile("nCol_angImgTest.fits")
  call checkFile("intensity_angImgTest.fits")
  call checkFile("intensity_neg_angImgTest.fits")
  call checkFile("intensity_pos_angImgTest.fits")
  call checkFile("tau_angImgTest.fits")

    if (fail) then 
       write(*,*) "TORUS: Test failed"
    else
       write(*,*) "TORUS: Test successful"
    endif

contains

  subroutine checkFile(thisFile)

    character(len=*), intent(in) :: thisFile
    logical :: foundFile

! Check that the file exists 
    inquire(file=thisFile, exist=foundFile)
    if (foundFile) then 
       write(*,*) "Found "//thisFile
    else
       write(*,*) "Did not find "//thisFile
       fail=.true.
    endif

! Add more checking here ...

  end subroutine checkFile

end program check

