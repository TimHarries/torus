subroutine readTioCrossSection(lamArray, nLambda, sigmaArray)
  use utils_mod
  use unix_mod
  implicit none
  integer :: nLambda
  real :: lamArray(nLambda), sigmaArray(nLambda)
  integer ::  i, j
  real :: xArray(20000)
  real :: sArray(20000)
  real, allocatable :: nArray(:)
  integer :: npts
  character(len=80) :: dataDirectory, filename

 call unixGetenv("TORUS_DATA",dataDirectory)
 dataDirectory = trim(dataDirectory)//"/"
 filename = trim(dataDirectory)//"tio.xsec"

  open(20, file=filename, status="old", form="formatted")

  npts = 1
  do
   read(20,*,end=10) xArray(npts), sArray(npts)
   npts=npts+1
  end do
10 continue
   close(20)
   npts=npts - 1

   allocate(nArray(1:nLambda))
   narray = 0.

   sigmaArray = 0.
   

   do i = 1, npts
    if ((xArray(i) >= lamArray(1)) .and. &
        (xArray(i)<=lamArray(nLambda))) then
           call locate(lamArray, nLambda, xArray(i), j)
           sigmaArray(j) = sigmaArray(j) + sArray(i)
           nArray(j) = nArray(j) + 1.
    endif
   enddo

   where(nArray > 0) 
    sigmaArray = sigmaArray/nArray
   endwhere

end subroutine readTioCrossSection

