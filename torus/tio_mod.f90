module tio_mod
implicit none

public

contains

  readLineList(filename, nLambda, vLow, vUpper, Jlow, gf, lambda, excitation)

    character(len=*) :: filename
    character(len=4) :: junk
    integer :: vLow(*), vUpper(*), Jlow(*)
    real :: gf(*)
    real :: lambda(*)
    real :: excitation(*)
    integer :: i, nLambda


    open(20,file=filename,status="old",form="formatted")

    nLambda = 1

10  continue
     read(20,'(i2,i2,i3,x,a4,e9.2,f10.3,f6.3)',end=30) &
      vUp(nLambda), vLow(nLambda), Jlow(nLambda), junk, gf(nLambda), &
      lambda(nLambda), excitation(nLambda)
     nLambda = nLambda + 1
    goto 10
30  continue
    nLambda = nLambda - 1
    close(20)
  end subroutine readLineList

end module tio_mod