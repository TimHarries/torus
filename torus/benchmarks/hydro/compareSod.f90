program compareSod
  implicit none

  character(len=*), parameter :: torus_file="sod.dat"
  character(len=*), parameter :: ref_file  ="sod_analytical.dat"
  integer, parameter :: nSodA = 19
  real :: xA(nSodA), rhoA(nSodA), testrho,fac
  integer :: nTorus
  real :: tot, percent
  integer :: i,j, status 
  integer :: nTorusMax = 100000
  real, allocatable :: xTorus(:), rhoTorus(:)
  logical :: found

  open(20,file=ref_file, status="old", form="formatted",iostat=status)
  if ( status /= 0) call file_read_error(ref_file)
  do i = 1, nSodA
     read(20,*) xA(i), rhoA(i)
  enddo
  close(20)

  inquire(file=torus_file, exist=found )
  if ( found ) then 
     open(20,file=torus_file, status="old", form="formatted")
     if ( status /= 0) call file_read_error(torus_file) 
  else
     write(*,*) torus_file//" not found"
     write(*,*) "TORUS: Test failed"
     STOP
  end if

  allocate(xTorus(1:nTorusMax), rhoTorus(1:nTorusMax))
  nTorus = 1
10 continue
  read(20,*,end=20) xTorus(nTorus), rhoTorus(nTorus)
  nTorus = nTorus + 1
  goto 10
20 continue
  nTorus = ntorus - 1
  close(20)
  tot = 0.d0
  percent = 0.d0
  do i = 1, nTorus
     call locate(xA, nsodA, xTorus(i), j)
     fac = (xTorus(i) - xA(j))/(xa(j+1)-xA(j))
     testRho = rhoa(j) + fac * (rhoa(j+1)-rhoa(j))
     percent = max(percent, abs(100.*(testRho-RhoTorus(i))/rhoTorus(i)))
     tot = tot + (testRho-rhoTorus(i))**2
  enddo
  tot = sqrt(tot/(real(nTorus-1)))
  write(*,*) "Sigma: ",tot
  if (tot < 1.d-2) then
     write(*,*) "TORUS: Test successful"
     
  else
     write(*,*) "TORUS: Test failed"
  endif
end program compareSod
SUBROUTINE LOCATE(XX,N,X,J)
  implicit none
  real XX(*), x
  integer :: n, j, jm, jl, ju
  JL=0
  JU=N+1
10 IF(JU-JL.GT.1)THEN
     JM=(JU+JL)/2
     IF((XX(N).GT.XX(1)).EQV.(X.GE.XX(JM)))THEN
        JL=JM
     ELSE
        JU=JM
     ENDIF
     GO TO 10
  ENDIF
  J=JL
END SUBROUTINE LOCATE

subroutine file_read_error(filename)

  implicit none

  character(len=*), intent(in) :: filename

  write(*,*) "Error reading from "//trim(filename)
  write(*,*) "TORUS: Test failed"
  STOP

end subroutine file_read_error

