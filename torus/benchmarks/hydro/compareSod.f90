program compareSod
  implicit none
  integer, parameter :: nSodA = 19
  real :: xA(nSodA), rhoA(nSodA), testrho,fac
  integer :: nTorus
  real :: tot, percent
  integer :: i,j 
  real, allocatable :: xTorus(:), rhoTorus(:)
  open(20,file="sod_analytical.dat", status="old", form="formatted")
  do i = 1, nSodA
     read(20,*) xA(i), rhoA(i)
  enddo
  close(20)
  open(20,file="sod.dat", status="old", form="formatted")
  read(20,*) nTorus
  allocate(xTorus(1:nTorus), rhoTorus(1:nTorus))
  do i = 1, nTorus
     read(20,*) xTorus(i), rhoTorus(i)
  enddo
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
