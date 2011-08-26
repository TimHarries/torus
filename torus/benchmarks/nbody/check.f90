program check

  implicit none
  integer :: i,j
  real(kind=8) :: pe(100), ke(100), etot(100),t
  character(len=10) :: junk


  open(20,file="energy.dat", status="old", form="formatted")
  read(20,'(a)') junk

  do i = 1, 100
     read(20,*) j,t,pe(i),ke(i),etot(i)
  enddo
  close(20)
  if (abs((etot(1)-etot(100))/etot(100)) > 1.d-5) then
     write(*,*) "Torus nbody test failed"
  else
     write(*,*) "Torus nbody test successful"
  endif
end program check
