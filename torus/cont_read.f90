!
! very old f77 subroutine which needs updating -

subroutine contread(filename,line_freq,hnu)
  !
  ! this subroutine reads in the continuum flux distribution from a list
  ! formatted file.
  !
  use math_mod
  implicit none
  real nu(2000)         ! frequency grid
  real hnugrid(2000)        ! intensity grid
  real line_freq,hnu
  character(len=80) :: filename
  integer n                        ! number of points
  integer j                        ! loop counter
  !

  if (filename .eq. "none") then
     hnu = 0.
     goto 999
  endif

  write(*,*) "! NB TORUS now expects files in fluxes (erg/s/cm^2/hz)"
  write(*,*) "! i.e. in dipso  (atlasrd, tofnu,xsort)"


  open(unit=11,file=filename,status='old',err=666)
  n=1
30 read(11,*,err=667,end=40) nu(n),hnugrid(n)
  n=n+1
  goto 30
40 continue
  n=n-1
  close(unit=11)
  !
  call locate(nu,n,line_freq,j)
  !
  hnu=logint(line_freq,nu(j),nu(j+1),hnugrid(j),hnugrid(j+1))      
!  write(*,*) '!!!!!!!continuum flux: ',hnu," at ",line_freq
  !
  goto 999
666 write(*,*) 'error opening file: ',trim(filename)
  stop
  !
667 write(*,*) 'error reading file: ',trim(filename)
  stop
  !
999 continue
end subroutine contread
