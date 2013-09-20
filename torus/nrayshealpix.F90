#ifdef PDR
module nrayshealpix
use constants_mod
use messages_mod
use parallel_mod
use gridio_mod
use source_mod
use timing
use grid_mod
use vtk_mod
use amr_mod
use mpi_amr_mod
use mpi_global_mod


implicit none


contains

subroutine donrayshealpix()
#ifdef MPI
  use mpi
  
#endif
!  use definitions
  use healpix
  use constants_mod
!  use healpix_types
  use healpix_guts
!  use healpix_module
  use inputs_mod, only : hlevel
  
  implicit none

  integer :: i

  real(kind=dp) :: vector(1:3)
  integer(kind=i4b) :: nside, nrays
!  integer(kind=i4b) ::pix2x(0:1023),pix2y(0:1023)
!  REAL(KIND=DP) :: vector(1:3)
  REAL(KIND=DP) :: vertex(1:3,1:4)
  integer :: ipix
  level=hlevel !<---- ONLY THIS IS NEEDED FOR INPUT PARAMETERS
  nside=2**level
  nrays = 12*nside**2
  ns_max=8192
  

!  allocate(vertex(1:3,1:4))
#ifdef MPI
  allocate(vectors(1:3,0:nrays-1))
  
  
  call mk_xy2pix()
  
  if(myrankglobal == 1) then
     write(6,*) 'Building HEALPix vectors...'
     open(unit=77,file='HEALPix_vectors.dat',status='replace')
  end if

  do i=1,nrays
     ipix=i-1 !ipix is the ID of a HEALPix ray. Runs with values 0:nrays-1
     call pix2vec_nest(nside,ipix,pix2x,pix2y,vector,vertex)
     vectors(1:3,ipix)=vector(1:3) !Store in memory
     if(myrankglobal == 1) then
        write(77,'(3E15.7,I7)') vectors(1:3,ipix),ipix
     endif
  enddo
  if(myrankglobal == 1) then
     Write(6,*) 'Done!';write(6,*) ''
  endif

#else

  allocate(vectors(1:3,0:nrays-1))  
  call mk_xy2pix()
  
  write(6,*) 'Building HEALPix vectors...'
  open(unit=77,file='HEALPix_vectors.dat',status='replace')
  do i=1,nrays
     ipix=i-1 !ipix is the ID of a HEALPix ray. Runs with values 0:nrays-1
     call pix2vec_nest(nside,ipix,pix2x,pix2y,vector,vertex)
     vectors(1:3,ipix)=vector(1:3) !Store in memory
     write(77,'(3E15.7,I7)') vectors(1:3,ipix),ipix
  enddo
  
  Write(6,*) 'Done!';write(6,*) ''
#endif

end subroutine donrayshealpix
end module
#endif
