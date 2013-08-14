#ifdef PDR
module nrayshealpix



implicit none


contains

subroutine donrayshealpix()
  use definitions
  use healpix
  use healpix_types
  use healpix_module
  use inputs_mod, only : hlevel
  
  implicit none

  integer :: i
  
  level=hlevel !<---- ONLY THIS IS NEEDED FOR INPUT PARAMETERS
  nside=2**level
  nrays = 12*nside**2
  ns_max=8192
  
  allocate(vector(1:3))
  allocate(vertex(1:3,1:4))
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
  

end subroutine donrayshealpix
end module
#endif
