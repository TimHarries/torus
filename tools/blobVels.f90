program blobVels
  use blob_mod
  use vector_mod
  implicit none
  integer :: i,j,k
  real, parameter :: t1 = 0.,t2 = 24. * 60. * 60.
  real :: thisTime
  character(len=80) :: blobfile, outfile
  integer, parameter:: maxBlobs=100
  type(BLOBTYPE) :: blobs(maxBlobs)
  real :: projVel



    do j=100,80,-1
       write(outfile,'(a,i3.3,a)') "blob",j,".mod"
       open(20,file=outfile,form='formatted',status='unknown')
       do i = 1, 50
          thisTime = t1 + real(i-1)*(t2-t1)/49.

          write(blobfile,'(a,i3.3,a)') "run",i,".blob"
          call readBlobs(blobfile, maxBlobs, blobs)
          
          projVel = blobs(j)%velocity .dot. VECTOR(0.,-1.,0.)
          if (blobs(j)%inUse) then
           write(20,*) thisTime,projVel*1.e5
          endif
       enddo
       close(20)
    enddo
end program blobVels

      
  
  
