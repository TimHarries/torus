program imconv
  implicit none
  
  character(LEN=80) ::  filename, datfilename, sdffilename
  integer :: i1,nSize, nPix, i, j, nphase,iphase
  real :: scale
  !
  real, allocatable :: thisImageI(:,:)   
  real, allocatable :: thisImageQ(:,:)   
  real, allocatable :: thisImageU(:,:)   
  real, allocatable :: thisImageV(:,:)   
  real, allocatable :: thisImageVel(:,:) 
  
  write(*,*) "filename?"
  read(*,'(a)') filename
  
  write(*,*) "nphase?"
  read(*,*) nphase
  
  i1 = index(filename," ")
  do iphase = 1, nphase
     write(datfilename,'(a,i3.3,a)')  filename(1:i1-1),iPhase,".dat"
     write(sdffilename,'(a,i3.3)')  filename(1:i1-1),iPhase
     
     write(*,*) datfilename,sdffilename
     open(20, file=datfilename, status="old", form="formatted")
     
     read(20,'(2I4)') nSize
     
     if (iphase==1) then 
        allocate(thisImageI(-nsize:nsize, -nsize:nsize))
        allocate(thisImageQ(-nsize:nsize, -nsize:nsize))
        allocate(thisImageU(-nsize:nsize, -nsize:nsize))
        allocate(thisImageV(-nsize:nsize, -nsize:nsize))
        allocate(thisImageVel(-nsize:nsize, -nsize:nsize))
     end if
     
     
     nPix = 2*nSize + 1
         
     read(20, '(1p,e13.5)') scale
     
     do i = -nSize, nSize
        do j = -nSize, nSize
           read(20, '(1p,5e13.5)') thisImageI(i,j), &
           &                       thisImageQ(i,j), &
           &                       thisImageU(i,j), &
           &                       thisImageV(i,j), &
           &                       thisImageVel(i,j) 
        enddo
     enddo
     close(20)

     
     call writeImagef77(thisImageI, thisImageQ, thisImageU, &
          thisImageV, thisImageVel, scale, nSize, nPix, &
          sdffilename)
  enddo

end program imconv


