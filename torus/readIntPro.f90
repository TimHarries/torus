subroutine rdintpro(int_pro_file,int_pro,x_int_pro,n_int_pro)
  !
  ! reads in an ascii format intrinsic profile tjh
  !
  implicit none
  character(80) int_pro_file
  real int_pro(*),x_int_pro(*),tot
  real xstart,xend
  integer n_int_pro,i

  if (int_pro_file(1:4).ne.'none') then

     open(20,file=int_pro_file,status='old',form='formatted')

     i=1
10   continue 
     read(20,*,end=30) x_int_pro(i),int_pro(i)
     i=i+1
     goto 10
30   continue
     write(*,*) "closed"
     n_int_pro=i-1
     close(20)

  else

     n_int_pro=100
     xstart=100.
     xend=100000.
     do i=1,n_int_pro
        x_int_pro(i)=xstart+(xend-xstart)*real(i-1)/real(n_int_pro-1)
        int_pro(i)=1.
     enddo
  end if
end subroutine rdintpro
  


