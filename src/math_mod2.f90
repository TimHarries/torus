module math_mod2

  ! some routines from numerical recipes. Kept in a 
  !   second math_mod file to avoid circular-dependency 
  !   problems.

  use kind_mod
  use random_mod
  implicit none

  INTERFACE ludcmp_f77
    MODULE PROCEDURE ludcmp_f77_single, ludcmp_f77_double
  END INTERFACE

  INTERFACE lubksb_f77
    MODULE PROCEDURE lubksb_f77_single, lubksb_f77_double
  END INTERFACE

contains
  
     subroutine ludcmp_f77_single(a,n,np,indx,d)
      real, parameter :: tiny=1.0e-20
      integer, parameter :: nmax=200
      integer, intent(in) :: n, np
      real, intent(out) :: d
      real,intent(inout) ::  a(np,np)
      real ::  vv(nmax)
      integer, intent(out) :: indx(n)
      integer :: i
      real(double) :: aamax
      real(double) :: sum, dum
      integer :: j, k
      integer :: imax
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) then
           write(*,*) 'Singular matrix in ludcmp_f77_single.'
           write(*,*) 'Correction forced with random number!'
           call randomNumberGenerator(getRealArray2d=a)
           aamax=abs(MAXVAL(a(:,:)))
!           stop
        endif
        vv(i)=real(1./aamax)
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=real(sum)
            endif
14        continue
        endif
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=real(sum)
          endif
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=real(dum)
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n)then
          if(a(j,j).eq.0.)a(j,j)=tiny
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=real(a(i,j)*dum)
18        continue
        endif
19    continue
      if(a(n,n).eq.0.)a(n,n)=tiny
      end subroutine ludcmp_f77_single

    subroutine ludcmp_f77_double(a,n,np,indx,d)
      real(double), parameter :: tiny=1.0e-20_db
      integer, parameter :: nmax=200
      integer, intent(in) :: n, np
      real(double), intent(out) :: d
      real(double),intent(inout) ::  a(np,np)
      real(double) ::  vv(nmax)
      integer, intent(out) :: indx(n)
      integer :: i
      real(double) :: aamax
      real(double) :: sum, dum
      integer :: j, k
      integer :: imax

      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.d0) then
           write(*,*) 'Singular matrix in ludcmp_f77_double.'
           write(*,*) 'Correction forced with random number!'
           call randomNumberGenerator(getDoubleArray2d=a)
           aamax=abs(MAXVAL(a(:,:)))
!           stop
        endif
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n)then
          if(a(j,j).eq.0.d0)a(j,j)=tiny
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      if(a(n,n).eq.0.)a(n,n)=tiny
      end subroutine ludcmp_f77_double

      pure subroutine lubksb_f77_single(a,n,np,indx,b2)
        integer, intent(in) :: np,n
        real,intent(in) ::  a(np,np)
        integer, intent(in) :: indx(n)
        real, intent(inout) :: b2(n)
        integer :: ii, i, ll, j
        real :: sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b2(ll)
        b2(ll)=b2(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b2(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b2(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b2(i)
        if(i.lt.n)then
          do 13 j=i+1,n
            sum=sum-a(i,j)*b2(j)
13        continue
        endif
        b2(i)=sum/a(i,i)
14    continue
      end subroutine lubksb_f77_single
      
      pure subroutine lubksb_f77_double(a,n,np,indx,b2)
        integer, intent(in) :: np,n
        real(double),intent(in) ::  a(np,np)
        integer, intent(in) :: indx(n)
        real(double), intent(inout) :: b2(n)
        integer :: ii, i, ll, j
        real(double) :: sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b2(ll)
        b2(ll)=b2(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b2(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b2(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b2(i)
        if(i.lt.n)then
          do 13 j=i+1,n
            sum=sum-a(i,j)*b2(j)
13        continue
        endif
        b2(i)=sum/a(i,i)
14    continue
      end subroutine lubksb_f77_double

end module math_mod2
