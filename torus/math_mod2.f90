module math_mod2

  ! some routines from numerical recipes. Kept in a 
  !   second math_mod file to avoid circular-dependency 
  !   problems.

  use kind_mod
  
  implicit none

  interface zriddr
    module procedure zriddr_single
    module procedure zriddr_double
  end interface 

  INTERFACE ludcmp_f77
    MODULE PROCEDURE ludcmp_f77_single, ludcmp_f77_double
  END INTERFACE

  INTERFACE lubksb_f77
    MODULE PROCEDURE lubksb_f77_single, lubksb_f77_double
  END INTERFACE

  INTERFACE mnewt
    MODULE PROCEDURE mnewt_single, mnewt_double
  END INTERFACE

contains
  
  REAL FUNCTION zriddr_single(func,x1,x2,xacc)

    REAL, INTENT(IN) :: x1, x2, xacc

    INTERFACE
      FUNCTION func(x)
        IMPLICIT NONE
        REAL, INTENT(IN) :: x
        REAL :: func
      END FUNCTION func
    END INTERFACE

    INTEGER, PARAMETER :: MAXIT = 60
    ! Using Ridders' method, return the root of a function func known to lie between x1 and
    ! x2. The root, returned as zriddr, will be refined to an approximate accuracy xacc.
    REAL, PARAMETER :: UNUSED = -1.11e30
    INTEGER :: j
    REAL :: fh, fl, fm, fnew, s, xh, xl, xm, xnew
    fl = func(x1)
    fh = func(x2)
    IF ((fl >0.0 .AND. fh <0.0).OR.(fl <0.0 .and.fh >0.0)) THEN
      xl = x1
      xh = x2
      zriddr_single = UNUSED ! Any highly unlikely value, to simplify logic below.
      DO j =1, MAXIT
        xm = 0.5 * (xl+xh)
        fm = func(xm) ! First of two function evaluations per iteration.
        s = SQRT(fm**2-fl*fh)
        IF (s ==0.0) RETURN
        xnew = xm + (xm-xl) * (sign(1.0,fl-fh)*fm/s) ! Updating formula.
        IF (ABS(xnew-zriddr_single)<=xacc )RETURN
        zriddr_single = xnew
        fnew = func(zriddr_single) ! Second of two function evaluations per
        ! iteration.
        IF (fnew ==0.0) RETURN
        IF (SIGN(fm,fnew)/=fm) THEN ! Bookkeeping to keep the root bracketed
          ! on next iteration.
          xl = xm
          fl = fm
          xh = zriddr_single
          fh = fnew
        ELSE IF (SIGN(fl,fnew)/=fl) THEN
          xh = zriddr_single
          fh = fnew
        ELSE IF (SIGN(fh,fnew)/=fh) THEN
          xl = zriddr_single
          fl = fnew
        ELSE
          write (*,*) 'zriddr:never get here '
        END IF
        IF (ABS(xh-xl)<=xacc) RETURN
      END DO
      write (*,*) 'zriddr:exceeded maximum iterations '
    ELSE IF (fl ==0.0) THEN
      zriddr_single = x1
    ELSE IF (fh ==0.0) THEN
      zriddr_single = x2
    ELSE
      write (*,*) 'zriddr:root must be bracketed '
    END IF
  END FUNCTION zriddr_single
  
  real(double) FUNCTION zriddr_double(func,x1,x2,xacc)

    real(double), INTENT(IN) :: x1, x2, xacc

    INTERFACE
      FUNCTION func(x)
        USE kind_mod
        IMPLICIT NONE
        real(double), INTENT(IN) :: x
        real(double) :: func
      END FUNCTION func
    END INTERFACE

    INTEGER, PARAMETER :: MAXIT = 60
    ! Using Ridders' method, return the root of a function func known to lie between x1 and
    ! x2. The root, returned as zriddr, will be refined to an approximate accuracy xacc.
    real(double), PARAMETER :: UNUSED = -1.11e299_db
    INTEGER :: j
    real(double) :: fh, fl, fm, fnew, s, xh, xl, xm, xnew
    fl = func(x1)
    fh = func(x2)
    IF ((fl >0.0 .AND. fh <0.0_db).OR.(fl <0.0_db .and.fh >0.0_db)) THEN
      xl = x1
      xh = x2
      zriddr_double = UNUSED ! Any highly unlikely value, to simplify logic below.
      DO j =1, MAXIT
        xm = 0.5_db * (xl+xh)
        fm = func(xm) ! First of two function evaluations per iteration.
        s = SQRT(fm**2-fl*fh)
        IF (s ==0.0_db) RETURN
        xnew = xm + (xm-xl) * (sign(1.0_db,fl-fh)*fm/s) ! Updating formula.
        IF (ABS(xnew-zriddr_double)<=xacc )RETURN
        zriddr_double = xnew
        fnew = func(zriddr_double) ! Second of two function evaluations per
        ! iteration.
        IF (fnew ==0.0_db) RETURN
        IF (SIGN(fm,fnew)/=fm) THEN ! Bookkeeping to keep the root bracketed
          ! on next iteration.
          xl = xm
          fl = fm
          xh = zriddr_double
          fh = fnew
        ELSE IF (SIGN(fl,fnew)/=fl) THEN
          xh = zriddr_double
          fh = fnew
        ELSE IF (SIGN(fh,fnew)/=fh) THEN
          xl = zriddr_double
          fl = fnew
        ELSE
          write (*,*) 'zriddr:never get here '
        END IF
        IF (ABS(xh-xl)<=xacc) RETURN
      END DO
      write (*,*) 'zriddr:exceeded maximum iterations '
    ELSE IF (fl ==0.0_db) THEN
      zriddr_double = x1
    ELSE IF (fh ==0.0_db) THEN
      zriddr_double = x2
    ELSE
      write (*,*) 'zriddr:root must be bracketed '
    END IF
  END FUNCTION zriddr_double

    SUBROUTINE mnewt_single(ntrial,x,tolx,tolf,usrfun)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ntrial 
    REAL, DIMENSION(:), INTENT(INOUT) :: x
    REAL, INTENT(IN) :: tolx,tolf 
  
    INTERFACE 
      SUBROUTINE usrfun(x,fvec,fjac)
        IMPLICIT NONE 
        REAL, DIMENSION(:), INTENT(IN) :: x 
        REAL, DIMENSION(:), INTENT(OUT) :: fvec 
        REAL, DIMENSION(:,:), INTENT(OUT) :: fjac 
      END SUBROUTINE usrfun 
    END INTERFACE 
    !Given an initial guess x for a root in N dimensions, take ntrial Newton-Raphson 
    !  steps to improve the root. Stop if the root converges in either summed
    !  absolute variable increments tolx or summed absolute function values tolf. 
    
    INTEGER :: i 
    INTEGER, DIMENSION(size(x)) :: indx 
    REAL :: d 
    REAL, DIMENSION(size(x)) :: fvec,p 
    REAL, DIMENSION(size(x),size(x)) :: fjac 
  
    do i=1,ntrial 
    
      call usrfun(x,fvec,fjac) 
      ! User subroutine supplies function values at x in fvec and Jacobian matrix
      !   in fjac. If your usrfun calls fdjac to get the Jacobian matrix 
      !   numerically, you must change its interface so that x is INTENT(INOUT),
      !   not INTENT(IN). 
      
      if (sum(abs(fvec)) <= tolf) RETURN 
      !Check function convergence. 
      p=-fvec ! Right-hand side of linear equations. 
      call ludcmp_f77(fjac,SIZE(fjac,1),SIZE(fjac,1),indx,d) !Solve linear equations using LU decomposition. 
      call lubksb_f77(fjac,SIZE(fjac,1),SIZE(fjac,1),indx,p) 
      x = x + p !Update solution. 
      if (sum(abs(p)) <= tolx) RETURN !Check root convergence.
      
    end do 
  
  END SUBROUTINE mnewt_single
  
  SUBROUTINE mnewt_double(ntrial,x,tolx,tolf,usrfun)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ntrial 
    real(double), DIMENSION(:), INTENT(INOUT) :: x
    real(double), INTENT(IN) :: tolx,tolf 
    INTERFACE 
      SUBROUTINE usrfun(x,fvec,fjac)
        use kind_mod
        IMPLICIT NONE 
        real(double), DIMENSION(:), INTENT(IN) :: x 
        real(double), DIMENSION(:), INTENT(OUT) :: fvec 
        real(double), DIMENSION(:,:), INTENT(OUT) :: fjac 
      END SUBROUTINE usrfun 
    END INTERFACE 
    !Given an initial guess x for a root in N dimensions, take ntrial Newton-Raphson 
    !  steps to improve the root. Stop if the root converges in either summed
    !  absolute variable increments tolx or summed absolute function values tolf. 
    
    INTEGER :: i 
    INTEGER, DIMENSION(size(x)) :: indx 
    real(double) :: d 
    real(double), DIMENSION(size(x)) :: fvec,p 
    real(double), DIMENSION(size(x),size(x)) :: fjac 
  
    do i=1,ntrial 

      call usrfun(x,fvec,fjac) 
      ! User subroutine supplies function values at x in fvec and Jacobian matrix
      !   in fjac. If your usrfun calls fdjac to get the Jacobian matrix 
      !   numerically, you must change its interface so that x is INTENT(INOUT),
      !   not INTENT(IN). 
      
      if (sum(abs(fvec)) <= tolf) RETURN 
      !Check function convergence. 
      p=-fvec ! Right-hand side of linear equations. 
      call ludcmp_f77(fjac,SIZE(fjac,1),SIZE(fjac,1),indx,d) !Solve linear equations using LU decomposition. 
      call lubksb_f77(fjac,SIZE(fjac,1),SIZE(fjac,1),indx,p) 
      x = x + p !Update solution. 
      if (sum(abs(p)) <= tolx) RETURN !Check root convergence.
      
    end do 
  
  END SUBROUTINE mnewt_double
  
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
      real(double) :: sum, dum, r
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
           call random_number(a(:,:))
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
        aamax=0.
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
          if(a(j,j).eq.0.)a(j,j)=tiny
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
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
      real(double) :: sum, dum, r
      integer :: j, k
      integer :: imax
      logical :: error
      error = .false.
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.d0) then
           write(*,*) 'Singular matrix in ludcmp_f77_double.'
           write(*,*) 'Correction forced with random number!'
           call random_number(a(:,:))
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
