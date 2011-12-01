  interface zriddr
    module procedure zriddr_single
    module procedure zriddr_double
  end interface 

  INTERFACE mnewt
    MODULE PROCEDURE mnewt_single, mnewt_double
  END INTERFACE


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
  
