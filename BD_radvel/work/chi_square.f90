module chi_square

  ! Collection of functions for computing chi_square and chi-quare 
  ! probability

  
  public :: chi_sq_dist, chi_sq_prob, chi_sq
  private :: gamma, simpson

contains
  

  !
  !
  ! 
  double precision function chi_sq(y, sigma, yfit, n)
    implicit none
    integer         , intent(in) :: n          ! number of data points
    double precision, intent(in) :: y(n)       ! data points
    double precision, intent(in) :: sigma(n)
    ! fitted values model evaluated at the data x values.
    double precision, intent(in) :: yfit(n) 
    !
    double precision :: tmp
    integer :: i 
    
    chi_sq = 0.0d0
    do i = 1, n
       tmp = (y(i)-yfit(i))/sigma(i)
       chi_sq = chi_sq + tmp*tmp
    end do
    
  end function chi_sq
  
  !
  ! The probability distribution function for chi-sq
  ! with "nfree" degrees of freedom.
  function chi_sq_dist(chi_sq,nfree) RESULT(out)
    implicit none
    double precision :: out 
    integer, intent(in) ::  nfree
    double precision, intent(in) ::  chi_sq
    !
    double precision :: nu
    !
    nu = nfree
    out  = ( chi_sq**((nu-2.0d0)/2.0d0) )  &
         /                                        &
         ( (2.0d0**(nu-2.0d0))*gamma(nu/2.0d0)*dexp(chi_sq/2.0d0) )
    !
    return
  end function chi_sq_dist
  
  
  
  
  !====================================================================
  !  This function calclulats the chi-squre probability.
  !  (c.f. "Datareduction and Error Analysis for the Physical
  !  Sciences" by P. Bevington, page 193.
  !===================================================================
  ! Created: 18-aug-1998 : Kurosawa  
  !====================================================================
  double precision function chi_sq_prob(chi_sq,nfree)
    implicit none
    !
    double precision, intent(in) :: chi_sq
    integer, intent(in) :: nfree 
    !
    double precision ::  z, term, sum,  f1
    integer :: imax, i
    logical :: nfree_is_even
    !Integration stuff
    integer, parameter :: nint=1000     ! # of integration pts
    double precision  ::  xi, xf
    !
    ! cheking the validity of the inputs
    if ( chi_sq.lt.0.0d0) then 
       print*, 'Error in ''chi_sq_prob'' : chi_sq must be positive.'
       stop
    end if
!    if (nfree.le.1) then 
    if (nfree.lt.1) then 
       print*, 'Error in ''chi_sq_prob'' : invalid value of nfree'
       stop
    end if
    ! cheking if nfree is even
    if (mod(nfree,2).eq.0) then
       nfree_is_even = .true.
    else
       nfree_is_even = .false.
    end if
    !
!    z = chi_sq/2.0d0   ! (Bug! Found on 18-Apr-2006, R. Kurosawa)
    z = chi_sq*dble(nfree)/2.0d0
     !
    !  Case for nfree even: We can use the summation formula
    if (nfree_is_even) then
       imax = nfree/2
       term = 1.0d0
       sum = 0.0d0
       do i = 1, imax
          f1 = i
          sum = sum + term
          term = term*z/f1
       end do
       chi_sq_prob = sum*exp(-z)
       ! Case when nfree is odd but z is a large number.
       !Perform numerical integration like an idiot.
    else
       xi=chi_sq           ! Lower integration limit
       xf=3.d0*chi_sq    ! Upper integration limit
       call simpson(nfree,xi,xf,nint,chi_sq_prob)
    end if
    !
    return
  end function chi_sq_prob


  !*********************************************************************
  ! Calculates gamma function using polynominal approximation formula
  ! See 6.1.35. on page 257 of Abromwitz and Stegun.
  ! For simplicity x must be "real" number and .ge. 1.
  ! The size of error: abs(err(x)) =< 3.0e-07
  !*********************************************************************
  double precision function gamma(z)
    implicit none
    double precision, intent(in) :: z
    double precision :: x,y,x_orig
    double precision :: a1,a2,a3,a4,a5,a6,a7,a8
    parameter(a1=-0.577191652d0, a2=0.988205891d0,  &
         &    a3=-0.897056937d0, a4=0.918206857d0,  &
         &    a5=-0.756704078d0, a6=0.482199394d0,  &
         &    a7=-0.193527871d0, a8=0.035868347d0)
    !
    x = z 
    x_orig = x
    !
    if (x.le.0.0d0) then
       write(*,*) 'Argument must be greater than 1. [function gamma]'
    elseif (x.gt.0.0d0.and.x.lt.1.0d0) then
       x = x + 1  ! then use the relation gamma(x) = gamma(x+1)/x
       gamma = 1.0d0
    elseif (x.gt.2.0d0) then
       gamma = 1.0d0
       do while (x.ge.2.0d0)
          gamma = gamma*(x-1.0d0) 
          x = x - 1.0d0
       end do
    elseif (x.ge.1.0d0.and.x.le.2.0d0) then
       gamma = 1.0d0 
    else
       print *, 'Wrong argument [function gamma]'
    endif
    !	
    ! Approximation formula valid for 0 <= y <= 1. |err(x)| <= 5e-7
    ! Warning: gamma = Gamma(y + 1). Also note: y = x - 1
    
    y = x - 1.0d0
    gamma = gamma * (1.0d0 + a1*y + a2*(y**2) + a3*(y**3) +   &
         a4*(y**4) + a5*(y**5) + a6*(y**6) + &
         a7*(y**7) + a8*(y**8)  )

    if (x_orig.gt.0.0d0.and.x_orig.lt.1.0d0) then
       gamma = gamma/x_orig
    end if
    !
    return
  end function gamma




  !*************************************************************************
  !  This subroutine calcurates the integral of a given function, F defined 
  !  in a main program.  The number of lattice spacings (N) a user want use 
  !  must
  !  be specified as the argument of this subroutine.  
  !  The result will be stored
  !  as OUT. The domain of function is [XI,XF]. NOTE:N must be a even number.
  !--------------------------------------------------------
  !  Modefied: 19-aug-1998: this version uses the function 
  !     externally defined.  kurosawa
  !
  !************************************************************************
  SUBROUTINE SIMPSON(nfree, XI,XF,N_in,OUT)
    IMPLICIT NONE
    double precision, intent(in) :: XI,XF
    integer, intent(in)          :: n_in
    integer, intent(in)          :: nfree
    double precision, intent(out) :: OUT
    !
    double precision ::H,SUM,FAC,X
    INTEGER I, n
    !
    !double precision ::  chi_sq
    !
    n = n_in
    IF (N.LT.1) THEN
       PRINT*,'N in Simpson must be greater than 2!'
       STOP
    ELSE IF (MOD(N,2).NE.0) THEN
       N=N+1
    END IF
    H=(XF-XI)/N
    SUM=CHI_SQ_DIST(XI,nfree)          ! cotribution from the initial point
    FAC=2.0d0                		! factor for Simpson's rule
    X=XI
    DO I=1,N-1
       IF (FAC.EQ.2.d0) THEN		! factor switching
          FAC=4.d0
       ELSE
          FAC=2.d0
       END IF
       X=X+H
       SUM=SUM+FAC*CHI_SQ_DIST(X,nfree)! contribution to the integral
    END DO
    SUM=SUM+CHI_SQ_DIST(XF,nfree)	! contribution from the end point
    OUT=SUM*H/3.d0			! the result
    RETURN
  END SUBROUTINE SIMPSON
  


  
  
end module chi_square
