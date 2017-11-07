module qshep3d_mod

use kind_mod
implicit none

contains

subroutine qshep3 ( n, x, y, z, f, nq, nw, nr, lcell, lnext, xyzmin, &
  xyzdel, rmax, rsq, a, ier )
!
!***********************************************************************
!
!! QSHEP3 defines a smooth trivariate interpolant of scattered 3D data.
!
!
!  Discussion:
!
!    This subroutine computes a set of parameters A and RSQ
!    defining a smooth (once continuously differentiable) trivariate 
!    function Q(X,Y,Z) which interpolates data values
!    F at scattered nodes (X,Y,Z).  The interpolant Q may be
!    evaluated at an arbitrary point by function QS3VAL, and
!    its first derivatives are computed by subroutine QS3GRD.
!
!    The interpolation scheme is a modified quadratic Shepard
!    method --
!
!      Q = (W(1)*Q(1)+W(2)*Q(2)+..+W(N)*Q(N))/(W(1)+W(2)+..+W(N))
!
!    for trivariate functions W(K) and Q(K).  The nodal functions are 
!    given by
!
!      Q(K)(X,Y,Z) = 
!          A(1,K) * DX**2 
!        + A(2,K) * DX * DY 
!        + A(3,K) * DY**2
!        + A(4,K) * DX * DZ
!        + A(5,K) * DY * DZ 
!        + A(6,K) * DZ**2
!        + A(7,K) * DX 
!        + A(8,K) * DY 
!        + A(9,K) * DZ 
!        + F(K)
!
!    where DX = (X-X(K)), DY = (Y-Y(K)), and DZ = (Z-Z(K)).
!
!    Thus, Q(K) is a quadratic function which interpolates the
!    data value at node K.  Its coefficients A(*,K) are obtained
!    by a weighted least squares fit to the closest NQ data
!    points with weights similar to W(K).  Note that the radius
!    of influence for the least squares fit is fixed for each
!    K, but varies with K.
!
!    The weights are taken to be
!
!      W(K)(X,Y,Z) = ( (R(K)-D(K))+ / R(K)*D(K) )**2
!
!    where (R(K)-D(K))+ = 0 if R(K) <= D(K), and D(K)(X,Y,Z)
!    is the euclidean distance between (X,Y,Z) and node K.  The
!    radius of influence R(K) varies with K and is chosen so
!    that NW nodes are within the radius.  Note that W(K) is
!    not defined at node (X(K),Y(K),Z(K)), but Q(X,Y,Z) has
!    limit F(K) as (X,Y,Z) approaches (X(K),Y(K),Z(K)).
!
!  Author:
!
!    Robert Renka,
!    University of North Texas,
!    (817) 565-2767.
!
!  Parameters:
!
!    Input, integer N, the number of nodes and associated data values.
!    N >= 10.
!
!    Input, real X(N), Y(N), Z(N), the coordinates of the nodes.
!
!    Input, real F(N), the data values at the nodes.
!
!    Input, integer NQ, the number of data points to be used in the least
!    squares fit for coefficients defining the nodal functions Q(K).  
!    A recommended value is NQ = 17.  9 <= NQ <= MIN ( 40, N-1 ).
!
!    Input, integer NW, the number of nodes within (and defining) the radii
!    of influence R(K) which enter into the weights W(K).  For N sufficiently
!    large, a recommended value is NW = 32.  1 <= NW <= min(40,N-1).
!
!    Input, integer NR, the number of rows, columns, and planes in the cell
!    grid defined in subroutine store3.  a box containing the nodes is
!    partitioned into cells in order to increase search efficiency.  
!    nr = (n/3)**(1/3) is recommended.  nr >= 1.!
!    Output, integer LCELL(NR,NR,NR), nodal indices asso-
!    ciated with cells.  refer to store3.
!
!    Output, integer LNEXT(N), next-node indices.  refer to store3.
!
!    Output, real xyzmin(3), xyzdel(3) = arrays of length 3 containing 
!    minimum nodal coordinates and cell dim-
!    ensions, respectively.  refer to
!    store3.
!
!    Output, rmax = square root of the largest element in rsq --
!    maximum radius r(k).
!
!    Output, real RSQ(N) = array containing the squares of the radii r(k)
!    which enter into the weights w(k).
!
!    Output, real A(9,N), the coefficients for
!    quadratic nodal function q(k) in column k.
!
!    Output, integer IER, error indicator.
!    0, if no errors were encountered.
!    1, if n, nq, nw, or nr is out of range.
!    2, if duplicate nodes were encountered.
!    3, if all nodes are coplanar.
!
!  Local parameters:
!
! av =         root-mean-square distance between k and the
!     nodes in the least squares system (unless
!     additional nodes are introduced for stabil-
!     ity).  the first 6 columns of the matrix
!     are scaled by 1/avsq, the last 3 by 1/av
! avsq =       av*av
! b =         transpose of the augmented regression matrix
! c =         first component of the plane rotation used to
!     zero the lower triangle of b**t -- computed
!     by subroutine givens
! dmin =       minimum of the magnitudes of the diagonal
!     elements of the regression matrix after
!     zeros are introduced below the diagonal
! dtol =       tolerance for detecting an ill-conditioned
!     system.  the system is accepted when dmin
!     >= dtol
! fk =         data value at node k -- f(k)
! i =         index for a, b, npts, xyzmin, xyzmn, xyzdel,
!     and xyzdl
! ib =         do-loop index for back solve
! ierr =       error flag for the call to store3
! ip1 =        i+1
! irm1 =       irow-1
! irow =       row index for b
! j =         index for a and b
! jp1 =        j+1
! k =         nodal function index and column index for a
! lmax =       maximum number of npts elements (must be con-
!     sistent with the dimension statement above)
! lnp =        current length of npts
! neq =        number of equations in the least squares fit
! nn,nnq,nnr = local copies of n, nq, and nr
! nnw =        local copy of nw
! np =         npts element
! npts =       array containing the indices of a sequence of
!     nodes to be used in the least squares fit
!     or to compute rsq.  the nodes are ordered
!     by distance from k and the last element
!     (usually indexed by lnp) is used only to
!     determine rq, or rsq(k) if nw > nq
! nqwmax =     max(nq,nw)
! rq =         radius of influence which enters into the
!     weights for q(k) (see subroutine setup3)
! rs =         squared distance between k and npts(lnp) --
!     used to compute rq and rsq(k)
! rsmx =       maximum rsq element encountered
! rsold =      squared distance between k and npts(lnp-1) --
!     used to compute a relative change in rs
!     between succeeding npts elements
! rtol =       tolerance for detecting a sufficiently large
!     relative change in rs.  if the change is
!     not greater than rtol, the nodes are
!     treated as being the same distance from k
! rws =        current value of rsq(k)
! s =         second component of the plane givens rotation
! sf =         marquardt stabilization factor used to damp
!     out the first 6 solution components (second
!     partials of the quadratic) when the system
!     is ill-conditioned.  as sf increases, the
!     fitting function approaches a linear
! sum2 =        sum of squared euclidean distances between
!     node k and the nodes used in the least
!     squares fit (unless additional nodes are
!     added for stability)
! t =         temporary variable for accumulating a scalar
!     product in the back solve
! xk,yk,zk =   coordinates of node k -- x(k), y(k), z(k)
! xyzdl =      local variables for xyzdel
! xyzmn =      local variables for xyzmin
!
  implicit none
!
  integer n
  integer nr
!
  real(double) a(9,n)
  real(double) av
  real(double) avsq
  real(double) b(10,10)
  real(double) c
  real(double) dmin
  real(double), parameter :: dtol = 0.01E+00
  real(double) f(n)
  real(double) fk
  integer i
  integer ib
  integer ier
  integer ierr
  integer ip1
  integer irm1
  integer irow
  integer j
  integer jp1
  integer k
  integer lcell(nr,nr,nr)
  integer lmax
  integer lnext(n)
  integer lnp
  integer neq
  integer nn
  integer nnq
  integer nnr
  integer nnw
  integer np
  integer npts(40)
  integer nq
  integer nqwmax
  integer nw
  real(double) rmax
  real(double) rq
  real(double) rs
  real(double) rsmx
  real(double) rsold
  real(double) rsq(n)
  real(double), parameter :: rtol = 1.0E-05
  real(double) rws
  real(double) s
  real(double), parameter :: sf = 1.0E+00
  real(double) sum2
  real(double) t
  real(double) x(n)
  real(double) xk
  real(double) xyzdel(3)
  real(double) xyzdl(3)
  real(double) xyzmin(3)
  real(double) xyzmn(3)
  real(double) y(n)
  real(double) yk
  real(double) z(n)
  real(double) zk
!
  nn = n
  nnq = nq
  nnw = nw
  nnr = nr
  nqwmax = max(nnq,nnw)
  lmax = min(40,nn-1)

  if ( 9 > nnq .or.  1 > nnw  .or.  nqwmax > &
    lmax  .or.  nnr < 1 ) then
    ier = 1
    return
  end if
!
!  Create the cell data structure, and initialize RSMX.
!
  call store3 ( nn, x, y, z, nnr, lcell, lnext, xyzmn, xyzdl, ierr )

  if ( ierr /= 0 ) then
    xyzmin(1:3) = xyzmn(1:3)
    xyzdel(1:3) = xyzdl(1:3)
    ier = 3
    return
  end if

  rsmx = 0.0E+00
!
!  Outer loop on node K.
!
  do k = 1, nn

    xk = x(k)
    yk = y(k)
    zk = z(k)
    fk = f(k)
!
!  Mark node K to exclude it from the search for nearest neighbors.
!
    lnext(k) = -lnext(k)
!
!  Initialize for loop on NPTS.
!
    rs = 0.0E+00
    sum2 = 0.0E+00
    rws = 0.0E+00
    rq = 0.0E+00
    lnp = 0
!
!  Compute NPTS, LNP, rws, neq, rq, and avsq.
!
1   continue

    sum2 = sum2 + rs
    if ( lnp == lmax ) go to 3
    lnp = lnp + 1
    rsold = rs
    call getnp3 ( xk, yk, zk, x, y, z, nnr, lcell, lnext, xyzmn, xyzdl, np, rs )
    if ( rs == 0.0E+00 ) go to 21
    npts(lnp) = np
    if ( (rs-rsold)/rs <  rtol ) go to 1
    if ( rws == 0.0E+00 .and.  lnp > nnw ) rws = rs
    if ( rq /= 0.0E+00  .or.  lnp <= nnq ) go to 2
!
!  RQ = 0 (not yet computed) and LNP > nq.  rq =
!  sqrt(rs) is sufficiently large to (strictly) include
!  nq nodes.  The least squares fit will include neq =
!  lnp-1 equations for 9 <= nq <= neq < lmax <= n-1.
!
    neq = lnp - 1
    rq = sqrt(rs)
    avsq = sum2 / dble( neq )
!
!  Bottom of loop -- test for termination.
!
2   continue

    if ( lnp > nqwmax ) go to 4
    go to 1
!
!  All LMAX nodes are included in npts.  rws and/or rq**2 is
!  (arbitrarily) taken to be 10 percent larger than the
!  distance RS to the last node included.
!
3   continue

    if ( rws == 0.0E+00 ) rws = 1.1E+00 * rs

    if ( rq == 0.0E+00 ) then
      neq = lmax
      rq = sqrt ( 1.1 * rs )
      avsq = sum2 / dble( neq )
    end if
!
!  Store RSQ(K), update RSMX if necessary, and compute av.
!
4   continue

    rsq(k) = rws
    if ( rws > rsmx ) rsmx = rws
    av = sqrt(avsq)
!
!  Set up the augmented regression matrix (transposed) as the
!  columns of B, and zero out the lower triangle (upper
!  triangle of b) with givens rotations -- qr decomposition
!  with orthogonal matrix q not stored.
!
    i = 0

5   continue

    i = i + 1
    np = npts(i)
    irow = min0(i,10)
    call setup3 (xk,yk,zk,fk,x(np),y(np),z(np),f(np),av,avsq,rq, b(1,irow))
    if ( i == 1 ) go to 5

    irm1 = irow-1

    do j = 1, irow-1
      jp1 = j + 1
      call givens (b(j,j),b(j,irow),c,s)
      call rotate (10-j,c,s,b(jp1,j),b(jp1,irow))
    end do

    if ( i < neq ) go to 5
!
!  Test the system for ill-conditioning.
!
    dmin = min ( abs(b(1,1)),abs(b(2,2)),abs(b(3,3)), &
               abs(b(4,4)),abs(b(5,5)),abs(b(6,6)), &
               abs(b(7,7)),abs(b(8,8)),abs(b(9,9)) )

    if ( dmin * rq >= dtol ) go to 13
    if ( neq == lmax ) go to 10
!
!  Increase RQ and add another equation to the system to
!  improve the conditioning.  The number of NPTS elements
!  is also increased if necessary.
!
7   continue

    rsold = rs
    neq = neq + 1
    if ( neq == lmax ) go to 9
    if ( neq == lnp ) go to 8
!
!  NEQ < LNP
!
    np = npts(neq+1)
    rs = (x(np)-xk)**2 + (y(np)-yk)**2 + (z(np)-zk)**2
    if ( (rs-rsold)/rs < rtol ) go to 7
    rq = sqrt(rs)
    go to 5
!
!  Add an element to NPTS.
!
8   continue

    lnp = lnp + 1
    call getnp3 (xk,yk,zk,x,y,z,nnr,lcell,lnext,xyzmn,xyzdl, np,rs)
    if ( np == 0 ) go to 21
    npts(lnp) = np
    if ( (rs-rsold)/rs < rtol ) go to 7
    rq = sqrt ( rs )
    go to 5

9   continue

    rq = sqrt ( 1.1E+00 * rs )
    go to 5
!
!  Stabilize the system by damping second partials.  Add
!  multiples of the first six unit vectors to the first
!  six equations.
!
10  continue

    do i = 1, 6

      b(i,10) = sf
      ip1 = i + 1

      do j = ip1, 10
        b(j,10) = 0.0E+00
      end do

      do j = i, 9
        jp1 = j + 1
        call givens (b(j,j),b(j,10),c,s)
        call rotate (10-j,c,s,b(jp1,j),b(jp1,10))
      end do

    end do
!
!  Test the stabilized system for ill-conditioning.
!
    dmin = min ( abs(b(1,1)),abs(b(2,2)),abs(b(3,3)), &
               abs(b(4,4)),abs(b(5,5)),abs(b(6,6)), &
               abs(b(7,7)),abs(b(8,8)),abs(b(9,9)) )

    if ( dmin * rq < dtol ) then
      go to 22
    end if
!
!  Solve the 9 by 9 triangular system for the coefficients
!
13  continue

    do ib = 1, 9
      i = 10-ib
      t = 0.0E+00
      do j = i+1, 9
        t = t + b(j,i)*a(j,k)
      end do
      a(i,k) = (b(10,i)-t)/b(i,i)
    end do
!
!  Scale the coefficients to adjust for the column scaling.
!
    a(1:6,k) = a(1:6,k) / avsq
    a(7,k) = a(7,k)/av
    a(8,k) = a(8,k)/av
    a(9,k) = a(9,k)/av
!
!  Unmark K and the elements of NPTS.
!
    lnext(k) = -lnext(k)
    do i = 1, lnp
      np = npts(i)
      lnext(np) = -lnext(np)
    end do

  end do
!
!  No errors encountered.
!
  xyzmin(1:3) = xyzmn(1:3)
  xyzdel(1:3) = xyzdl(1:3)
  rmax = sqrt ( rsmx )
  ier = 0
  return
!
!  N, NQ, NW, or NR is out of range.
!
!20 continue

  ier = 1
  return
!
!  Duplicate nodes were encountered by GETNP3.
!
21 ier = 2
  return
!
!  No unique solution due to collinear nodes.
!
22 continue

  xyzmin(1:3) = xyzmn(1:3)
  xyzdel(1:3) = xyzdl(1:3)

  ier = 3
end subroutine qshep3
function qs3val ( px, py, pz, n, x, y, z, f, nr, lcell, lnext, xyzmin, &
  xyzdel, rmax, rsq, a )
!
!***********************************************************************
!
!! QS3VAL evaluates the interpolant function Q(X,Y,Z) created by QSHEP3.
!
!
!  Discussion:
!
!    This function returns the value Q(PX,PY,PZ) where Q is
!    the weighted sum of quadratic nodal functions defined in
!    subroutine QSHEP3.  QS3GRD may be called to compute a
!    gradient of Q along with the value, or to test for errors.
!
!    This function should not be called if a nonzero error flag was
!    returned by QSHEP3.
!
!  Author:
!
!    Robert Renka,
!    University of North Texas,
!    (817) 565-2767.
!
!  Parameters:
!
!    Input, real(double) PX, PY, PZ, the point P at which Q is to be evaluated.
!
!    Input, integer N, the number of nodes and data values defining Q.
!    N >= 10.
!
!    Input, real(double) X(N), Y(N), Z(N), F(N), the node coordinates
!    and data values interpolated by Q.
!
!    Input, integer NR, the number of rows, columns and planes in the cell
!    grid.  Refer to STORE3.  NR >= 1.
!
!    Input, integer LCELL(NR,NR,NR), nodal indices associated with cells.  
!    Refer to STORE3.
!
!    Input, integer LNEXT(N), the next-node indices.  Refer to STORE3.
!
!    Input, real(double) XYZMIN(3), XYZDEL(3), the minimum nodal coordinates and 
!    cell dimensions, respectively.  XYZDEL elements must be positive.  
!    Refer to STORE3.
!
!    Input, real(double) RMAX, the square root of the largest element in RSQ,
!    the maximum radius.
!
!    Input, real(double) RSQ(N), the squared radii which enter into the weights 
!    defining Q.
!
!    Input, real(double) A(9,N), the coefficients for the nodal functions defining Q.
!
!    Output, real(double) QS3VAL, the function value Q(PX,PY,PZ) unless N, NR,
!    XYZDEL, or RMAX is invalid, in which case the value 0 is returned.
!
  implicit none
!
  integer n
  integer nr
!
  real(double) a(9,n)
  real(double) delx
  real(double) dely
  real(double) delz
  real(double) dxsq
  real(double) dysq
  real(double) dzsq
  real(double) ds
  real(double) dx
  real(double) dy
  real(double) dz
  real(double) f(n)
  integer i
  integer imax
  integer imin
  integer j
  integer jmax
  integer jmin
  integer k
  integer kmax
  integer kmin
  integer l
  integer lcell(nr,nr,nr)
  integer lnext(n)
  integer lp
  real(double) px
  real(double) py
  real(double) pz
  real(double) qs3val
  real(double) rd
  real(double) rds
  real(double) rmax
  real(double) rs
  real(double) rsq(n)
  real(double) sw
  real(double) swq
  real(double) w
  real(double) x(n)
!  real(double) xmax
  real(double) xmin
  real(double) xp
  real(double) xyzdel(3)
  real(double) xyzmin(3)
  real(double) y(n)
!  real(double) ymax
  real(double) ymin
  real(double) yp
  real(double) z(n)
!  real(double) zmax
  real(double) zmin
  real(double) zp
!
  xp = px
  yp = py
  zp = pz
  xmin = xyzmin(1)
  ymin = xyzmin(2)
  zmin = xyzmin(3)
  dx = xyzdel(1)
  dy = xyzdel(2)
  dz = xyzdel(3)

  if ( n < 10  .or.  nr < 1  .or.  dx <= 0.0 &
         .or.  dy <= 0.0  .or.  dz <= 0.0  .or. &
         rmax < 0.0 ) then
    qs3val = 0.0E+00
    return
  end if
!
!  Set IMIN, imax, jmin, jmax, kmin, and kmax to cell indices
!  defining the range of the search for nodes whose radii
!  include P.  The cells which must be searched are those
!  intersected by (or contained in) a sphere of radius rmax
!  centered at P.
!
  imin = int((xp-xmin-rmax)/dx) + 1
  imin = max ( imin, 1 )

  imax = int((xp-xmin+rmax)/dx) + 1
  imax = min ( imax, nr )

  jmin = int((yp-ymin-rmax)/dy) + 1
  jmin = max ( jmin, 1 )

  jmax = int((yp-ymin+rmax)/dy) + 1
  jmax = min ( jmax, nr )

  kmin = int((zp-zmin-rmax)/dz) + 1
  kmin = max ( kmin, 1 )

  kmax = int((zp-zmin+rmax)/dz) + 1
  kmax = min ( kmax, nr )
!
!  Test for no cells within the sphere of radius RMAX.
!
  if ( imin > imax .or. jmin > jmax .or. kmin > kmax ) then
    qs3val = 0.0E+00
    return
  end if
!
!  Accumulate weight values in SW and weighted nodal function
!  values in SWQ.  The weights are w(l) = ((r-d)+/(r*d))**2
!  for r**2 = rsq(l) and d = distance between P and node L.
!
  sw = 0.0E+00
  swq = 0.0E+00
!
!  Outer loop on cells (i,j,k).
!
  do k = kmin, kmax
    do j = jmin, jmax
      do i = imin, imax

        l = lcell(i,j,k)

        if ( l == 0 ) then
          cycle
        end if
!
!  Inner loop on nodes L.
!
        do

          delx = xp - x(l)
          dely = yp - y(l)
          delz = zp - z(l)

          dxsq = delx*delx
          dysq = dely*dely
          dzsq = delz*delz

          ds = dxsq + dysq + dzsq
          rs = rsq(l)

          if ( ds < rs ) then

            if ( ds == 0.0E+00 ) then
              qs3val = f(l)
              return
            end if

            rds = rs*ds
            rd = sqrt(rds)
            w = (rs+ds-rd-rd)/rds
            sw = sw + w

            swq = swq + w *( a(1,l)*dxsq + a(2,l)*delx*dely + &
              a(3,l)*dysq + a(4,l)*delx*delz + &
              a(5,l)*dely*delz + a(6,l)*dzsq + &
              a(7,l)*delx + a(8,l)*dely + &
              a(9,l)*delz + f(l) )

          end if

          lp = l
          l = lnext(lp)

          if ( l == lp ) then
            exit
          end if

        end do

      end do
    end do
  end do
!
!  SW = 0 iff P is not within the radius R(L) for any node L.
!
  if ( sw == 0.0E+00 ) then
    qs3val = 0.0E+00
  else
    qs3val = swq / sw
  end if

end function qs3val
subroutine qs3grd ( px, py, pz, n, x, y, z, f, nr, lcell, lnext, xyzmin, &
  xyzdel, rmax, rsq, a, q, qx, qy, qz, ier )
!
!***********************************************************************
!
!! QS3GRD computes the value and gradient of the interpolant function.
!
!
!  Discussion:
!
!    This subroutine computes the value and gradient at (PX,PY,PZ) of 
!    the interpolatory function Q defined in subroutine QSHEP3.  
!
!    Q(X,Y,Z) is a weighted sum of quadratic nodal functions.
!
!  Author:
!
!    Robert Renka,
!    University of North Texas,
!    (817) 565-2767.
!
!  Parameters:
!
!    Input, real(double) PX, PY, PZ, the point P at which Q and its partials are 
!    to be evaluated.
!
!    Input, integer N, the number of nodes and data values defining Q.
!    N >= 10.
!
!    Input, real(double) X(N), Y(N), Z(N), F(N), the node coordinates and
!    data values interpolated by Q.
!
!    Input, integer NR, the number of rows, columns and planes in the cell
!    grid.  Refer to STORE3.  NR >= 1.
!
!    Input, integer LCELL(NR,NR,NR), nodal indices associated with cells.  
!    Refer to STORE3.
!
!    Input, integer LNEXT(N), the next-node indices.  Refer to STORE3.
!
!    Input, real(double) XYZMIN(3), XYZDEL(3), the minimum nodal coordinates and 
!    cell dimensions, respectively.  XYZDEL elements must be positive.  
!    Refer to STORE3.
!
!    Input, real(double) RMAX, the square root of the largest element in RSQ,
!    the maximum radius.
!
!    Input, real(double) RSQ(N), the squared radii which enter into the weights 
!    defining Q.
!
!    Input, real(double) A(9,N), the coefficients for the nodal functions defining Q.
!
!    Output, real(double) Q, the value of Q at (PX,PY,PZ) unless IER == 1, in
!    which case no values are returned.
!
!    Output, real(double) QX, QY, QZ, the first partial derivatives of Q at
!    (PX,PY,PZ) unless IER == 1.
!
!    Output, integer IER, error indicator
!    0, if no errors were encountered.
!    1, if N, NR, XYZDEL, or RMAX is invalid.
!    2, if no errors were encountered but (PX.PY.PZ) is not within the 
!       radius R(K) for any node K (and thus Q = QX = QY = QZ = 0).
!
  implicit none
!
  integer n
  integer nr
!
  real(double) a(9,n)
  real(double) delx
  real(double) dely
  real(double) delz
  real(double) ds
  real(double) dx
  real(double) dxsq
  real(double) dy
  real(double) dysq
  real(double) dz
  real(double) dzsq
  real(double) f(n)
  integer i
  integer ier
  integer imax
  integer imin
  integer j
  integer jmax
  integer jmin
  integer k
  integer kmax
  integer kmin
  integer l
  integer lcell(nr,nr,nr)
  integer lnext(n)
  integer lp
  real(double) px
  real(double) py
  real(double) pz
  real(double) q
  real(double) ql
  real(double) qlx
  real(double) qly
  real(double) qlz
  real(double) qx
  real(double) qy
  real(double) qz
  real(double) rd
  real(double) rds
  real(double) rmax
  real(double) rs
  real(double) rsq(n)
  real(double) sw
  real(double) swq
  real(double) swqx
  real(double) swqy
  real(double) swqz
  real(double) sws
  real(double) swx
  real(double) swy
  real(double) swz
  real(double) t
  real(double) w
  real(double) wx
  real(double) wy
  real(double) wz
  real(double) x(n)
!  real(double) xmax
  real(double) xmin
  real(double) xp
  real(double) xyzdel(3)
  real(double) xyzmin(3)
  real(double) y(n)
!  real(double) ymax
  real(double) ymin
  real(double) yp
  real(double) z(n)
!  real(double) zmax
  real(double) zmin
  real(double) zp
!
  xp = px
  yp = py
  zp = pz
  xmin = xyzmin(1)
  ymin = xyzmin(2)
  zmin = xyzmin(3)
  dx = xyzdel(1)
  dy = xyzdel(2)
  dz = xyzdel(3)

  if ( n < 10  .or.  nr < 1  .or.  dx <= 0. &
         .or.  dy <= 0.0E+00  .or.  dz <= 0.0E+00  .or. &
         rmax < 0.0E+00 ) then
    ier = 1
    return
  end if
!
!  Set IMIN, imax, jmin, jmax, kmin, and kmax to cell indices
!  defining the range of the search for nodes whose radii
!  include P.  The cells which must be searched are those
!  intersected by (or contained in) a sphere of radius RMAX
!  centered at P.
!
  imin = int((xp-xmin-rmax)/dx) + 1
  imin = max ( imin, 1 )

  imax = int((xp-xmin+rmax)/dx) + 1
  imax = min ( imax, nr )

  jmin = int((yp-ymin-rmax)/dy) + 1
  jmin = max ( jmin, 1 )

  jmax = int((yp-ymin+rmax)/dy) + 1
  jmax = min ( jmax, nr )

  kmin = int((zp-zmin-rmax)/dz) + 1
  kmin = max ( kmin, 1 )

  kmax = int((zp-zmin+rmax)/dz) + 1
  kmax = min ( kmax, nr )
!
!  Test for no cells within the sphere of radius RMAX.
!
  if ( imin > imax .or. jmin > jmax .or. kmin > kmax ) then
    q = 0.0E+00
    qx = 0.0E+00
    qy = 0.0E+00
    qz = 0.0E+00
    ier = 2
    return
  end if
!
!  Q = swq/sw = sum(w(l)*q(l))/sum(w(l)) where the sum is
!  from l = 1 to N, q(l) is the quadratic nodal function,
!  and w(l) = ((r-d)+/(r*d))**2 for radius r(l) and dist-
!  ance d(l).  thus
!
!    qx = (swqx*sw - swq*swx)/sw**2
!    qy = (swqy*sw - swq*swy)/sw**2
!    qz = (swqz*sw - swq*swz)/sw**2
!
!  where swqx and swx are partial derivatives with respect
!  to X of swq and sw, respectively.  swqy, swy, swqz, and
!  swz are defined similarly.
!
  sw = 0.0E+00
  swx = 0.0E+00
  swy = 0.0E+00
  swz = 0.0E+00
  swq = 0.0E+00
  swqx = 0.0E+00
  swqy = 0.0E+00
  swqz = 0.0E+00
!
!  Outer loop on cells (i,j,k).
!
  do k = kmin, kmax
    do j = jmin, jmax
      do i = imin, imax

        l = lcell(i,j,k)
        if ( l == 0 ) then
          cycle
        end if
!
!  Inner loop on nodes L.
!
        do

          delx = xp - x(l)
          dely = yp - y(l)
          delz = zp - z(l)
          dxsq = delx*delx
          dysq = dely*dely
          dzsq = delz*delz
          ds = dxsq + dysq + dzsq
          rs = rsq(l)

          if ( ds < rs ) then

            if ( ds == 0.0E+00 ) then
              q = f(l)
              qx = a(7,l)
              qy = a(8,l)
              qz = a(9,l)
              ier = 0
              return
            end if

            rds = rs*ds
            rd = sqrt(rds)
            w = (rs+ds-rd-rd)/rds
            t = 2.0E+00 *(rd-rs)/(ds*rds)
            wx = delx*t
            wy = dely*t
            wz = delz*t

            qlx = 2.0E+00 *a(1,l)*delx + a(2,l)*dely + a(4,l)*delz
            qly = a(2,l)*delx + 2.0E+00 * a(3,l)*dely + a(5,l)*delz
            qlz = a(4,l)*delx + a(5,l)*dely + 2.0E+00 * a(6,l)*delz
            ql = (qlx*delx + qly*dely + qlz*delz) / 2.0E+00 + &
              a(7,l)*delx + a(8,l)*dely + a(9,l)*delz + f(l)
            qlx = qlx + a(7,l)
            qly = qly + a(8,l)
            qlz = qlz + a(9,l)

            sw = sw + w
            swx = swx + wx
            swy = swy + wy
            swz = swz + wz
            swq = swq + w*ql
            swqx = swqx + wx*ql + w*qlx
            swqy = swqy + wy*ql + w*qly
            swqz = swqz + wz*ql + w*qlz

          end if

          lp = l
          l = lnext(lp)

          if ( l == lp ) then
            exit
          end if

        end do

      end do
    end do
  end do
!
!  SW = 0 iff P is not within the radius R(L) for any node L.
!
  if ( sw /= 0.0E+00 ) then

    q = swq/sw
    sws = sw*sw
    qx = (swqx*sw - swq*swx)/sws
    qy = (swqy*sw - swq*swy)/sws
    qz = (swqz*sw - swq*swz)/sws
    ier = 0
!
!  No cells contain a point within RMAX of P, or
!  SW = 0 and thus DS >= RSQ(L) for all L.
!
  else

    q = 0.0E+00
    qx = 0.0E+00
    qy = 0.0E+00
    qz = 0.0E+00
    ier = 2

  end if

end subroutine qs3grd

subroutine getnp3 ( px, py, pz, x, y, z, nr, lcell, lnext, xyzmin, &
  xyzdel, np, dsq )
!
!***********************************************************************
!
!! GETNP3 finds the closest node to a given point.
!
!
!  Discussion:
!
!    Given a set of N nodes and the data structure defined in
!    subroutine STORE3, this subroutine uses the cell method to
!    find the closest unmarked node NP to a specified point P.
!
!    NP is then marked by setting LNEXT(NP) to -LNEXT(NP).  (A
!    node is marked if and only if the corresponding lnext element 
!    is negative.  The absolute values of LNEXT elements,
!    however, must be preserved.)  Thus, the closest M nodes to
!    P may be determined by a sequence of M calls to this routine.  
!    Note that if the nearest neighbor to node K is to
!    be determined (PX = X(K), PY = Y(K), and PZ = Z(K)), then
!    K should be marked before the call to this routine.
!
!    The search is begun in the cell containing (or closest
!    to) P and proceeds outward in box-shaped layers until all
!    cells which contain points within distance R of P have
!    been searched, where R is the distance from P to the first
!    unmarked node encountered (infinite if no unmarked nodes
!    are present).
!
!  Author:
!
!    Robert Renka,
!    University of North Texas,
!    (817) 565-2767.
!
!  Parameters:
!
!    Input, real(double) PX, PY, PZ, the coordinates of the point P whose nearest
!    unmarked neighbor is to be found.
!
!    Input, real(double) X(N), Y(N), Z(N), the coordinates of the nodes.
!
!    Input, integer NR, the number of rows, columns, and planes in the cell
!    grid.  NR >= 1.
!
!    Input, integer LCELL(NR,NR,NR), nodal indices associated with cells.
!
!    Input/output, integer LNEXT(N), next-node indices (or their negatives).
!
!    Input, real(double) XYZMIN(3), XYZDEL(3), minimum nodal coordinates and cell 
!    dimensions, respectively.  XYZEL elements must be positive.
!
!    Output, integer NP, index of the nearest unmarked node to P, or 0 
!    if all nodes are marked or NR < 1 or an element of XYZDEL is not 
!    positive.  LNEXT(NP) < 0 if NP /= 0.
!
!    Output, real(double) DSQ, squared euclidean distance between P and node
!    NP, or 0 if NP = 0.
!
  implicit none
!
  integer nr
!
  real(double) delx
  real(double) dely
  real(double) delz
  real(double) dsq
  real(double) dx
  real(double) dy
  real(double) dz
  logical first
  integer i
  integer i0
  integer i1
  integer i2
  integer imax
  integer imin
  integer j
  integer j0
  integer j1
  integer j2
  integer jmax
  integer jmin
  integer k
  integer k0
  integer k1
  integer k2
  integer kmax
  integer kmin
  integer l
  integer lcell(nr,nr,nr)
  integer lmin
  integer ln
  integer lnext(*)
  integer np
  real(double) px
  real(double) py
  real(double) pz
  real(double) r
  real(double) rsmin
  real(double) rsq
  real(double) x(*)
  real(double) xp
  real(double) xyzdel(3)
  real(double) xyzmin(3)
  real(double) y(*)
  real(double) yp
  real(double) z(*)
  real(double) zp
!
  xp = px
  yp = py
  zp = pz
  dx = xyzdel(1)
  dy = xyzdel(2)
  dz = xyzdel(3)
!
!  Test for invalid input parameters.
!
  if ( nr < 1 .or. dx <= 0.0E+00 .or. dy <= 0.0E+00 .or. dz <= 0.0E+00 ) then
    np = 0
    dsq = 0.0E+00
    return
  end if
!
!  Initialize parameters --
!
!  first = true iff the first unmarked node has yet to be encountered,
!  imin,...,kmax = cell indices defining the range of the search,
!  delx,dely,delz = px-xyzmin(1), py-xyzmin(2), and pz-xyzmin(3),
!  i0,j0,k0 = cell containing or closest to p,
!  i1,...,k2 = cell indices of the layer whose intersection
!  with the range defined by imin,...,kmax is
!  currently being searched.
!
  first = .true.
  imin = 1
  imax = nr
  jmin = 1
  jmax = nr
  kmin = 1
  kmax = nr
  delx = xp - xyzmin(1)
  dely = yp - xyzmin(2)
  delz = zp - xyzmin(3)
  i0 = int(delx/dx) + 1
  if ( i0 < 1 ) i0 = 1
  if ( i0 > nr ) i0 = nr
  j0 = int ( dely / dy ) + 1
  if ( j0 < 1 ) j0 = 1
  if ( j0 > nr ) j0 = nr
  k0 = int(delz/dz) + 1
  if ( k0 < 1 ) k0 = 1
  if ( k0 > nr ) k0 = nr
  i1 = i0
  i2 = i0
  j1 = j0
  j2 = j0
  k1 = k0
  k2 = k0
!
!  Outer loop on layers, inner loop on layer cells, excluding
!  those outside the range (imin,imax) x (jmin,jmax) x (kmin,kmax).
!
1 continue

  do 7 k = k1, k2

    if ( k > kmax ) go to 8
    if ( k < kmin ) go to 7

    do 6 j = j1, j2

      if ( j > jmax ) go to 7
      if ( j < jmin ) go to 6

      do 5 i = i1, i2

        if ( i > imax ) go to 6
        if ( i < imin ) go to 5
        if ( k /= k1  .and.  k /= k2  .and. j /= j1  .and.  &
          j /= j2  .and. i /= i1  .and.  i /= i2) go to 5
!
!  Search cell (i,j,k) for unmarked nodes L.
!
        l = lcell(i,j,k)
        if ( l == 0 ) go to 5
!
!  Loop on nodes in cell (i,j,k).
!
2       continue

        ln = lnext(l)
        if ( ln < 0 ) go to 4
!
!  Node L is not marked.
!
        rsq = (x(l)-xp)**2 + (y(l)-yp)**2 + (z(l)-zp)**2
        if ( .not. first ) go to 3
!
!  Node L is the first unmarked neighbor of P encountered.
!  initialize LMIN to the current candidate for NP, and
!  rsmin to the squared distance from P to lmin.  imin,
!  imax, jmin, jmax, kmin, and kmax are updated to define
!  the smallest rectangle containing a sphere of radius
!  r = sqrt(rsmin) centered at P, and contained in (1,nr)
!  x (1,nr) x (1,nr) (except that, if P is outside the
!  box defined by the nodes, it is possible that imin
!  > nr or imax < 1, etc.).  FIRST is reset to
!  false.
!
  lmin = l
  rsmin = rsq
  r = sqrt(rsmin)
  imin = int((delx-r)/dx) + 1
  if ( imin < 1 ) imin = 1
  imax = int((delx+r)/dx) + 1
  if ( imax > nr ) imax = nr
  jmin = int((dely-r)/dy) + 1
  if ( jmin < 1 ) jmin = 1
  jmax = int((dely+r)/dy) + 1
  if ( jmax > nr ) jmax = nr
  kmin = int((delz-r)/dz) + 1
  if ( kmin < 1 ) kmin = 1
  kmax = int((delz+r)/dz) + 1
  if ( kmax > nr ) kmax = nr
  first = .false.
  go to 4
!
!  Test for node L closer than LMIN to P.
!
3      continue
!
!  Update LMIN and RSMIN.
!
  if ( rsq < rsmin ) then
    lmin = l
    rsmin = rsq
  end if
!
!  Test for termination of loop on nodes in cell (i,j,k).
!
4      continue

  if ( abs(ln) == l ) go to 5
  l = abs(ln)
  go to 2

5      continue
6    continue
7  continue
!
! Test for termination of loop on cell layers.
!
8 continue

  if ( i1 <= imin  .and.  i2 >= imax  .and. &
         j1 <= jmin  .and.  j2 >= jmax  .and. &
         k1 <= kmin  .and.  k2 >= kmax ) go to 9
  i1 = i1 - 1
  i2 = i2 + 1
  j1 = j1 - 1
  j2 = j2 + 1
  k1 = k1 - 1
  k2 = k2 + 1
  go to 1
!
!  Unless no unmarked nodes were encountered, LMIN is the
!  closest unmarked node to P.
!
9 continue

  if ( .not. first ) then
    np = lmin
    dsq = rsmin
    lnext(lmin) = -lnext(lmin)
  else
    np = 0
    dsq = 0.0E+00
  end if

end subroutine getnp3

subroutine givens ( a, b, c, s )
!
!***********************************************************************
!
!! GIVENS constructs a Givens plane rotation.
!
!
!  Discussion:
!
!    This routine constructs the Givens plane rotation
!
!          ( c  s)
!      g = (     ) 
!          (-s  c)
!
!    where c*c + s*s = 1, which zeros the second entry of the 2-vector 
!    (a b)-transpose.  A call to GIVENS is normally followed by a call 
!    to ROTATE which applies the transformation to a 2 by N matrix.  
!    This routine was taken from LINPACK.
!
!  Author:
!
!    Robert Renka,
!    University of North Texas,
!    (817) 565-2767.
!
!  Parameters:
!
! on input --
!
!  a,b = components of the 2-vector to be rotated.
!
! on output --
!
!  a = value overwritten by r = +/-sqrt(a*a + b*b)
!
!  b = value overwritten by a value z which allows c
!      and s to be recovered as follows --
!        c = sqrt(1-z*z), s=z     if abs(z) <= 1.
!        c = 1/z, s = sqrt(1-c*c) if abs(z) > 1.
!
!  c = +/-(a/r)
!
!  s = +/-(b/r)
!
!  Local parameters:
!
!    aa,bb = local copies of a and b
!    r =    c*a + s*b = +/-sqrt(a*a+b*b)
!    u,v =   variables used to scale a and b for computing r
!
  implicit none
!
  real(double) a
  real(double) aa
  real(double) b
  real(double) bb
  real(double) c
  real(double) r
  real(double) s
  real(double) u
  real(double) v
!
  aa = a
  bb = b
!
!  Abs(a) > abs(b)
!  Note that R has the sign of A, C > 0, and S has sign(a)*sign(b).
!
  if ( abs ( aa ) > abs ( bb ) ) then
    u = aa + aa
    v = bb/u
    r = sqrt ( 0.25E+00 + v*v ) * u
    c = aa/r
    s = v * (c + c)
    b = s
    a = r
    return
  end if
!
!  abs(a) <= abs(b)
!
  if ( bb == 0.0E+00 ) go to 2
  u = bb + bb
  v = aa/u
!
!  Store R in A.
!
  a = sqrt ( 0.25E+00 + v*v ) * u
  s = bb / a
  c = v * (s + s)
!
!  Note that R has the sign of b, s > 0, and c has sign(a)*sign(b).
!
  b = 1.0E+00
  if ( c /= 0.0E+00 ) then
    b = 1.0E+00 / c
  end if

  return
!
!  A = B = 0.
!
2 continue

  c = 1.0E+00
  s = 0.0E+00

end subroutine givens
subroutine rotate ( n, c, s, x, y )
!
!***********************************************************************
!
!! ROTATE applies a Givens rotation to two vectors.
!
!
!  Discussion:
!
!   This routine applies the Givens rotation
!
!     ( c  s)
!     (-s  c)
!
!   to the 2 by n matrix
!
!    (x(1) ... x(n))
!    (y(1) ... y(n))
!
!  Author:
!
!    Robert Renka,
!    University of North Texas,
!    (817) 565-2767.
!
!  Parameters:
!
! on input --
!
!    Input, integer N, the number of columns to be rotated.
!
!  c, s = elements of the givens rotation.  These may be
!        determined by subroutine GIVENS.
!
!  x, y = arrays of length >= n containing the vectors
!        to be rotated.
!
! parameters n, c, and s are not altered by this routine.
!
! on output --
!
!  x,y = rotated vectors.
!
  implicit none
!
  integer n
!
  real(double) c
  integer i
  real(double) s
  real(double) x(n)
  real(double) xi
  real(double) y(n)
  real(double) yi
!
  if ( c == 1.0E+00 .and. s == 0.0E+00 ) then
    return
  end if

  do i = 1, n
    xi = x(i)
    yi = y(i)
    x(i) = c*xi + s*yi
    y(i) = -s*xi + c*yi
  end do

end subroutine rotate
subroutine setup3 ( xk, yk, zk, fk, xi, yi, zi, fi, s1, s2, r, row )
!
!***********************************************************************
!
!! SETUP3 sets up the weighted least-squares fit of the data.
!
!
!  Discussion:
!
!    This routine sets up the I-th row of an augmented regression matrix 
!    for a weighted least-squares fit of a quadratic function Q(X,Y,Z) 
!    to a set of data values F, where Q(XK,YK,ZK) = FK.  
!
!    The first 6 columns (quadratic terms) are scaled by 1/S2, and columns 
!    7, 8, and 9 (linear terms) are scaled by 1/S1.  The weight is 
!    (R-D)/(R*D) if R > D, and 0 if R <= D, where D is the distance 
!    between nodes I and K.
!
!  Author:
!
!    Robert Renka,
!    University of North Texas,
!    (817) 565-2767.
!
!  Parameters:
!
!    Input, real(double) XK, YK, ZK, FK = coordinates and data value at node K
!    (interpolated by q).
!
!    Input, real(double) XI, YI, ZI, FI = coordinates and data value at node I.
!
!    Input, real(double) S1, S2 = reciprocals of the scale factors.
!
!    Input, real(double) R = radius of influence about node K defining the
!    weight.
!
!    Output, real(double) ROW(10), a row of the augmented
!    regression matrix.
!
!  Local parameters:
!
!    d =   distance between nodes k and i
!    w =   weight associated with the row
!    w1 =   w/s1
!    w2 =   w/s2
!
  implicit none
!
  real(double) d
  real(double) dx
  real(double) dxsq
  real(double) dy
  real(double) dysq
  real(double) dz
  real(double) dzsq
  real(double) fi
  real(double) fk
!  integer i
  real(double) r
  real(double) row(10)
  real(double) s1
  real(double) s2
  real(double) w
  real(double) w1
  real(double) w2
  real(double) xi
  real(double) xk
  real(double) yi
  real(double) yk
  real(double) zi
  real(double) zk
!
  dx = xi - xk
  dy = yi - yk
  dz = zi - zk
  dxsq = dx*dx
  dysq = dy*dy
  dzsq = dz*dz
  d = sqrt ( dxsq + dysq + dzsq )

  if ( d <= 0.0E+00 .or. d >= r ) then
    row(1:10) = 0.0E+00
    return
  end if

  w = ( r - d ) / r / d
  w1 = w/s1
  w2 = w/s2

  row(1) = dxsq*w2
  row(2) = dx*dy*w2
  row(3) = dysq*w2
  row(4) = dx*dz*w2
  row(5) = dy*dz*w2
  row(6) = dzsq*w2
  row(7) = dx*w1
  row(8) = dy*w1
  row(9) = dz*w1
  row(10) = (fi - fk)*w

end subroutine setup3
subroutine store3 ( n, x, y, z, nr, lcell, lnext, xyzmin, xyzdel, ier )
!
!***********************************************************************
!
!! STORE3 sets up a data structure for N scattered nodes in 3D.
!
!
!  Discussion:
!
!    Given a set of N arbitrarily distributed nodes in three-space, 
!    this subroutine creates a data structure for a cell-based method of 
!    solving closest-point problems.  The smallest box containing the nodes 
!    is partitioned into an NR by NR by NR uniform grid of cells, and 
!    nodes are associated with cells.  In particular, the data structure
!    stores the indices of the nodes contained in each cell.  For a 
!    uniform random distribution of nodes, the nearest node to an 
!    arbitrary point can be determined in constant expected time.
!
!  Author:
!
!    Robert Renka,
!    University of North Texas,
!    (817) 565-2767.
!
!  Parameters:
!
!    Input, integer N, the number of nodes.  N >= 2.
!
!    Input, real(double) X(N), Y(N), Z(N), the coordinates of the nodes.
!
!    Input, integer NR, the number of rows, columns, and planes in the
!    grid.  The cell density (average number of nodes per cell) is 
!    D = N/(NR**3).  A recommended value, based on empirical evidence, 
!    is D = 3, so NR = (N/3)**(1/3).  NR >= 1.
!
!    Output, integer LCELL(NR,NR,NR), a cell array such that
!    lcell(i,j,k) contains the index (for x, y,
!    and z) of the first node (node with smallest
!    index) in cell (i,j,k), or lcell(i,j,k) = 0
!    if no nodes are contained in the cell.  the
!    corner of cell (i,j,k) which is farthest
!    from the box corner defined by xyzmin has
!    coordinates (xmin+i*dx,ymin+j*dy,zmin+k*dz),
!    where (xmin,ymin,zmin) are the elements of
!    xyzmin.  lcell is not defined if ier /= 0.
!
!    Output, integer LNEXT(N),  next-node indices such that
!    lnext(l) contains the index of the next node
!    in the cell which contains node l, or
!    lnext(l) = l if l is the last node in the
!    cell for l = 1,...,n.  (the nodes contained
!    in a cell are ordered by their indices.)
!    if, for example, cell (i,j,k) contains nodes
!    2, 3, and 5 (and no others), then
!    lcell(i,j,k) = 2, lnext(2) = 3, lnext(3) =
!    5, and lnext(5) = 5.  lnext is not defined
!    if ier /= 0.
!
!    Output, real(double) XYZMIN(3), the minimum
!    nodal coordinates xmin, ymin, and zmin (in
!    that order) unless ier = 1.  the opposite
!    corner of the box defined by the nodes is
!    (xmin+nr*dx,ymin+nr*dy,zmin+nr*dz).
!
!    Output, real(double) XYZDEL(3), the dimensions
!    of the cells unless ier = 1.  xyzdel(1) =
!    (xmax-xmin)/nr, xyzdel(2) = (ymax-ymin)/nr,
!    and xyzdel(3) = (zmax-zmin)/nr, where xmin,
!    xmax, ymin, ymax, zmin, and zmax are the
!    extrema of x, y, and z.
!
!    Output, integer IER,  = error indicator --
!    0, if no errors were encountered.
!    1, if n < 2 or nr < 1.
!    2, if a component of xyzdel is not positive.
!
  implicit none
!
  integer n
  integer nr
!
  real(double) delx
  real(double) dely
  real(double) delz
  integer i
  integer ier
  integer j
  integer k
  integer l
  integer lb
  integer lcell(nr,nr,nr)
  integer ll
  integer lnext(n)
  integer nn
  integer nnr
  integer np1
  real(double) x(n)
  real(double) xmn
  real(double) xmx
  real(double) xyzdel(3)
  real(double) xyzmin(3)
  real(double) y(n)
  real(double) ymx
  real(double) ymn
  real(double) z(n)
  real(double) zmx
  real(double) zmn
!
  ier = 0
  nn = n
  nnr = nr

  if ( nn < 2 .or. nnr < 1 ) then
    ier = 1
    return
  end if
!
!  Compute the dimensions of the box containing the nodes.
!
  xmn = minval ( x(1:nn) )
  xmx = maxval ( x(1:nn) )
  ymn = minval ( y(1:nn) )
  ymx = maxval ( y(1:nn) )
  zmn = minval ( z(1:nn) )
  zmx = maxval ( z(1:nn) )

  xyzmin(1) = xmn
  xyzmin(2) = ymn
  xyzmin(3) = zmn
!
!  Compute cell dimensions and test for zero area.
!
  delx = (xmx-xmn)/dble(nnr)
  dely = (ymx-ymn)/dble(nnr)
  delz = (zmx-zmn)/dble(nnr)

  xyzdel(1) = delx
  xyzdel(2) = dely
  xyzdel(3) = delz

  if ( delx == 0.0E+00 .or. dely == 0.0E+00 .or. delz == 0.0E+00 ) then
    ier = 2
    return
  end if
!
!  Initialize LCELL.
!
  lcell(1:nnr,1:nnr,1:nnr) = 0
!
!  Loop on nodes, storing indices in LCELL and LNEXT.
!
  np1 = nn + 1

  do ll = 1, nn
    lb = np1 - ll
    i = int((x(lb)-xmn)/delx) + 1
    if ( i > nnr ) i = nnr
    j = int((y(lb)-ymn)/dely) + 1
    if ( j > nnr ) j = nnr
    k = int((z(lb)-zmn)/delz) + 1
    if ( k > nnr ) k = nnr
    l = lcell(i,j,k)
    lnext(lb) = l
    if ( l == 0 ) lnext(lb) = lb
    lcell(i,j,k) = lb
  end do

end subroutine store3
subroutine timestamp ( )
!
!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none
!
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

end subroutine timestamp

end module qshep3d_mod
