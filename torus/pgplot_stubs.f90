! pgplot stubs for compilation of torus without
! pgplot - th 12/12/05

subroutine pgctab(x1, x2, x3, x4,i, x5, x6)
  integer :: i
  real :: x1, x2, x3, x4, x5, x6
end subroutine pgctab

subroutine pglab(s1, s2, s3)
  character(len=*) :: s1, s2, s3
end subroutine pglab

subroutine pglabel(s1, s2, s3)
  character(len=*) :: s1, s2, s3
end subroutine pglabel

subroutine pgqwin(x1,x2,x3,x4)
  real :: x1, x2, x3, x4
end subroutine pgqwin

subroutine pgqvp(x1,x2,x3,x4)
  real :: x1, x2, x3, x4
end subroutine pgqvp

subroutine pgsvp(x1,x2,x3,x4)
  real :: x1, x2, x3, x4
end subroutine pgsvp

subroutine pgswin(x1,x2,x3,x4)
  real :: x1, x2, x3, x4
end subroutine pgswin

subroutine pgsls(i)
  integer :: i
end subroutine pgsls

subroutine pgslw(f)
  real :: f
end subroutine pgslw

subroutine pgline(n,x,y)
  integer :: n
  real :: x(n), y(n)
end subroutine pgline

subroutine pgmtxt(s1, f, g, h, s2)
  character(len=*) :: s1, s2
  real :: f, g, h
end subroutine pgmtxt

subroutine pgcont(x1, i, j, k, l, m, n, f, n1, tr)
  integer :: i, j, k, l, m, n, n1
  real :: f(*), tr(*)
  real :: x1(i,j)
end subroutine pgcont
  

subroutine pgmove(x, y)
  real :: x, y
end subroutine pgmove

subroutine pgdraw(x, y)
  real :: x, y
end subroutine pgdraw

subroutine pgend()
end subroutine pgend

subroutine pgbegin(i, fname, j, k)
  integer :: i, j, k
  character(len=*) :: fname
end subroutine pgbegin

function pgbeg(i, fname, j, k) result (n)
  integer :: i, j, k, n
  character(len=*) :: fname
  n = 1
end function pgbeg

subroutine pgpaper(x, y)
  real :: x, y
end subroutine pgpaper

subroutine pgsci(i)
  integer :: i
end subroutine pgsci

subroutine pgwnad(x1, x2, x3, x4)
  real :: x1,x2,x3,x4
end subroutine pgwnad

subroutine pgsitf(i)
  integer :: i
end subroutine pgsitf

function pgopen(s1) result (i)
  integer :: i
  character(len=*) :: s1
  i = 0
end function pgopen

function pgcurs(x, y, s1) result (i)
  integer :: i
  real :: x, y
  character(len=*) :: s1
  i = 0
end function pgcurs

subroutine pgvsiz(x1, x2, x3, x4)
  real :: x1, x2, x3, x4
end subroutine pgvsiz


subroutine pgwedg(s1, x1, x2, x3, x4, s2)
  character(len=*) :: s1, s2
  real :: x1, x2, x3, x4
end subroutine pgwedg


subroutine pgqinf(s1, f, i)
  character(len=*) :: s1
  real :: f
  integer ::i 
end subroutine pgqinf

subroutine pgpt1(x, y, i)
  integer :: i
  real :: x, y
end subroutine pgpt1

subroutine pgvect(x, y, i, j, k, l, m, n, f, n1, tr, f2)
  real :: x(*), y(*)
  integer :: i, j, k, l, m, n, n1
  real :: f, f2
end subroutine pgvect

subroutine pgpoint(i, x, y, j)
  integer :: i, j
  real :: x, y
end subroutine pgpoint

subroutine pgqcir(i, j)
  integer :: i, j
end subroutine pgqcir

subroutine pgscir(i, j)
  integer :: i, j
end subroutine pgscir

subroutine pgbin(n, x, y, j)
  integer :: n
  real :: x(*), y(*)
  logical :: j
end subroutine pgbin

subroutine pgvport(x1, x2, x3, x4)
  real :: x1, x2, x3, x4
end subroutine pgvport

subroutine pgbox(s1, x1, i, s2, x2, j)
  integer :: i,j
  character(len=*) :: s1, s2
  real :: x1, x2
end subroutine pgbox

subroutine pgenv(x1, x2, x3, x4, i, j)
  integer :: i, j
  real :: x1, x2, x3, x4
end subroutine pgenv

subroutine pgpage()
end subroutine pgpage

subroutine pgmtext(s1,x,y,z,s2)
  character(len=*) :: s1, s2
  real :: x, y, z
end subroutine pgmtext

subroutine pgsfs(i)
  integer :: i
end subroutine pgsfs

subroutine pgsch(f)
  real :: f
end subroutine pgsch

subroutine pgask(f)
  logical :: f
end subroutine pgask

subroutine pgimag(f, i, j, k, l, m, f1, f2, tr)
  integer :: i, j, k, l, m
  real :: f(i,j),f1,f2
  real :: tr(*)
end subroutine pgimag

subroutine pgrect(x1,x2,x3,x4)
  real :: x1, x2, x3, x4
end subroutine pgrect

subroutine pgqcol(i, j)
  integer :: i, j
  i = 100
  j = 1000
end subroutine pgqcol

subroutine pgtext(x, y, s1)
  integer :: x, y
  character(len=*) :: s1
end subroutine pgtext

subroutine pgwindow(x1, x2, x3, x4)
  real :: x1, x2, x3, x4
end subroutine pgwindow

