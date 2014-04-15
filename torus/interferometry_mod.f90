module interferometry_mod

  use kind_mod
  use constants_mod
  use random_mod
  use messages_mod
  implicit none

  public

  contains

    subroutine testDFT()
      real(double),pointer :: image(:,:)
      real(double), pointer :: xAxis(:), yAxis(:)
      integer :: nx, ny, i, j

      nx = 100
      ny = 100
      allocate(image(1:nx,1:ny), xAxis(1:nx), yAxis(1:ny))
      do i = 1, nx
         xAxis(i) = -1.d0 + 2.d0*dble(i-1)/dble(nx-1)
         yAxis(i) = -1.d0 + 2.d0*dble(i-1)/dble(nx-1)
      enddo
      image = 0.d0

      do i = 1, nx
         do j = 1, ny
            if ((xAxis(i)**2 + yAxis(j)**2) < 0.9d0) image(i,j) = 1.d0
            if ((xAxis(i)**2 + yAxis(j)**2) < 0.5d0) image(i,j) = 0.d0
         enddo
      enddo
      xAxis = xAxis * autocm/(140.d0*pctocm)
      yAxis = yAxis * autocm/(140.d0*pctocm)

!      call visibilityCurve(filename, image, nx, ny, xAxis, yAxis, 2.2d-6, 300.d0)

    end subroutine testDFT


    subroutine discreteFourierTransform(image, nx, ny, xAxis, yAxis, dft, uAxis, vAxis)
      real(kind=double) :: image(:,:)
      integer :: nx, ny
      real(kind=double) :: xAxis(:), yAxis(:)
      complex(double), pointer :: dft(:,:)
      real(double),pointer :: uAxis(:), vAxis(:)
      integer :: i, j, m, n

      allocate(uAxis(1:nx), vAxis(1:nx), dft(1:nx, 1:ny))

      uAxis = xAxis*25.d0
      vAxis = yAxis*25.d0
      dft = 0.d0
      do i = 1, nx
         do j = 1, ny

            do m = 1, nx
               do n = 1, ny
                  dft(i,j) = dft(i,j) + &
                       cmplx(image(m,n), 0.d0,kind=db) &
                       * exp(cmplx(-twoPi,0.d0,kind=db) * cmplx(0.d0, 1.d0,kind=db) &
                             * cmplx((xAxis(m)*uAxis(i) + yAxis(n)*vAxis(j)),0.d0,kind=db))
               enddo
            enddo
         enddo
      enddo
      dft = dft / SUM(image)
    end subroutine discreteFourierTransform

    function visibility(image, nx, ny, xAxis, yAxis, u, v) result(nu)
      real(kind=double) :: image(:,:)
      integer :: nx, ny
      real(kind=double) :: xAxis(:), yAxis(:)
      real(double) :: u, v, nu
      integer :: m, n

      nu = 0.d0
      do m = 1, nx
         do n = 1, ny
            nu = nu + &
            dble(cmplx(image(m,n), 0.d0,kind=db) &
                 * exp(cmplx(-twoPi,0.d0,kind=db) * cmplx(0.d0, 1.d0,kind=db) * cmplx((xAxis(m)*u + yAxis(n)*v),0.d0,kind=db)))
         enddo
      enddo
      nu = nu / SUM(image)
    end function visibility

    subroutine visibilityCurve(vfile, image, nx, ny, xAxis, yAxis, wavelength, baseline)
      real(double) :: image(:,:), xaxis(:), yAxis(:), wavelength, baseLine, thisBaseline
      integer :: nx, ny
      real(double) :: nu, u, v, ang
      integer :: i, j
      character(len=*) :: vFile
      if (writeoutput) then

         open(34, file=vfile, status="unknown", form="formatted")
         do i = 1, 100
            thisBaseline = baseline *  real(i-1)/99.d0
            
            nu = 0.d0
            do j = 1, 100
               call randomNumberGenerator(getdouble = ang)
               ang = ang * twoPi
               u = thisBaseline*cos(ang) / wavelength
               v = thisBaseline*sin(ang) / wavelength
               nu = nu + visibility(image, nx, ny, xAxis, yAxis, u, v)
            enddo
            nu = nu / 100.d0
            write(34,*) thisBaseline,nu**2
         enddo
         close(34)
      endif
    end subroutine visibilityCurve


end module interferometry_mod
