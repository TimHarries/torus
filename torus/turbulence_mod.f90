module turbulence_mod
  use kind_mod
  use messages_mod
  use vector_mod
  use constants_mod
  use random_mod
  use utils_mod
  implicit none

contains

  subroutine createBox(box, n)
    integer :: n, i, j, k, ik
    type(VECTOR) :: box(n,n,n)
    real(double) :: kMod, r, phase
    type(VECTOR) :: kVec, xVec,rHat, meanVec
    type(VECTOR), allocatable :: tBox(:,:,:)
    real(double) :: power
    allocate(tBox(1:n,1:n,1:n))
    tBox = VECTOR(0.d0, 0.d0, 0.d0)

    box = VECTOR(0.d0, 0.d0, 0.d0)
    do ik = 1, n

       kmod = (twoPi/dble(ik))
       kVec = kmod*randomUnitVector()
       
       call randomNumberGenerator(getdouble=r)
       phase = r * twoPi
       do i = 1, n
          do j = 1, n
             do k = 1, n

                xVec = VECTOR(dble(i-1)/dble(n-1), dble(j-1)/dble(n-1), dble(k-1)/dble(n-1))
                rHat = kVec
                call normalize(rHat)
                r = 1.d-9 * gasdev() / kmod**4
                power = r 

                box(i,j,k) = box(i,j,k) + rhat*(power*cos((kVec.dot.xVec)+phase))
             enddo
          enddo
       enddo
    enddo

       do i = 1, n
          do j = 1, n
             do k = 1, n

                tbox(i,j,k) = curl(box,n,i,j,k)
             enddo
          enddo
       enddo
       meanVec = VECTOR(0.d0, 0.d0, 0.d0)
       do i = 1, n
          do j = 1, n
             do k = 1, n

                box(i,j,k) = tbox(i,j,k)
                meanVec = meanVec + box(i,j,k)
             enddo
          enddo
       enddo
       meanVec = (1.d0/dble(n**3))*meanVec 
       do i = 1, n
          do j = 1, n
             do k = 1, n
                box(i,j,k) = box(i,j,k) - meanVec
             enddo
          enddo
       enddo
       deallocate(tbox)
          

  end subroutine createBox

  type(VECTOR) function curl(box,n,i,j,k)
    type(VECTOR) :: box(:,:,:)
    integer :: n, i, j, k, ip1, im1, jp1, jm1, kp1, km1
    real(double) :: dFzdy, dFydz, dFxdz, dFzdx, dFydx, dFxdy

    ip1 = i + 1
    if (ip1 > n) then
       curl = VECTOR(0.d0,0.d0,0.d0)
       goto 666
    endif
    im1 = i - 1
    if (im1 < 1) then
       curl = VECTOR(0.d0,0.d0,0.d0)
       goto 666
    endif
       

    jp1 = j + 1
    if (jp1 > n) then
       curl = VECTOR(0.d0,0.d0,0.d0)
       goto 666
    endif
    jm1 = j - 1
    if (jm1 < 1) then
       curl = VECTOR(0.d0,0.d0,0.d0)
       goto 666
    endif

    kp1 = k + 1
    if (kp1 > n) then
       curl = VECTOR(0.d0,0.d0,0.d0)
       goto 666
    endif
    km1 = k - 1
    if (km1 < 1) then
       curl = VECTOR(0.d0,0.d0,0.d0)
       goto 666
    endif


    dFzdy = box(i,jp1,k)%z - box(i,jm1, k)%z
    dFydz = box(i,j,kp1)%y - box(i,j,km1)%y
    dFxdz = box(i,j,kp1)%x - box(i,j,km1)%x
    dFzdx = box(ip1,j,k)%z - box(im1,j,k)%z
    dFydx =  box(ip1,j,k)%y - box(im1,j,k)%y
    dFxdy =  box(i,jp1,k)%x - box(i,jm1,k)%x
!    write(*,*) "box ",box(i,j,k)
!    write(*,*) dfZdy,dFydz,dFxdz,dFzdx, dFydx,dFxdy
    curl%x = dfZdy - dFydz
    curl%y = dFxdz - dFzdx
    curl%z = dFydx - dFxdy
666 continue
  end function curl

  real(double) function div(box,n,i,j,k)
    type(VECTOR) :: box(:,:,:)
    integer :: n, i, j, k, ip1, im1, jp1, jm1, kp1, km1

    ip1 = min(n, i + 1)
    im1 = max(1, i - 1)

    jp1 = min(n, j + 1)
    jm1 = max(1, j - 1)

    kp1 = min(n, k + 1)
    km1 = max(1, k - 1)
    div = 0.
    div = div + box(ip1,j,k)%x - box(im1,j,k)%x
    div = div + box(i,jp1,k)%y - box(i,jm1,k)%y
    div = div + box(i,j,kp1)%z - box(i,j,km1)%z

  end function div
end module turbulence_mod
