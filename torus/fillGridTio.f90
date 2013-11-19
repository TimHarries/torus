module fillGridTio_mod

  public  :: fillGridTio
  private :: readTioCrossSection 

contains

  subroutine fillGridTio(grid,  scale)

    use gridtype_mod, only: gridtype
    use constants_mod

    implicit none
    real :: scale
    type(GRIDTYPE) :: grid
    real :: abundance
    real :: massTio
    real, allocatable :: xArray(:), sigmaArray(:)
    integer :: i, j, k, m

    massTio = real(15.9994*mHydrogen + 47.867*mHydrogen)
    abundance = 1.e-7


    allocate(xArray(1:grid%nlambda))
    allocate(sigmaArray(1:grid%nlambda))

    xArray = grid%lamarray


    call readTioCrossSection(xArray, grid%nLambda, sigmaArray)

    write(*,*) "Filling grid with TiO opacity"

    if (.not.grid%cartesian) then
       do i = 1, grid%nr
          do j = 1, grid%nmu
             do k = 1, grid%nPhi
                do m = 1, grid%nLambda
                   grid%kappaAbs(i,j,k,m) = grid%kappaAbs(i,j,k,m) + &
                        sigmaArray(m)*abundance*grid%rho(i,j,k)/massTio * scale
                enddo
             enddo
          enddo
       enddo
    else
       do i = 1, grid%nx
          do j = 1, grid%ny
             do k = 1, grid%nz
                do m = 1, grid%nLambda
                   grid%kappaAbs(i,j,k,m) = grid%kappaAbs(i,j,k,m) + &
                        sigmaArray(m)*abundance*grid%rho(i,j,k)/massTio * scale
                enddo
             enddo
          enddo
       enddo
    endif
  end subroutine fillGridTio


  subroutine readTioCrossSection(lamArray, nLambda, sigmaArray)
    use utils_mod
    use unix_mod, only: unixGetenv
    implicit none
    integer :: nLambda
    real :: lamArray(nLambda), sigmaArray(nLambda)
    integer ::  i, j
    real :: xArray(20000)
    real :: sArray(20000)
    real, allocatable :: nArray(:)
    integer :: npts
    character(len=80) :: dataDirectory, filename
    
    call unixGetenv("TORUS_DATA",dataDirectory)
    dataDirectory = trim(dataDirectory)//"/"
    filename = trim(dataDirectory)//"tio.xsec"

    open(20, file=filename, status="old", form="formatted")

    npts = 1
    do
       read(20,*,end=10) xArray(npts), sArray(npts)
       npts=npts+1
    end do
10  continue
    close(20)
    npts=npts - 1

    allocate(nArray(1:nLambda))
    narray = 0.

    sigmaArray = 0.
   

    do i = 1, npts
       if ((xArray(i) >= lamArray(1)) .and. &
            (xArray(i)<=lamArray(nLambda))) then
          call locate(lamArray, nLambda, xArray(i), j)
          sigmaArray(j) = sigmaArray(j) + sArray(i)
          nArray(j) = nArray(j) + 1.
       endif
    enddo

    where(nArray > 0) 
       sigmaArray = sigmaArray/nArray
    endwhere

  end subroutine readTioCrossSection

end module fillGridTio_mod
