module vh1_mod

  use kind_mod

  implicit none

  public :: read_vh1, assign_from_vh1
  
  integer, parameter, private :: nx=800
  integer, parameter, private :: ny=800
  real(db), private, save :: rho(nx,ny)
  real(db), private, save :: xaxis(nx)
  real(db), private, save :: yaxis(ny)

contains


  subroutine read_vh1
    
    ! Read data from VH-1 output file written by printold.f
    
    use messages_mod

    implicit none

    character(len=*), parameter :: infile="bowsc.1008"
    character :: head1*7, head2*6, head3*4
    integer :: vh1_cycle
    real :: time, dt

    real(db) :: pres, vx, vy
    integer :: i, j

    character(len=200) message

    open(unit=10, status="old", file=infile)
    
! Open file
    write(message,*) "Opening file "//infile
    call writeInfo(message, FORINFO)

! Read header
    read(10,*) head1, vh1_cycle, head2, time, head3, dt
    write(message,*) "Header: ", head1, vh1_cycle, head2, time, head3, dt
    call writeInfo(message, TRIVIAL)
    
! Read vh1 data values
    do j=1, ny
       do i=1, nx
          read(10,*) rho(i,j), pres, vy, vx
       end do
    end do

! Read x-axis values
    do i=1, nx
       read(10,*) xaxis(i)
    end do

! Read y-axis values
    do j=1, ny
       read(10,*) yaxis(j)
    end do

! Convert axes to Torus units
    xaxis(:) = xaxis(:) / 1.0e10_db
    yaxis(:) = yaxis(:) / 1.0e10_db

    close(10)

    call writeInfo("Finished reading VH-1 data", FORINFO)

  end subroutine read_vh1

  subroutine assign_from_vh1(thisOctal, subcell)

    use octal_mod
    use vector_mod
    implicit none

    TYPE(OCTAL), intent(inout) :: thisOctal
    integer, intent(in) :: subcell

    TYPE(vector) :: thisCentre
    integer :: this_i, this_j, i, j

    thisCentre = subcellcentre(thisOctal, subcell)

    if ( thisCentre%x < xaxis(1)  .or. &
         thisCentre%x > xaxis(nx) .or. &
         thisCentre%z < yaxis(1)  .or. &
         thisCentre%z > yaxis(ny) ) then 
       
       thisOctal%rho(subcell) = 1.0e-33_db

    else

       do i=2, nx
          if ( thisCentre%x < xaxis(i) ) then 
             this_i = i
             exit
          end if
       end do

! In a 2D grid the octal y value doesn't change
       do j=2, ny
          if ( thisCentre%z < yaxis(j) ) then 
             this_j = j
             exit
          end if
       end do

       thisOctal%rho(subcell) = rho(this_i, this_j)
       thisOctal%temperature(subcell) = 10.0

    end if

  end subroutine assign_from_vh1

end module vh1_mod
