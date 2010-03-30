module vh1_mod

  use kind_mod

  implicit none

  public :: read_vh1, assign_from_vh1, get_density_vh1
  
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

    real(db), parameter :: rho_bg=1.0e-23_db

    TYPE(OCTAL), intent(inout) :: thisOctal
    integer, intent(in) :: subcell

    TYPE(vector) :: thisCentre
    integer :: this_i, this_j, i, j

    thisCentre = subcellcentre(thisOctal, subcell)

! Remeber that the axes are flipped. 
    if ( thisCentre%x < yaxis(1)  .or. &
         thisCentre%x > yaxis(nx) .or. &
         thisCentre%z < xaxis(1)  .or. &
         thisCentre%z > xaxis(ny) ) then 
       
       thisOctal%rho(subcell) = rho_bg

    else

       do i=2, nx
          if ( thisCentre%x < yaxis(i) ) then 
             this_i = i
             exit
          end if
       end do

! In a 2D grid the octal y value doesn't change
       do j=2, ny
          if ( thisCentre%z < xaxis(j) ) then 
             this_j = j
             exit
          end if
       end do

! this_j is for VH-1 x-axis and this_i is for VH-1 y-axis
       thisOctal%rho(subcell) = rho(this_j, this_i)

    end if

  end subroutine assign_from_vh1

  subroutine get_density_vh1(thisOctal, subcell, mean_rho, min_rho, max_rho, ncells)

    use octal_mod

    implicit none 

    TYPE(OCTAL), intent(in) :: thisOctal
    integer, intent(in)     :: subcell
    real(db), intent(out)   :: mean_rho
    real(db), intent(out)   :: min_rho
    real(db), intent(out)   :: max_rho
    integer, intent(out)    :: ncells

    integer :: imin, imax, jmin, jmax
    integer :: i, j
    TYPE(VECTOR) :: cellcentre
    real(db) :: this_loc, sum_rho

! Find cells which are in the octal
! Not efficient for a uniform grid but will work for the non-uniform case
! Remember that: Torus x-axis -> VH-1 y-axis
!                Torus z-axis -> VH-1 x-axis

    cellCentre = subcellCentre(thisOctal,subCell)

    this_loc = cellCentre%z - (0.5 * thisOctal%subcellsize)
    imin = -99
    do i=1, nx
       ! find the first vh1 point within this subcell
       if ( xaxis(i) > this_loc ) then
          imin = i 
          exit
       end if
    end do
    if ( imin == -99 ) then 
       mean_rho = 0.0
       min_rho  = 0.0
       max_rho  = 0.0
       ncells   = -1
       return
    end if

    this_loc = cellCentre%z + (0.5 * thisOctal%subcellsize)
    imax = -99
    do i=nx, 1, -1
       ! find the first vh1 point within this subcell
       if ( xaxis(i) < this_loc ) then
          imax = i 
          exit
       end if
    end do
    if ( imax == -99 ) then 
       mean_rho = 0.0
       min_rho  = 0.0
       max_rho  = 0.0
       ncells   = -2
       return
    end if

    this_loc = cellCentre%x - (0.5 * thisOctal%subcellsize)
    jmin = -99
    do j=1, ny
       ! find the first vh1 point within this subcell
       if ( yaxis(j) > this_loc ) then
          jmin = j 
          exit
       end if
    end do
    if ( jmin == -99 ) then 
       mean_rho = 0.0
       min_rho  = 0.0
       max_rho  = 0.0
       ncells   = -3
       return
    end if

    this_loc = cellCentre%x + (0.5 * thisOctal%subcellsize)
    jmax = -99
    do j=ny, 1, -1
       ! find the first vh1 point within this subcell
       if ( yaxis(j) < this_loc ) then
          jmax = j 
          exit
       end if
    end do
    if ( jmax == -99 ) then 
       mean_rho = 0.0
       min_rho  = 0.0
       max_rho  = 0.0
       ncells   = -4
       return
    end if

    max_rho = -1.0e30_db
    min_rho =  1.0e30_db
    sum_rho =  0.0 
    ncells  =  0 
    do j=jmin, jmax
       do i=imin, imax
          max_rho = max(max_rho, rho(i,j) )
          min_rho = min(min_rho, rho(i,j) )
          sum_rho = sum_rho + rho(i,j)
          ncells  = ncells + 1 
       end do
    end do
    mean_rho = sum_rho / real(ncells,db)

  end subroutine get_density_vh1

end module vh1_mod
