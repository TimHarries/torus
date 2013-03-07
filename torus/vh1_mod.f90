module vh1_mod

  use kind_mod

  implicit none

! Public components
  public :: read_vh1, assign_from_vh1, get_density_vh1, vh1FileRequired, &
       setGridFromVh1Parameters
  
  private

! Private components
  integer, private :: nx
  integer, private :: ny
  real(db), private, save, allocatable :: rho(:,:)
  real(db), private, save, allocatable :: xaxis(:)
  real(db), private, save, allocatable :: yaxis(:)
  real(db), private, save, allocatable :: vx(:,:), vy(:,:)

  real(db), private, parameter :: xSource = 3.3e7_db
  logical, private, parameter  :: centreOnSource=.true.

  logical, private  :: isRequired=.false.
  character(len=80), private :: infile

contains

! Set module variables from values the parameters file. Called from inputs_mod.
! Could add xSource and centreOnSource 
  subroutine setGridFromVh1Parameters(vh1filename)
    character(len=80), intent(in) :: vh1filename

    infile=vh1filename

    isRequired=.true.

  end subroutine setGridFromVh1Parameters

! Accessor function
  logical function  vh1FileRequired()
     vh1FileRequired=isRequired
  end function vh1FileRequired

  subroutine read_vh1
    
    ! Read data from VH-1 output file written by printold.f
    
    use messages_mod

    implicit none

    character(len=7) :: head1
    character(len=6) :: head2
    character(len=4) :: head3
    integer :: vh1_cycle
    real :: time, dt

    real(db) :: pres
    integer :: i, j
    integer :: nlines

    character(len=200) message

! Work out the size of the grid and allocate storage. 
! nlines includes one header line i.e.
! nx = sqrt(l+1) - 1 
! if l is lines of data only

    nlines = file_line_count(infile)
    nx = int(sqrt(real(nlines))) - 1
    ny = int(sqrt(real(nlines))) - 1

    write(message,*) "Found ", nlines, " lines in ", infile
    call writeInfo(message, FORINFO)
    write(message,*) "Assuming grid is square: nx, ny= ", nx, ny
    call writeInfo(message, FORINFO)

    allocate(rho(nx,ny))
    allocate(vx(nx,ny))
    allocate(vy(nx,ny))
    allocate(xaxis(nx))
    allocate(yaxis(ny))

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
          read(10,*) rho(i,j), pres, vy(i,j), vx(i,j)
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

! Apply offset to put star at the origin
    if (centreOnSource) then 
       xaxis(:) = xaxis(:) - xSource
       write(message,*) "X-axis offset by ", xSource, "to centre grid on source"
       call writeInfo(message, FORINFO)
    end if

    close(10)

    call writeInfo("Finished reading VH-1 data", FORINFO)
    write(message,*) "x-axis (min/max) =", xaxis(1), xaxis(nx)
    call writeInfo(message, FORINFO)
    write(message,*) "y-axis (min/max) =", yaxis(1), yaxis(nx)
    call writeInfo(message, FORINFO)

  contains

! Count the number of lines in a file. 
! Should return the same answer as wc -l i.e. includes blank lines in the total
! D. Acreman, June 2010
    integer function file_line_count(filename)

      character(len=*) :: filename
      character(len=1) :: dummy
      integer :: status 

      file_line_count = 0
      open(unit=30, status="old", file=filename)
      do
         read(30,'(a1)',iostat=status) dummy
         if ( status /= 0 ) exit
         file_line_count = file_line_count + 1 
      end do
      close(30)

    end function file_line_count

  end subroutine read_vh1

  subroutine assign_from_vh1(thisOctal, subcell)

    use octal_mod
    use vector_mod
    implicit none

    real(db), parameter :: rho_bg=1.0e-23_db

    TYPE(OCTAL), intent(inout) :: thisOctal
    integer, intent(in) :: subcell

    TYPE(vector) :: thisCentre, thisVel, sourceToCell
    real(db) :: vOutflow
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

! Extract cell velocity from the VH-1 grid
       thisVel=VECTOR(vy(this_j, this_i), 0.0, vx(this_j, this_i))
! Store velocity for writing into VTK file. Needs to be as fraction of c. 
       thisOctal%velocity(subcell) = thisVel/cspeed
! Vector from source position to this cell. Remember that VH-1 x-axis is Torus z-axis.
       if (centreOnSource) then
          sourceToCell = thisCentre - VECTOR(0.0, 0.0,0.0)
       else
          sourceToCell = thisCentre - VECTOR(0.0, 0.0, xSource)
       endif
! Project velocity onto this vector to get the outflow velocity
       vOutflow = (thisVel.dot.sourceToCell)/modulus(sourceToCell)

! Set up initial dust distribution
       if (vOutflow > 2000.0*1e5) then 
! No dust in outflow
          thisOctal%dustTypeFraction(subcell,:) = 0.0
       else
          thisOctal%dustTypeFraction(subcell,:) = 0.01
       endif

    end if

    thisOctal%temperature(subcell) = 10000.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-8
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%inFlow(subcell) = .true.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

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
