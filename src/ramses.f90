module ramses

  ! To do: Fill with missing data outside ramses grid
  !        Report number of points in/outside ramses grid
  !        Set up HI
  !        Set up corner velocities
  !        Split grid according to original grid resolution
  !        Deallocate arrays once grid has been set up
  !        Use constants mod and double instead of kind=8
  !        Write output correctly with MPI
  !        Smarter cell search
  !        Change module name to ramses_mod
  
  implicit none

  private

  public rd_gas, fillRamses, splitRamses

  real(kind=8), dimension(:),   allocatable :: dx,HI,temp,rho,mg,ratio,nH
  real(kind=8), dimension(:,:), allocatable :: xg,vg
  integer :: nleaf
  
  contains

    subroutine rd_gas(fname)
      use messages_mod
      
      implicit none
      
      real(kind=8), parameter :: Msol = 1.9891d+33
      real(kind=8), parameter :: kpc = 3.086d+21
      real(kind=8),parameter ::mH = 1.6600000d-24
      real(kind=8),parameter ::X = 0.76
      real(kind=8) :: scale_l,scale_d,scale_t,aexp
      character(len=*) :: fname

      call writeInfo("Reading SeleneTORUSinp.dat", TRIVIAL)
      
      open(13, file=fname, status='unknown', form='unformatted')
      rewind 13
      read(13) scale_l,scale_d,scale_t,aexp
      read(13) nleaf
      print *,'nleaf, aexp, scale_l,scale_d,scale_t ',nleaf,aexp,scale_l,scale_d,scale_t
      allocate(xg(nleaf,3),vg(nleaf,3),mg(nleaf),rho(nleaf),dx(nleaf),temp(nleaf),HI(nleaf),ratio(nleaf),nH(nleaf))
      read(13) xg(:,:)
      read(13) vg(:,:)
      read(13) mg(:)
      read(13) rho(:)
      read(13) dx(:)
      read(13) temp(:)
      read(13) HI(:)
      close(13)
  
      print *,'Min/Max values in code units'
      print *,'x ',minval(xg(:,1)),maxval(xg(:,1))
      print *,'y ',minval(xg(:,2)),maxval(xg(:,2))
      print *,'z ',minval(xg(:,3)),maxval(xg(:,3))
      print *,'vx ',minval(vg(:,1)),maxval(vg(:,1))
      print *,'vy ',minval(vg(:,2)),maxval(vg(:,2))
      print *,'vz ',minval(vg(:,3)),maxval(vg(:,3))
      print *,'mg ',minval(mg),maxval(mg)
      print *,'rho ',minval(rho),maxval(rho)
      print *,'dx ',minval(dx),maxval(dx)
      print *,'T ',minval(temp),maxval(temp)
      print *,'HI ',minval(HI),maxval(HI)
      print *
      
      xg=xg*scale_l/kpc             ! convert to kpc
      vg=vg*scale_l/scale_t/1d5     ! convert to km/s
      mg=mg*scale_d*scale_l**3/Msol ! convert to Msol
      nH = rho*X/mH                 ! convert to nH/cm^3
      rho=rho*scale_d               ! convert to g/cc
      dx=dx*scale_l/kpc             ! convert to kpc
      HI = 10**HI                   ! convert to nHI/cm^3
      
      print *,'Min/Max values in useful units'
      print *,'x (kpc)   ', minval(xg(:,1)),maxval(xg(:,1))
      print *,'y (kpc)   ', minval(xg(:,2)),maxval(xg(:,2))
      print *,'z (kpc)   ', minval(xg(:,3)),maxval(xg(:,3))
      print *,'vx (km/s) ',minval(vg(:,1)),maxval(vg(:,1))
      print *,'vy (km/s) ',minval(vg(:,2)),maxval(vg(:,2))
      print *,'vz (km/s) ',minval(vg(:,3)),maxval(vg(:,3))
      print *,'mg (Msol) ',minval(mg),maxval(mg)
      print *,'rho (g/cc)',minval(rho),maxval(rho)
      print *,'dx (kpc)  ',minval(dx),maxval(dx)
      print *,'T         ',minval(temp),maxval(temp)
      print *,'nHI       ',minval(HI),maxval(HI)
      print *,'nH        ',minval(nH),maxval(nH)
      
    end subroutine rd_gas

    subroutine fillRamses(thisOctal, subcell)
      use octal_mod
      implicit none

      type(OCTAL) :: thisOctal
      type(VECTOR) :: position
      integer, intent(in) :: subcell
      integer :: i, iup, idown
      integer, save :: prevIndex=1

      position = subcellCentre(thisOctal, subcell)
      
      ! See if this point is in the same cell as last time
      if (incell(prevIndex, position)) then
         call fillFromCell(prevIndex)
      else
         
      ! Otherwise search for the cell assuming it is close to the previous cell
         iup   = prevIndex+1
         idown = prevIndex-1
         cellLoop: do i=1, nleaf

            if (idown > 0) then
               if (incell(idown, position)) then
                  call fillFromCell(idown)
                  prevIndex=idown
                  exit cellLoop
               end if
            end if

            if (iup < nleaf+1) then
               if (incell(iup, position)) then
                  call fillFromCell(iup)
                  prevIndex=iup
                  exit cellLoop
               end if
            endif

            idown = idown-1
            iup   = iup+1
         
         end do cellLoop

      endif

    contains

      subroutine fillFromCell(i)
        integer, intent(in) :: i
        
        thisOctal%rho(subcell)         = rho(i)
        thisOctal%temperature(subcell) = temp(i)
        thisOctal%velocity(subcell)    = VECTOR(vg(i,1), vg(i,2), vg(i,3))
      end subroutine fillFromCell
      
    end subroutine fillRamses

    logical function splitRamses(size, position)
      use constants_mod
      use vector_mod
      implicit none
      integer :: index
      
      real(double), intent(in) :: size
      type(VECTOR), intent(in) :: position
      real(double) :: size_kpc
      
      size_kpc = size * 1.0e10_db/kpcToCm
      splitRamses=.false. ! In case the point is not within the Ramses grid
      do index=1, nleaf
         if (incell(index, position)) then
            splitRamses=splitThisCell()
            exit
         end if
      end do

    contains

      logical function splitThisCell()
        if ( size_kpc > dx(index) ) then
           splitThisCell=.true.
        else
           splitThisCell=.false.
           end if
      end function splitThisCell
      
    end function splitRamses
    
    logical function inCell(index, position)
      use constants_mod
      use vector_mod
      implicit none
      
      integer, intent(in)      :: index
      type(VECTOR), intent(in) :: position
      real(double) ::  x,  y,  z
      real(double) :: diffx, diffy, diffz
      real (double) :: halfCell

      halfCell = dx(index)/2.0
      
      ! Convert position of cell centre from Torus units to kpc
      x = position%x * 1.0e10_db/kpcToCm
      y = position%y * 1.0e10_db/kpcToCm
      z = position%z * 1.0e10_db/kpcToCm
      
      diffx = xg(index,1) - x
      diffy = xg(index,2) - y
      diffz = xg(index,3) - z

      if ( abs(diffx) <= halfcell .and. abs(diffy) <= halfcell .and. abs(diffz) <= halfcell) then
         incell = .true.
      else
         incell = .false.
      endif
      
    end function inCell
    
  end module ramses

  
