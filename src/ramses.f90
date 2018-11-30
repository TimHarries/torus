module ramses

  ! To do: Get fillRamses to assign correct values
  !        Fill with missing data outside ramses grid
  !        Report number of points in/outside ramses grid
  !        Set up HI
  !        Set up corner velocities
  !        Split grid according to original grid resolution
  !        Deallocate arrays once grid has been set up
  !        Use constants mod and double instead of kind=8
  !        Write output correctly with MPI
  
  implicit none

  private

  public rd_gas, fillRamses

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
      integer, intent(in) :: subcell
      integer :: index

      do index=1, nleaf
         if (incell(index, subcellCentre(thisOctal, subcell))) then 
            thisOctal%rho(subcell)         = rho(index)
            thisOctal%temperature(subcell) = temp(index)
            thisOctal%velocity(subcell)    = VECTOR(vg(index,1), vg(index,2), vg(index,3))
            exit
         end if
      end do
      
    end subroutine fillRamses

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

  
