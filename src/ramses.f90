module ramses

  implicit none

  private

  public rd_gas, fillRamses
  
  contains

    subroutine rd_gas(fname)
      use messages_mod
      
      implicit none
      
      real(kind=8), parameter :: Msol = 1.9891d+33
      real(kind=8), parameter :: kpc = 3.086d+21
      real(kind=8),parameter ::mH = 1.6600000d-24
      real(kind=8),parameter ::X = 0.76
      integer :: nleaf
      real(kind=8) :: scale_l,scale_d,scale_t,aexp
      real(kind=8),dimension(:),allocatable::dx,HI,temp,rho,mg,ratio,nH
      real(kind=8),dimension(:,:),allocatable::xg,vg
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
  
      print *,'code units'
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
      
      print *,'useful units'
      print *,'x ',minval(xg(:,1)),maxval(xg(:,1))
      print *,'y ',minval(xg(:,2)),maxval(xg(:,2))
      print *,'z ',minval(xg(:,3)),maxval(xg(:,3))
      print *,'vx ',minval(vg(:,1)),maxval(vg(:,1))
      print *,'vy ',minval(vg(:,2)),maxval(vg(:,2))
      print *,'vz ',minval(vg(:,3)),maxval(vg(:,3))
      print *,'mg (Msol) ',minval(mg),maxval(mg)
      print *,'rho ',minval(rho),maxval(rho)
      print *,'dx ',minval(dx),maxval(dx)
      print *,'T ',minval(temp),maxval(temp)
      print *,'nHI ',minval(HI),maxval(HI)
      print *, 'nH',minval(nH),maxval(nH)
      
    end subroutine rd_gas

    subroutine fillRamses(thisOctal, subcell)
      use octal_mod
      implicit none

      type(OCTAL) :: thisOctal
      integer, intent(in) :: subcell

      ! WRITE ME !!!
      ! Fill with dummy values for now
      thisOctal%rho(subcell)         = 1.0e-20
      thisOctal%temperature(subcell) = 100.0
      thisOctal%velocity             = VECTOR(0.,0.,0.)
      
    end subroutine fillRamses
      
  end module ramses

  
