module ramses_mod
  use kind_mod
  
  ! Module for working with Gareth Few's Ramses data
  ! The rd_gas subroutine is based on code contributed by Gareth
  ! D. Acreman, December 2018
  
  implicit none

  private

  public setGridFromRamsesParameters, rd_gas, fillRamses, splitRamses, finishRamses

  character(len=80) :: fname
  real(double), dimension(:),   allocatable :: dx,HI,temp,rho,mg,ratio,nH
  real(double), dimension(:,:), allocatable :: xg,vg
  integer :: nleaf
  integer, save :: num_outside=0
  
  contains

    subroutine setGridFromRamsesParameters(ramsesfilename)
      use messages_mod

      character(len=*), intent(in) :: ramsesfilename
      logical :: foundFile
      
      fname=ramsesfilename

      inquire(file=fname, exist=foundFile)
      if (.not. foundFile) then
         call writeFatal("File "//trim(fname)//" not found")
      end if
      
    end subroutine setGridFromRamsesParameters
      
    subroutine rd_gas
      use messages_mod
      use constants_mod
      implicit none
      
      real(double), parameter :: kpc = kpcToCm
      real(double), parameter :: mH = mHydrogen
      real(double), parameter :: X = 0.76
      real(double) :: scale_l,scale_d,scale_t,aexp
      character(len=120) :: message

      call writeInfo("Reading "//fname, TRIVIAL)
      
      open(13, file=fname, status='old', form='unformatted')
      rewind 13
      read(13) scale_l,scale_d,scale_t,aexp
      read(13) nleaf
      write(message,*) 'nleaf, aexp: ', nleaf, aexp
      call writeInfo(message, TRIVIAL)
      write(message,*) 'scale_l, scale_d, scale_t: ', scale_l, scale_d, scale_t
      call writeInfo(message, TRIVIAL)
      allocate(xg(nleaf,3),vg(nleaf,3),mg(nleaf),rho(nleaf),dx(nleaf),temp(nleaf),HI(nleaf),ratio(nleaf),nH(nleaf))
      read(13) xg(:,:)
      read(13) vg(:,:)
      read(13) mg(:)
      read(13) rho(:)
      read(13) dx(:)
      read(13) temp(:)
      read(13) HI(:)
      close(13)
  
      call writeInfo('Min/Max values in code units:', TRIVIAL)
      write(message,*) 'x ',minval(xg(:,1)),maxval(xg(:,1))
      call writeInfo(message, TRIVIAL)
      write(message,*) 'y ',minval(xg(:,2)),maxval(xg(:,2))
      call writeInfo(message, TRIVIAL)
      write(message,*) 'z ',minval(xg(:,3)),maxval(xg(:,3))
      call writeInfo(message, TRIVIAL)
      write(message,*) 'vx ',minval(vg(:,1)),maxval(vg(:,1))
      call writeInfo(message, TRIVIAL)
      write(message,*) 'vy ',minval(vg(:,2)),maxval(vg(:,2))
      call writeInfo(message, TRIVIAL)
      write(message,*) 'vz ',minval(vg(:,3)),maxval(vg(:,3))
      call writeInfo(message, TRIVIAL)
      write(message,*) 'mg ',minval(mg),maxval(mg)
      call writeInfo(message, TRIVIAL)
      write(message,*) 'rho ',minval(rho),maxval(rho)
      call writeInfo(message, TRIVIAL)
      write(message,*) 'dx ',minval(dx),maxval(dx)
      call writeInfo(message, TRIVIAL)
      write(message,*) 'T ',minval(temp),maxval(temp)
      call writeInfo(message, TRIVIAL)
      write(message,*) 'HI ',minval(HI),maxval(HI)
      call writeInfo(message, TRIVIAL)
      call writeInfo("", TRIVIAL)
      
      xg=xg*scale_l/kpc             ! convert to kpc
      vg=vg*scale_l/scale_t/1d5     ! convert to km/s
      mg=mg*scale_d*scale_l**3/Msol ! convert to Msol
      nH = rho*X/mH                 ! convert to nH/cm^3
      rho=rho*scale_d               ! convert to g/cc
      dx=dx*scale_l/kpc             ! convert to kpc
      HI = 10**HI                   ! convert to nHI/cm^3
      
      call writeInfo('Min/Max values in useful units:', TRIVIAL)
      write(message,*) 'x (kpc)   ', minval(xg(:,1)),maxval(xg(:,1))
      call writeInfo(message, TRIVIAL)
      write(message,*) 'y (kpc)   ', minval(xg(:,2)),maxval(xg(:,2))
      call writeInfo(message, TRIVIAL)
      write(message,*) 'z (kpc)   ', minval(xg(:,3)),maxval(xg(:,3))
      call writeInfo(message, TRIVIAL)
      write(message,*) 'vx (km/s) ',minval(vg(:,1)),maxval(vg(:,1))
      call writeInfo(message, TRIVIAL)
      write(message,*) 'vy (km/s) ',minval(vg(:,2)),maxval(vg(:,2))
      call writeInfo(message, TRIVIAL)
      write(message,*) 'vz (km/s) ',minval(vg(:,3)),maxval(vg(:,3))
      call writeInfo(message, TRIVIAL)
      write(message,*) 'mg (Msol) ',minval(mg),maxval(mg)
      call writeInfo(message, TRIVIAL)
      write(message,*) 'rho (g/cc)',minval(rho),maxval(rho)
      call writeInfo(message, TRIVIAL)
      write(message,*) 'dx (kpc)  ',minval(dx),maxval(dx)
      call writeInfo(message, TRIVIAL)
      write(message,*) 'T         ',minval(temp),maxval(temp)
      call writeInfo(message, TRIVIAL)
      write(message,*) 'nHI       ',minval(HI),maxval(HI)
      call writeInfo(message, TRIVIAL)
      write(message,*) 'nH        ',minval(nH),maxval(nH)
      call writeInfo(message, TRIVIAL)
      call writeInfo("Finished reading "//fname, TRIVIAL)
      
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

            ! If we're off the array both ends then this point is not in the grid
            ! so set missing data values
            if (idown<1 .and. iup>nleaf) then
               thisOctal%rho(subcell)         = 1.0e-33_db
               thisOctal%temperature(subcell) = 10.0_db
               thisOctal%velocity(subcell)    = VECTOR(2.0_db,2.0_db,2.0_db)
               num_outside=num_outside+1
               exit cellLoop
            end if
            
            idown = idown-1
            iup   = iup+1
         
         end do cellLoop

      endif

    contains

      subroutine fillFromCell(i)
        use constants_mod
        integer, intent(in) :: i

        ! Total gas density
        !thisOctal%rho(subcell)         = rho(i)
        ! HI density
        thisOctal%rho(subcell)         = HI(i) * mHydrogen
        thisOctal%temperature(subcell) = temp(i)
        thisOctal%velocity(subcell)    = VECTOR(vg(i,1)* 1.0e5/cSpeed, vg(i,2)* 1.0e5/cSpeed, vg(i,3)* 1.0e5/cSpeed) 
      end subroutine fillFromCell
      
    end subroutine fillRamses

    logical function splitRamses(size, position)
      use constants_mod
      use vector_mod
      implicit none
      
      real(double), intent(in) :: size
      type(VECTOR), intent(in) :: position
      real(double) :: size_kpc
      integer :: i, iup, idown
      integer, save :: prevIndex=1

      ! In case the point is not within the Ramses grid
      splitRamses=.false.
      size_kpc = size * 1.0e10_db/kpcToCm

      ! Is this point the same cell as last time?
      if (incell(prevIndex, position)) then
         splitRamses=splitThisCell(prevIndex)
      else

      ! Otherwise search for the cell assuming it is close to the previous cell
         iup   = prevIndex+1
         idown = prevIndex-1
         cellLoop: do i=1, nleaf

            if (idown > 0) then
               if (incell(idown, position)) then
                  splitRamses=splitThisCell(idown)
                  prevIndex=idown
                  exit cellLoop
               end if
            end if

            if (iup < nleaf+1) then
               if (incell(iup, position)) then
                  splitRamses=splitThisCell(iup)
                  prevIndex=iup
                  exit cellLoop
               end if
            endif

            idown = idown-1
            iup   = iup+1
         
         end do cellLoop
      end if
      
    contains

      logical function splitThisCell(i)
        integer, intent(in) :: i
        if ( size_kpc > dx(i) ) then
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

    ! Deallocate arrays and report information
    subroutine finishRamses
      use messages_mod
      implicit none
      character(len=80) :: message

      write(message,*) "Number of points allocated missing data values: ", num_outside
      call writeInfo(trim(message), TRIVIAL)
      deallocate(xg, vg, mg, rho, dx, temp, HI, ratio, nH)
      
    end subroutine finishRamses
    
  end module ramses_mod
