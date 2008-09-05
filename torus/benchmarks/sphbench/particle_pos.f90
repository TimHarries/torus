module particle_pos_mod

  implicit none

  public :: particle_pos

  contains

    subroutine particle_pos(npart)

      implicit none

      integer, parameter :: db = selected_real_kind(15,307)

      real, parameter :: r_out = 1000.0      ! AU
      real, parameter :: r_d   = r_out / 2.0 ! AU
      real, parameter :: z_d   = r_out / 8.0 ! AU 
      real(db), parameter :: pi = 3.1415926535897932_db
      real(db), parameter :: rho_zero = 0.81614E-17_db
      real(db), parameter :: auToCm = 1.495979d13
      real(db), parameter :: mSol = 1.9891e33_db

      real(db) :: this_mass
      integer  :: i, part 
      integer, parameter  :: imax=10000 ! Radial sampling
      integer, intent(in) :: npart
      real :: r, pdf(imax)
      real :: ran_num, part_r, part_z

      integer, parameter :: npts=10000 ! ! z sampling
      real :: gaus_pdf(0:npts)
      real :: z_sig(0:npts)

      ! Calculate look up table of mass out to given radius
      do i=1, imax
         r = ( real(i) / real(imax) ) * r_out
         call mass_vs_r(1.0, r, this_mass)
         pdf(i) = this_mass
      end do

      ! Normalise by total mass to get probability
      pdf(:) = pdf(:) / pdf(imax)

      call setup_gaussian(z_sig, gaus_pdf, npts)

      open (unit=61, file="part.dat")
      do part=1, npart

         call random_number(ran_num)

         do i=1, imax
            if ( pdf(i) > ran_num ) then 
               part_r = ( real(i) / real(imax) ) * r_out
               exit
            end if
         end do
     
         call get_z(part_r, part_z)

         write(61,*) part_r, part_z

      end do
      close(61)

    contains

      subroutine mass_vs_r(r_in, this_r, mass)
        
        implicit none
    
        real, intent(in)      :: r_in, this_r ! AU
        real(db), intent(out) :: mass

        real :: u_in, this_u

        this_u = this_r / r_d 
        u_in   = r_in   / r_d 

        mass = ( (this_u ** 2.125 ) / 2.125 ) - &
             ( (u_in   ** 2.125 ) / 2.125 )

        mass = 4.0 * pi * rho_zero * r_d**2 * z_d * mass * auToCm**3 

        mass = mass / msol

      end subroutine mass_vs_r

! Sample z distribution for a given r value. 
      subroutine get_z(part_r, part_z)

        implicit none

        real, intent(in)  :: part_r
        real, intent(out) :: part_z

        real :: h, sigma
        real :: z_ran_num

        integer :: k 

        h     = z_d * ( (part_r/r_d)**1.125 )
        sigma = h * sqrt(2.0/pi)

        call random_number(z_ran_num)

        do k=1, npts
           if ( gaus_pdf(k) > z_ran_num ) then 
              part_z = z_sig(k) * sigma
              exit
           end if
        end do

      end subroutine get_z

      subroutine setup_gaussian(x, gaus_pdf, npts)

        implicit none

! Arguments
        integer, intent(in) :: npts
        real, intent(out)   :: gaus_pdf(0:npts)
        real, intent(out)   :: x(0:npts)

! Local
        real, parameter :: xmin = -10.0
        real, parameter :: xmax =  10.0
        real    :: dx
        real    :: gaus
        integer :: ng

        dx = (xmax - xmin) / real(npts)
        gaus_pdf(0) = 0.0
        x(0)        = xmin

        do ng = 1, npts    
           x(ng) = x(ng-1) + dx
           gaus =  exp(-1.0 * x(ng)**2 / 2.0 ) / ( sqrt(2.0 * pi) )
           gaus_pdf(ng) = gaus_pdf(ng-1) + gaus * dx
        end do
        
        if ( abs(gaus_pdf(npts) - 1.0_db)  > 1e-5_db ) then
           print *, "Error: gaussian integral incorrect"
           print *, gaus_pdf(npts)
           stop
        end if

      end subroutine setup_gaussian

    end subroutine particle_pos

  end module particle_pos_mod


