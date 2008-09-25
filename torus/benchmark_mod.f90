module benchmark_mod

  implicit none

  public :: check_benchmark_values

  private

  contains

! Check a list of temperature and density values read from a file 
! against the values on the torus grid. The input file format is 
! 
! x (cm), y(cm), z(cm), density(g/cm^3), temperature(K)
!
! D. Acreman, September 2008

    subroutine check_benchmark_values(grid, filename)

      use kind_mod
      use messages_mod
      use gridtype_mod, only: gridtype
      use vector_mod,   only: vector
      use amr_mod,      only: amrGridDensity, amrGridTemperature

! Arguments
      type(GRIDTYPE), intent(in)   :: grid
      character(len=*), intent(in) :: filename

! Local variables
      character(len=100) :: message
      integer            :: num_lines, iline
      integer, parameter :: lun_in   = 60
      integer, parameter :: lun_out  = 61
      integer, parameter :: lun_out2 = 62
      real :: temperature, temperature_grid, temperature_diff
      real(double) :: density, density_grid, density_diff
      TYPE(vector)  :: point
      logical :: do_check

      inquire(file=filename, exist=do_check)
      if ( do_check ) then 

         open(unit=lun_in,   status='old',     file=filename                )
         open(unit=lun_out,  status='replace', file='temperature_stats.dat' )
         open(unit=lun_out2, status='replace', file='density_stats.dat'     )

         read(lun_in,*) num_lines      
         write(message,*) "Reading ", num_lines, "lines" 
         call writeInfo(trim(message), TRIVIAL)

         do iline = 1, num_lines

            read(60,*) point%x, point%y, point%z, density, temperature

            ! Convert from cm to torus units
            point%x = point%x * 1.0e-10
            point%y = point%y * 1.0e-10
            point%z = point%z * 1.0e-10

            density_grid     = amrGridDensity(grid%octreeRoot, point)
            temperature_grid = amrGridTemperature(grid%octreeRoot, point)

            temperature_diff = temperature - temperature_grid
            density_diff     = density     - density_grid

            write(lun_out,*)  temperature, temperature_grid, temperature_diff
            write(lun_out2,*) density,     density_grid,     density_diff  

         end do

         close( lun_in   )
         close( lun_out  )
         close( lun_out2 )

      end if

    end subroutine check_benchmark_values

end module benchmark_mod

