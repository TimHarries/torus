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
      use amr_mod,      only: amrGridValues
      use constants_mod, only: auToCm

! Arguments
      type(GRIDTYPE), intent(in)   :: grid
      character(len=*), intent(in) :: filename

! Local variables
      character(len=100) :: message
      integer            :: num_lines, iline
      integer, parameter :: lun_in   = 60
      integer, parameter :: lun_out  = 61
      integer, parameter :: lun_out2 = 62
      integer, parameter :: lun_out3 = 63
      integer, parameter :: lun_out4 = 64
      real :: temperature, temperature_grid, temperature_diff
      real(double) :: density, density_grid, density_diff
      real(double) :: part_mass, h, x_au, y_au, z_au
      TYPE(vector) :: point
      logical :: do_check

      inquire(file=filename, exist=do_check)
      if ( do_check ) then 

         open(unit=lun_in,   status='old',     file=filename                )
         open(unit=lun_out,  status='replace', file='temperature_stats.dat' )
         open(unit=lun_out2, status='replace', file='density_stats.dat'     )
         open(unit=lun_out3, status='replace', file='part_diffs.dat'        )
         open(unit=lun_out4, status='replace', file='part_frac_diffs.dat'   )

         read(lun_in,*) num_lines      
         write(message,*) "Reading ", num_lines, "lines" 
         call writeInfo(trim(message), TRIVIAL)

         do iline = 1, num_lines

            read(60,*) point%x, point%y, point%z, part_mass, h, density, temperature

            ! Distance in AU
            x_au = point%x / auToCm
            y_au = point%y / auToCm
            z_au = point%z / auToCm

            ! Convert from cm to torus units
            point%x = point%x * 1.0e-10
            point%y = point%y * 1.0e-10
            point%z = point%z * 1.0e-10

            call amrGridValues(grid%octreeRoot, point, rho=density_grid, temperature=temperature_grid)

            temperature_diff = temperature - temperature_grid
            density_diff     = density     - density_grid

            write(lun_out,*)  temperature, temperature_grid, temperature_diff
            write(lun_out2,*) density,     density_grid,     density_diff  
            write(lun_out3,'(7(e15.8,2x))') x_au, y_au, z_au, part_mass, h, density_diff, temperature_diff
            write(lun_out4,'(7(e15.8,2x))') x_au, y_au, z_au, part_mass, h, &
                 (density_diff/density_grid), (temperature_diff/temperature_grid)

         end do

         close( lun_in   )
         close( lun_out  )
         close( lun_out2 )
         close( lun_out3 )
         close( lun_out4 )

      end if

    end subroutine check_benchmark_values

end module benchmark_mod

