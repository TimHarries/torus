! Module relating to H 21cm images
! 
! Dave Acreman, January 2009


module h21cm_mod

  use kind_mod

  implicit none

  real(double), public, parameter :: h21cm_lambda = 21.12

  public :: hi_emop

  private

contains

!------------------------------------------------------------------------------

  subroutine  hi_emop(dens,temp,emis,opac)

!**************************************************
!
!! HI_EMOP Returns the HI emissivity and opacity
!!  when given the density and temperature of a cell.
!
! Version 1 of subroutine Nov 7 2008
! Kevin Douglas
!
    implicit none

    real, intent(in)         :: temp
    real(double), intent(in) :: dens
    real(double), intent(out):: emis, opac

! minimum density threshold for non-zero emissivity and opacity 
    real(double), parameter :: min_dens = 1.0e-99_db 

    if ( dens < min_dens ) then 
       opac = 0.0
       emis = 0.0
    else
       opac = 1.5472d+09 * dens / temp
       emis = 9.586d-10 * dens
    end if

  end subroutine hi_emop

!------------------------------------------------------------------------------

end module h21cm_mod

