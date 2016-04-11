module gridanalysis_mod
  use kind_mod
  use gridtype_mod
  use amr_utils_mod
  use random_mod
  use vector_mod

contains

  subroutine analysis(grid)
    type(GRIDTYPE) :: grid

    flux = fluxThroughSphericalSurface(grid, 1000.d0*autocm/1.d10)
    flux = flux / msol * 365.25d0*24.d0*3600.d0
    write(*,*) "flux through 1000 AU sphere is ",flux, " solar masses/year"
  end subroutine analysis



  function fluxThroughSphericalSurface(grid, radius) result(totflux)
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: radius
    integer :: i, n
    real(double) :: flux, totflux
    type(VECTOR) :: rVec, inVec
    n = 10000
    da = fourPi*(radius*1.d10)**2/dble(n)
    totflux = 0.d0
    do i = 1, n

       rVec = randomUnitVector()

       inVec = (-1.d0)*rVec
       call normalize(inVec)

       call findSubcellLocal(rVec, thisOctal, subcell)
       flux = invec%x * thisOctal%rhou(subcell) + &
            invec%y * thisOctal%rhov(subcell) + &
            invec%z * thisOctal%rhow(subcell)
       totFlux = totFlux + flux * da
    enddo
  end function fluxThroughSphericalSurface

end module gridanalysis_mod
