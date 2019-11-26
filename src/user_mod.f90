module user_mod

  use vector_mod
  use constants_mod
  use octal_mod
  implicit none

  private
  public userSplit, userSubroutine, userDensity

  real(double) :: rhoFloor = 1.d-24 ! lowest density in the grid in g per cm^3
  real(double) :: rhoMax = 1.d-17 ! maximum density of 
  real(double), parameter :: rSphere = 0.1d0 * pcTocm   ! constant is from constants mod
  real(double), parameter :: rHole = 200.d0 * auTocm     ! constant is from constants mod

contains

  subroutine userSubroutine(thisOctal, subcell)
    type(OCTAL) :: thisOctal
    integer :: subcell
    type(VECTOR) :: centre

    centre = subcellCentre(thisOctal, subcell)
    thisOctal%rho(subcell) = userDensity(centre*1.d10)
    thisOctal%temperature(subcell) = 10.
  end subroutine userSubroutine


! function to return a density at a particular position
  
  real(double) function userDensity(position) result(density) ! returns density in g/cm^3
    type(VECTOR) :: position ! position vector in units of cm

! you can change anything below this to produce the density structure that you like


    real(double) :: r ! radius
    real(double) :: rCylindrical, rCav
    

    r = modulus(position)
    rCylindrical  = sqrt(position%x**2 + position%y**2)/rSphere
    rCav = 0.3d0*sqrt(abs(position%z/rSphere))

    density = rhoFloor

    if ( (r < rSphere).and.(r > rHole) ) then
       density = rhoMax 
    endif

    if (rCylindrical < rCav) then
       density = rhoFloor
    endif
  end function userDensity

  function userSplit(position, size) result(split)
    logical :: split
    type(VECTOR) :: position
    real(double) :: size
    real(double) :: r
    real(double) :: rCav, rCylindrical

    split = .false.

    r = modulus(position)

    rCav = 0.3d0*sqrt(abs(position%z/rSphere))
    rCylindrical = sqrt(position%x**2 + position%y**2)/rSphere

    if ((r < 1.1d0*rSphere).and.(size/rSphere > 0.05)) split = .true.

    if (r > rHole) then
       if (abs(rCav - rCylindrical) < size/rSphere) then
          if (size/rSphere > 0.01d0*rCav) then
             split = .true.
          endif
       endif
    endif
    


    if ( (abs(r-rHole) < size).and.(size > 0.01d0*rHole) ) then
       split = .true.
    endif



  

  end function userSplit

end module user_mod
