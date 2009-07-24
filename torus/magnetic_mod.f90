! Magnetic field stuff added in order to compute zeeman profiles
! TJH in December 2003
! This is for the AMR grid style only

module magnetic_mod

  implicit none

  public

contains

  recursive subroutine addDipoleField(thisOctal, theta, rStar, bStrength, dipoleAxis)
    
  use vector_mod
  use octal_mod, only: OCTAL

! this subroutine adds a dipole field to all the grid points

    real, intent=in :: bStrength, rStar
    type(OCTAL), POINTER :: thisOctal
    type(OCTAL), POINTER :: child
    integer :: i




    if ( thisOctal%nChildren > 0) then
       do i = 1, thisOctal%nChildren
          child => thisOctal%child(i)
          call addDipoleField(child, theta, rStar, bStrength, dipoleAxis)
       end do
    else
       thisCellCentre = subcellCentre(thisOctal, subcell)
       r  = modulus(thisCellCentre)
       rHat = thisCellCentre / r
       theta = acos(rHat .dot. dipoleAxis)
       phi = atan2(thisCellCentre%y, thisCellCentre%x)
       bz = bStrength * (3.*cos(theta)**2-1.) / (r/rStar)**3
       bx = bStrength * (3.*sin(theta)*cos(theta))/(r/rStar)**3
       bField = VECTOR(bx, 0., bz)
       call rotateX(bField, dipoleOffset)
       call arbitraryRotate(bField, phi, dipoleAxis)
       thisOctal%bField(subcell) = bField
    endif
  end subroutine addDipoleField


end module magnetic_mod
