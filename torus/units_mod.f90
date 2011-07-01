module units_mod

use constants_mod

implicit none

contains

  subroutine convertToTorusUnits(unitString, unitType, torusUnit)

    implicit none
    
    character(len=*) :: unitString
    real(double) :: torusUnit


    select case(unitString)

!Distance units - TORUS uses cm 
    case("cm")
       torusUnit = 1.d0
       
    case("m")
       torusUnit = 1.d2
       
    case("au")
       torusUnit = auToCm
       
    case("pc")
       torusUnit = pcToCm
       
    case("rSol")
       torusUnit = rSol


!Wavelength/dust units - TORUS uses Angstroms (A) and microns (um)
    case("A")
       torusUnit = 1.d0 
       
    case("nm")
       torusUnit = 10.d0
       
    case("um")
       if(unitType == "wavelength") then
          torusUnit = 1.d4
       elseif(unitType == "dust") then
          torusUnit = 1.d0
       end if

    case("mm")
       if(unitType == "wavelength") then
          torusUnit = 1.d7
       elseif(unitType == "dust") then
          torusUnit = 1.d3
       end if


!Angular units - TORUS uses radians
       
    case("rad")
       torusUnit = 1.d0
       
    case("deg")
       torusUnit = degToRad
       
    case("arcmin")
       torusUnit = degToRad/60.d0
       
    case("arcsec")
       torusUnit = degToRad/3600.d0        

       
!Mass units - torus uses g

    case("g")
       torusUnit = 1.d0
       
    case("mSol")
       torusUnit = mSol
       
    case("kg")
       torusUnit = 1.d3

       
!Time units - TORUS uses s
    case("s")
       torusUnit = 1.d0
       
    case("yr")
       torusUnit = 1.d0/secsToYears
       
    case("kyr")
       torusUnit = 1.d0/(secsToYears*1.d3)
       
    case("Myr")
       torusUnit = 1.d0/(secsToYears*1.d6)
       

!Temperature units - TORUS uses K
    case("K")
       torusUnit = 1.d0


!Luminosity units - TORUS uses erg s^-1
    case("ergsec")
       torusUnit = 1.d0

    case("lSol")
       torusUnit = lsol
       
    case(default)
       write(*,*) "Unit not recognized """, unitString """
       write(*,*) "Using TORUS defaults"
       torusUnit = 1.d0

    end select
    
    
  end subroutine convertToTorusUnits

end module
