module units_mod

use constants_mod

implicit none

contains

  subroutine convertToTorusUnits(unitString, unitType, torusUnit, inputValue)

    implicit none
    
    character(len=*) :: unitString, unitType
    real(double) :: torusUnit, inputValue
    
    
    if(unitType == "distance") then
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
          
       case default
          write(*,*) "Unrecognized ",unitType ," unit '", unitString, "'"
             torusUnit = 1.d0
       end select

    else if(unitType == "wavelength" .or. unitType == "dust") then
       !Wavelength/dust units - TORUS uses Angstroms (A) and microns (um)
       select case(unitString)
       case("A")
          if(unitType == "wavelength") then
             torusUnit = 1.d0 
          elseif(unitType == "dust") then
             torusUnit = 1.d-4 
          end if
       

       case("nm")
          if(unitType == "wavelength") then
             torusUnit = 10.d0
          elseif(unitType == "dust") then
             torusUnit = 1.d-3
          end if
       
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

       case default
          write(*,*) "Unrecognized ",unitType ," unit '", unitString, "'"
             torusUnit = 1.d0
       end select

    else if(unitType == "angle") then
       !Angular units - TORUS uses radians
       select case(unitString)
       case("rad")
          torusUnit = 1.d0
          
       case("deg")
          torusUnit = degToRad
          
       case("arcmin")
          torusUnit = degToRad/60.d0
          
       case("arcsec")
          torusUnit = degToRad/3600.d0        

       case default
          write(*,*) "Unrecognized ",unitType ," unit '", unitString, "'"
             torusUnit = 1.d0
       end select
       
    else if (unitType == "mass") then
       !Mass units - torus uses g
       select case(unitString) 
       case("g")
          torusUnit = 1.d0
          
       case("mSol")
          torusUnit = mSol
          
       case("kg")
          torusUnit = 1.d3

       case default
          write(*,*) "Unrecognized ",unitType ," unit '", unitString, "'"
             torusUnit = 1.d0
       end select
       

    else if(unitType == "time") then
       !Time units - TORUS uses s
       select case(unitString)
       case("s")
          torusUnit = 1.d0
          
       case("yr")
          torusUnit = 1.d0/secsToYears
          
       case("kyr")
          torusUnit = 1.d0/(secsToYears*1.d3)
          
       case("Myr")
          torusUnit = 1.d0/(secsToYears*1.d6)

       case default
          write(*,*) "Unrecognized ",unitType ," unit '", unitString, "'"
             torusUnit = 1.d0
       end select


       else if(unitType == "temperature") then

          !Temperature units - TORUS uses K
          if(unitString == "K") then
             torusUnit = 1.d0
          else
             write(*,*) "Unrecognized ",unitType ," unit '", unitString, "'"
             torusUnit = 1.d0
          end if


       else if(unitType == "luminosity") then
          !Luminosity units - TORUS uses erg s^-1
          select case(unitString)
          case("ergsec")
             torusUnit = 1.d0
             
          case("lSol")
             torusUnit = lsol
             
          case default
             write(*,*) "Unrecognized ",unitType ," unit '", unitString, "'"
             torusUnit = 1.d0
          end select
       end if

       inputvalue = inputValue * torusUnit
    
  end subroutine convertToTorusUnits

end module
