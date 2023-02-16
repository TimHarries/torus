module units_mod
use messages_mod
use constants_mod
implicit none

type UNITTYPE
   character(len=20) :: unitName
   character(len=10) :: unitDimension
   character(len=80) :: longUnitName
   real(double) :: unitInCGS
end type UNITTYPE

integer :: nUnits
type(UNITTYPE) :: unitList(100)


contains

  subroutine addUnit(newUnit)
    type(UNITTYPE) :: newUnit
    nUnits = nUnits + 1
    unitList(nUnits) = newUnit
  end subroutine addUnit


  subroutine setupUnits()
    nUnits = 0

    call addUnit(UNITTYPE("codeunits","L","10^10 centimetre",1.d10))


! lengths

    call addUnit(UNITTYPE("cm","L","centimetre",1.d0))
    call addUnit(UNITTYPE("m","L","metre",100.d0))
    call addUnit(UNITTYPE("km","L","kilometre",1e5))
    call addUnit(UNITTYPE("mm","L","millimetre",0.1d0))
    call addUnit(UNITTYPE("micron","L","micrometre",1.d-4))
    call addUnit(UNITTYPE("angstrom","L","angstrom",1.d-8))
    call addUnit(UNITTYPE("rsol","L","solar radii",rsol))
    call addUnit(UNITTYPE("au","L","AU",autocm))
    call addUnit(UNITTYPE("pc","L","parsec",pctocm))

! masses

    call addUnit(UNITTYPE("gram","M","gram",1.d0))
    call addUnit(UNITTYPE("msol","M","solar masses",mSol))
    call addUnit(UNITTYPE("mmoon","M","lunar masses",mMoon))
    call addUnit(UNITTYPE("mhydrogen","M","Hydrogen atom mass",mHydrogen))
    call addUnit(UNITTYPE("mearth","M","earth masses",mEarth))

! Times

    call addUnit(UNITTYPE("year","T","year",yearsTosecs))

! Speeds

    call addUnit(UNITTYPE("c","L/T","Speed of light",cSpeed))
    call addUnit(UNITTYPE("km/s","L/T","km/s",1.d5))

! Temperatures
    call addUnit(UNITTYPE("K","L/T","kelvin",1.d0))
    
  end subroutine setupUnits

  logical function correctDimensions(thisUnit, testDimensions)
    type(UNITTYPE) :: thisUnit
    character(len=*) :: testDimensions

    correctDimensions = .true.
    if (thisUnit%unitDimension .ne. testDimensions) then
       correctDimensions = .false.
    endif
  end function correctDimensions
 
  
  integer function unitNumber(unitString)
    character(len=*) unitString
    logical :: found
    integer :: i
    found = .false.
    unitNumber = 0
    do i = 1, nUnits
       if (unitString == unitList(i)%unitName) then
          found = .true.
          unitNumber = i
          cycle
       endif
    enddo
    if (.not.found) then
       call writeFatal("Unit "//trim(unitString)//" not recognised")
    endif
  end function unitNumber
  
  subroutine convertToTorusUnits(unitString, unitType, inputValue)

    implicit none
    
    character(len=*) :: unitString
    character(len=*) :: unitType
    real(double) :: torusUnit
    real(double), intent(inout) :: inputValue
   

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
!          write(*,*) "Unrecognized ",unitType ," unit '", unitString, "'"
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
!          write(*,*) "Unrecognized ",unitType ," unit '", unitString, "'"
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
!          write(*,*) "Unrecognized ",unitType ," unit '", unitString, "'"
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
!          write(*,*) "Unrecognized ",unitType ," unit '", unitString, "'"
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
!          write(*,*) "Unrecognized ",unitType ," unit '", unitString, "'"
             torusUnit = 1.d0
       end select


       else if(unitType == "temperature") then

          !Temperature units - TORUS uses K
          if(unitString == "K") then
             torusUnit = 1.d0
          else
 !            write(*,*) "Unrecognized ",unitType ," unit '", unitString, "'"
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
!             write(*,*) "Unrecognized ",unitType ," unit '", unitString, "'"
             torusUnit = 1.d0
          end select


       else if(unitType == "velocity") then
          select case(unitString)
          
          case("cms")
             torusUnit = 1.d0
          case("ms")
             torusUnit = 1.d2
          case("kms")
             torusUnit = 1.d5
          case("c")
             torusUnit = cspeed
          case default
             torusUnit = 1.d0
          end select

       else if(unitType == "density") then
          select case(unitString)
            
          case default
             torusUnit = 1.d0
          end select
       end if
    
       inputvalue = inputValue * torusUnit
           
  end subroutine convertToTorusUnits

end module
