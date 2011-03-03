module dimensionality_mod

  use kind_mod
  use messages_mod
  implicit none
  real(double) :: lengthToCodeUnits, lengthToPhysicalUnits
  real(double) :: timeToCodeUnits, timeToPhysicalUnits
  real(double) :: massToCodeUnits, massToPhysicalUnits
  real(double) :: chargeToCodeUnits, chargeToPhysicalUnits
  real(double) :: temperatureToCodeUnits, temperatureToPhysicalUnits

contains

  subroutine initializeCodeUnits()

    lengthToCodeUnits = 1.d0
    lengthToPhysicalUnits = 1.d0

    timeToCodeUnits = 1.d0
    timeToPhysicalUnits = 1.d0

    massToCodeUnits = 1.d0
    massToPhysicalUnits = 1.d0

    chargeToCodeUnits = 1.d0
    chargeToPhysicalUnits = 1.d0

    temperatureToCodeUnits = 1.d0
    temperatureToPhysicalUnits = 1.d0

  end subroutine initializeCodeUnits


  subroutine writeCodeUnits()
    
    if (writeoutput) then

       call writeBanner("Code Units","-")    
       write(*,'(a,1pe12.5)') "Length code unit: ",lengthToPhysicalUnits
       write(*,'(a,1pe12.5)') "Time code unit: ", timeToPhysicalUnits
       write(*,'(a,1pe12.5)') "Mass code unit: ", massToPhysicalUnits
       write(*,'(a,1pe12.5)') "Temperature code unit: ", temperatureToPhysicalUnits
       write(*,'(a,1pe12.5)') "Charge code unit: ", chargeToPhysicalUnits
       write(*,*) " "
    endif
  end subroutine writeCodeUnits

  subroutine setCodeUnit(length, time, mass, charge, temperature)
    real(double), optional :: length, time, mass, charge, temperature


    if (PRESENT(length)) then
       lengthtoCodeUnits = 1.d0 / length
       lengthtoPhysicalUnits = length
    endif

    if (PRESENT(time)) then
       timetoCodeUnits = 1.d0 / time
       timetoPhysicalUnits = time
    endif

    if (PRESENT(mass)) then
       masstoCodeUnits = 1.d0 / mass
       masstoPhysicalUnits = mass
    endif

    if (PRESENT(charge)) then
       chargetoCodeUnits = 1.d0 / charge
       chargetoPhysicalUnits = charge
    endif

    if (PRESENT(temperature)) then
       temperaturetoCodeUnits = 1.d0 / temperature
       temperaturetoPhysicalUnits = temperature
    endif

  end subroutine setCodeUnit

  real(double) function returnCodeUnitDensity(rho) result (codeRho)
    real(double) :: rho
    codeRho = rho * massToCodeUnits / (lengthToCodeUnits**3)
  end function returnCodeUnitDensity


  real(double) function returnPhysicalUnitDensity(rho) result (physicalRho)
    real(double) :: rho
    physicalRho = rho * massToPhysicalUnits / (lengthToPhysicalUnits**3)
  end function returnPhysicalUnitDensity

  real(double) function returnPhysicalUnitTime(time) result (physicalTime)
    real(double) :: time
    physicalTime = time * timeToPhysicalUnits
  end function returnPhysicalUnitTime

  real(double) function returnCodeUnitSpeed(speed) result (codeSpeed)
    real(double) :: speed
    codeSpeed = speed * lengthToCodeUnits / timeToCodeUnits
  end function returnCodeUnitSpeed


  real(double) function returnPhysicalUnitMass(mass) result (physicalMass)
    real(double) :: mass
    physicalmass = mass * massToPhysicalUnits
  end function returnPhysicalUnitMass

  real(double) function returnPhysicalUnitLength(length) result (physicalLength)
    real(double) :: length
    physicalLength = length * lengthToPhysicalUnits
  end function returnPhysicalUnitLength

  real(double) function returnPhysicalUnitSpeed(speed) result (physicalspeed)
    real(double) :: speed
    physicalSpeed = speed * lengthToPhysicalUnits/timeToPhysicalUnits
  end function returnPhysicalUnitSpeed



  real(double) function returnCodeUnitRhoV(rhov) result (coderhov)
    real(double) :: rhoV

    codeRhoV = rhoV * massToCodeUnits / (lengthToCodeUnits**2  * timeToCodeUnits)

  end function returnCodeUnitRhoV

  real(double) function returnCodeUnitRhoE(rhoe) result (coderhoE)
    real(double) :: rhoE

    codeRhoE = rhoE * massToCodeUnits**2 / (lengthToCodeUnits * timeToCodeUnits**2)

  end function returnCodeUnitRhoE

  real(double) function returnCodeUnitPressure(pressure) result (codePressure)
    real(double) :: pressure

    codePressure = pressure * massToCodeUnits / (lengthToCodeUnits*timeToCodeUnits**2)

  end function returnCodeUnitPressure

  real(double) function returnCodeUnitEnergy(energy) result (codeEnergy)
    real(double) :: energy

    codeEnergy = energy * masstoCodeUnits * lengthtoCodeUnits**2 / timeToCodeUnits**2

  end function returnCodeUnitEnergy

  real(double) function returnCodeUnitLength(length) result (codeLength)
    real(double) :: length

    codeLength = length * lengthToCodeUnits

  end function returnCodeUnitLength


  real(double) function returnCodeUnitTime(time) result (codeTime)
    real(double) :: time

    codetime = time * timetoCodeUnits

  end function returnCodeUnitTime

  

end module dimensionality_mod
