module memory_mod
  use kind_mod
  use messages_mod
  use gridtype_mod
  

  implicit none

  public
  integer(kind=bigInt) :: globalMemory
  
  contains

    subroutine findTotalMemory(grid, totalMemory)
      type(GRIDTYPE) :: grid
      integer(kind=bigInt), intent(out) :: totalMemory

      totalMemory = 0 
      call findTotalMemoryPrivate(grid%octreeRoot, totalMemory)

      contains
        
        recursive subroutine findTotalMemoryPrivate(thisOctal, totalMemory)
      
        TYPE(OCTAL), pointer  :: thisOctal 
        TYPE(OCTAL), pointer  :: child
        integer(kind=bigInt) :: totalMemory
        integer:: i

        totalMemory = totalMemory + octalMemory(thisOctal)

        
        IF ( thisOctal%nChildren > 0 ) THEN
          ! call this subroutine recursively on each of its children
          DO i = 1, thisOctal%nChildren, 1
            child => thisOctal%child(i)
            CALL findTotalMemoryPrivate(child, totalMemory)
          END DO
        END IF

        
      END SUBROUTINE findTotalMemoryPrivate


    end subroutine findTotalMemory


    function octalMemory(thisOctal) result(tot)
      type(OCTAL), pointer :: thisOctal
      integer(kind=bigInt) tot

      tot = 0



#ifdef MEMCHECK
      tot = sizeof(thisOctal)



      tot = tot + sizeof(thisOctal%etaLine)
      tot = tot + sizeof(thisOctal%chiLine)
      tot = tot + sizeof(thisOctal%NH2)

      tot = tot + sizeof(thisOctal%oldFrac)
      tot = tot + sizeof(thisOctal%fixedTemperature)
      tot = tot + sizeof(thisOctal%dustType)
      tot = tot + sizeof(thisOctal%dusttypefraction)

      tot = tot + sizeof(thisOctal%meanIntensity)
      tot = tot + sizeof(thisOctal%diffusionApprox)
      tot = tot + sizeof(thisOctal%changed)
      tot = tot + sizeof(thisOctal%nDiffusion)
      tot = tot + sizeof(thisOctal%eDens)
      tot = tot + sizeof(thisOctal%diffusionCoeff)
      tot = tot + sizeof(thisOctal%oldeDens)
      tot = tot + sizeof(thisOctal%nDirectPhotons)
       
      tot = tot + sizeof(thisOctal%underSampled)
      tot = tot + sizeof(thisOctal%oldTemperature)
      tot = tot + sizeof(thisOctal%kappaRoss)
      tot = tot + sizeof(thisOctal%distanceGrid)
      tot = tot + sizeof(thisOctal%nCrossings)
      tot = tot + sizeof(thisOctal%nTot)
      tot = tot + sizeof(thisOctal%etaCont)
      tot = tot + sizeof(thisOctal%biasCont3D)
      tot = tot + sizeof(thisOctal%biasLine3D)
      tot = tot + sizeof(thisOctal%probDistLine)
      tot = tot + sizeof(thisOctal%probDistCont)

      tot = tot + sizeof(thisOctal%nTot)
      tot = tot + sizeof(thisOctal%ne)
      tot = tot + sizeof(thisOctal%N)
      tot = tot + sizeof(thisOctal%kappaAbs)
      tot = tot + sizeof(thisOctal%kappaSca)

      tot = tot + sizeof(thisOctal%molAbundance)
      tot = tot + sizeof(thisOctal%temperatureGas)
      tot = tot + sizeof(thisOctal%temperatureDust)
      tot = tot + sizeof(thisOctal%microturb)
      tot = tot + sizeof(thisOctal%molmicroturb)

      tot = tot + sizeof(thisOctal%nh)
      tot = tot + sizeof(thisOctal%nhi)
      tot = tot + sizeof(thisOctal%nhei)
      tot = tot + sizeof(thisOctal%nhii)
      tot = tot + sizeof(thisOctal%HHeating)
      tot = tot + sizeof(thisOctal%HeHeating)

      tot = tot + sizeof(thisOctal%ionFrac)
      tot = tot + sizeof(thisOctal%photoionCoeff)


      tot = tot + sizeof(thisOctal%uDens)
      tot = tot + sizeof(thisOctal%aDot)
      tot = tot + sizeof(thisOctal%distancegridaDot)
      tot = tot + sizeof(thisOctal%distanceGridPhotonFromGas)
      tot = tot + sizeof(thisOctal%distanceGridPhotonFromSource)
      tot = tot + sizeof(thisOctal%photonEnergyDensityFromGas)
      tot = tot + sizeof(thisOctal%photonEnergyDensityFromSource)
      tot = tot + sizeof(thisOctal%photonEnergyDensity)
      tot = tot + sizeof(thisOctal%oldphotonEnergyDensity)


      tot = tot + sizeof(thisOctal%q_i)
      tot = tot + sizeof(thisOctal%q_i_plus_1)
      tot = tot + sizeof(thisOctal%q_i_minus_1)
      tot = tot + sizeof(thisOctal%q_i_minus_2)

      tot = tot + sizeof(thisOctal%x_i)
      tot = tot + sizeof(thisOctal%x_i_plus_1)
      tot = tot + sizeof(thisOctal%x_i_minus_1)

      tot = tot + sizeof(thisOctal%u_interface)
      tot = tot + sizeof(thisOctal%u_i_plus_1)
      tot = tot + sizeof(thisOctal%u_i_minus_1)
      tot = tot + sizeof(thisOctal%flux_i)
      tot = tot + sizeof(thisOctal%flux_i_plus_1)
      tot = tot + sizeof(thisOctal%flux_i_minus_1)

      tot = tot + sizeof(thisOctal%phiLimit)
      tot = tot + sizeof(thisOctal%ghostCell)
      tot = tot + sizeof(thisOctal%feederCell)
      tot = tot + sizeof(thisOctal%edgeCell)
      tot = tot + sizeof(thisOctal%refinedLastTime)
      tot = tot + sizeof(thisOctal%pressure_i)
      tot = tot + sizeof(thisOctal%pressure_i_plus_1)
      tot = tot + sizeof(thisOctal%pressure_i_minus_1)
      tot = tot + sizeof(thisOctal%rhou)
      tot = tot + sizeof(thisOctal%rhov)
      tot = tot + sizeof(thisOctal%rhow)
      tot = tot + sizeof(thisOctal%rhoe)
      tot = tot + sizeof(thisOctal%energy)
      tot = tot + sizeof(thisOctal%phi_i)
      tot = tot + sizeof(thisOctal%phi_i_plus_1)
      tot = tot + sizeof(thisOctal%phi_i_minus_1)
      tot = tot + sizeof(thisOctal%rho_i_plus_1)
      tot = tot + sizeof(thisOctal%rho_i_minus_1)
      tot = tot + sizeof(thisOctal%boundaryCondition)
      tot = tot + sizeof(thisOctal%boundaryPartner)
      tot = tot + sizeof(thisOctal%gravboundaryPartner)
      tot = tot + sizeof(thisOctal%rLimit)
      tot = tot + sizeof(thisOctal%iEquationOfState)
      tot = tot + sizeof(thisOctal%gamma)
      tot = tot + sizeof(thisOctal%radiationMomentum)

      tot = tot + sizeof(thisOctal%molecularLevel)
      tot = tot + sizeof(thisOctal%molcellparam)
      tot = tot + sizeof(thisOctal%newmolecularLevel)
      tot = tot + sizeof(thisOctal%oldmolecularLevel)
      tot = tot + sizeof(thisOctal%oldestmolecularLevel)

#endif
    end function octalMemory


    subroutine reportMemory(iBytes)
      integer(kind=bigInt) :: iBytes, i
      character(len=80) :: message
      real :: f
      call writeInfo(" ",TRIVIAL)
      if ( iBytes < 1024 ) then
         i = iBytes 
         write(message,'(a,i4,a)') "Memory footprint of grid is ",i, " bytes"
         call writeInfo(message,TRIVIAL)
      else if ( iBytes < 1000000000 ) then
         f = real(iBytes)/1000000.
         write(message,'(a,f7.2,a)') "Memory footprint of grid is ",f, " Mbytes"
         call writeInfo(message,TRIVIAL)
      else
         f = real(iBytes)/1000000000.
         write(message,'(a,f7.2,a)') "Memory footprint of grid is ",f, " Gbytes"
         call writeInfo(message,TRIVIAL)
      endif
      call writeInfo(" ",TRIVIAL)

    end subroutine reportMemory

  end module memory_mod
