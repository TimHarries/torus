module memory_mod
  use kind_mod
  use vector_mod
  use messages_mod
  use gridtype_mod

  implicit none

#ifdef MEMCHECK
  interface mySizeOf
     module procedure mySizeOfDouble
     module procedure mySizeOfDouble2D
     module procedure mySizeOfDouble3D
     module procedure mySizeOfReal
     module procedure mySizeOfInteger
     module procedure mySizeOfLogical
     module procedure mySizeOfVector
  end interface
#endif  

  public
  integer(kind=bigInt) :: globalMemoryFootprint
  
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
      integer :: i
      tot = 0
      i = thisOctal%nChildren


#ifdef MEMCHECK
      tot = sizeOf(thisOctal)



      tot = tot + mySizeOf(thisOctal%etaLine)
      tot = tot + mySizeOf(thisOctal%chiLine)
      tot = tot + mySizeOf(thisOctal%NH2)

      tot = tot + mySizeOf(thisOctal%oldFrac)
      tot = tot + mySizeOf(thisOctal%fixedTemperature)
      tot = tot + mySizeOf(thisOctal%dustType)
      tot = tot + mySizeOf(thisOctal%dusttypefraction)

      tot = tot + mySizeOf(thisOctal%meanIntensity)
      tot = tot + mySizeOf(thisOctal%diffusionApprox)
      tot = tot + mySizeOf(thisOctal%changed)
      tot = tot + mySizeOf(thisOctal%nDiffusion)
      tot = tot + mySizeOf(thisOctal%eDens)
      tot = tot + mySizeOf(thisOctal%diffusionCoeff)
      tot = tot + mySizeOf(thisOctal%oldeDens)
      tot = tot + mySizeOf(thisOctal%nDirectPhotons)
       
      tot = tot + mySizeOf(thisOctal%underSampled)
      tot = tot + mySizeOf(thisOctal%oldTemperature)
      tot = tot + mySizeOf(thisOctal%kappaRoss)
      tot = tot + mySizeOf(thisOctal%distanceGrid)
      tot = tot + mySizeOf(thisOctal%nCrossings)
      tot = tot + mySizeOf(thisOctal%nTot)
      tot = tot + mySizeOf(thisOctal%etaCont)
      tot = tot + mySizeOf(thisOctal%biasCont3D)
      tot = tot + mySizeOf(thisOctal%biasLine3D)
      tot = tot + mySizeOf(thisOctal%probDistLine)
      tot = tot + mySizeOf(thisOctal%probDistCont)

      tot = tot + mySizeOf(thisOctal%nTot)
      tot = tot + mySizeOf(thisOctal%ne)
      tot = tot + mySizeOf(thisOctal%N)
      tot = tot + mySizeOf(thisOctal%kappaAbs)
      tot = tot + mySizeOf(thisOctal%kappaSca)

      tot = tot + mySizeOf(thisOctal%molAbundance)
      tot = tot + mySizeOf(thisOctal%temperatureGas)
      tot = tot + mySizeOf(thisOctal%temperatureDust)
      tot = tot + mySizeOf(thisOctal%microturb)
      tot = tot + mySizeOf(thisOctal%molmicroturb)

      tot = tot + mySizeOf(thisOctal%nh)
      tot = tot + mySizeOf(thisOctal%nhi)
      tot = tot + mySizeOf(thisOctal%nhei)
      tot = tot + mySizeOf(thisOctal%nhii)
      tot = tot + mySizeOf(thisOctal%HHeating)
      tot = tot + mySizeOf(thisOctal%HeHeating)

      tot = tot + mySizeOf(thisOctal%ionFrac)
      tot = tot + mySizeOf(thisOctal%photoionCoeff)


      tot = tot + mySizeOf(thisOctal%uDens)
      tot = tot + mySizeOf(thisOctal%aDot)
      tot = tot + mySizeOf(thisOctal%distancegridaDot)
      tot = tot + mySizeOf(thisOctal%distanceGridPhotonFromGas)
      tot = tot + mySizeOf(thisOctal%distanceGridPhotonFromSource)
      tot = tot + mySizeOf(thisOctal%photonEnergyDensityFromGas)
      tot = tot + mySizeOf(thisOctal%photonEnergyDensityFromSource)
      tot = tot + mySizeOf(thisOctal%photonEnergyDensity)
      tot = tot + mySizeOf(thisOctal%oldphotonEnergyDensity)


      tot = tot + mySizeOf(thisOctal%q_i)
      tot = tot + mySizeOf(thisOctal%q_i_plus_1)
      tot = tot + mySizeOf(thisOctal%q_i_minus_1)
      tot = tot + mySizeOf(thisOctal%q_i_minus_2)

      tot = tot + mySizeOf(thisOctal%x_i)
      tot = tot + mySizeOf(thisOctal%x_i_plus_1)
      tot = tot + mySizeOf(thisOctal%x_i_minus_1)

      tot = tot + mySizeOf(thisOctal%u_interface)
      tot = tot + mySizeOf(thisOctal%u_i_plus_1)
      tot = tot + mySizeOf(thisOctal%u_i_minus_1)
      tot = tot + mySizeOf(thisOctal%flux_i)
      tot = tot + mySizeOf(thisOctal%flux_i_plus_1)
      tot = tot + mySizeOf(thisOctal%flux_i_minus_1)

      tot = tot + mySizeOf(thisOctal%phiLimit)
      tot = tot + mySizeOf(thisOctal%ghostCell)
      tot = tot + mySizeOf(thisOctal%feederCell)
      tot = tot + mySizeOf(thisOctal%edgeCell)
      tot = tot + mySizeOf(thisOctal%refinedLastTime)
      tot = tot + mySizeOf(thisOctal%pressure_i)
      tot = tot + mySizeOf(thisOctal%pressure_i_plus_1)
      tot = tot + mySizeOf(thisOctal%pressure_i_minus_1)
      tot = tot + mySizeOf(thisOctal%rhou)
      tot = tot + mySizeOf(thisOctal%rhov)
      tot = tot + mySizeOf(thisOctal%rhow)
      tot = tot + mySizeOf(thisOctal%rhoe)
      tot = tot + mySizeOf(thisOctal%energy)
      tot = tot + mySizeOf(thisOctal%phi_i)
      tot = tot + mySizeOf(thisOctal%phi_i_plus_1)
      tot = tot + mySizeOf(thisOctal%phi_i_minus_1)
      tot = tot + mySizeOf(thisOctal%rho_i_plus_1)
      tot = tot + mySizeOf(thisOctal%rho_i_minus_1)
      tot = tot + mySizeOf(thisOctal%boundaryCondition)
      tot = tot + mySizeOf(thisOctal%boundaryPartner)
      tot = tot + mySizeOf(thisOctal%gravboundaryPartner)
      tot = tot + mySizeOf(thisOctal%rLimit)
      tot = tot + mySizeOf(thisOctal%iEquationOfState)
      tot = tot + mySizeOf(thisOctal%gamma)
      tot = tot + mySizeOf(thisOctal%radiationMomentum)

      tot = tot + mySizeOf(thisOctal%molecularLevel)
      tot = tot + mySizeOf(thisOctal%molcellparam)
      tot = tot + mySizeOf(thisOctal%newmolecularLevel)
      tot = tot + mySizeOf(thisOctal%oldmolecularLevel)
      tot = tot + mySizeOf(thisOctal%oldestmolecularLevel)

#endif

    end function octalMemory


    subroutine reportMemory(iBytes)
      integer(kind=bigInt) :: iBytes, i
      character(len=80) :: message
      real :: f
      call writeInfo(" ",TRIVIAL)
      if ( iBytes < 1024 ) then
         i = iBytes 
         write(message,'(a,i8,a)') "Memory footprint of grid is ",i, " bytes"
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


    function  humanReadableMemory(iBytes) result (cString)
      integer(kind=bigInt) :: iBytes, i
      character(len=20) :: cString
      real :: f
      if ( iBytes < 1024 ) then
         i = iBytes 
         write(cstring,'(i8,a)') i, " bytes"
      else if ( iBytes < 1000000000 ) then
         f = real(iBytes)/1000000.
         write(cString,'(f7.2,a)') f, " Mbytes"
      else
         f = real(iBytes)/1000000000.
         write(cString,'(f7.2,a)') f, " Gbytes"
      endif

    end function humanReadableMemory

#ifdef MEMCHECK
    function mySizeOfInteger(p) result (i)
      integer, pointer :: p(:)
      integer :: i
      i = 0
      if (associated(p)) i = sizeof(p)
    end function mySizeOfInteger

    function mySizeOfReal(p) result (i)
      real, pointer :: p(:)
      integer :: i
      i = 0
      if (associated(p)) i = sizeof(p)
    end function mySizeOfReal

    function mySizeOfDouble(p) result (i)
      real(double), pointer :: p(:)
      integer :: i
      i = 0
      if (associated(p)) i = sizeof(p)
    end function mySizeOfDouble

    function mySizeOfDouble2d(p) result (i)
      real(double), pointer :: p(:,:)
      integer :: i
      i = 0
      if (associated(p)) i = sizeof(p)
    end function mySizeOfDouble2d

    function mySizeOfDouble3d(p) result (i)
      real(double), pointer :: p(:,:,:)
      integer :: i
      i = 0
      if (associated(p)) i = sizeof(p)
    end function mySizeOfDouble3d

    function mySizeOfLogical(p) result (i)
      logical, pointer :: p(:)
      integer :: i
      i = 0
      if (associated(p)) i = sizeof(p)
    end function mySizeOfLogical

    function mySizeOfVector(p) result (i)
      type(VECTOR), pointer :: p(:)
      integer :: i
      i = 0
      if (associated(p)) i = sizeof(p)
    end function mySizeOfVector
#endif

  end module memory_mod
