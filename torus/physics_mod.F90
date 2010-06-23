module physics_mod

  use constants_mod
  use messages_mod
  use mpi_global_mod
  use utils_mod
  use timing
  use grid_mod
  use setupamr_mod
  use vector_mod

  implicit none


contains

  subroutine setupMicrophysics(grid)
    use input_variables, only : atomicPhysics, photoionPhysics, nAtom, photoionization, molecular
    use input_variables, only : molecularPhysics, moleculeFile
    use molecular_mod
    use modelatom_mod
    use source_mod
    type(GRIDTYPE) :: grid

    if (molecularPhysics) then
       molecular = .true.
       call readMolecule(globalMolecule, moleculefile)
    endif


    if (atomicPhysics) then
       if (associated(globalAtomArray)) deallocate(globalAtomArray)
       allocate(globalAtomArray(1:nAtom))
       call setupAtoms(nAtom, globalAtomArray)
    endif

  if (photoionPhysics) then
     call addIons(grid%ion, grid%nion)
     photoionization = .true.
  endif



  end subroutine setupMicrophysics

  subroutine setupAtoms(nAtom, atomarray)
    use modelatom_mod
    use input_variables, only : atomFilename
    integer :: nAtom, iAtom
    type(MODELATOM) :: atomArray(:)

    do iatom = 1, nAtom
       call readAtom(atomArray(iatom), atomFilename(iatom))
    end do
  end subroutine setupAtoms

  subroutine setupSources(nSource, source, grid)
    use spectrum_mod
    use surface_mod
    use source_mod
    use input_variables
    type(GRIDTYPE) :: grid
    integer :: nSource, iSource
    type(SOURCETYPE) :: source(:)
    logical :: ok
    real(double) :: distToEdge
    
    do iSource = 1, nSource
       
       source(iSource)%teff = sourceTeff(iSource)
       source(iSource)%mass = sourceMass(iSource)
       source(iSource)%radius = sourceRadius(iSource)
       source(iSource)%position = sourcePos(iSource)
       source(isource)%luminosity = fourPi * stefanBoltz * &
            (source(isource)%radius*1.d10)**2 * (source(isource)%teff)**4

       source(isource)%outsideGrid = .false.
       source(isource)%onEdge      = .false. 
       distToEdge = source(iSource)%position%z - (grid%octreeRoot%centre%z - grid%octreeRoot%subcellSize)
       if (.not.inOctal(grid%octreeRoot, source(iSource)%position)) then
          source(isource)%outsideGrid = .true.
          source(isource)%distance = modulus(source(isource)%position)*1.d10
          source(isource)%luminosity = source(isource)%luminosity * (2.d0*grid%octreeRoot%subcellSize*1.d10)**2 / &
               (fourPi*source(isource)%distance**2)
       else if ( distToEdge < grid%halfSmallestSubcell) then
          source%onEdge = .true.
          ! only half the photons will end up on the grid
          source(isource)%luminosity = source(isource)%luminosity / 2.0
       endif

       if (inputcontfluxfile(isource) /= "blackbody") then
          call readSpectrum(source(isource)%spectrum, inputcontfluxfile(isource), ok)
       else
          call fillSpectrumBB(source(isource)%spectrum, source(isource)%teff, 10.d0, 1000.d4,1000)
       endif
       call normalizedSpectrum(source(isource)%spectrum)
       lamStart = 10.d0
       lamEnd = 1000.d4
       nlambda = 1000
       call buildSphere(source(isource)%position, dble(source(isource)%radius), &
            source(isource)%surface, 100, inputcontFluxFile(isource), source(isource)%teff)

    end do
  end subroutine setupSources

  subroutine doPhysics(grid)
    use phasematrix_mod
    use dust_mod
    use modelatom_mod, only : globalAtomArray
    use input_variables, only : atomicPhysics, photoionPhysics, photoionEquilibrium
    use input_variables, only : dustPhysics, lowmemory, radiativeEquilibrium
    use input_variables, only : statisticalEquilibrium, nAtom, nDustType, nLucy, &
         lucy_undersampled, molecularPhysics
    use input_variables, only : useDust, realDust, readlucy, writelucy
    use input_variables, only : lucyfilenameOut, lucyFilenamein
    use cmf_mod, only : atomloop
    use photoionAMR_mod, only: photoionizationLoopAMR, ionizeGrid
    use photoion_mod, only : refineLambdaArray, photoionizationLoop
    use source_mod, only : globalNsource, globalSourceArray
    use molecular_mod, only : molecularLoop, globalMolecule
    use lucy_mod, only : lucyRadiativeEquilibriumAMR
#ifdef MPI
    use input_variables, only : hydrovelocityconv
    use mpi_amr_mod, only : fillVelocityCornersFromHydro
#endif
    use amr_mod, only : hydroVelocityConvert
    use setupamr_mod, only: doSmoothOnTau
    real, pointer :: xArray(:) => null()
    integer :: nLambda 
    type(PHASEMATRIX), pointer :: miePhase(:,:,:) => null()
    integer, parameter :: nMuMie = 20
    type(GRIDTYPE) :: grid
    real :: massEnvelope


     if (dustPhysics.and.radiativeEquilibrium) then
        call setupXarray(grid, xarray, nLambda)
        call setupDust(grid, xArray, nLambda, miePhase, nMumie)
!        call fillDustUniform(grid, grid%octreeRoot)
        call doSmoothOnTau(grid)
        call lucyRadiativeEquilibriumAMR(grid, miePhase, nDustType, nMuMie, nLambda, xArray, &
             globalsourcearray, globalnSource, nLucy, massEnvelope, lucy_undersampled, finalPass=.true.)
     endif

     if (molecularPhysics.and.statisticalEquilibrium) then
#ifdef MPI
        if (grid%splitOverMPI.and.hydrovelocityconv) then
           call setallUnchanged(grid%octreeRoot)
           call hydroVelocityConvert(grid%octreeRoot)
           if (myrankGlobal /= 0) then
              call fillVelocityCornersFromHydro(grid)
           endif
        endif
#endif
        if (dustPhysics) then
           call setupXarray(grid, xarray, nLambda)
           call setupDust(grid, xArray, nLambda, miePhase, nMumie)
           usedust = .true.
           realdust = .true.
        endif
        lowMemory = .false.
        call molecularLoop(grid, globalMolecule)
     endif
     
     if (atomicPhysics.and.statisticalEquilibrium) then
        call atomLoop(grid, nAtom, globalAtomArray, globalnsource, globalsourcearray)
     endif

     if (photoionPhysics.and.photoionEquilibrium) then 

        call setupXarray(grid, xArray, nLambda)
        if (dustPhysics) call setupDust(grid, xArray, nLambda, miePhase, nMumie)

        if (.not.grid%splitOverMPI) then
           call photoIonizationloop(grid, globalsourceArray, globalnSource, nLambda, xArray, readlucy, writelucy, &
             lucyfileNameout, lucyfileNamein)
        else

           call photoIonizationloopAMR(grid, globalsourceArray, globalnSource, nLambda, xArray, .false., .false., &
             " ", " ", 5, 1.d30, sublimate=.false.)
        endif

     end if

!     if (hydrodynamics) then
!        call hydrodynamics(grid)
!     endif

!     if (radiationHydrodynamics) then
!        call radiationHydro(grid, source, nSource, nLambda, xArray, readlucy, writelucy, &
!             lucyfilenameout, lucyfilenamein)
!        call torus_mpi_barrier
!     endif

   end subroutine doPhysics

   subroutine setupXarray(grid, xArray, nLambda, lamMin, lamMax, wavLin)
     use input_variables, only : photoionPhysics, dustPhysics, molecularPhysics
     use photoion_mod, only : refineLambdaArray
     type(GRIDTYPE) :: grid
     real, pointer :: xArray(:)
     integer :: nLambda
     real, optional, intent(in) :: lamMin, lamMax
     logical, optional, intent(in) :: wavLin
     integer :: nCurrent, nt, i
     real(double) :: fac, logLamStart, logLamEnd, lamStart,lamEnd
     if (associated(xarray)) deallocate(xarray)
     if (molecularPhysics) then
        nLambda = 100
        allocate(xarray(1:nLambda))
        lamStart = 1200.d0
        lamEnd = 1.d7
        logLamStart = log10(lamStart)
        logLamEnd = log10(lamEnd)
        do i = 1, nlambda
           fac = logLamStart + real(i-1)/real(nLambda-1)*(logLamEnd - logLamStart)
           fac = 10.**fac
           xArray(i) = fac
        enddo
     endif

     if (dustPhysics) then

        nLambda = 1000
        allocate(xarray(1:nLambda))

        if ( present(lamMin) ) then 
           lamStart = lamMin
        else
           lamStart = 10.d0
        end if
        if ( present(lamMax) ) then 
           lamEnd = lamMax
        else
           lamEnd = 1.d7
        end if

        if ( present(wavLin) ) then 
           if (wavLin) then 
              call setupLinSpacing
           else
              call setupLogSpacing
           endif
        else
           call setupLogSpacing
        end if
     endif

     if (photoionPhysics) then
        nLambda = 1000
        allocate(xarray(1:nLambda))
        lamStart = 10.d0
        lamEnd = 1.d7
        logLamStart = log10(lamStart)
        logLamEnd = log10(lamEnd)
        xArray(1) = lamStart
        xArray(2) = lamEnd
        nCurrent = 2
        call refineLambdaArray(xArray, nCurrent, grid)
        nt = nLambda - nCurrent
        do i = 1, nt
           fac = logLamStart + real(i)/real(nt+1)*(logLamEnd - logLamStart)
           fac = 10.**fac
           nCurrent=nCurrent + 1
           xArray(nCurrent) = fac
           call sort(nCurrent, xArray)
        enddo
     endif




     grid%nLambda = nLambda
     if (associated(grid%lamArray)) deallocate(grid%lamArray)
     allocate(grid%lamArray(1:nLambda))
     grid%lamArray = xarray

     contains

       subroutine setupLogSpacing

         logLamStart = log10(lamStart)
         logLamEnd = log10(lamEnd)
         do i = 1, nlambda
            fac = logLamStart + real(i-1)/real(nLambda-1)*(logLamEnd - logLamStart)
            fac = 10.**fac
            xArray(i) = fac
        enddo
        
       end subroutine setupLogSpacing

       subroutine setupLinSpacing

         do i = 1, nlambda
            xArray(i) = LamStart + real(i-1)/real(nLambda-1)*(LamEnd - LamStart)
        enddo
        
       end subroutine setupLinSpacing

   end subroutine setupXarray

   subroutine setupGlobalSources(grid)     
     use parallel_mod
     use starburst_mod
     use source_mod, only : globalNsource, globalSourceArray
     use input_variables, only : inputNsource
     type(GRIDTYPE) :: grid
     real(double) :: coreContinuumFlux, lAccretion
     real :: fAccretion


     if (associated(globalsourceArray)) then
        deallocate(globalSourceArray)
     endif

     if (inputNsource > 0 ) call writeBanner("Source setup","-",TRIVIAL)
     if (inputNsource > 0) then
        globalnSource = inputNSource
        allocate(globalsourceArray(1:globalnSource))
        call setupSources(globalnSource, globalsourceArray, grid)
     endif

     if (grid%geometry == "starburst") then
#ifdef MPI
        call sync_random_seed()
#endif
        allocate(globalsourcearray(1:10000))
        globalsourceArray(:)%outsideGrid = .false.
        globalnSource = 0
        call createSources(globalnSource,globalsourcearray, "instantaneous", 1.d6, 1.d3, 1.d0)
        call init_random_seed()
    endif
    

    if (grid%geometry(1:6) == "ttauri") then
       coreContinuumFlux = 0.d0
       call buildSphere(globalsourceArray(1)%position, globalSourceArray(1)%radius, &
            globalsourcearray(1)%surface, 1000, "blackbody", &
            globalsourcearray(1)%teff)
       call genericAccretionSurface(globalsourcearray(1)%surface, grid, 1.e16, coreContinuumFlux,fAccretion, lAccretion) 
       globalsourcearray(1)%luminosity = globalsourcearray(1)%luminosity + lAccretion
    endif



end subroutine setupGlobalSources

subroutine setupDust(grid, xArray, nLambda, miePhase, nMumie)
  use input_variables, only : mie, nDustType
  use phasematrix_mod
  use dust_mod
  type(GRIDTYPE) :: grid
  real, pointer :: xArray(:)
  integer :: nLambda
  type(PHASEMATRIX), pointer :: miePhase(:,:,:)
  integer :: nMuMie

  grid%oneKappa = .true.
  mie = .true.
  grid%nTempRossArray = 1000
  if (associated(grid%kappaRossArray)) deallocate(grid%kappaRossArray)
  if (associated(grid%tempRossArray)) deallocate(grid%tempRossArray)
  allocate(grid%kappaRossArray(nDustType,1:grid%nTempRossArray))
  allocate(grid%tempRossArray(1:grid%nTempRossArray))


  call  createDustCrossSectionPhaseMatrix(grid, xArray, nLambda, miePhase, nMuMie)
  call allocateMemoryForDust(grid%octreeRoot)
end subroutine setupDust

end module physics_mod
