module physics_mod

  use constants_mod
  use messages_mod
  use mpi_global_mod
  use utils_mod
  use inputs_mod
  use timing
  use grid_mod
  use setupamr_mod
  implicit none

contains

  subroutine setupMicrophysics(grid)
    use input_variables, only : atomicPhysics, photoionPhysics, nAtom
    use modelatom_mod
    use source_mod
    type(GRIDTYPE) :: grid


    if (atomicPhysics) then
       if (associated(globalAtomArray)) deallocate(globalAtomArray)
       allocate(globalAtomArray(1:nAtom))
       call setupAtoms(nAtom, globalAtomArray)
    endif

  if (photoionPhysics) then
     call addIons(grid%ion, grid%nion)
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
    
    do iSource = 1, nSource
       
       source(iSource)%teff = sourceTeff(iSource)
       source(iSource)%mass = sourceMass(iSource)
       source(iSource)%radius = sourceRadius(iSource)
       source(iSource)%position = sourcePos(iSource)
       source(isource)%luminosity = fourPi * stefanBoltz * &
            (source(isource)%radius*1.d10)**2 * (source(isource)%teff)**4

       source(isource)%outsideGrid = .false.
       if (.not.inOctal(grid%octreeRoot, source(iSource)%position)) then
          source(isource)%outsideGrid = .true.
          source(isource)%distance = modulus(source(isource)%position)*1.d10
          source(isource)%luminosity = source(isource)%luminosity * (2.d0*grid%octreeRoot%subcellSize*1.d10)**2 / &
               (fourPi*source(isource)%distance**2)
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
    use input_variables, only : atomicPhysics, mie, photoionPhysics, photoionEquilibrium
    use input_variables, only : statisticalEquilibrium, nAtom, nDustType
    use cmf_mod, only : atomloop
    use photoionAMR_mod, only: photoionizationLoopAMR, ionizeGrid
    use photoion_mod, only : refineLambdaArray
    use source_mod, only : globalNsource, globalSourceArray
    use photoion_mod, only : refineLambdaArray

    real, allocatable :: xArray(:)
    integer :: nLambda
    type(PHASEMATRIX), pointer :: miePhase(:,:,:) => null()
    integer, parameter :: nMuMie = 20
    integer :: nCurrent, nt, i
    real(double) :: fac, logLamStart, logLamEnd, lamStart,lamEnd
    type(GRIDTYPE) :: grid

!    if (molecularPhysics) then
!       call setupMolecule(thisMolecule)
!    endif


!     if (dustPhysics.and.radiativeEquilibrium) then
!        call lucyRadiativeEquilibriumAMR(grid, miePhase, nDustType, nMuMie, nLambda, xArray, &
!             source, nSource, nLucy, massEnvelope, lucy_undersampled, finalPass=.true.)
!     endif

!     if (molecularPhysics.and.statisticalEquilibrium) then
!        call molecularLoop(grid, co)
!     endif
     
     if (atomicPhysics.and.statisticalEquilibrium) then
        call atomLoop(grid, nAtom, globalAtomArray, globalnsource, globalsourcearray)
     endif

     if (photoionPhysics.and.photoionEquilibrium) then 
        if (.not.grid%splitOverMPI) then
!           call photoIonizationloop(grid, source, nSource, nLambda, xArray, readlucy, writelucy, &
!             lucyfileNameout, lucyfileNamein)
        else
           nLambda = SIZE(globalsourcearray(1)%spectrum%flux)
           allocate(xArray(1:nLambda))
           xArray = globalsourcearray(1)%spectrum%lambda
           mie = .true.

           lamStart = 10.d0
           lamEnd = 1.d7
           logLamStart = log10(lamStart)
           logLamEnd = log10(lamEnd)
           nLambda = 1000
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

           grid%nLambda = nLambda
           if (associated(grid%lamArray)) deallocate(grid%lamArray)
           allocate(grid%lamArray(1:nLambda))
           grid%lamArray = xarray



           grid%nTempRossArray = 1000
           if (associated(grid%kappaRossArray)) deallocate(grid%kappaRossArray)
           if (associated(grid%tempRossArray)) deallocate(grid%tempRossArray)
           allocate(grid%kappaRossArray(nDustType,1:grid%nTempRossArray))
           allocate(grid%tempRossArray(1:grid%nTempRossArray))


           call  createDustCrossSectionPhaseMatrix(grid, xArray, nLambda, miePhase, nMuMie)
!           call refineLambdaArray(xArray, nLambda, grid)
!           call ionizeGrid(grid%octreeRoot)
           call photoIonizationloopAMR(grid, globalsourceArray, globalnSource, nLambda, xArray, .false., .false., &
             " ", " ", 5, 1.d30)
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

   subroutine setupGlobalSources(grid)     
     use source_mod, only : globalNsource, globalSourceArray
     use input_variables, only : inputNsource
     type(GRIDTYPE) :: grid


     call writeBanner("Source setup","-",TRIVIAL)
     if (inputNsource > 0) then
        globalnSource = inputNSource
        allocate(globalsourceArray(1:globalnSource))
        call setupSources(globalnSource, globalsourceArray, grid)
     endif
end subroutine setupGlobalSources


end module physics_mod
