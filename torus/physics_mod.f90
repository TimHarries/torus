module physics_mod

  use constants_mod
  use messages_mod
  use mpi_global_mod
  use utils_mod
  use inputs_mod
  use timing
  use grid_mod
  use setupamr_mod

contains

  subroutine setupMicrophysics()
    use input_variables, only : atomicPhysics
    use input_variables, only : inputNsource
    use modelatom_mod
    use source_mod


    if (inputNsource > 0) then
       globalnSource = inputNSource
       allocate(globalsourceArray(1:globalnSource))
       call setupSources(globalnSource, globalsourceArray)
    endif

    if (atomicPhysics) then
       if (associated(globalAtomArray)) deallocate(globalAtomArray)
       allocate(globalAtomArray(1:nAtom))
       call setupAtoms(nAtom, globalAtomArray)
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

  subroutine setupSources(nSource, source)
    use spectrum_mod
    use surface_mod
    use input_variables, only : lamStart, lamEnd, nLambda
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
    use modelatom_mod, only : globalAtomArray
    use source_mod, only : globalNsource, globalSourceArray
    use input_variables, only : atomicPhysics
    use input_variables, only : statisticalEquilibrium
    use cmf_mod, only : atomloop
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

!     if (photoionPhysics.and.photoionEquilibrium) then 
!        call photoIonizationloop(grid, source, nSource, nLambda, xArray, readlucy, writelucy, &
!             lucyfileNameout, lucyfileNamein)
!     end if

!     if (hydrodynamics) then
!        call hydrodynamics(grid)
!     endif

!     if (radiationHydrodynamics) then
!        call radiationHydro(grid, source, nSource, nLambda, xArray, readlucy, writelucy, &
!             lucyfilenameout, lucyfilenamein)
!        call torus_mpi_barrier
!     endif

   end subroutine doPhysics



end module physics_mod
