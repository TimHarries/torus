module physics_mod
  
  use constants_mod
  use messages_mod
  use mpi_global_mod
  use utils_mod
  use timing
  use grid_mod
  use setupamr_mod
  use vector_mod
  use parallel_mod
  use random_mod
  implicit none


contains

  subroutine setupMicrophysics(grid)
#ifdef CMFATOM
    use inputs_mod, only : atomicPhysics, nAtom
    use modelatom_mod
#endif
#ifdef MOLECULAR
    use inputs_mod, only : molecularPhysics, moleculeFile, molecular
    use molecular_mod, only: readMolecule, globalMolecule
#endif
#ifdef PHOTOION
    use ion_mod, only: addIons, globalIonArray, nGlobalIon
    use inputs_mod, only : photoionization, photoionPhysics, usemetals, hOnly, usexraymetals
!    use inputs_mod, only : startFromNeutral
#endif


    type(GRIDTYPE) :: grid

#ifdef MOLECULAR
    if (molecularPhysics) then
       molecular = .true.
       call readMolecule(globalMolecule, moleculefile)
    endif
#endif

#ifdef CMFATOM
    if (atomicPhysics) then
       if (associated(globalAtomArray)) deallocate(globalAtomArray)
       allocate(globalAtomArray(1:nAtom))
       call setupAtoms(nAtom, globalAtomArray)
    endif
#endif

#ifdef PHOTOION
  if (photoionPhysics) then
     call addIons(grid%ion, grid%nion, usemetals, hOnly, usexraymetals)
     call addIons(globalIonArray, nGlobalIon, usemetals, hOnly, usexraymetals)
     photoionization = .true.
  endif
#endif


  end subroutine setupMicrophysics

#ifdef CMFATOM
  subroutine setupAtoms(nAtom, atomarray)
    use modelatom_mod
    use inputs_mod, only : atomFilename
    integer :: nAtom, iAtom
    type(MODELATOM) :: atomArray(:)

    do iatom = 1, nAtom
       call readAtom(atomArray(iatom), atomFilename(iatom))
!       call stripAtomLevels(atomArray(iAtom), 5)
    end do

  end subroutine setupAtoms
#endif

  subroutine setupSources(nSource, source, grid)
    use spectrum_mod
    use surface_mod
    use source_mod
    use inputs_mod
    use vtk_mod
    type(GRIDTYPE) :: grid
    integer :: nSource, iSource
    type(SOURCETYPE) :: source(:)
    logical :: ok
    real(double) :: distToEdge, fac, sumSurfaceLuminosity
    character(len=120) :: message
    integer :: i
    i = 0

    if (readSources) then

#ifdef MPI
       do i = 0, nThreadsGlobal-1
          if (myrankWorldGlobal == i) then
#endif
             call readSourceArray(nsource, source, sourceFilename)

#ifdef MPI
          endif
          call torus_mpi_barrier
       enddo
#endif
       call writeSourceList(source, nSource)
    else

       do iSource = 1, nSource

          if (stellarSource(isource)) then


             if (doNbodyOnly) then
                source(iSource)%mass = sourceMass(iSource)
                source(iSource)%position = sourcePos(iSource)
                source(iSource)%velocity = sourceVel(iSource)
             else
                source(isource)%stellar = .true.
                source(isource)%diffuse = .false.

                source(isource)%limbDark = 0.d0
                source(iSource)%teff = sourceTeff(iSource)
                source(iSource)%mass = sourceMass(iSource)
                source(iSource)%mdot = sourceMdot(iSource)
                source(iSource)%radius = sourceRadius(iSource)
                source(iSource)%position = sourcePos(iSource)
                source(iSource)%velocity = sourceVel(iSource)
                source(iSource)%accretionRadius = 1.d10*accretionRadius*smallestCellSize
                source(isource)%luminosity = fourPi * stefanBoltz * &
                     (source(isource)%radius*1.d10)**2 * (source(isource)%teff)**4
                source(iSource)%prob = sourceProb(iSource)
                source(iSource)%pointSource = pointSourceArray(iSource)
                source(isource)%time = 0.d0
                source(isource)%outsideGrid = .false.
                source(isource)%onEdge      = .false. 
                source(isource)%onCorner    = .false. 
                !       distToEdge = abs(source(iSource)%position%z) - abs((grid%octreeRoot%centre%z - grid%octreeRoot%subcellSize))
                !       if (.not.inOctal(grid%octreeRoot, source(iSource)%position)) then
                !          source(isource)%outsideGrid = .true.
                !          source(isource)%distance = modulus(source(isource)%position)*1.d10
                !          source(isource)%luminosity = source(isource)%luminosity * (2.d0*grid%octreeRoot%subcellSize*1.d10)**2 / &
                !               (fourPi*source(isource)%distance**2)
                !       else if ( distToEdge < grid%halfSmallestSubcell) then
                !          source(iSource)%onEdge = .true.
                !          source(iSource)%onCorner = .false.
                !          !Thaw - accomodating corner sources
                !          if(grid%octreeRoot%twoD) then
                !             distToEdge = abs(source(iSource)%position%x) - abs((grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize))
                !             if ( distToEdge < grid%halfSmallestSubcell) then
                !                source(iSource)%onCorner = .true.
                !             end if
                !          else if(grid%octreeRoot%threeD) then
                !             distToEdge = abs(source(iSource)%position%x) - abs((grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize))
                !             if ( distToEdge < grid%halfSmallestSubcell) then
                !                distToEdge = abs(source(iSource)%position%y) - abs((grid%octreeRoot%centre%y - grid%octreeRoot%subcellSize))
                !                if ( distToEdge < grid%halfSmallestSubcell) then
                !                   source(iSource)%onCorner = .true.
                !                end if
                !             else
                !                source(iSource)%onCorner = .false.
                !             end if
                !          end if

                source(iSource)%onCorner = .false.


                !This currently only captures corners, edges will come later
                distToEdge = abs(abs((grid%octreeRoot%centre%z - grid%octreeRoot%subcellSize))-abs(source(iSource)%position%z))
                !       print *, "distToEdge ", distToEdge
                !       print *, "smallestSubcell", grid%halfSmallestSubcell
                !       print *, "abs(source(iSource)%position%z)", abs(source(iSource)%position%z)       
                !       print *, "grid%octreeRoot%centre%z", grid%octreeRoot%centre%z




                ! Find distance of source from upper or lower z boundary, whichever is smaller.


                distToedge = min ( abs( source(iSource)%position%z - &
                     (grid%octreeRoot%centre%z - grid%octreeRoot%subcellSize) ), &
                     abs( source(iSource)%position%z - (grid%octreeRoot%centre%z + grid%octreeRoot%subcellSize) ) )


                if (.not.inOctal(grid%octreeRoot, source(iSource)%position)) then
                   source(isource)%outsideGrid = .true.
                   source(isource)%distance = modulus(source(isource)%position)*1.d10


                   !source(isource)%luminosity = source(isource)%luminosity * (2.d0*grid%octreeRoot%subcellSize*1.d10)**2 / &
                   !     (fourPi*source(isource)%distance**2)
                else if ( distToEdge < grid%halfSmallestSubcell) then          
                   source(iSource)%onEdge = .true.
                   source(iSource)%onCorner = .false.
                   !Thaw - accomodating corner sources

                   if(grid%octreeRoot%twoD) then
                      distToEdge = abs(abs((grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize)) &
                           -abs(source(iSource)%position%x))
                      !             distToEdge = abs(abs(source(iSource)%position%x) - abs((grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize)))

                      if ( grid%octreeRoot%cylindrical ) then 
                         source(iSource)%onCorner = .false.
                      else if(grid%octreeRoot%twoD) then
                         distToEdge = abs(source(iSource)%position%x) &
                              - abs((grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize))
                      end if
                      if ( distToEdge < grid%halfSmallestSubcell) then
                         source(iSource)%onCorner = .true.
                      end if
                   else if(grid%octreeRoot%threeD) then
                      distToEdge = abs(abs((grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize)) &
                           -abs(source(iSource)%position%x))
                      !             distToEdge = abs(abs(source(iSource)%position%x) - abs((grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize)))
                      if ( distToEdge < grid%halfSmallestSubcell) then
                         distToEdge = abs(abs((grid%octreeRoot%centre%y - grid%octreeRoot%subcellSize)) &
                              -abs(source(iSource)%position%y))
                         !                distToEdge = abs(abs(source(iSource)%position%y) - abs((grid%octreeRoot%centre%y - grid%octreeRoot%subcellSize)))
                         if ( distToEdge < grid%halfSmallestSubcell) then
                            source(iSource)%onCorner = .true.
                         end if
                      else
                         source(iSource)%onCorner = .false.
                      end if
                   end if

                   if(source(iSource)%onCorner .and. myRankGlobal == 1) then
                      write(*,*) "Source on corner!"
                   else  if(source(iSource)%outsideGrid) then
                      write(*,*) "Source outside grid!"
                   end if

                   if(source(iSource)%onCorner) call writeinfo("Source is on corner", TRIVIAL)
                   if(source(iSource)%onedge)   call writeinfo("Source is on edge",   TRIVIAL)

                endif

!                if (cylindricalHydro) then
!                   source(iSource)%onCorner = .true.
!                endif

                if(source(iSource)%outsideGrid .and. myRankGlobal == 1) then
                   write(*,*) "Source outside grid!"
                end if

                if (cylindricalHydro) then
                   source(isource)%onEdge = .false.
                   source(isource)%onCorner = .false.
                endif

                select case(inputContFluxFile(isource))
                case("blackbody")
                   !          if(biasToLyman) then
                   !!             call fillSpectrumBB(source(isource)%spectrum, source(isource)%teff, 10.d0, 1000.d4, & 
                   !                  1000)
                   !          else
                   call fillSpectrumBB(source(isource)%spectrum, source(isource)%teff, 10.d0, 1000.d4,1000)
                   !          end if
                case("kurucz")
                   call fillSpectrumKurucz(source(isource)%spectrum, source(isource)%teff, source(isource)%mass, &
                        source(isource)%radius*1.d10)
                case("tlusty")
                   call fillSpectrumTlusty(source(isource)%spectrum, source(isource)%teff, source(isource)%mass, &
                        source(isource)%radius*1.d10)
                case DEFAULT
                   call readSpectrum(source(isource)%spectrum, inputcontfluxfile(isource), ok)
                end select

                call normalizedSpectrum(source(isource)%spectrum)
                !       lamStart = 10.d0
                !       lamEnd = 1000.d4
                !       nlambda = 1000
                call buildSphere(source(isource)%position, dble(source(isource)%radius), &
                     source(isource)%surface, 100, source(isource)%teff, & 
                     source(isource)%spectrum)
                call sumSurface(source(isource)%surface, sumSurfaceluminosity)
                fac = fourPi * stefanBoltz * (source(1)%radius*1.d10)**2 * (source(1)%teff)**4
                !       if (abs(fac-source(1)%luminosity)/source(1)%luminosity > 0.01d0) then
                !          if (myrankGlobal==0) then
                !             write(*,*) "WARNING: luminosity from effective temperature and that from the SED differ by >1%"
                !             write(*,*) "Lum from Teff (lSol): ",fac/lSol
                !             write(*,*) "Lum from SED  (lSol): ",source(1)%luminosity/lSol
                !             write(*,*) "Implied source of accretion luminosity of: ",(source(1)%luminosity-fac)/lSol
                !             tmp = 1.d0/(1.d0/(source(1)%radius * 1.e10) - 1.d0/(rInner*1.d10)) ! [cm]
                !             fac = (source(1)%luminosity-fac) * tmp / (bigG * mCore)
                !             write(*,*) "Mass accretion rate could be ",fac/mSol/ secstoYears, " solar masses/year"
                !          endif
                !       endif

                fac = fourPi * stefanBoltz * (source(isource)%radius*1.d10)**2 * (source(isource)%teff)**4
                write(message,*) "Lum from spectrum / lum from teff ",sumSurfaceLuminosity/fac
                call writeInfo(message, TRIVIAL)
                write(message,*) "Setting source luminosity to luminosity from spectrum: ",sumSurfaceLuminosity/lsol, " lsol"
                call writeInfo(message, TRIVIAL)
                source(iSource)%luminosity = sumSurfaceLuminosity
             endif
             else

                source(isource)%stellar = .false.
                source(isource)%diffuse = .true.
                select case (diffuseType(iSource))

                case("isrf")
                   call createSourceISRF(source(isource), grid)
                case("cmb")
                   call createSourceCMB(source(isource), grid)
                case DEFAULT
                   call writeFatal("Diffuse source type not recognized")
                end select
          endif
       end do
       if (.not.donBodyOnly) then
          if (Writeoutput) call testRandomSource(source, nsource)
          if (Writeoutput) call testSourceSpectrum(source, nsource)
       endif

    endif
!    call writeVtkFile(nsource, source, "sources_at_setup.vtk")

  end subroutine setupSources

  subroutine doPhysics(grid)
    use phasematrix_mod
    use dust_mod
    use inputs_mod, only : atomicPhysics, photoionPhysics, photoionEquilibrium, cmf, nBodyPhysics
    use inputs_mod, only : dustPhysics, lowmemory, radiativeEquilibrium, gasOpacityPhysics
    use inputs_mod, only : statisticalEquilibrium, nAtom, nDustType, nLucy, &
         lucy_undersampled, molecularPhysics, hydrodynamics!, UV_vector
    use inputs_mod, only : useDust, realDust, variableDustSublimation, massEnvelope
    use inputs_mod, only : mCore, solveVerticalHydro, sigma0!, scatteredLightWavelength,  storeScattered
    use inputs_mod, only : tEnd, tDump
    use gas_opacity_mod
#ifdef CMFATOM
    use modelatom_mod, only : globalAtomArray
    use cmf_mod, only : atomloop
#endif
#ifdef STATEQ
    use stateq_mod, only: amrstateqNew
#endif
#ifdef HYDRO
    use nbody_mod, only : donBodyonly
#endif
    use source_mod, only : globalNsource, globalSourceArray, randomSource
    use lucy_mod, only : lucyRadiativeEquilibriumAMR, setFixedtemperatureOnTau
    use setupamr_mod, only: doSmoothOnTau
    use disc_hydro_mod, only: verticalHydrostatic

#ifdef PHOTOION
    use photoion_mod, only : photoionizationLoop
#endif

#ifdef HYDRO
    use inputs_mod, only : hydrodynamics
#ifdef MPI
    use hydrodynamics_mod, only : doHydrodynamics, setupevenuparray
#endif
#endif

#ifdef PDR
!!    use nrayshealpix, only : donrayshealpix
    use pdr_mod, only : PDR_MAIN!, castAllRaysOverGrid
    use inputs_mod, only pdrcalc
#endif

#ifdef MPI
#ifdef PHOTOION
    use photoionAMR_mod, only: photoionizationLoopAMR, ionizegrid

    use photoion_utils_mod, only: setupphotogrid
!    use inputs_mod, only : optimizeStack

#ifdef HYDRO
    use photoionAMR_mod, only: radiationHydro
#endif
#endif
#endif

#ifdef MOLECULAR
    use molecular_mod, only : molecularLoop, globalMolecule
    use inputs_mod, only : lowmemory, molecularPhysics,  useDust, realDust
#ifdef MPI
    use inputs_mod, only : hydrovelocityconv
!    use mpi_amr_mod, only : fillVelocityCornersFromHydro
    use amr_mod, only : hydroVelocityConvert
#endif
#endif

    real, pointer :: xArray(:) => null()
    real(double), pointer :: xArrayDouble(:) => null()
    integer :: nLambda 
    real(double) :: packetWeight
    real(double) :: temp, dustMass, ksca, kabs
    real :: temp2
    integer :: i
    type(PHASEMATRIX), pointer :: miePhase(:,:,:) => null()
    integer, parameter :: nMuMie = 200
    integer :: nlower, nupper, sign
    type(GRIDTYPE) :: grid
#ifdef MPI
#ifdef PHOTOION
    integer :: evenuparray(nthreadsGlobal-1)
    real :: iterTime(3)
    integer :: iterStack(3)
    integer :: optID
#endif
#endif
    nLower = 2
    nUpper = 3

!#ifdef MPI
!#ifdef PHOTOION
!    if(optimizeStack .and. photoionPhysics .and. photoionEquilibrium) then
!       call writeInfo("Optimizing photon stack size.", TRIVIAL)
!       call setupevenuparray(grid, evenuparray)
!       iterTime = 1.e20
!       call photoIonizationloopAMR(grid, globalsourceArray, globalnSource, nLambda, xArray, 200, 1.d40, 1.d40, .false., &
!            iterTime, .true., evenuparray, sublimate=.false.)
!    end if
!#endif
!#endif

    if (gasOpacityPhysics) then
       call setupXarray(grid, xarray, nLambda, dustRadeq=.true.)
       allocate(xArrayDouble(1:nLambda))
       xArrayDouble = dble(xArray)
       call createAllMolecularTables(nLambda, xArrayDouble)
       deallocate(xArrayDouble)
       if (writeoutput) then
          do i = 1, nLambda
             call returnGasKappaValue(grid, 2000., 1.d0, kappaAbs=kabs, kappaSca=ksca, lambda=xarray(i), ilambda=i)
             write(*,'(1p,5e13.4)') xArray(i),kabs+ksca,kabs,ksca, ksca/(ksca+kabs)
          enddo
       endif
    endif




     if (dustPhysics.and.radiativeEquilibrium) then
        call setupXarray(grid, xarray, nLambda, dustRadeq=.true.)
        if (globalnSource > 0) then
           call randomSource(globalsourcearray, globalnSource, &
                i, packetWeight, grid%lamArray, grid%nLambda, initialize=.true.)  
        endif
        call setupDust(grid, xArray, nLambda, miePhase, nMumie, filestart="dust")

        if (grid%geometry == "shakara") then
           call fillDustShakara(grid, grid%octreeRoot, dustmass)
        endif
!        call fillDustUniform(grid, grid%octreeRoot)
#ifdef MPI
        call randomNumberGenerator(syncIseed=.true.)
#endif

!        call fillDustUniform(grid, grid%octreeRoot)
        if (.not.variableDustSublimation) call doSmoothOnTau(grid)

        
!        scatteredlightWavelength = 2.2d4 ! 2.2 microns
!        storeScattered = .true.
#ifdef MPI
        call randomNumberGenerator(randomSeed=.true.)
#endif

     call returnKappa(grid, grid%octreeRoot, 1, atthistemperature=1500., rosselandKappa = temp)
     call returnKappa(grid, grid%octreeRoot, 1, atthistemperature=10000., kappap = temp2)
!     write(*,*) "Ross (1500), Planck (10000): ",temp,temp2/grid%octreeRoot%rho(1)

!     iLambda = findIlambda(1.e5, grid%lamArray, nLambda, ok)
!     call setFixedTemperatureOnTau(grid, iLambda)

        if (solveVerticalHydro) then

           call verticalHydrostatic(grid, mCore, sigma0, miePhase, nDustType, nMuMie, nLambda, xArray, &
                globalsourcearray, globalnSource, nLucy, massEnvelope)
        else
           call lucyRadiativeEquilibriumAMR(grid, miePhase, nDustType, nMuMie, nLambda, xArray, &
                globalsourcearray, globalnSource, nLucy, massEnvelope, lucy_undersampled, finalPass=.true.)
        endif
     endif

#ifdef MOLECULAR
     if (molecularPhysics.and.statisticalEquilibrium) then
#ifdef MPI
        if (grid%splitOverMPI.and.hydrovelocityconv) then
           call setallUnchanged(grid%octreeRoot)
           call hydroVelocityConvert(grid%octreeRoot)
!           if (myrankGlobal /= 0) then
!              call fillVelocityCornersFromHydro(grid)
!           endif
        endif
#endif
        if (dustPhysics) then
           call setupXarray(grid, xarray, nLambda,dustRadeq=.true.)
           call setupDust(grid, xArray, nLambda, miePhase, nMumie)
           usedust = .true.
           realdust = .true.
        endif
        lowMemory = .false.
        call molecularLoop(grid, globalMolecule)
        call writeVTKfile(grid, "molResults.vtk", valueTypeString=(/"J=0       ",&
             "J=1       ", "J=2       ", "J=3       ", "J=4       ", "J=5       "/))
     endif
#endif

#ifdef CMFATOM
     if (atomicPhysics.and.statisticalEquilibrium.and.cmf) then
        call atomLoop(grid, nAtom, globalAtomArray, globalnsource, globalsourcearray)
     endif
#endif
#ifdef STATEQ
     if (atomicPhysics.and.statisticalEquilibrium.and.(.not.cmf)) then
        call amrStateqnew(grid, .false., nLower, nUpper, globalSourceArray(1)%surface,&
                       recalcPrevious=.false.)
        call writeVTKfile(grid, "eta.vtk", valueTypeString=(/"etaline   ",&
             "sourceline"/))
     endif
#endif

#ifdef PHOTOION
     if (photoionPhysics.and.photoionEquilibrium.and. .not. hydrodynamics) then
        sign = 1
        call setupXarray(grid, xArray, nLambda,photoion=.true.)
        if (dustPhysics) call setupDust(grid, xArray, nLambda, miePhase, nMumie)
        if (dustPhysics) call fillDustUniform(grid, grid%octreeRoot)
        if (dustPhysics) call setupOrigDustFraction(grid%octreeRoot)

        if (.not.grid%splitOverMPI) then
           call photoIonizationloop(grid, globalsourceArray, globalnSource, nLambda, xArray )
        else
#ifdef MPI
           call setupevenuparray(grid, evenuparray)
           
!           if(.not. startFromNeutral) then
           print *, "ionizing grid"
           call ionizeGrid(grid%octreeRoot)
           call setupPhotoGrid(grid%octreeRoot)
 !          endif
           call photoIonizationloopAMR(grid, globalsourceArray, globalnSource, nLambda, xArray, 20, 1.d40, &
                1.d40, .false.,iterTime,.true., evenuparray, optID, iterStack, miePhase, nMuMie, sublimate=.false.)

#else
           call writeFatal("Domain decomposed grid requires MPI")
#endif
        endif
     end if
#endif

#ifdef PDR
     if(pdrcalc) then
        call PDR_MAIN(grid, globalsourceArray, globalnSource)     
     end if
#endif

#ifdef HYDRO
     if (hydrodynamics) then
        if (.not.photoionPhysics) then
#ifdef MPI 
        call dohydrodynamics(grid)
#else
        call writeFatal("hydrodynamics not available in single processor version")
        stop
#endif
     else 
        call setupXarray(grid, xArray, nLambda,photoion=.true.)
        if (dustPhysics) call setupDust(grid, xArray, nLambda, miePhase, nMumie)
#ifdef PHOTOION
#ifdef MPI 
        call radiationHydro(grid, globalSourceArray, globalNSource, nLambda, xArray, miePhase, nMuMie)
#else
        call writeFatal("hydrodynamics not available in single processor version")
        stop
#endif
#endif
        call torus_mpi_barrier
     endif
  endif
#endif

#ifdef HYDRO
  if (nbodyPhysics.and.(.not.hydrodynamics)) then
     call  donBodyOnly(tEnd, tdump, grid)
  endif
#endif

! Free memory allocated in this subroutine
  if (associated(xArray)) then
     deallocate(xArray)
     nullify(xArray)
  end if
  if (associated(miePhase)) then
     deallocate(miePhase)
     nullify(miePhase)
  end if


   end subroutine doPhysics

   subroutine setupXarray(grid, xArray, nLambda, lamMin, lamMax, wavLin, numLam, dustRadeq, photoion, atomicDataCube)
     use inputs_mod, only : lamFile, lamFilename, lamLine, vMinSpec, vMaxSpec, nv, resolveSilicateFeature, writepolar, &
          polarWavelength
#ifdef PHOTOION
     use photoion_utils_mod, only : refineLambdaArray
#endif

     logical, optional :: dustRadeq, photoion, atomicDataCube
     type(GRIDTYPE) :: grid
     real, pointer :: xArray(:)
     integer :: nLambda
     real, optional, intent(in) :: lamMin, lamMax
     logical, optional, intent(in) :: wavLin
     integer, optional, intent(in) :: numLam
     integer :: nCurrent, nt, i
     real(double) :: fac, logLamStart, logLamEnd, lamStart,lamEnd, junk

     if (associated(xarray)) deallocate(xarray)

     if (writePolar) then
        nLambda = numLam
        allocate(xarray(1:nLambda))
        xArray(1) = real(polarWavelength)
     endif
     
     if (PRESENT(dustRadEq)) then
        if (dustRadEq) then
              
           if ( present(numLam) ) then 
              nLambda = numLam
           else
              nLambda = 200
           end if

           if (.not.resolveSilicateFeature) then
              allocate(xarray(1:nLambda))
           else
              allocate(xarray(1:nLambda+40))
           endif
           
           if ( present(lamMin) ) then 
              lamStart = lamMin
           else
              lamStart = 100.d0
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
           
           if (lamFile) then
              call setupLamFile
           endif
           
           if (resolveSilicateFeature) then
              nt = 20
              do i = 1, nt
                 xArray(i+nLambda) = real((7.d0 + 7.d0*(dble(i-1)/dble(nt-1)))*1.d4)
              enddo
              nLambda = nLambda + nt
              do i = 1, nt
                 xArray(i+nlambda) = real((15.d0 + 20.d0*(dble(i-1)/dble(nt-1)))*1.d4)
              enddo
              nLambda = nLambda + nt
              call sort(nLambda, xArray)
           endif

        endif
     endif
     
     if (PRESENT(atomicDataCube)) then
        if (atomicDataCube) then
           nLambda = nv
           allocate(xArray(1:nLambda))
           lamStart = lamLine*(1.d0 + (vMinSpec*1.d5)/cSpeed)
           lamend =  lamLine*(1.d0 + (vMaxSpec*1.d5)/cSpeed)
           do i = 1, nLambda
              xArray(i) = real(lamStart+(lamEnd-lamStart)*dble(i-1)/dble(nLambda-1))
           enddo
        endif
     endif

    if (PRESENT(photoion)) then
       if (photoion) then
          nLambda = 1000
          allocate(xarray(1:nLambda))
          lamStart = 10.d0
          lamEnd = 1.d7
          logLamStart = log10(lamStart)
          logLamEnd = log10(lamEnd)
          xArray(1) = real(lamStart)
          xArray(2) = real(lamEnd)
          nCurrent = 2
#ifdef PHOTOION
          call refineLambdaArray(xArray, nCurrent, grid)
#endif
          nt = nLambda - nCurrent
          do i = 1, nt
             fac = logLamStart + real(i)/real(nt+1)*(logLamEnd - logLamStart)
             fac = 10.**fac
             nCurrent=nCurrent + 1
             xArray(nCurrent) = real(fac)
             call sort(nCurrent, xArray)
          enddo
       endif
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
            xArray(i) =real(fac)
        enddo
        
       end subroutine setupLogSpacing

       subroutine setupLamfile

         if (associated(xArray)) then
            deallocate(xArray)
         endif


         call writeInfo("Reading wavelength points from file.", TRIVIAL)
         open(77, file=lamfilename, status="old", form="formatted")
         nLambda = 1
         ! Count the number of entries
333      continue
         read(77,*,end=334) junk
         nLambda = nLambda + 1
         goto 333
334      continue
         nlambda = nlambda - 1
         allocate(xArray(nlambda))
         ! Rewind the file and read them in
         rewind(77)
         do i = 1, nLambda
            read(77,*) xArray(i)
         enddo
         close(77)

       end subroutine setupLamfile

       subroutine setupLinSpacing

         if (nLambda > 1) then
            do i = 1, nlambda
               xArray(i) = real(LamStart + real(i-1)/real(nLambda-1)*(LamEnd - LamStart))
            enddo
         else 
            xArray(1) = real(lamStart)
         endif
       end subroutine setupLinSpacing

   end subroutine setupXarray

   subroutine setupGlobalSources(grid)     
     use parallel_mod
#ifdef SPH
     use sph_data_class
#endif
     use starburst_mod
#ifdef MPI
#ifdef HYDRO
     use hydrodynamics_mod, only : gatherSinks
     use inputs_mod, only: splitOverMPI
#endif
#endif
     use source_mod, only : globalNsource, globalSourceArray
     use inputs_mod, only : inputNsource, mstarburst, lxoverlbol, readsources, &
          hosokawaTracks, nbodyPhysics, nSphereSurface, discardSinks, hotSpot, starburst
#ifdef MPI
     use mpi
#endif

!     integer, parameter :: maxSources =1000
     integer, parameter :: maxSources =10000
     integer(bigInt) :: itest
     integer :: i
     type(GRIDTYPE) :: grid
     real(double) :: coreContinuumFlux, lAccretion, xRayFlux
     real :: fAccretion
     character(len=80) :: message

     if (associated(globalsourceArray)) then
        deallocate(globalSourceArray)
     endif

 !    print *, "setting up global sources"

     if (readSources.or.(inputNsource > 0 )) call writeBanner("Source setup","-",TRIVIAL)
     globalNSource = 0
     allocate(globalsourceArray(1:maxSources))
     if (readSources.or.(inputNsource > 0)) then
        globalnSource = inputNSource
        call setupSources(globalnSource, globalsourceArray, grid)
     endif

!     print *, "delta"
     if (grid%geometry == "theGalaxy" .or. grid%geometry == "sphfile") then
!        if(.not. discardsinks .and. 0 == 1) then
        if(.not. discardsinks) then
           !        print *, "echo"
#ifdef SPH
           globalnSource = get_nptmass()
           !        print *, "using ", globalnsource, "sources"
           if ( globalnSource > size(globalSourceArray)) then 
              write(message,*) "Number of sources exceeds size of source array", globalnSource, size(globalSourceArray)
              call writeFatal(message)
           endif
           do i = 1, globalnSource
              globalSourceArray(i)%stellar = .true.
              globalSourceArray(i)%mass = get_pt_mass(i) * get_umass()
              globalSourceArray(i)%position = get_pt_position(i) * (get_udist()/1.d10)
              globalSourceArray(i)%velocity = get_pt_velocity(i) * get_udist() / get_utime()
           enddo
#else
           call writeFatal("This geometry requires SPH functionality which is not built.")
#endif
#ifdef MPI
#ifdef HYDRO
           if (splitOverMPI) call gatherSinks()
#endif 
#endif
           !        print *, "setting up source array properties"
           if (nbodyPhysics.and.hosokawaTracks)  then
              call setSourceArrayProperties(globalsourceArray, globalnSource, 1.d0)
              print *, "writing source list"
              !        call writeSourceList(globalsourceArray, globalnSource)
#ifdef MPI
              !        call MPI_BARRIER(MPI_COMM_WORLD, ier)
              if(myrankglobal == 0) then
              !           print *, " " 
                 call writeIonizingFLuxes(globalsourceArray, globalnSource)
           
                 !           call MPI_BARRIER(MPI_COMM_WORLD, ier)
              endif
#endif
           endif
        endif
     endif

     if (starburst.and.(.not.readSources)) then
#ifdef MPI
        call randomNumberGenerator(randomSeed = .true.)
        call randomNumberGenerator(syncIseed=.true.)
#endif
        call randomNumberGenerator(getIseed=itest)
        allocate(globalsourcearray(1:10000))
        globalsourceArray(:)%outsideGrid = .false.
        globalnSource = 0
        call createSources(globalnSource,globalsourcearray, "instantaneous", 1.d6, mStarburst, 1.d0)
        call randomNumberGenerator(randomSeed = .true.)
    endif
    

    if (grid%geometry(1:6) == "ttauri") then
       coreContinuumFlux = 0.d0
       call addXray(globalSourceArray(1)%spectrum, lxoverlbol)
       xRayFlux = integrateSpectrumOverBand(globalSourceArray(1)%spectrum, &
            (cSpeed / (10.d0 * 1000.d0 * evtoerg/hcgs))*1.d8, &
            (cSpeed / (0.1d0 * 1000.d0 * evtoerg/hcgs))*1.d8)
       if (writeoutput) write(*,*) "EUV/X-ray flux (erg/s)", &
            xRayFlux*fourPi*globalSourceArray(1)%radius**2*1.d20

       call buildSphere(globalsourceArray(1)%position, globalSourceArray(1)%radius, &
            globalsourcearray(1)%surface, nSphereSurface, &
            globalsourcearray(1)%teff, globalsourceArray(1)%spectrum)
       call genericAccretionSurface(globalsourcearray(1)%surface, 1.e16, coreContinuumFlux,fAccretion, lAccretion) 
       call writeVTKfileSource(1, globalSourceArray(1:1), "source.vtk")
       if (writeoutput) write(*,*) "Added accretion luminosity of ",lAccretion/lsol, " lsol"
       if (writeoutput) write(*,*) "Accretion luminosity in stellar units ",lAccretion/globalsourceArray(1)%luminosity
       globalsourcearray(1)%luminosity = globalsourcearray(1)%luminosity + lAccretion
       globalNSource = 1
    endif

    if (hotSpot) then
       coreContinuumFlux = 0.d0

       call buildSphere(globalsourceArray(1)%position, globalSourceArray(1)%radius, &
            globalsourcearray(1)%surface, nSphereSurface, &
            globalsourcearray(1)%teff, globalsourceArray(1)%spectrum)
       call hotSpotSurface(globalsourcearray(1)%surface, real(cspeed/(6562.8d-8)), coreContinuumFlux,fAccretion, lAccretion) 
       call writeVTKfileSource(1, globalSourceArray(1:1), "source.vtk")
       if (writeoutput) write(*,*) "Added accretion luminosity of ",lAccretion/lsol, " lsol"
       if (writeoutput) write(*,*) "Accretion luminosity in stellar units ",lAccretion/globalsourceArray(1)%luminosity
       globalsourcearray(1)%luminosity = globalsourcearray(1)%luminosity + lAccretion
       globalNSource = 1
    endif

!    if ((myrankGlobal==0).and.(globalnSource > 0)) call writeVtkfileSource(globalnSource, globalsourcearray, "source.vtk")


end subroutine setupGlobalSources

subroutine setupDust(grid, xArray, nLambda, miePhase, nMumie, fileStart)
  use inputs_mod, only : mie, nDustType, readDustFromFile, writeDustToFile
  use phasematrix_mod
  use dust_mod
  use gas_opacity_mod
  type(GRIDTYPE) :: grid
  real, pointer :: xArray(:)
  character(len=*), optional :: fileStart
  integer :: nLambda, i
  type(PHASEMATRIX), pointer :: miePhase(:,:,:)
  integer :: nMuMie
  character(len=80) :: dustFilename

  grid%oneKappa = .true.
  mie = .true.
  grid%nTempRossArray = 1000
  if (associated(grid%kappaRossArray)) deallocate(grid%kappaRossArray)
  if (associated(grid%tempRossArray)) deallocate(grid%tempRossArray)
  allocate(grid%kappaRossArray(nDustType,1:grid%nTempRossArray))
  allocate(grid%tempRossArray(1:grid%nTempRossArray))


  if (.not.readDustFromFile) then
     call  createDustCrossSectionPhaseMatrix(grid, xArray, nLambda, miePhase, nMuMie)
     if (writeDustToFile) then
        do i = 1, nDustType
           write(dustFilename, '(a,i2.2,a)') "dust_",i,".dat"
           if (PRESENT(fileStart)) write(dustFilename, '(a,i2.2,a)') trim(fileStart),i,".dat"
           call writeDust(dustfilename, i, grid, xArray, nLambda, miePhase, nMuMie)
        enddo
     endif
  else
     if (associated(grid%onekappaAbs)) deallocate(grid%onekappaAbs)
     if (associated(grid%onekappaSca)) deallocate(grid%onekappaSca)
     if (associated(miePhase)) deallocate(miePhase)
     allocate(grid%oneKappaAbs(1:nDustType, 1:nLambda))
     allocate(grid%oneKappaSca(1:nDustType, 1:nLambda))
     allocate(miePhase(1:nDustType,1:nLambda,1:nMuMie))
     do i = 1, nDustType
        write(dustFilename, '(a,i2.2,a)') "dust_",i,".dat"
        if (PRESENT(fileStart)) write(dustFilename, '(a,i2.2,a)') trim(fileStart),i,".dat"
        call readDust(dustfilename, i, grid, xArray, nLambda, miePhase, nMuMie)
     enddo

     call writeInfo("Creating Rosseland opacity lu table",TRIVIAL)
     call createRossArray(grid)
     call writeInfo("Done.",TRIVIAL)
     call returnKappa(grid, grid%OctreeRoot, 1, reset_kappa=.true.)
  end if

  if (associated(grid%octreeRoot)) call allocateMemoryForDust(grid%octreeRoot)



end subroutine setupDust

subroutine testSuiteRandom()
  integer :: i
  integer(bigint) :: iseed
  real :: r
  real(double) :: rd
#ifdef _OPENMP
  integer :: np, omp_get_num_threads, omp_get_thread_num, j
  real :: array(100,10)
#endif
  character(len=80) :: message

  call randomNumberGenerator(randomSeed=.true.)
  call randomNumberGenerator(getIseed=iseed)
  write(message,*) "Random number seed is initialized from clock with ",iseed
  call writeInfo(message)
  if (writeoutput) write(*,*) "Sequence of 10 random numbers (real , double)"
  
  do i = 1, 10
     call randomNumberGenerator(getReal=r)
     call randomNumberGenerator(getIseed=iseed)
     call randomNumberGenerator(getDouble=rd)
     if (writeoutput) write(*,*) i, r, rd
  enddo

#ifdef _OPENMP
  if (writeoutput) write(*,*) "OMP thread test - random seed"
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(I,j) &
  !$OMP SHARED(ARRAY,np)
  call randomNumberGenerator(randomSeed=.true.)
  np = omp_get_num_threads()
  do i = 1, 10
     j = omp_get_thread_num()+1
     call randomNumberGenerator(getReal=array(i,j))
  enddo
  !$OMP END PARALLEL
  do i = 1, 10
     write(*,*) array(i,1:np)
  enddo
#endif


#ifdef _OPENMP
  if (writeoutput) write(*,*) "OMP thread test - fixed seed"
  call randomNumberGenerator(randomSeed=.true.)
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(I,j) &
  !$OMP SHARED(ARRAY,np)
  call randomNumberGenerator(syncISeed = .true.)
  np = omp_get_num_threads()
  do i = 1, 10
     j = omp_get_thread_num()+1
     call randomNumberGenerator(getReal=array(i,j))
  enddo
  !$OMP END PARALLEL
  do i = 1, 10
     write(*,*) array(i,1:np)
  enddo
#endif

end subroutine testSuiteRandom



end module physics_mod
