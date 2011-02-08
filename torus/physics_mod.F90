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
     call addIons(globalIonArray, nGlobalIon)
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
    real(double) :: distToEdge, fac, sumSurfaceLuminosity
    character(len=100) :: message

    do iSource = 1, nSource
   
       source(iSource)%teff = sourceTeff(iSource)
       source(iSource)%mass = sourceMass(iSource)
       source(iSource)%radius = sourceRadius(iSource)
       source(iSource)%position = sourcePos(iSource)
       source(isource)%luminosity = fourPi * stefanBoltz * &
            (source(isource)%radius*1.d10)**2 * (source(isource)%teff)**4
       source(iSource)%prob = sourceProb(iSource)
       source(iSource)%pointSource = pointSourceArray(iSource)

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


       distToedge = min ( abs( source(iSource)%position%z - (grid%octreeRoot%centre%z - grid%octreeRoot%subcellSize) ), &
                          abs( source(iSource)%position%z - (grid%octreeRoot%centre%z + grid%octreeRoot%subcellSize) ) )


       if (.not.inOctal(grid%octreeRoot, source(iSource)%position)) then
          source(isource)%outsideGrid = .true.
          source(isource)%distance = modulus(source(isource)%position)*1.d10
          source(isource)%luminosity = source(isource)%luminosity * (2.d0*grid%octreeRoot%subcellSize*1.d10)**2 / &
               (fourPi*source(isource)%distance**2)
       else if ( distToEdge < grid%halfSmallestSubcell) then          
          source(iSource)%onEdge = .true.
          source(iSource)%onCorner = .false.
          !Thaw - accomodating corner sources

          if(grid%octreeRoot%twoD) then
             distToEdge = abs(abs((grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize))-abs(source(iSource)%position%x))
!             distToEdge = abs(abs(source(iSource)%position%x) - abs((grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize)))

          if ( grid%octreeRoot%cylindrical ) then 
             source(iSource)%onCorner = .false.
          else if(grid%octreeRoot%twoD) then
             distToEdge = abs(source(iSource)%position%x) - abs((grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize))
          end if
             if ( distToEdge < grid%halfSmallestSubcell) then
                source(iSource)%onCorner = .true.
             end if
          else if(grid%octreeRoot%threeD) then
             distToEdge = abs(abs((grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize))-abs(source(iSource)%position%x))
!             distToEdge = abs(abs(source(iSource)%position%x) - abs((grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize)))
             if ( distToEdge < grid%halfSmallestSubcell) then
                distToEdge = abs(abs((grid%octreeRoot%centre%y - grid%octreeRoot%subcellSize))-abs(source(iSource)%position%y))
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

       if(source(iSource)%outsideGrid .and. myRankGlobal == 1) then
         write(*,*) "Source outside grid!"
       end if

       select case(inputContFluxFile(isource))
          case("blackbody")
             call fillSpectrumBB(source(isource)%spectrum, source(isource)%teff, 10.d0, 1000.d4,1000)
          case("kurucz")
             call fillSpectrumKurucz(source(isource)%spectrum, source(isource)%teff, source(isource)%mass, &
                  source(isource)%radius*1.d10)
          case DEFAULT
             call readSpectrum(source(isource)%spectrum, inputcontfluxfile(isource), ok)
       end select

       call normalizedSpectrum(source(isource)%spectrum)
!       lamStart = 10.d0
!       lamEnd = 1000.d4
!       nlambda = 1000
       call buildSphere(source(isource)%position, dble(source(isource)%radius), &
            source(isource)%surface, 400, source(isource)%teff, &
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
    end do
  end subroutine setupSources

  subroutine doPhysics(grid)
    use phasematrix_mod
    use dust_mod
    use modelatom_mod, only : globalAtomArray
    use input_variables, only : atomicPhysics, photoionPhysics, photoionEquilibrium
    use input_variables, only : dustPhysics, lowmemory, radiativeEquilibrium
    use input_variables, only : statisticalEquilibrium, nAtom, nDustType, nLucy, &
         lucy_undersampled, molecularPhysics, hydrodynamics, cmf, lte, nlower, nupper
    use input_variables, only : useDust, realDust, readlucy, writelucy, variableDustSublimation
    use input_variables, only : lucyfilenameOut, lucyFilenamein, massEnvelope
    use input_variables, only : mCore, solveVerticalHydro, sigma0, scatteredLightWavelength,  storeScattered
    !use input_variables, only : radiationHydrodynamics
    use cmf_mod, only : atomloop
    use photoion_mod, only : photoionizationLoop
#ifdef MPI
    use photoionAMR_mod, only: photoionizationLoopAMR, radiationHydro !
    use hydrodynamics_mod, only : doHydrodynamics
#endif
    use source_mod, only : globalNsource, globalSourceArray, randomSource
    use molecular_mod, only : molecularLoop, globalMolecule
    use lucy_mod, only : lucyRadiativeEquilibriumAMR
#ifdef MPI
    use input_variables, only : hydrovelocityconv
!    use mpi_amr_mod, only : fillVelocityCornersFromHydro
    use amr_mod, only : hydroVelocityConvert
#endif
    use setupamr_mod, only: doSmoothOnTau
    use disc_hydro_mod, only: verticalHydrostatic

    real, pointer :: xArray(:) => null()
    integer :: nLambda 
    real(double) :: packetWeight
    integer :: i
    type(PHASEMATRIX), pointer :: miePhase(:,:,:) => null()
    integer, parameter :: nMuMie = 20
    type(GRIDTYPE) :: grid


     if (dustPhysics.and.radiativeEquilibrium) then
        call setupXarray(grid, xarray, nLambda)
        if (globalnSource > 0) then
           call randomSource(globalsourcearray, globalnSource, &
                i, packetWeight, grid%lamArray, grid%nLambda, initialize=.true.)  
        endif
        call setupDust(grid, xArray, nLambda, miePhase, nMumie)
!        call fillDustUniform(grid, grid%octreeRoot)
#ifdef MPI
        call randomNumberGenerator(syncIseed=.true.)
#endif

        call fillDustUniform(grid, grid%octreeRoot)
        if (.not.variableDustSublimation) call doSmoothOnTau(grid)

        
        scatteredlightWavelength = 2.2d4 ! 2.2 microns
        storeScattered = .true.
#ifdef MPI
        call randomNumberGenerator(randomSeed=.true.)
#endif

        if (solveVerticalHydro) then

           call verticalHydrostatic(grid, mCore, sigma0, miePhase, nDustType, nMuMie, nLambda, xArray, &
                globalsourcearray, globalnSource, nLucy, massEnvelope)
        else
           call lucyRadiativeEquilibriumAMR(grid, miePhase, nDustType, nMuMie, nLambda, xArray, &
                globalsourcearray, globalnSource, nLucy, massEnvelope, lucy_undersampled, finalPass=.true.)
        endif
     endif

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
           call setupXarray(grid, xarray, nLambda)
           call setupDust(grid, xArray, nLambda, miePhase, nMumie)
           usedust = .true.
           realdust = .true.
        endif
        lowMemory = .false.
        call molecularLoop(grid, globalMolecule)
     endif
     
     if (atomicPhysics.and.statisticalEquilibrium.and.cmf) then
        call atomLoop(grid, nAtom, globalAtomArray, globalnsource, globalsourcearray)
     endif

     if (atomicPhysics.and.statisticalEquilibrium.and.(.not.cmf)) then
        call amrStateqnew(grid, lte, nLower, nUpper, globalSourceArray(1)%surface,&
                       recalcPrevious=.false.)
     endif


     if (photoionPhysics.and.photoionEquilibrium) then 

        call setupXarray(grid, xArray, nLambda)
        if (dustPhysics) call setupDust(grid, xArray, nLambda, miePhase, nMumie)

        if (.not.grid%splitOverMPI) then
           call photoIonizationloop(grid, globalsourceArray, globalnSource, nLambda, xArray, readlucy, writelucy, &
             lucyfileNameout, lucyfileNamein)
        else
#ifdef MPI
           call photoIonizationloopAMR(grid, globalsourceArray, globalnSource, nLambda, xArray, 100, 1.d40, 1.d40, .false., &
                .true., sublimate=.false.)
#else
           call writeFatal("Domain decomposed grid requires MPI")
#endif
        endif

     end if

     if (hydrodynamics) then
        if (.not.photoionPhysics) then
#ifdef MPI 
        call dohydrodynamics(grid)
#else
        call writeFatal("hydrodynamics not available in single processor version")
        stop
#endif
     else 
        call setupXarray(grid, xArray, nLambda)
        if (dustPhysics) call setupDust(grid, xArray, nLambda, miePhase, nMumie)

#ifdef MPI 
        call radiationHydro(grid, globalSourceArray, globalNSource, nLambda, xArray)
#else
        call writeFatal("hydrodynamics not available in single processor version")
        stop
#endif
        call torus_mpi_barrier
     endif
  endif

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

   subroutine setupXarray(grid, xArray, nLambda, lamMin, lamMax, wavLin, numLam)
     use input_variables, only : photoionPhysics, dustPhysics, molecularPhysics, atomicPhysics
     use input_variables, only : lamFile, lamFilename, lamLine, vMinSpec, vMaxSpec, nv, calcDataCube
     use modelatom_mod, only : globalAtomArray
     use photoion_mod, only : refineLambdaArray
     type(GRIDTYPE) :: grid
     real, pointer :: xArray(:)
     integer :: nLambda
     real, optional, intent(in) :: lamMin, lamMax
     logical, optional, intent(in) :: wavLin
     integer, optional, intent(in) :: numLam
     integer :: nCurrent, nt, i
     real(double) :: fac, logLamStart, logLamEnd, lamStart,lamEnd, junk

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

        if ( present(numLam) ) then 
           nLambda = numLam
        else
           nLambda = 1000
        end if
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

        if (lamFile) then
           call setupLamFile
        endif


     endif

     if (atomicPhysics.and.calcDataCube) then
        nLambda = nv
       allocate(xArray(1:nLambda))
       lamStart = lamLine*(1.d0 + (vMinSpec*1.d5)/cSpeed)
       lamend =  lamLine*(1.d0 + (vMaxSpec*1.d5)/cSpeed)
       do i = 1, nLambda
          xArray(i) = lamStart+(lamEnd-lamStart)*dble(i-1)/dble(nLambda-1)
       enddo
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




     if (dustPhysics.or.photoIonPhysics.or.atomicPhysics.or.molecularPhysics) then
        grid%nLambda = nLambda
        if (associated(grid%lamArray)) deallocate(grid%lamArray)
        allocate(grid%lamArray(1:nLambda))
        grid%lamArray = xarray
     endif

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
               xArray(i) = LamStart + real(i-1)/real(nLambda-1)*(LamEnd - LamStart)
            enddo
         else 
            xArray(1) = lamStart
         endif
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
        call randomNumberGenerator(syncIseed=.true.)
#endif
        allocate(globalsourcearray(1:10000))
        globalsourceArray(:)%outsideGrid = .false.
        globalnSource = 0
        call createSources(globalnSource,globalsourcearray, "instantaneous", 1.d6, 1.d3, 1.d0)
        call randomNumberGenerator(randomSeed = .true.)
    endif
    

    if (grid%geometry(1:6) == "ttauri") then
       coreContinuumFlux = 0.d0
       call buildSphere(globalsourceArray(1)%position, globalSourceArray(1)%radius, &
            globalsourcearray(1)%surface, 1000, &
            globalsourcearray(1)%teff, globalsourceArray(1)%spectrum)
       call genericAccretionSurface(globalsourcearray(1)%surface, grid, 1.e16, coreContinuumFlux,fAccretion, lAccretion) 
       globalsourcearray(1)%luminosity = globalsourcearray(1)%luminosity + lAccretion
       globalNSource = 1
    endif

    if ((myrankGlobal==0).and.(globalnSource > 0)) call writeVtkfileSource(globalnSource, globalsourcearray, "source.vtk")


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
