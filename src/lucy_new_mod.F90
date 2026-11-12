module lucy_new_mod

  use constants_mod
  use messages_mod
  use pah_mod, only: getPAHfreqFromAdot, getKappaAbsPAH, getKappaScaPAH, PAHemissivityFromAdot
  use vector_mod
  use amr_mod, only: addNewChild, inOctal, distanceToCellBoundary, returnKappa, amrGridValues, &
       countvoxels, findsubcelllocal, findsubcelltd
  use utils_mod, only: hunt, locate, solvequaddble
  use gridtype_mod, only: GRIDTYPE
  use octal_mod, only: OCTAL, octalWrapper, subcellCentre, cellVolume
  use spectrum_mod, only: getwavelength
  use atom_mod, only: bnu, blambda
  use timing, only: tune
  use vtk_mod, only: writeVtkFile
  use mpi_global_mod, only: myRankGlobal
  use memory_mod

  implicit none

#ifdef USEMKL
   include 'mkl_vml.fi'
#endif

  ! Cells below this density are treated as empty (e.g. the rho = 1e-30 placeholder inside
  ! an inner radius) and are ignored when checking for undersampled cells.
  real(double), parameter, private :: emptyCellRho = 1.d-29

  ! state for the sublimation ray sweep (sublimationmethod sweep)
  integer, parameter, private :: nSweepBins = 360        ! angular bins for the backwarming factor
  real(double), parameter, private :: sweepTauCrit = 1.d0 ! cells with tau from the star below this are "directly lit"
  integer, parameter, private :: nKPTable = 300           ! Planck-mean absorption table, log T from 1 to 1e5 K
  ! B(theta) per dust type is constrained to be smooth and symmetric about the midplane: each
  ! measure iteration, the measured ratios q = B_measured/B_trial of all lit fronts are fitted
  ! with an even polynomial in theta, and B is multiplied by the (under-relaxed) fit. Per-bin
  ! or per-group B let protrusions of the front lower their own B (they see less of the wall),
  ! which grew fingers and dips.
  real(double), allocatable, save, private :: sweepB(:,:)        ! (nSweepBins,nDustType) trial B
  real(double), allocatable, save, private :: fitRelax(:)        ! (nDustType) relaxation factor of the B update
  real(double), allocatable, save, private :: fitLastSign(:)     ! (nDustType) sign of the last midplane correction
  character(len=240), allocatable, save, private :: sweepFitText(:) ! (nDustType) convergence diagnostic of the last measure sweep
  real(double), allocatable, save, private :: sweepFrontR(:,:)   ! front radius at the previous move sweep
  real(double), allocatable, save, private :: kPAbsTable(:,:)    ! (nKPTable, nDustType), includes 1e10 factor
  real(double), allocatable, save, private :: kPAbsSource(:,:)   ! (nDustType, nSource) Planck mean abs at source Teff
  real(double), allocatable, save, private :: kPExtSource(:)     ! (nDustType) Planck mean ext at source 1 Teff

contains

   function relative_q(source, grid,tdust) result(q)
      use source_mod, only: SOURCETYPE
      use utils_mod, only: loginterp_dble
      use spectrum_mod, only: fillSpectrumBB, SPECTRUMTYPE, normalizedSpectrum
      type(SOURCETYPE) :: source(:)
      type(SPECTRUMTYPE) :: bbSpectrum
      type(GRIDTYPE) :: grid
      real(double) ::  tdust
      real(double) :: q, qstar, qdust, flux, lam, dlam
      integer :: i
      qstar = 0.d0
      qdust = 0.d0
      call fillSpectrumBB(bbSpectrum, tdust, 0.d0, 0.d0, grid%nLambda, grid%lamArray)
      call normalizedSpectrum(bbSpectrum)
      do i = 1, grid%nLambda-1
         dlam = grid%lamArray(i+1) - grid%lamArray(i)
         lam = grid%lamArray(i)
         flux = loginterp_dble(source(1)%spectrum%normflux, source(1)%spectrum%nlambda, source(1)%spectrum%lambda, lam)
         
         qstar = qstar + flux * grid%oneKappaAbs(1,i)* dlam
         qdust = qdust + grid%oneKappaAbs(1,i)*bbSpectrum%normflux(i)*dlam
         q = qstar / qdust
      enddo
      end function



  subroutine lucyRadiativeEquilibriumAMR(grid, miePhase, nDustType, nMuMie, nLambda, lamArray, &
       source, nSource, nLucy, massEnvelope,  percent_undersampled_min, iHydro, finalPass)
    use inputs_mod, only : usemultidust, grainfrac, mdisc, tsub
    use inputs_mod, only : variableDustSublimation, iterlucy, solveVerticalHydro, dustSettling, restartLucy
    use inputs_mod, only : sublimationMethod
    use inputs_mod, only : smoothFactor, lambdasmooth, taudiff, forceLucyConv, multiLucyFiles, doSmoothGridTau
    use inputs_mod, only : object, convergeOnUndersampled, storeScattered, scatteredLightWavelength
    use inputs_mod, only : writelucyTmpfile, discWind, mincrossings, maxiterLucy, solveDiffusionZone, quickSublimate, usePAH
    use source_mod, only: SOURCETYPE, randomSource, getPhotonPositionDirection, distanceToSource
    use phasematrix_mod, only: PHASEMATRIX, newDirectionMie
    use diffusion_mod, only: solvearbitrarydiffusionzones, defineDiffusionOnRosseland, defineDiffusionOnUndersampled, randomwalk, &
         unsetDiffusion
    use amr_mod, only: myScaleSmooth, myTauSmooth, findtotalmass, scaledensityamr
!    use intensity_storage_mod
    use dust_mod, only: filldustuniform, stripdustaway, sublimatedust, sublimatedustwr104, fillDustShakara, &
         normalizeDustFractions, findDustMass, reportMasses, fillDustSettled, &
         getDustMasses, setupDustFracOriginal, adaptiveSublimateDust
    use random_mod
    use gas_opacity_mod, only: atomhydrogenRayXsection
    use gridio_mod, only: writeAMRgrid
#ifdef CHEMISTRY
    use inputs_mod, only : doChemistry
    use chemistry_mod
#endif
#ifdef MPI
    use mpi_global_mod, only: myRankGlobal, nThreadsGlobal
    use mpi
#endif
    implicit none

    type(GRIDTYPE) :: grid
    integer :: nSource
    integer :: iSmoothLam
    logical :: gridConverged
    integer :: nDustType
    logical, optional :: finalPass
    type(SOURCETYPE) :: source(:)
    integer :: nLucy
    ! threshold value for undersampled cell in percent (for stopping iteration).
    real, intent(in) :: percent_undersampled_min  
    type(VECTOR) ::  uHat, uNew, rVec, rHat, olduHat
    type(VECTOR) :: octVec
    type(OCTAL), pointer :: thisOctal, sOctal, tempOctal, foundOctal
    integer :: foundSubcell
    integer :: tempSubcell
    integer, intent(in)  :: nLambda, nMuMie
    type(PHASEMATRIX):: miePhase(:,:,:)
    real  :: lamArray(:)
    real(oct) :: r
    real :: massEnvelope, scaleFac
    real :: tauMax
    real(oct) :: totalMass
    integer :: nFreq
    integer :: i, j
    real(oct) :: freq(nLambda), dnu(nLambda), probDistJnu(nLambda+1)
    !    real(oct) :: probDistPlanck(nFreq)
    real(double) :: kappaScadb, kappaAbsdb, totalLumInPackets, kappaAbsPAH, kappaScaPAH
    integer(double) :: nDiffusion
    integer(double) :: nScat, nAbs
    integer(bigInt) :: nMonte, iMonte
    integer :: thisPhotonAbs, thisPhotonSca
    real(oct) :: thisFreq
    real(oct) :: albedo
    logical :: escaped
    real(oct) :: t1, hrecip_kt
    real(oct) :: thisLam,wavelength
    integer :: iSource
    integer :: iLam
    integer(double) :: nInf
    integer(double) :: nStarHit            ! packets absorbed (and removed) on hitting a stellar surface
    real(double) :: dStar                  ! distance along the packet direction to a stellar surface
    logical :: hitsStarSurface, hitStar
    integer :: iIter, nIter
    real(oct) ::   epsOverDeltaT
    integer :: nDT, nUndersampled
    real(oct) :: totalEmission, oldTotalEmission
    integer :: subcell
    logical :: leftHandBoundary
    real :: treal
    real(oct) :: lCore
    type(vector) :: vec_tmp
    real(oct)::  dT_sum ! [Kelvins]  the sum of the temperature change
    real(oct)::  dT_min ! [kelvins]  the minimum change of temperature
    real(oct)::  dT_max ! [kelvins]  the maximum change of temperature
    real(oct)::  dT_over_T_max ! [kelvins]  the maximum fractional change of temperature
    character(len=80) :: tfilename
    character(len=80) :: message
    logical, save :: first_time_to_open_file = .true.
    real :: totFrac
    integer :: nFrac
    integer, parameter :: LU_OUT = 42
    integer :: iIter_grand, iMultiplier, iHydro
    real(oct)::  dT_mean_new, dT_mean_old ! [Kelvins]
    logical :: converged
    real:: percent_undersampled
    integer :: nKilled
    integer(bigInt)       :: imonte_beg, imonte_end   ! the index boundary for the photon loops.
    real(double)  :: kAbsArray(nlambda),kAbsArray2(nlambda)
    logical :: ok                             ! Did call to randomWalk complete OK?
    logical :: photonInDiffusionZone
    real :: diffusionZoneTemp, temp
    logical :: directPhoton !, smoothconverged
    logical :: thisSmooth
    logical :: thermalPhoton, scatteredPhoton
    integer :: nCellsInDiffusion
    integer :: nColdDiff
    logical :: adaptiveSublimation, sweepSub, sweepMeasure
    integer :: nSubCells, nSubChanged, nSubDamped, nDustSettled
    integer :: nSweepChanged, nSweepExposed, nSweepFronts, nSweepFrontMoved, nSweepFrontConv, nSweepMidUnconv
    real(double) :: sweepFrontErr
    real(double) :: maxSubChange
    character(len=300) :: subMessage
    real(double) :: fac1(nLambda), fac1dnu(nlambda)
    integer :: nVoxels, nOctals, nInFlow
    real(double), allocatable :: dustMassArray(:)
    integer :: icritupper, icritlower
    real(double) :: logt, weight
    real(double) :: logNu1, logNuN, dlogNu, scaleNu
    real(double) :: loglam1, loglamN, scalelam
    real(double) :: q
    real(double) :: logNucritUpper, logNucritLower
    real(double) :: this_bnu(nlambda), fac2(nlambda), hNuOverkT(nlambda)
    real(double) :: packetWeight, wavelengthWeight
    real(double) :: subRadius, dustMass, PAHprob, subradius_monnier
    real :: lamSmoothArray(5)
    logical :: thisIsFinalPass
    real(double) :: farray(10)
    integer(bigInt) :: totMem
   character(len=10) :: stringArray(20)
   character(len=12), allocatable :: lucyValueTypeString(:)
   integer :: nLucyValueTypes
#ifdef USEMKL
    integer :: oldmode
    real(oct) :: hrecip_ktarray(nlambda)
    real(double) :: OneArray(nlambda)
#endif

#ifdef MPI
    ! For MPI implementations =====================================================
    integer ::   ierr           ! error flag
    !  integer ::   n_rmdr, m      !
    integer, parameter :: tag = 0
    integer(bigInt) ::  np, n_rmdr, m  
    real :: buffer_real(nThreadsGlobal)     
    real(double) :: buffer_double(nThreadsGlobal)     

    ! FOR MPI IMPLEMENTATION=======================================================
    writeoutput = .false.
    if (myRankIsZero) then
       writeoutput = .true.

!       print *, ' '
!       print *, 'Lucy radiative equilibrium routine computed by ', nThreadsGlobal, ' processors.'
!       print *, ' '
    endif

    ! ============================================================================
#endif

#ifdef USEMKL
    call writeinfo('Using Intel MKL', TRIVIAL)
    OneArray(1:nlambda) = 1.d0

    oldmode = vmlgetmode()
    oldmode = vmlsetmode(VML_LA)
#endif

    thisIsFinalPass = .false.
    if (PRESENT(finalPass)) thisIsFinalPass = finalPass


    do i = 1, 1000
       lognu1 = log10(1200.d0) + real(i-1)*(log10(1e7)-log10(1200.))/999.
       lognu1 = 10.d0**lognu1
       logt = atomhydrogenRayXsection(lognu1)
!       if (myrankglobal == 1) write(105,*) lognu1, logt/sigmaE
    enddo

    lamSmoothArray = (/5500., 1.e4, 2.e4, 5.e4, 10.e4/)

    oldTotalEmission = 1.d30
    kappaScaDb = 0.d0; kappaAbsDb = 0.d0
    diffusionZoneTemp = 0.d0; foundSubcell = 0; iSource = 0
    kAbsArray = 0.; kAbsArray2 = 0.; leftHandBoundary = .true.; ok = .true.
    photonInDiffusionZone = .false.; rHat = VECTOR(0.d0, 0.d0, 0.d0);temp = 0.
    wavelength = 0.;  nfreq = 0


    nFreq = nLambda
    do i = 1, nFreq
       freq(nFreq-i+1) = cSpeed / (lamArray(i)*1.e-8)
    enddo

    ! Set up factor which is used in bnu calculation
    fac1(:) = (2.d0*hCgs*freq(:)**3)/(cSpeed**2)

    do i = 2, nFreq-1
       dnu(i) = 0.5*((freq(i+1)+freq(i))-(freq(i)+freq(i-1)))
    enddo

    dnu(1) = freq(2)-freq(1)
    dnu(nFreq) = freq(nFreq)-freq(nFreq-1)

    fac1dnu(:) = fac1(:) * dnu(:)

    lognu1 = log(freq(1))
    logNuN = log(freq(nFreq))
    dlogNu = logNuN - lognu1
    scaleNu  = dble(nFreq) / dlogNu 

    loglam1 = log(lamarray(1))
    loglamN = log(lamarray(nFreq))
    scalelam  = dble(nLambda - 1) / dlogNu 



    lCore = 0.d0
    do i = 1, nSource
       lCore = lCore + source(i)%luminosity
    enddo

    write(message,'(a,1pe12.5)') "Total souce luminosity (lsol): ",lCore/lSol
    call writeInfo(message, TRIVIAL)

    ! sublimation radius from equation 1 of Monnier & Millan-Gabet 2002
   q = relative_q(source, grid, tsub(1)) 
   subradius_monnier = 1.1d0 *sqrt(q) *  sqrt(source(1)%luminosity/(1000.d0*lSol))*(1500./tsub(1))**2


    if (.not.restartLucy) then
       if (grid%geometry.eq."wr104") then
          totalMass = 0.
          call findTotalMass(grid%octreeRoot, totalMass)
          scaleFac = real(massEnvelope / totalMass)
          if (writeoutput) write(*,'(a,1pe12.5)') "Density scale factor: ",scaleFac
          call scaleDensityAMR(grid%octreeRoot, dble(scaleFac))
       endif

       call unsetDiffusion(grid%octreeRoot)

     

       !    call setupFreqProb(temperature, freq, dnu, nFreq, ProbDistPlanck)
       
       call writeInfo("Computing lucy radiative equilibrium in AMR...",TRIVIAL)


       if ((.not.discWind).and.(grid%geometry == "shakara")) then
          dustMass = 0.d0
          call fillDustShakara(grid, grid%octreeRoot, dustMass)
          if (dustSettling) call fillDustSettled(grid)
          if (usemultidust) then
             allocate(dustMassArray(1:nDustType))
             dustMassArray = 0.d0 
             call getDustMasses(grid, dustMassArray)
             farray = dble(grainFrac(1:10)*mDisc)
             if (nDustType > 1) call normalizeDustFractions(grid, grid%octreeRoot, dustMassArray(1:10), &
                  farray, 10)
          endif
          call reportMasses(grid)
       endif


       
       if (object == "ab_aur") then
          call writeInfo("Filling dust with large dust in midplane", FORINFO)
          call fillDustUniform(grid, grid%octreeRoot)
       endif
       

       if (nDustType >= 1) then
! Lots of dust types signifies using sphNG dump file with multigrain
          if (nDustType >= 12) then
             do i = 1, nDustType
                write(stringArray(i),'(a,i2.2)') "dustN",i
             enddo
             call writeVTKfile(grid,"dust.vtk",valueTypeString=stringArray(1:nDustType))
          else
             do i = 1, nDustType
                write(stringArray(i),'(a,i2.2)') "dust",i
             enddo
             call writeVTKfile(grid,"dust.vtk",valueTypeString=stringArray(1:nDustType))
          endif
       endif


       if (grid%geometry == "wr104") then
          call fillDustUniform(grid, grid%octreeRoot)
          call stripDustAway(grid%octreeRoot, 1.d-2, 1.d30)
       endif

       if (variableDustSublimation.or.quicksublimate) then
          call setupDustFracOriginal(grid%octreeRoot, nDustType)
          tauMax = 1.e-20
          call sublimateDust(grid, grid%octreeRoot, totFrac, nFrac, tauMax)
       endif


       nCellsInDiffusion = 0
       if (solveDiffusionZone) then
          call defineDiffusionOnRosseland(grid,grid%octreeRoot, taudiff, ndiff=nCellsInDiffusion)
          write(message,*) "Number of cells in diffusion zone: ", nCellsInDiffusion
          call writeInfo(message,IMPORTANT)
       endif
    endif

    iIter_grand = 0
    sweepSub = variableDustSublimation .and. (trim(sublimationMethod) == "sweep")
    adaptiveSublimation = variableDustSublimation .and. &
         ((trim(sublimationMethod) == "adaptive") .or. sweepSub)
    nDustSettled = 0
    converged = .false.
    dT_mean_new = 10.0d0
    dT_mean_old = 10.0d0
    iMultiplier = 1


    call zeroAdot(grid%octreeRoot)

    if (restartLucy) then
       open(33, file="restart.dat", status="unknown",form="formatted")
       read(33,'(a)') tfilename
       read(33, *) iHydro, iiter_grand, iMultiplier
       close(33)
    endif

    do while (.not.converged)
       ! ensure we do at least three iterations
       ! before removing the high temperature cells.

       !       if (variableDustSublimation) then
       !          nMonte = 10000000
       !          if (iIter_grand == 0) then
       !             nIter = 3
       !          else
       !             nIter = 1 !2
       !          endif
       !       endif

       nIter = 1


       do iIter = 1, nIter

          iIter_grand =  iIter_grand + 1  ! total number of iterations so far

          thisSmooth = .false.

          !          if (iIter_grand <= 3) thisSmooth = .true.

          !          if (variableDustSublimation) thisSmooth = .false.


          if (thisSmooth) then
             call locate(grid%lamArray, nLambda,lambdaSmooth, ismoothlam)
             call writeInfo("Smoothing adaptive grid structure for optical depth...", TRIVIAL)
             do j = iSmoothLam, nLambda, 2
                write(message,*) "Smoothing at lam = ",grid%lamArray(j), " angs"
                call writeInfo(message, TRIVIAL)
                if (grid%octreeRoot%twoD) then
                   do
                      gridConverged = .true.
                      call putTau(grid, grid%lamArray(j))
                      call myTauSmooth(grid%octreeRoot, grid, j, gridConverged, &
                           inheritProps = .false., interpProps = .true.)!, photosphereSplit = .true.)
                      if (gridConverged) exit
                   end do
                else
                   do
                      gridConverged = .true.
                      call myTauSmooth(grid%octreeRoot, grid, j, gridConverged, &
                           inheritProps = .false., interpProps = .true.)
                      if (gridConverged) exit
                   end do
                endif
             enddo
             call writeInfo("...grid smoothing complete", TRIVIAL)

             call writeInfo("Smoothing adaptive grid structure (again)...", TRIVIAL)
             do
                gridConverged = .true.
                call myScaleSmooth(smoothFactor, grid, &
                     gridConverged,  inheritProps = .false., interpProps = .true.)
                if (gridConverged) exit
             end do
             call countVoxels(grid%OctreeRoot,nOctals,nVoxels)  
             call writeInfo("...grid smoothing complete", TRIVIAL)

          endif

          nCellsInDiffusion = 0
          if (solveDiffusionZone) then
             call defineDiffusionOnRosseland(grid,grid%octreeRoot,tauDiff,  nDiff=nCellsInDiffusion)
          endif



          if (doTuning) call tune(6, "One Lucy Rad Eq Itr")  ! start a stopwatch

          call zeroDistanceGrid(grid%octreeRoot)

#ifdef CHEMISTRY
          if (doChemistry) then
             call zeroPathLength(grid%octreeRoot)
          endif
#endif


          call resetDirectGrid(grid%octreeRoot)

          nInf = 0
          nStarHit = 0
          nScat = 0
          nAbs = 0
          nDiffusion = 0
          nKilled = 0
          totalLumInPackets = 0.d0

          call countVoxels(grid%OctreeRoot,nOctals,nVoxels)  
          nInflow = 0
          call countInflow(grid%OctreeRoot,nInFlow)
          if (nLucy /= 0) then
             nMonte = nLucy
          else
             if (.not.variableDustSublimation) then
                nMonte = nVoxels * 10
             else
                nMonte = nVoxels * 50
             endif
          endif
          nMonte = nMonte * iMultiplier
          epsOverDeltaT = lCore / dble(nMonte)

          write(message,*) "Iteration",iIter_grand,",",nmonte," photons"
          call writeInfo(message, TRIVIAL)

          imonte_beg=1; imonte_end=nMonte  ! default value

#ifdef MPI
                 ! Set the range of index for a photon loop used later.     
                 np = nThreadsGlobal
                 n_rmdr = MOD(nMonte,np)
                 m = nMonte/np
          
                 if (myRankGlobal .lt. n_rmdr ) then
                    imonte_beg = (m+1)*myRankGlobal + 1
                    imonte_end = imonte_beg + m
                 else
                    imonte_beg = m*myRankGlobal + 1 + n_rmdr
                    imonte_end = imonte_beg + m -1
                 end if
             !    print *, ' '
             !    print *, 'imonte_beg = ', imonte_beg
             !    print *, 'imonte_end = ', imonte_end
            
                
                 !  Just for safety.
                 if (imonte_end .gt. nMonte .or. imonte_beg < 1) then
                    print *, 'Index out of range: i_beg and i_end must be ' 
                    print *, ' 0< index < ', nMonte , '    ... [lucy_mod::lucyRadiativeEquilibriumAMR]'
                    print *, 'imonte_beg = ', imonte_beg
                    print *, 'imonte_end = ', imonte_end
                    stop
                 end if



#endif    

                 if (.not. TorusSerial) call test_random_hybrid()

                 if (doTuning) call tune(6, "Photon loop")

                !$OMP PARALLEL DEFAULT(NONE) &
                !$OMP PRIVATE(iMonte, iSource, rVec, uHat, rHat) &
                !$OMP PRIVATE(escaped, wavelength, thisFreq, thisLam, iLam, octVec) &
                !$OMP PRIVATE(thisOctal, albedo, r) &
                !$OMP PRIVATE(vec_tmp, uNew, Treal, subcell, probDistJnu) &
                !$OMP PRIVATE(i, j, T1) &
                !$OMP PRIVATE(thisPhotonAbs, thisPhotonSca) &
                !$OMP PRIVATE( photonInDiffusionZone, leftHandBoundary, directPhoton) &
                !$OMP PRIVATE(diffusionZoneTemp, kappaAbsdb, sOctal, kappaScadb, kappaAbsPAH, kappaScaPAH, kAbsArray) &
                !$OMP PRIVATE(oldUHat, packetWeight, wavelengthWeight) &
                !$OMP PRIVATE(tempOctal, tempSubCell, temp, ok) &
                !$OMP PRIVATE(foundOctal, foundSubcell, hrecip_kt, logt, logNucritUpper, logNucritLower) &
                !$OMP PRIVATE(icritupper, icritlower,  kAbsArray2, hNuOverkT, fac2, this_bnu ) &
                !$OMP PRIVATE(thermalPhoton, scatteredPhoton, weight, PAHprob) &
                !$OMP PRIVATE(dStar, hitsStarSurface, hitStar) &
                !$OMP SHARED(logNu1, fac1dnu, loglam1, scalelam, scalenu)  &
                !$OMP SHARED(grid, nLambda, lamArray,miePhase, nMuMie, nDustType) &
                !$OMP SHARED(imonte_beg, imonte_end, source, nsource, usePAH) &
                !$OMP SHARED(dnu, nFreq, freq, nMonte, epsOverDeltaT) &
                !$OMP REDUCTION (+: nAbs, nScat, nInf, totalLumInPackets, nDiffusion, nKilled, nStarHit)

                 call returnKappa(grid, grid%OctreeRoot, 1, reset_kappa=.true.)

                 if (nSource > 0) then
                    call randomSource(source, nSource, &
                         i, packetWeight, grid%lamArray, grid%nLambda, initialize=.true.)  
                 endif


                !$OMP DO SCHEDULE(DYNAMIC,10)
                photonloop: do iMonte = imonte_beg, imonte_end
!                    if (mod(iMonte,imonte_end/10) == 0) write(*,*) "imonte ",imonte
#ifdef MPI
                   !  if (MOD(i,nThreadsGlobal) /= myRankGlobal) cycle photonLoop
#endif
                   

                   thisPhotonAbs = 0
                   thisPhotonSca = 0
                   call randomSource(source, nSource, iSource, packetWeight)
                   call getPhotonPositionDirection(Source(isource), rVec, uHat, rHat,grid, weight=weight)
                   packetWeight = packetWeight * weight
!                   write(*,*) isource, rVec%x*1.d10/rsol, rVec%y*1.d10/rsol, rvec%z*1.d10/rsol,packetweight
!                   write(*,*) "pos ",rVec,packetweight
                   thermalphoton = .true.
                   directPhoton = .true.
                   call amrGridValues(grid%octreeRoot, rVec, foundOctal=tempOctal, &
                        foundSubcell=tempsubcell)
                   thisOctal => tempOctal
                   subcell = tempSubcell

                   if (tempOctal%diffusionApprox(tempsubcell)) then

                      call randomWalk(grid, tempOctal, tempSubcell, thisOctal, Subcell, temp, ok)
                      if (.not.ok) cycle photonLoop ! abort photon if random walk has failed
                      directPhoton = .false.
                      rVec = subcellCentre(thisOctal, subcell)
                   endif
                   sOctal => thisOctal

                   escaped = .false.
                   call getWavelength(source(isource)%spectrum, wavelength, wavelengthWeight)
                   packetWeight = packetWeight * wavelengthWeight

                   thisFreq = cSpeed/(wavelength / 1.e8)
                   thislam = wavelength


                   do while(.not.escaped)


                      ilam = min(floor((log(thislam) - loglam1) * scalelam) + 1, nfreq)
                      ilam = max(ilam, 1)

                      ! distance to a (non-point) stellar surface along the current direction;
                      ! 1e30 if the packet will not hit one
                      call distanceToSource(source, nSource, rVec, uHat, hitsStarSurface, dStar)

                      call toNextEventAMR(grid, rVec, uHat, packetWeight, escaped, thisFreq, nLambda, lamArray,  &
                           photonInDiffusionZone, diffusionZoneTemp,  &
                           directPhoton, scatteredPhoton,  &
                           sOctal, foundOctal, foundSubcell, iLamIn=ilam, kappaAbsOut = kappaAbsdb, kappaScaOut = kappaScadb ,&
                           kappaAbsPAHout = kappaAbsPAH, kappaScaPAHout = kappaScaPAH, &
                           distToStar=dStar, hitStar=hitStar)

                      ! packets reaching a stellar surface are absorbed and removed
                      if (hitStar) then
                         nStarHit = nStarHit + 1
                         cycle photonloop
                      endif

                      If (escaped) then
                         !$OMP ATOMIC
                         nInf = nInf + 1
                         totalLumInPackets = totalLumInPackets + packetWeight*epsOverDeltaT
                      endif

                      if (photonInDiffusionZone) then
                         nDiffusion = nDiffusion + 1
                         cycle photonloop
                      endif

                      if (.not. escaped) then

                         thisOctal => foundOctal
                         subcell = foundSubcell
                         !                thisLam = (cSpeed / thisFreq) * 1.e8
                         !                call locate(lamArray, nLambda, real(thisLam), iLam)
                         octVec = rVec 

                         ! This call is dead because kappaabs and kappasca come from tonextevent amr

                         !                call amrGridValues(grid%octreeRoot, octVec, startOctal=thisOctal, actualsubcell=subcell,iLambda=iLam, &
                         !                     kappaSca=kappaScadb, kappaAbs=kappaAbsdb, grid=grid)
                         sOctal => thisOctal

!                         if (thisOctal%diffusionApprox(subcell)) then
!                            write(*,*) "photon in diffusion zone",photonindiffusionzone
!                         endif

                         if (kappaScadb+kappaAbsdb /= 0.0d0) then
                            albedo = (kappaScadb + kappaScaPAH) / (kappaScadb + kappaAbsdb + kappaAbsPAH)
                         else
                            albedo = 0.5
                         end if

                         directPhoton = .false.

                         ! photon is always absorbed/reprocessed if it has come from diffusion zone

                         if (PhotonInDiffusionZone) albedo = 0. 

                         albedo = min(albedo,0.9999d0)

                         call randomNumberGenerator(getDouble=r)

                         ! scattering case
                         if (r < albedo) then 

                            thermalPhoton = .false.
                            scatteredPhoton = .true.
                            vec_tmp = uhat
                            uNew = newDirectionMie(grid, thisOctal, subcell, vec_tmp, real(thisLam), lamArray, nLambda, miePhase, &
                                 nDustType, nMuMie)

                            nScat = nScat + 1
                            thisPhotonSca = thisPhotonSca + 1
                            if (thisPhotonSca > 50000) then
                               nKilled = nKilled + 1
                               cycle photonLoop
                            endif
                            uHat = uNew

                         else

                            nAbs = nAbs + 1
                            thisPhotonAbs = thisPhotonAbs + 1
                            if (thisPhotonAbs > 50000) then
                               nKilled = nKilled + 1
                               cycle photonLoop
                            endif

                            if (usepah) then
                               PAHprob = kappaAbsPAH / (kappaAbsdb + kappaAbsPAH)
                            else
                               PAHprob = 0.d0
                            endif
                            call randomNumberGenerator(getDouble=r)
                            if (r < PAHprob) then
                               thisFreq = getPAHfreqFromAdot(thisOctal%adot(subcell))
                               thisLam = (cSpeed / thisFreq) * 1.e8
                            else
                                  


                               call amrGridValues(grid%octreeRoot, octVec, startOctal=thisOctal, &
                                    actualSubcell=subcell, temperature=treal,grid=grid, kappaAbsArray=kAbsArray)
                               t1 = dble(treal)
                               
                               ! if the photon has come from the diffusion zone then it has to be
                               ! reprocessed using the appropriate temperature

                               if (photonInDiffusionZone) then
                                  t1 = dble(diffusionZoneTemp)
                               endif
                               
                               hrecip_kt = hcgs / (kErg * t1 )
                               
                               logt = log(t1)
                               logNucritUpper = 27.9500d0 + logt ! 23.76 is log(k) - log(h) ! 66.0
                               logNucritLower = 26.0626d0 + logt ! 23.76 is log(k) - log(h) ! 10.0
                               icritupper = min( floor((logNucritUpper - logNu1) * scalenu) + 1, nfreq        )
                               icritLower = min( floor((logNucritLower - logNu1) * scalenu) + 1, icritupper-1 )
                               !			    write(*,*) "icrit upper,lower", icritupper,icritlower,nfreq,t1
                               
                               do i = 1, icritupper
                                  iLam = nfreq - i + 1
                                  kAbsArray2(i) = kabsArray(ilam)
                               enddo
                               
                               !                  do i = 1, icritupper
                               !                     iLam = nfreq - i + 1
                               !                     write(*,*) kabsarray2(i), kabsarray(ilam)
                               !                  enddo
                               !stop
                               probDistJnu(1)  = 1.d-50
#ifdef USEMKL
                               hrecip_ktarray(1:icritupper) = hrecip_kt
                               
                               call vdmul(icritupper, freq(1:icritupper), hrecip_ktarray(1:icritupper) , hNuOverkT(1:icritupper)) ! start from second frequency 
                               call vdexp(icritlower, hNuOverkT(1:icritlower), hNuOverKt(1:icritlower)) ! hnuoverkt now e^hv/kT
                               call vdsub(icritlower, hnuoverkt(1:icritlower), OneArray(1:icritlower), hNuOverKt(1:icritlower))! hnuoverkt now e^nv/KT - 1
                               call vdinv(icritlower, hnuoverkt(1:icritlower), fac2(1:icritlower))
                               
                               call vdexp(icritupper - icritlower, -hNuOverkT(icritlower + 1 : icritupper), fac2(icritlower+1 : icritupper))
                               
                               call vdmul(icritupper, fac1dnu(1:icritupper), fac2(1:icritupper), hNuOverkT(1:icritupper)) ! hnuoverkt is temp array now
                               call vdmul(icritupper, hnuoverkt(1:icritupper), kabsarray2(1:icritupper), this_bnu(1:icritupper))
#else             
                               hNuOverkT(1:icritupper) = freq(1:icritupper) * hrecip_kt
                               
                               fac2(1:icritlower) = 1.d0 / (exp(hNuOverkT(1:icritlower)) - 1.d0)
                               fac2(icritlower+1:icritupper) = exp(-hNuOverkT(icritlower+1:icritupper))
                               
                               this_bnu(1:icritupper) = fac1dnu(1:icritupper) * fac2(1:icritupper) * kabsarray2(1:icritupper) 
                               !this_bnu here is actually bnu * dnu * kabs
#endif
                               
                               do i = 2, icritupper
                                  probDistJnu(i) = probDistJnu(i-1) + this_bnu(i) ! should be this_bnu(i-1)?
                               enddo
                               
                               if (probDistJnu(icritupper) /= 0.d0) then
                                  call randomNumberGenerator(getDouble=r)
                                  r = r * probdistjnu(icritupper) ! this is equivalent to dividing probdist (normalising)
                                  call locate(probDistJnu(1:icritupper), icritupper, r, j)
                                  if (j == nFreq) j = nFreq -1
                                  thisFreq = freq(j) + (freq(j+1) - freq(j))* &
                                       (r - probDistJnu(j))/(probDistJnu(j+1)-probDistJnu(j))
                                  ! note that probdistjnu(nfreq) cancels here so it's fine
                                  thisLam = (cSpeed / thisFreq) * 1.e8
                                  !                      write(*,*) cSpeed/thisFreq/angstromtocm
                               endif
                            endif
                            oldUhat = uHat
                            uHat = randomUnitVector()
                            thermalPhoton = .true.
                            scatteredPhoton = .false.
                            ! make sure diffused photon is moving out of diffusion zone
                               
                            if (photonInDiffusionZone) then
                               if ((oldUhat.dot.uHat)  > 0.) then
                                  uHat = (-1.d0) * uHat
                               endif
                            endif

                         endif

                      endif
                   enddo

                enddo photonLoop
                !$OMP END DO

          !$OMP END PARALLEL

                if (doTuning) call tune(6, "Photon loop")

#ifdef MPI
          ! Summing the value in octal computed by each processors.
          ! This is a recursive function which involves the many 
          ! comunications between processors. It may take a while....
          ! (Maybe faster to pack the values in 1D arrays and distribute.)

          if(doTuning) call tune(6, "  Lucy Loop Update ")  ! start a stopwatch
!          if(myRankIsZero) write(*,*) "Calling update_octal_MPI"
          !   call update_octal_MPI(grid%octreeRoot, grid)

          call MPI_BARRIER(MPI_COMM_WORLD, ierr)


          call updateGridMPI(grid)
!          if(myRankGlobal == 0) write(*,*) "Done update."

          if(doTuning) call tune(6, "  Lucy Loop Update ")  ! stop a stopwatch

          ! collect some statical info from each node.
          call MPI_ALLGATHER(REAL(nInf), 1, MPI_REAL, buffer_real, 1, MPI_REAL, MPI_COMM_WORLD, ierr)
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          nInf = INT(SUM(buffer_real), double)

          call MPI_ALLGATHER(REAL(nStarHit), 1, MPI_REAL, buffer_real, 1, MPI_REAL, MPI_COMM_WORLD, ierr)
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          nStarHit = INT(SUM(buffer_real), double)

          call MPI_ALLGATHER(REAL(nScat), 1, MPI_REAL, buffer_real, 1, MPI_REAL, MPI_COMM_WORLD, ierr)
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          nScat = INT(SUM(buffer_real), double)

          call MPI_ALLGATHER(REAL(nAbs), 1, MPI_REAL, buffer_real, 1, MPI_REAL, MPI_COMM_WORLD, ierr)
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          nAbs = INT(SUM(buffer_real), double)

          call MPI_ALLGATHER(totalLumInPackets, 1, MPI_DOUBLE, buffer_double, 1, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          totalLumInPackets = SUM(buffer_double)

          call MPI_ALLGATHER(REAL(nKilled), 1, MPI_REAL, buffer_real, 1, MPI_REAL, MPI_COMM_WORLD, ierr)
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          nKilled = INT(SUM(buffer_real))

          call MPI_ALLGATHER(REAL(nDiffusion), 1, MPI_REAL, buffer_real, 1, MPI_REAL, MPI_COMM_WORLD, ierr)
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          nDiffusion = INT(SUM(buffer_real), double)
#endif

          write(message,'(a,f7.2)') "Photons done.",real(ninf)/real(nmonte)
          call writeInfo(message,TRIVIAL)
          write(message,'(a,f13.3)') "Mean number of scatters per photon: ",real(nScat)/real(nMonte)
          call writeInfo(message,IMPORTANT)
          write(message,'(a,f13.3)') "Mean number of absorbs  per photon: ",real(nAbs)/real(nMonte)
          call writeInfo(message,IMPORTANT)
          write(message,'(a,f13.3)') "Fraction of photons killed: ",real(nKilled)/real(nMonte)
          call writeInfo(message,IMPORTANT)
          write(message,'(a,f13.3)') "Fraction of photons in diffusion zone: ",real(nDiffusion)/real(nMonte)
          call writeInfo(message,IMPORTANT)
          write(message,'(a,1pe13.3)') "Total luminosity in packets: ",totalLumInPackets/lsol
          call writeInfo(message,IMPORTANT)


          ! Packets absorbed by a star are removed: their energy is lost rather than shared out
          ! over the escaping packets, so they are included in the normalisation.
          write(message,'(a,f8.4)') "Fraction of packets absorbed by a star: ", &
               real(nStarHit)/real(max(nInf+nStarHit,1_bigint))
          call writeInfo(message,IMPORTANT)
          epsOverDeltaT = lCore / dble(nInf + nStarHit)

          nDT = 0
          nUndersampled = 0


          totalEmission = 0.

#ifdef CHEMISTRY
          if (doChemistry) then
             call pathLengthToIntensity(grid%octreeRoot, epsOverDeltat)
          endif
#endif

             call calculateMeanIntensity(grid%octreeRoot, epsOverDeltaT)
             call calculateAdotPAH(grid%octreeRoot, epsOverDeltaT)


          if (storeScattered) then
             call locate(freq, nFreq, cSpeed/(scatteredLightWavelength*angstromtocm),i)
             call calculateMeanIntensityOverdNu(grid%octreeRoot, epsOverDeltaT,dnu(i))
          endif

          

          call calculateTemperatureCorrections(.true., grid%octreeRoot, totalEmission, epsOverDeltaT, &
               nFreq, freq, dnu, lamarray, nLambda, grid, nDt, nUndersampled,  &
               dT_sum, dT_min, dT_max, dT_over_T_max)


          if (solveDiffusionZone) then
             if (doTuning) call tune(6, "Gauss-Seidel sweeps")
             call defineDiffusionOnRosseland(grid,grid%octreeRoot, taudiff)
          endif




          nCellsInDiffusion = 0
          if (solveDiffusionZone) then
             call defineDiffusionOnUndersampled(grid%octreeroot, nDiff=nCellsInDiffusion)
          endif

          if (solveDiffusionZone) then
             if ((.not.variableDustSublimation).or.adaptiveSublimation) then
                call solveArbitraryDiffusionZones(grid)
             else
                if (iIter_grand >= 5) call solveArbitraryDiffusionZones(grid)
             endif
             ! The diffusion solver sets only the cell temperature. Photon packets deposit little
             ! or no energy in diffusion cells, so the per-species temperatures computed in
             ! calculateTemperatureCorrections collapse towards TMinGlobal there. Those stale
             ! values are then used for the emissivity (calculateEtaCont) and for sublimation.
             ! In the optically thick diffusion limit all grain species share the local
             ! radiation-field temperature, so set them to the cell temperature.
             nColdDiff = 0
             call syncDustTempsInDiffusionZone(grid%octreeRoot, nColdDiff)
             write(message,'(a,i8)') "Diffusion cells with a dust species <50% of cell T (reset to cell T): ", nColdDiff
             call writeInfo(message, FORINFO)
          endif
       call writeVtkFile(grid, "diff.vtk", &
            valueTypeString=(/"rho        ", "temperature", "tau        ", "crossings  ", "etacont    " , &
            "dust       ", "deltaT     ", "etaline    ","fixedtemp  ",     "inflow     ", "diff       ", &
            "chiline    ", &
            "adot       ", "under      "/))


          nCellsInDiffusion = 0
          nUndersampled = 0
          if (writeoutput) write(*,*) "Checking for undersampled cells with mincrossings ",mincrossings
          call checkUndersampled(grid%octreeRoot, nUndersampled, nCellsInDiffusion)
          percent_undersampled  = 100.*real(nUndersampled)/real(nInflow-nCellsInDiffusion)

          call calculateEtaCont(grid, grid%octreeRoot, nFreq, freq, dnu, lamarray, nLambda, kAbsArray)
          dt_min = 1.d10
          dt_max = -1.d10
          dt_sum = 0.d0
          totalEmission = 0.d0
          call calculateDeltaTStats(grid%octreeRoot, dt_min, dt_max, dt_sum, dt_over_t_max, totalEmission, nDt)
          dT_mean_old = dT_mean_new  ! save this for the next iteration
          dT_mean_new = dT_sum/real(nDt)

          if (myRankIsZero) then

             write(message,*) "Number of photons entering diffusion zone", ndiffusion
             call writeInfo(message, TRIVIAL)
             write(message,*) "Mean change in temperature    : ", dT_mean_new , " Kelvins"
             call writeInfo(message, TRIVIAL)
             write(message,*) "Minimum change in temperature : ", dT_min, " Kelvins"
             call writeInfo(message, TRIVIAL)
             write(message,*) "Maximum change in temperature : ", dT_max, " Kelvins"
             call writeInfo(message, TRIVIAL)
             write(message,*) "Maximum fractional change in temperature : ", dT_over_T_max
             call writeInfo(message, TRIVIAL)
             write(message,*) "Fractional change in emissivity : ", abs(totalEmission-oldTotalEmission)/totalEmission
             call writeInfo(message, TRIVIAL)
             write(message,'(a,i8)') "Number of undersampled cells: ",nUndersampled
             call writeInfo(message, TRIVIAL)
             write(message,'(a,f8.2)') "Percentage of undersampled cells: ",percent_undersampled
             call writeInfo(message, TRIVIAL)
             write(message,*) "Emissivity of dust / core", totalEmission /lCore * 1.e30
             call writeInfo(message, TRIVIAL)
             write(message,*) "Total emission",totalEmission
             call writeInfo(message, TRIVIAL)


             !
             ! writing the info above to a file.
10           format(a3, a12, 3(2x, a12),     4(2x, a12),     (2x, a12))
11           format(3x, i12, 6(2x, f12.2), (2x,f12.4), (2x, i12))
             if (first_time_to_open_file) then
                open(unit=LU_OUT, file='convergence_lucy.dat', status='replace')
                first_time_to_open_file=.false.
                write(LU_OUT, 10) "#", "iteration", "Mean dT [K]", "Min dT [K]", "Max dT [K]", &
                     "% bad cells", "Eta dust/*", "Eta Total", "Max dT/T", "nMonte"
             else
                open(unit=LU_OUT, file='convergence_lucy.dat', status='old', position='append')
             end if

             write(LU_OUT, 11) iIter_grand, dT_sum/real(nDt), dT_min, dT_max, & 
                  percent_undersampled, totalEmission /lCore *1.e30, &
                  totalEmission, dT_over_T_max, nMonte

             close(LU_OUT)

          end if


          if (doTuning) call tune(6, "One Lucy Rad Eq Itr")  ! stop a stopwatch



          nCellsInDiffusion = 0
          if (solvediffusionZone) then
                call defineDiffusionOnRosseland(grid,grid%octreeRoot,tauDiff,  nDiff=nCellsInDiffusion)
             !       call unsetOnDirect(grid%octreeRoot)
             write(message,*) "Number of cells in diffusion zone: ", nCellsInDiffusion
             call writeInfo(message,IMPORTANT)
             if (doTuning) call tune(6, "Gauss-Seidel sweeps")
          endif



          if (grid%geometry.eq."wr104") then
             totalMass = 0.
             call findTotalMass(grid%octreeRoot, totalMass)
             scaleFac = real(massEnvelope / totalMass)
             if (writeoutput) write(*,'(a,1pe12.5)') "Density scale factor: ",scaleFac
             call scaleDensityAMR(grid%octreeRoot, dble(scaleFac))
             call sublimateDustWR104(grid%octreeRoot)
          endif



       enddo

       if (((grid%geometry == "ppdisk").or.(grid%geometry=="warpeddisc")).and.(nDustType > 1)) then
          call fillDustUniform(grid, grid%octreeRoot)
       endif

       if (quickSublimate) call quickSublimateLucy(grid%octreeRoot, minLevel=1.d-6)




       if (adaptiveSublimation) then
          ! every iteration: move each cell's surviving dust fraction towards equilibrium
          ! with an adaptive per-cell step (see adaptiveSublimateDust in dust_mod)
          nSubCells = 0
          nSubChanged = 0
          nSubDamped = 0
          maxSubChange = 0.d0
          ! With the sweep, iterations alternate after the first: odd iterations "move"
          ! (the sweep places the fronts with the current backwarming factor B(theta), and
          ! the Monte Carlo update adjusts shadowed cells), even iterations "measure"
          ! (nothing changes the dust, so the temperatures are consistent with it, and the
          ! sweep measures B and updates the fit of B(theta)). See sweepSublimation.
          sweepMeasure = sweepSub .and. (iIter_grand > 1) .and. (mod(iIter_grand,2) == 0)
          if ((.not.sweepSub) .or. ((iIter_grand >= 3).and.(.not.sweepMeasure))) then
             call adaptiveSublimateDust(grid%octreeRoot, (iIter_grand == 1).and.(.not.sweepSub), &
                  nSubCells, nSubChanged, nSubDamped, maxSubChange)
          endif
          if (sweepSub) then
             call sweepSublimation(grid, source, nSource, iIter_grand == 1, sweepMeasure, &
                  nSweepChanged, nSweepExposed, nSweepFronts, nSweepFrontMoved, &
                  nSweepFrontConv, sweepFrontErr, nSweepMidUnconv)
             nSubChanged = nSubChanged + nSweepChanged
             ! settled when, on a measure iteration, at least 95% of the lit fronts, and all
             ! of those within 10 degrees of the midplane, have their radius determined to
             ! within sublimationTolerance; converged after two consecutive settled measure
             ! iterations
             if (sweepMeasure) then
                if ((nSweepFronts > 0).and.(real(nSweepFrontConv) >= 0.95*real(nSweepFronts)) .and. &
                     (nSweepMidUnconv == 0)) then
                   nDustSettled = nDustSettled + 1
                else
                   nDustSettled = 0
                endif
                write(subMessage,'(a,i3,a,i4,a,i4,a,i3,a,f6.2,a,i2,a)') "Sublimation (measure) iteration ",iIter_grand, &
                     ": fronts converged ",nSweepFrontConv," of ",nSweepFronts,", unconverged within 10 deg ", &
                     nSweepMidUnconv,", mean front error ",100.d0*sweepFrontErr,"%, settled ",nDustSettled,"/2"
             else
                write(subMessage,'(a,i3,a,i4,a,i4,a,i8,a,i8,a)') "Sublimation (move) iteration ",iIter_grand, &
                     ": fronts moved ",nSweepFrontMoved," of ",nSweepFronts,", cells changed by sweep ", &
                     nSweepChanged,", by Monte Carlo update ",nSubChanged-nSweepChanged
             endif
             call writeInfo(subMessage, FORINFO)
             if (sweepMeasure) call writeSweepFitInfo()
          else
             if ((iIter_grand > 1).and.(nSubCells > 0).and. &
                  (real(nSubChanged) <= 1.e-3*real(max(nSubCells,1)))) then
                nDustSettled = nDustSettled + 1
             else
                nDustSettled = 0
             endif
             write(subMessage,'(a,i4,a,i9,a,i9,a,i9,a,1pe10.3,a,i3)') "Adaptive sublimation iteration ",iIter_grand, &
                  ": cells ",nSubCells,"  changed(>0.01) ",nSubChanged,"  damped ",nSubDamped, &
                  "  max |ds| ",maxSubChange,"  settled count ",nDustSettled
             call writeInfo(subMessage, FORINFO)
          endif
       else if (variableDustSublimation) then
          totFrac = 0.
          nFrac = 0


          if (iiter_grand == 3) then
             taumax = 0.1
             call sublimateDust(grid, grid%octreeRoot, totFrac, nFrac, taumax, undercorrectopt=1d0)
          endif
          if ((iIter_grand > 3).and.(iIter_grand <=10)) then
            if (iIter_grand <= 6) then
               taumax = 0.1
            else  
               taumax = 1.e30
            endif
            if (mod(iIter_grand-4,2)==0) call sublimateDust(grid, grid%octreeRoot, totFrac, nFrac, tauMax, undercorrectopt=0.1d0)
         endif
         if (iiter_grand==12) call sublimateDust(grid, grid%octreeRoot, totFrac, nFrac, tauMax, undercorrectopt=1.d0)
          
!
!             if (iIter_grand == 6) then
!                tauMax = 1.e-3
!                call sublimateDust(grid, grid%octreeRoot, totFrac, nFrac, tauMax)
!             endif
!
!             if (iIter_grand == 9) then
!                tauMax = 1.e-2
!                call sublimateDust(grid, grid%octreeRoot, totFrac, nFrac, tauMax)
!             endif
!
!             if (iIter_grand == 12) then
!                tauMax = 1.e-1
!                call sublimateDust(grid, grid%octreeRoot, totFrac, nFrac, tauMax)
!             endif
!
!             
!             if (iIter_grand > 15) then
!                tauMax = 1.e30
!!                if (dustSettling) call fillDustSettled(grid)
!!                call setupOrigDustFraction(grid%octreeRoot)
!                call sublimateDust(grid, grid%octreeRoot, totFrac, nFrac, tauMax)
!             endif

             if ((iiter_grand == 11).and.doSmoothGridTau) then
                call locate(grid%lamArray, nLambda,lambdasmooth,ismoothlam)

                call writeInfo("Smoothing adaptive grid structure for optical depth...", TRIVIAL)
                do j = iSmoothLam, iSmoothLam
                   write(message,*) "Smoothing at lam = ",grid%lamArray(j), " angs"
                   call putTau(grid, grid%lamArray(j))
!                      call writeVtkFile(grid, "puttau.vtk", &
!                           valueTypeString=(/"rho        ", "temperature", "dust       ", &
!                                    "etaline    "/))
                   call writeInfo(message, TRIVIAL)
                   do
                      gridConverged = .true.

                      if (solveVerticalHydro) then
                         call myTauSmooth(grid%octreeRoot, grid, j, gridConverged, &
                              inheritProps = .false., interpProps = .true.)!, photosphereSplit = .true.)
                      else
                         call getSublimationRadius(grid, subRadius)
                         call myTauSmooth(grid%octreeRoot, grid, j, gridConverged, &
                              inheritProps = .false., interpProps = .true., photosphereSplit = .true., splitToRadius=2.d0*subRadius)
                      endif

                      if (gridConverged) exit
                   end do

!                   do
!                      gridConverged = .true.
!                      call unrefineThinCells(grid%octreeRoot, grid, j, nUnrefine, gridconverged)
!                      if (gridConverged) exit
!                   end do
                   
                enddo

                
                call countVoxels(grid%OctreeRoot,nOctals,nVoxels)  
                call writeInfo("...grid smoothing complete", TRIVIAL)

                call writeInfo("Smoothing adaptive grid structure (again)...", TRIVIAL)
                do
                   gridConverged = .true.
                   call myScaleSmooth(smoothFactor, grid, &
                        gridConverged,  inheritProps = .false., interpProps = .true.)
                   if (gridConverged) exit
                end do
                call writeInfo("...grid smoothing complete", TRIVIAL)

                taumax = 1.d30
                call fixDust(grid,grid%octreeRoot)                
!
             if (writeoutput) write(*,*) "Global Memory has ",humanReadableMemory(globalMemoryFootprint)
             call findTotalMemory(grid, totMem)
             if (writeoutput) write(*,*) "Full check  has ", humanReadableMemory(totMem)
!             if (writeoutput) write(*,*) "Max memory available ", humanReadableMemory(maxMemoryAvailable)
!



             call writeVtkFile(grid, "afterphotorefine.vtk", &
                  valueTypeString=(/"rho        ", "temperature", "dust       ", &
                                    "etaline    "/))


          endif

          if (iIter_grand >= 5) then
             nCellsInDiffusion = 0
             if (solveDiffusionZone) then
                call defineDiffusionOnUndersampled(grid%octreeroot, nDiff=nCellsInDiffusion, reset=.true.)
                call defineDiffusionOnRosseland(grid,grid%octreeRoot, taudiff, ndiff=nCellsInDiffusion)
                write(message,*) "Number of cells in diffusion zone: ", nCellsInDiffusion
                call writeInfo(message,IMPORTANT)
             endif
          endif
       endif


       if (grid%geometry.eq."wr104") then
          totalMass = 0.
          call findTotalMass(grid%octreeRoot, totalMass)
          scaleFac = real(massEnvelope / totalMass)
          if (writeoutput) write(*,'(a,1pe12.5)') "Density scale factor: ",scaleFac
          call scaleDensityAMR(grid%octreeRoot, dble(scaleFac))
          call sublimateDustWR104(grid%octreeRoot)
       endif


       if (abs(totalEmission-oldTotalEmission)/totalEmission < 1.d-2) converged = .true.



       oldTotalEmission = totalEmission

       if (percent_undersampled > percent_undersampled_min) then
          iMultiplier  = iMultiplier * 2
          if ( convergeOnUndersampled ) converged = .false.
       endif
       

       if ((grid%geometry == "shakara").and.variableDustSublimation) then
          call getSublimationRadius(grid, subRadius)
          write(message, '(a, i3, a, f8.2,a )') "End of lucy iteration ",iIter_grand,&
               ": Dust Sublimation radius is: ",(1.d10*subRadius/autocm)," au"
          call writeInfo(message, FORINFO)
          
          
          write(message, '(a, i3, a, f8.2,a )') "End of lucy iteration ",iIter_grand,&
               ": Theoretical Dust Sublimation radius is: ",(subRadius_monnier), " au"
          call writeInfo(message, FORINFO)

      endif



       if (sweepSub) then
          ! Sweep method: stop only on a measure iteration (its temperatures are
          ! consistent with the dust) once the sublimation front is determined. The
          ! emission criterion is not used: the emission always changes between a move and
          ! a measure iteration, and the front criterion already requires the dust to be
          ! settled.
          converged = sweepMeasure .and. (nDustSettled >= 2) .and. &
               (.not.(convergeOnUndersampled .and. (percent_undersampled > percent_undersampled_min)))
       endif
       if (iIter_grand < iterlucy) converged = .false.
       if (adaptiveSublimation) then
          ! dust must have settled (<0.1% of cells changing) for two consecutive iterations
          if (nDustSettled < 2) converged = .false.
       else
          if (variableDustSublimation.and.(iIter_grand < 11)) converged = .false.
       endif

       if (iIter_grand >= maxiterLucy) then
          write(message,'(a)') "Lucy loop exceeded max iterations. Forcing convergence"
          call writeInfo(message, IMPORTANT)
          converged = .true.
       endif

       ! forceLucyConv is set in the parameters.dat file if required
       if ( forceLucyConv ) then
          converged = .true.
          if (myRankIsZero) write(*,*) "FORCING CONVERGENCE FOR TESTS!!!!"
       end if

       if (multiLucyFiles) then
          write(tfilename, '(a,i3.3,i3.3,i3.3,a)') "lucy_",iHydro,iIter_grand,imultiplier,".vtk"
       else
          write(tfilename,'(a,i3.3,a)') "lucy_",iIter_grand,".vtk"
       endif


        nLucyValueTypes = 13 + nDustType
        allocate(lucyValueTypeString(1:nLucyValueTypes))
        lucyValueTypeString(1:13) = (/"rho        ", "temperature", "tau        ", "crossings  ", "etacont    " , &
           "dust       ", "deltaT     ", "etaline    ","fixedtemp  ", "inflow     ", "diff       ", &
           "chiline    ", "adot       " /)
        do i = 1, nDustType
         write(lucyValueTypeString(13+i), '(a8,i2.2)') 'dusttemp', i
        enddo

        call writeVtkFile(grid, tfilename, valueTypeString=lucyValueTypeString)
        deallocate(lucyValueTypeString)
       !    !
       !    ! Write grid structure to a tmp file.
       !    !


       if (myRankIsZero.and.writelucytmpfile) then
          write(tfilename, '(a,i3.3,a,i3.3,a,i3.3,a)') "lucy_",iHydro,"_",iIter_grand,"_",imultiplier,".dat"
          call writeAMRgrid(tfilename, .false., grid)
!          call writeAMRgrid("lucy_grid_tmp.dat", .false., grid)
       endif
       if (myrankIsZero) then
          open(33, file="restart.dat", status="unknown",form="formatted")
          write(33,'(a)') tfilename
          write(33, *) iHydro,iiter_grand, iMultiplier
          close(33)
       endif
    enddo

!    if (variableDustSublimation.and.thisIsFinalPass.and.doSmoothGridTau) then
!       call sublimateDust(grid, grid%octreeRoot, totFrac, nFrac, tauMax)
!       if ((nfrac /= 0).and.(writeoutput)) then
!          write(*,*) "Average absolute change in sublimation fraction: ",totFrac/real(nfrac)
!       endif
!       call locate(grid%lamArray, nLambda,lambdasmooth,ismoothlam)
!       
!       call writeInfo("Smoothing adaptive grid structure for optical depth...", TRIVIAL)
!       do j = iSmoothLam, nLambda, 2
!          write(message,*) "Smoothing at lam = ",grid%lamArray(j), " angs"
!          call writeInfo(message, TRIVIAL)
!          do
!             gridConverged = .true.
!             call putTau(grid, grid%lamArray(j))
!             call myTauSmooth(grid%octreeRoot, grid, j, gridConverged, &
!                  inheritProps = .false., interpProps = .true., photosphereSplit = thisIsFinalPass)
!             
!             if (gridConverged) exit
!          end do
!       enddo
!       call countVoxels(grid%OctreeRoot,nOctals,nVoxels)  
!       call writeInfo("...grid smoothing complete", TRIVIAL)
!
!       call writeInfo("Smoothing adaptive grid structure (again)...", TRIVIAL)
!       do
!          gridConverged = .true.
!          call myScaleSmooth(smoothFactor, grid, &
!               gridConverged,  inheritProps = .false., interpProps = .true.)
!          if (gridConverged) exit
!       end do
!       call writeInfo("...grid smoothing complete", TRIVIAL)
!       
!    endif

 


    if (grid%geometry == "shakara") then
       
       write(message, '(a, f8.2,a )') "End of lucy loop: Dust Sublimation radius is: ",(subRadius*1.e10/autocm), " au"
       call writeInfo(message, FORINFO)
       write(message, '(a,f8.2,a)') "Theoretical dust sublimation radius is ",subRadius_monnier, " au"
       call writeInfo(message,TRIVIAL)
    endif

        if (storescattered) then 
           call locate(freq, nFreq, cSpeed/(scatteredlightwavelength*angstromtocm),i)
           call calcIntensityFromGrid(grid%octreeRoot, epsOverDeltaT, dnu(i))
           if (writeoutput) call writeVTKfile(grid, "scattered.vtk", valueTypeString = (/"scattered"/))
        endif

  end subroutine lucyRadiativeEquilibriumAMR

  subroutine getSublimationRadius(grid, subRadius, temperature, sublimationTemp, density)
    use amr_mod, only: tauAlongPath
    type(GRIDTYPE) :: grid
    real(double), intent(out) :: subRadius
    real(double), optional :: temperature, density, sublimationTemp
    real(double) :: tau
    type(VECTOR) :: point
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double), allocatable :: tauArray(:), xArray(:)

    allocate(tauArray(1:100000), xArray(1:100000))

   call tauAlongPath(1, grid, VECTOR(0.d0, 0.d0, 1.d-3*grid%halfSmallestSubcell), VECTOR(1.d0, 0.d0, 0.d0), &
         tau,  tauMax=0.667d0, subRadius=subRadius)

   if (tau >= 0.667d0) then
       if (present(temperature)) then
          point = VECTOR(subRadius, 0.d0, 0.d0)
          call findSubcellTD(point, grid%octreeRoot, thisOctal, subcell)
          temperature = thisOctal%temperature(subcell)
          density = thisOctal%rho(subcell)
          sublimationTemp = max(700.d0,2000.d0 *density**(1.95d-2))
          if (temperature > (sublimationTemp + 100.0)) &
               call writeWarning("Inner edge too hot by >100K") 
       endif
    else
       call writeWarning("No tau_ross = 2/3 boundary")
    endif
    deallocate(tauArray, xArray)
  end subroutine getSublimationRadius

  subroutine getTauinMidPlane(grid)
    use amr_mod, only: tauAlongPath
    type(GRIDTYPE) :: grid
    real(double) :: tau
    integer :: i


    do i = 1, grid%nLambda
       call tauAlongPath(i, grid, VECTOR(0.d0, 0.d0, 1.d-3*grid%halfSmallestSubcell), VECTOR(1.d0, 0.d0, 0.d0), &
            tau)
       if (writeoutput) write(*,*) "Midplane tau at ",grid%lamArray(i),": ",tau
    enddo
  end subroutine getTauinMidPlane

!!$  subroutine lucyRadiativeEquilibrium(grid, miePhase, nDustType, nMuMie, nLambda, lamArray, temperature, nLucy)
!!$    use source_mod, only: random_direction_from_sphere
!!$    use phasematrix_mod, only: PHASEMATRIX, newDirectionMie
!!$    use grid_mod, only: getindices
!!$
!!$    type(GRIDTYPE) :: grid
!!$    type(VECTOR) :: uHat, uNew
!!$    type(VECTOR) :: rVec
!!$    integer :: nlucy
!!$    real(double) :: packetWeight
!!$    integer :: nDustType
!!$    integer,intent(in) :: nLambda, nMuMie
!!$    type(PHASEMATRIX):: miePhase(:,:,:)
!!$    real :: lamArray(:)
!!$    real(oct), allocatable :: distanceGrid(:, :, :)
!!$    real(oct) :: r
!!$    integer, parameter :: nFreq = 100
!!$    integer :: i, j, k
!!$    real(oct) :: freq(nFreq), dnu(nFreq), probDistPlanck(nFreq), probDistJnu(nFreq)
!!$    real(oct) :: temperature
!!$    integer :: nMonte, iMonte, nScat, nAbs
!!$    real(oct) :: thisFreq,  logFreqStart, logFreqEnd
!!$    real(oct) :: albedo
!!$    logical :: escaped
!!$    integer :: i1, i2, i3
!!$    real(oct) :: t1, t2, t3
!!$    real(oct) :: thisLam
!!$    integer :: iLam
!!$    integer :: nInf
!!$    real(oct) :: t
!!$    integer :: nt
!!$    real(oct) :: ang
!!$    integer :: iIter, nIter = 5
!!$    real(oct) :: kabs
!!$    real(oct) ::  kappaP, norm
!!$    real(oct) :: adot, V, epsOverDeltaT
!!$    real(oct) :: newT, deltaT
!!$    real(oct) :: meanDeltaT
!!$    real(double) :: dummy(1)
!!$    real(oct), allocatable :: tempImage(:,:)
!!$    integer :: nDT
!!$    real(oct) :: totalEmission
!!$    type(vector) :: vec_tmp
!!$
!!$    allocate(tempImage(1:grid%nx,1:grid%ny))
!!$    dummy = 0.d0; probDistPlanck = 0.
!!$
!!$    logFreqStart = log10((cSpeed / (lamArray(nLambda)*1.e-8)))
!!$    logFreqEnd =  log10((cSpeed / (lamArray(1)*1.e-8)))
!!$    write(*,*) logFreqStart, logFreqEnd
!!$    do i = 1, nFreq
!!$       freq(i) = logFreqStart + (dble(i-1)/dble(nFreq-1))*(logFreqEnd-logFreqStart)
!!$       freq(i) = 10.**freq(i)
!!$    enddo
!!$    write(*,*) "Lam",(cSpeed/freq(1))*1.e8,(cSpeed/freq(nFreq))*1.e8
!!$
!!$    do i = 2, nFreq-1
!!$       dnu(i) = 0.5*((freq(i+1)+freq(i))-(freq(i)+freq(i-1)))
!!$    enddo
!!$    dnu(1) = freq(2)-freq(1)
!!$    dnu(nFreq) = freq(nFreq)-freq(nFreq-1)
!!$
!!$
!!$
!!$    call setupFreqProb(temperature, freq, dnu, nFreq, ProbDistPlanck)
!!$
!!$
!!$    write(*,'(a)') "Computing lucy radiative equilibrium..."
!!$
!!$
!!$    allocate(distanceGrid(1:grid%nx, 1:grid%ny, 1:grid%nz))
!!$
!!$
!!$    do iIter = 1, nIter
!!$       distanceGrid = 0.d0
!!$       write(*,*) "Iteration",iIter
!!$       nInf = 0
!!$       nScat = 0
!!$       nAbs = 0
!!$       nMonte = nLucy
!!$
!!$       do iMonte = 1, nMonte
!!$          escaped = .false.
!!$          call randomNumberGenerator(getDouble=r)
!!$          call locate(probDistPlanck, nFreq, r, j)
!!$          thisFreq = freq(j) + (freq(j+1) - freq(j))* &
!!$               (r - probDistPlanck(j))/(probDistPlanck(j+1)-probDistPlanck(j))
!!$
!!$          rVec  = real(grid%rCore, kind=oct) * randomUnitVECTOR()
!!$!          uHat = fromPhotosphereVector(rVec)
!!$          ! -- using a new routine in source_mod.f90 (RK)
!!$          uHat = random_direction_from_sphere(rVec)
!!$
!!$          if ( (rVec .dot. uHat) < 0.) uHat = (-1.d0) * uHat
!!$
!!$          do while(.not.escaped)
!!$
!!$             call toNextEvent(grid, rVec, uHat, packetWeight, escaped, distanceGrid, thisFreq, nLambda, lamArray)
!!$             if (escaped) then
!!$!$OMP ATOMIC
!!$                nInf = nInf + 1
!!$             endif
!!$             if (.not. escaped) then
!!$
!!$                thisLam = (cSpeed / thisFreq) * 1.e8
!!$                call locate(lamArray, nLambda, real(thisLam), iLam)
!!$
!!$                call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
!!$                if (.not.grid%oneKappa) then
!!$                   albedo = grid%kappaSca(i1,i2,i3,iLam)/(grid%kappaSca(i1,i2,i3,iLam)+grid%kappaAbs(i1,i2,i3,iLam))
!!$                else
!!$                   albedo = grid%oneKappaSca(1,iLam) / (grid%oneKappaSca(1,iLam)+grid%oneKappaAbs(1,iLam))
!!$                endif
!!$
!!$                call randomNumberGenerator(getDouble=r)
!!$                if (r < albedo) then
!!$
!!$                   vec_tmp=uhat 
!!$!                   uNew = newDirectionMie(grid, thisOctal, subcell, vec_tmp, real(thisLam), lamArray, nLambda, miePhase, nDustType, nMuMie)
!!$                   nScat = nScat + 1
!!$                   uHat = uNew
!!$                else
!!$
!!$                   nAbs = nAbs + 1
!!$                   probDistJnu(1) = 0.d0
!!$                   do i = 2, nFreq
!!$
!!$                      thisLam = (cSpeed / freq(i)) * 1.e8
!!$                      call hunt(lamArray, nLambda, real(thisLam), iLam)
!!$
!!$                      if ((ilam >=1).and.(ilam <= nlambda)) then
!!$                         if (.not.grid%oneKappa) then
!!$                            kabs = grid%kappaAbs(i1,i2,i3,iLam)
!!$                         else
!!$                            kabs = grid%oneKappaAbs(1,iLam) * grid%rho(i1,i2,i3)
!!$                         endif
!!$                         probDistJnu(i) = probDistJnu(i-1) + &
!!$                              (bnu(freq(i),dble(grid%temperature(i1,i2,i3)))) &
!!$                              * kabs * dnu(i)
!!$                      endif
!!$                   enddo
!!$
!!$                   probDistJnu(1:nFreq) = probDistJnu(1:nFreq) / probDistJnu(nFreq)
!!$                   call randomNumberGenerator(getDouble=r)
!!$                   call locate(probDistJnu, nFreq, r, j)
!!$                   if (j == nFreq) j = nFreq -1
!!$                   thisFreq = freq(j) + (freq(j+1) - freq(j))* &
!!$                        (r - probDistJnu(j))/(probDistJnu(j+1)-probDistJnu(j))
!!$                   uHat = randomUnitVector()
!!$                endif
!!$
!!$             endif
!!$          enddo
!!$       enddo
!!$
!!$       write(*,'(a,f7.2)') "Photons done.",real(ninf)/real(nmonte)
!!$       write(*,'(a,f7.3)') "Mean number of scatters per photon: ",real(nScat)/real(nMonte)
!!$       write(*,'(a,f7.3)') "Mean number of absorbs  per photon: ",real(nAbs)/real(nMonte)
!!$
!!$       V = (dble(grid%xAxis(2)-grid%xAxis(1))*dble(grid%yAxis(2)-grid%yAxis(1))*dble(grid%zAxis(2)-grid%zAxis(1)))
!!$
!!$       epsOverDeltaT = (dble(grid%lCore)) / dble(nInf)
!!$
!!$       meanDeltaT = 0.
!!$       nDT = 0
!!$
!!$       totalEmission = 0.
!!$       do i1 = 1, grid%nx
!!$          do i2 = 1, grid%ny
!!$             do i3 = 1, grid%nz
!!$                if (grid%inUse(i1,i2,i3)) then
!!$                   adot = epsoverDeltaT * (1.d0 / v) * distancegrid(i1,i2,i3) / 1.d30
!!$                   kappaP = 0.d0
!!$                   norm = 0.d0
!!$                   do i = 1, nFreq
!!$                      thisLam = (cSpeed / freq(i)) * 1.e8
!!$                      call hunt(lamArray, nLambda, real(thisLam), iLam)
!!$                      if ((iLam >=1) .and. (iLam <= nLambda)) then
!!$                         if (.not.grid%oneKappa) then
!!$                            kabs = grid%kappaAbs(i1,i2,i3,iLam)
!!$                         else
!!$                            kabs = grid%oneKappaAbs(1,iLam) * grid%rho(i1,i2,i3)
!!$                         endif
!!$
!!$                         kappaP = kappaP + kabs * &
!!$                              bnu(freq(i),dble(grid%temperature(i1,i2,i3)))  * dnu(i)
!!$                         norm = norm + bnu(freq(i),dble(grid%temperature(i1,i2,i3)))  * dnu(i)
!!$                      endif
!!$                   enddo
!!$                   kappaP = kappaP / norm /1.d10
!!$
!!$
!!$                   newT = (pi / stefanBoltz) * aDot / (fourPi * kappaP)
!!$                   newT = newT**0.25
!!$
!!$                   deltaT = newT - grid%temperature(i1,i2,i3)
!!$                   grid%temperature(i1,i2,i3) = max(1.e-3,grid%temperature(i1,i2,i3) + 0.8 * real(deltaT))
!!$                   nDT = nDT  + 1
!!$                   meanDeltaT = meanDeltaT + deltaT
!!$                   kappaP = 0.d0
!!$                   norm = 0.d0
!!$                   do i = 1, nFreq
!!$                      thisLam = (cSpeed / freq(i)) * 1.e8
!!$                      call hunt(lamArray, nLambda, real(thisLam), iLam)
!!$                      if ((iLam >=1) .and. (iLam <= nLambda)) then
!!$                         if (.not.grid%oneKappa) then
!!$                            kabs = grid%kappaAbs(i1,i2,i3,iLam)
!!$                         else
!!$                            kabs = grid%oneKappaAbs(1,iLam) * grid%rho(i1,i2,i3)
!!$                         endif
!!$                         kappaP = kappaP + kabs * &
!!$                              bnu(freq(i),dble(grid%temperature(i1,i2,i3))) * dnu(i)
!!$                         norm = norm + bnu(freq(i),dble(grid%temperature(i1,i2,i3))) * dnu(i)
!!$                      endif
!!$                   enddo
!!$                   kappaP = kappaP / norm /1.e10
!!$                   grid%etaCont(i1,i2,i3) = real(fourPi * kappaP * (stefanBoltz/pi) * (grid%temperature(i1,i2,i3)**4))
!!$                   totalEmission = totalEmission + grid%etaCont(i1,i2,i3) * V
!!$
!!$                endif
!!$             enddo
!!$          enddo
!!$       enddo
!!$
!!$       write(*,*) meanDeltaT / real(nDT)
!!$
!!$       write(*,*) "Emissivity of dust / core", totalEmission /grid%lCore * 1.e30
!!$
!!$       do i = 1, grid%nx
!!$          do j = 1, grid%ny
!!$             tempImage(i,j) = grid%temperature(i,j,grid%nz/2)
!!$          enddo
!!$       enddo
!!$    enddo
!!$    open(21,file="r.dat",status="unknown",form="formatted")
!!$    do i = 1, 100
!!$       r = grid%xAxis(grid%nx) * real(i-1)/99.
!!$       t = 0
!!$       nT = 0
!!$       do j = 1, 100
!!$          ang = twoPi * real(j-1)/100.
!!$          rVec = VECTOR(r*cos(ang), r*sin(ang),0.)
!!$          call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
!!$          if (grid%inUse(i1,i2,i3)) then
!!$             nT = nT + 1
!!$             t = t + grid%temperature(i1,i2,i3)
!!$          endif
!!$       enddo
!!$       if (nt > 0) then
!!$          t = t / real(nt)
!!$          write(21,*) r / grid%rCore, t
!!$       endif
!!$    enddo
!!$    close(21)
!!$
!!$    open(22,file="temps.dat",status="unknown",form="unformatted")
!!$    do i = 1, grid%nx
!!$       do j = 1, grid%ny
!!$          do k = 1, grid%nz
!!$             write(22) grid%temperature(i,j,k), grid%etaCont(i,j,k)
!!$          enddo
!!$       enddo
!!$    enddo
!!$    close(22)
!!$
!!$
!!$  end subroutine lucyRadiativeEquilibrium

  subroutine setupFreqProb(temperature, freq, dnu, nFreq, probDist)

    ! Lucy 1999, A&A, 344, 282 Equation 3

    real(oct) :: temperature
    integer :: nFreq
    real(oct) :: freq(:)
    real(oct) :: probDist(:)
    real(oct) :: dnu(:)

    integer :: i
    
    probDist(1:nFreq) = 0.
    do i = 2, nFreq
       probDist(i) = probDist(i-1) + bnu(dble(freq(i)),dble(temperature)) * dnu(i)
    enddo

    probDist(1:nFreq) = probDist(1:nFreq) / probDist(nFreq)

  end subroutine setupFreqProb


 subroutine toNextEvent(grid, rVec, uHat, packetWeight,  escaped, distanceGrid, thisFreq, nLambda, lamArray)
   use grid_mod, only: getindices

   type(GRIDTYPE) :: grid
   type(VECTOR) :: rVec, uHat
   integer :: i1, i2, i3
   real(oct) :: t1, t2, t3
   real(oct) :: tval, tau, r
   real :: lamArray(:)
   integer :: nLambda
   real(double) :: packetWeight
   logical :: stillinGrid
   logical :: escaped
   real(oct) :: thisTau
   real(oct) :: kabs, ksca
   real(oct) :: thisFreq
   real(oct) :: distanceGrid(1:grid%nx, 1:grid%ny, 1:grid%nz)
   real(oct) :: thisLam
   integer :: iLam

   tVal = 0.d0
    stillinGrid = .true.
    escaped = .false.


    thisLam = (cSpeed / thisFreq) * 1.e8
    call hunt(lamArray, nLambda, real(thisLam), iLam)
    if ((ilam < 1).or.(ilam > nlambda)) then
       write(*,*) "ilam errro",ilam
    endif

    call randomNumberGenerator(getDouble=r)
    tau = -log(1.0-r)
    call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
    call intersectCubeCart(grid, rVec, i1, i2, i3, uHat, tVal)

    if (.not.grid%oneKappa) then
       kabs = grid%kappaAbs(i1,i2,i3,iLam)
       ksca = grid%kappaSca(i1,i2,i3,iLam)
    else
       kabs = grid%oneKappaAbs(1,iLam) * grid%rho(i1,i2,i3)
       ksca = grid%oneKappaSca(1,iLam) * grid%rho(i1,i2,i3)
    endif

    thisTau = tVal * (kabs + ksca)

    do while(stillinGrid .and. (tau > thisTau)) 
       rVec = rVec + tVal * uHat


       if ((rVec%x > grid%xAxis(grid%nx)).or. &
           (rVec%y > grid%yAxis(grid%ny)).or. &
           (rVec%z > grid%zAxis(grid%nz)).or. &
           (rVec%x < grid%xAxis(1)).or. &
           (rVec%y < grid%yAxis(1)).or. &
           (rVec%z < grid%zAxis(1))) then
          stillinGrid = .false.
          escaped = .true.
          if (.not.grid%oneKappa) then
             kabs = grid%kappaAbs(i1,i2,i3,iLam)
          else
             kabs = grid%oneKappaAbs(1,iLam) * grid%rho(i1,i2,i3)
          endif
          distanceGrid(i1,i2,i3) = distanceGrid(i1,i2,i3) + tVal * kabs
       endif


       if (stillinGrid) then
          call randomNumberGenerator(getDouble=r)
          tau = -log(1.0-r)
          call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
          call intersectCubeCart(grid, rVec, i1, i2, i3, uHat, tVal)
          if (.not.grid%oneKappa) then
             kabs = grid%kappaAbs(i1,i2,i3,iLam)
             ksca = grid%kappaSca(i1,i2,i3,iLam)
          else
             kabs = grid%oneKappaAbs(1,iLam) * grid%rho(i1,i2,i3)
             ksca = grid%oneKappaSca(1,iLam) * grid%rho(i1,i2,i3)
          endif
          thisTau = tVal * (kabs+ksca)
          if (tVal == 0.) then
             escaped = .true.
             stillingrid = .false.
          endif
       endif
    enddo

    if (.not.escaped) then
       call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
       call intersectCubeCart(grid, rVec, i1, i2, i3, uHat, tVal)

       if ((tVal*tau/thisTau) > 2.*(grid%xAxis(2)-grid%xAxis(1))) then
          write(*,*) "tval*tau/thistau too big", tval*tau/thisTau
       else
          if (.not.grid%oneKappa) then
             kabs = grid%kappaAbs(i1,i2,i3,iLam)
          else
             kabs = grid%oneKappaAbs(1,iLam) * grid%rho(i1,i2,i3)
          endif
          distanceGrid(i1,i2,i3) = distanceGrid(i1,i2,i3) + (tVal*tau/thisTau) * kabs * packetWeight
       endif
       

       if (tau > thisTau) then
          write(*,*) "tau > thistau"
          stop
       endif
       
       rVec = rVec + (tVal*tau/thisTau) * uHat
    endif




 end subroutine toNextEvent


  subroutine intersectCubeCart(grid, posVec, i1,i2,i3,direction, tval)
   use vector_mod
   use grid_mod
   implicit none
   type(GRIDTYPE) :: grid
   type(VECTOR) :: direction
   type(VECTOR) :: posVec, norm(6), p3(6)
   real(oct) :: t(6),tval,denom(6)
   integer :: i,j
   logical :: ok, thisOk(6)
   integer :: i1, i2, i3

   ok = .true.

   norm(1) = VECTOR(1., 0., 0.)
   norm(2) = VECTOR(0., 1., 0.)
   norm(3) = VECTOR(0., 0., 1.)
   norm(4) = VECTOR(-1., 0., 0.)
   norm(5) = VECTOR(0., -1., 0.)
   norm(6) = VECTOR(0., 0., -1.)

   p3(1) = VECTOR(grid%xAxis(i1+1), 0., 0.)
   p3(2) = VECTOR(0.,grid%yAxis(i2+1),0.)
   p3(3) = VECTOR(0.,0.,grid%zAxis(i3+1))
   p3(4) = VECTOR(grid%xAxis(i1), 0., 0.)
   p3(5) = VECTOR(0.,grid%yAxis(i2),0.)
   p3(6) = VECTOR(0.,0.,grid%zAxis(i3))

   thisOk = .true.
   
   do i = 1, 6

      denom(i) = norm(i) .dot. direction
      if (denom(i) /= 0.) then
         t(i) = (norm(i) .dot. (p3(i)-posVec))/denom(i)
      else
         thisOk(i) = .false.
         t(i) = 0.
      endif
      if (t(i) < 0.) thisOk(i) = .false.
!      if (denom > 0.) thisOK(i) = .false.
 enddo




  
  j = 0
  do i = 1, 6
    if (thisOk(i)) j=j+1
  enddo

  if (j == 0) ok = .false.
   
  if (.not.ok) then
     write(*,*) i1, i2, i3
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  tval = minval(t, mask=thisOk)
  tval = max(tval * 1.001d0,dble((grid%xAxis(2)-grid%xAxis(1))/1000.))


  if (tval == 0.) then
     write(*,*) i1, i2, i3,tval
     write(*,*) posVec
     write(*,*) grid%xAxis(i1),grid%yAxis(i2),grid%zAxis(i3)
     write(*,*) grid%xAxis(i1+1),grid%yAxis(i2+1),grid%zAxis(i3+1)
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  if (tval > 2.*(grid%xAxis(2)-grid%xAxis(1))) then
!     write(*,*) "tval too big",tval,i1,i2,i3,posvec
!     write(*,*) "direction",direction
!     write(*,*) t(1:6)
!     write(*,*) denom(1:6)
  endif


  end subroutine intersectCubeCart 




  ! =====================================================================================
  ! Ray sweep for variable dust sublimation (sublimationmethod sweep)
  !
  ! The sublimation front for direct starlight is causal along rays from the star: the
  ! first dust along a ray is where a directly lit grain is cooler than T_sub, and dust
  ! behind it is shadowed. Updating dust only from the Monte Carlo temperatures moves the
  ! front ~1 optical depth per iteration (a shadowed cell only learns it is too hot once
  ! the dust in front of it has gone). The sweep instead visits cells in order of distance
  ! from the star, integrates the optical depth from the star using the dust already
  ! decided for closer cells, and for directly lit cells (tau < sweepTauCrit) sets the
  ! dust from the directly lit temperature B*T_thin, where T_thin is the analytic optically
  ! thin dust temperature and B a backwarming factor. Shadowed cells are left to the
  ! Monte Carlo (adaptiveSublimateDust), but keep at least the dust that B*T_thin allows.
  !
  ! B is measured along each ray from the Monte Carlo temperature of the last dust-free
  ! directly lit cell in front of the dust (T_MC/T_thin), which samples the radiation field
  ! at the face of the dust. (The cell average temperature of an unresolved, optically thick
  ! front cell is not used: it is dominated by the shadowed interior.) B(theta) is fitted,
  ! per dust type, as a smooth function symmetric about the midplane (see sweepSublimation).
  !
  ! The first source is taken as the centre for rays and angular bins.
  ! =====================================================================================

  real(double) function planckMeanDustKappa(grid, iDust, temperature, extinction)
    ! Planck mean of grid%oneKappaAbs (or abs+sca) of dust type iDust (includes 1e10 factor)
    type(GRIDTYPE) :: grid
    integer, intent(in) :: iDust
    real(double), intent(in) :: temperature
    logical, intent(in) :: extinction
    integer :: i
    real(double) :: nu1, nu2, numid, dnu, b, k, num, den
    num = 0.d0
    den = 0.d0
    do i = 2, grid%nLambda
       nu1 = cSpeed / (dble(grid%lamArray(i-1))*1.d-8)
       nu2 = cSpeed / (dble(grid%lamArray(i))*1.d-8)
       numid = 0.5d0*(nu1 + nu2)
       dnu = abs(nu1 - nu2)
       b = bnu(numid, temperature)
       k = 0.5d0*dble(grid%oneKappaAbs(iDust,i-1) + grid%oneKappaAbs(iDust,i))
       if (extinction) k = k + 0.5d0*dble(grid%oneKappaSca(iDust,i-1) + grid%oneKappaSca(iDust,i))
       num = num + k*b*dnu
       den = den + b*dnu
    enddo
    planckMeanDustKappa = 0.d0
    if (den > 0.d0) planckMeanDustKappa = num/den
  end function planckMeanDustKappa

  subroutine setupKPAbsTable(grid)
    ! Planck-mean absorption opacity of each dust type tabulated against temperature
    use inputs_mod, only : nDustType
    type(GRIDTYPE) :: grid
    integer :: i, j
    real(double) :: t
    if (allocated(kPAbsTable)) return
    allocate(kPAbsTable(nKPTable, nDustType))
    do j = 1, nDustType
       do i = 1, nKPTable
          t = 10.d0**(5.d0*dble(i-1)/dble(nKPTable-1))
          kPAbsTable(i,j) = planckMeanDustKappa(grid, j, t, .false.)
       enddo
    enddo
  end subroutine setupKPAbsTable

  subroutine setupSweepTables(grid, source, nSource)
    use inputs_mod, only : nDustType
    use source_mod, only : SOURCETYPE
    type(GRIDTYPE) :: grid
    integer :: nSource
    type(SOURCETYPE) :: source(:)
    integer :: j, k
    call setupKPAbsTable(grid)
    if (allocated(kPAbsSource)) return
    allocate(kPAbsSource(nDustType, nSource), kPExtSource(nDustType))
    allocate(sweepB(nSweepBins,nDustType), sweepFrontR(nSweepBins,nDustType))
    allocate(fitRelax(nDustType), fitLastSign(nDustType), sweepFitText(nDustType))
    sweepFitText = " "
    sweepB = 1.d0
    fitRelax = 1.d0
    fitLastSign = 0.d0
    sweepFrontR = -1.d0
    do j = 1, nDustType
       do k = 1, nSource
          kPAbsSource(j,k) = planckMeanDustKappa(grid, j, source(k)%teff, .false.)
       enddo
       kPExtSource(j) = planckMeanDustKappa(grid, j, source(1)%teff, .true.)
    enddo
  end subroutine setupSweepTables

  real(double) function kPAbsAt(iDust, temperature)
    integer, intent(in) :: iDust
    real(double), intent(in) :: temperature
    real(double) :: x
    integer :: i
    x = log10(max(1.d0, min(1.d5, temperature)))/5.d0*dble(nKPTable-1) + 1.d0
    i = min(nKPTable-1, max(1, int(x)))
    kPAbsAt = kPAbsTable(i,iDust) + (x - dble(i))*(kPAbsTable(i+1,iDust) - kPAbsTable(i,iDust))
  end function kPAbsAt

  real(double) function thinDustTemperature(position, iDust, source, nSource)
    ! Optically thin equilibrium temperature of dust type iDust at position:
    ! T^4 kP(T) = sum over sources of W T_eff^4 kP(T_eff), W the geometric dilution factor
    ! (stellar radius from luminosity and T_eff, so point sources are handled too).
    use source_mod, only : SOURCETYPE
    type(VECTOR), intent(in) :: position
    integer, intent(in) :: iDust, nSource
    type(SOURCETYPE) :: source(:)
    integer :: k, it
    real(double) :: d, rEff, w, rhs, lo, hi, mid, t
    rhs = 0.d0
    do k = 1, nSource
       d = modulus(position - source(k)%position)
       rEff = sqrt(source(k)%luminosity / (fourPi*stefanBoltz*source(k)%teff**4)) / 1.d10
       if (d > rEff) then
          w = 0.5d0*(1.d0 - sqrt(max(0.d0, 1.d0 - (rEff/d)**2)))
       else
          w = 0.5d0
       endif
       rhs = rhs + w * source(k)%teff**4 * kPAbsSource(iDust,k)
    enddo
    lo = 0.d0
    hi = 5.d0
    do it = 1, 60
       mid = 0.5d0*(lo + hi)
       t = 10.d0**mid
       if (t**4 * kPAbsAt(iDust, t) > rhs) then
          hi = mid
       else
          lo = mid
       endif
    enddo
    thinDustTemperature = 10.d0**(0.5d0*(lo + hi))
  end function thinDustTemperature

  recursive subroutine collectSweepCells(thisOctal, n, octs, subs, store)
    type(octal), pointer :: thisOctal
    type(octal), pointer :: child
    integer :: n
    type(octalWrapper) :: octs(:)
    integer :: subs(:)
    logical, intent(in) :: store
    integer :: subcell, i
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call collectSweepCells(child, n, octs, subs, store)
                exit
             end if
          end do
       else
          if (thisOctal%rho(subcell) < emptyCellRho) cycle
          if (.not.associated(thisOctal%origDustTypeFraction)) cycle
          n = n + 1
          if (store) then
             octs(n)%content => thisOctal
             subs(n) = subcell
          endif
       endif
    enddo
  end subroutine collectSweepCells

  real(double) function tauFromSource(grid, src, targetOctal, targetSub, tauStop)
    ! Optical depth (Planck mean extinction at the source T_eff) along the straight line
    ! from the source to where it enters the target cell. Stops early at tauStop.
    use inputs_mod, only : nDustType
    use source_mod, only : SOURCETYPE
    type(GRIDTYPE) :: grid
    type(SOURCETYPE) :: src
    type(octal), pointer :: targetOctal
    integer, intent(in) :: targetSub
    real(double), intent(in) :: tauStop
    type(VECTOR) :: cen, dir, pos
    type(octal), pointer :: thisOctal
    integer :: subcell, nStep, j
    real(oct) :: tVal
    real(double) :: dist, travelled, seg, kext, rStart

    tauFromSource = 0.d0
    cen = subcellCentre(targetOctal, targetSub)
    dir = cen - src%position
    dist = modulus(dir)
    if (dist <= 0.d0) return
    dir = dir / dist
    rStart = max(src%radius, 1.d-6*dist)
    if (rStart >= dist) return
    pos = src%position + rStart*dir
    travelled = rStart
    thisOctal => grid%octreeRoot
    nStep = 0
    do
       nStep = nStep + 1
       if (nStep > 1000000) exit
       call findSubcellLocal(pos, thisOctal, subcell)
       if (associated(thisOctal, targetOctal) .and. (subcell == targetSub)) exit
       call distanceToCellBoundary(grid, pos, dir, tVal, thisOctal)
       seg = dble(tVal) + 1.d-2*grid%halfSmallestSubcell
       if (thisOctal%rho(subcell) >= emptyCellRho) then
          kext = 0.d0
          do j = 1, nDustType
             kext = kext + thisOctal%dustTypeFraction(subcell,j) * kPExtSource(j)
          enddo
          tauFromSource = tauFromSource + seg * kext * thisOctal%rho(subcell)
          if (tauFromSource >= tauStop) exit
       endif
       travelled = travelled + seg
       if (travelled >= dist) exit
       pos = pos + seg*dir
    enddo
  end function tauFromSource

  subroutine sweepSublimation(grid, source, nSource, firstSweep, measureOnly, nChanged, nExposed, &
       nFronts, nFrontMoved, nFrontConverged, meanFrontErr, nMidUnconverged)
    ! Each dust type d has, in each angular bin, its own sublimation front and its own
    ! backwarming factor B(bin,d).
    ! Move iteration (measureOnly = .false.): in directly lit cells (optical depth from the
    !   star, all dust types, below sweepTauCrit) dust type d survives according to
    !   B(bin,d)*T_thin(d) < T_sub(d); shadowed cells keep at least that much dust.
    ! Measure iteration (measureOnly = .true.): no dust is changed, so the Monte Carlo
    !   temperatures are consistent with the dust. B(bin,d) is measured as T_MC(d)/T_thin(d)
    !   in the last unattenuated cell in front of dust type d's front. The ratios
    !   q = B_measured/B_trial of all lit fronts of dust type d are fitted with an even
    !   polynomial in theta (B smooth and symmetric about the midplane), and B is multiplied
    !   by the under-relaxed fit, with B >= 1.
    use inputs_mod, only : nDustType, sublimationTolerance
    use source_mod, only : SOURCETYPE
    use dust_mod, only : sublimationSurvivalFraction, localSublimationTemperature
    use utils_mod, only : sortdouble2index
    type(GRIDTYPE) :: grid
    integer :: nSource
    type(SOURCETYPE) :: source(:)
    logical, intent(in) :: firstSweep    ! first call: shadowed cells start with full dust
    logical, intent(in) :: measureOnly   ! measure iteration: do not change the dust
    integer, intent(out) :: nChanged, nExposed
    integer, intent(out) :: nFronts          ! (bin, dust type) pairs with a directly lit front
    integer, intent(out) :: nFrontMoved      ! move iteration: fronts that moved by more than a cell
    integer, intent(out) :: nFrontConverged  ! fronts whose B (and so front radius) is determined
    real(double), intent(out) :: meanFrontErr ! mean estimated relative uncertainty of the front radius
    integer, intent(out) :: nMidUnconverged  ! lit fronts within midplaneZone of the midplane not yet converged
    real(double), parameter :: midplaneZone = 10.d0  ! degrees
    real(double), parameter :: tauMeasure = 0.1d0    ! B is only measured where the starlight is unattenuated
    real(double), parameter :: dustyCellTau = 0.1d0  ! a dust type with optical depth above this in a cell is "present"
    integer, parameter :: nFitMin = 6                ! fewest measured fronts for a fit of B(theta)
    real(double), parameter :: relaxMin = 0.25d0
    type(octalWrapper), allocatable :: octs(:)
    integer, allocatable :: subs(:), idx(:)
    real(double), allocatable :: rr(:)
    real(double), allocatable :: bMeasR(:,:), bMeas(:,:), frontR(:,:), frontSize(:,:), bTrialUsed(:,:)
    logical, allocatable :: frontFound(:,:), anyDust(:,:)
    real(double), allocatable :: tThinCell(:), ratioCell(:)
    logical, allocatable :: presentCell(:)
    type(octal), pointer :: thisOctal
    type(VECTOR) :: cen, rel, posF, posF2
    integer :: n, i, ii, j, subcell, ib, nErr, k
    real(double) :: theta, tau, orig, sOld, sNew, slope, frontErr, errSum, thetaBin
    logical :: exposed, haveRatio
    real(double) :: xPt(nSweepBins), yPt(nSweepBins), coef(4), sigma, thetaMax, resid, sgn
    logical :: usePt(nSweepBins), fitOK(nDustType)
    real(double) :: fitQ(nSweepBins,nDustType), fitSig(nDustType), fitTmax(nDustType), fitCoef(4,nDustType)
    integer :: nPts, nUse, nCoef, iPass, fitNC(nDustType), fitNUse(nDustType)
    character(len=12) :: roughText

    call setupSweepTables(grid, source, nSource)

    n = 0
    allocate(octs(1), subs(1))
    call collectSweepCells(grid%octreeRoot, n, octs, subs, .false.)
    deallocate(octs, subs)
    allocate(octs(max(n,1)), subs(max(n,1)), rr(max(n,1)), idx(max(n,1)))
    i = 0
    call collectSweepCells(grid%octreeRoot, i, octs, subs, .true.)
    do i = 1, n
       rr(i) = modulus(subcellCentre(octs(i)%content, subs(i)) - source(1)%position)
       idx(i) = i
    enddo
    if (n > 1) call sortdouble2index(rr(1:n), idx(1:n))

    allocate(bMeasR(nSweepBins,nDustType), bMeas(nSweepBins,nDustType), frontR(nSweepBins,nDustType))
    allocate(frontSize(nSweepBins,nDustType), bTrialUsed(nSweepBins,nDustType))
    allocate(frontFound(nSweepBins,nDustType), anyDust(nSweepBins,nDustType))
    allocate(tThinCell(nDustType), ratioCell(nDustType), presentCell(nDustType))

    nChanged = 0
    nExposed = 0
    bMeasR = -1.d0
    bMeas = 0.d0
    frontFound = .false.
    anyDust = .false.
    frontR = -1.d0
    frontSize = 0.d0
    do ii = 1, n
       i = idx(ii)
       thisOctal => octs(i)%content
       subcell = subs(i)
       cen = subcellCentre(thisOctal, subcell)
       rel = cen - source(1)%position
       theta = atan2(rel%z, sqrt(rel%x**2 + rel%y**2))
       ib = min(nSweepBins, max(1, 1 + int((theta + 0.5d0*pi)/pi*dble(nSweepBins))))

       tau = tauFromSource(grid, source(1), thisOctal, subcell, sweepTauCrit)
       exposed = (tau < sweepTauCrit)
       if (exposed) nExposed = nExposed + 1

       haveRatio = .false.
       if (associated(thisOctal%tempStorage)) haveRatio = (size(thisOctal%tempStorage,2) >= 2*nDustType)
       presentCell = .false.
       do j = 1, nDustType
          orig = thisOctal%origDustTypeFraction(subcell,j)
          if (orig <= 0.d0) cycle
          sOld = thisOctal%dustTypeFraction(subcell,j) / orig
          sNew = sOld
          tThinCell(j) = 0.d0
          ratioCell(j) = 0.d0
          if (exposed) then
             tThinCell(j) = thinDustTemperature(cen, j, source, nSource)
             if (haveRatio) ratioCell(j) = thisOctal%tempStorage(subcell,nDustType+j)/tThinCell(j)
             if (.not.measureOnly) then
                ! directly lit temperature: thin temperature times backwarming. (No exp(-tau/4)
                ! attenuation: the diffuse field largely makes up for the attenuated starlight
                ! at tau < 1, and the attenuation put dust in lit cells whose Monte Carlo
                ! temperature was above T_sub, ahead of the front.) B*T_thin falls monotonically
                ! along the ray, so the front is the first cell where B*T_thin < T_sub.
                sNew = sublimationSurvivalFraction(sweepB(ib,j)*tThinCell(j), &
                     localSublimationTemperature(thisOctal, subcell, j))
                ! hysteresis: a cell the sweep empties is marked (subDir = 2) so that the
                ! Monte Carlo update cannot put dust back into it once it is shadowed
                if (associated(thisOctal%subStep)) then
                   thisOctal%subStep(subcell,j) = 1.d0
                   if (sNew < 0.5d0) then
                      thisOctal%subDir(subcell,j) = 2.d0
                   else
                      thisOctal%subDir(subcell,j) = 0.d0
                   endif
                endif
             endif
          else if (.not.measureOnly) then
             if (firstSweep) then
                sNew = 1.d0
             else if (sOld < 1.d0 - 1.d-6) then
                ! A shadowed cell is no hotter than it would be if directly lit, B*T_thin, so it
                ! keeps at least the dust that allows (the same rule as the lit cells). This fills, in the same move, the cells left empty behind a front that has
                ! moved inwards (the new front cell is optically thick even when nearly
                ! empty, so they are shadowed, and their Monte Carlo temperatures are stale).
                sNew = max(sOld, sublimationSurvivalFraction(sweepB(ib,j)*thinDustTemperature(cen, j, source, nSource), &
                     localSublimationTemperature(thisOctal, subcell, j)))
             endif
             ! a cell that is no longer directly lit is handed back to the Monte Carlo update,
             ! which may put dust back if it is cool enough
             if (associated(thisOctal%subDir)) then
                if (thisOctal%subDir(subcell,j) == 2.d0) then
                   thisOctal%subDir(subcell,j) = 0.d0
                   thisOctal%subStep(subcell,j) = 1.d0
                endif
             endif
          endif
          if (.not.measureOnly) then
             sNew = max(1.d-20, min(1.d0, sNew))
             thisOctal%dustTypeFraction(subcell,j) = max(1.d-30, orig*sNew)
             if (abs(sNew - sOld) > 1.d-2) nChanged = nChanged + 1
          endif
          ! dust type j is present in this cell if most of it survives or if its own optical
          ! depth is significant (near the front cells are so optically thick at full dust that
          ! even a small surviving fraction blocks the starlight)
          presentCell(j) = (sNew > 0.5d0) .or. &
               (thisOctal%dustTypeFraction(subcell,j)*kPExtSource(j)*thisOctal%rho(subcell)* &
               thisOctal%subcellSize > dustyCellTau)
       enddo

       ! Front of dust type j in this bin: its first present cell (lit or shadowed). B for
       ! dust type j is measured from the last unattenuated cell in front of it.
       do j = 1, nDustType
          if (anyDust(ib,j)) cycle
          if (thisOctal%origDustTypeFraction(subcell,j) <= 0.d0) cycle
          if (.not.presentCell(j)) then
             if (exposed .and. (tau < tauMeasure) .and. haveRatio) then
                bMeasR(ib,j) = rr(ii)
                bMeas(ib,j) = ratioCell(j)
             endif
          else
             anyDust(ib,j) = .true.
             if (exposed) frontFound(ib,j) = .true.
             frontR(ib,j) = rr(ii)
             frontSize(ib,j) = thisOctal%subcellSize
          endif
       enddo
    enddo

    nFronts = count(frontFound)
    nFrontMoved = 0
    if (.not.measureOnly) then
       do j = 1, nDustType
          do ib = 1, nSweepBins
             if (frontFound(ib,j) .neqv. (sweepFrontR(ib,j) > 0.d0)) then
                nFrontMoved = nFrontMoved + 1
             else if (frontFound(ib,j)) then
                if (abs(frontR(ib,j) - sweepFrontR(ib,j)) > frontSize(ib,j)) nFrontMoved = nFrontMoved + 1
             endif
             sweepFrontR(ib,j) = frontR(ib,j)
          enddo
       enddo
    endif

    ! Fit q = B_measured/B_trial over the lit fronts of each dust type as an even polynomial in
    ! theta (symmetric about the midplane), with one pass of 3-sigma clipping.
    bTrialUsed = sweepB
    fitQ = 1.d0
    fitOK = .false.
    fitSig = 0.d0
    fitTmax = 0.d0
    fitCoef = 0.d0
    fitNC = 0
    fitNUse = 0
    do j = 1, nDustType
       nPts = 0
       thetaMax = 0.d0
       do ib = 1, nSweepBins
          if (frontFound(ib,j) .and. (bMeasR(ib,j) > 0.d0)) then
             nPts = nPts + 1
             xPt(nPts) = abs(-0.5d0*pi + (dble(ib) - 0.5d0)*pi/dble(nSweepBins))
             yPt(nPts) = bMeas(ib,j)/bTrialUsed(ib,j)
             thetaMax = max(thetaMax, xPt(nPts))
          endif
       enddo
       if (nPts < nFitMin) cycle
       thetaMax = max(thetaMax, pi/180.d0)
       xPt(1:nPts) = xPt(1:nPts)/thetaMax
       nCoef = max(1, min(4, nPts/8))
       usePt(1:nPts) = .true.
       do iPass = 1, 2
          call fitEvenPoly(nPts, xPt, yPt, usePt, nCoef, coef)
          nUse = count(usePt(1:nPts))
          sigma = 0.d0
          do i = 1, nPts
             if (usePt(i)) sigma = sigma + (yPt(i) - evalEvenPoly(nCoef, coef, xPt(i)))**2
          enddo
          sigma = sqrt(sigma/dble(max(1, nUse - nCoef)))
          if (iPass == 1) then
             do i = 1, nPts
                resid = yPt(i) - evalEvenPoly(nCoef, coef, xPt(i))
                if (abs(resid) > 3.d0*sigma) usePt(i) = .false.
             enddo
          endif
       enddo
       fitOK(j) = .true.
       fitNC(j) = nCoef
       fitNUse(j) = count(usePt(1:nPts))
       fitCoef(1:nCoef,j) = coef(1:nCoef)
       fitSig(j) = sigma/sqrt(dble(max(1, fitNUse(j))))
       fitTmax(j) = thetaMax
       do ib = 1, nSweepBins
          thetaBin = abs(-0.5d0*pi + (dble(ib) - 0.5d0)*pi/dble(nSweepBins))
          fitQ(ib,j) = evalEvenPoly(nCoef, coef, min(thetaBin, thetaMax)/thetaMax)
       enddo
    enddo

    if (measureOnly) then
       do j = 1, nDustType
          if (.not.fitOK(j)) cycle
          ! under-relax when the midplane correction reverses (overshoot)
          sgn = sign(1.d0, fitCoef(1,j) - 1.d0)
          if (sgn*fitLastSign(j) < 0.d0) fitRelax(j) = max(relaxMin, 0.5d0*fitRelax(j))
          fitLastSign(j) = sgn
          do ib = 1, nSweepBins
             sweepB(ib,j) = max(1.d0, sweepB(ib,j)*(1.d0 + fitRelax(j)*(fitQ(ib,j) - 1.d0)))
          enddo
       enddo
    endif

    ! convergence diagnostics (written by the Lucy loop), one line per dust type: B at 0, 10 and 20 degrees, the fit
    ! (number of lit fronts used, angular range, midplane correction q(0)-1 and its
    ! statistical error, relaxation factor) and the front roughness within midplaneZone
    ! (rms relative radius jump between adjacent lit bins)
    if (measureOnly) then
       do j = 1, nDustType
          errSum = 0.d0
          nErr = 0
          do ib = 1, nSweepBins-1
             thetaBin = -90.d0 + dble(ib)*180.d0/dble(nSweepBins)
             if (abs(thetaBin) > midplaneZone) cycle
             if (frontFound(ib,j) .and. frontFound(ib+1,j)) then
                errSum = errSum + (2.d0*(frontR(ib+1,j) - frontR(ib,j))/(frontR(ib+1,j) + frontR(ib,j)))**2
                nErr = nErr + 1
             endif
          enddo
          roughText = "     -"
          if (nErr > 0) write(roughText,'(f5.1,a)') 100.d0*sqrt(errSum/dble(nErr)),"%"
          if (fitOK(j)) then
             write(sweepFitText(j),'(a,i2,a,3f7.3,a,i4,a,f5.1,a,f6.2,a,f5.2,a,f5.2,a,a)') "  dust ",j, &
                  ": B(0,10,20 deg)",(sweepB(nSweepBins/2 + 1 + k*nSweepBins/18, j), k = 0, 2), &
                  "  fit to ",fitNUse(j)," fronts within ",fitTmax(j)*180.d0/pi," deg, midplane correction ", &
                  100.d0*(fitCoef(1,j) - 1.d0),"% +/-",100.d0*fitSig(j),"%, step ",fitRelax(j), &
                  ", front roughness ",trim(roughText)
          else
             write(sweepFitText(j),'(a,i2,a,a)') "  dust ",j,": too few lit fronts to fit B, front roughness ", &
                  trim(roughText)
          endif
       enddo
    endif

    ! convergence: the fitted correction |q - 1| (plus its statistical error), times the
    ! sensitivity of the front radius to B (1/|dlnT_thin/dlnr| at the front), estimates the
    ! remaining error in the front radius
    nFrontConverged = 0
    nMidUnconverged = 0
    nErr = 0
    errSum = 0.d0
    do j = 1, nDustType
       do ib = 1, nSweepBins
          if (.not.frontFound(ib,j)) cycle
          thetaBin = -0.5d0*pi + (dble(ib) - 0.5d0)*pi/dble(nSweepBins)
          frontErr = 1.d0
          if (fitOK(j)) then
             posF = source(1)%position + frontR(ib,j)*VECTOR(cos(thetaBin), 0.d0, sin(thetaBin))
             posF2 = source(1)%position + 1.02d0*frontR(ib,j)*VECTOR(cos(thetaBin), 0.d0, sin(thetaBin))
             slope = log(thinDustTemperature(posF2, j, source, nSource) / &
                  thinDustTemperature(posF, j, source, nSource)) / log(1.02d0)
             frontErr = (abs(fitQ(ib,j) - 1.d0) + fitSig(j)) / max(0.05d0, abs(slope))
             errSum = errSum + frontErr
             nErr = nErr + 1
          endif
          if (frontErr < sublimationTolerance) then
             nFrontConverged = nFrontConverged + 1
          else if (abs(thetaBin)*180.d0/pi <= midplaneZone) then
             nMidUnconverged = nMidUnconverged + 1
          endif
       enddo
    enddo
    meanFrontErr = 0.d0
    if (nErr > 0) meanFrontErr = errSum/dble(nErr)
    deallocate(octs, subs, rr, idx, bMeasR, bMeas, frontR, frontSize, bTrialUsed)
    deallocate(frontFound, anyDust, tThinCell, ratioCell, presentCell)
  contains
    subroutine fitEvenPoly(np, x, y, use, nc, c)
      ! least squares fit of y = sum_k c(k) x^(2(k-1)), k = 1..nc, to the points with use
      integer, intent(in) :: np, nc
      real(double), intent(in) :: x(:), y(:)
      logical, intent(in) :: use(:)
      real(double), intent(out) :: c(4)
      real(double) :: a(4,4), rhs(4), f(4), piv
      integer :: ip, k1, k2, kr
      a = 0.d0
      rhs = 0.d0
      do ip = 1, np
         if (.not.use(ip)) cycle
         do k1 = 1, nc
            f(k1) = x(ip)**(2*(k1-1))
         enddo
         do k1 = 1, nc
            rhs(k1) = rhs(k1) + f(k1)*y(ip)
            do k2 = 1, nc
               a(k1,k2) = a(k1,k2) + f(k1)*f(k2)
            enddo
         enddo
      enddo
      ! Gaussian elimination (the normal matrix is symmetric positive definite)
      do k1 = 1, nc
         piv = a(k1,k1)
         if (abs(piv) < 1.d-30) piv = 1.d-30
         do kr = k1+1, nc
            f(1) = a(kr,k1)/piv
            a(kr,k1:nc) = a(kr,k1:nc) - f(1)*a(k1,k1:nc)
            rhs(kr) = rhs(kr) - f(1)*rhs(k1)
         enddo
      enddo
      c = 0.d0
      do k1 = nc, 1, -1
         c(k1) = (rhs(k1) - sum(a(k1,k1+1:nc)*c(k1+1:nc)))/max(1.d-30, a(k1,k1))
      enddo
    end subroutine fitEvenPoly

    real(double) function evalEvenPoly(nc, c, x)
      integer, intent(in) :: nc
      real(double), intent(in) :: c(4), x
      integer :: k1
      evalEvenPoly = 0.d0
      do k1 = 1, nc
         evalEvenPoly = evalEvenPoly + c(k1)*x**(2*(k1-1))
      enddo
    end function evalEvenPoly
  end subroutine sweepSublimation

  subroutine writeSweepFitInfo()
    ! per dust type convergence diagnostics of the last measure sweep
    integer :: j
    if (.not.allocated(sweepFitText)) return
    do j = 1, size(sweepFitText)
       call writeInfo(trim(sweepFitText(j)), FORINFO)
    enddo
  end subroutine writeSweepFitInfo

  ! Set the per-species dust temperatures (tempStorage columns nDustType+1:2*nDustType)
  ! and tDust to the cell temperature in cells treated with the diffusion approximation.
  ! nCold counts cells where any species was below half the cell temperature beforehand.
  recursive subroutine syncDustTempsInDiffusionZone(thisOctal, nCold)
    use inputs_mod, only : nDustType
    type(octal), pointer :: thisOctal
    type(octal), pointer :: child
    integer :: nCold
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call syncDustTempsInDiffusionZone(child, nCold)
                exit
             end if
          end do
       else
          if (.not.associated(thisOctal%diffusionApprox)) cycle
          if (.not.thisOctal%diffusionApprox(subcell)) cycle
          if (associated(thisOctal%tempStorage)) then
             if (size(thisOctal%tempStorage,2) >= 2*nDustType) then
                if (any(thisOctal%tempStorage(subcell,nDustType+1:2*nDustType) < &
                     0.5d0*dble(thisOctal%temperature(subcell)))) nCold = nCold + 1
                thisOctal%tempStorage(subcell,nDustType+1:2*nDustType) = dble(thisOctal%temperature(subcell))
             endif
          endif
          if (associated(thisOctal%tDust)) thisOctal%tDust(subcell) = dble(thisOctal%temperature(subcell))
       endif
    enddo
  end subroutine syncDustTempsInDiffusionZone

   recursive subroutine zeroDistanceGrid(thisOctal)
   use inputs_mod, only : nDustType
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call zeroDistanceGrid(child)
                exit
             end if
          end do
       else
          if ((.not.associated(thisOctal%tempStorage)).or.(size(thisOctal%tempStorage,1) /= thisOctal%maxChildren).or.&
               (size(thisOctal%tempStorage,2) < 2*nDustType)) then
             if (associated(thisOctal%tempStorage)) deallocate(thisOctal%tempStorage)
             allocate(thisOctal%tempStorage(1:thisOctal%maxChildren, 1:2*nDustType))
          endif
          thisOctal%distanceGrid(subcell) = 0.d0
          thisOctal%tempStorage(subcell,1:nDustType) = 0.d0
          thisOctal%tempStorage(subcell,nDustType+1:2*nDustType) = dble(thisOctal%temperature(subcell))
          thisOctal%nCrossings(subcell) = 0
          thisOctal%undersampled(subcell) = .false.
          thisOCtal%nDiffusion(subcell) = 0.
          if (associated(thisOctal%scatteredIntensity)) &
               thisOctal%scatteredIntensity(subcell,:,:) = 0.d0

          if (associated(thisOctal%meanIntensity)) &
               thisOctal%meanIntensity(subcell) = 0.d0

          if (associated(thisOctal%AdotPAH)) &
               thisOctal%aDotPAH(subcell) = 0.d0

       endif
    enddo
  end subroutine zeroDistanceGrid

  recursive subroutine fixDust(grid,thisOctal)
    use inputs_mod, only : nDusttype, grainfrac
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real(double) :: thisTau, tauMax, frac, kappaAbs, kappaSca
    integer :: subcell, i, j
    integer, save :: ilambda
    logical, save :: firstTime = .true.
    !$OMP THREADPRIVATE (ilambda, FirstTime)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fixDust(grid,child)
                exit
             end if
          end do
       else

          if (firstTime) then
             call locate(grid%lamArray, grid%nLambda, 5500., iLambda)
             firstTime = .false.
          endif
          do j = 1,nDustType
             if  (thisOctal%dustTypeFraction(subcell, j) < 0.999d0*grainFrac(j))  then
                call returnKappa(grid, thisOctal, subcell, ilambda, kappaSca=kappaSca, kappaAbs=kappaAbs)
                thisTau = (kappaAbs+kappaSca)*thisOctal%subcellSize
                tauMax = 1.d0
                if (thisTau > tauMax) then
                   write(*,*) "fixing ",thisOctal%dustTypeFraction(subcell,j),thisTau
                   frac = tauMax / thisTau 
                   thisOctal%dustTypeFraction(subcell,j) = thisOctal%dustTypeFraction(subcell,j) * frac
                endif
             endif
          enddo

       endif
    enddo
  end subroutine fixDust

  recursive subroutine zeroAdot(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call zeroAdot(child)
                exit
             end if
          end do
       else

               thisOctal%aDot(subcell) = 0.d0

       endif
    enddo
  end subroutine zeroAdot

  recursive subroutine countInflow(thisOctal,n)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i, n
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call countInflow(child, n)
                exit
             end if
          end do
       else
          if (thisOctal%inFlow(subcell).and.(thisOctal%rho(subcell) >= emptyCellRho)) n = n + 1
       endif
    enddo
  end subroutine countInflow

  recursive subroutine calculateMeanIntensity(thisOctal, epsOverDt)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: epsOverDt, dV
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateMeanIntensity(child, epsOVerDt)
                exit
             end if
          end do
       else

          if (associated(thisOctal%meanIntensity)) then
             dv = cellVolume(thisOctal, subcell)*1.d30
             thisOctal%meanIntensity(subcell) = (1.d0/fourPi) * (1.d0/dv) * (epsOverDt) * thisOctal%meanIntensity(subcell)
          endif
       endif
    enddo
  end subroutine calculateMeanIntensity

  recursive subroutine calculateAdotPAH(thisOctal, epsOverDt)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: epsOverDt, dV
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateadotPAH(child, epsOVerDt)
                exit
             end if
          end do
       else
          dv = cellVolume(thisOctal, subcell)*1.d30
          thisOctal%adotPAH(subcell) = (1.d0/dv) * (epsOverDt) * thisOctal%adotPAH(subcell)
       endif
    enddo
  end subroutine calculateAdotPAH

  recursive subroutine calculateMeanIntensityOverdNu(thisOctal, epsOverDt, dnu)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: epsOverDt, dV, dnu
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateMeanIntensityoverdNu(child, epsOVerDt, dnu)
                exit
             end if
          end do
       else

          if (associated(thisOctal%meanIntensity)) then
             dv = cellVolume(thisOctal, subcell)*1.d30
             thisOctal%meanIntensity(subcell) = (1.d0/fourPi) * (1.d0/dv) * (epsOverDt) * thisOctal%meanIntensity(subcell)/dnu
          endif
       endif
    enddo
  end subroutine calculateMeanIntensityOverdNu

  recursive subroutine checkUndersampled(thisOctal, nUndersampled, nCellsInDiffusion)
    use inputs_mod, only : minCrossings
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: nUndersampled, nCellsinDiffusion
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call checkUndersampled(child, nUndersampled, nCellsInDiffusion)
                exit
             end if
          end do
       else
          if (thisOctal%rho(subcell) < emptyCellRho) cycle   ! empty cells don't count either way
          if (.not.thisOctal%fixedTemperature(subcell)) then
             if (thisOctal%diffusionApprox(subcell)) then
                nCellsInDiffusion = nCellsInDiffusion + 1
             else
                if (thisOctal%inflow(subcell).and.(thisOctal%nCrossings(subcell) < minCrossings)) then
                   nUnderSampled = nUndersampled + 1
                endif
             endif
          endif
       endif
    enddo
  end subroutine checkUndersampled


  recursive subroutine calculateTemperatureCorrections(this_is_root, thisOctal, totalEmission, &
       epsOverDeltaT, nFreq, freq, dnu, lamarray, nLambda, grid, nDt, nUndersampled,  &
       dT_sum, dT_min, dT_max, dT_over_T_max)

   use inputs_mod, only : minCrossings, TMinGlobal, nDustType
    logical, intent(in) :: this_is_root    ! T if thisOctal is a root node.
    real(oct) :: totalEmission
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(gridtype) :: grid
    integer :: subcell, i
    real(oct) :: V, epsOverDeltaT
    real(oct) :: adot
    real(oct) :: kappaP, norm
    real(oct) :: thisLam
    integer :: ilam, nLambda
    real(oct) :: newT, deltaT
    integer, intent(inout) :: nDt
    real(oct), intent(inout) :: dT_sum ! [Kelvins]  the sum of the temperature change
    real(oct), intent(inout) :: dT_min ! [kelvins]  the minimum change of temperature
    real(oct), intent(inout) :: dT_max ! [kelvins]  the maximum change of temperature
    ! [kelvins]  the maximum fractional change of temperature
    real(oct), intent(inout) :: dT_over_T_max 

!    real :: kabs
    real(double) :: kabsArray(1000)
    real :: lamArray(:)
    integer :: nFreq
    real(oct) :: freq(:)
    real(oct) :: dnu(:)
    integer :: nUndersampled
    integer, save  :: nwarning = 0
    integer, parameter :: nmaxwarning = 20
    real, parameter :: underCorrect = 0.9
   integer :: iDust
   real(oct) :: fracDust, fracNorm
   real(oct) :: aDotDust, kappaPDust, normDust, newTDust, etaDust
   real(double) :: eqTarget, logTlo, logThi, logTmid, tMid
   integer :: iBisect
   real(oct) :: tDustEff, oldTDustEff

    if(this_is_root) then !initialize some values
       dT_sum = 0.0
       dT_min = 1.0e10
       dT_max = 0.0
       nDT = 0
       dT_over_T_max = 0.0
    end if
    

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateTemperatureCorrections(.false., child, totalEmission, epsOverDeltaT, &
                     nFreq, freq, dnu, lamarray, nLambda, grid, nDt, nUndersampled, &
                     dT_sum, dT_min, dT_max, dT_over_T_max)
                exit
             end if
          end do
       else

          if (.not.associated(thisOctal%oldTemperature)) then
             allocate(thisOctal%oldTemperature(1:thisOctal%maxChildren))
          endif
          if (.not.associated(thisOctal%tDust)) then
             allocate(thisOctal%tDust(1:thisOctal%maxChildren))
             thisOctal%tDust(1:thisOctal%maxChildren) = dble(thisOctal%temperature(1:thisOctal%maxChildren))
          endif
          thisOctal%oldTemperature(subcell) = thisOctal%temperature(subcell)

          if (thisOctal%inFlow(subcell)) then
             v = cellVolume(thisOctal, subcell)
             adot = epsoverDeltaT * (thisOctal%distancegrid(subcell)/v) / 1.d30
             kappaP = 0.d0
             norm = 0.d0

             call amrGridValues(grid%octreeRoot, subcellCentre(thisOctal,subcell), startOctal=thisOctal, &
                  actualSubcell=subcell, kappaAbsArray=kabsArray, grid=grid)

             do i = 1, nFreq
                thisLam = (cSpeed / freq(i)) * 1.d8
                call hunt(lamArray, nLambda, real(thisLam), iLam)
                if ((iLam >=1) .and. (iLam <= nLambda)) then

                   ! tjh added on 7/4/05 - stop direct addressing of octal kappa arrays
                   ! to allow implemntation of gas opacities...
                   
!                   call amrGridValues(grid%octreeRoot, subcellCentre(thisOctal,subcell), startOctal=thisOctal, &
!                        actualSubcell=subcell, kappaAbs=kabs, ilambda=ilam, grid=grid)

!                   if (.not.grid%oneKappa) then
!                      kabs = thisOctal%kappaAbs(subcell,iLam)
!                   else
!                      kabs = grid%oneKappaAbs(thisOctal%dustType(subcell),iLam) * thisOctal%rho(subcell)
!                   endif
                   kappaP = kappaP + dble(kabsArray(ilam)) * &
                        dble(bnu(dble(freq(i)),dble(thisOctal%temperature(subcell))))  * dble(dnu(i))
                   norm = norm + dble(bnu(dble(freq(i)),dble(thisOctal%temperature(subcell))))  * dble(dnu(i))
                endif
             enddo
             if (norm /= 0.d0) then
                kappaP = kappaP / norm /1.d10
             else
                kappaP = 0.d0
             endif
             
             
             if (kappaP /= 0.d0) then
                newT = dble(pi / stefanBoltz) * aDot / (dble(fourPi) * kappaP)
                newT = newT**0.25d0
             else
                newT = thisOctal%temperature(subcell)
             endif

             if (grid%oneKappa) then
                if ((.not.associated(thisOctal%tempStorage)).or.(size(thisOctal%tempStorage,2) < 2*nDustType)) then
                   if (associated(thisOctal%tempStorage)) deallocate(thisOctal%tempStorage)
                   allocate(thisOctal%tempStorage(1:thisOctal%maxChildren, 1:2*nDustType))
                   do iDust = 1, nDustType
                      thisOctal%tempStorage(1:thisOctal%maxChildren,nDustType+iDust) = &
                           dble(thisOctal%temperature(1:thisOctal%maxChildren))
                      thisOctal%tempStorage(1:thisOctal%maxChildren,iDust) = 0.d0
                   enddo
                endif

                ! Per-species dust temperature from radiative equilibrium,
                !    4 pi kappa_P(T) (sigma/pi) T^4 = adot,
                ! solved self-consistently in T (bisection in log T on the tabulated Planck
                ! mean). Evaluating kappa_P at the previous temperature instead makes the
                ! temperature lag the radiation field by several iterations whenever a cell's
                ! illumination changes (e.g. 10 K -> 10^4 K -> 1300 K -> 2000 K at the start).
                if (.not.allocated(kPAbsTable)) call setupKPAbsTable(grid)
                do iDust = 1, nDustType
                   fracDust = max(1.d-30, thisOctal%dustTypeFraction(subcell,iDust))
                   aDotDust = epsoverDeltaT * (thisOctal%tempStorage(subcell,iDust)/v) / 1.d30
                   ! target value of T^4 * kPAbsTable(T) (the table includes the 1e10 factor)
                   eqTarget = dble(pi/stefanBoltz) * aDotDust / (dble(fourPi) * fracDust * &
                        dble(thisOctal%rho(subcell)) / 1.d10)
                   if ((aDotDust > 0.d0) .and. (thisOctal%rho(subcell) > 0.d0)) then
                      logTlo = log10(max(1.d0, dble(TMinGlobal)))
                      logThi = 4.d0
                      do iBisect = 1, 50
                         logTmid = 0.5d0*(logTlo + logThi)
                         tMid = 10.d0**logTmid
                         if (tMid**4 * kPAbsAt(iDust, tMid) > eqTarget) then
                            logThi = logTmid
                         else
                            logTlo = logTmid
                         endif
                      enddo
                      newTDust = 10.d0**(0.5d0*(logTlo + logThi))
                   else
                      newTDust = dble(TMinGlobal)
                   endif
                   thisOctal%tempStorage(subcell,nDustType+iDust) = max(dble(TMinGlobal), min(newTDust, 10000.d0))
                enddo

                fracNorm = 0.d0
                tDustEff = 0.d0
                oldTDustEff = 0.d0
                do iDust = 1, nDustType
                   fracDust = max(1.d-30, thisOctal%dustTypeFraction(subcell,iDust))
                   tDustEff = tDustEff + fracDust * thisOctal%tempStorage(subcell,nDustType+iDust)
                   oldTDustEff = oldTDustEff + fracDust * dble(thisOctal%temperature(subcell))
                   fracNorm = fracNorm + fracDust
                enddo
                if (fracNorm > 0.d0) then
                   newT = tDustEff / fracNorm
                else
                   newT = thisOctal%temperature(subcell)
                endif
             endif

             deltaT = newT - thisOctal%temperature(subcell)
             if ((deltaT > 20000.).and.(thisOctal%nCrossings(subcell) .ge. minCrossings)) then
                if (nwarning < nmaxwarning) then
                   write(*,*) "Warning:: deltaT > 20000 kelvins calculateTemperatureCorrections.", &
                        deltaT, thisOctal%nCrossings(subcell), &
                        modulus(subcellcentre(thisOctal,subcell))/grid%rinner,thisOctal%rho(subcell), &
                        dble(1./grid%rinner)*subcellcentre(thisOctal,subcell) 
                   write(*,*) "kappap",kappap,"adot",adot
                   nwarning = nwarning + 1
                elseif (nwarning == nmaxwarning) then
                   write(*,*) " "
                   write(*,*) " ==> Suppressing the further delatT > 1000 Kelvins  wanrning .... "
                   write(*,*) " "
                   nwarning = nwarning + 1
                else
                   continue
                end if                   
                !                write(*,*) deltaT, thisOctal%temperature(subcell), newT
                !                write(*,*) adot,kappap,thisOctal%rho(subcell)
             endif


             if (.not.thisOctal%fixedTemperature(subcell)) then
                if (thisOctal%nCrossings(subcell) > 0) then
                   thisOctal%temperature(subcell) = max(TMinGlobal,thisOctal%temperature(subcell) + underCorrect*real(deltaT))
                else
                   thisOctal%temperature(subcell) = TminGlobal
                endif
                if (thisOctal%inflow(subcell).and.(thisOctal%nCrossings(subcell) .lt. minCrossings)) then
                   nUnderSampled = nUndersampled + 1
                   thisOctal%undersampled(subcell) = .true.
                endif

                thisOctal%temperature(subcell) = min(thisOctal%temperature(subcell), 10000.)
             endif


             kappaP = 0.d0
             norm = 0.d0
             if ((grid%oneKappa).and.associated(thisOctal%tempStorage)) then
                thisOctal%etaCont(subcell) = tiny(thisOctal%etaCont(subcell))
                do iDust = 1, nDustType
                   fracDust = max(1.d-30, thisOctal%dustTypeFraction(subcell,iDust))
                   kappaPDust = 0.d0
                   normDust = 0.d0
                   do i = 1, nFreq
                      thisLam = (cSpeed / freq(i)) * 1.e8
                      call hunt(lamArray, nLambda, real(thisLam), iLam)
                      if ((iLam >=1) .and. (iLam <= nLambda)) then
                         kappaPDust = kappaPDust + fracDust * dble(grid%oneKappaAbs(iDust,iLam)) * dble(thisOctal%rho(subcell)) * &
                              dble(bnu(dble(freq(i)),thisOctal%tempStorage(subcell,nDustType+iDust))) * dble(dnu(i))
                         normDust = normDust + &
                              dble(bnu(dble(freq(i)),thisOctal%tempStorage(subcell,nDustType+iDust))) * dble(dnu(i))
                      endif
                   enddo
                   if (normDust /= 0.d0) then
                      kappaPDust = kappaPDust / normDust / 1.d10
                   else
                      kappaPDust = 0.d0
                   endif
                   etaDust = fourPi * kappaPDust * (stefanBoltz/pi) * (thisOctal%tempStorage(subcell,nDustType+iDust)**4)
                   thisOctal%etaCont(subcell) = thisOctal%etaCont(subcell) + etaDust
                enddo
                thisOctal%tDust(subcell) = dble(thisOctal%temperature(subcell))
             else if (.not.grid%oneKappa) then
                kabsArray(1:nlambda) = thisOctal%kappaAbs(subcell,1:nlambda)
             else
                call amrGridValues(grid%octreeRoot, subcellCentre(thisOctal,subcell), startOctal=thisOctal, &
                     actualSubcell=subcell, kappaAbsArray=kabsArray, grid=grid)
             endif

             if (.not.((grid%oneKappa).and.associated(thisOctal%tempStorage))) then
                do i = 1, nFreq
                   thisLam = (cSpeed / freq(i)) * 1.e8
                   call hunt(lamArray, nLambda, real(thisLam), iLam)
                   if ((iLam >=1) .and. (iLam <= nLambda)) then
                      kappaP = kappaP + dble(kabsArray(ilam)) * &
                           dble(bnu(dble(freq(i)),dble(thisOctal%temperature(subcell)))) * dble(dnu(i)) 
                      norm = norm + dble(bnu(dble(freq(i)),dble(thisOctal%temperature(subcell)))) * dble(dnu(i)) 
                   endif
                enddo
                if (norm /= 0.d0) then
                   kappaP = kappaP / norm / 1.d10
                else
                   kappaP = TINY(kappaP)
                endif
                thisOctal%etaCont(subcell) = fourPi * kappaP * (stefanBoltz/pi) * &
                     (thisOctal%temperature(subcell)**4)
             endif
             totalEmission = totalEmission + thisOctal%etaCont(subcell) * V
             nDT = nDT  + 1
             dT_sum = dT_sum + ABS(deltaT)
             dT_min = MIN(deltaT, dT_min)
             dT_max = MAX(deltaT, dT_max)
             if (newT>0.0) dT_over_T_max = MAX(ABS(deltaT)/newT, dT_over_T_max)
          else 
             thisOctal%etaCont(subcell) = 0.
          endif
       endif
    enddo
  end subroutine calculateTemperatureCorrections

subroutine toNextEventAMR(grid, rVec, uHat, packetWeight,  escaped,  thisFreq, nLambda, lamArray,  &
      photonInDiffusionZone, diffusionZoneTemp,  directPhoton, scatteredPhoton, &
       startOctal, foundOctal, foundSubcell, ilamIn, kappaAbsOut, kappaScaOut, kappaAbsPAHout, kappaScaPAHout, &
       distToStar, hitStar)
   use diffusion_mod, only: randomwalk
   use inputs_mod, only : scatteredLightWavelength, storeScattered, usePAH, nDustType

#ifdef CHEMISTRY
   use chemistry_mod
   use inputs_mod, only : doChemistry
#endif

   type(GRIDTYPE) :: grid
   type(VECTOR) :: rVec,uHat, octVec
   type(OCTAL), pointer :: thisOctal, tempOctal !, sourceOctal
   type(OCTAL),pointer :: oldOctal, sOctal, topOctal
   type(OCTAL),pointer :: foundOctal, endOctal, startOctal
   real(double) :: packetWeight
   logical :: scatteredPhoton
   integer :: foundSubcell
   integer :: endSubcell
   integer :: subcell, tempSubcell!, sourceSubcell
   logical :: directPhoton
   real(oct) :: tval, tau, r
   real(double) :: tval_db
   real :: lamArray(:)
   integer :: nLambda
   logical :: stillinGrid, ok
   logical :: escaped
   real(double) :: kappaScaDb, kappaAbsDb, kappaAbsPAH, kappaScaPAH
   real(oct) :: thisTau
   real(double) :: tauRatio
   real(oct) :: thisFreq
   real(oct) :: thisLam
   integer :: iLam
   logical ::inFlow
   real :: diffusionZoneTemp
   logical :: photonInDiffusionZone
!   real(double) :: prob
   integer :: i
   integer :: iDust
   real(double) :: kappaAbsDustType
   real(double), parameter :: fudgeFac = 1.d-2
   integer, optional :: ilamIn
   real(double), optional :: kappaScaOut, kappaAbsOut, kappaAbsPAHout, kappaScaPAHout
   integer, save :: iLamScat
   logical, save :: firstTime = .true.
   ! distToStar: distance along uHat to a stellar surface (>= 1e30 if none). If the packet
   ! reaches it before interacting, the path is cut there and hitStar is returned .true.
   real(double), optional, intent(in) :: distToStar
   logical, optional, intent(out) :: hitStar
   real(double) :: pathSoFar
   logical :: reachesStar
    real :: test
    logical :: test2

!$OMP THREADPRIVATE(firstTime, iLamScat)


    test = lamArray(1)
    test2 = scatteredPhoton
   endSubcell = 0
   stillinGrid = .true.
   escaped = .false.
   ok = .true.
   photonInDiffusionZone = .false.
   kappaAbsDb = 0.d0; kappaScaDb = 0.d0
   topOctal => grid%octreeRoot
   pathSoFar = 0.d0
   reachesStar = .false.
   if (present(hitStar)) hitStar = .false.


   if (firstTime) then
      call hunt(lamArray, nLambda, scatteredLightWavelength, iLamScat)
      firstTime = .false.
   endif


   if(.not. present(ilamIn)) then
      thisLam = (cSpeed / thisFreq) * 1.e8
      if ((ilam < 1).or.(ilam > nlambda)) then
         write(*,*) "ilam error in tonexteventamr",ilam,thislam
      endif
   else
      ilam = ilamin
   endif

! select an initial random tau and find distance to next cell
   call randomNumberGenerator(getDouble=r)
   tau = -log(1.0-r)

    octVec = rVec

    call amrGridValues(grid%octreeRoot, octVec, iLambda=iLam,  startOctal=startOctal, foundOctal=thisOctal, &
         foundSubcell=subcell, kappaSca=kappaScadb, kappaAbs=kappaAbsdb, &
         grid=grid, inFlow=inFlow, direction=uHat)

    kappaAbsPAH = 0.d0
    kappaScaPAH = 0.d0
    if (usePAH)  then
       kappaAbsPAH = getKappaAbsPAH(thisFreq) * thisOctal%rho(subcell) * 1.d10
       kappaAbsPAH = getKappaScaPAH(thisFreq) * thisOctal%rho(subcell) * 1.d10
    endif

    if(present(kappaAbsOut)) kappaAbsOut = kappaAbsdb
    if(present(kappaScaOut)) kappaScaOut = kappaScadb
    if(present(kappaAbsPAHOut)) kappaAbsPAHOut = kappaAbsPAH
    if(present(kappaScaPAHOut)) kappaScaPAHOut = kappaScaPAH

    oldOctal => thisOctal
    sOctal => thisOctal

! moved from before call to amrgridvalues - th 16/11/05

    call distanceToCellBoundary(grid, rVec, uHat, tVal, sOctal)

    tval = tval + fudgeFac * grid%halfSmallestSubcell

    ! stop the step at a stellar surface if it lies within this cell
    if (present(distToStar)) then
       if (pathSoFar + dble(tVal) >= distToStar) then
          tVal = real(max(0.d0, distToStar - pathSoFar), oct)
          reachesStar = .true.
       endif
    endif

    if (inFlow) then
       thisTau = dble(tVal) * (kappaAbsdb + kappaScadb + kappaAbsPAH + kappaScaPAH)
    else
       thisTau = 1.0e-28
    end if

! if tau > thisTau then the photon has traversed a cell with no interactions

    do while(stillinGrid .and. (tau > thisTau)) 

! add on the distance to the next cell
       rVec = rVec + tVal * uHat
       octVec = rVec

       if (.not.inOctal(topOctal, octVec)) then
          stillinGrid = .false.
          escaped = .true.
       endif

! It is here that the photon path has potentially entered a cell
! in the diffusion approximation zone. In this case we have to return
! a position _before_ the photon enters the cell, and we need to
! tell the main photon loop that the photon has been absorbed
! and needs to re-emitted with a temperature corresponding to
! the top layer of the diffusion zone (this will conserve
! energy, and will create the appropriate boundary condition
! at the surface of the diffusion zone.
!       write(*,*) "octvec",octvec
!       write(*,*) sOctal%centre,sOctal%subcellSize
!       call amrGridValues(grid%octreeRoot, octVec,  startOctal=sOctal, foundOctal=tempOctal, &
!            foundSubcell=tempsubcell)

       if (.not.escaped) then

          if (.not.inOctal(topOctal, octVec)) then
             write(*,*) "photon outside grid but escaped is set to ",escaped,octVec
             write(*,*) grid%octreeRoot%subcellsize,grid%octreeRoot%xMax, inOctal(topOctal, octVec), &
                  inOctal(topOctal,rVec)
          endif

          tempOctal => sOctal
          call findSubcellLocal(octVec, tempOctal, tempSubcell)
!          call amrGridValues(grid%octreeRoot, octVec,  startOctal=sOctal, foundOctal=tempOctal, &
!            foundSubcell=tempsubcell)

       sOctal => tempOctal

       if (tempOctal%diffusionApprox(tempsubcell)) then
!XXXXXXXXX
          photoninDiffusionZone = .true.
          escaped = .true.
          goto 666

          call randomWalk(grid, tempOctal, tempSubcell,  endOctal, endSubcell, diffusionZoneTemp, ok)

          if (.not.ok) goto 666
          photonInDiffusionZone = .true.
          directPhoton = .false.
          rVec = subcellCentre(endOctal,endSubcell)
          octVec = rVec

          tempOctal => sOctal
          call findSubcellLocal(rVec, tempOctal, tempSubcell)
!          call amrGridValues(grid%octreeRoot, rVec,  startOctal=sOctal, foundOctal=tempOctal, &
!            foundSubcell=tempsubcell)

          sOctal => tempOctal
          if (tempOctal%diffusionApprox(tempsubcell)) then
             write(*,*) "moved to cell with diffapprox",tempOctal%diffusionApprox(tempsubcell)
             write(*,*) "endOctal",endOctal%diffusionApprox(endSubcell)
             write(*,*) "temp octal",subcellCentre(tempOctal, tempSubcell),tempSubcell, tempOctal%rho(tempSubcell)
             write(*,*) "end octal",subcellCentre(endOctal, endSubcell), endSubcell,endOctal%rho(endSubcell)
             write(*,*) tempOctal%chiLine(tempSubcell),endOctal%chiLine(endSubcell)
             write(*,*) tempOctal%chiLine(1:tempOctal%maxChildren)
             write(*,*) tempOctal%splitAzimuthally, endOctal%splitAzimuthally
             write(*,*) tempOctal%haschild(tempSubcell), endOctal%hasChild(endSubcell)
             write(*,*) tempOctal%nDepth,endOctal%nDepth
             write(*,*) tempOctal%phi*radtodeg,tempOctal%dphi*radtodeg
             write(*,*) endOctal%phi*radtodeg,endOctal%dphi*radtodeg
             write(*,*) tempOctal%r, endOctal%r
             write(*,*) tempOctal%parent%chiline(1:tempOctal%parent%maxChildren)
             write(*,*) endOctal%parent%chiline(1:endOctal%parent%maxChildren)
             write(*,*) "temp"
             do i = 1, tempOctal%parent%maxChildren
                write(*,*) subcellCentre(tempOctal,i)
             enddo
             write(*,*) "end"
             do i = 1, tempOctal%parent%maxChildren
                write(*,*) subcellCentre(tempOctal,i)
             enddo

!             do;enddo
          endif

       endif


    endif

! check whether the photon has escaped from the grid


! two cases here now. in the 2D case we only update the distance grid in a single plane (y=0, x>=0)

! update the distance grid

! ifort (v11 and earlier) has problems with type conversion in atomic statements.
    tval_db = real(tval,db)
!$OMP ATOMIC
    thisOctal%distanceGrid(subcell) = thisOctal%distanceGrid(subcell) + tVal_db * kappaAbsdb * packetWeight
   if (associated(thisOctal%tempStorage).and.grid%oneKappa) then
      do iDust = 1, nDustType
        kappaAbsDustType = thisOctal%dustTypeFraction(subcell, iDust) * grid%oneKappaAbs(iDust, iLam) * thisOctal%rho(subcell)
!$OMP ATOMIC
          thisOctal%tempStorage(subcell,iDust) = thisOctal%tempStorage(subcell,iDust) + tVal_db * kappaAbsDustType * packetWeight
      enddo
   endif
    if (usePAH) then
!$OMP ATOMIC
       thisOctal%aDotPAH(subcell) = thisOctal%adotPAH(subcell) + tVal_db * getkappaAbsPAH(thisFreq) * packetWeight
    endif


!$OMP ATOMIC
    thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1

#ifdef CHEMISTRY
    if (doChemistry) call storePathLength(thisOctal, subcell, thisFreq, tval_db * packetWeight)
#endif

    if (directPhoton) then
!$OMP ATOMIC
       thisOctal%nDirectPhotons(subcell) = thisOctal%nDirectPhotons(subcell)+1
    end if

! These lines were previously in an openmp critical section. If they are reinstated with OpenMP
! in use then access to shared memory may need to be protected from simutaneous updates. 
          if (storeScattered.and.(iLam==iLamScat)) &
               call addToScatteredIntensity(octVec, thisOctal, subcell, uHat, tVal)

    call addMeanIntensity(thisOctal, subcell, tVal*1.d10)


          if (storeScattered.and.(iLam==iLamScat)) &
               call addMeanIntensity(thisOctal, subcell, tVal*1.d10)




       ! the packet has now travelled this (possibly shortened) step; if it ended on a
       ! stellar surface the packet is absorbed by the star
       pathSoFar = pathSoFar + dble(tVal)
       if (reachesStar) then
          if (present(hitStar)) hitStar = .true.
          escaped = .false.
          goto 666
       endif

! now we need to return if the photon is in the diffusionzone

       if (photonInDiffusionZone) then
          goto 666
       endif
         

! now if the photon is in the grid choose a new random tau

       if (stillinGrid) then
          call randomNumberGenerator(getDouble=r)
          tau = -log(1.0-r)
!          tau = -log(r)  ! modified as above

          sOctal => thisOctal

          call amrGridValues(topOctal, octVec, iLambda=iLam,  foundOctal=thisOctal, startOctal=sOctal,&
               foundSubcell=subcell, kappaSca=kappaScadb, kappaAbs=kappaAbsdb, &
               grid=grid, inFlow=inFlow)

          if(present(kappaAbsOut)) kappaAbsOut = kappaAbsdb
          if(present(kappaScaOut)) kappaScaOut = kappaScadb

          kappaAbsPAH = 0.d0
          kappaScaPAH = 0.d0
          if (usePAH) then
             kappaAbsPAH = getKappaAbsPAH(thisFreq) * thisOctal%rho(subcell) * 1.d10
             kappaScaPAH = getKappaScaPAH(thisFreq) * thisOctal%rho(subcell) * 1.d10
          endif
          sOctal => thisOctal
          oldOctal => thisOctal

! calculate the distance to the next cell

          call distanceToCellBoundary(grid, rVec, uHat, tVal, sOctal)
          tval = tval + fudgeFac*grid%halfSmallestSubcell
          octVec = rVec

          ! stop the step at a stellar surface if it lies within this cell
          if (present(distToStar)) then
             if (pathSoFar + dble(tVal) >= distToStar) then
                tVal = real(max(0.d0, distToStar - pathSoFar), oct)
                reachesStar = .true.
             endif
          endif

! calculate the optical depth to the next cell boundary

          if (inFlow) then
             thisTau = dble(tVal) * (kappaAbsdb + kappaScadb + kappaAbsPAH)
          else
             thisTau = 1.0e-28
          end if

          if ((tVal == 0.0d0).and.(.not.reachesStar)) then
             escaped = .true.
             stillingrid = .false.
          endif
       endif ! still in grid

! if photon is still in grid and  tau > tau_to_the_next_cell then loop round again
! choosing a new tau

       
    enddo
    

! the photon may have escaped the grid...

    if (.not.inOctal(topOctal, octVec))  escaped = .true.

 ! if not the photon must interact in this cell
    if (.not.escaped) then
       octVec = rVec
!       if (.not.inOctal(grid%octreeRoot, octVec)) then
!          write(*,*) "Error:: Photon location is out of boundaries, but its status is not ESCAPED."
!          write(*,*) "        .... [lucy_mod::toNextEventAMR]"
!          write(*,*) "octVec-centre = ",octVec-grid%octreeRoot%centre
!          write(*,*) "cell size = ",grid%octreeRoot%subcellsize
!          stop
!       endif
       if (dble(tau)/thisTau .gt. 1.d0) then
          write(*,*) "tau prob",tau,thisTau
       endif

       call amrGridValues(topOctal, octVec, startOctal=oldOctal,iLambda=iLam, &
            foundOctal=thisOctal, foundSubcell=subcell, & 
            kappaAbs=kappaAbsdb,kappaSca=kappaScadb, grid=grid, inFlow=inFlow)

       if(present(kappaAbsOut)) kappaAbsOut = kappaAbsdb
       if(present(kappaScaOut)) kappaScaOut = kappaScadb

       kappaAbsPAH = 0.d0
       kappaScaPAH = 0.d0
       if (usePAH) then
          kappaAbsPAH = getKappaAbsPAH(thisFreq) * thisOctal%rho(subcell) * 1.d10
          kappaScaPAH = getKappaScaPAH(thisFreq) * thisOctal%rho(subcell) * 1.d10
       endif

       if (thisOctal%diffusionApprox(subcell)) then
          write(*,*) "Photon in diff zone but not escaped"
       endif

       if (.not.inFlow) then
          kappaAbsdb =0.0d0
          if(present(kappaAbsOut)) kappaAbsOut = kappaAbsdb
       endif


! update the distance grid

       if (thisTau > 0.d0) then


! ifort (v11 and earlier) has problems with type conversion in atomic statements.
    tval_db = real(tval,db)
    tauRatio = real( (tau/thisTau) ,db)
!$OMP ATOMIC
    thisOctal%distanceGrid(subcell) = thisOctal%distanceGrid(subcell) + (tVal_db*tauRatio) * kappaAbsdb * packetWeight
   if (associated(thisOctal%tempStorage).and.grid%oneKappa) then
      do iDust = 1, nDustType
        kappaAbsDustType = thisOctal%dustTypeFraction(subcell, iDust) * &
             grid%oneKappaAbs(iDust, iLam) * thisOctal%rho(subcell)
!$OMP ATOMIC
          thisOctal%tempStorage(subcell,iDust) = thisOctal%tempStorage(subcell,iDust) + &
               (tVal_db*tauRatio) * kappaAbsDustType * packetWeight
      enddo
   endif
    if (usePAH) then
!$OMP ATOMIC
       thisOctal%adotPAH(subcell) = thisOctal%aDotPAH(subcell) + (tVal_db*tauRatio) * getKappaAbsPAH(thisFreq) * packetWeight
    endif

#ifdef CHEMISTRY
    if (doChemistry) call storePathLength(thisOctal, subcell, thisFreq, (tval_db*tauRatio) * packetWeight)
#endif

!$OMP ATOMIC
    thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1
    if (directPhoton) then
!$OMP ATOMIC
       thisOctal%nDirectPhotons(subcell) = thisOctal%nDirectPhotons(subcell)+1
    end if

! These lines were previously in an openmp critical section. If they are reinstated with OpenMP
! in use then access to shared memory may need to be protected from simutaneous updates. 
          if (storeScattered.and.(iLam==iLamScat)) &
               call addToScatteredIntensity(octVec, thisOctal, subcell, uHat, dble(tVal)*dble(tau)/thisTau)

               call addMeanIntensity(thisOctal, subcell, 1.d10*dble(tVal)*dble(tau)/thisTau)

          if (storeScattered.and.(iLam==iLamScat)) &
               call addMeanIntensity(thisOctal, subcell, 1.d10*dble(tVal)*dble(tau)/thisTau)





          oldOctal => thisOctal
          
       endif

       if (tau > thisTau) then
          write(*,*) "tau > thistau"
          stop
       endif

! move the requisite distance within the cell and return. Reduce tval slightly to ensure
! event is within the grid if we are in the root octal.
       tVal = tVal - 2.0 * fudgeFac * grid%halfSmallestSubcell
       rVec = rVec + (dble(tVal)*dble(tau)/thisTau) * uHat


! this is a workaround due to numerical problems with long pathlengths
! and small cells. needs fixing.

       tempOctal => sOctal
       call findSubcellLocal(rVec, tempOctal, tempSubcell)
!       call amrGridValues(grid%octreeRoot, rVec,  startOctal=sOctal, foundOctal=tempOctal, &
!            foundSubcell=tempsubcell)
       if (tempOctal%diffusionApprox(tempsubcell)) then
          photonInDiffusionZone = .true.
          call randomWalk(grid, tempOctal, tempSubcell,  endOctal, endSubcell, diffusionZoneTemp, ok)
          if (.not.ok) goto 666
          photonInDiffusionZone = .true.
          rVec = subcellCentre(endOctal,endSubcell)
          octVec = rVec
          call amrGridValues(topOctal, rVec,  startOctal=sOctal,foundOctal=tempOctal, &
            foundSubcell=tempsubcell)
          if (tempOctal%diffusionApprox(tempsubcell)) then
             write(*,*) "! error - photon produced in diffusion zone 2..."
             write(*,*) Tau,thistau,tval
          endif
       endif
       foundSubcell = tempSubcell
       foundOctal => tempOctal

    endif
    
666 continue

  end subroutine toNextEventAMR


#ifdef MPI
  !
  ! For MPI implementation.
  ! Updates distanceGrid and nCrossings from each octal
  !
  recursive subroutine update_octal_MPI(thisOctal, grid)
    use mpi_global_mod, only: nThreadsGlobal

#ifdef CHEMISTRY
    use krome_main
    use krome_user
    use inputs_mod, only : doChemistry
#endif

    use mpi
    implicit none

    type(gridtype) :: grid
    type(octal), pointer  :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    ! data space to store values from all processors
    real(double)  :: buffer_ncrossings(nThreadsGlobal)     
    real(double) :: buffer_distanceGrid(nThreadsGlobal) 
    real(double) :: buffer_adotPah(nThreadsGlobal) 
#ifdef CHEMISTRY
    real(double) :: buffer_pathlength(nThreadsGlobal) 
    integer :: j
#endif
    integer  :: ierr


    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call update_octal_MPI(child, grid)
                exit
             end if
          end do
       else
          ! 
          ! collecting the data from all the processors including itself.
          call MPI_ALLGATHER(thisOctal%distanceGrid(subcell), 1, MPI_REAL, &
               buffer_distanceGrid, 1, MPI_REAL, MPI_COMM_WORLD, ierr)  
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)

          ! collecting the data from all the processors including itself.
          call MPI_ALLGATHER(thisOctal%aDotPAH(subcell), 1, MPI_REAL, &
               buffer_adotPah, 1, MPI_REAL, MPI_COMM_WORLD, ierr)  
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)

#ifdef CHEMISTRY
          if (doChemistry) then
             do j = 1, krome_nPhotoBins
                call MPI_ALLGATHER(thisOctal%kromeIntensity(subcell,j), 1, MPI_REAL, &
                     buffer_pathlength, 1, MPI_REAL, MPI_COMM_WORLD, ierr)  
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                thisOctal%kromeIntensity(subcell,j) = sum(buffer_pathlength)
             enddo
          endif
#endif


          
          
          call MPI_ALLGATHER(DBLE(thisOctal%ncrossings(subcell)), 1, MPI_DOUBLE_PRECISION, &
               buffer_ncrossings, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)  
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          
          thisOctal%distanceGrid(subcell) = SUM(buffer_distanceGrid)
          thisOctal%aDotPAH(subcell) = SUM(buffer_adotpah)
          thisOctal%nCrossings(subcell) = INT(SUM(buffer_ncrossings), kind=bigint)
          
       endif
    enddo
  end subroutine update_octal_MPI


   subroutine updateGridMPI(grid)
      use inputs_mod, only : usePAH, nDustType
    use mpi
    implicit none

    type(gridtype) :: grid
    integer :: nOctals, nVoxels
    real(double), allocatable :: nCrossings(:)
    real, allocatable :: tempRealArray(:)
    real(double), allocatable :: distanceGrid(:),adotPAH(:), tempDoubleArray(:), meanIntensity(:)
   real(double), allocatable :: tempStorage(:,:), tempStorageReduced(:,:)
    real, allocatable :: nDiffusion(:)
    integer, parameter :: nTheta = 11, nPhi = 10
   integer :: ierr, nIndex, nIndexScattered, nTempStorage, nTempStorageAlloc

    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
    nOctals = 0
    nVoxels = 0
    call countVoxels(grid%octreeRoot,nOctals,nVoxels)
    allocate(nCrossings(1:nVoxels))
    allocate(nDiffusion(1:nVoxels))
    allocate(distanceGrid(1:nVoxels))
    allocate(adotPAH(1:nVoxels))
    allocate(meanIntensity(1:nVoxels))

    nIndex = 0
    nIndexScattered = 0
    nTempStorage = 2*nDustType
    nTempStorageAlloc = max(1, nTempStorage)
    adotPAH = 0.d0
    allocate(tempStorage(nVoxels, nTempStorageAlloc))
    allocate(tempStorageReduced(nVoxels, nTempStorageAlloc))
    tempStorage = 0.d0
    call packValues(grid%octreeRoot,nIndex,nIndexScattered, &
          distanceGrid,adotPAH,nCrossings,ndiffusion, meanIntensity, tempStorage, nTempStorage)



!    if (storeScattered) then
!       allocate(tempDoubleArray(nVoxels*nTheta*nPhi))
!       tempDoubleArray = 0.d0
!       call MPI_ALLREDUCE(scatteredIntensity,tempDoubleArray, &
!            nVoxels*nTheta*nPhi,MPI_DOUBLE_PRECISION,&
!            MPI_SUM,MPI_COMM_WORLD,ierr)
!       scatteredIntensity= tempDoubleArray 
!       deallocate(tempDoubleArray)
!    endif


    allocate(tempRealArray(nVoxels))
    allocate(tempDoubleArray(nVoxels))


    tempDoubleArray = 0.d0
    call MPI_ALLREDUCE(distanceGrid,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    distanceGrid = tempDoubleArray 

    if (usePAH) then
       tempDoubleArray = 0.d0
       call MPI_ALLREDUCE(adotPAH,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
            MPI_SUM,MPI_COMM_WORLD,ierr)
       adotPAH = tempDoubleArray 
    endif

    tempDoubleArray = 0.d0
    call MPI_ALLREDUCE(meanIntensity,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    meanIntensity = tempDoubleArray 

    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(nCrossings,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    nCrossings = tempDoubleArray 

    tempRealArray = 0.0
    call MPI_ALLREDUCE(nDiffusion,tempRealArray,nVoxels,MPI_REAL,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    nDiffusion = tempRealArray 

    if (nTempStorage > 0) then
       tempStorageReduced = 0.d0
       ! Only columns 1:nDustType (absorbed energy per dust type) are tallies to be summed.
       ! Columns nDustType+1:2*nDustType hold per-dust temperatures, which are identical on
       ! every rank; summing them would scale the temperatures by the number of ranks.
       ! The first nDustType columns are contiguous, so reduce just those elements.
       call MPI_ALLREDUCE(tempStorage,tempStorageReduced,nVoxels*nDustType,MPI_DOUBLE_PRECISION,&
            MPI_SUM,MPI_COMM_WORLD,ierr)
       tempStorage(:,1:nDustType) = tempStorageReduced(:,1:nDustType)
    endif
    
      deallocate(tempRealArray, tempDoubleArray)
     
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
    
    nIndex = 0
    nIndexScattered = 0
    call unpackValues(grid%octreeRoot,nIndex,nIndexScattered, &
          distanceGrid,adotPAH,nCrossings,nDiffusion, meanIntensity, tempStorage, nTempStorage)

    deallocate(nCrossings, nDiffusion, distanceGrid, adotPAH, meanIntensity)
      deallocate(tempStorage, tempStorageReduced)

  end subroutine updateGridMPI
#endif


  subroutine find2Doctal(rVec, grid, foundOctal, subcell)
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec, rotVec
    type(OCTAL), pointer :: foundOctal, resultOctal
    integer :: subcell
    real(oct) :: r

! This subroutine finds the location of an octal in the y=0, x>=0 plane
! that corresponds to the position vector rVec when it is rotated
! into the y=0 plane 


    rotVec%z = rVec%z     ! the new vector has the same "z" value

    r = sqrt(rVec%x**2 + rVec%y**2) ! distance from z-axis

    rotVec%y = 0.

    rotVec%x = r

! now we find the octal and subcell that rotVec lies in

    call findSubcellTD(rotVec, grid%ocTreeRoot, resultOctal, subcell)
    foundOctal => resultOctal

  end subroutine find2Doctal


  recursive subroutine resetDistanceGrid(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call resetDistanceGrid(child)
                exit
             end if
          end do
       else
          thisOctal%distanceGrid(subcell) = thisOctal%chiline(subcell) 
       endif
    enddo
  end subroutine resetDistanceGrid

  recursive subroutine resetDirectGrid(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call resetDirectGrid(child)
                exit
             end if
          end do
       else
          if (.not. associated(thisOctal%nDirectPhotons)) allocate (thisOctal%nDirectPhotons(thisOctal%maxChildren))
          thisOctal%nDirectPhotons(subcell) = 0
       endif
    enddo
  end subroutine resetDirectGrid

  recursive subroutine recountDiffusionCells(thisOctal, ncells)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  integer :: ncells
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call recountDiffusionCells(child, ncells)
                exit
             end if
          end do
       else
          if (thisOctal%diffusionApprox(subcell)) then
             ncells = ncells + 1
          endif
       endif
    enddo
  end subroutine recountDiffusionCells

  recursive subroutine calculateDeltaTStats(thisOctal, dt_min, dt_max, dt_sum, dt_over_t_max, totalEmission, nDt)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    integer :: nDt
    real(double) :: dt_min, dt_max, dt_sum, dt_over_t_max, dt, totalEmission

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call  calculateDeltaTStats(child, dt_min, dt_max, dt_sum, dt_over_t_max, totalEmission, nDt)
                exit
             end if
          end do
       else
          dt = thisOctal%temperature(subcell) - thisOctal%oldTemperature(subcell)
          dt_min = MIN(dt_min, dt)
          dt_max = MAX(dt_max, dt)
          dt_sum = dt_sum + abs(dt)
          if (thisOctal%temperature(subcell)>0.0) dT_over_T_max = &
               MAX(ABS(dT)/thisOctal%temperature(subcell), dT_over_T_max)
          totalEmission = totalEmission + thisOctal%etacont(subcell) * cellVolume(thisOctal, subcell)
          ndt = ndt + 1
       endif
    enddo
  end subroutine calculateDeltaTStats

  recursive subroutine calculateEtaCont(grid, thisOctal, nFreq, freq, dnu, lamarray, nLambda, kabsArray)
      use inputs_mod, only : nDustType
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real :: lamArray(:)
    real(double) :: kappap, norm, freq(:), dnu(:), kabsArray(:), thisLam
    integer :: j, iLam, nLambda, nFreq
   integer :: iDust
   real(double) :: kappaPDust, normDust, fracDust

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call  calculateEtaCont(grid, child, nFreq, freq, dnu, lamarray, nLambda, kabsArray)
                exit
             end if
          end do
       else
             if ((grid%oneKappa).and.associated(thisOctal%tempStorage)) then
                thisOctal%etaCont(subcell) = tiny(thisOctal%etaCont(subcell))
                do iDust = 1, nDustType
                   fracDust = max(1.d-30, thisOctal%dustTypeFraction(subcell,iDust))
                   kappaPDust = 0.d0
                   normDust = 0.d0
                   do j = 1, nFreq
                      thisLam = (cSpeed / freq(j)) * 1.e8
                      call hunt(lamArray, nLambda, real(thisLam), iLam)
                      if ((iLam >=1) .and. (iLam <= nLambda)) then
                         kappaPDust = kappaPDust + fracDust * dble(grid%oneKappaAbs(iDust,iLam)) * dble(thisOctal%rho(subcell)) * &
                              dble(bnu(dble(freq(j)),thisOctal%tempStorage(subcell,nDustType+iDust))) * dble(dnu(j))
                         normDust = normDust + &
                              dble(bnu(dble(freq(j)),thisOctal%tempStorage(subcell,nDustType+iDust))) * dble(dnu(j))
                      endif
                   enddo
                   if (normDust /= 0.d0) then
                      kappaPDust = kappaPDust / normDust / 1.d10
                   else
                      kappaPDust = tiny(kappaPDust)
                   endif
                   thisOctal%etaCont(subcell) = thisOctal%etaCont(subcell) + fourPi * kappaPDust * (stefanBoltz/pi) * &
                        (thisOctal%tempStorage(subcell,nDustType+iDust)**4)
                enddo
             else
                kappaP = 0.d0
                norm = 0.d0
                if (.not.grid%oneKappa) then
                   kabsArray(1:nlambda) = thisOctal%kappaAbs(subcell,1:nlambda)
                else
                   call amrGridValues(grid%octreeRoot, subcellCentre(thisOctal,subcell), startOctal=thisOctal, &
                        actualSubcell=subcell, kappaAbsArray=kabsArray, grid=grid)
                endif

                do j = 1, nFreq
                   thisLam = (cSpeed / freq(j)) * 1.e8
                   call hunt(lamArray, nLambda, real(thisLam), iLam)
                   if ((iLam >=1) .and. (iLam <= nLambda)) then
                      kappaP = kappaP + dble(kabsArray(ilam)) * &
                           dble(bnu(dble(freq(j)),dble(thisOctal%temperature(subcell)))) * dble(dnu(j)) 
                      norm = norm + dble(bnu(dble(freq(j)),dble(thisOctal%temperature(subcell)))) * dble(dnu(j)) 
                   endif
                enddo
                if (norm /= 0.d0) then
                   kappaP = kappaP / norm / 1.d10
                else
                   kappaP = TINY(kappaP)
                endif
                thisOctal%etaCont(subcell) = fourPi * kappaP * (stefanBoltz/pi) * &
                     (thisOctal%temperature(subcell)**4)
             endif
       endif
    enddo
  end subroutine calculateEtaCont

  recursive subroutine packvalues(thisOctal,nIndex, nIndexScattered,&
     distanceGrid,adotPAH,nCrossings, nDiffusion, meanIntensity, tempStorage, nTempStorage)
    use inputs_mod, only : usePAH
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: distanceGrid(:)
  real(double) :: aDotPAH(:)
  real(double) :: meanIntensity(:)
  real(double) :: nCrossings(:)
  real(double) :: tempStorage(:,:)
  real :: nDiffusion(:)
  integer :: nIndex, nIndexScattered
  integer :: nTempStorage
  integer :: subcell, i !, j , k
!  integer, parameter :: nTheta = 10, nPhi = 10
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call packvalues(child,nIndex,nIndexScattered,distanceGrid,adotPAH, nCrossings,nDiffusion, &
                     meanIntensity, tempStorage, nTempStorage)
                exit
             end if
          end do
       else
          nIndex = nIndex + 1
          distanceGrid(nIndex) = thisOctal%distanceGrid(subcell)
          if (usePAH) adotPAH(nIndex) = thisOctal%adotPAH(subcell)
          meanIntensity(nIndex) = thisOctal%meanIntensity(subcell)
          nCrossings(nIndex) = dble(thisOctal%nCrossings(subcell))
          nDiffusion(nIndex) = thisOctal%nDiffusion(subcell)
          if (nTempStorage > 0) then
             if (associated(thisOctal%tempStorage).and.(size(thisOctal%tempStorage,2) >= nTempStorage)) then
                tempStorage(nIndex,1:nTempStorage) = thisOctal%tempStorage(subcell,1:nTempStorage)
             else
                tempStorage(nIndex,1:nTempStorage) = 0.d0
             endif
          endif
!          if (storeScattered) then
!             do j = 1, nTheta
!                do k = 1, nPhi
!                   nIndexScattered = nIndexScattered+1
!                   scatteredIntensity(nIndexScattered) = thisOctal%scatteredIntensity(subcell, j ,k)
!                enddo
!             enddo
!          endif

       endif
    enddo
  end subroutine packvalues

   recursive subroutine unpackvalues(thisOctal,nIndex,nIndexScattered,distanceGrid,adotPAH,nCrossings, nDiffusion, meanIntensity, &
          tempStorage, nTempStorage)
    use inputs_mod, only : usePAH
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: distanceGrid(:), adotPAH(:), meanIntensity(:)
  real(double) :: ncrossings(:)
   real(double) :: tempStorage(:,:)
  real :: ndiffusion(:)
  integer :: nIndex
  integer :: subcell, i !, j, k
  integer :: nIndexScattered
   integer :: nTempStorage
!  integer, parameter :: nTheta = 10, nPhi = 10
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unpackvalues(child,nIndex,nIndexScattered,distanceGrid,adotPAH,nCrossings, &
                   nDiffusion, meanIntensity, tempStorage, nTempStorage)
                exit
             end if
          end do
       else
          nIndex = nIndex + 1
          thisOctal%distanceGrid(subcell) = distanceGrid(nIndex)
          if (usePAH) thisOctal%adotPAH(subcell) = adotPAH(nIndex)
          thisOctal%meanIntensity(subcell) = meanIntensity(nIndex)
          thisOctal%nCrossings(subcell) = int(nCrossings(nIndex), kind=bigint)
          thisOctal%nDiffusion(subcell) = nDiffusion(nIndex)
          if (nTempStorage > 0) then
             if ((.not.associated(thisOctal%tempStorage)).or.(size(thisOctal%tempStorage,2) < nTempStorage)) then
                if (associated(thisOctal%tempStorage)) deallocate(thisOctal%tempStorage)
                allocate(thisOctal%tempStorage(1:thisOctal%maxChildren, 1:nTempStorage))
                thisOctal%tempStorage = 0.d0
             endif
             thisOctal%tempStorage(subcell,1:nTempStorage) = tempStorage(nIndex,1:nTempStorage)
          endif
!          if (storescattered) then
!             do j = 1, nTheta
!                do k = 1, nPhi
!                   nIndexScattered = nIndexScattered + 1
!                   thisOctal%scatteredIntensity(subcell,j,k) = scatteredIntensity(nIndexScattered)
!                enddo
!             enddo
!          endif
       endif
    enddo
  end subroutine unpackvalues

  recursive subroutine directPhotonSmooth(thisOctal, grid, converged, inheritProps, interpProps)
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    logical, optional :: inheritProps, interpProps
    integer :: subcell, i
    logical :: converged
    real(double) :: r
    type(VECTOR) :: dirVec(6), centre, octVec
    real(double) :: thisFac, neighbourFac
    integer :: neighbourSubcell, j
    dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
    dirVec(2) = VECTOR( 0.d0,+1.d0,  0.d0)
    dirVec(3) = VECTOR(+1.d0, 0.d0,  0.d0)
    dirVec(4) = VECTOR(-1.d0, 0.d0,  0.d0)
    dirVec(5) = VECTOR( 0.d0,-1.d0,  0.d0)
    dirVec(6) = VECTOR( 0.d0, 0.d0, -1.d0)

!    do subcell = 1, thisOctal%maxChildren
    subcell = 1
    do while (subcell < thisOctal%maxChildren)
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call directPhotonSmooth(child, grid, converged, inheritProps, interpProps)
                exit
             end if
          end do
       else

          if (thisOctal%nCrossings(subcell) > 0)  then
             thisFac = dble(thisOctal%nDirectPhotons(subcell)) / dble(thisOctal%ncrossings(subcell))
          else
             thisFac = 0.
          endif

          r = thisOctal%subcellSize/2. + grid%halfSmallestSubcell 
          centre = subcellCentre(thisOctal, subcell)
          do j = 1, 6
             octVec = centre + r * dirvec(j)
             if (inOctal(grid%octreeRoot, octVec)) then
                startOctal => thisOctal
                call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                     foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                if (neighbourOctal%nCrossings(neighboursubcell) > 0)  then
                   neighbourFac = dble(neighbourOctal%nDirectPhotons(neighboursubcell)) &
                        / dble(neighbourOctal%ncrossings(neighboursubcell))
                else
                   neighbourFac = 0.
                endif
                if ((min(thisFac, neighbourFac) == 0.0).and.(max(thisFac, neighbourFac) > 0.5)) then
                   if ((thisFac >= 0.).and.(thisFac < neighbourFac)) then
                      call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                           inherit=inheritProps, interp=interpProps)
                      write(*,*) "split",thisFac,neighbourFac
                      converged = .false.
                      subcell = 0
                      exit
                   endif
                endif
             endif
          enddo
       endif
       subcell = subcell + 1
    end do

  end subroutine directPhotonSmooth


  recursive subroutine  calcContinuumEmissivityLucy(grid, thisOctal, nlambda, lamArray)
    type(GRIDTYPE) :: grid
    integer :: nLambda
    real :: lamArray(:)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calcContinuumEmissivityLucy(grid, child, nlambda, lamArray)
                exit
             end if
          end do
       else
          thisOctal%etaCont(subcell) = 1.d-40
          if (thisOctal%temperature(subcell) > 1.d-3) then

             call addDustContinuumLucy(thisOctal, subcell, grid, nlambda, lamArray)
             
          endif

       endif
    enddo
  end subroutine calcContinuumEmissivityLucy

  recursive subroutine  calcContinuumEmissivityLucyMono(grid, thisOctal, lamArray, lambda, iPhotonLambda)
    type(GRIDTYPE) :: grid
    integer :: iPhotonLambda
    real :: lamArray(:), lambda
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calcContinuumEmissivityLucyMono(grid, child, lamArray, lambda, iPhotonLambda)
                exit
             end if
          end do
       else
          thisOctal%etaCont(subcell) = tiny(thisOctal%etaCont)
          if ((thisOctal%temperature(subcell) > 1.d-3).and.(thisOctal%rho(subcell) > 1.d-30)) then

             call addDustContinuumLucyMono(thisOctal, subcell, grid, lambda, iPhotonLambda)
             
          endif

       endif
    enddo
  end subroutine calcContinuumEmissivityLucyMono

  recursive subroutine  setDustTemperatureIfZero(thisOctal)
      use inputs_mod, only : nDustType
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
      integer :: subcell, i, iDust
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setDustTemperatureIfZero(child)
                exit
             end if
          end do
       else
          if ((.not.associated(thisOctal%tempStorage)).or.(size(thisOctal%tempStorage,2) < 2*nDustType)) then
             if (associated(thisOctal%tempStorage)) deallocate(thisOctal%tempStorage)
             allocate(thisOctal%tempStorage(1:thisOctal%maxChildren, 1:2*nDustType))
             do iDust = 1, nDustType
                thisOctal%tempStorage(1:thisOctal%maxChildren,nDustType+iDust) = &
                     dble(thisOctal%temperature(1:thisOctal%maxChildren))
                thisOctal%tempStorage(1:thisOctal%maxChildren,iDust) = 0.d0
             enddo
          endif
          if (thisOctal%tDust(subcell) < 1.d-10) then
             thisOctal%tDust(subcell) = dble(thisOctal%temperature(subcell))
          endif
          do iDust = 1, nDustType
             if (thisOctal%tempStorage(subcell,nDustType+iDust) < 1.d-10) then
                thisOctal%tempStorage(subcell,nDustType+iDust) = dble(thisOctal%temperature(subcell))
             endif
          enddo
       endif
    enddo
  end subroutine setDustTemperatureIfZero

  recursive subroutine  calcContinuumEmissivityLucyMonoAtDustTemp(grid, thisOctal, lamArray, lambda, iPhotonLambda)
    type(GRIDTYPE) :: grid
    integer :: iPhotonLambda
    real :: lamArray(:), lambda
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calcContinuumEmissivityLucyMonoAtDustTemp(grid, child, lamArray, lambda, iPhotonLambda)
                exit
             end if
          end do
       else
          thisOctal%etaCont(subcell) = tiny(thisOctal%etaCont)
          if ((thisOctal%temperature(subcell) > 1.d-3).and.(thisOctal%tdust(subcell) > 1.d-3) & 
             .and.(thisOctal%rho(subcell) > 1.d-30)) then

             call addDustContinuumLucyMonoAtDustTemp(thisOctal, subcell, grid, lambda, iPhotonLambda)
             
          endif

       endif
    enddo
  end subroutine calcContinuumEmissivityLucyMonoAtDustTemp



subroutine addDustContinuumLucy(thisOctal, subcell, grid, nlambda, lamArray)
   use inputs_mod, only : usePAH, nDustType
  type(OCTAL), pointer :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  integer :: nLambda
  real :: lamArray(:)
  integer :: i
   integer :: iDust
  real :: dlam
  real(double), allocatable :: kabsArray(:)


  allocate(kAbsArray(1:nlambda))

  call returnKappa(grid, thisOctal, subcell, kappaAbsArray=kAbsArray)

  thisOctal%etaCont(subcell) = tiny(thisOctal%etaCont(subcell))

  do i = 2, nLambda
     dlam = lamArray(i)-lamArray(i-1)
     if (associated(thisOctal%tempStorage).and.grid%oneKappa) then
        do iDust = 1, nDustType
           thisOctal%etaCont(subcell) = thisOctal%etaCont(subcell) + &
                bLambda(dble(lamArray(i)), thisOctal%tempStorage(subcell,nDustType+iDust)) * &
                thisOctal%dustTypeFraction(subcell,iDust) * grid%oneKappaAbs(iDust,i) * thisOctal%rho(subcell) * &
                1.d-10* dlam * fourPi * 1.d-8
        enddo
     else
        thisOctal%etaCont(subcell) = thisOctal%etaCont(subcell) + &
             bLambda(dble(lamArray(i)), dble(thisOctal%temperature(subcell))) * &
                kAbsArray(i) *1.d-10* dlam * fourPi * 1.d-8 ! conversion from per cm to per A
     endif


     if (usePAH) then
        thisOctal%pahEmissivity(subcell) = PAHemissivityFromAdot(dble(lamArray(i)), &
             thisOctal%adotPAH(subcell), &
             thisOctal%rho(subcell)) & 
             *cSpeed/(lamArray(i)*angstromtocm)**2 *dlam * fourPi * 1.d-8

        thisOctal%etaCont(subcell) = thisOctal%etaCont(subcell) + thisOctal%pahEmissivity(subcell)
     endif
  enddo

  deallocate(kAbsArray)

end subroutine addDustContinuumLucy

!-------------------------------------------------------------------------------

subroutine addDustContinuumLucyMono(thisOctal, subcell, grid,  lambda, iPhotonLambda)
  use pah_mod
   use inputs_mod, only : usePAH, nDustType
  type(OCTAL), pointer :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  integer :: iPhotonLambda
   integer :: iDust
  real ::  lambda
  real(double) :: kappaAbs
  kappaAbs = 0.d0
  thisOctal%etaCont(subcell) = tiny(thisOctal%etaCont(subcell))

  call returnKappa(grid, thisOctal, subcell, lambda=lambda, iLambda=iPhotonLambda, kappaAbs=kappaAbs)

  if (associated(thisOctal%tempStorage).and.grid%oneKappa) then
     thisOctal%etaCont(subcell) = tiny(thisOctal%etaCont(subcell))
     do iDust = 1, nDustType
        thisOctal%etaCont(subcell) = thisOctal%etaCont(subcell) + &
             bLambda(dble(lambda), thisOctal%tempStorage(subcell,nDustType+iDust)) * &
             thisOctal%dustTypeFraction(subcell,iDust) * grid%oneKappaAbs(iDust,iPhotonLambda) * thisOctal%rho(subcell) * &
             1.d-10 * fourPi * 1.d-8
     enddo
  else
     thisOctal%etaCont(subcell) =  bLambda(dble(lambda), dble(thisOctal%temperature(subcell))) * &
                kappaAbs * 1.d-10 * fourPi * 1.d-8 ! conversion from per cm to per A
  endif

  if (usePAH) then
     thisOctal%pahEmissivity(subcell) =  PAHemissivityFromAdot(dble(lambda), &
          thisOctal%adotPAH(subcell), &
          thisOctal%rho(subcell)) &
          *cSpeed/(lambda*angstromtocm)**2 * fourPi * 1.d-8


     thisOctal%etaCont(subcell) = thisOctal%etaCont(subcell) + thisOctal%pahEmissivity(subcell)
  endif

  if (.not.thisOctal%inFlow(subcell)) thisOctal%etaCont(subcell) = 0.d0

end subroutine addDustContinuumLucyMono

!-------------------------------------------------------------------------------

subroutine addDustContinuumLucyMonoAtDustTemp(thisOctal, subcell, grid,  lambda, iPhotonLambda)
   use inputs_mod, only : tminGlobal, decoupleGasDustTemperature, usePAH, nDustType
  type(OCTAL), pointer :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  integer :: iPhotonLambda
   integer :: iDust
  real ::  lambda
  real(double) :: kappaAbs
  logical, save :: firstTime = .true.

  kappaAbs = 0.d0
  thisOctal%etaCont(subcell) = tiny(thisOctal%etaCont(subcell))

  if (.not.decoupleGasDustTemperature) then
     if (thisOctal%tDust(subcell) < tMinglobal) then
        if (firstTime.and.writeoutput) write(*,*) "Looks like tdust not set up, setting tdust to temperature"
        thisOctal%tDust(subcell) = dble(thisOctal%temperature(subcell))
        firstTime = .false.
     endif
  endif

  if (associated(thisOctal%tempStorage)) then
     do iDust = 1, nDustType
        if (thisOctal%tempStorage(subcell,nDustType+iDust) < dble(tMinglobal)) then
           thisOctal%tempStorage(subcell,nDustType+iDust) = dble(thisOctal%temperature(subcell))
        endif
     enddo
  endif


  call returnKappa(grid, thisOctal, subcell, lambda=lambda, iLambda=iPhotonLambda, kappaAbs=kappaAbs)
  if (associated(thisOctal%tempStorage).and.grid%oneKappa) then
     thisOctal%etaCont(subcell) = tiny(thisOctal%etaCont(subcell))
     do iDust = 1, nDustType
        thisOctal%etaCont(subcell) = thisOctal%etaCont(subcell) + &
             bLambda(dble(lambda), thisOctal%tempStorage(subcell,nDustType+iDust)) * &
             thisOctal%dustTypeFraction(subcell,iDust) * grid%oneKappaAbs(iDust,iPhotonLambda) * thisOctal%rho(subcell) * &
             1.d-10 * fourPi * 1.d-8
     enddo
  else
     thisOctal%etaCont(subcell) =  bLambda(dble(lambda), thisOctal%tdust(subcell)) * &
                kappaAbs * 1.d-10 * fourPi * 1.d-8 ! conversion from per cm to per A
  endif
  if (usePAH) then
     thisOctal%pahEmissivity(subcell) =  PAHemissivityFromAdot(dble(lambda), &
          thisOctal%adotPAH(subcell), &
          thisOctal%rho(subcell)) &
          *cSpeed/(lambda*angstromtocm)**2 * fourPi * 1.d-8


     thisOctal%etaCont(subcell) = thisOctal%etaCont(subcell) + thisOctal%pahEmissivity(subcell)
  endif


  if (.not.thisOctal%inFlow(subcell)) thisOctal%etaCont(subcell) = 0.d0

end subroutine addDustContinuumLucyMonoAtDustTemp

!-------------------------------------------------------------------------------

subroutine setBiasOnTau(grid, iLambda)
    use inputs_mod, only : cylindrical, amr3d, amr1d, smallestCellSize
    use amr_mod, only: tauAlongPath, getOctalArray
#ifdef MPI
    use mpi_global_mod,  only : myRankGlobal, nThreadsGlobal
    use mpi
#endif
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    real(double) :: tau, thisTau
  type(VECTOR) :: rVec, direction
  integer :: i
    integer :: subcell
    integer :: iLambda
    integer                     :: nOctal        ! number of octals in grid
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: iOctal
    integer :: iOctal_beg, iOctal_end
    real(double) :: kappaSca, kappaAbs, kappaExt
    type(VECTOR) :: arrayVec(6)
     integer :: nDir
#ifdef MPI
! Only declared in MPI case
     real(double), allocatable :: eArray(:), tArray(:)
     integer :: nVoxels, ierr
     integer :: nBias
     integer :: np, n_rmdr, m
#endif


     kappaAbs = 0.d0; kappasca = 0.d0; thisTau = 0.d0

     call writeInfo("Computing bias on tau...",TRIVIAL)

    allocate(octalArray(grid%nOctals))
    nOctal = 0
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)
    if (nOctal /= grid%nOctals) then
       write(*,*) "Screw up in get octal array", nOctal,grid%nOctals
       stop
    endif

    
    ! default loop indices
    ioctal_beg = 1
    ioctal_end = nOctal

#ifdef MPI
    
            ! Set the range of index for octal loop used later.     
            np = nThreadsGlobal
            n_rmdr = MOD(SIZE(octalArray),np)
            m = SIZE(octalArray)/np
            
            if (myRankGlobal .lt. n_rmdr ) then
               ioctal_beg = (m+1)*myRankGlobal + 1
               ioctal_end = ioctal_beg + m
            else
               ioctal_beg = m*myRankGlobal + 1 + n_rmdr
               ioctal_end = ioctal_beg + m - 1
            end if
            
#endif


!$OMP PARALLEL DEFAULT (NONE) &
!$OMP PRIVATE (iOctal, subcell,  kappaExt, kappaAbs, KappaSca, tau,  thisOctal, direction, thisTau, ndir, arrayvec, rvec) &
!$OMP SHARED (iOctal_beg, iOctal_end, octalArray, grid, cylindrical, ilambda, amr3d, amr1d, smallestCellSize)
     call returnKappa(grid, grid%OctreeRoot, 1, reset_kappa=.true.)

     if (amr3d) then
        nDir = 6
     else
        nDir = 4
     endif

     if (amr1d) ndir = 2
     if (nDir == 4) then
        arrayVec(1) = VECTOR(1.d0, 1.d-10, 1.d-10)
        arrayVec(2) = VECTOR(-1.d0, 1.d-10, 1.d-10)
        arrayVec(3) = VECTOR(1.d-10, 1.d-10, 1.d0)
        arrayVec(4) = VECTOR(1.d-10, 1.d-10,-1.d0)
     endif

     if (nDir == 2) then
        arrayVec(1) = VECTOR(1.d0, 0.d0, 0.d0)
        arrayVec(2) = VECTOR(-1.d0, 0.d0, 0.d0)
     endif


!$OMP DO SCHEDULE (DYNAMIC, 1)
    do iOctal =  iOctal_beg, iOctal_end
       thisOctal => octalArray(iOctal)%content

       do subcell = 1, thisOctal%maxChildren

          if (.not.thisOctal%hasChild(subcell)) then

             thisOctal%biasCont3d(subcell) = 1.d0
!             cycle

             rVec = subcellCentre(thisOctal, subcell)
             if (thisOctal%threed) then
                rVec = rVec + 0.01d0*smallestCellSize*randomUnitVector()
             endif
              
             call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, kappaSca=kappaSca, kappaAbs=kappaAbs)
             kappaExt = kappaAbs + kappaSca
             if (thisOctal%subcellSize*kappaExt < 0.1d0) then
                thisOctal%biasCont3D(subcell) = 1.d0
             else

                tau = 1.d30
                if (amr3d) then
                   if (cylindrical) then
                      ndir = 6
                      arrayVec(1) = VECTOR(0.d0, 0.d0, 1.d0)
                      arrayVec(2) = VECTOR(0.d0, 0.d0,-1.d0)
                      arrayVec(3) = VECTOR(rVec%x, rVec%y,0.d0)
                      call normalize(arrayVec(3))
                      arrayVec(4) = (-1.d0)*arrayVec(3)
                      arrayVec(3) = rotateZ(arrayVec(3), 0.1d0*degtorad)
                      arrayVec(4) = rotateZ(arrayVec(4), 0.1d0*degtorad)
                      arrayVec(5) = arrayVec(3).cross.arrayVec(1)
                      call normalize(arrayVec(5))
                      arrayVec(6) = (-1.d0)*arrayVec(5)
                   else
                      ndir = 6
                      arrayVec(1) = VECTOR(0.d0, 0.d0, 1.d0)
                      arrayVec(2) = VECTOR(0.d0, 0.d0,-1.d0)
                      arrayVec(3) = VECTOR(1.d0, 0.d0, 0.d0)
                      arrayVec(4) = VECTOR(-1.d0, 0.d0,0.d0)
                      arrayVec(5) = VECTOR(0.d0, 1.d0, 0.d0)
                      arrayVec(6) = VECTOR(0.d0,-1.d0, 0.d0)
                   endif
                endif


                do i = 1, ndir 
                  direction = arrayVec(i)
                   call tauAlongPath(ilambda, grid, rVec, direction, thistau, 10.d0, startOctal=thisOctal, startSubcell=subcell)
                   tau = min(tau, thisTau)
                enddo
!                if (tau < 5.) then
!                   thisOctal%biasCont3D(subcell) = 1.d0
!                else
                thisOctal%biasCont3D(subcell) = max(exp(-tau),1.d-6)
                
!                endif
                 

             endif

            
          endif

       enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

#ifdef MPI
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 


     call countVoxels(grid%octreeRoot,nOctal,nVoxels)
     allocate(eArray(1:nVoxels))
     allocate(tArray(1:nVoxels))
     eArray = 0.d0
     call packBias(octalArray, nBias, eArray, ioctal_beg, ioctal_end)
     call MPI_ALLREDUCE(eArray,tArray,nBias,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
     eArray = tArray
     call unpackBias(octalArray, nBias, eArray)
     deallocate(eArray, tArray)
#endif

    deallocate(octalArray)
!    call writeVtkFile(grid, "bias.vtk", &
!            valueTypeString=(/"bias"/))
     call writeInfo("Done.",TRIVIAL)

  end subroutine setBiasOnTau

subroutine setFixedTemperatureOnTau(grid, iLambda)
    use inputs_mod, only : cylindrical, amr3d, amr1d, smallestCellSize
    use amr_mod, only: tauAlongPath, getOctalArray
#ifdef MPI
    use mpi_global_mod,  only : myRankGlobal, nThreadsGlobal
    use mpi
#endif
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    real(double) :: tau, thisTau
  type(VECTOR) :: rVec, direction
  integer :: i
    integer :: subcell
    integer :: iLambda
    integer                     :: nOctal        ! number of octals in grid
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: iOctal
    integer :: iOctal_beg, iOctal_end
    real(double) :: kappaSca, kappaAbs, kappaExt, r
    type(VECTOR) :: arrayVec(6)
     integer :: nDir
#ifdef MPI
! Only declared in MPI case
     logical, allocatable :: eArray(:), tArray(:)
     integer :: nVoxels, ierr
     integer :: nFixed
     integer :: np, n_rmdr, m
#endif


     kappaAbs = 0.d0; kappasca = 0.d0; thisTau = 0.d0

     call writeInfo("Fixing temperature on tau...",TRIVIAL)

    allocate(octalArray(grid%nOctals))
    nOctal = 0
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)
    if (nOctal /= grid%nOctals) then
       write(*,*) "Screw up in get octal array", nOctal,grid%nOctals
       stop
    endif

    
    ! default loop indices
    ioctal_beg = 1
    ioctal_end = nOctal

#ifdef MPI
    
            ! Set the range of index for octal loop used later.     
            np = nThreadsGlobal
            n_rmdr = MOD(SIZE(octalArray),np)
            m = SIZE(octalArray)/np
            
            if (myRankGlobal .lt. n_rmdr ) then
               ioctal_beg = (m+1)*myRankGlobal + 1
               ioctal_end = ioctal_beg + m
            else
               ioctal_beg = m*myRankGlobal + 1 + n_rmdr
               ioctal_end = ioctal_beg + m - 1
            end if
            
#endif


!$OMP PARALLEL DEFAULT (NONE) &
!$OMP PRIVATE (iOctal, subcell,  kappaExt, kappaAbs, KappaSca, tau,  thisOctal, direction, thisTau, ndir, arrayvec, rvec, r) &
!$OMP SHARED (iOctal_beg, iOctal_end, octalArray, grid, cylindrical, ilambda, amr3d, amr1d, smallestCellSize)
     call returnKappa(grid, grid%OctreeRoot, 1, reset_kappa=.true.)

     if (amr3d) then
        nDir = 6
     else
        nDir = 4
     endif

     if (amr1d) ndir = 2
     if (nDir == 4) then
        arrayVec(1) = VECTOR(1.d0, 1.d-10, 1.d-10)
        arrayVec(2) = VECTOR(-1.d0, 1.d-10, 1.d-10)
        arrayVec(3) = VECTOR(1.d-10, 1.d-10, 1.d0)
        arrayVec(4) = VECTOR(1.d-10, 1.d-10,-1.d0)
     endif

     if (nDir == 2) then
        arrayVec(1) = VECTOR(1.d0, 0.d0, 0.d0)
        arrayVec(2) = VECTOR(-1.d0, 0.d0, 0.d0)
     endif


!$OMP DO SCHEDULE (DYNAMIC, 1)
    do iOctal =  iOctal_beg, iOctal_end
       thisOctal => octalArray(iOctal)%content

       do subcell = 1, thisOctal%maxChildren

          if (.not.thisOctal%hasChild(subcell)) then

             rVec = subcellCentre(thisOctal, subcell)
             if (thisOctal%threed) then
                rVec = rVec + 0.01d0*smallestCellSize*randomUnitVector()
             endif
              
             call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, kappaSca=kappaSca, kappaAbs=kappaAbs)
             kappaExt = kappaAbs + kappaSca
             if (thisOctal%subcellSize*kappaExt < 0.1d0) then
                thisOctal%biasCont3D(subcell) = 1.d0
             else

                tau = 1.d30
                if (amr3d) then
                   if (cylindrical) then
                      ndir = 6
                      arrayVec(1) = VECTOR(0.d0, 0.d0, 1.d0)
                      arrayVec(2) = VECTOR(0.d0, 0.d0,-1.d0)
                      arrayVec(3) = VECTOR(rVec%x, rVec%y,0.d0)
                      call normalize(arrayVec(3))
                      arrayVec(4) = (-1.d0)*arrayVec(3)
                      arrayVec(3) = rotateZ(arrayVec(3), 0.1d0*degtorad)
                      arrayVec(4) = rotateZ(arrayVec(4), 0.1d0*degtorad)
                      arrayVec(5) = arrayVec(3).cross.arrayVec(1)
                      call normalize(arrayVec(5))
                      arrayVec(6) = (-1.d0)*arrayVec(5)
                   else
                      ndir = 6
                      arrayVec(1) = VECTOR(0.d0, 0.d0, 1.d0)
                      arrayVec(2) = VECTOR(0.d0, 0.d0,-1.d0)
                      arrayVec(3) = VECTOR(1.d0, 0.d0, 0.d0)
                      arrayVec(4) = VECTOR(-1.d0, 0.d0,0.d0)
                      arrayVec(5) = VECTOR(0.d0, 1.d0, 0.d0)
                      arrayVec(6) = VECTOR(0.d0,-1.d0, 0.d0)
                   endif
                endif


                do i = 1, ndir 
                  direction = arrayVec(i)
                   call tauAlongPath(ilambda, grid, rVec, direction, thistau, 10.d0, startOctal=thisOctal, startSubcell=subcell)
                   tau = min(tau, thisTau)
                enddo

                rVec = subcellCentre(thisOctal, subcell)
                r = sqrt(rVec%x**2 + rVec%y**2)

                if ((tau > 5.).and.(r > 2.*autocm/1.d10).and.(r < 12.d0*autocm)) thisOctal%fixedTemperature(subcell) = .true.
            
             endif
          endif
       enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

#ifdef MPI
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 


     call countVoxels(grid%octreeRoot,nOctal,nVoxels)
     allocate(eArray(1:nVoxels))
     allocate(tArray(1:nVoxels))
     eArray = .false.
     call packFixed(octalArray, nFixed, eArray, ioctal_beg, ioctal_end)
     call MPI_ALLREDUCE(eArray,tArray,nFixed,MPI_LOGICAL,&
         MPI_LOR,MPI_COMM_WORLD,ierr)
     eArray = tArray
     call unpackFixed(octalArray, nFixed, eArray)
     deallocate(eArray, tArray)
#endif

    deallocate(octalArray)
!    call writeVtkFile(grid, "bias.vtk", &
!            valueTypeString=(/"bias"/))
     call writeInfo("Done.",TRIVIAL)

  end subroutine setFixedTemperatureOnTau

  subroutine packBias(octalArray, nBias, eArray, iOctal_beg, iOctal_end)
    type(OCTALWRAPPER) :: octalArray(:)
    integer :: iOctal_beg, iOctal_end
    integer :: nBias
    real(double) :: eArray(:)
    integer :: iOctal, iSubcell
    type(OCTAL), pointer :: thisOctal

       !
       ! Update the bias values of grid computed by all processors.
       !
    nBias = 0
    do iOctal = 1, SIZE(octalArray)
       
       thisOctal => octalArray(iOctal)%content
          
       do iSubcell = 1, thisOctal%maxChildren
                
          if (.not.thisOctal%hasChild(iSubcell)) then
             nBias = nBias + 1
             if ((ioctal >= ioctal_beg).and.(iOctal<=ioctal_end)) then
                eArray(nBias) = octalArray(iOctal)%content%Biascont3d(iSubcell)
             else 
                eArray(nBias) = 0.d0
             endif
          endif
          
       end do
    end do
  end subroutine packBias

  subroutine unpackBias(octalArray, nBias, eArray)
    type(OCTALWRAPPER) :: octalArray(:)
    integer :: nBias
    real(double) :: eArray(:)
    integer :: iOctal, iSubcell
    type(OCTAL), pointer :: thisOctal

    !
    ! Update the bias values of grid computed by all processors.
    !
    nBias = 0
    do iOctal = 1, SIZE(octalArray)

       thisOctal => octalArray(iOctal)%content
          
       do iSubcell = 1, thisOctal%maxChildren
          
          if (.not.thisOctal%hasChild(iSubcell)) then
             nBias = nBias + 1
             octalArray(iOctal)%content%biasCont3d(iSubcell) = eArray(nBias)
          endif
       end do
    end do
  end subroutine unpackBias

  subroutine packFixed(octalArray, nFixed, eArray, iOctal_beg, iOctal_end)
    type(OCTALWRAPPER) :: octalArray(:)
    integer :: iOctal_beg, iOctal_end
    integer :: nFixed
    logical :: eArray(:)
    integer :: iOctal, iSubcell
    type(OCTAL), pointer :: thisOctal

       !
       ! Update the bias values of grid computed by all processors.
       !
    nFixed = 0
    do iOctal = 1, SIZE(octalArray)
       
       thisOctal => octalArray(iOctal)%content
          
       do iSubcell = 1, thisOctal%maxChildren
                
          if (.not.thisOctal%hasChild(iSubcell)) then
             nFixed = nFixed + 1
             if ((ioctal >= ioctal_beg).and.(iOctal<=ioctal_end)) then
                eArray(nFixed) = octalArray(iOctal)%content%fixedTemperature(iSubcell)
             else 
                eArray(nFixed) = .false.
             endif
          endif
          
       end do
    end do
  end subroutine packFixed

  subroutine unpackFixed(octalArray, nFixed, eArray)
    type(OCTALWRAPPER) :: octalArray(:)
    integer :: nFixed
    logical :: eArray(:)
    integer :: iOctal, iSubcell
    type(OCTAL), pointer :: thisOctal

    !
    ! Update the bias values of grid computed by all processors.
    !
    nFixed = 0
    do iOctal = 1, SIZE(octalArray)

       thisOctal => octalArray(iOctal)%content
          
       do iSubcell = 1, thisOctal%maxChildren
          
          if (.not.thisOctal%hasChild(iSubcell)) then
             nFixed = nFixed + 1
             octalArray(iOctal)%content%fixedTemperature(iSubcell) = eArray(nFixed)
          endif
       end do
    end do
  end subroutine unpackFixed

  recursive subroutine allocateMemoryForLucy(thisOctal)
    use octal_mod, only: allocateattribute
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call allocateMemoryForLucy(child)
                exit
             end if
          end do
       else
          call allocateAttribute(thisOctal%diffusionApprox, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%nDiffusion, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%eDens, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%diffusionCoeff, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%oldeDens, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%nDirectPhotons, thisOctal%maxChildren)

          call allocateAttribute(thisOctal%underSampled, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%oldTemperature, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%kappaRoss, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%distanceGrid, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%adotPAH, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%adot, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%nCrossings, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%nTot, thisOctal%maxChildren)


       endif
    enddo
  end subroutine allocateMemoryForLucy

  subroutine addToScatteredIntensity(position, thisOctal, subcell, uHat, tVal)

!    use intensity_storage_mod

    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    type(VECTOR) :: uHat, thisVec, position
    real(double) :: tval, ang
    integer :: ipix

    ang = atan2(position%y, position%x)

    thisVec = uHat
    
    if (thisOctal%twoD) then
       thisVec = rotateZ(uHat, -ang)
    endif

    ipix = 1!returnIpix(thisVec)
    thisOctal%scatteredIntensity(subcell,1,ipix) = thisOctal%scatteredIntensity(subcell,1, ipix) + tval*1.d10
  end subroutine addToScatteredIntensity

  subroutine addMeanIntensity(thisOctal, subcell, tVal)

    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: tVal
!$OMP ATOMIC
    thisOctal%meanIntensity(subcell) = thisOctal%meanIntensity(subcell) + tVal
    
  end subroutine addMeanIntensity


  recursive subroutine calcIntensityFromGrid(thisOctal, epsOverDt, dnu)
!    use intensity_storage_mod
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real(double) :: epsOverDt, dnu
    real(double) :: v
    integer :: subcell, i
    real(double), allocatable :: iNuDOmega(:,:)
!    integer :: npix

!    npix = returnHealpixPixelNumber()
!    allocate(iNuDOmega(1:1,1:nPix))




    Do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calcIntensityFromGrid(child, epsOverDt, dnu)
                exit
             end if
          end do
       else

          V = cellVolume(thisOctal, subcell)*1.d30
!          InuDomega(:,:) = (epsOverDt / V) * thisOctal%scatteredIntensity(subcell, :,:) / dnu

!          thisOctal%scatteredIntensity(subcell, :, :) = inuDomega(:,:) / (fourPi/dble(nPix))

          write(*,*) "scat ",thisOctal%scatteredIntensity(subcell,1,:)

       endif
       
    enddo

    deallocate(iNuDOmega)
  end subroutine calcIntensityFromGrid

  recursive subroutine unrefineThinCells(thisOctal, grid, ilambda, nUnrefine, converged)
    use inputs_mod, only : minDepthAMR
    use amr_mod, only: deleteChild
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: ilambda
    real(double) :: kappaAbs, kappaSca, tau
    integer :: subcell, i, j
    logical :: unrefine, converged
    integer :: nc
    integer :: oldNChildren
    integer :: nUnrefine
    integer :: nChild

    

    kappaAbs = 0.d0; kappaSca = 0.d0
    nc = 0

    if ( thisOctal%nChildren > 0) then
       do i = 1, thisOctal%nChildren
          child => thisOctal%child(i)
          oldNChildren = thisOctal%nChildren
          call unrefineThinCells(child, grid, ilambda, nUnrefine, converged)
          if (thisOctal%NChildren /= oldNChildren) return
       end do
    else
       unrefine = .true.
       do subcell = 1, thisOctal%maxChildren
          call returnKappa(grid, thisOctal, subcell, ilambda, kappaAbs=kappaAbs,kappaSca=kappaSca)
          tau = thisOctal%subcellSize*(kappaAbs+kappaSca)
          if (tau > 1.e-5) then
             unrefine = .false.
             exit
          endif
       enddo
       if (.not.associated(thisOctal%parent).and.writeoutput) then
          write(*,*) "Parent of thisoctal not assocaited ",thisOctal%nDepth
       endif
       if (associated(thisOctal%parent).and.unrefine) then       
          if (thisOctal%nDepth > minDepthAMR) then
             if ((thisOctal%nChildren == 0).and.unrefine.and.converged) then
                nChild = 0
                do j = 1, thisOctal%parent%nChildren
                   if (thisOctal%parent%indexChild(j) == thisOctal%parentSubcell) then
                      nChild = j
                      exit
                   endif
                enddo
                if (nChild == 0) then
                   write(*,*) "fatal bug"
                   write(*,*) thisOctal%parent%nChildren,thisOctal%parent%indexchild
                   stop
                endif
                   
                call deleteChild(thisOctal%parent, thisOctal%parentSubcell, adjustParent = .true., &
                     grid = grid, adjustGridInfo = .true.)
                converged = .false.
                nUnrefine = nUnrefine + 1
             endif
          endif
       endif
    endif

  end subroutine unrefineThinCells

  recursive subroutine unrefineBack(thisOctal, grid, beta, height, rSub, nUnrefine, converged)
    use inputs_mod, only : rOuter, rInner, minDepthAMR, heightSplitFac, maxDepthAMR, smoothinneredge, rSublimation
    use amr_mod, only: deleteChild
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: i
    logical :: unrefine, converged
    logical :: split
    real(double) :: cellSize, r, hr, beta, height, rSub, thisHeightSplitFac
    integer :: oldNChildren
    integer :: nUnrefine

    type(VECTOR) :: cellCentre
    

    unrefine = .true.

    if ( thisOctal%nChildren > 0) then
       do i = 1, thisOctal%nChildren
          child => thisOctal%child(i)
          oldNChildren = thisOctal%nChildren
          call unrefineBack(child, grid, beta, height, rSub, nUnrefine, converged)
          if (thisOctal%NChildren /= oldNChildren) return
       end do
    else
       split = .true.
       if (.not.associated(thisOctal%parent).and.(myrankGlobal==1)) then
          write(*,*) "Parent of thisoctal not assocaited ",thisOctal%nDepth
       endif
       if (associated(thisOctal%parent)) then
          cellSize = thisOctal%parent%subcellSize 
          cellCentre = subcellCentre(thisOctal%parent,thisOctal%parentSubcell)
          r = sqrt(cellcentre%x**2 + cellcentre%y**2)
          hr = height * (r / (100.d0*autocm/1.d10))**beta
          thisHeightSplitfac = heightSplitFac
          if (r < rSublimation) thisHeightSplitFac = 1.
          split = .false.
          

          if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > thisheightSplitFac)) split = .true.
          
          if (.not.split) then
             if ((abs(cellcentre%z)/hr > 2.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.
          endif
          
          if (.not.split) then
             if (.not.smoothInnerEdge) then
                if (((r-cellsize/2.d0) < rinner).and. ((r+cellsize/2.d0) > rInner) .and. &
                     (thisOctal%nDepth < maxdepthamr) .and. (abs(cellCentre%z/hr) < 3.d0) ) split=.true.
             endif
          endif

          if (.not.split) then
             if ((abs(cellcentre%z)/hr > 2.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.
          endif
          if (.not.split) then
             if (((r-cellsize/2.d0) < rSub).and. ((r+cellsize/2.d0) > rSub) .and. &
               (thisOctal%nDepth < maxdepthamr) .and. (abs(cellCentre%z/hr) < 3.d0) ) split=.true.
          endif

          if (.not.split) then
             if (((r-cellsize/2.d0) < rOuter).and. ((r+cellsize/2.d0) > rOuter) .and. &
                  (thisOctal%subcellSize/rOuter > 0.01) .and. (abs(cellCentre%z/hr) < 7.d0) ) split=.true.
          endif


          if (.not.split) then
             if ((r+cellsize/2.d0) < grid%rinner*1.) split = .false.
          endif

          if (.not.split) then
             if ((r-cellsize/2.d0) > grid%router*1.) split = .false.
          endif
       
          if ((.not.split).and.(thisOctal%nDepth > minDepthAMR)) then
             nUnrefine = nUnrefine + thisOctal%parent%nChildren
             call deleteChild(thisOctal%parent, thisOctal%parentSubcell, adjustParent = .true., &
               grid = grid, adjustGridInfo = .true.)
             converged = .false.
          endif
       endif
    endif
  end subroutine unrefineBack


  subroutine putTau(grid, wavelength)
    use amr_mod, only: getxValues, getzValues
    use utils_mod, only: stripSimilarValues
    type(GRIDTYPE) :: grid
    real :: wavelength
    integer :: nx, nz
    real(double), allocatable :: xAxis(:),zAxis(:)
    integer :: iLambda, i
    integer :: nOctals,nVoxels
    call countVoxels(grid%octreeRoot, nOctals, nVoxels)
    call locate(grid%lamArray, grid%nLambda, wavelength, ilambda)
    allocate(xAxis(1:nVoxels))
    nx = 0
    call getxValues(grid%octreeRoot,nx,xAxis)
    call stripSimilarValues(xAxis,nx,1.d-5*grid%halfSmallestSubcell)
    xAxis(1:nx) = xAxis(1:nx) + 1.d-5*grid%halfSmallestSubcell
    do i = 1, nx
       call integrateDownWards(grid, xAxis(i), iLambda)
       call integrateUpWards(grid, xAxis(i), iLambda)
    end do
    deallocate(xAxis)

    allocate(zAxis(1:nVoxels))
    nz = 0
    call getzValues(grid%octreeRoot,nz,zAxis)
    call stripSimilarValues(zAxis,nz,1.d-5*grid%halfSmallestSubcell)
    zAxis(1:nz) = zAxis(1:nz) + 1.d-5*grid%halfSmallestSubcell
    do i = 1, nz
       call integrateRightwards(grid, zAxis(i), iLambda)
    end do
    deallocate(zAxis)

  end subroutine putTau

  subroutine integrateDownwards(grid, x, ilambda)
    type(GRIDTYPE) :: grid
    real(double) :: x, tau
    integer :: iLambda
    real(double) ::  i0, snu, dtau, jnu
    type(VECTOR) :: position, viewVec
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: currentDistance, tVal
    real(double) :: kappaAbs, kappaSca
    real(double) :: kappaAbsDust, kappaScaDust
    viewVec = VECTOR(0.d0, 0.d0, -1.d0)

    position = VECTOR(x, 0.d0, &
         grid%octreeRoot%subcellSize-0.01d0*grid%halfSmallestSubcell)
    thisOctal => grid%octreeRoot
    subcell = 1
    i0 = 0.d0
    tau = 0.d0
    currentDistance  = 0.d0
    do while(inOctal(grid%octreeRoot, position).and.(position%z > 0.))
       call findSubcelllocal(position, thisOctal, subcell)
       call distanceToCellBoundary(grid, position, viewVec, tVal, sOctal=thisOctal)
       call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, kappaAbs = kappaAbs, kappaSca = kappaSca, &
            kappaAbsDust = kappaAbsDust, kappaScaDust=kappaScaDust)
       dTau = (kappaAbsDust + kappaScaDust) * tVal
       
       jnu = kappaAbs * bnu(cspeed/(grid%lamArray(ilambda)*angstromTocm), dble(thisOctal%temperature(subcell)))
       
       if (kappaAbs .ne. 0.d0) then
          snu = jnu/kappaAbs
          i0 = i0 +  exp(-tau) * (1.d0-exp(-dtau))*snu
       else
          snu = tiny(snu)
          i0 = i0 + tiny(i0)
       endif
       
       tau = tau + dtau
       position = position + (tval+1.d-3*grid%halfSmallestSubcell) * viewVec
       thisOctal%etaLine(subcell) =  max(tau,1.d-30)
    end do
  end subroutine integrateDownwards

  subroutine integrateRightwards(grid, z, ilambda)
    type(GRIDTYPE) :: grid
    real(double) :: z, tau
    integer :: iLambda
    real(double) ::  i0, snu, dtau, jnu
    type(VECTOR) :: position, viewVec
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: currentDistance, tVal
    real(double) :: kappaAbs, kappaSca
    real(double) :: kappaAbsDust, kappaScaDust
    viewVec = VECTOR(1.d0, 0.d0, 0.d0)

    position = VECTOR(0.01d0*grid%halfSmallestSubcell, 0.d0, z)
    thisOctal => grid%octreeRoot
    subcell = 1
    i0 = 0.d0
    tau = 0.d0
    currentDistance  = 0.d0
    do while(inOctal(grid%octreeRoot, position).and.(position%x > 0.))
       call findSubcelllocal(position, thisOctal, subcell)
       call distanceToCellBoundary(grid, position, viewVec, tVal, sOctal=thisOctal)
       call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, kappaAbs = kappaAbs, kappaSca = kappaSca, &
            kappaAbsDust = kappaAbsDust, kappaScaDust=kappaScaDust)
       dTau = (kappaAbsDust + kappaScaDust) * tVal
       
       jnu = kappaAbs * bnu(cspeed/(grid%lamArray(ilambda)*angstromTocm), dble(thisOctal%temperature(subcell)))
       
       if (kappaAbs .ne. 0.d0) then
          snu = jnu/kappaAbs
          i0 = i0 +  exp(-tau) * (1.d0-exp(-dtau))*snu
       else
          snu = tiny(snu)
          i0 = i0 + tiny(i0)
       endif
       
       tau = tau + dtau
       position = position + (tval+1.d-3*grid%halfSmallestSubcell) * viewVec
       thisOctal%etaLine(subcell) =  max(tau,thisOctal%etaLine(subcell))
    end do
  end subroutine integrateRightwards


  subroutine integrateUpwards(grid, x, ilambda)
    type(GRIDTYPE) :: grid
    real(double) :: x, tau
    integer :: iLambda
    real(double) ::  i0, snu, dtau, jnu
    type(VECTOR) :: position, viewVec
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: currentDistance, tVal
    real(double) :: kappaAbs, kappaSca
    real(double) :: kappaAbsDust, kappaScaDust
    viewVec = VECTOR(0.d0, 0.d0, 1.d0)

    position = VECTOR(x, 0.d0, -grid%octreeRoot%subcellSize+0.01d0*grid%halfSmallestSubcell)
    thisOctal => grid%octreeRoot
    subcell = 1
    i0 = 0.d0
    tau = 0.d0
    currentDistance  = 0.d0
    do while(inOctal(grid%octreeRoot, position).and.(position%z < 0.))
       call findSubcelllocal(position, thisOctal, subcell)
       call distanceToCellBoundary(grid, position, viewVec, tVal, sOctal=thisOctal)
       call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, kappaAbs = kappaAbs, kappaSca = kappaSca, &
            kappaAbsDust = kappaAbsDust, kappaScaDust=kappaScaDust)
       dTau = (kappaAbsDust + kappaScaDust) * tVal
       
       jnu = kappaAbs * bnu(cspeed/(grid%lamArray(ilambda)*angstromTocm), dble(thisOctal%temperature(subcell)))
       
       if (kappaAbs .ne. 0.d0) then
          snu = jnu/kappaAbs
          i0 = i0 +  exp(-tau) * (1.d0-exp(-dtau))*snu
       else
          snu = tiny(snu)
          i0 = i0 + tiny(i0)
       endif
       
       tau = tau + dtau
       position = position + (tval+1.d-3*grid%halfSmallestSubcell) * viewVec
       thisOctal%etaLine(subcell) = max(tau,1.d-30)
    end do
  end subroutine integrateUpwards

  recursive subroutine  refineDiscGrid(thisOctal, grid, beta, height, rSub, gridconverged, inheritProps, interpProps)
    use inputs_mod, only : rOuter, heightsplitfac, maxDepthAMR, rInner, maxMemoryAvailable
    use memory_mod, only : globalMemoryFootprint
    logical :: gridConverged
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    type(VECTOR) :: cellCentre
    logical, optional :: inheritProps, interpProps
    integer :: subcell, i
    real(double) :: rSub, cellSize, hr, r, beta, height
    logical, save :: firstTimeMem = .true.
    logical :: split
    character(len=80) :: message

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call refineDiscGrid(child, grid, beta, height, rSub, gridconverged, inheritProps, interpProps)
                exit
             end if
          end do
       else

          split = .false.


          cellSize = thisOctal%subcellSize 
          cellCentre = subcellCentre(thisOctal,subCell)
          r = sqrt(cellcentre%x**2 + cellcentre%y**2)
          hr = height * (r / (100.d0*autocm/1.d10))**beta


          if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > heightSplitFac)) split = .true.
          
          if ((abs(cellcentre%z)/hr > 2.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.
          
          if (((r-cellsize/2.d0) < rSub).and. ((r+cellsize/2.d0) > rSub) .and. &
               (thisOctal%nDepth < maxdepthamr) .and. (abs(cellCentre%z/hr) < 3.d0) ) split=.true.
          
          if (((r-cellsize/2.d0) < rInner).and. ((r+cellsize/2.d0) > rInner) .and. &
               (thisOctal%nDepth < maxdepthamr) .and. (abs(cellCentre%z/hr) < 3.d0) ) split=.true.


          if (((r-cellsize/2.d0) < rOuter).and. ((r+cellsize/2.d0) > rOuter) .and. &
               (thisOctal%subcellSize/rOuter > 0.01) .and. (abs(cellCentre%z/hr) < 7.d0) ) split=.true.

          if ((r+cellsize/2.d0) < grid%rinner*1.) split = .false.
          if ((r-cellsize/2.d0) > grid%router*1.) split = .false.

          if (globalMemoryFootprint > maxMemoryAvailable) then
             split = .false.
             if (firstTimeMem) then
                write(message,'(a)') "Maxmimum memory exceeded for grid :"//humanReadableMemory(globalMemoryFootprint)
                call writeWarning(message)
                firstTimeMem = .false.
             endif
          endif


          if (split) then
             call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                  inherit=inheritProps, interp=interpProps)
             gridconverged = .false.
             return
          endif


       endif
    enddo

  end subroutine refineDiscGrid


recursive subroutine quickSublimateLucy(thisOctal, minLevel)
  use inputs_mod, only : grainFrac, nDustType, tsub
  type(octal), pointer   :: thisOctal
  real(double), optional :: minLevel
  real(double) :: thisMinLevel
  type(octal), pointer  :: child 
  ! Where dust is present set dustTypeFraction to this value. 
  integer :: subcell, i, idust
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call quickSublimateLucy(child, minLevel)
                exit
             end if
          end do
       else

          if (present(minLevel)) then
             thisMinLevel = minLevel
          else
             thisMinLevel = 1.d-20
          endif
          do idust = 1, nDustType
             if (associated(thisOctal%tempStorage)) then
                if (thisOctal%tempStorage(subcell,nDustType+idust) > tsub(idust)) then
                   thisOctal%dustTypeFraction(subcell,idust) = thisMinLevel
                else
                   thisOctal%dustTypeFraction(subcell,iDust) = grainFrac(iDust)
                endif
             else if (thisOctal%temperature(subcell) > tsub(idust)) then
                thisOctal%dustTypeFraction(subcell,idust) = thisMinLevel
             else
                thisOctal%dustTypeFraction(subcell,iDust) = grainFrac(iDust)
             endif
          enddo
       endif
    enddo
  end subroutine quickSublimateLucy

end module lucy_new_mod

