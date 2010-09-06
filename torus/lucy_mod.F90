module lucy_mod

  use constants_mod
  use messages_mod
  use vector_mod
  use amr_mod, only: addNewChild, inOctal, distanceToCellBoundary, returnKappa, amrGridValues, &
       countvoxels, findsubcelllocal, findsubcelltd
  use utils_mod, only: hunt, locate, solvequaddble
  use gridtype_mod, only: GRIDTYPE
  use octal_mod, only: OCTAL, octalWrapper, subcellCentre, cellVolume
  use spectrum_mod, only: bnu, blambda, getwavelength
  use timing, only: tune
  use vtk_mod, only: writeVtkFile
  use mpi_global_mod, only: myRankGlobal
  use memory_mod
  implicit none

#ifdef USEMKL
   include 'mkl_vml.fi'
#endif

contains


  subroutine lucyRadiativeEquilibriumAMR(grid, miePhase, nDustType, nMuMie, nLambda, lamArray, &
       source, nSource, nLucy, massEnvelope,  percent_undersampled_min, finalPass)
    use input_variables, only : variableDustSublimation, iterlucy, rCore, scatteredLightWavelength, solveVerticalHydro
    use input_variables, only : smoothFactor, lambdasmooth, taudiff, forceLucyConv, multiLucyFiles
    use input_variables, only : object, maxMemoryAvailable
    use source_mod, only: SOURCETYPE, randomSource, getPhotonPositionDirection
    use phasematrix_mod, only: PHASEMATRIX, newDirectionMie
    use diffusion_mod, only: solvearbitrarydiffusionzones, defineDiffusionOnRosseland, defineDiffusionOnUndersampled, randomwalk
    use amr_mod, only: myScaleSmooth, myTauSmooth, findtotalmass, scaledensityamr
    use dust_mod, only: filldustuniform, stripdustaway, sublimatedust, sublimatedustwr104
    use random_mod
    use gas_opacity_mod, only: atomhydrogenRayXsection
    use gridio_mod, only: writeAMRgrid
#ifdef MPI
    use mpi_global_mod, only: myRankGlobal, nThreadsGlobal
    use parallel_mod, only: mpiBlockHandout, mpiGetBlock
#endif
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    type(GRIDTYPE) :: grid
    integer :: nSource
    integer :: iSmoothLam
    logical :: gridConverged
    integer :: nDustType
    logical, optional :: finalPass
    type(SOURCETYPE) :: source(:), thisSource
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
    real(double) :: kappaScadb, kappaAbsdb
    integer(double) :: nDiffusion
    integer(double) :: nDiffusion_sub
    integer(double) :: nScat, nAbs
    integer(bigInt) :: nMonte, iMonte
    integer :: thisPhotonAbs
    real(oct) :: thisFreq
    real(oct) :: albedo
    logical :: escaped
    real(oct) :: t1, hrecip_kt
    real(oct) :: thisLam,wavelength
    integer :: iSource
    integer :: iLam
    integer(double) :: nInf
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
    integer :: iIter_grand, iMultiplier
    real(oct)::  dT_mean_new, dT_mean_old ! [Kelvins]
    logical :: converged
    real:: percent_undersampled
    integer(double) :: nAbs_sub, nScat_sub, nInf_sub  ! For OpenMP
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
    real(double) :: fac1(nLambda), fac1dnu(nlambda)
    integer :: nVoxels, nOctals

    integer :: icritupper, icritlower
    real(double) :: logt
    real(double) :: logNu1, logNuN, dlogNu, scaleNu
    real(double) :: loglam1, loglamN, scalelam
    real(double) :: logNucritUpper, logNucritLower
    real(double) :: this_bnu(nlambda), fac2(nlambda), hNuOverkT(nlambda)
    real(double) :: subRadius
    real :: lamSmoothArray(5)
    logical :: thisIsFinalPass
    integer(bigInt) :: totMem
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

    ! FOR MPI IMPLEMENTATION=======================================================
    writeoutput = .false.
    if (myRankIsZero) then
       writeoutput = .true.

       print *, ' '
       print *, 'Lucy radiative equilibrium routine computed by ', nThreadsGlobal, ' processors.'
       print *, ' '
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
       if (myrankglobal == 1) write(105,*) lognu1, logt/sigmaE
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


    if (grid%geometry.eq."wr104") then
       totalMass = 0.
       call findTotalMass(grid%octreeRoot, totalMass)
       scaleFac = massEnvelope / totalMass
       if (writeoutput) write(*,'(a,1pe12.5)') "Density scale factor: ",scaleFac
       call scaleDensityAMR(grid%octreeRoot, scaleFac)
    endif



    !    call setupFreqProb(temperature, freq, dnu, nFreq, ProbDistPlanck)

    call writeInfo("Computing lucy radiative equilibrium in AMR...",TRIVIAL)

    lCore = 0.d0
    do i = 1, nSource
       lCore = lCore + source(i)%luminosity
    enddo

    write(message,'(a,1pe12.5)') "Total souce luminosity (lsol): ",lCore/lSol
    call writeInfo(message, TRIVIAL)

    !    if ((grid%geometry == "shakara").or.(grid%geometry == "circumbin")) then
    !       if ((nDustType ==2).and.(aMax(1) <  aMax(2))) then
    !          call writeInfo("Filling dust with large dust in midplane", FORINFO)
    !          call fillDustShakara(grid, grid%octreeRoot)
    !       else
    !          call writeInfo("Filling disc with uniform dust fractions", FORINFO)
    !          call fillDustUniform(grid, grid%octreeRoot)
    !       endif
    !    endif
    !
    !    if (((grid%geometry == "ppdisk").or.(grid%geometry=="warpeddisc")).and.(nDustType > 1)) then
    !       call fillDustUniform(grid, grid%octreeRoot)
    !    endif


    if (object == "ab_aur") then
       call writeInfo("Filling dust with large dust in midplane", FORINFO)
       call fillDustUniform(grid, grid%octreeRoot)
    endif

    if (grid%geometry == "wr104") then
       call fillDustUniform(grid, grid%octreeRoot)
       call stripDustAway(grid%octreeRoot, 1.d-2, 1.d30)
    endif

    if (variableDustSublimation) then
       call stripDustAway(grid%octreeRoot, 1.d-7, 1.d30)
    endif


    nCellsInDiffusion = 0
    call defineDiffusionOnRosseland(grid,grid%octreeRoot, taudiff, ndiff=nCellsInDiffusion)
    write(message,*) "Number of cells in diffusion zone: ", nCellsInDiffusion
    call writeInfo(message,IMPORTANT)

!    call getTauinMidPlane(grid)

    iIter_grand = 0
    converged = .false.
    dT_mean_new = 10.0d0
    dT_mean_old = 10.0d0
    iMultiplier = 1

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

          !          if ( suppressLucySmooth ) thisSmooth = .false.

          if (thisSmooth) then
             call locate(grid%lamArray, nLambda,lambdaSmooth, ismoothlam)
             call writeInfo("Smoothing adaptive grid structure for optical depth...", TRIVIAL)
             do j = iSmoothLam, iSmoothLam !nLambda, 2
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

          call defineDiffusionOnRosseland(grid,grid%octreeRoot,tauDiff,  nDiff=nCellsInDiffusion)


          if (doTuning) call tune(6, "One Lucy Rad Eq Itr")  ! start a stopwatch

          call zeroDistanceGrid(grid%octreeRoot)



          call resetDirectGrid(grid%octreeRoot)

          nInf = 0
          nScat = 0
          nAbs = 0
          nDiffusion = 0
          nKilled = 0

          call countVoxels(grid%OctreeRoot,nOctals,nVoxels)  
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

#ifdef _OPENMP
                 call testRandomOMP()
#endif
                !$OMP PARALLEL DEFAULT(NONE) &
                !$OMP PRIVATE(iMonte, iSource, thisSource, rVec, uHat, rHat) &
                !$OMP PRIVATE(escaped, wavelength, thisFreq, thisLam, iLam, octVec) &
                !$OMP PRIVATE(thisOctal, albedo, r) &
                !$OMP PRIVATE(vec_tmp, uNew, Treal, subcell, probDistJnu) &
                !$OMP PRIVATE(i, j, T1) &
                !$OMP PRIVATE(nAbs_sub, nScat_sub, nInf_sub, nDiffusion_sub, thisPhotonAbs) &
                !$OMP PRIVATE( photonInDiffusionZone, leftHandBoundary, directPhoton) &
                !$OMP PRIVATE(diffusionZoneTemp, kappaAbsdb, sOctal, kappaScadb, kAbsArray) &
                !$OMP PRIVATE(oldUHat) &
                !$OMP PRIVATE(tempOctal, tempSubCell, temp, ok) &
                !$OMP PRIVATE(foundOctal, foundSubcell, hrecip_kt, logt, logNucritUpper, logNucritLower) &
                !$OMP PRIVATE(icritupper, icritlower,  kAbsArray2, hNuOverkT, fac2, this_bnu ) &
                !$OMP PRIVATE(thermalPhoton, scatteredPhoton) &
                !$OMP SHARED(logNu1, fac1dnu, loglam1, scalelam, scalenu)  &
                !$OMP SHARED(grid, nLambda, lamArray,miePhase, nMuMie, nDustType) &
                !$OMP SHARED(imonte_beg, imonte_end, source, nsource) &
                !$OMP SHARED(dnu, nFreq, freq, nMonte) &
                !$OMP SHARED(nAbs, nScat, nInf, nDiffusion, nKilled) 


                 nInf_sub = 0
                 nScat_sub = 0
                 nAbs_sub = 0
                 nDiffusion_sub = 0

                !$OMP DO SCHEDULE(static)
                photonloop: do iMonte = imonte_beg, imonte_end
!                    if (mod(iMonte,imonte_end/10) == 0) write(*,*) "imonte ",imonte
#ifdef MPI
                   !  if (MOD(i,nThreadsGlobal) /= myRankGlobal) cycle photonLoop
#endif
                   

                   thisPhotonAbs = 0
                   call randomSource(source, nSource, iSource)
                   thisSource = source(iSource)
                   call getPhotonPositionDirection(thisSource, rVec, uHat, rHat,grid)
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
                   call getWavelength(thisSource%spectrum, wavelength)
                   thisFreq = cSpeed/(wavelength / 1.e8)
                   thislam = wavelength

                   do while(.not.escaped)

                      ilam = min(floor((log(thislam) - loglam1) * scalelam) + 1, nfreq)
                      ilam = max(ilam, 1)

                      call toNextEventAMR(grid, rVec, uHat, escaped, thisFreq, nLambda, lamArray,  &
                           photonInDiffusionZone, diffusionZoneTemp,  &
                           directPhoton, scatteredPhoton,  &
                           sOctal, foundOctal, foundSubcell, iLamIn=ilam, kappaAbsOut = kappaAbsdb, kappaScaOut = kappaScadb)

                      If (escaped) nInf_sub = nInf_sub + 1

                      if (photonInDiffusionZone) then
                         nDiffusion_sub = nDiffusion_sub + 1
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
                            albedo = kappaScadb / (kappaScadb + kappaAbsdb)
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
                            uNew = newDirectionMie(vec_tmp, real(thisLam), lamArray, nLambda, miePhase, nDustType, nMuMie, &
                                 thisOctal%dustTypeFraction(subcell, 1:nDusttype))

                            nScat_sub = nScat_sub + 1
                            uHat = uNew

                         else

                            nAbs_sub = nAbs_sub + 1
                            thisPhotonAbs = thisPhotonAbs + 1
                            if (thisPhotonAbs > 50000) then
                               nKilled = nKilled + 1
                               cycle photonLoop
                            endif

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
                !$OMP CRITICAL (update)

          nScat = nScat_sub  + nScat   ! sum from each thread for OpenMP
          nInf = nInf_sub    + nInf    ! sum from each thread for OpenMP
          nAbs = nAbs_sub    + nAbs    ! sum from each thread for OpenMP
          nDiffusion = nDiffusion_sub   + nDiffusion    ! sum from each thread for OpenMP
          !$OMP END CRITICAL (update)

          !$OMP BARRIER
          !$OMP END PARALLEL



#ifdef MPI
          ! Summing the value in octal computed by each processors.
          ! This is a recursive function which involves the many 
          ! comunications between processors. It may take a while....
          ! (Maybe faster to pack the values in 1D arrays and distribute.)

          if(doTuning) call tune(6, "  Lucy Loop Update ")  ! start a stopwatch
          if(myRankIsZero) write(*,*) "Calling update_octal_MPI"
          !   call update_octal_MPI(grid%octreeRoot, grid)

          call MPI_BARRIER(MPI_COMM_WORLD, ierr)


          call updateGridMPI(grid)
          if(myRankGlobal == 0) write(*,*) "Done update."

          if(doTuning) call tune(6, "  Lucy Loop Update ")  ! stop a stopwatch

          ! collect some statical info from each node.
          call MPI_ALLGATHER(REAL(nInf), 1, MPI_REAL, buffer_real, 1, MPI_REAL, MPI_COMM_WORLD, ierr)
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          nInf = INT(SUM(buffer_real), double)

          call MPI_ALLGATHER(REAL(nScat), 1, MPI_REAL, buffer_real, 1, MPI_REAL, MPI_COMM_WORLD, ierr)
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          nScat = INT(SUM(buffer_real), double)

          call MPI_ALLGATHER(REAL(nAbs), 1, MPI_REAL, buffer_real, 1, MPI_REAL, MPI_COMM_WORLD, ierr)
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          nAbs = INT(SUM(buffer_real), double)

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


          epsOverDeltaT = lCore / dble(nInf)

          nDT = 0
          nUndersampled = 0


          totalEmission = 0.

          call locate(freq, nFreq, cSpeed/(scatteredLightWavelength*angstromtocm),i)
          call calculateMeanIntensity(grid%octreeRoot, epsOverDeltaT,dnu(i))

          call calculateTemperatureCorrections(.true., grid%octreeRoot, totalEmission, epsOverDeltaT, &
               nFreq, freq, dnu, lamarray, nLambda, grid, nDt, nUndersampled,  &
               dT_sum, dT_min, dT_max, dT_over_T_max)


          if (doTuning) call tune(6, "Gauss-Seidel sweeps")
          call defineDiffusionOnRosseland(grid,grid%octreeRoot, taudiff)


          nCellsInDiffusion = 0
          call defineDiffusionOnUndersampled(grid%octreeroot, nDiff=nCellsInDiffusion)

          if (.not.variableDustSublimation) then
             call solveArbitraryDiffusionZones(grid)
          else
             if (iITer_grand >= 6) then
                call solveArbitraryDiffusionZones(grid)
             endif
          endif

          nCellsInDiffusion = 0
          nUndersampled = 0
          call checkUndersampled(grid%octreeRoot, nUndersampled, nCellsInDiffusion)
          percent_undersampled  = 100.*real(nUndersampled)/real(nVoxels-nCellsInDiffusion)


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



          call defineDiffusionOnRosseland(grid,grid%octreeRoot,tauDiff,  nDiff=nCellsInDiffusion)
          !       call unsetOnDirect(grid%octreeRoot)
          write(message,*) "Number of cells in diffusion zone: ", nCellsInDiffusion
          call writeInfo(message,IMPORTANT)

          if (doTuning) call tune(6, "Gauss-Seidel sweeps")



          if (grid%geometry.eq."wr104") then
             totalMass = 0.
             call findTotalMass(grid%octreeRoot, totalMass)
             scaleFac = massEnvelope / totalMass
             if (writeoutput) write(*,'(a,1pe12.5)') "Density scale factor: ",scaleFac
             call scaleDensityAMR(grid%octreeRoot, scaleFac)
             call sublimateDustWR104(grid%octreeRoot)
          endif



       enddo

       if (((grid%geometry == "ppdisk").or.(grid%geometry=="warpeddisc")).and.(nDustType > 1)) then
          call fillDustUniform(grid, grid%octreeRoot)
       endif

       if (variableDustSublimation) then
          totFrac = 0.
          nFrac = 0
          tauMax = 1.e30


          if (mod(iIter_grand, 2) == 0) then

             if (iIter_grand == 2) tauMax = 0.1
             if (iIter_grand == 4) tauMax = 1.d0
             if (iIter_grand == 6) tauMax = 1.d0
             if (iIter_grand == 8) tauMax = 10.d0
             if (iIter_grand == 10) tauMax = 1.d30

             ! Sublimate the dust and smooth at the photosphere on the last pass
             if (iIter_Grand <= 10) &
                  call sublimateDust(grid, grid%octreeRoot, totFrac, nFrac, tauMax)
             if ((nfrac /= 0).and.(writeoutput)) then
                write(*,*) "Average absolute change in sublimation fraction: ",totFrac/real(nfrac)
             endif

             !             if (iiter_grand == 6) then
             !
             !                call locate(grid%lamArray, nLambda,lambdasmooth,ismoothlam)
             !                if (writeoutput) write(*,*) "Unrefining very optically thin octals..."
             !                gridconverged = .false.
             !                do while(.not.gridconverged)
             !                   gridconverged = .true.
             !                   nUnrefine = 0
             !                   call unrefineThinCells(grid%octreeRoot, grid, ismoothlam, nUnrefine, gridconverged, .false.)
             !                   if (writeoutput) write(*,*) "Unrefined ",nUnrefine, " cells on this pass"
             !                end do
             !                call countVoxels(grid%OctreeRoot,nOctals,nVoxels)  
             !                if (writeoutput) then
             !                   write(*,*) "done."
             !                endif
             !             endif

             if ((iiter_grand >= 6).and.(iiter_grand <= 10)) then
                call locate(grid%lamArray, nLambda,lambdasmooth,ismoothlam)

                call writeInfo("Smoothing adaptive grid structure for optical depth...", TRIVIAL)
                do j = iSmoothLam, nLambda, 20
                   write(message,*) "Smoothing at lam = ",grid%lamArray(j), " angs"
                   call writeInfo(message, TRIVIAL)
                   do
                      gridConverged = .true.
!                      call writeInfo("putting tau")
                      call putTau(grid, grid%lamArray(j))
!                      call writeInfo("done")
                      if (solveVerticalHydro) then
                         call myTauSmooth(grid%octreeRoot, grid, j, gridConverged, &
                              inheritProps = .false., interpProps = .true., photosphereSplit = thisIsFinalPass)
                      else
                         call myTauSmooth(grid%octreeRoot, grid, j, gridConverged, &
                              inheritProps = .false., interpProps = .true.)!, photosphereSplit = .true.)
                      endif

                      if (gridConverged) exit
                   end do
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

             endif


             if (writeoutput) write(*,*) "Global Memory has ",humanReadableMemory(globalMemoryFootprint)
             call findTotalMemory(grid, totMem)
             if (writeoutput) write(*,*) "Full check  has ", humanReadableMemory(totMem)
             if (writeoutput) write(*,*) "Max memory availables ", humanReadableMemory(maxMemoryAvailable)


             nCellsInDiffusion = 0
             call defineDiffusionOnRosseland(grid,grid%octreeRoot, taudiff, ndiff=nCellsInDiffusion)
             write(message,*) "Number of cells in diffusion zone: ", nCellsInDiffusion
             call writeInfo(message,IMPORTANT)



          endif

       endif



       if (grid%geometry == "wr104") then

          totalMass = 0.
          call findTotalMass(grid%octreeRoot, totalMass)
          scaleFac = massEnvelope / totalMass
          if (myRankIsZero) write(*,'(a,1pe12.5)') "Density scale factor: ",scaleFac
          call scaleDensityAMR(grid%octreeRoot, scaleFac)

       endif


       if (abs(totalEmission-oldTotalEmission)/totalEmission < 1.d-2) converged = .true.

       oldTotalEmission = totalEmission

       if (percent_undersampled > percent_undersampled_min) then
          iMultiplier  = iMultiplier * 2
          converged = .false.
       endif


       if ((grid%geometry == "shakara").and.variableDustSublimation.and.(iIter_grand > 2)) then
          call getSublimationRadius(grid, subRadius)
          write(message, '(a, i3, a, f7.3,a )') "End of lucy iteration ",iIter_grand,&
               ": Dust Sublimation radius is: ",(1.d10*subRadius/rSol), " solar radii"
          call writeInfo(message, FORINFO)
          write(message, '(a, i3, a, f7.3,a )') "End of lucy iteration ",iIter_grand,&
               ": Dust Sublimation radius is: ",(subRadius/rCore), " core radii"
          call writeInfo(message, FORINFO)
      endif



       if (iIter_grand < iterlucy) converged = .false.

       if (variableDustSublimation) then
          if (iIter_grand < 9) then
             converged = .false.
          endif
       endif

       if (iIter_grand > 20) then
          write(message,'(a)') "Lucy loop exceeded 20 iterations. Forcing convergence"
          converged = .true.
       endif

       if (variableDustSublimation.and. &
            solveVerticalHydro.and.(.not.thisisfinalPass).and.(iiter_grand == 12)) converged =.true.

       ! forceLucyConv is set in the parameters.dat file if required
       if ( forceLucyConv ) then
          converged = .true.
          if (myRankIsZero) write(*,*) "FORCING CONVERGENCE FOR TESTS!!!!"
       end if

       if (multiLucyFiles) then
          write(tfilename, '(a,i2.2,a)') "lucy",iIter_grand,".vtk"
       else
          tfilename = "lucy.vtk"
       endif


       call writeVtkFile(grid, tfilename, &
            valueTypeString=(/"rho        ", "temperature", "tau        ", "crossings  ", "etacont    " , &
            "dust1      ", "deltaT     ", "etaline    ","fixedtemp  "/))

       !    !
       !    ! Write grid structure to a tmp file.
       !    !


       if (myRankIsZero) call writeAMRgrid("lucy_grid_tmp.dat", .false., grid)



    enddo

    if (variableDustSublimation.and.thisIsFinalPass) then
       call sublimateDust(grid, grid%octreeRoot, totFrac, nFrac, tauMax)
       if ((nfrac /= 0).and.(writeoutput)) then
          write(*,*) "Average absolute change in sublimation fraction: ",totFrac/real(nfrac)
       endif
       call locate(grid%lamArray, nLambda,lambdasmooth,ismoothlam)
       
       call writeInfo("Smoothing adaptive grid structure for optical depth...", TRIVIAL)
       do j = iSmoothLam, nLambda, 2
          write(message,*) "Smoothing at lam = ",grid%lamArray(j), " angs"
          call writeInfo(message, TRIVIAL)
          do
             gridConverged = .true.
             call putTau(grid, grid%lamArray(j))
             call myTauSmooth(grid%octreeRoot, grid, j, gridConverged, &
                  inheritProps = .false., interpProps = .true., photosphereSplit = thisIsFinalPass)
             
             if (gridConverged) exit
          end do
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
       
    endif


    if (grid%geometry == "shakara") then
       call getSublimationRadius(grid, subRadius)
       write(message, '(a, f7.3,a )') "End of lucy loop: Dust Sublimation radius is: ",(1.d10*subRadius/rSol), " solar radii"
       call writeInfo(message, FORINFO)
       write(message, '(a, f7.3,a )') "End of lucy loop: Dust Sublimation radius is: ",(subRadius/rCore), " core radii"
       call writeInfo(message, FORINFO)
    endif

    !    if (storescattered) then 
    !       call locate(freq, nFreq, cSpeed/(1.e4*angstromtocm),i)
    !       call calcIntensityFromGrid(grid%octreeRoot, epsOverDeltaT, dnu(i))
    !       if (writeoutput) call writeVTKfile(grid, "scattered.vtk", valueTypeString = (/"scattered"/))
    !    endif

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
         tau,  ross = .true., tauMax=0.667d0, subRadius=subRadius)

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

  subroutine lucyRadiativeEquilibrium(grid, miePhase, nDustType, nMuMie, nLambda, lamArray, temperature, nLucy)
    use source_mod, only: random_direction_from_sphere
    use phasematrix_mod, only: PHASEMATRIX, newDirectionMie
    use grid_mod, only: getindices

    type(GRIDTYPE) :: grid
    type(VECTOR) :: uHat, uNew
    type(VECTOR) :: rVec
    integer :: nlucy
    integer :: nDustType
    integer,intent(in) :: nLambda, nMuMie
    type(PHASEMATRIX):: miePhase(:,:,:)
    real :: lamArray(:)
    real(oct), allocatable :: distanceGrid(:, :, :)
    real(oct) :: r
    integer, parameter :: nFreq = 100
    integer :: i, j, k
    real(oct) :: freq(nFreq), dnu(nFreq), probDistPlanck(nFreq), probDistJnu(nFreq)
    real(oct) :: temperature
    integer :: nMonte, iMonte, nScat, nAbs
    real(oct) :: thisFreq,  logFreqStart, logFreqEnd
    real(oct) :: albedo
    logical :: escaped
    integer :: i1, i2, i3
    real(oct) :: t1, t2, t3
    real(oct) :: thisLam
    integer :: iLam
    integer :: nInf
    real(oct) :: t
    integer :: nt
    real(oct) :: ang
    integer :: iIter, nIter = 5
    real(oct) :: kabs
    real(oct) ::  kappaP, norm
    real(oct) :: adot, V, epsOverDeltaT
    real(oct) :: newT, deltaT
    real(oct) :: meanDeltaT
    real(double) :: dummy(1)
    real(oct), allocatable :: tempImage(:,:)
    integer :: nDT
    real(oct) :: totalEmission
    type(vector) :: vec_tmp

    allocate(tempImage(1:grid%nx,1:grid%ny))
    dummy = 0.d0; probDistPlanck = 0.

    logFreqStart = log10((cSpeed / (lamArray(nLambda)*1.e-8)))
    logFreqEnd =  log10((cSpeed / (lamArray(1)*1.e-8)))
    write(*,*) logFreqStart, logFreqEnd
    do i = 1, nFreq
       freq(i) = logFreqStart + (dble(i-1)/dble(nFreq-1))*(logFreqEnd-logFreqStart)
       freq(i) = 10.**freq(i)
    enddo
    write(*,*) "Lam",(cSpeed/freq(1))*1.e8,(cSpeed/freq(nFreq))*1.e8

    do i = 2, nFreq-1
       dnu(i) = 0.5*((freq(i+1)+freq(i))-(freq(i)+freq(i-1)))
    enddo
    dnu(1) = freq(2)-freq(1)
    dnu(nFreq) = freq(nFreq)-freq(nFreq-1)



    call setupFreqProb(temperature, freq, dnu, nFreq, ProbDistPlanck)


    write(*,'(a)') "Computing lucy radiative equilibrium..."


    allocate(distanceGrid(1:grid%nx, 1:grid%ny, 1:grid%nz))


    do iIter = 1, nIter
       distanceGrid = 0.d0
       write(*,*) "Iteration",iIter
       nInf = 0
       nScat = 0
       nAbs = 0
       nMonte = nLucy

       do iMonte = 1, nMonte
          escaped = .false.
          call randomNumberGenerator(getDouble=r)
          call locate(probDistPlanck, nFreq, r, j)
          thisFreq = freq(j) + (freq(j+1) - freq(j))* &
               (r - probDistPlanck(j))/(probDistPlanck(j+1)-probDistPlanck(j))

          rVec  = real(grid%rCore, kind=oct) * randomUnitVECTOR()
!          uHat = fromPhotosphereVector(rVec)
          ! -- using a new routine in source_mod.f90 (RK)
          uHat = random_direction_from_sphere(rVec)

          if ( (rVec .dot. uHat) < 0.) uHat = (-1.d0) * uHat

          do while(.not.escaped)

             call toNextEvent(grid, rVec, uHat, escaped, distanceGrid, thisFreq, nLambda, lamArray)
             if (escaped) nInf = nInf + 1

             if (.not. escaped) then

                thisLam = (cSpeed / thisFreq) * 1.e8
                call locate(lamArray, nLambda, real(thisLam), iLam)

                call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
                if (.not.grid%oneKappa) then
                   albedo = grid%kappaSca(i1,i2,i3,iLam)/(grid%kappaSca(i1,i2,i3,iLam)+grid%kappaAbs(i1,i2,i3,iLam))
                else
                   albedo = grid%oneKappaSca(1,iLam) / (grid%oneKappaSca(1,iLam)+grid%oneKappaAbs(1,iLam))
                endif

                call randomNumberGenerator(getDouble=r)
                if (r < albedo) then

                   vec_tmp=uhat 
                   uNew = newDirectionMie(vec_tmp, real(thisLam), lamArray, nLambda, miePhase, nDustType, nMuMie, dummy)
                   nScat = nScat + 1
                   uHat = uNew
                else

                   nAbs = nAbs + 1
                   probDistJnu(1) = 0.d0
                   do i = 2, nFreq

                      thisLam = (cSpeed / freq(i)) * 1.e8
                      call hunt(lamArray, nLambda, real(thisLam), iLam)

                      if ((ilam >=1).and.(ilam <= nlambda)) then
                         if (.not.grid%oneKappa) then
                            kabs = grid%kappaAbs(i1,i2,i3,iLam)
                         else
                            kabs = grid%oneKappaAbs(1,iLam) * grid%rho(i1,i2,i3)
                         endif
                         probDistJnu(i) = probDistJnu(i-1) + &
                              (bnu(freq(i),dble(grid%temperature(i1,i2,i3)))) &
                              * kabs * dnu(i)
                      endif
                   enddo

                   probDistJnu(1:nFreq) = probDistJnu(1:nFreq) / probDistJnu(nFreq)
                   call randomNumberGenerator(getDouble=r)
                   call locate(probDistJnu, nFreq, r, j)
                   if (j == nFreq) j = nFreq -1
                   thisFreq = freq(j) + (freq(j+1) - freq(j))* &
                        (r - probDistJnu(j))/(probDistJnu(j+1)-probDistJnu(j))
                   uHat = randomUnitVector()
                endif

             endif
          enddo
       enddo

       write(*,'(a,f7.2)') "Photons done.",real(ninf)/real(nmonte)
       write(*,'(a,f7.3)') "Mean number of scatters per photon: ",real(nScat)/real(nMonte)
       write(*,'(a,f7.3)') "Mean number of absorbs  per photon: ",real(nAbs)/real(nMonte)

       V = (dble(grid%xAxis(2)-grid%xAxis(1))*dble(grid%yAxis(2)-grid%yAxis(1))*dble(grid%zAxis(2)-grid%zAxis(1)))

       epsOverDeltaT = (dble(grid%lCore)) / dble(nInf)

       meanDeltaT = 0.
       nDT = 0

       totalEmission = 0.
       do i1 = 1, grid%nx
          do i2 = 1, grid%ny
             do i3 = 1, grid%nz
                if (grid%inUse(i1,i2,i3)) then
                   adot = epsoverDeltaT * (1.d0 / v) * distancegrid(i1,i2,i3) / 1.d30
                   kappaP = 0.d0
                   norm = 0.d0
                   do i = 1, nFreq
                      thisLam = (cSpeed / freq(i)) * 1.e8
                      call hunt(lamArray, nLambda, real(thisLam), iLam)
                      if ((iLam >=1) .and. (iLam <= nLambda)) then
                         if (.not.grid%oneKappa) then
                            kabs = grid%kappaAbs(i1,i2,i3,iLam)
                         else
                            kabs = grid%oneKappaAbs(1,iLam) * grid%rho(i1,i2,i3)
                         endif

                         kappaP = kappaP + kabs * &
                              bnu(freq(i),dble(grid%temperature(i1,i2,i3)))  * dnu(i)
                         norm = norm + bnu(freq(i),dble(grid%temperature(i1,i2,i3)))  * dnu(i)
                      endif
                   enddo
                   kappaP = kappaP / norm /1.d10


                   newT = (pi / stefanBoltz) * aDot / (fourPi * kappaP)
                   newT = newT**0.25

                   deltaT = newT - grid%temperature(i1,i2,i3)
                   grid%temperature(i1,i2,i3) = max(1.e-3,grid%temperature(i1,i2,i3) + 0.8 * real(deltaT))
                   nDT = nDT  + 1
                   meanDeltaT = meanDeltaT + deltaT
                   kappaP = 0.d0
                   norm = 0.d0
                   do i = 1, nFreq
                      thisLam = (cSpeed / freq(i)) * 1.e8
                      call hunt(lamArray, nLambda, real(thisLam), iLam)
                      if ((iLam >=1) .and. (iLam <= nLambda)) then
                         if (.not.grid%oneKappa) then
                            kabs = grid%kappaAbs(i1,i2,i3,iLam)
                         else
                            kabs = grid%oneKappaAbs(1,iLam) * grid%rho(i1,i2,i3)
                         endif
                         kappaP = kappaP + kabs * &
                              bnu(freq(i),dble(grid%temperature(i1,i2,i3))) * dnu(i)
                         norm = norm + bnu(freq(i),dble(grid%temperature(i1,i2,i3))) * dnu(i)
                      endif
                   enddo
                   kappaP = kappaP / norm /1.e10
                   grid%etaCont(i1,i2,i3) = fourPi * kappaP * (stefanBoltz/pi) * (grid%temperature(i1,i2,i3)**4)
                   totalEmission = totalEmission + grid%etaCont(i1,i2,i3) * V

                endif
             enddo
          enddo
       enddo

       write(*,*) meanDeltaT / real(nDT)

       write(*,*) "Emissivity of dust / core", totalEmission /grid%lCore * 1.e30

       do i = 1, grid%nx
          do j = 1, grid%ny
             tempImage(i,j) = grid%temperature(i,j,grid%nz/2)
          enddo
       enddo
    enddo
    open(21,file="r.dat",status="unknown",form="formatted")
    do i = 1, 100
       r = grid%xAxis(grid%nx) * real(i-1)/99.
       t = 0
       nT = 0
       do j = 1, 100
          ang = twoPi * real(j-1)/100.
          rVec = VECTOR(r*cos(ang), r*sin(ang),0.)
          call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
          if (grid%inUse(i1,i2,i3)) then
             nT = nT + 1
             t = t + grid%temperature(i1,i2,i3)
          endif
       enddo
       if (nt > 0) then
          t = t / real(nt)
          write(21,*) r / grid%rCore, t
       endif
    enddo
    close(21)

    open(22,file="temps.dat",status="unknown",form="unformatted")
    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             write(22) grid%temperature(i,j,k), grid%etaCont(i,j,k)
          enddo
       enddo
    enddo
    close(22)


  end subroutine lucyRadiativeEquilibrium

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


 subroutine toNextEvent(grid, rVec, uHat,  escaped, distanceGrid, thisFreq, nLambda, lamArray)
   use grid_mod, only: getindices

   type(GRIDTYPE) :: grid
   type(VECTOR) :: rVec, uHat
   integer :: i1, i2, i3
   real(oct) :: t1, t2, t3
   real(oct) :: tval, tau, r
   real :: lamArray(:)
   integer :: nLambda
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
          distanceGrid(i1,i2,i3) = distanceGrid(i1,i2,i3) + (tVal*tau/thisTau) * kabs
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


  recursive subroutine zeroDistanceGrid(thisOctal)
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
          thisOctal%distanceGrid(subcell) = 0.d0
          thisOctal%nCrossings(subcell) = 0
          thisOctal%undersampled(subcell) = .false.
          thisOCtal%nDiffusion(subcell) = 0.
          if (associated(thisOctal%scatteredIntensity)) &
               thisOctal%scatteredIntensity(subcell,:,:) = 0.d0

          if (associated(thisOctal%meanIntensity)) &
               thisOctal%meanIntensity(subcell) = 0.d0
       endif
    enddo
  end subroutine zeroDistanceGrid

  recursive subroutine calculateMeanIntensity(thisOctal, epsOverDt, dnu)
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
                call calculateMeanIntensity(child, epsOVerDt, dnu)
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
  end subroutine calculateMeanIntensity

  recursive subroutine checkUndersampled(thisOctal, nUndersampled, nCellsInDiffusion)
    use input_variables, only : minCrossings
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

    use input_variables, only : minCrossings, TMinGlobal
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
                if (thisOctal%nCrossings(subcell) .ge. 10) then
                   thisOctal%temperature(subcell) = max(TMinGlobal,thisOctal%temperature(subcell) + real(deltaT))
                endif
 
                if (thisOctal%inflow(subcell).and.(thisOctal%nCrossings(subcell) .lt. minCrossings)) then
                   nUnderSampled = nUndersampled + 1
                   thisOctal%undersampled(subcell) = .true.
                endif
             endif


             kappaP = 0.d0
             norm = 0.d0
             if (.not.grid%oneKappa) then
                kabsArray(1:nlambda) = thisOctal%kappaAbs(subcell,1:nlambda)
             else
                call amrGridValues(grid%octreeRoot, subcellCentre(thisOctal,subcell), startOctal=thisOctal, &
                     actualSubcell=subcell, kappaAbsArray=kabsArray, grid=grid)
             endif

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

  subroutine intersectCubeAMR(grid, posVec, direction, tval)
   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(VECTOR), intent(in) :: posVec
   type(VECTOR), intent(in) :: direction
   real(oct), intent(out) :: tval
   !
   type(VECTOR) :: norm(6), p3(6)
   type(OCTAL),pointer :: thisOctal
   type(VECTOR) :: subcen, point
   integer :: subcell
   
   real(oct) :: t(6),denom(6)
   integer :: i,j
   logical :: ok, thisOk(6)


   point = posVec
   write(*,*) "in"
   call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
   write(*,*) "out"
   subcen =  subcellCentre(thisOctal,subcell)
   ok = .true.

   norm(1) = VECTOR(1.0d0, 0.d0, 0.0d0)
   norm(2) = VECTOR(0.0d0, 1.0d0, 0.0d0)
   norm(3) = VECTOR(0.0d0, 0.0d0, 1.0d0)
   norm(4) = VECTOR(-1.0d0, 0.0d0, 0.0d0)
   norm(5) = VECTOR(0.0d0, -1.0d0, 0.0d0)
   norm(6) = VECTOR(0.0d0, 0.0d0, -1.0d0)

   p3(1) = VECTOR(subcen%x+thisOctal%subcellsize/2.0d0, subcen%y, subcen%z)
   p3(2) = VECTOR(subcen%x, subcen%y+thisOctal%subcellsize/2.0d0 ,subcen%z)
   p3(3) = VECTOR(subcen%x,subcen%y,subcen%z+thisOctal%subcellsize/2.0d0)
   p3(4) = VECTOR(subcen%x-thisOctal%subcellsize/2.0d0, subcen%y,  subcen%z)
   p3(5) = VECTOR(subcen%x,subcen%y-thisOctal%subcellsize/2.0d0, subcen%z)
   p3(6) = VECTOR(subcen%x,subcen%y,subcen%z-thisOctal%subcellsize/2.0d0)

   thisOk = .true.

   do i = 1, 6

      denom(i) = norm(i) .dot. direction
      if (denom(i) /= 0.0d0) then
         t(i) = (norm(i) .dot. (p3(i)-posVec))/denom(i)
      else
         thisOk(i) = .false.
         t(i) = 0.0d0
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
     write(*,*) "Error: j=0 (no intersection) in lucy_mod::intersectCubeAMR. "
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  tval = minval(t, mask=thisOk)
  tval = max(tval * 1.001d0,dble(thisOctal%subCellSize/1000.))


  if (tval == 0.) then
     write(*,*) posVec
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  if (tval > sqrt(3.)*thisOctal%subcellsize) then
!     write(*,*) "tval too big",tval/(sqrt(3.)*thisOctal%subcellSize)
!     write(*,*) "direction",direction
!     write(*,*) t(1:6)
!     write(*,*) denom(1:6)
  endif


  end subroutine intersectCubeAMR

  subroutine intersectCubeAMR2D(grid, posVec, direction, tval, sOctal)

! this is to find a cell intersection for a 2D AMR grid
! which is essentially a 2D-grid that projects into
! a 3D cylindrical coordinate system


   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(OCTAL), pointer :: sOctal
   type(VECTOR), intent(in) :: posVec
   type(VECTOR), intent(inout) :: direction
   real(oct), intent(out) :: tval
   !
   type(OCTAL),pointer :: thisOctal
   type(VECTOR) :: subcen, point
   integer :: subcell
   real(double) :: compZ,currentZ
   real(double) :: distToZBoundary, distToXboundary
   real(oct) :: r1,r2,d,cosmu,x1,x2,distTor1,distTor2, theta, mu
   logical :: ok
   type(VECTOR) :: xHat, zHAt

   point = posVec
   write(*,*) 
   call amrGridValues(grid%octreeRoot, point, startOctal=sOctal,foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
   subcen =  subcellCentre(thisOctal,subcell)
   ok = .true.

   r1 = subcen%x - thisOctal%subcellSize/2.d0
   r2 = subcen%x + thisOctal%subcellSize/2.d0
   d = sqrt(point%x**2+point%y**2)
   xHat = VECTOR(point%x, point%y,0.d0)
   call normalize(xHat)

   cosmu =((-1.d0)*xHat).dot.direction
   call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
   if (.not.ok) then
      write(*,*) "Quad solver failed in intersectcubeamr2d"
      direction = randomUnitVector()
      x1 = thisoctal%subcellSize/2.d0
      x2 = 0.d0
   endif
   distTor2 = max(x1,x2)

   theta = asin(max(-1.d0,min(1.d0,r1 / d)))
   cosmu = xHat.dot.direction
   mu = acos(max(-1.d0,min(1.d0,cosmu)))
   distTor1 = 1.e30
   if (mu  < theta ) then
      call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
      if (.not.ok) then
         write(*,*) "Quad solver failed in intersectcubeamr2d"
         direction = randomUnitVector()
         x1 = thisoctal%subcellSize/2.d0
         x2 = 0.d0
      endif
      distTor1 = max(x1,x2)
   endif
         
   distToXboundary = min(distTor1, distTor2)


   zHat = VECTOR(0.d0, 0.d0, 1.d0)
   compZ = zHat.dot.direction
   currentZ = point%z

   if (compZ /= 0.d0 ) then
      if (compZ > 0.d0) then
         distToZboundary = (subcen%z + thisOctal%subcellsize/2.d0 - currentZ ) / compZ
      else
         distToZboundary = abs((subcen%z - thisOctal%subcellsize/2.d0 - currentZ ) / compZ)
      endif
   else
      disttoZboundary = 1.e30
   endif

   tVal = min(distToZboundary, distToXboundary) +0.0001d0*grid%halfsmallestsubcell
   if (tVal > 1.e29) then
      write(*,*) tVal,compZ, distToZboundary,disttoxboundary
      write(*,*) "subcen",subcen
      write(*,*) "z",currentZ
   endif
   if (tval < 0.) then
      write(*,*) tVal,compZ, distToZboundary,disttoxboundary
      write(*,*) "subcen",subcen
      write(*,*) "z",currentZ
   endif

  end subroutine intersectCubeAMR2D

  subroutine intersectCubeCylindrical(grid, posVec, direction, tval, sOctal)

! this is to find a cell intersection for a 3D AMR grid
! which is a cylindrical amr mesh


   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(OCTAL), pointer :: sOctal
   type(VECTOR), intent(in) :: posVec
   type(VECTOR), intent(inout) :: direction
   real(oct), intent(out) :: tval
   !
   type(OCTAL),pointer :: thisOctal
   type(VECTOR) :: subcen, point
   integer :: subcell
   real(double) ::  compZ,currentZ
   real(double) :: distToZBoundary, distToXboundary
   real(oct) :: r1,r2,d,cosmu,x1,x2,distTor1,distTor2, theta, mu
   logical :: ok
   type(VECTOR) :: xHat, zHAt

   point = posVec

   call amrGridValues(grid%octreeRoot, point, startOctal=sOctal,foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
   subcen =  subcellCentre(thisOctal,subcell)
   ok = .true.


   r1 = subcen%x - thisOctal%subcellSize/2.d0
   r2 = subcen%x + thisOctal%subcellSize/2.d0
   d = sqrt(point%x**2+point%y**2)
   xHat = VECTOR(point%x, point%y,0.d0)
   call normalize(xHat)

   cosmu =((-1.d0)*xHat).dot.direction
   call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
   if (.not.ok) then
      write(*,*) "Quad solver failed in intersectcubeamr2d"
      direction = randomUnitVector()
      x1 = thisoctal%subcellSize/2.d0
      x2 = 0.d0
   endif
   distTor2 = max(x1,x2)

   theta = asin(max(-1.d0,min(1.d0,r1 / d)))
   cosmu = xHat.dot.direction
   mu = acos(max(-1.d0,min(1.d0,cosmu)))
   distTor1 = 1.e30
   if (mu  < theta ) then
      call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
      if (.not.ok) then
         write(*,*) "Quad solver failed in intersectcubeamr2d"
         direction = randomUnitVector()
         x1 = thisoctal%subcellSize/2.d0
         x2 = 0.d0
      endif
      distTor1 = max(x1,x2)
   endif
         
   distToXboundary = min(distTor1, distTor2)


   zHat = VECTOR(0.d0, 0.d0, 1.d0)
   compZ = zHat.dot.direction
   currentZ = point%z

   if (compZ /= 0.d0 ) then
      if (compZ > 0.d0) then
         distToZboundary = (subcen%z + thisOctal%subcellsize/2.d0 - currentZ ) / compZ
      else
         distToZboundary = abs((subcen%z - thisOctal%subcellsize/2.d0 - currentZ ) / compZ)
      endif
   else
      disttoZboundary = 1.e30
   endif

   tVal = min(distToZboundary, distToXboundary) +0.0001d0*grid%halfsmallestsubcell
   if (tVal > 1.e29) then
      write(*,*) tVal,compZ, distToZboundary,disttoxboundary
      write(*,*) "subcen",subcen
      write(*,*) "z",currentZ
   endif
   if (tval < 0.) then
      write(*,*) tVal,compZ, distToZboundary,disttoxboundary
      write(*,*) "subcen",subcen
      write(*,*) "z",currentZ
   endif

 end subroutine intersectCubeCylindrical



subroutine toNextEventAMR(grid, rVec, uHat,  escaped,  thisFreq, nLambda, lamArray,  &
      photonInDiffusionZone, diffusionZoneTemp,  directPhoton, scatteredPhoton, &
       startOctal, foundOctal, foundSubcell, ilamIn, kappaAbsOut, kappaScaOut)

   use diffusion_mod, only: randomwalk

   type(GRIDTYPE) :: grid
   type(VECTOR) :: rVec,uHat, octVec
   type(OCTAL), pointer :: thisOctal, tempOctal !, sourceOctal
   type(OCTAL),pointer :: oldOctal, sOctal, topOctal
   type(OCTAL),pointer :: foundOctal, endOctal, startOctal
   logical :: scatteredPhoton
   integer :: foundSubcell
   integer :: endSubcell
   integer :: subcell, tempSubcell!, sourceSubcell
   logical :: directPhoton
   real(oct) :: tval, tau, r
   real :: lamArray(:)
   integer :: nLambda
   logical :: stillinGrid, ok
   logical :: escaped
   real(double) :: kappaScaDb, kappaAbsDb
   real(oct) :: thisTau
   real(oct) :: thisFreq
   real(oct) :: thisLam
   integer :: iLam
   logical ::inFlow
   real :: diffusionZoneTemp
   logical :: photonInDiffusionZone
!   real(double) :: prob
   integer :: i
   real(double), parameter :: fudgeFac = 1.d-2
   integer, optional :: ilamIn
   real(double), optional :: kappaScaOut, kappaAbsOut
!   integer, save :: iLamScat
!   logical, save :: firstTime = .true.
    real :: test  
    logical :: test2
    test = lamArray(1)
    test2 = scatteredPhoton
   endSubcell = 0
   stillinGrid = .true.
   escaped = .false.
   ok = .true.
   photonInDiffusionZone = .false.
   kappaAbsDb = 0.d0; kappaScaDb = 0.d0
   topOctal => grid%octreeRoot

!   if (firstTime) then
!      call hunt(lamArray, nLambda, scatteredLightWavelength, iLamScat)
!      firstTime = .false.
!   endif


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
         grid=grid, inFlow=inFlow)

    if(present(kappaAbsOut)) kappaAbsOut = kappaAbsdb
    if(present(kappaScaOut)) kappaScaOut = kappaScadb

    oldOctal => thisOctal
    sOctal => thisOctal

! moved from before call to amrgridvalues - th 16/11/05

    call distanceToCellBoundary(grid, rVec, uHat, tVal, sOctal)

    tval = tval + fudgeFac * grid%halfSmallestSubcell

!    if (grid%octreeRoot%threed) then
!       call intersectCubeAMR(grid, rVec, uHat, tVal)
!    else
!       call intersectCubeAMR2D(grid, rVec, uHat, tVal, sOctal)
!    endif

    if (inFlow) then
       thisTau = dble(tVal) * (kappaAbsdb + kappaScadb)
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

             do;enddo
          endif

       endif


    endif

! check whether the photon has escaped from the grid


! two cases here now. in the 2D case we only update the distance grid in a single plane (y=0, x>=0)

! update the distance grid

!$OMP CRITICAL (changegrid)

          thisOctal%distanceGrid(subcell) = thisOctal%distanceGrid(subcell) &
               + tVal * dble(kappaAbsdb)
          thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1
          if (directPhoton) thisOctal%nDirectPhotons(subcell) = thisOctal%nDirectPhotons(subcell)+1

!          if (scatteredPhoton.and.storeScattered.and.(iLam==iLamScat)) &
!               call addToScatteredIntensity(octVec, thisOctal, subcell, uHat, tVal)

!          if (storeScattered.and.(iLam==iLamScat)) &
!               call addMeanIntensity(octVec, thisOctal, subcell, uHat, tVal*1.d10)


!$OMP END CRITICAL (changegrid)


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

          sOctal => thisOctal
          oldOctal => thisOctal

! calculate the distance to the next cell

          call distanceToCellBoundary(grid, rVec, uHat, tVal, sOctal)
          tval = tval + fudgeFac*grid%halfSmallestSubcell

!          if (grid%octreeRoot%threed) then
!             call intersectCubeAMR(grid, rVec, uHat, tVal)
!          else
!             call intersectCubeAMR2D(grid, rVec, uHat, tVal, thisOctal)
!          endif
          octVec = rVec

! calculate the optical depth to the next cell boundary

          if (inFlow) then
             thisTau = dble(tVal) * (kappaAbsdb + kappaScadb)
          else
             thisTau = 1.0e-28
          end if

          if (tVal == 0.0d0) then
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

!       if (grid%octreeRoot%threed) then
!          call intersectCubeAMR(grid, rVec, uHat, tVal)
!       else
!          call intersectCubeAMR2D(grid, rVec, uHat, tVal)
!       endif

       call amrGridValues(topOctal, octVec, startOctal=oldOctal,iLambda=iLam, &
            foundOctal=thisOctal, foundSubcell=subcell, & 
            kappaAbs=kappaAbsdb,kappaSca=kappaScadb, grid=grid, inFlow=inFlow)

       if(present(kappaAbsOut)) kappaAbsOut = kappaAbsdb
       if(present(kappaScaOut)) kappaScaOut = kappaScadb



       if (thisOctal%diffusionApprox(subcell)) then
          write(*,*) "Photon in diff zone but not escaped"
       endif

       if (.not.inFlow) then
          kappaAbsdb =0.0d0
          if(present(kappaAbsOut)) kappaAbsOut = kappaAbsdb
       endif


! update the distance grid

       if (thisTau > 0.d0) then

!$OMP CRITICAL (changegrid2)

          

             thisOctal%distanceGrid(subcell) = thisOctal%distanceGrid(subcell) &
                  + (dble(tVal)*dble(tau)/thisTau) * dble(kappaAbsdb)
             thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1
          if (directPhoton) thisOctal%nDirectPhotons(subcell) = thisOctal%nDirectPhotons(subcell)+1
!          if (scatteredPhoton.and.storeScattered.and.(iLam==iLamScat)) &
!               call addToScatteredIntensity(octVec, thisOctal, subcell, uHat, dble(tVal)*dble(tau)/thisTau)

!          if (storeScattered.and.(iLam==iLamScat)) &
!               call addMeanIntensity(octVec, thisOctal, subcell, uHat, 1.d10*dble(tVal)*dble(tau)/thisTau)



!$OMP END CRITICAL (changegrid2)


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
    implicit none
    include 'mpif.h'

    type(gridtype) :: grid
    type(octal), pointer  :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    ! data space to store values from all processors
    real  :: buffer_ncrossings(nThreadsGlobal)     
    real(double) :: buffer_distanceGrid(nThreadsGlobal) 
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
          
          
          call MPI_ALLGATHER(REAL(thisOctal%ncrossings(subcell)), 1, MPI_REAL, &
               buffer_ncrossings, 1, MPI_REAL, MPI_COMM_WORLD, ierr)  
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          
          thisOctal%distanceGrid(subcell) = SUM(buffer_distanceGrid)
          thisOctal%nCrossings(subcell) = INT(SUM(buffer_ncrossings))
          
       endif
    enddo
  end subroutine update_octal_MPI


  subroutine updateGridMPI(grid)
!    use input_variables, only : storeScattered
    implicit none
    include 'mpif.h'
    type(gridtype) :: grid
    integer :: nOctals, nVoxels
    real, allocatable :: nCrossings(:)
    real, allocatable :: tempRealArray(:)
    real(double), allocatable :: distanceGrid(:),tempDoubleArray(:), scatteredIntensity(:)
    real, allocatable :: nDiffusion(:)
    integer, parameter :: nTheta = 11, nPhi = 10
    integer :: ierr, nIndex, nIndexScattered

    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
    nOctals = 0
    nVoxels = 0
    call countVoxels(grid%octreeRoot,nOctals,nVoxels)
    allocate(nCrossings(1:nVoxels))
    allocate(nDiffusion(1:nVoxels))
    allocate(distanceGrid(1:nVoxels))

    allocate(scatteredIntensity(1:(nVoxels*nTheta*nPhi)))

    nIndex = 0
    nIndexScattered = 0 
    call packValues(grid%octreeRoot,nIndex,nIndexScattered, &
         distanceGrid,nCrossings,ndiffusion, scatteredIntensity)



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

    tempRealArray = 0.0
    call MPI_ALLREDUCE(nCrossings,tempRealArray,nVoxels,MPI_REAL,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    nCrossings = tempRealArray 

    tempRealArray = 0.0
    call MPI_ALLREDUCE(nDiffusion,tempRealArray,nVoxels,MPI_REAL,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    nDiffusion = tempRealArray 
    
    deallocate(tempRealArray, tempDoubleArray)
     
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
    
    nIndex = 0
    nIndexScattered = 0
    call unpackValues(grid%octreeRoot,nIndex,nIndexScattered, &
         distanceGrid,nCrossings,nDiffusion,scatteredIntensity)

    deallocate(nCrossings, nDiffusion, distanceGrid, scatteredIntensity)

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
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real :: lamArray(:)
    real(double) :: kappap, norm, freq(:), dnu(:), kabsArray(:), thisLam
    integer :: j, iLam, nLambda, nFreq

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
    enddo
  end subroutine calculateEtaCont

  recursive subroutine packvalues(thisOctal,nIndex, nIndexScattered,&
       distanceGrid,nCrossings, nDiffusion, scatteredIntensity)
!    use input_variables, only : storeScattered
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: distanceGrid(:)
  real(double) :: scatteredIntensity(:)
  real :: nCrossings(:)
  real :: nDiffusion(:)
  integer :: nIndex, nIndexScattered
  integer :: subcell, i !, j , k
!  integer, parameter :: nTheta = 10, nPhi = 10
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call packvalues(child,nIndex,nIndexScattered,distanceGrid,nCrossings,nDiffusion, scatteredIntensity)
                exit
             end if
          end do
       else
          nIndex = nIndex + 1
          distanceGrid(nIndex) = thisOctal%distanceGrid(subcell)
          nCrossings(nIndex) = real(thisOctal%nCrossings(subcell))
          nDiffusion(nIndex) = thisOctal%nDiffusion(subcell)
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

  recursive subroutine unpackvalues(thisOctal,nIndex,nIndexScattered,distanceGrid,nCrossings, nDiffusion, scatteredIntensity)
!    use input_variables, only : storeScattered
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: distanceGrid(:), scatteredIntensity(:)
  real :: ncrossings(:)
  real :: ndiffusion(:)
  integer :: nIndex
  integer :: subcell, i !, j, k
  integer :: nIndexScattered
!  integer, parameter :: nTheta = 10, nPhi = 10
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unpackvalues(child,nIndex,nIndexScattered,distanceGrid,nCrossings, &
                     nDiffusion, scatteredIntensity)
                exit
             end if
          end do
       else
          nIndex = nIndex + 1
          thisOctal%distanceGrid(subcell) = distanceGrid(nIndex)
          thisOctal%nCrossings(subcell) = int(nCrossings(nIndex))
          thisOctal%nDiffusion(subcell) = nDiffusion(nIndex)
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
    real :: thisFac, neighbourFac
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
             thisFac = real(thisOctal%nDirectPhotons(subcell)) / real(thisOctal%ncrossings(subcell))
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
                   neighbourFac = real(neighbourOctal%nDirectPhotons(neighboursubcell)) &
                        / real(neighbourOctal%ncrossings(neighboursubcell))
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

  recursive subroutine  calcContinuumEmissivityLucyMono(grid, thisOctal, nlambda, lamArray, lambda, iPhotonLambda)
    type(GRIDTYPE) :: grid
    integer :: nLambda
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
                call calcContinuumEmissivityLucyMono(grid, child, nlambda, lamArray, lambda, iPhotonLambda)
                exit
             end if
          end do
       else
          thisOctal%etaCont(subcell) = 1.d-40
          if (thisOctal%temperature(subcell) > 1.d-3) then

             call addDustContinuumLucyMono(thisOctal, subcell, grid, lambda, iPhotonLambda)
             
          endif

       endif
    enddo
  end subroutine calcContinuumEmissivityLucyMono



subroutine addDustContinuumLucy(thisOctal, subcell, grid, nlambda, lamArray)

  type(OCTAL), pointer :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  integer :: nLambda
  real :: lamArray(:)
  integer :: i
  real :: dlam
  real(double), allocatable :: kabsArray(:)


  allocate(kAbsArray(1:nlambda))

  call returnKappa(grid, thisOctal, subcell, kappaAbsArray=kAbsArray)

  thisOctal%etaCont(subcell) = tiny(thisOctal%etaCont(subcell))

  do i = 2, nLambda
     dlam = lamArray(i)-lamArray(i-1)
     thisOctal%etaCont(subcell) = thisOctal%etaCont(subcell) + &
          bLambda(dble(lamArray(i)), dble(thisOctal%temperature(subcell))) * &
             kAbsArray(i) *1.d-10* dlam * fourPi * 1.d-8 ! conversion from per cm to per A
  enddo

  deallocate(kAbsArray)

end subroutine addDustContinuumLucy

!-------------------------------------------------------------------------------

subroutine addDustContinuumLucyMono(thisOctal, subcell, grid,  lambda, iPhotonLambda)

  type(OCTAL), pointer :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  integer :: iPhotonLambda
  real ::  lambda
  real(double) :: kappaAbs
  kappaAbs = 0.d0
  thisOctal%etaCont(subcell) = tiny(thisOctal%etaCont(subcell))

  call returnKappa(grid, thisOctal, subcell, lambda=lambda, iLambda=iPhotonLambda, kappaAbs=kappaAbs)

  thisOctal%etaCont(subcell) =  bLambda(dble(lambda), dble(thisOctal%temperature(subcell))) * &
             kappaAbs * 1.d-10 * fourPi * 1.d-8 ! conversion from per cm to per A
  if (.not.thisOctal%inFlow(subcell)) thisOctal%etaCont(subcell) = 0.d0

end subroutine addDustContinuumLucyMono

!-------------------------------------------------------------------------------

subroutine setBiasOnTau(grid, iLambda)
    use input_variables, only : cylindrical
    use amr_mod, only: tauAlongPath, getOctalArray
#ifdef MPI
    use mpi_global_mod,  only : myRankGlobal, nThreadsGlobal
    use input_variables, only : blockHandout
    use parallel_mod, only: mpiBlockHandout, mpiGetBlock
    include 'mpif.h'
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
     integer, dimension(:), allocatable :: octalsBelongRank
     logical :: rankComplete
     integer :: tag = 0
     real(double), allocatable :: eArray(:), tArray(:)
     integer :: nVoxels, ierr
     integer :: nBias
#endif
     kappaAbs = 0.d0; kappasca = 0.d0; thisTau = 0.d0

     if (cylindrical) then
        nDir = 6
     else
        nDir = 4
     endif

     if (nDir == 4) then
        arrayVec(1) = VECTOR(1.d0, 1.d-10, 1.d-10)
        arrayVec(2) = VECTOR(-1.d0, 1.d-10, 1.d-10)
        arrayVec(3) = VECTOR(1.d-10, 1.d-10, 1.d0)
        arrayVec(4) = VECTOR(1.d-10, 1.d-10,-1.d0)
     endif

    allocate(octalArray(grid%nOctals))
    nOctal = 0
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)
    if (nOctal /= grid%nOctals) then
       write(*,*) "Screw up in get octal array", nOctal,grid%nOctals
       stop
    endif

#ifdef MPI
    ! FOR MPI IMPLEMENTATION=======================================================
    
    ! we will use an array to store the rank of the process
    !   which will calculate each octal's variables
    allocate(octalsBelongRank(size(octalArray)))
    
    if (myRankGlobal == 0) then
       print *, ' '
       print *, 'Tau bias  computed by ', nThreadsGlobal-1, ' processors.'
       print *, ' '
       call mpiBlockHandout(nThreadsGlobal,octalsBelongRank,blockDivFactor=1,tag=tag,&
                            setDebug=.false.)
    
    endif
    ! ============================================================================
#endif
    
    ! default loop indices
    ioctal_beg = 1
    ioctal_end = nOctal

#ifdef MPI
 if (myRankGlobal /= 0) then
  blockLoop: do     
 call mpiGetBlock(myRankGlobal,iOctal_beg,iOctal_end,rankComplete,tag,setDebug=.false.)
   if (rankComplete) exit blockLoop 
#endif
    do iOctal =  iOctal_beg, iOctal_end
       thisOctal => octalArray(iOctal)%content

       do subcell = 1, thisOctal%maxChildren

          if (.not.thisOctal%hasChild(subcell)) then

             rVec = subcellCentre(thisOctal, subcell)
             if (thisOctal%threed) then
                rVec = rVec + 0.01d0*grid%halfSmallestSubcell*randomUnitVector()
             endif
              
             call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, kappaSca=kappaSca, kappaAbs=kappaAbs)
             kappaExt = kappaAbs + kappaSca
             if (thisOctal%subcellSize*kappaExt < 0.1d0) then
                thisOctal%biasCont3D(subcell) = 1.d0
             else

                tau = 1.d30
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
                endif



                do i = 1, ndir
                   direction = arrayVec(i)
                   call tauAlongPath(ilambda, grid, rVec, direction, thistau, 20.d0, startOctal=thisOctal, startSubcell=subcell)
                   tau = min(tau, thisTau)
                enddo
                if (tau < 5.) then
                   thisOctal%biasCont3D(subcell) = 1.d0
                else
                   thisOctal%biasCont3D(subcell) = 1.e-20
                endif
!                max(exp(-tau),1.e-5)

             endif

            
          endif

       enddo
    enddo

#ifdef MPI
 if (.not.blockHandout) exit blockloop
 end do blockLoop       
 end if ! (myRankGlobal /= 0)


     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     ! have to send out the 'octalsBelongRank' array
     call MPI_BCAST(octalsBelongRank,SIZE(octalsBelongRank),  &
                    MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     call countVoxels(grid%octreeRoot,nOctal,nVoxels)
     allocate(eArray(1:nVoxels))
     allocate(tArray(1:nVoxels))
     eArray = 0.d0
     call packBias(octalArray, nBias, eArray,octalsBelongRank)
     call MPI_ALLREDUCE(eArray,tArray,nBias,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
     eArray = tArray
     call unpackBias(octalArray, nBias, eArray)
     deallocate(eArray, tArray, octalsBelongRank)
#endif

    deallocate(octalArray)
!    call writeVtkFile(grid, "bias.vtk", &
!            valueTypeString=(/"bias"/))

  end subroutine setBiasOnTau

  subroutine packBias(octalArray, nBias, eArray, octalsBelongRank)
    USE mpi_global_mod, ONLY: myRankGlobal
    type(OCTALWRAPPER) :: octalArray(:)
    integer :: octalsBelongRank(:)
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
             if (octalsBelongRank(iOctal) == myRankGlobal) then
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
          call allocateAttribute(thisOctal%nCrossings, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%nTot, thisOctal%maxChildren)


       endif
    enddo
  end subroutine allocateMemoryForLucy

  subroutine addToScatteredIntensity(position, thisOctal, subcell, uHat, tVal)

    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    type(VECTOR) :: uHat, thisVec, position
    real(double) :: tVal, thisTheta, thisPhi, ang
    integer :: nTheta, nPhi
    integer :: iTheta, iPhi

    nTheta = SIZE(thisOctal%scatteredIntensity,2)
    nPhi = SIZE(thisOctal%scatteredIntensity,3)


    ang = atan2(position%y, position%x)

    thisVec = uHat
    
    if (thisOctal%twoD) then
       thisVec = rotateZ(uHat, -ang)
    endif

    thisTheta = acos(thisvec%z)
    thisPhi = atan2(thisVec%y,thisVec%x)

    if (thisPhi < 0.d0) thisPhi = thisPhi + twoPi 
    iTheta = nint((thisTheta / pi) * dble(nTheta-1))+1
    iphi = nint((thisPhi / twoPi) * dble(nPhi-1))+1

    thisOctal%scatteredIntensity(subcell,iTheta,iPhi) = thisOctal%scatteredIntensity(subcell,iTheta, iPhi) + tval*1.d10
!    write(*,*) thisOctal%scatteredIntensity(subcell,itheta,iphi)
  end subroutine addToScatteredIntensity

  subroutine addMeanIntensity(thisOctal, subcell, tVal)

    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: tVal
    thisOctal%meanIntensity(subcell) = thisOctal%meanIntensity(subcell) + tVal
    
  end subroutine addMeanIntensity


  recursive subroutine calcIntensityFromGrid(thisOctal, epsOverDt, dnu)

    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real(double) :: epsOverDt, dnu
    real(double) :: v
    integer :: subcell, i, j, k
    integer :: nTheta, nphi
    real(double), allocatable :: iNuDomega(:,:)
    real(double) :: dTheta, dphi, theta, dOmega




    nTheta = SIZE(thisOctal%scatteredIntensity,2)
    nPhi = SIZE(thisOctal%scatteredIntensity,3)

    allocate(iNudOmega(nTheta, nPhi))

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
          InuDomega(:,:) = (epsOverDt / V) * thisOctal%scatteredIntensity(subcell, :,:) / dnu
          dTheta = pi / dble(nTheta-1)
          dPhi = twoPi / dble(nPhi)
          do j = 1, nTheta
             theta = pi * dble(j-1)/dble(nTheta-1)
             do k = 1, nPhi
                dOmega = dTheta * dphi * sin(theta)
                thisOctal%scatteredIntensity(subcell, j, k) = inuDomega(j,k) / max(dOmega,1.d-10)
             enddo
          enddo
       endif
       
    enddo

    deallocate(iNuDOmega)
  end subroutine calcIntensityFromGrid

  recursive subroutine unrefineThinCells(thisOctal, grid, ilambda, nUnrefine, converged)
    use input_variables, only : minDepthAMR
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
    use input_variables, only : rOuter, rInner, minDepthAMR, heightSplitFac, maxDepthAMR, smoothinneredge
    use amr_mod, only: deleteChild
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: i
    logical :: unrefine, converged
    logical :: split
    real(double) :: cellSize, r, hr, beta, height, rSub
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

          split = .false.
          
          if (r < 10.*grid%rinner) then
             if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > heightSplitFac/2.)) split = .true.
          else
             if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > heightSplitFac)) split = .true.
          endif
          
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
    use amr_mod, only: getxValues
    use utils_mod, only: stripSimilarValues
    type(GRIDTYPE) :: grid
    real :: wavelength
    integer :: nx
    real(double), allocatable :: xAxis(:)
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

    position = VECTOR(x, 0.d0, grid%octreeRoot%subcellSize-0.01d0*grid%halfSmallestSubcell)
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
    use input_variables, only : rOuter, heightsplitfac, maxDepthAMR, rInner, smoothinneredge, maxMemoryAvailable
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


          if (r < 10.*grid%rinner) then
             if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > heightSplitFac/2.)) split = .true.
          else
             if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > heightSplitFac)) split = .true.
          endif
          
          if ((abs(cellcentre%z)/hr > 2.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.
          
          if (.not.smoothInnerEdge) then
             if (((r-cellsize/2.d0) < grid%rinner).and. ((r+cellsize/2.d0) > grid%rInner) .and. &
                  (thisOctal%nDepth < maxdepthamr) .and. (abs(cellCentre%z/hr) < 3.d0) ) split=.true.
          endif
          
          if ((abs(cellcentre%z)/hr > 2.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.

!          if (((r-cellsize/2.d0) < rSub).and. ((r+cellsize/2.d0) > rSub) .and. &
!               (thisOctal%nDepth < maxdepthamr) .and. (abs(cellCentre%z/hr) < 3.d0) ) split=.true.

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


end module lucy_mod

