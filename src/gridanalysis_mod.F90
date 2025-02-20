module gridanalysis_mod
  use kind_mod
  use gridtype_mod
  use amr_utils_mod
  use random_mod
  use vector_mod
  use source_mod
  implicit none

contains

  subroutine analysis(grid)
    use inputs_mod, only : imodel
    use utils_mod, only : findMultifilename
    type(GRIDTYPE) :: grid
    real(double) :: mass, mass14, mass15, mass16, mdisc
    real(double) :: flux
    character(len=80) :: thisFile

    flux = fluxThroughSphericalSurface(grid, 1500.d0*autocm/1.d10)
    flux = flux / msol * 365.25d0*24.d0*3600.d0
    mass = 0.d0
    call findMassWithBounds(grid%octreeRoot, mass, maxRadius=1500.d0*autocm/1.d10)

    mass14 = 0.d0
    mass15 = 0.d0
    mass16 = 0.d0
    mDisc = 0.d0
    call findMassWithBounds(grid%octreeRoot, mass14, minRho=1.d-14)
    call findMassWithBounds(grid%octreeRoot, mass15, minRho=1.d-15)
    call findMassWithBounds(grid%octreeRoot, mass16, minRho=1.d-16)
    call findDiscMass(grid%octreeRoot, globalSourceArray(1)%mass, mDisc, 2000.d0*autocm/1.d10)
    write(*,*) "disc mass ",mdisc/msol
    call averageRadialProfile(grid, 10000.d0*autocm/1.d10)
    call findMultifilename("vel****.dat",iModel,thisfile)
    open(69, file=thisFile, status="unknown", form="formatted")
    call writeVelocityWithBounds(grid%octreeRoot, maxRadius=1500.d0*autocm/1.d10)
    close(69)
    call findMultifilename("kep****.dat",iModel,thisfile)
    call  writeKeplerianButterfly(thisfile,globalSourceArray(1)%mass, 26.d0*autocm, 1500.d0*autocm)
    write(*,'(i4,1p,11e12.3)') iModel, globalSourceArray(1)%mass/msol, globalSourceArray(1)%mdot/msol*356.25d0*24.d0*3600.d0, &
flux, mass/msol, mass14/msol, mass15/msol, mass16/msol, mdisc/msol
    write(48,'(i4,1p,11e12.3)') iModel, globalSourceArray(1)%mass/msol, globalSourceArray(1)%mdot/msol*356.25d0*24.d0*3600.d0, &
flux, mass/msol, mass14/msol, mass15/msol, mass16/msol, mdisc/msol


  end subroutine analysis

  subroutine calculateColumn(columnDensity, cs, inputOctal, inputSubcell, grid)
    type(VECTOR) :: rVec, direction
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, inputOctal
    real(double) :: ds,cs,totweight,columndensity,mu,pressure
    integer :: subcell, inputsubcell


    columnDensity = 0.d0
    cs = 0.d0
    totweight = 0.d0

    direction = VECTOR(0.d0, 0.d0, 1.d0)
    thisOctal => inputOctal
    subcell = inputSubcell
    rVec = subcellCentre(inputOctal, inputSubcell)
    do while(inOctal(grid%octreeRoot, rVec))
       call findSubcellLocal(rVec, thisOctal, subcell)
       call distanceToCellBoundary(grid, rVec, direction, ds, sOctal=thisOctal)
       columnDensity = columnDensity + thisOctal%subcellSize * thisOctal%rho(subcell) * 1.d10
       mu = 1.27
       pressure = thisOctal%temperature(subcell) * kerg * thisOctal%rho(subcell) / (mu * mHydrogen)
       cs = cs + sqrt(pressure / thisOctal%rho(subcell)) * thisOctal%rho(subcell)
       totweight = totweight + thisOctal%rho(subcell)
       rVec = rVec + (ds + 1.d-6)*direction
    end do
    direction = VECTOR(0.d0, 0.d0, -1.d0)
    thisOctal => inputOctal
    subcell = inputSubcell
    rVec = subcellCentre(inputOctal, inputSubcell)
    do while(inOctal(grid%octreeRoot, rVec))
       call findSubcellLocal(rVec, thisOctal, subcell)
       call distanceToCellBoundary(grid, rVec, direction, ds, sOctal=thisOctal)
       columnDensity = columnDensity + thisOctal%subcellSize * thisOctal%rho(subcell) * 1.d10
       mu = 1.27
       pressure = thisOctal%temperature(subcell) * kerg * thisOctal%rho(subcell) / (mu * mHydrogen)
       cs = cs + sqrt(pressure / thisOctal%rho(subcell)) * thisOctal%rho(subcell)
       totweight = totweight + thisOctal%rho(subcell)
       rVec = rVec + (ds + 1.d-6)*direction
    end do
    cs = cs / totweight
  end subroutine calculateColumn

  real(double) function massWithinR(thisR, grid)
    real(double) :: thisR
    type(GRIDTYPE) :: grid

    massWithinR = 0.d0
    call findMassWithBounds(grid%octreeRoot, massWithinR, maxRadius=thisR)
  end function massWithinR

  recursive subroutine calculateToomreQ(thisOctal, grid)
    use ion_mod
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(VECTOR) :: rVec
    real(double) :: r, cs, omegaK, mass, sigma, toomreQ, mu, pressure
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateToomreQ(child, grid)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal, subcell)
          if (abs(abs(rVec%z) - thisOCtal%subcellSize/2.d0) < 1496.) then
             call calculateColumn(sigma, cs, thisOctal, subcell, grid)
             mass = massWithinR(modulus(rVec), grid)
             mass = mass + SUM(globalSourceArray(1:GlobalnSource)%mass)
             r = modulus(rVec)
             omegaK = sqrt(bigG * mass / (r*1d10)**3)

             mu = returnMu(thisOctal, subcell, globalIonArray, nGlobalIon)
             pressure = thisOctal%temperature(subcell) * kerg * thisOctal%rho(subcell) / (mu * mHydrogen)
             cs = sqrt(pressure / thisOctal%rho(subcell))
             toomreQ = (cs * omegaK) / (pi * bigG * sigma)
             
             if (.not.associated(thisOctal%etaCont)) allocate(thisOctal%etaCont(1:thisOctal%maxChildren))
             thisOctal%etaCont(subcell) = toomreQ
             write(*,*) "toomre Q " , toomreQ
          endif
       endif
    enddo
  end subroutine calculateToomreQ

          



  function fluxThroughSphericalSurface(grid, radius) result(totflux)
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: radius,da,flux
    integer :: i, n
    real(double) ::  totflux
    type(VECTOR) :: rVec, inVec
    n = 100000
    da = fourPi*(radius*1.d10)**2/dble(n)
    totflux = 0.d0
    thisOctal => grid%octreeRoot
    do i = 1, n

       rVec = radius*randomUnitVector()

       inVec = (-1.d0)*rVec
       call normalize(inVec)

       call findSubcellLocal(rVec, thisOctal, subcell)
       flux = invec%x * thisOctal%rhou(subcell) + &
            invec%y * thisOctal%rhov(subcell) + &
            invec%z * thisOctal%rhow(subcell)
       totFlux = totFlux + flux * da
    enddo
  end function fluxThroughSphericalSurface

  recursive subroutine findMassWithBounds(thisOctal, totalMass, minRho, maxRho, minRadius, maxRadius, highestRho, totalVolume, &
            ionizedMass, ionizedVolume, totalKE, totalTE)
    use ion_mod
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(VECTOR) :: rVec, velocity
    real(double) :: totalMass
    real(double),optional :: ionizedMass, ionizedVolume, totalKE, totalTE, totalVolume
    real(double),optional :: minRho, maxRho, minRadius, maxRadius
    real(double),optional :: highestRho
    real(double) :: dv,r, mu
    integer :: subcell, i
    logical :: includeThisCell

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findMassWithBounds(child, totalMass, minRho, maxRho, minRadius, maxRadius, highestRho, totalVolume, &
                  ionizedMass, ionizedVolume, totalKE, totalTE)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
             dv = cellVolume(thisOctal, subcell)*1.d30
             rVec = subcellCentre(thisOctal, subcell) 
             r = modulus(rVec)

             includeThisCell = .true.

             if (PRESENT(minRadius)) then
                if (r < minRadius) includeThisCell = .false.
             endif

             if (PRESENT(maxRadius)) then
                if (r > maxRadius) includeThisCell = .false.
             endif

             if (PRESENT(minRho)) then
                if (thisOctal%rho(subcell) < minRho) includeThisCell = .false.
             endif

             if (PRESENT(maxRho)) then
                if (thisOctal%rho(subcell) > maxRho) includeThisCell = .false.
             endif

             
             if (includethiscell) then
                totalMass = totalMass + thisOctal%rho(subcell) * dv
                if (present(highestRho)) highestRho = max(highestRho, thisOctal%rho(subcell))
                if (present(totalVolume)) totalVolume = totalVolume + dv
                if (present(ionizedMass)) then
                   if (thisOctal%ionfrac(subcell, 1) < 0.1d0 .and. .not.thisOctal%undersampled(subcell)) then
                      ionizedMass = ionizedMass + thisOctal%rho(subcell) * dv
                      ionizedVolume = ionizedVolume + dv
                   endif
                endif
                if (present(totalKE)) then
                   mu = returnMu(thisOctal, subcell, globalIonArray, nGlobalIon)
                   totalTE = totalTE + (1.5d0 * thisOctal%temperature(subcell)*kerg/(mu*mHydrogen))
                   velocity = VECTOR(thisOctal%rhou(subcell) / thisOctal%rho(subcell), &
                        thisOctal%rhov(subcell) / thisOctal%rho(subcell), &
                        thisOctal%rhow(subcell) / thisOctal%rho(subcell))
                   totalKE = totalKE + 0.5d0 * (thisOctal%rho(subcell) * dv) * (modulus(velocity))**2

                endif
             endif

          endif
       endif
    enddo
  end subroutine findMassWithBounds


  subroutine averageRadialProfile(grid, maxRadius)
    type(GRIDTYPE) :: grid
    real(double) :: maxRadius
    integer, parameter :: nr = 30
    integer :: i
    real(double), allocatable :: r(:), rho(:), t(:),n(:)
    allocate(r(1:nr), rho(1:nr), t(1:nr), n(1:nr))
    r = 0.d0; rho = 0.d0; t = 0.d0; n = 0.d0
    do i = 1, nr
!       r(i) = log10(0.01d0*maxradius) + (log10(maxRadius)-log10(0.01d0*maxradius))*(dble(i-1)/dble(nr-1))
!       r(i) = 10.d0**r(i)
       r(i) = maxRadius*dble(i-1)/dble(nr-1)
    enddo
    
    call fillRadialMidplaneArray(grid%octreeRoot, nr, r, rho, t, n, maxRadius)
    open(20,file="midplane.dat",status="unknown",form="formatted")
    write(20,*) "radius (au)    density(g cm^-3)    temperature (K)"
    do i = 1, nr
       write(20,'(f7.0,1pe12.2,f8.0,i5)') r(i)*1.d10/autocm, rho(i)/n(i), t(i)/n(i),int(n(i))
    enddo
    close(20)
  end subroutine averageRadialProfile


  recursive subroutine fillradialMidplaneArray(thisOctal, nr, r, rho, t, n ,maxRadius)

    type(OCTAL), pointer :: thisOctal, child
    real(double) :: r(:), rho(:), t(:), n(:), maxRadius
    integer :: subcell, i, j, nr
    type(VECTOR):: cen
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fillradialMidplaneArray(child, nr, r, rho, t, n ,maxRadius)
                exit
             end if
          end do
       else
          cen = subcellCentre(thisOctal, subcell)
          if (sqrt(cen%x**2 + cen%y**2)  < maxRadius+0.5d0*(r(nr)-r(nr-1))) then
             if ( (cen%z + thisOctal%subcellSize/2.d0*1.001d0)* (cen%z - thisOctal%subcellSize/2.d0*1.001d0) < 0.d0) then
                call locate(r,nr,sqrt(cen%x**2+cen%y**2),j)
                if (sqrt(cen%x**2+cen%y**2) > maxRadius) j = nr
                rho(j) = rho(j) + thisOCtal%rho(subcell)
                t(j) = t(j) + thisOctal%temperature(Subcell)
                n(j) = n(j) + 1.d0
             endif
          endif
       endif
    enddo
  end subroutine fillradialMidplaneArray

  recursive subroutine findDiscMass(thisOctal, sourceMass, totalMass,  maxRadius)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(VECTOR) :: rVec, rHat, vHat, vel, vkep
    real(double) :: totalMass
    real(double) :: maxRadius, sourceMass
    real(double) :: dv,r
    integer :: subcell, i
    logical :: includeThisCell

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findDiscMass(child, sourceMass, totalMass, maxRadius)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
             dv = cellVolume(thisOctal, subcell)*1.d30
             rVec = subcellCentre(thisOctal, subcell) 
             r = modulus(rVec)

             includeThisCell = .true.

             if (r > maxRadius) includeThisCell = .false.
             vel = VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell),thisOctal%rhov(subcell)/thisOctal%rho(subcell), 0.d0)
             rHat = rVec
             call normalize(rHat)
             vHat = rHat.cross.zhat
             call normalize(vHat)
             vKep = sqrt(bigG*sourceMass/(r*1.d10))*vHat
             if (modulus(vel-vKep)/modulus(vKep) > 0.2d0) then
                includeThisCell = .false.
             endif
             if (includethiscell) totalMass = totalMass + thisOctal%rho(subcell) * dv
          endif
       endif
    enddo
  end subroutine findDiscMass

  recursive subroutine writeVelocityWithBounds(thisOctal, minRho, maxRho, minRadius, maxRadius)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(VECTOR) :: rVec
    real(double) :: r
    real(double),optional :: minRho, maxRho, minRadius, maxRadius
    real(double) :: dv
    integer :: subcell, i
    logical :: includeThisCell

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call writeVelocityWithBounds(child, minRho, maxRho, minRadius, maxRadius)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
             dv = cellVolume(thisOctal, subcell)*1.d30
             rVec = subcellCentre(thisOctal, subcell) 

             r = modulus(rVec)

             includeThisCell = .false.

             if ((abs(rVec%z)-thisOctal%subcellSize/2.d0)  < 1.d0) then
!                if ((abs(rVec%y)-thisOctal%subcellSize/2.d0)  < 1.d0) then
                   includeThisCell = .true.
!                endif
             endif


             if (PRESENT(minRadius)) then
                if (r < minRadius) includeThisCell = .false.
             endif

             if (PRESENT(maxRadius)) then
                if (r > maxRadius) includeThisCell = .false.
             endif

             if (PRESENT(minRho)) then
                if (thisOctal%rho(subcell) < minRho) includeThisCell = .false.
             endif

             if (PRESENT(maxRho)) then
                if (thisOctal%rho(subcell) > maxRho) includeThisCell = .false.
             endif

             if (includethiscell) write(69,*) rVec%x/1496.d0, &
                  thisOctal%rhov(subcell)/thisOctal%rho(subcell)/1.d5
          endif
       endif
    enddo
  end subroutine writeVelocityWithBounds

  subroutine writeKeplerianButterfly(thisfile,mass, rMin, rMax)
    real(double) :: mass, rMin, rMax, r, v
    character(len=*) :: thisFile
    integer :: i
    open(66, file=thisFile, status="unknown", form="formatted")
    do i = 1, 1000
       r = rMax + (rMin-rMax)*dble(i-1)/999.d0
       v = sqrt(bigG * mass/r)
       write(66,*) -r/autocm,v/1.d5
    enddo
    do i = 1, 1000
       r = rMin + (rMax-rMin)*dble(i-1)/999.d0
       v = -sqrt(bigG * mass/r)
       write(66,*) r/autocm,v/1.d5
    enddo
    v = sqrt(bigG * mass/rMax)
    write(66,*) -rMax/autocm,v/1.d5

    close(66)
  end subroutine writeKeplerianButterfly

#ifdef MPI
  subroutine clusterAnalysis(grid, source, nSource, nLambda, lamArray, miePhase, nMuMie)
#else
  subroutine clusterAnalysis(grid, source, nsource)
#endif
    use inputs_mod, only : splitOverMPI, burstTime, findHabing, findHabingCluster, findHIIradius, clustersinks
    use inputs_mod, only:  clusterRadial, radiusTotals
    use inputs_mod, only : smallestCellSize, calculateLymanFlux, calculateTauFUV, primarySource
    use inputs_mod, only : iModel
    use vtk_mod
!    use utils_mod, only : findMultifilename
    use phasematrix_mod
#ifdef MPI
    use mpi
    use mpi_amr_mod, only : dumpValuesAlongLine
    use hydrodynamics_mod, only : setupevenuparray
    use inputs_mod, only: columnImageDirection, findNUndersampled, calculateGlobalAvgTdust, calculateGlobalAvgTemp 
    use inputs_mod, only: calculateEmissionMeasure, nClusterIonLoops, calculateAvgPressure, plotAvgTemp, plotAvgTdust, writeLums
#ifdef PHOTOION
    use photoionAMR_mod, only : createTemperatureColumnImage, writeLuminosities, calculateAverageTemperature, photoionizationLoopAMR
    use photoionAMR_mod, only : createTdustColumnImage, calculateAverageTdust, createEmissionMeasureImage, calculateAveragePressure
    use photoionAMR_mod, only : createPressureImage
    use mpi_amr_mod, only : createIonizationImage, createColumnDensityImage, createhiimage, findNumberUndersampled
#endif
#endif
#ifdef USECFITSIO
    use image_mod, only : writeFitsColumnDensityImage 
#endif

! Arguments
    type(GRIDTYPE)    :: grid
    type(SOURCETYPE)  :: source(:)
    integer           :: nSource
#ifdef MPI
    integer           :: nLambda    
    integer           :: nMuMie
    type(PHASEMATRIX) :: miePhase(:,:,:)
    real              :: lamArray(:)
#endif

    logical, save :: firstTime=.true.
    logical, save :: firstTimeStars(1:1000,0:1000) = .true.
    type(VECTOR) :: rcom, rVec
    real(double)  :: weightedFluxInRadius, massInRadius, meang0, g0inCell, distance
    integer :: iSource
    real(double) :: nlyA, nlyB, nlyR, pgas, prad!, sigmapgas, sigmaprad
    character(len=80) :: thisFile
    real(double) :: tauAbs, tauSca, tauExt, columnDensity
    real(double) :: total
#ifdef MPI
    real(double) :: n0, t0, trms, nrms, sigma, sigmane, tempDouble
    real(double) :: loopLimitTime
    character(len=20) :: weighting, comp
    integer :: ierr
    real :: iterTime(3)
    integer :: nSampled, nUndersampled, nCells, nPhotoIter
#ifdef USECFITSIO
    real(double), pointer :: image(:,:)=>null(), rmsImage(:,:)=>null()
#endif
#endif
    real(double) :: vmean(7), vms(7), sigmav(7), mion, mneu
    integer :: i, j, imax, jmax
    real(double) :: radii(5), flux, trad0, mmax!, sums(12) 
    real(double), allocatable :: lum(:)
#ifdef MPI
    integer :: evenUpArray(nHydroThreadsGlobal), iterStack(3), optID
#endif
!    type(VECTOR) :: startPoint, endPoint, firststartpoint, thisDir
!    real(double) :: maxRho, totalMass 
!    character(len=80) :: thisFileGrid, thisFileRadius!, rootFilename, fm
!    real(double), save :: radius 
!    real(double) :: mass, highestRho, volume, ionizedMass, ionizedVolume, totalKE, totalTE, massflux

#ifdef PHOTOION
#ifdef MPI
    if (nClusterIonLoops > 0) then
!       nPhotoIter = 1
       nPhotoIter = nClusterIonLoops
       loopLimitTime = 1.d40 
       iterTime = 1.e30
       call setupevenuparray(grid, evenuparray)
       call photoIonizationloopAMR(grid, source, nsource, nLambda, &
            lamArray, nPhotoIter, loopLimitTime, looplimittime, .false.,iterTime,.true., &
            evenuparray, optID, iterStack, miePhase, nMuMie) 
    endif
#endif
#endif

!    write(rootfilename, '(a)') "habingRadial****.dat"
!    call findMultiFilename(rootfilename, iModel, thisFile)
!    startPoint = source(1)%position 
!    endPoint = VECTOR(grid%octreeRoot%subcellSize, 0.d0, 0.d0) 
!    call dumpValuesAlongLine(grid, thisFile, startPoint, endPoint, 1000)



!    if (firsttime) then
!       firstStartPoint = source(1)%position
!    endif
!    ! values in the cell containing the most massive star
!    write(thisFile, '(a)') "firststartpoint.dat"
!    ! write header
!    if (firstTime) then
!       if (writeoutput) then
!          open(69, file=thisFile, status="replace", form="formatted")
!!             write(69, '(a6,a20, 3(1x,a12), 2(1x,a9))') "# dump", "t (s)", "rho (g/cc)", "UV (G0)", "tval (TU)", & 
!             write(69, '(a6,a20, 3(1x,a12), 2(1x,a9))') "# dump", "t (s)", "rho (g/cc)", "UV (G0)", "rtostar(TU)", & 
!                   "Tgas (K)", "Tdust (K)" 
!          close(69)
!       endif
!    endif
    ! write data
!    open(69, file=thisFile, status="old", position="append", form="formatted")
!       call writeCellValuesForSource(69, grid, source(1), firstStartPoint)
!    close(69)

!    ! get total mass and max rho
!    call findMassOverAllThreads(grid, totalMass, maxRho=maxRho)
!    if (writeoutput) then
!       write(thisFile, '(a)') "gridvalues.dat"
!       ! write header
!       if (firstTime) then
!          open(69, file=thisFile, status="replace", form="formatted")
!             write(69, '(a6,a20, 2(1x,a12))') "# dump", "t (s)", "maxrho(g/cc)", "tot M(msol)"
!          close(69)
!       endif
!       ! write data
!       open(69, file=thisFile, status="old", position="append", form="formatted")
!           write(69, '(i6.4,f20.2, 2(1x,es12.5)') grid%idump, time, maxrho, totalMass/msol 
!       close(69)
!    endif

!   temperature weighted-average along los 
#ifdef MPI
    if (plotAvgTemp .or. calculateGlobalAvgTemp) then
       do i=3,7
          if (i==1) then
             write(weighting, '(a)') "emissivity" 
          elseif (i==2) then
             write(weighting, '(a)') "oiiidensity" 
          elseif (i==3) then
             write(weighting, '(a)') "mass" 
          elseif (i==4) then
             !write(weighting, '(a)') "hiidensity" 
             cycle
          elseif (i==5) then
             write(weighting, '(a)') "ne2" 
          elseif (i==6) then
             write(weighting, '(a)') "none" 
          elseif (i==7) then
             write(weighting, '(a)') "ionNone" 
          endif
#ifdef USECFITSIO
          if (plotAvgTemp) then
             ! pixel-by-pixel average
             if (writeoutput) write(*,*) "Making temperature image with weighting ", weighting
             call createTemperatureColumnImage(grid, columnImageDirection, image, trim(weighting))
             write(thisFile, '(i4.4,a,a,a)') grid%idump, "_avgtemp_", trim(weighting), ".fits"
             if (writeoutput) call writeFitsColumnDensityImage(image, thisFile)
          endif
#else
          call writeInfo("FITS not enabled, not writing",FORINFO)
#endif

          if (calculateGlobalAvgTemp) then
             ! global average
             if (writeoutput) write(*,*) "Calculating avg temperature with weighting ", weighting
             sigma = 0.d0
             total = 0.d0
             sigmaNe = 0.d0
             call calculateAverageTemperature(grid%octreeRoot, grid,sigma, total, trim(weighting), .false., 0.d0, sigmaNe, 0.d0)
#ifdef MPI
             call MPI_ALLREDUCE(sigma, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             sigma = tempDouble ! sum(T w dV)
             call MPI_ALLREDUCE(total, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             total = tempDouble ! sum(w dV)
             call MPI_ALLREDUCE(sigmaNe, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             sigmaNe = tempDouble ! sum(ne w dV)
#endif
             t0 = sigma/total ! mean temperature
             n0 = sigmaNe/total ! mean electron density

             ! rms 
             sigma = 0.d0
             total = 0.d0
             sigmaNe = 0.d0
             call calculateAverageTemperature(grid%octreeRoot, grid, sigma, total, trim(weighting), .true., t0, sigmaNe, n0)
#ifdef MPI
             call MPI_ALLREDUCE(sigma, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             sigma = tempDouble ! sum( (T-T0)^2 w dV)
             call MPI_ALLREDUCE(total, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             total = tempDouble ! sum( w dV)
             call MPI_ALLREDUCE(sigmaNe, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             sigmaNe = tempDouble ! sum( (ne-n0)^2 w dV)
#endif
             ! rms = sqrt(<(T-T0)^2>)
             trms = sqrt(sigma/total)
             nrms = sqrt(sigmaNe/total)
#ifdef MPI
             if (myrankglobal == 1) then
#endif
                write(thisFile, '(a)') "avgtempGlobalRms.dat"
                open(55, file=thisFile, status="unknown", position="append")
                if (firstTime) then
                   write(55,'(a5,a20,a13,4(a12,1x))') "#dump", "t(s)", "weight", "T0(K)", "Trms(K)", "ne_0(cm-3)", "neRms(cm-3)"
                endif
                write(55,'(i5.4,f20.1, a13, 4(es12.5,1x))') grid%idump, grid%currenttime, trim(weighting), t0, trms, n0, nrms
                close(55)
#ifdef MPI
             endif
#endif
          endif
       enddo
    endif

#ifdef MPI
    if (plotAvgTdust .or. calculateGlobalAvgTdust) then
       do i=1,2
          if (i==1) then
             write(weighting, '(a)') "dustmass"
          elseif (i==2) then
             write(weighting, '(a)') "none" 
          endif
#ifdef USECFITSIO
          if (plotAvgTdust) then
             ! pixel-by-pixel average
             if (writeoutput) write(*,*) "Making tdust image with weighting ", weighting
             ! calculate mean, save to image
             call createTdustColumnImage(grid, columnImageDirection, image, trim(weighting))
             write(thisFile, '(i4.4,a,a,a)') grid%idump, "_avgtdust_", trim(weighting), ".fits"
             if (writeoutput) call writeFitsColumnDensityImage(image, thisFile)
             ! calculate normalised rms
             call createTdustColumnImage(grid, columnImageDirection, rmsImage, trim(weighting), image)
             rmsImage = rmsImage / image**2 ! <(T-Tavg)^2> / Tavg^2
             write(thisFile, '(i4.4,a,a,a)') grid%idump, "_normvarianceTdust_", trim(weighting), ".fits"
             if (writeoutput) call writeFitsColumnDensityImage(rmsImage, thisFile)
          endif
#else
          call writeInfo("FITS not enabled, not writing",FORINFO)
#endif

          if (calculateGlobalAvgTdust) then
             ! global average
             if (writeoutput) write(*,*) "Calculating avg tdust with weighting ", weighting
             sigma = 0.d0
             total = 0.d0
             sigmaNe = 0.d0
             call calculateAverageTdust(grid%octreeRoot, grid,sigma, total, trim(weighting), .false., 0.d0, sigmaNe)
             call MPI_ALLREDUCE(sigma, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             sigma = tempDouble
             call MPI_ALLREDUCE(total, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             total = tempDouble
             t0 = sigma/total ! mean

             ! rms 
             sigma = 0.d0
             total = 0.d0
             sigmaNe = 0.d0
             call calculateAverageTdust(grid%octreeRoot, grid, sigma, total, trim(weighting), .true., t0, sigmaNe)
             call MPI_ALLREDUCE(sigma, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             sigma = tempDouble
             call MPI_ALLREDUCE(total, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             total = tempDouble
             ! trms  = sqrt(<(T-Tavg)^2>)
             trms = sqrt(sigma/total)
             if (myrankglobal == 1) then
                write(thisFile, '(a)') "avgtdustGlobalRms.dat"
                open(55, file=thisFile, status="unknown", position="append")
                if (firstTime) then
                   write(55,'(a5,a20,a13,2(a12,1x))') "#dump", "t(s)", "weight", "Td0(K)", "TdRms(K)"
                endif
                write(55,'(i5.4,f20.1, a13, 2(es12.5,1x))') grid%idump, grid%currenttime, trim(weighting), t0, trms
                close(55)
             endif
          endif
       enddo
    endif

! emission measure image
#ifdef USECFITSIO
    if (calculateEmissionMeasure) then
!       write(thisFile,'(a, i4.4, a)') "emissionMeasurez_", grid%iDump,".fits"
!       call createEmissionMeasureImage(grid, VECTOR(0.d0, 0.d0, 1.d0), image)
!       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))
!
!       write(thisFile,'(a, i4.4, a)') "emissionMeasurey_", grid%iDump,".fits"
!       call createEmissionMeasureImage(grid, VECTOR(0.d0, 1.d0, 0.d0), image)
!       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))
!
!       write(thisFile,'(a, i4.4, a)') "emissionMeasurex_", grid%iDump,".fits"
!       call createEmissionMeasureImage(grid, VECTOR(1.d0, 0.d0, 0.d0), image)
!       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))

! integral(P) (nb need to divide by ncells for mean along column)
! Fig 8 in Ali 2020

!       write(thisFile,'(a, i4.4, a)') "prad_colz_", grid%iDump,".fits"
!       call createPressureImage(grid, VECTOR(0.d0, 0.d0, 1.d0), image, 'prad')
!       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))
!
!       write(thisFile,'(a, i4.4, a)') "prad_coly_", grid%iDump,".fits"
!       call createPressureImage(grid, VECTOR(0.d0, 1.d0, 0.d0), image, 'prad')
!       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))
!
!       write(thisFile,'(a, i4.4, a)') "prad_colx_", grid%iDump,".fits"
!       call createPressureImage(grid, VECTOR(1.d0, 0.d0, 0.d0), image, 'prad')
!       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))
!
!       write(thisFile,'(a, i4.4, a)') "pgas_colz_", grid%iDump,".fits"
!       call createPressureImage(grid, VECTOR(0.d0, 0.d0, 1.d0), image, 'pgas')
!       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))
!
!       write(thisFile,'(a, i4.4, a)') "pgas_coly_", grid%iDump,".fits"
!       call createPressureImage(grid, VECTOR(0.d0, 1.d0, 0.d0), image, 'pgas')
!       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))
!
!       write(thisFile,'(a, i4.4, a)') "pgas_colx_", grid%iDump,".fits"
!       call createPressureImage(grid, VECTOR(1.d0, 0.d0, 0.d0), image, 'pgas')
!       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))
!
!       write(thisFile,'(a, i4.4, a)') "pradOverPgas_colz_", grid%iDump,".fits"
!       call createPressureImage(grid, VECTOR(0.d0, 0.d0, 1.d0), image, 'ratio')
!       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))
!
!       write(thisFile,'(a, i4.4, a)') "pradOverPgas_coly_", grid%iDump,".fits"
!       call createPressureImage(grid, VECTOR(0.d0, 1.d0, 0.d0), image, 'ratio')
!       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))
!
!       write(thisFile,'(a, i4.4, a)') "pradOverPgas_colx_", grid%iDump,".fits"
!       call createPressureImage(grid, VECTOR(1.d0, 0.d0, 0.d0), image, 'ratio')
!       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))

       ! pdir from lbol (lopez)
       write(thisFile,'(a, i4.4, a)') "pbol_colz_", grid%iDump,".fits"
       call createPressureImage(grid, VECTOR(0.d0, 0.d0, 1.d0), image, 'pbol')
       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))

       write(thisFile,'(a, i4.4, a)') "pbol_coly_", grid%iDump,".fits"
       call createPressureImage(grid, VECTOR(0.d0, 1.d0, 0.d0), image, 'pbol')
       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))

       write(thisFile,'(a, i4.4, a)') "pbol_colx_", grid%iDump,".fits"
       call createPressureImage(grid, VECTOR(1.d0, 0.d0, 0.d0), image, 'pbol')
       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))

       write(thisFile,'(a, i4.4, a)') "pradOverPbol_colz_", grid%iDump,".fits"
       call createPressureImage(grid, VECTOR(0.d0, 0.d0, 1.d0), image, 'pradOverPbol')
       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))

       write(thisFile,'(a, i4.4, a)') "pradOverPbol_coly_", grid%iDump,".fits"
       call createPressureImage(grid, VECTOR(0.d0, 1.d0, 0.d0), image, 'pradOverPbol')
       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))

       write(thisFile,'(a, i4.4, a)') "pradOverPbol_colx_", grid%iDump,".fits"
       call createPressureImage(grid, VECTOR(1.d0, 0.d0, 0.d0), image, 'pradOverPbol')
       if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))
   endif
#else
    call writeInfo("FITS not enabled, not writing",FORINFO)
#endif


!! luminosities
    if (writeLums) then
       call writeLuminosities(grid)
    endif

    if (findNUndersampled) then
       call findNumberUndersampled(grid, nSampled, nUndersampled, nCells)
          if (myrankglobal == 1) then
             write(thisFile, '(a)') "nsampled.dat"
             open(56, file=thisFile, status="unknown", position="append")
             if (firstTime) then
                write(56,'(a5,a20,3(1x,a10))') "#dump", "t(s)", "nSampled", "nUnderSamp", "nCells"
             endif
             write(56,'(i5.4,f20.1,3(1x,i10))') grid%idump, grid%currenttime, nSampled, nUndersampled, nCells 
             close(56)
          endif
    endif
#endif
#endif
! end mpi
    if (radiusTotals.and..not.splitOverMPI) then
!       ! calculate totals inside specified radii
       radii = (/ 3.d0, 5.d0, 10.d0, 15.d0, 20.d0 /) ! pc
       do i = 1, size(radii)
! this can be uncommented and it works
!          sums(:) = 0.d0
!          call calculateTotalInRadius(grid%octreeRoot, grid, sums, radii(i)*pctocm/1.d10, VECTOR(0.d0,0.d0,0.d0))
!          do iSource = 1, nSource
!             if (modulus(globalSourceArray(iSource)%position) <= radii(i)*pctocm/1.d10) then
!                ! sink mass (msol)
!                sums(7) = sums(7) + globalSourceArray(iSource)%mass/msol
!                ! star mass (msol)
!                do j = 1, globalSourceArray(iSource)%nSubsource
!                   sums(8) = sums(8) + globalsourceArray(iSource)%subsourceArray(j)%mass/msol
!                enddo
!                ! luminosity (lsol)
!                sums(9) = sums(9) + globalSourceArray(iSource)%luminosity/lsol
!             endif
!          enddo
!          flux = fluxThroughSphericalSurface(grid, radii(i)*pctocm/1.d10)
!          flux = flux / msol * 365.25d0*24.d0*3600.d0
!          write(thisFile, '(a,i2.2,a)') "totinradius_", int(radii(i)), "pc.dat"
!          open(56, file=thisFile, status="unknown", position="append")
!          if (firstTime) then
!             write(56,'(a5,a20,10a13)') "#dump", "t(s)", "mass(msol)", "pgas(d/cm2)", "prad(d/cm2)", "mom(gcm/s)", &
!             "ionmass", "ionvol(pc3)", "msink", "mstars", "lsink(lsol)", "flux(sol/yr)"
!          endif
!          write(56,'(i5.4,f20.1,10es13.5)') grid%idump, grid%currenttime, sums(1:9), -flux
!          close(56)

          ! VELOCITY DISPERSION
          vmean = 0.d0
          vms = 0.d0
          total = 0.d0
          mion = 0.d0
          mneu = 0.d0
          call calculateMeanVelocity(grid%octreeRoot, grid, vmean, vms, total, mion, mneu, radii(i)*pctocm/1.d10)
          ! ionized and neutral
          vmean(1:2) = vmean(1:2)/total
          vms(1:2) = vms(1:2)/total 
          ! ionised only
          mion = max(mion, 1.d-50)
          vmean(3:4) = vmean(3:4)/mion
          vms(3:4) = vms(3:4)/mion
          ! neutral only
          mneu = max(mneu, 1.d-50)
          vmean(5:6) = vmean(5:6)/mneu
          vms(5:6) = vms(5:6)/mneu
          ! ionised temperature
          vmean(7) = vmean(7)/mion
          vms(7) = vms(7)/mion
          ! vrms is sqrt(vms) 
          ! standard deviation
          sigmav = sqrt(vms - vmean**2)
          ! write to file
!          write(thisFile, '(a,i2.2,a)') "vdispersion_", int(radii(i)), "pc.dat"
          write(thisFile, '(a,i2.2,a)') "velMeanStd_", int(radii(i)), "pc.dat"
          open(56, file=thisFile, status="unknown", position="append")
          if (firstTime) then
             write(56,'(a5,a20,17a13)') "#dump", "t(s)", &
             "vModMean", "vModStd",&
             "vrMean",   "vrStd",&
             "vModMean_II", "vModStd_II",&
             "vrMean_II",   "vrStd_II",&
             "vModMean_I", "vModStd_I",&
             "vrMean_I",   "vrStd_I",&
             "TMean_II",   "TStd_II",&
             "mtot", "mion", "mneu"
          endif
          write(56,'(i5.4,f20.1, 17es13.5)') grid%idump, grid%currenttime, &
          ! velocities in km/s
          vmean(1)/1e5, sigmav(1)/1e5, &
          vmean(2)/1e5, sigmav(2)/1e5, &
          vmean(3)/1e5, sigmav(3)/1e5, &
          vmean(4)/1e5, sigmav(4)/1e5, &
          vmean(5)/1e5, sigmav(5)/1e5, &
          vmean(6)/1e5, sigmav(6)/1e5, &
          ! temperature in K
          vmean(7), sigmav(7), &
          total/msol, mion/msol, mneu/msol
          close(56)
       enddo
    endif

    if (clusterRadial.and..not.splitOverMPI) then
       ! radial averages from cluster centre of mass
       if (writeoutput) write(*,*) "Calculating radial profiles from COM of nsource ", nsource
       rcom = VECTOR(0.d0,0.d0,0.d0)
       if (nsource > 0) then
          do i = 1, nsource
             rcom = rcom + source(i)%mass * source(i)%position
          enddo
          rcom = rcom / sum(source(1:nsource)%mass)
       endif
       call shellAverage(grid, rcom, .false.)

       ! radial averages from source positions 
       do isource = 1, nsource
          if (source(isource)%luminosity > 0.d0) then
             if (writeoutput) write(*,*) "Calculating radial profiles from source ", isource, source(isource)%luminosity/lsol
             call shellAverage(grid, source(isource)%position, .true., isource=isource)
          endif
       enddo

       ! total pressure on grid (counts pgas in ionized cells only)
       pgas = 0.d0; prad = 0.d0
       call calculateTotalPressure(grid%octreeRoot, grid, pgas, prad)
       write(thisFile, '(a)') "totpressure.dat"
       open(56, file=thisFile, status="unknown", position="append")
       if (firstTime) then
          write(56,'(a5,a20,2a13)') "#dump", "t(s)", "pgas", "prad"
       endif
       write(56,'(i5.4,f20.1,2es13.5)') grid%idump, grid%currenttime, pgas, prad
       close(56)
    endif

#ifdef MPI
    if (calculateAvgPressure) then
       do j = 1,2
          if (j == 1) then
             write(comp, '(a)') "gas" 
          elseif (j == 2) then
             write(comp, '(a)') "rad" 
          endif
          do i = 2,4
             if (i==1) then
                write(weighting, '(a)') "mass" 
             elseif (i==2) then
                write(weighting, '(a)') "none" 
             elseif (i==3) then
                write(weighting, '(a)') "ionNone" 
             elseif (i==4) then
                write(weighting, '(a)') "neuNone" 
             endif

             ! global average
             if (writeoutput) write(*,*) "Calculating avg pressure with weighting ", weighting
             sigma = 0.d0
             total = 0.d0
             sigmaNe = 0.d0
             call calculateAveragePressure(grid%octreeRoot, grid,sigma, total, trim(weighting), .false., 0.d0, sigmaNe, 0.d0, comp)
#ifdef MPI
             call MPI_ALLREDUCE(sigma, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             sigma = tempDouble ! sum(T w dV)
             call MPI_ALLREDUCE(total, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             total = tempDouble ! sum(w dV)
             call MPI_ALLREDUCE(sigmaNe, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             sigmaNe = tempDouble ! sum(ne w dV)
#endif
             t0 = sigma/total ! mean temperature
             n0 = sigmaNe/total ! mean electron density

             ! std dev / rms deviation
             sigma = 0.d0
             total = 0.d0
             sigmaNe = 0.d0
#ifdef MPI
             call calculateAveragePressure(grid%octreeRoot, grid, sigma, total, trim(weighting), .true., t0, sigmaNe, n0, comp)
             call MPI_ALLREDUCE(sigma, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             sigma = tempDouble ! sum( (T-T0)^2 w dV)
             call MPI_ALLREDUCE(total, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             total = tempDouble ! sum( w dV)
             call MPI_ALLREDUCE(sigmaNe, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             sigmaNe = tempDouble ! sum( (ne-n0)^2 w dV)
#endif
             ! std dev = sqrt(<(T-T0)^2>)
             trms = sqrt(sigma/total)
             nrms = sqrt(sigmaNe/total)

             ! TODO calculate total pressure
#ifdef MPI
             if (myrankglobal == 1) then
#endif
                write(thisFile, '(3a)') "avgpressure_",trim(comp),".dat"
                open(55, file=thisFile, status="unknown", position="append")
                if (firstTime) then
                   write(55,'(a5,a20,a9,4(a13,1x))') "#dump", "t(s)", "weight","P_0(dy/cm2)","Prmsd","ne_0(cm-3)","neRmsd(cm-3)"
                endif
                write(55,'(i5.4,f20.1,a9,4(es13.5,1x))') grid%idump, grid%currenttime, trim(weighting),t0,trms,n0,nrms
                close(55)
#ifdef MPI
             endif
#endif
          enddo
       enddo
    endif
#endif


! assumes 1 thread only
!   mass flux through spherical surface
!    if (.not. splitOverMPI) then
!       write(thisFileGrid, '(a)') "valuesOnGrid.dat"
!       write(thisFileRadius, '(a)') "valuesWithinRadius.dat"
!
!       if (edgeRadius > 1.d0) then
!          ! use input variable
!          radius = edgeRadius
!       else
!          ! calculate the radius from the bursttime dump and save for future dumps
!          if (firstTime) then
!             if (grid%currentTime*secsToYears < burstTime*1.1d0) then
!                call findCloudEdge(grid%octreeRoot, radius)
!                ! replace file with header 
!                ! entire grid
!                open(69, file=thisFileGrid, status="replace")
!                   fm = '(a6,1x,a9,1x,a20,1x,a20,1x,7(a16,1x))'
!                   write(69, fm) "# dump", "[1]L(pc)", "[2]t(s)", "[3]Mtot(msol)", "[4]flux(msol/yr)", & 
!                     "[5]vol(pc3)", "[6]ionMass(msol)", "[7]ionVol(pc3)", "[8]rhomax(g/cc)", "[9]KE(erg)", &
!                     "[10]Eth(erg)"
!                close(69)
!                ! within radius
!                open(69, file=thisFileRadius, status="replace")
!                   fm = '(a6,1x,a9,1x,a20,1x,a20,1x,7(a16,1x))'
!                   write(69, fm) "# dump", "[1]R(pc)", "[2]t(s)", "[3]Mtot(msol)", "[4]flux(msol/yr)", & 
!                     "[5]vol(pc3)", "[6]ionMass(msol)", "[7]ionVol(pc3)", "[8]rhomax(g/cc)", "[9]KE(erg)", &
!                     "[10]Eth(erg)"
!                close(69)
!             else
!                write(*,*) "Input variable edgeradius not defined, and first dump isn't near burstTime so haven't calculated edgeradius"
!                stop
!             endif
!          endif
!       endif
!       write(*,*) "radius ", radius, " 1e10 cm / ", radius*1.d10/pcToCm, " pc"
!
!       ! entire grid
!       mass = 0.d0
!       highestRho = 1.d-30
!       volume  = 0.d0
!       ionizedMass = 0.d0
!       ionizedVolume = 0.d0
!       totalKE = 0.d0
!       totalTE = 0.d0
!       call findMassWithBounds(grid%octreeRoot, mass, highestRho=highestRho, totalVolume=volume,&
!                ionizedMass=ionizedMass, ionizedVolume=ionizedVolume, totalKE=totalKE, totalTE=totalTE)
!
!       massflux = fluxThroughSphericalSurface(grid, amrgridsize/2.d0)
!       massflux = massflux / msol * 365.25d0*24.d0*3600.d0 ! Msol/yr
!       ! write data
!       open(69, file=thisFileGrid, status="old", position="append")
!          fm = '(i6.4,1x,f9.5,1x,f20.1,1x,f20.8,1x,7(es16.8,1x))'
!          write(69, fm) grid%idump, amrgridsize/2.d0*1.d10/pcToCm, grid%currenttime, mass/msol, massflux, &
!               volume/(pctocm**3), ionizedMass/msol, ionizedVolume/(pctocm**3), highestRho, totalKE, totalTE 
!       close(69)
!
!       ! within radius 
!       mass = 0.d0
!       highestRho = 1.d-30
!       volume  = 0.d0
!       ionizedMass = 0.d0
!       ionizedVolume = 0.d0
!       totalKE = 0.d0
!       totalTE = 0.d0
!       call findMassWithBounds(grid%octreeRoot, mass, highestRho=highestRho, totalVolume=volume,&
!                ionizedMass=ionizedMass, ionizedVolume=ionizedVolume, totalKE=totalKE, totalTE=totalTE,&
!                maxRadius=radius)
!
!       massflux = fluxThroughSphericalSurface(grid, radius)
!       massflux = massflux / msol * 365.25d0*24.d0*3600.d0 ! Msol/yr
!
!       ! write data
!       open(69, file=thisFileRadius, status="old", position="append")
!          fm = '(i6.4,1x,f9.5,1x,f20.1,1x,f20.8,1x,7(es16.8,1x))'
!          write(69, fm) grid%idump, radius*1.d10/pcToCm, grid%currenttime, mass/msol, massflux, &
!               volume/(pctocm**3), ionizedMass/msol, ionizedVolume/(pctocm**3), highestRho, totalKE, totalTE 
!       close(69)
!    else
!       if (writeoutput) write(*,*) "clusterAnalysis requires not splitting over mpi"
!       stop
!    endif

! average habing flux over accretionRadius of first source
    if (findHabing) then
       if (.not. splitOverMPI) then
          write(thisFile, '(a,i3.3,a)') "averageHabingFlux_allStars_primarySource",primarySource,".dat"
          ! write header
          if (firstTime) then
             open(69, file=thisFile, status="unknown", form="formatted", position="append")
                write(69, '(a6,a20,a5,a9,3(1x,a12))') "#dump", "t(s)", "i*", "M*(Msol)", "dToM*p(pc)", "G0inCell(H)", "avgG0(H)"
             close(69)
          endif
          open(69, file=thisFile, status="old", position="append", form="formatted")
          do iSource = 1, globalnSource
             weightedFluxInRadius = 0.d0
             massInRadius = 0.d0
             call calculateAverageHabingFlux(grid%octreeRoot, globalSourceArray(iSource), 2.5d0*smallestCellSize, &
                     weightedFluxInRadius, massInRadius, g0inCell)
             meanG0 = weightedFluxInRadius/massInRadius
             ! distance to primary star
             if (iSource /= primarySource) then
                distance = modulus(globalSourceArray(iSource)%position - globalSourceArray(primarySource)%position)*1.d10/pctocm
             else
                distance = 0.d0
             endif

             ! write data
              write(69, '(i6.4,f20.2, i5.3, f9.4, 3(1x,es12.5))') grid%idump, grid%currentTime-burstTime, iSource, & 
                        globalSourceArray(iSource)%mass/msol, distance, g0InCell, meanG0 
          enddo
          close(69)
       endif
    endif

    ! tables for tom haworth (FUV vs time for clustersinks)
    if (findHabingCluster .and. clustersinks) then
       if (.not. splitOverMPI) then
          trad0 = 25425924913214.7 ! fixme t=0 hardcoded for metallicity models (Ali 2020). Time rad starts.
          ! calculate emitted luminosity
          allocate(lum(1:globalnSource))
          lum(:) = 0.d0
          do i = 1, globalnSource
             if (globalSourceArray(i)%luminosity > 0.d0) then
                lum(i) = integrateSpectrumOverBand(globalSourceArray(i)%spectrum, 912.d0, 2400.d0) &
                      * (fourPi*(1.d10*globalSourceArray(i)%radius)**2)
             endif
          enddo

          ! calculate flux received at each sink i from all other sinks j
          do i = 1, globalnSource
             flux = 0.d0
             do j = 1, globalnSource
                if (i == j) then
                   distance = globalSourceArray(i)%accretionRadius / 2.5d0 / 2.d0
                else
                   distance = modulus(globalSourceArray(i)%position - globalSourceArray(j)%position) * 1.d10
                endif
                flux = flux + lum(j) / (fourPi * distance**2)
             enddo
             flux = flux / habing

             write(thisFile, '(a,i3.3,a)') "fuv_dilution_",i,".dat"
             ! write header
             if (firstTimeStars(i,0)) then
                open(69, file=thisFile, status="unknown", form="formatted", position="append")
                   write(69, '(a5, 2(1x,a10), 5(1x,a12),a6)') "dump", "time", "mass", "lum", "x", "y", "z", &
                   "g0_dil","nstar"
                close(69)
                firstTimeStars(i,0) = .false.
             endif

             open(69, file=thisFile, status="old", position="append", form="formatted")
                write(69, '(i5, 2(1x,f10.5), 5(1x,es12.5), i6)') iModel, (globalSourceArray(i)%time-trad0)*secsToYears/1e6, &
                globalSourceArray(i)%mass/msol, & 
                globalSourceArray(i)%luminosity/lsol, & 
                globalSourceArray(i)%position%x * 1.d10/pctocm, & 
                globalSourceArray(i)%position%y * 1.d10/pctocm, & 
                globalSourceArray(i)%position%z * 1.d10/pctocm, & 
                flux, globalSourceArray(i)%nsubsource 
             close(69)
          enddo
          deallocate(lum)

          

! stuff below works, just commented out
! get formation time relative to time rad starts (only need to do this for end dump)
!          do i = 1, globalnSource
!             if (globalSourceArray(i)%nsubsource > 0) then
!                write(thisFile, '(a,i3.3,a)') "birth_",i,".dat"
!                open(69, file=thisFile, status="unknown", form="formatted", position="append")
!                write(69, '(2a11)') "mini", "tborn"
!                do j = 1, globalSourceArray(i)%nsubsource
!                   write(69, '(2f11.5)') globalSourceArray(i)%subsourceArray(j)%initialMass, &
!                   (grid%currentTime - globalSourceArray(i)%subsourceArray(j)%age/secsToYears - trad0)*secsToYears/1e6
!                enddo
!                close(69)
!             endif
!          enddo

! stuff below works, just commented out
! write out g0 at/near every clustersink
!          do i = 1, globalnSource
!             write(thisFile, '(a,i3.3,a)') "fuv_",i,".dat"
!             ! write header
!             if (firstTimeStars(i,0)) then
!                open(69, file=thisFile, status="unknown", form="formatted", position="append")
!                   write(69, '(a5, 2(1x,a10), 6(1x,a12),a6)') "dump", "time", "mass", "lum", "x", "y", "z", "g0_cell", "g0_avg", &
!                   "nstar"
!                close(69)
!                firstTimeStars(i,0) = .false.
!             endif
!
!             weightedFluxInRadius = 0.d0
!             massInRadius = 0.d0
!             call calculateAverageHabingFlux(grid%octreeRoot, globalSourceArray(i), globalSourceArray(i)%accretionRadius/1.d10, &
!                     weightedFluxInRadius, massInRadius, g0inCell)
!             meanG0 = weightedFluxInRadius/massInRadius
!
!             open(69, file=thisFile, status="old", position="append", form="formatted")
!                write(69, '(i5, 2(1x,f10.5), 6(1x,es12.5), i6)') grid%idump, (grid%currentTime-trad0)*secsToYears/1e6, &
!                globalSourceArray(i)%mass/msol, & 
!                globalSourceArray(i)%luminosity/lsol, & 
!                globalSourceArray(i)%position%x * 1.d10/pctocm, & 
!                globalSourceArray(i)%position%y * 1.d10/pctocm, & 
!                globalSourceArray(i)%position%z * 1.d10/pctocm, & 
!                g0InCell, meanG0, globalSourceArray(i)%nsubsource 
!             close(69)
!
!             ! massive subsources
!             do j = 1, globalSourceArray(i)%nsubsource
!                if (globalSourceArray(i)%subsourceArray(j)%mass > 18.d0*msol) then
!                   write(thisFile, '(a,i3.3,a,i3.3,a)') "fuv_",i,"_",j,".dat"
!                   ! write header
!                   if (firstTimeStars(i,j)) then
!                      open(69, file=thisFile, status="unknown", form="formatted", position="append")
!                         write(69, '(a5, 2(1x,a10), 2a12)') "dump", "time", "mass", "lum", "teff"
!                      close(69)
!                      firstTimeStars(i,j) = .false.
!                   endif
!                   open(69, file=thisFile, status="old", position="append", form="formatted")
!                      write(69, '(i5, 2(1x,f10.5), 2es12.5)') grid%idump, (grid%currentTime-trad0)*secsToYears/1e6, &
!                      globalSourceArray(i)%subsourceArray(j)%mass/msol, & 
!                      globalSourceArray(i)%subsourceArray(j)%luminosity/lsol, &
!                      globalSourceArray(i)%subsourceArray(j)%teff
!                   close(69)
!                endif
!             enddo
!         enddo
      endif
   endif

    ! calculate the mean HII region radius, and then calculate mean values inside this radius
    if (findHIIRadius) then
       if (.not. splitOverMPI) then
          ! find most massive star
          mmax = 1e-30
          imax = 0
          jmax = 0
          do i = 1, globalnSource
             if (globalSourceArray(i)%luminosity > 0.) then
                if (clustersinks) then
                   do j = 1, globalSourceArray(i)%nsubsource
                      if (globalSourceArray(i)%subsourceArray(j)%mass > mmax) then
                         mmax = globalSourceArray(i)%subsourceArray(j)%mass
                         imax = i
                         jmax = j
                         rvec = globalSourceArray(i)%position
                      endif
                   enddo
                else
                   if (globalSourceArray(i)%mass > mmax) then
                      mmax = globalSourceArray(i)%mass
                      imax = i
                      rvec = globalSourceArray(i)%position + VECTOR(grid%halfSmallestSubcell*1.d-3, 0.d0, 0.d0)
                   endif
                endif
             endif
          enddo
          if (imax > 0) then
             write(*,*) "most massive star i,j", imax, jmax, " mass ", mmax/msol, " position ", rvec

             ! tag HII region cells and find its mean radius (with the most massive star as the origin)
             call identifyHiiRegion(grid%octreeRoot, grid, rvec)

             ! calculate median and IQR in the tagged hii region
             call hiiMedianPressure(grid, firstTime, imax)

! the below commented out code calculates the mean and std dev in the tagged HII region
!             ! calculate the mean pressure using cells which have been tagged as belonging to the HII region
!             sums(:) = 0.d0
!             call calculateTotalInRadius(grid%octreeRoot, grid, sums, 1.d50, rvec, tagged=.true.)
!             ! mean = sum(x)/Ncells
!             pgas = sums(2)/sums(10)
!             prad = sums(3)/sums(10)
!             ! stddev**2= mean of squares - square of mean 
!             sigmapgas = sqrt(sums(11)/sums(10) - pgas**2)
!             sigmaprad = sqrt(sums(12)/sums(10) - prad**2)
!
!
!
!             write(thisFile, '(a)') "hiiregion_mean.dat"
!             open(69, file=thisFile, status="unknown", form="formatted", position="append")
!             if (firstTime) then
!                write(69, '(a)') "# pgas and prad (dyn/cm2) are mean = sum(P)/N. ionmass (msol) and ionvol (pc3) are total values"
!                write(69,'(a5,a20,6a13,a5)') "dump", "t", "pgas_mean", "pgas_std", "prad_mean", "prad_std", &
!                "ionmass", "ionvol", "sink"
!             endif
!             write(69,'(i5.4,f20.1,6es13.5,i5)') grid%idump, grid%currenttime, pgas, sigmapgas, prad, sigmaprad, & 
!             sums(5), sums(6), imax
!             close(69)
!
!             write(thisFile, '(a,i4.4,a)') "hiiregion_", grid%idump, ".vtk"
!             call writeVtkFile(grid, thisFile, &
!                  valueTypeString=(/"HI           ", "etacont      " /))
          endif
       endif
    endif
    
    
#ifdef USECFITSIO
!    ! column along x, ionfrac along x and y 
!    ! x dir
!    thisDir = VECTOR(1.d0, 0.d0, 0.d0)
!!    write(thisFile, '(a,i4.4,a)') "columnx_", grid%idump, ".fits"
!!    call createColumnDensityImage(grid, thisDir, image)
!!    if (writeoutput) call writeFitsColumnDensityImage(image, thisFile)
!    write(thisFile, '(a,i4.4,a)') "ionx_", grid%idump, ".fits"
!    call createIonizationImage(grid, thisDir, image)
!    if (writeoutput) call writeFitsColumnDensityImage(image, thisFile)
!    write(thisFile, '(a,i4.4,a)') "hix_", grid%idump, ".fits"
!    call createHiImage(grid, thisDir, image)
!    if (writeoutput) call writeFitsColumnDensityImage(image, thisFile)
!
!    ! z dir
!    thisDir = VECTOR(0.d0, 0.d0, 1.d0)
!!    write(thisFile, '(a,i4.4,a)') "columnz_", grid%idump, ".fits"
!!    call createColumnDensityImage(grid, thisDir, image)
!!    if (writeoutput) call writeFitsColumnDensityImage(image, thisFile)
!    write(thisFile, '(a,i4.4,a)') "ionz_", grid%idump, ".fits"
!    call createIonizationImage(grid, thisDir, image)
!    if (writeoutput) call writeFitsColumnDensityImage(image, thisFile)
!    write(thisFile, '(a,i4.4,a)') "hiz_", grid%idump, ".fits"
!    call createHiImage(grid, thisDir, image)
!    if (writeoutput) call writeFitsColumnDensityImage(image, thisFile)
#endif

    if (calculateLymanFlux) then
       if (.not. splitOverMPI) then
          nlyA = 0.d0 
          nlyB = 0.d0 
          nlyR = 0.d0
          call estimateLymanFlux(grid%octreeRoot, grid, nlyA, nlyB, nlyR) 

          write(thisFile, '(a)') "estimatedLymanFlux.dat"
          ! write header
          if (firstTime) then
             open(69, file=thisFile, status="replace", form="formatted")
                write(69, '(a6,a20, 3(1x,a12))') "#dump", "t(s)", "nly(s-1)(A)", "nly(s-1)(B)", "nly(s-1)(R)"
             close(69)
          endif
          ! write data
          open(69, file=thisFile, status="old", position="append", form="formatted")
              write(69, '(i6.4,f20.2, 3(1x,es12.5))') grid%idump, grid%currentTime, nlyA, nlyB, nlyR 
          close(69)
       endif
    endif

    ! calculate tau in FUV band (due to dust) between star 1 and all other stars
    if (calculateTauFUV) then
       if (.not. splitOverMPI) then


          write(thisFile, '(a,i3.3,a)') "tauFUVBetweenStars_primarySource",primarySource,".dat"
          if (firstTime) then
             open(69, file=thisFile, status="unknown", form="formatted", position="append")
                write(69, '(a6,a20,a5,4(1x,a12))') "#dump", "t(s)", "i*", "tauAbs", "tauSca", "tauExt", "colDen(g/cm2)"
             close(69)
          endif
          do iSource = 1, globalnSource
             if (iSource /= primarySource) then

                tauAbs = 0.d0; tauSca = 0.d0; tauExt = 0.d0; columnDensity = 0.d0
                call meanTauBetweenPoints(912., 2400., grid, globalSourceArray(primarySource)%position, &
                   globalSourceArray(iSource)%position, tauAbs, tauSca, tauExt, columnDensity)

                open(69, file=thisFile, status="old", position="append", form="formatted")
                   write(69, '(i6.4,f20.2, i5.3, 4(1x,es12.5))') grid%iDump, grid%currentTime, iSource, tauAbs, tauSca, &
                     tauExt, columnDensity
                close(69)

             endif
          enddo
       endif
    endif

    ! sum PAH emission and absorption
    if (.not. splitOverMPI) then
       write(*,*) "doing sumPAH"
       total = 0.d0 ! absorption rate
       flux = 0.d0 ! emission rate
       distance = 0.d0 ! sum of Adot for which emission==0
       call sumPAH(grid%octreeRoot, grid, total, flux, distance)
       write(*,*) "emi ", flux
       write(*,*) "abs ", total
       write(*,*) "emi/abs", flux/total
       write(*,*) "Adot > 0 but emi > 0 ", distance
       write(*,*) "   frac (/abs)", distance/total

    endif

    firstTime = .false.
  end subroutine clusterAnalysis
  
  recursive subroutine findCloudEdge(thisOctal, radius)
    use inputs_mod, only : surfacedensity, sphereRadius, beta
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: i, subcell
    real(double) :: radius, r, rhoOutside, rhoCore
    type(VECTOR) :: rVec

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findCloudEdge(child, radius) 
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then

             ! from geometry krumholz
             rhoCore = 6.d0 * surfaceDensity / ((2.d0**2.5d0 - 1.d0) * sphereRadius*1.d10)
             rhoOutside = 1.d-2 * 2.d0**beta * rhoCore
             if (thisOctal%rho(subcell) > 1.1d0*rhoOutside) then
                rVec = subcellCentre(thisOctal, subcell) 
                r = modulus(rVec)
                radius = max(radius, r)
             endif 

          endif
       endif
    enddo
  
  end subroutine findCloudEdge

  subroutine writeCellValuesForSource(grid, thisSource, position) 
     use inputs_mod, only : iModel
     use amr_mod, only : octalOnThread
     type(GRIDTYPE) :: grid
     type(OCTAL), pointer :: thisOctal
     type(SOURCETYPE) :: thisSource
     type(VECTOR) :: position
     integer :: subcell
     real(double) :: time, rho, habingFlux, tgas, tdust, rMod !, tval 
!     type(VECTOR) :: vHat
    character(len=80) :: thisFile

     if (myrankGlobal /= 0 .and..not. loadBalancingThreadGlobal) then
        thisOctal => grid%octreeRoot
        call findSubcellTD(position, grid%octreeRoot, thisOctal, subcell) 
        if (octalOnThread(thisOctal, subcell, myRankGlobal)) then

           time = grid%currentTime
           rho = thisOctal%rho(subcell)
           habingFlux = thisOctal%habingFlux(subcell)
           tgas = thisOctal%temperature(subcell)
           tdust = thisOctal%tdust(subcell)

!           vHat = VECTOR(thisSource%velocity%x, thisSource%velocity%y, thisSource%velocity%z)  
!           call normalize(vHat)
!           call distanceToCellBoundary(grid, thisSource%position, vHat, tval, thisOctal, subcell)
           rMod = modulus(thisSource%position - position)

           write(thisFile, '(a)') "firststartpoint.dat"
           open(69, file=thisFile, status="old", position="append", form="formatted")
           write(69, '(i6.4, f20.2, 3(1x,es12.5), 2(1x,f9.2))') iModel, time, rho, habingFlux*1.d10, rmod, tgas, tdust 
           close(69)

        endif
     endif
     
  end subroutine writeCellValuesForSource 
  recursive subroutine calculateAverageHabingFlux(thisOctal, source, radius, totalFlux, totalMass,cellFlux)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    type(vector) :: rcell
    real(double) :: totalFlux, totalMass, thisMass, radius, cellFlux
     type(SOURCETYPE) :: source


    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateAverageHabingFlux(child, source, radius, totalFlux, totalMass, cellFlux)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then

             rcell = subcellCentre(thisOctal, subcell)
             if (modulus(rcell - source%position) < radius) then
                thisMass = thisOctal%rho(subcell) * cellVolume(thisOctal, subcell)*1.d30
                totalFlux = totalFlux + thisMass*thisOctal%habingFlux(subcell)*1.d10 ! habingFlux is in 1e10 G0
                totalMass = totalMass + thisMass
             endif
             if (inSubcell(thisOctal, subcell, source%position)) then
                cellFlux = thisOctal%habingFlux(subcell)*1.d10
             endif 

          endif
       endif
    enddo


  end subroutine calculateAverageHabingFlux

  recursive subroutine estimateLymanFlux(thisOctal, grid, nlyA, nlyB, nlyR) 
    use photoion_utils_mod, only : recombRate
    use inputs_mod, only : honly
    type(octal), pointer   :: thisOctal
    type(GRIDTYPE) :: grid
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: nlyA, nlyB, ne, nhii, nheii, alphaAHi, alphaAHei, alphaB, dV, alphaR, nlyR


    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call estimateLymanFlux(child, grid, nlyA, nlyB, nlyR)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
             dv = cellVolume(thisOctal, subcell)*1.d30

             ne = thisOctal%ne(subcell)
             nHii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,2) * grid%ion(2)%abundance
             if(.not. hOnly) then
                nHeii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,4) * grid%ion(4)%abundance
             else
                nHeii = 0.d0
             end if

             ! case A recombination
             alphaAHi = recombrate(grid%ion(1), thisOctal%temperature(subcell))
             alphaAHei = recombrate(grid%ion(3), thisOctal%temperature(subcell))
             ! case B recombination
             alphaB = 2.7d-13
             ! rubin's
             alphaR = 4.1d-10 * thisOctal%temperature(subcell)**(-0.8d0)
             
             nlyA = nlyA + ne*(nHii*alphaAHi + nHeii*alphaAHei)*dV
             nlyB = nlyB + ne*(nHii + nHeii)*alphaB*dV
             nlyR = nlyR + ne*(nHii + nHeii)*alphaR*dV



          endif
       endif
    enddo
    
  end subroutine estimateLymanFlux 

  subroutine meanTauBetweenPoints(startLambda, endLambda, grid, startVec, endVec, tauAbs, tauSca, tauExt, columnDensity)
    use amr_mod, only : returnKappa
    type(GRIDTYPE) :: grid
    type(VECTOR) :: startVec, endVec, currentPosition, direction
    real(double), intent(out) :: tauAbs, tauSca, tauExt, columnDensity
    real(double) :: distToNextCell, totalDistance
    real :: startLambda, endLambda
    type(OCTAL), pointer :: thisOctal, sOctal
    real(double) :: fudgeFac = 1.d-3
    real(double) :: meankappaSca, meankappaAbs, meankappaExt !, rTot
    real(double), allocatable :: kappaAbsArray(:), kappaScaArray(:), dLambda(:)
    integer :: subcell
    integer :: iLamStart, iLamEnd, i
    logical :: endLoop
    logical, save :: firstTime=.true.

    tauAbs = 0.d0
    tauSca = 0.d0
    tauExt = 0.d0
    columnDensity = 0.d0
    totalDistance = 0.d0
    allocate(kappaAbsArray(1:grid%nlambda), kappaScaArray(1:grid%nlambda),dLambda(1:grid%nLambda))
    call locate(grid%lamArray, SIZE(grid%lamArray), startLambda, iLamStart)
    call locate(grid%lamArray, SIZE(grid%lamArray), endLambda, iLamEnd)
    iLamEnd=iLamEnd+1
    do i = 2, grid%nLambda-1
       dLambda(i) = (grid%lamArray(i+1)-grid%lamArray(i-1))/2.d0
    enddo
    dLambda(1) = 2.d0*(grid%lamArray(2)-grid%lamArray(1))
    dLambda(grid%nlambda) = 2.d0*(grid%lamArray(grid%nLambda)-grid%lamArray(grid%nLambda-1))

    currentPosition = startVec
    direction = endVec - startVec 
!    rTot = modulus(direction)
    call normalize(direction)

    CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)

    endLoop = .false.
    do while (inOctal(grid%octreeRoot, currentPosition) .and. .not.endLoop)

       call findSubcellLocal(currentPosition, thisOctal,subcell)

       call returnKappa(grid, thisOctal, subcell, kappaAbsArray=kappaAbsArray, kappaScaArray=kappaScaArray)

       meanKappaAbs = 0.d0; meanKappaSca=0.d0; meanKappaExt=0.d0
       do i = iLamStart, iLamEnd
          meanKappaAbs = meanKappaAbs + kappaAbsArray(i) * dLambda(i)
          meanKappaSca = meanKappaSca + kappaScaArray(i) * dLambda(i)
          meanKappaExt = meanKappaExt + (kappaAbsArray(i) + kappaScaArray(i)) * dLambda(i)
       enddo

       meanKappaAbs = meanKappaAbs / (grid%lamArray(iLamEnd) - grid%lamArray(iLamStart))
       meanKappaSca = meanKappaSca / (grid%lamArray(iLamEnd) - grid%lamArray(iLamStart))
       meanKappaExt = meanKappaExt / (grid%lamArray(iLamEnd) - grid%lamArray(iLamStart))

       sOctal => thisOctal
       call distanceToCellBoundary(grid, currentPosition, direction, DisttoNextCell, sOctal)

!       if (modulus(currentPosition + distToNextCell*direction) >= modulus(endVec - currentPosition)) then
       if (inSubcell(thisOctal, subcell, endVec)) then
          distToNextCell = modulus(endVec - currentPosition) 
          currentPosition = currentPosition + (distToNextCell)*direction
          endLoop = .true.
       else
          currentPosition = currentPosition + (distToNextCell+fudgeFac*grid%halfSmallestSubcell)*direction
       endif
       totalDistance = totalDistance + distToNextCell

       tauAbs = tauAbs + distToNextCell*meankappaAbs
       tauSca = tauSca + distToNextCell*meankappaSca
       tauExt = tauExt + distToNextCell*meankappaExt
       columnDensity = columnDensity + distToNextCell*1.d10*thisOctal%rho(subcell)

    end do
    deallocate(kappaAbsArray, kappaScaArray)

    if (firstTime) then
       write(*,*) "lamstart, lamend ", grid%lamArray(iLamStart), grid%lamArray(iLamEnd)
       write(*,*) "mean kappaAbs ", meanKappaAbs*1.d-10/thisOctal%rho(subcell)
       write(*,*) "mean kappaSca ", meanKappaSca*1.d-10/thisOctal%rho(subcell)
       write(*,*) "mean kappaExt ", meanKappaExt*1.d-10/thisOctal%rho(subcell)
       firstTime = .false.
    endif

!    write(*,*) "did ", startVec, " to ", currentPosition
    write(*,'(a,3(es12.3))') "Ended at ", currentPosition%x, currentPosition%y, currentPosition%z
    write(*,'(a,es12.3)') "distance travelled ", modulus(currentPosition-startVec)
    write(*,'(a,es12.3)') "distance betw star ", modulus(endVec-startVec)
  end subroutine meanTauBetweenPoints

  ! calculate mean and median pressure in radial bins extending out from the origin
  recursive subroutine averagePressureRadius(thisOctal, grid, weighting, origin, rbin, nbin, mean, vtot, values, ncells)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, nbin, ibin
    real(double) :: mean(:,:), vtot(:), rbin(:), values(:,:,:)
    integer :: ncells(:)
    character(len=*) :: weighting
    real(double) :: dx, r, fudgeFac
    type(VECTOR) :: origin, rvec


    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call averagePressureRadius(child, grid, trim(weighting), origin, rbin, nbin, mean, vtot, values, ncells)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
! not split over mpi
!             if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
                rVec = subcellCentre(thisOctal, subcell)
                dx = thisoctal%subcellsize ! torus unit
                fudgeFac = 1.d-1*grid%halfSmallestSubcell  

                r = modulus(rVec - origin)
                if (r > dx*rbin(nbin)+fudgefac) then
                   cycle
                elseif (r <= dx*rbin(1)+fudgeFac) then
                   call storePressureInR(thisOctal, subcell, 1, weighting, mean, vtot, values, ncells)
                else
                   do ibin = 2, nbin
                      if (r > dx*rbin(ibin-1)+fudgeFac .and. r <= dx*rbin(ibin)+fudgeFac) then
                         call storePressureInR(thisOctal, subcell, ibin, weighting, mean, vtot, values, ncells)
                         exit
                      endif
                   enddo
                endif

!             endif
          endif
       endif
    enddo
    
  end subroutine averagePressureRadius

  subroutine storePressureInR(thisOctal, subcell, ibin, weighting, mean, vtot, values, ncells)
    use ion_mod
    type(OCTAL), pointer :: thisOctal
    integer :: subcell, ibin
    real(double) :: mean(:,:), vtot(:), values(:,:,:)
    integer :: ncells(:)
    character(len=*) :: weighting
    real(double) :: dx, pgas,prad, mu
    real(double) :: weight, dv

    dv = cellVolume(thisOctal, subcell) * 1.d30 ! cm
    dx = thisoctal%subcellsize * 1.d10 ! cm

    mu = returnMu(thisOctal, subcell, globalIonArray, nGlobalIon)
    pgas = thisOctal%rho(subcell)*kerg*thisOctal%temperature(subcell)/(mu*mHydrogen)

    prad = dx * sqrt(thisOctal%radiationMomentum(subcell)%x**2 + &
                         thisOctal%radiationMomentum(subcell)%y**2 + &
                         thisOctal%radiationMomentum(subcell)%z**2)

    weight = 0.d0
    if (weighting=='mass') then
       weight = thisOctal%rho(subcell) * dv 

    ! V is ionized volume 
    elseif (weighting=='ionNone') then
       if (thisOctal%ionfrac(subcell, 2) > 0.99d0) then
          weight = 1.d0 
       else
          weight = 0.d0
       endif

    ! V is neutral volume
    elseif (weighting=='neuNone') then
       if (thisOctal%ionfrac(subcell, 1) > 0.99d0) then
          weight = 1.d0 
       else
          weight = 0.d0
       endif

    elseif (weighting=='none') then
       weight = 1.d0
    endif

    mean(1,ibin) = mean(1,ibin) + pgas*weight*dv
    mean(2,ibin) = mean(2,ibin) + prad*weight*dv
    mean(3,ibin) = mean(3,ibin) + (prad/pgas)*weight*dv
    mean(4,ibin) = mean(4,ibin) + thisOctal%rho(subcell)*weight*dv
    mean(5,ibin) = mean(5,ibin) + thisOctal%ionfrac(subcell, 2)*weight*dv

    vtot(ibin) = vtot(ibin) + weight*dv

    ncells(ibin) = ncells(ibin) + 1
    if (ncells(ibin) > size(values,dim=1)) then
       write(*,*) "ncells(ibin) > size(values)" 
       write(*,*) ibin, ncells(ibin), size(values,dim=1)
       stop
    endif
    if (weight >= 0.99d0) then
       values(ncells(ibin), 1, ibin) = pgas 
       values(ncells(ibin), 2, ibin) = prad
       values(ncells(ibin), 3, ibin) = prad/pgas 
       values(ncells(ibin), 4, ibin) = thisOctal%rho(subcell)
       values(ncells(ibin), 5, ibin) = thisOctal%ionfrac(subcell,2)
    endif
  end subroutine storePressureInR


  ! calculate mean and median pressure in all cells tagged with etacont == 1
  recursive subroutine averagePressure(thisOctal, grid, ivar, mean, vtot, values, ncells)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, ivar, ncells
    real(double) :: mean, vtot, values(:)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call averagePressure(child, grid, ivar, mean, vtot, values, ncells)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
! not split over mpi
!             if (octalOnThread(thisOctal, subcell, myRankGlobal)) then

                ! only take into account cells tagged with etacont = 1
                if (thisOctal%etaCont(subcell) == 1.d0) then
                   if (thisOctal%ionfrac(subcell, 1) > 0.9d0) then
                      write(*,*) "error - found neutral cell in averagePressure"
                      write(*,*) "etacont ", thisOctal%etacont(subcell)
                      write(*,*) "temperature ", thisOctal%temperature(subcell)
                      write(*,*) "radmom ", thisOctal%radiationMomentum(subcell)
                      stop
                   endif
                   call storePressure(thisOctal, subcell, ivar, mean, vtot, values, ncells)
                endif

!             endif
          endif
       endif
    enddo
    
  end subroutine averagePressure

  subroutine storePressure(thisOctal, subcell, ivar, mean, vtot, values, ncells)
    use ion_mod
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: mean, vtot, values(:)
    integer :: ncells, ivar
    real(double) :: dx, pgas,prad, mu
    real(double) :: weight, dv

    weight = 1.d0 

    dv = cellVolume(thisOctal, subcell) * 1.d30 ! cm
    dx = thisoctal%subcellsize * 1.d10 ! cm

    mu = returnMu(thisOctal, subcell, globalIonArray, nGlobalIon)
    pgas = thisOctal%rho(subcell)*kerg*thisOctal%temperature(subcell)/(mu*mHydrogen)

    prad = dx * sqrt(thisOctal%radiationMomentum(subcell)%x**2 + &
                         thisOctal%radiationMomentum(subcell)%y**2 + &
                         thisOctal%radiationMomentum(subcell)%z**2)

    if (prad == 0.d0 .or. pgas == 0.d0) then
       return
    endif

    select case (ivar)
       case(1)
          mean = mean + pgas*weight*dv
       case(2)
          mean = mean + prad*weight*dv
       case(3)
          mean = mean + (prad/pgas)*weight*dv
       case(4)
          mean = mean + thisOctal%rho(subcell)*weight*dv
       case(5)
          mean = mean + thisOctal%ionfrac(subcell, 2)*weight*dv
    end select

    vtot = vtot + weight*dv

    ncells = ncells + 1
    if (ncells > size(values)) then
       write(*,*) "ncells > size(values)" 
       write(*,*) ncells, size(values)
       stop
    endif
    if (ncells < 0 .or. ncells == huge(ncells)) then
       write(*,*) "ncells overflow ", ncells
       stop
    endif
    select case (ivar)
       case(1)
          values(ncells) = pgas 
       case(2)
          values(ncells) = prad
       case(3)
          values(ncells) = prad/pgas 
       case(4)
          values(ncells) = thisOctal%rho(subcell)
       case(5)
          values(ncells) = thisOctal%ionfrac(subcell,2)
       case default 
          write(*,*) "unknown ivar in store pressure ", ivar
          stop
     end select
  end subroutine storePressure

  subroutine quartiles(temp, n, q1, q2, q3) 
     integer :: n, m
     real(double) :: temp(:)
     real(double) :: q1, q2, q3
     q2 = medianDouble(temp,n)
     ! excludes median value
     m = n/2
     q1 = medianDouble(temp(1:m), m)
     if (mod(n,2) == 0) then
        ! even
        q3 = medianDouble(temp(m+1:n), m)
     else
        ! odd
        q3 = medianDouble(temp(m+2:n), m)
     endif
  end subroutine quartiles
  
  real(double) function medianDouble(temp,n)
     integer :: n
     real(double) :: temp(:) 
     if (mod(n,2) == 0) then
        ! even 
        medianDouble = (temp(n/2) + temp(n/2+1)) / 2.
     else
        ! odd
        medianDouble = temp(n/2+1)
     endif
  end function medianDouble

  recursive subroutine calculateTotalPressure(thisOctal, grid, pgas, prad)
    use ion_mod
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dx, mu, pgas, prad

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateTotalPressure(child, grid, pgas, prad)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
! not splitovermpi
!             if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
                dx = thisoctal%subcellsize * 1.d10 ! cm

                if (thisOctal%ionfrac(subcell, 2) > 0.99d0) then
                   mu = returnMu(thisOctal, subcell, globalIonArray, nGlobalIon)
                   pgas = pgas + thisOctal%rho(subcell)*kerg*thisOctal%temperature(subcell)/(mu*mHydrogen)
                endif

                prad = prad + dx * sqrt(thisOctal%radiationMomentum(subcell)%x**2 + &
                                     thisOctal%radiationMomentum(subcell)%y**2 + &
                                     thisOctal%radiationMomentum(subcell)%z**2)

!             endif
          endif
       endif
    enddo
    
  end subroutine calculateTotalPressure

  subroutine hiiMedianPressure(grid, firstTime, iSource)
    type(GRIDTYPE)    :: grid
    real(double), allocatable :: mean(:), vtot(:), values(:), temp(:)
    real(double), allocatable :: median(:), mins(:), maxs(:), q1(:), q3(:) 
    integer :: ncells, iSource
    integer :: ivar, nvar, n
    character(len=20) :: comp
    character(len=80) :: thisFile
    logical :: firsttime

    nvar = 3
    allocate(vtot(nvar))
    allocate(mean(nvar))
    allocate(mins(nvar))
    allocate(maxs(nvar))
    allocate(median(nvar))
    allocate(q1(nvar))
    allocate(q3(nvar))
    allocate(values(256**3)) 
    vtot = 0.d0; mean = 0.d0; mins = 0.d0; maxs = 0.d0; median = 0.d0; q1 = 0.d0; q3 = 0.d0
    do ivar = 1, nvar 
       ! ivar
       ! (1) pgas
       ! (2) prad
       ! (3) prad/pgas
       ! (4) rho
       ! (5) HII frac

       ncells = 0
       call averagePressure(grid%octreeRoot, grid, ivar, mean(ivar), vtot(ivar), values, ncells)

       ! volume-average
       if (vtot(ivar) > 0.d0) then
          mean(ivar) = mean(ivar) / vtot(ivar)
       else
          write(*,*) "ivar ", ivar, " vtot ", vtot(ivar)
          stop
       endif

       ! median
       n = ncells ! number of cells with non-zero values
       if (n > 3) then
          allocate(temp(1:n))
          temp = values(1:n)
          call sort(n, temp) ! sort in ascending order
          mins(ivar) = temp(1)
          maxs(ivar) = temp(n)
          call quartiles(temp, n, q1(ivar), median(ivar), q3(ivar))
          deallocate(temp)
       endif

    end do

    write(thisFile, '(a)') "hiiregion_median.dat"
    open(55, file=thisFile, status="unknown", form="formatted", position="append")
    if (firstTime) then
       write(55,'(a5,a6,7(a13,1x),a5)') "dump", "var", "volume", &
       "mean","median","q1","q3","min","max","sink"
    endif

    do ivar = 1, nvar
       if (ivar==1) then
          write(comp, '(a)') "pgas"
       elseif (ivar==2) then
          write(comp, '(a)') "prad" 
       elseif (ivar==3) then
          write(comp, '(a)') "pr/pg" 
       elseif (ivar==4) then
          write(comp, '(a)') "rho" 
       elseif (ivar==5) then
          write(comp, '(a)') "hii"
       endif

       write(55,'(i5.4, a6, 7(es13.5,1x), i5)') grid%idump, trim(comp),&
       vtot(ivar)/(pctocm**3),&
       mean(ivar), &
       median(ivar), &
       q1(ivar), &
       q3(ivar), &
       mins(ivar), &
       maxs(ivar), &
       iSource
    enddo
    deallocate(mean,vtot,values,median,mins,maxs,q1,q3)
    close(55)
  end subroutine hiiMedianPressure

  subroutine shellAverage(grid, rVec, fromSource, isource)

    type(GRIDTYPE)    :: grid
    type(OCTAL), pointer :: thisOctal
    type(VECTOR) :: rVec, cen
    real(double) :: dx
    real(double), allocatable :: mean(:,:), vtot(:), rbin(:), values(:,:,:), temp(:)
    real(double), allocatable :: median(:,:), mins(:,:), maxs(:,:), q1(:,:), q3(:,:) 
    integer, allocatable :: ncellsbin(:)
    integer, optional :: isource
    integer :: i, j, nbin, ibin, nvar, n, subcell
    character(len=20) :: weighting, comp
    character(len=80) :: thisFile
    logical :: fromSource

    if (fromSource) then
       write(thisFile, '(a,i4.4,a,i3.3,a)') "rbin_",grid%idump,"_",isource,".dat"
    else
       write(thisFile, '(a,i4.4,a)') "rbin_",grid%idump,"_com.dat"
    endif
    open(55, file=thisFile, status='replace')
    write(55,'(a9,a6,a13,a8,7(a13,1x))') "#weight","var","rbin(pc)","ncells","w*V(pc3)",&
    "mean","median","q1","q3","min","max"

    thisOctal => grid%octreeRoot
    call findSubcellTD(rVec, grid%octreeRoot, thisOctal, subcell) 
    ! define r=0 as the centre of the cell containing rVec
    cen = subcellCentre(thisOctal, subcell) 
    dx = thisoctal%subcellsize * 1.d10 / pctocm ! pc
    if (writeoutput) write(*,*) "vec ", rVec
    if (writeoutput) write(*,*) "cen ", cen

    nbin = 256  ! fixme
    nvar = 5
    allocate(mean(nvar,nbin))
    allocate(vtot(nbin))
    allocate(rbin(nbin))
    allocate(values(int(4*pi*nbin**2),nvar,nbin))
    allocate(ncellsbin(nbin))
    allocate(median(nvar,nbin))
    allocate(mins(nvar,nbin))
    allocate(maxs(nvar,nbin))
    allocate(q1(nvar,nbin))
    allocate(q3(nvar,nbin))
    do i = 1, nbin
       rbin(i) = dble(i) ! units of dx (torus units)
    enddo
    do i = 2,4
       if (i==1) then
          write(weighting, '(a)') "mass" 
       elseif (i==2) then
          write(weighting, '(a)') "none" 
       elseif (i==3) then
          write(weighting, '(a)') "ionNone" 
       elseif (i==4) then
          write(weighting, '(a)') "neuNone" 
       endif
       ! mean(1,ibin) pgas
       ! mean(2,ibin) prad
       ! mean(3,ibin) prad/pgas
       ! mean(4,ibin) rho
       ! mean(5,ibin) HII
       mean(:,:) = 0.d0
       vtot(:) = 0.d0
       ! values(ncellsbin(ibin), 1, ibin) = pgas
       ! values(ncellsbin(ibin), 2, ibin) = prad
       ! values(ncellsbin(ibin), 3, ibin) = prad/pgas
       ! values(ncellsbin(ibin), 4, ibin) = rho
       ! values(ncellsbin(ibin), 5, ibin) = HII
       values(:,:,:) = 0.d0
       ncellsbin(:) = 0
       call averagePressureRadius(grid%octreeRoot, grid, trim(weighting), cen, rbin, nbin, mean, vtot, values, ncellsbin)

       ! volume-average
       do ibin = 1, nbin
          if (vtot(ibin) > 0.d0) then
             mean(:,ibin) = mean(:,ibin) / vtot(ibin)
          endif
       enddo

       ! median
       mins = 0.d0; maxs = 0.d0; median = 0.d0; q1 = 0.d0; q3 = 0.d0
       do ibin = 1, nbin
          n = ncellsbin(ibin) ! number of cells with non-zero values
          if (n > 3) then
             allocate(temp(1:n))
             do j = 1, nvar
                temp = values(1:n, j, ibin)
                call sort(n, temp) ! sort in ascending order
                mins(j,ibin) = temp(1)
                maxs(j,ibin) = temp(n)
                call quartiles(temp, n, q1(j,ibin), median(j,ibin), q3(j,ibin))
             enddo
             deallocate(temp)
          endif
       enddo

       do j = 1, nvar
          if (j==1) then
             write(comp, '(a)') "pgas"
          elseif (j==2) then
             write(comp, '(a)') "prad" 
          elseif (j==3) then
             write(comp, '(a)') "pr/pg" 
          elseif (j==4) then
             write(comp, '(a)') "rho" 
          elseif (j==5) then
             write(comp, '(a)') "hii"
          endif

          do ibin = 1, nbin
             write(55,'(a9,a6,es13.5,i8,7(es13.5,1x))') trim(weighting),trim(comp),dx*rbin(ibin),ncellsbin(ibin),&
             vtot(ibin)/(pctocm**3),&
             mean(j,ibin), &
             median(j,ibin), &
             q1(j,ibin), &
             q3(j,ibin), &
             mins(j,ibin), &
             maxs(j,ibin)
          enddo
       enddo
    enddo
    deallocate(mean,vtot,rbin,values,ncellsbin,median,mins,maxs,q1,q3)
    close(55)
  end subroutine shellAverage

  recursive subroutine calculateTotalInRadius(thisOctal, grid, total, radius, origin, tagged)
    use ion_mod
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    logical, optional :: tagged
    integer :: subcell, i
    type(VECTOR) :: rVec, origin
    real(double) :: dx, mu, tot, radius, dv
    real(double) :: total(:), r, prad, pgas

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateTotalInRadius(child, grid, total, radius, origin, tagged)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
! not splitovermpi
!             if (octalOnThread(thisOctal, subcell, myRankGlobal)) then

                if (present(tagged)) then
                   if (tagged) then
                      ! cells have been tagged - only include cells with etacont > 0 in the totals
                      if (.not.associated(thisOctal%etaCont)) then 
                         cycle
                      elseif (thisOctal%etaCont(subcell) < 0.1d0) then
                         cycle
                      endif
                   endif
                endif

                rVec = subcellCentre(thisOctal, subcell) - origin
                r = modulus(rVec)
                if (r <= radius) then
                   dv = cellVolume(thisOctal, subcell) * 1.d30 ! cm
                   dx = thisoctal%subcellsize * 1.d10 ! cm
                   ! mass (msol)
                   tot = thisOctal%rho(subcell) * dv / msol
                   total(1) = total(1) + tot
                   ! gas thermal pressure (dyn cm-2)
                   mu = returnMu(thisOctal, subcell, globalIonArray, nGlobalIon)
                   pgas = thisOctal%rho(subcell)*kerg*thisOctal%temperature(subcell)/(mu*mHydrogen)
                   total(2) = total(2) + pgas
                   ! radiation pressure (dyn cm-2)
                   prad = dx * sqrt(thisOctal%radiationMomentum(subcell)%x**2 + &
                                   thisOctal%radiationMomentum(subcell)%y**2 + &
                                   thisOctal%radiationMomentum(subcell)%z**2)
                   total(3) = total(3) + prad
                   ! momentum (g cm s-1)
                   tot = dv * sqrt(thisOctal%rhou(subcell)**2 + &
                                   thisOctal%rhov(subcell)**2 + &
                                   thisOctal%rhow(subcell)**2)
                   total(4) = total(4) + tot
                   if (thisOctal%ionfrac(subcell, 1) < 0.1d0 .and. .not.thisOctal%undersampled(subcell)) then
                      ! ionized mass (msol)
                      tot = thisOctal%rho(subcell) * dv / msol
                      total(5) = total(5) + tot
                      ! ionized volume (pc3)
                      tot = dv / (pctocm**3)
                      total(6) = total(6) + tot
                   endif
                   ! no of cells in volume
                   tot = 1.d0 
                   total(10) = total(10) + tot

                   ! sum of squares
                   total(11) = total(11) + pgas**2
                   total(12) = total(12) + prad**2

                endif
!             endif
          endif
       endif
    enddo
    
  end subroutine calculateTotalInRadius

  recursive subroutine calculateMeanVelocity(thisOctal, grid, vmean, vms, mtot, mion, mneu, radius)
    use ion_mod
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    type(VECTOR) :: rVec, rHat, vel
    real(double) :: radius, dv, mass, mion, mneu, vmod, vr
    real(double) :: r, vms(:), vmean(:), mtot, vi(7)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateMeanVelocity(child, grid, vmean, vms, mtot, mion, mneu, radius)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
! not splitovermpi
!             if (octalOnThread(thisOctal, subcell, myRankGlobal)) then

                rVec = subcellCentre(thisOctal, subcell)
                rHat = rVec
                call normalize(rHat)
                r = modulus(rVec)
                if (r <= radius) then
                   dV = cellVolume(thisOctal, subcell) * 1.d30 ! cm
                   mass = thisOctal%rho(subcell) * dV
                   vel = VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell), &
                                thisOctal%rhov(subcell)/thisOctal%rho(subcell), &
                                thisOctal%rhow(subcell)/thisOctal%rho(subcell))
                   vmod = modulus(vel)
                   vr = vel .dot. rhat

                   ! ionized + neutral
                   vi(1) = vmod
                   vi(2) = vr 

                   if (thisOctal%ionfrac(subcell, 2) > 0.98d0) then
                      ! ionized
                      vi(3) = vmod
                      vi(4) = vr
                      ! cell doesn't count towards neutral vmean
                      vi(5) = 0.d0
                      vi(6) = 0.d0
                      mion = mion + mass
                      ! might as well get mean HII temperature while we're here 
                      vi(7) = thisOctal%temperature(subcell)
                   else
                      ! cell doesn't count towards ionized vmean
                      vi(3) = 0.d0
                      vi(4) = 0.d0
                      ! neutral
                      vi(5) = vmod
                      vi(6) = vr
                      mneu = mneu + mass
                      ! cell doesn't count towards ionized temperature 
                      vi(7) = 0.d0 
                   endif

                   vmean(:) = vmean(:) + mass * vi(:)
                   vms(:) = vms(:) + mass * vi(:)**2
                   mtot = mtot + mass 
                endif
!             endif
          endif
       endif
    enddo
  end subroutine calculateMeanVelocity

  recursive subroutine identifyHIIregion(thisOctal, grid, sourcePosition)
    type(GRIDTYPE) :: grid
    type(VECTOR) :: sourcePosition, rVec
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call identifyHIIregion(child, grid, sourcePosition)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
! not splitovermpi
!             if (octalOnThread(thisOctal, subcell, myRankGlobal)) then

                rVec = subcellCentre(thisOctal, subcell)
                call taglineOfSight(grid, rVec, sourcePosition)

!             endif
          endif
       endif
    enddo
  end subroutine identifyHiiRegion

  subroutine tagLineOfSight(grid, startVec, endVec) 
    type(GRIDTYPE) :: grid
    type(VECTOR) :: startVec, endVec, currentPosition, direction
    real(double) :: distToNextCell
    type(OCTAL), pointer :: thisOctal, sOctal
    real(double) :: fudgeFac = 1.d-3
    integer :: subcell
    logical :: endLoop, ionizedPath, blocked


    currentPosition = startVec
    direction = endVec - startVec 
    call normalize(direction)

    CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)

    endLoop = .false.
    ionizedPath = .false.
    blocked = .false.
    do while (inOctal(grid%octreeRoot, currentPosition) .and. .not.endLoop)
       call findSubcellLocal(currentPosition, thisOctal,subcell)
       
       if (thisOctal%ionfrac(subcell, 1) > 0.9d0) then
          ! neutral cell in between origin and source
          blocked = .true.
          endLoop = .true.
       endif

       if (inSubcell(thisOctal, subcell, endVec)) then
          endLoop = .true.
          if (.not. blocked) then
             ! managed to reach the source without hitting neutral block
             ionizedPath = .true.
          endif
       endif

       if (.not. endLoop) then
          sOctal => thisOctal
          call distanceToCellBoundary(grid, currentPosition, direction, DisttoNextCell, sOctal)
          currentPosition = currentPosition + (distToNextCell+fudgeFac*grid%halfSmallestSubcell)*direction
       endif
    end do

    ! reverse direction and tag
    direction = startVec - endVec 
    call normalize(direction)

    endLoop = .false.
    do while (inOctal(grid%octreeRoot, currentPosition) .and. .not.endLoop)
       call findSubcellLocal(currentPosition, thisOctal,subcell)

       if (.not.associated(thisOctal%etaCont)) allocate(thisOctal%etaCont(1:thisOctal%maxChildren))
       if (blocked) then
          thisOctal%etaCont(subcell) = 0.d0
       elseif (ionizedPath) then
          thisOctal%etaCont(subcell) = 1.d0
          ! sanity check
          if (thisOctal%ionfrac(subcell, 1) > 0.9d0) then
             write(*,*) "error - found neutral cell in ionized path"
             stop
          endif
       endif

       if (inSubcell(thisOctal, subcell, startVec)) then
          endLoop = .true.
       endif

       if (.not. endLoop) then
          sOctal => thisOctal
          call distanceToCellBoundary(grid, currentPosition, direction, DisttoNextCell, sOctal)
          currentPosition = currentPosition + (distToNextCell+fudgeFac*grid%halfSmallestSubcell)*direction
       endif
     end do


  end subroutine tagLineOfSight

  recursive subroutine sumPAH(thisOctal, grid, totAbsorption, totEmission, adotZeroEmission)
    use pah_mod, only : PAHemissivityFromAdot
    use atom_mod, only: blambda
    use amr_mod, only : returnKappa
#ifdef MPI
    use amr_mod, only : octalOnThread
#endif
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell, i, j
    real(double) :: totAbsorption, totEmission, emission, lambda, dlam, adotZeroEmission
    real(double) :: kappaAbs, kappaAbsDust, V, absorption

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call sumPAH(child, grid, totAbsorption, totEmission, adotZeroEmission)
                exit
             end if
          end do
       else
!          if(.not. thisoctal%ghostcell(subcell)) then
#ifdef MPI
             if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
#endif
                if (.not. associated(thisOctal%etacont)) allocate(thisOctal%etaCont(1:thisOctal%maxchildren))
                if (thisOctal%adotPAH(subcell) > 0.d0) then
                   ! sum emission spectrum
                   emission = 0.d0
                   do j = 2, grid%nlambda
                      lambda = dble(grid%lamArray(j))
                      dlam = dble(grid%lamArray(j)-grid%lamArray(j-1))
                      ! pah
                      emission = emission +  PAHemissivityFromAdot(lambda, &
                           thisOctal%adotPAH(subcell), &
                           thisOctal%rho(subcell), &
                           thisOctal%dustTypeFraction(subcell,1)) &
                           *cSpeed/(lambda*angstromtocm)**2 * fourPi * 1.d-8 &
                           *dlam
                      ! dust
                      call returnKappa(grid, thisOctal, subcell, lambda=real(lambda), iLambda=j, &
                       kappaAbs=kappaAbs, kappaAbsDust=kappaAbsDust)
                       emission = emission + bLambda(lambda, thisOctal%tdust(subcell)) &
                           *kappaAbsDust * 1.d-10 * fourPi * 1.d-8 & ! conversion from per cm to per A
                           *dlam
                   enddo
                   v = cellVolume(thisOctal, subcell) * 1.d30
                   absorption = thisOctal%adotPAH(subcell) + thisOctal%distanceGrid(subcell)*9.8935326161648629E+032/V

                   ! cell comparison
                   thisOctal%etaCont(subcell) = emission/absorption
                   if (emission == 0.d0) then
!                      write(*,*) "emission 0 for Adot", thisOctal%adotPAH(subcell) &
!                         / (thisOctal%dustTypeFraction(subcell,1)*thisOctal%rho(subcell))
                      adotZeroEmission = adotZeroEmission + thisOctal%adotPAH(subcell)
                   endif

                   ! tots over grid
                   totEmission = totEmission + emission
                   totAbsorption = totAbsorption + absorption
                else
                   thisOctal%etaCont(subcell) = 1.d-15
                endif
#ifdef MPI
             endif
#endif
!          endif
       endif
    enddo
  end subroutine sumPAH

end module gridanalysis_mod
