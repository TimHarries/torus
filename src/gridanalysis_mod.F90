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
    call findDiscMass(grid%octreeRoot, globalSourceArray(1)%mass, mDisc, 1500.d0*autocm/1.d10)
    write(*,*) "disc mass ",mdisc/msol
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
  
  subroutine clusterAnalysis(grid, source, nSource, nLambda, lamArray, miePhase, nMuMie)
    use inputs_mod, only : imodel, edgeRadius, splitOverMPI, burstTime, amrgridsize
    use inputs_mod, only : columnImageDirection, columnImageFilename, findNUndersampled, nClusterIonLoops, findHabing
    use utils_mod, only : findMultifilename
    use phasematrix_mod
#ifdef MPI
    use mpi
#ifdef PHOTOION
    use photoionAMR_mod, only : createTemperatureColumnImage, writeLuminosities, calculateAverageTemperature, photoionizationLoopAMR
    use photoionAMR_mod, only : createTdustColumnImage, calculateAverageTdust, createEmissionMeasureImage
    use mpi_amr_mod, only : createIonizationImage, createColumnDensityImage, createhiimage, findNumberUndersampled
#endif
    use hydrodynamics_mod, only : setupevenuparray
#endif
#ifdef CFITSIO
    use image_mod, only : writeFitsColumnDensityImage 
#endif
    use inputs_mod, only : calculateEmissionMeasure, calculateLymanFlux
    use inputs_mod, only : smallestCellSize, writeLums, plotAvgTemp, calculateGlobalAvgTemp, plotAvgTdust, calculateGlobalAvgTdust
    integer :: nMuMie
    type(PHASEMATRIX) :: miePhase(:,:,:)
    integer :: nLambda
    real :: lamArray(:)
    type(GRIDTYPE) :: grid
    type(SOURCETYPE) :: source(:)
    integer :: nSource, i
!    type(VECTOR) :: startPoint, endPoint, firststartpoint, thisDir
!    real(double) :: maxRho, totalMass 
    character(len=80) :: thisFile!, thisFileGrid, thisFileRadius!, rootFilename, fm
#ifdef CFITSIO
    real(double), pointer :: image(:,:), rmsImage(:,:)
#endif
    logical, save :: firstTime=.true.
!    real(double), save :: radius 
!    real(double) :: mass, highestRho, volume, ionizedMass, ionizedVolume, totalKE, totalTE, massflux
    character(len=20) :: weighting
    integer :: ierr
    real(double) :: weightedFluxInRadius, massInRadius, meang0, g0inCell, distance
    real(double) :: n0, t0, trms, nrms, sigma, sigmane, total, tempDouble
    integer :: nSampled, nUndersampled, nCells, nPhotoIter
    real(double) :: loopLimitTime
    real :: iterTime(3)
    integer :: evenUpArray(nHydroThreadsGlobal), iterStack(3), optID, iSource
    real(double) :: nlyA, nlyB, nlyR

#ifdef PHOTOION
#ifdef MPI
    if (nClusterIonLoops > 0) then
       nPhotoIter = 1
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
             write(weighting, '(a)') "hiidensity" 
          elseif (i==5) then
             write(weighting, '(a)') "ne2" 
          elseif (i==6) then
             write(weighting, '(a)') "none" 
          elseif (i==7) then
             write(weighting, '(a)') "ionNone" 
          endif
#ifdef CFITSIO
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
             call MPI_ALLREDUCE(sigma, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             sigma = tempDouble ! sum(T w dV)
             call MPI_ALLREDUCE(total, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             total = tempDouble ! sum(w dV)
             call MPI_ALLREDUCE(sigmaNe, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             sigmaNe = tempDouble ! sum(ne w dV)
             t0 = sigma/total ! mean temperature
             n0 = sigmaNe/total ! mean electron density

             ! rms 
             sigma = 0.d0
             total = 0.d0
             sigmaNe = 0.d0
             call calculateAverageTemperature(grid%octreeRoot, grid, sigma, total, trim(weighting), .true., t0, sigmaNe, n0)
             call MPI_ALLREDUCE(sigma, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             sigma = tempDouble ! sum( (T-T0)^2 w dV)
             call MPI_ALLREDUCE(total, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             total = tempDouble ! sum( w dV)
             call MPI_ALLREDUCE(sigmaNe, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroPlusAMRCommunicator, ierr)
             sigmaNe = tempDouble ! sum( (ne-n0)^2 w dV)
             ! rms = sqrt(<(T-T0)^2>)
             trms = sqrt(sigma/total)
             nrms = sqrt(sigmaNe/total)
             if (myrankglobal == 1) then
                write(thisFile, '(a)') "avgtempGlobalRms.dat"
                open(55, file=thisFile, status="unknown", position="append")
                if (firstTime) then
                   write(55,'(a5,a20,a13,4(a12,1x))') "#dump", "t(s)", "weight", "T0(K)", "Trms(K)", "ne_0(cm-3)", "neRms(cm-3)"
                endif
                write(55,'(i5.4,f20.1, a13, 4(es12.5,1x))') grid%idump, grid%currenttime, trim(weighting), t0, trms, n0, nrms
                close(55)
             endif
          endif
       enddo
    endif

    if (plotAvgTdust .or. calculateGlobalAvgTdust) then
       do i=1,2
          if (i==1) then
             write(weighting, '(a)') "dustmass"
          elseif (i==2) then
             write(weighting, '(a)') "none" 
          endif
#ifdef CFITSIO
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
#ifdef CFITSIO
    if (calculateEmissionMeasure) then
!       call resetNe(grid%octreeRoot)
       call createEmissionMeasureImage(grid, columnImageDirection, image)
       write(thisFile, '(a,i4.4,a)') "emissionMeasure_", grid%idump, ".fits"
       if (writeoutput) call writeFitsColumnDensityImage(image, thisFile)
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
          write(thisFile, '(a)') "averageHabingFlux_allStars.dat"
          ! write header
          if (firstTime) then
             open(69, file=thisFile, status="unknown", form="formatted", position="append")
                write(69, '(a6,a20,a5,a9,3(1x,a12))') "#dump", "t(s)", "i*", "M*(Msol)", "dToM*1(pc)", "G0inCell(H)", "avgG0(H)"
             close(69)
          endif
          open(69, file=thisFile, status="old", position="append", form="formatted")
          do iSource = 1, globalnSource
             weightedFluxInRadius = 0.d0
             massInRadius = 0.d0
             call calculateAverageHabingFlux(grid%octreeRoot, globalSourceArray(iSource), 2.5d0*smallestCellSize, &
                     weightedFluxInRadius, massInRadius, g0inCell)
             meanG0 = weightedFluxInRadius/massInRadius
             ! distance to most massive star
             if (iSource > 1) then
                distance = modulus(globalSourceArray(iSource)%position - globalSourceArray(1)%position)*1.d10/pctocm
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

    
    
#ifdef CFITSIO
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

end module gridanalysis_mod
