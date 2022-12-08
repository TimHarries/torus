module outputs_mod


  use gridio_mod
  implicit none

contains

SUBROUTINE doOutputs(grid)
  USE gas_opacity_mod
#ifdef CMFATOM
  USE cmf_mod, ONLY : calculateAtomSpectrum
  USE modelatom_mod, ONLY : globalAtomArray
#endif
    use source_mod, only : globalNSource, globalSourceArray, writeSourceHistory
    use inputs_mod, only : gridOutputFilename, writegrid, calcPhotometry, amr2d
    use inputs_mod, only : calcDataCube, atomicPhysics, nAtom, sourceHistory, calcDustCube, doAnalysis, doClusterAnalysis
    use inputs_mod, only : iTransLine, iTransAtom, gridDistance, gasOpacityPhysics
    use inputs_mod, only : writePolar
    use inputs_mod, only : calcImage, calcSpectrum, calcBenchmark, calcMovie, calcColumnImage
    use inputs_mod, only : photoionPhysics, splitoverMpi, dustPhysics, thisinclination
    use inputs_mod, only : mie, gridDistance, nLambda, ncubes
    use inputs_mod, only : postsublimate, lineEmission, nv
    use inputs_mod, only : dowriteradialfile, radialfilename
    use inputs_mod, only : scatteringSurface, ScatteringSurfaceFilename
    use inputs_mod, only : sourcelimbaB, sourcelimbbB ,sourcelimbaV, sourcelimbbV
    use sed_mod, only : SEDlamMin, SEDlamMax, SEDwavLin, SEDnumLam
    use image_mod
    use image_utils_mod, only: getImageWavelength, getnImage, getimageFilename, getImagenPixelsX, &
         getImageNPixelsY,getImageSizeX, &
         getImageSizeY, getFluxUnits
    use vtk_mod, only : writeVtkFilenBody
#ifdef MPI
#ifdef HYDRO
  USE hydrodynamics_mod, ONLY : checkMaclaurinBenchmark
#endif
#endif
  !    use inputs_mod, only : rotateViewAboutX, rotateViewAboutY, rotateViewAboutZ
  USE inputs_mod, ONLY : cmf, lamline, ttauriRouter,amrgridsize
  USE h21cm_mod, ONLY: h21cm
  USE physics_mod, ONLY : setupXarray, setupDust
  USE gridanalysis_mod
#ifdef MOLECULAR
  USE molecular_mod
  USE angularImage
  USE angularImage_utils, ONLY: internalView
  USE inputs_mod, ONLY : molecularPhysics, useDust, realdust
  USE datacube_mod, ONLY: dataCubeFilename
#endif
  USE phasematrix_mod
  USE phaseloop_mod, ONLY : do_phaseloop
  USE surface_mod, ONLY : surfacetype
  USE disc_class, ONLY : alpha_disc
  USE blob_mod, ONLY : blobtype
  USE dust_mod, ONLY : readLambdaFile, dumpPolarizability
  USE setupamr_mod, ONLY : writegridkengo, writeFogel
  USE lucy_mod, ONLY : getSublimationRadius
  USE inputs_mod, ONLY : fastIntegrate, geometry, intextfilename, outtextfilename, sourceHistoryFilename, lambdatau, itrans
  USE inputs_mod, ONLY : lambdaFilename, polarWavelength, polarFilename, nPhotSpec, nPhotImage, nPhotons
  USE inputs_mod, ONLY : nDataCubeInclinations, datacubeInclinations, nLamLine, lamLineArray !, rgapinner1
  USE formal_solutions, ONLY :compute_obs_line_flux
#ifdef PHOTOION
  USE photoion_utils_mod, ONLY: quickSublimate
  USE photoion_mod, ONLY: createImagePhotoion
#ifdef MPI
  USE inputs_mod, ONLY : columnimagedirection, imodel, columnImageFilename
  USE photoionAMR_mod, ONLY : createImageSplitGrid
  USE mpi_global_mod, ONLY : loadBalancingThreadGlobal
#endif
#endif

  !    real(double) :: rSub
  TYPE(GRIDTYPE) :: grid
  TYPE(DATACUBE) :: thisCube
  REAL, POINTER :: xArray(:)=>null()
  CHARACTER(len=80) :: tempChar, tempFilename
  TYPE(PHASEMATRIX), POINTER :: miePhase(:,:,:) => NULL()
  INTEGER, PARAMETER :: nMuMie = 180
  INTEGER :: i, j
  CHARACTER(len=80) :: message
  REAL(DOUBLE) :: lambdaArray(2000), dx
  INTEGER :: nimage, nCubeLambda
  TYPE(IMAGETYPE) :: imageSlice

#ifdef MPI
  REAL(DOUBLE), POINTER :: image(:,:)
  CHARACTER(len=80) :: thisFile
#endif
  REAL :: lambdaImage
  REAL, ALLOCATABLE :: tarray(:,:)
  REAL(DOUBLE), ALLOCATABLE :: xArrayDouble(:)
  REAL(DOUBLE) :: kabs, ksca, oldMass, oldAge
  INTEGER :: iInc
  LOGICAL, SAVE :: firstSource = .TRUE.
#ifdef MOLECULAR
  !    integer :: nAng
  !    type(VECTOR) :: thisVec,  axis
  !    real(double) :: ang
#endif
#ifdef CMFATOM
  REAL(DOUBLE) :: totalFlux, ang, vFlux, bFlux
  TYPE(ALPHA_DISC) :: tdisc
  TYPE(VECTOR) :: viewVec
  REAL(DOUBLE), ALLOCATABLE :: flux(:)
#endif
  IF (FirstSource) THEN
     oldMAss = 0.d0
     oldAge = 0.d0
     firstSource = .FALSE.
  ENDIF

  CALL writeBanner("Creating outputs","-",TRIVIAL)

  IF (writegrid.AND.(.NOT.loadBalancingThreadGlobal)) THEN
     CALL writeAMRgrid(gridOutputFilename,.FALSE.,grid)

     IF (geometry == "kengo") THEN
        CALL writegridkengo(grid)
     ENDIF
     IF (geometry == "fogel") THEN
        CALL writeFogel(grid, intextFilename, outtextFilename)
     ENDIF


  ENDIF

    if (scatteringSurface) then
       call writeInfo("Writing scattering surface file")
       call writeScatteringSurfaceFile(scatteringSurfaceFilename, grid)
    endif


    if (geometry == "envelope") then
       if (writeoutput.and.amr2d) call writeEduard(grid)
       goto 666
    endif

  IF (geometry == "mgascii") THEN
     CALL writeVtkFilenBody(globalnSource, globalsourceArray, "nbody_rosette.vtk")
     CALL writeVtkFile(grid, "dump_rosette.vtk", &
          valueTypeString=(/"rho          ","HI           ","temperature  ", &
          "tdust        ", "scatters     ","ioncross     ","tempconv     ","crossings    "/))
  ENDIF


  IF (dustPhysics.AND.writepolar) THEN
     CALL setupXarray(grid, xarray, nLambda, numLam = 1)
     CALL setupDust(grid, xArray, nLambda, miePhase, nMumie)
     CALL dumpPolarizability(miePhase, nMuMie, xarray, nLambda)
  ENDIF


  IF ( (.NOT.calcImage).AND. &
       (.NOT.calcSpectrum).AND. &
       (.NOT.calcColumnImage).AND. &
       (.NOT.calcMovie).AND. &
       (.NOT.doAnalysis).AND. &
       (.NOT.doClusterAnalysis).AND. &
       (.NOT.calcDataCube).AND. &
       (.NOT.calcPhotometry).AND. &
       (.NOT.calcBenchmark) .AND. &
       (.NOT. sourceHistory).AND. &
       (.NOT.dowriteRadialfile) &
       ) THEN
     GOTO 666
  ENDIF

  nimage = getnImage()

  IF (calcColumnImage) THEN
#ifdef MPI
     IF (.NOT.loadBalancingThreadGlobal) THEN
        CALL createColumnDensityImage(grid, columnImageDirection,image)
        CALL findmultifilename(TRIM(columnImageFilename), iModel, thisFile)
#ifdef USECFITSIO
        IF (writeoutput) CALL writeFitsColumnDensityImage(image, TRIM(thisFile))
#else
        CALL writeWarning("Torus was build without FITS support. No image written.")
#endif
     ENDIF
#endif
  ENDIF

  IF (doClusterAnalysis) THEN
#ifdef PHOTOION
     IF (photoionPhysics.OR.dustPhysics) THEN
        CALL setupXarray(grid, xArray, nLambda, photoion=.TRUE.)
        IF (dustPhysics) CALL setupDust(grid, xArray, nLambda, miePhase, nMumie)
     ENDIF
#endif
#ifdef MPI
     CALL clusterAnalysis(grid, globalSourceArray, globalnSource, nLambda, xArray, miePhase, nMuMie)
#else
       call clusterAnalysis(grid, globalSourceArray, globalnSource)
#endif
  ENDIF

  !    if (calcColumnImage) then
#ifdef MPI
  !       if (.not.loadBalancingThreadGlobal) then
  !          call createColumnDensityImage(grid, columnImageDirection,image)
  !          call findmultifilename(columnImageFilename, iModel, thisFile)
  !          if (writeoutput) call writeFitsColumnDensityImage(image, thisFile)
  !       endif
#endif
  !    endif

  IF (doAnalysis) CALL analysis(grid)
  !    if (doAnalysis) call calculateToomreQ(grid%octreeRoot, grid)
  !    if (doAnalysis) call writeVtkFile(grid, "toomreq.vtk",  valueTypeString=(/"etacont"/))

  IF (calcBenchmark) THEN
     SELECT CASE (grid%geometry)
#ifdef MPI
#ifdef HYDRO
     CASE("maclaurin")
        CALL checkMaclaurinBenchmark(grid)
#endif
#endif
     CASE DEFAULT
        WRITE(message,'(a)') "No benchmark calculation for: "//TRIM(grid%geometry)
     END SELECT
  ENDIF



#ifdef CMFATOM
  IF (atomicPhysics.AND.calcDataCube) THEN
     CALL setupXarray(grid, xArray, nLambda, atomicDataCube=.TRUE.)
     IF (dustPhysics) CALL setupDust(grid, xArray, nLambda, miePhase, nMumie)
     DO i = 1, ncubes
        DO ilambda = 1, nLamLine
           DO iInc = 1, nDataCubeInclinations

              viewVec = VECTOR(SIN(dataCubeInclinations(iInc)), 0.d0, -COS(datacubeInclinations(iInc)))
              ang = 2.*pi * DBLE(i-1)/DBLE(ncubes)
              viewVec =  rotatez(viewVec, ang)

              lamLine = lamLineArray(iLambda)

              tempChar = TRIM(dataCubeFilename)
              IF (INDEX(datacubefilename,".fits")/=0) THEN
                 tempChar = datacubeFilename(1:(INDEX(datacubefilename,".fits")-1))
              ENDIF
              IF (ncubes==1) THEN
                 WRITE(tempFilename,'(a,a,i3.3,a,i5.5)') TRIM(tempChar),"_", &
                      NINT(radtodeg*dataCubeInclinations(iInc)),"_", &
                      NINT(lamLine)
              ELSE
                 WRITE(tempFilename,'(a,a,i3.3,a,i5.5,a,i3.3)') TRIM(tempChar),"_",&
                      NINT(radtodeg*dataCubeInclinations(iInc)),"_", &
                      NINT(lamLine),"_",i
              ENDIF
              IF (writeoutput) WRITE(*,*) "Calculating spectrum: ",TRIM(tempFilename)
              CALL calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
                   viewVec, DBLE(gridDistance), &
                   globalSourceArray, globalnsource, totalflux, occultingDisc=.FALSE., prefix=tempFilename)
           ENDDO
        ENDDO
     ENDDO
  ENDIF

  IF (atomicPhysics.AND.calcPhotometry) THEN
     CALL setupXarray(grid, xArray, nLambda, atomicDataCube=.TRUE.)
     IF (myrankGlobal == 1) THEN
        OPEN(28, file="phase.dat", status="unknown", form="formatted")
     ENDIF
     DO i = 1, 50
        viewVec = VECTOR(SIN(thisInclination), 0.d0, -COS(thisinclination))
        !       gridDistance = 140.d0* pctocm/1.d10
        ang = twoPi * DBLE(i-1)/50.d0
        viewVec =  rotatez(viewVec, ang)

        !          globalSourceArray(1)%limbDark(1) = +1.05395E+00
        !          globalSourceArray(1)%limbDark(2) = -1.64891E-01
        globalSourceArray(1)%limbDark(1) = sourcelimbaB
        globalSourceArray(1)%limbDark(2) = sourcelimbbB
        CALL calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
             viewVec, DBLE(gridDistance), &
             globalSourceArray, globalnsource,  bflux, " ",forceLambda=4400.d0, occultingDisc=.TRUE.)
        !          globalSourceArray(1)%limbDark(1) = +8.29919E-01
        !          globalSourceArray(1)%limbDark(2) = +1.62937E-02
        globalSourceArray(1)%limbDark(1) = sourcelimbaV
        globalSourceArray(1)%limbDark(2) = sourcelimbbV


        CALL calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
             viewVec, DBLE(gridDistance), &
             globalSourceArray, globalnsource, vflux, " ",forceLambda=5500.d0, occultingDisc=.TRUE.)
        IF (myrankGlobal == 1) WRITE(28,'(4f13.4)') ang/twoPi, returnMagnitude(bFlux, "B"), returnMagnitude(vFlux,"V"), &
             returnMagnitude(bFlux, "B") -  returnMagnitude(vFlux,"V")
        !          write(*,*)  myrankGlobal,writeoutput, " v, b flux",vflux, bflux
     ENDDO
     IF (myrankGlobal==1) CLOSE(28)
  ENDIF

  IF (.FALSE.) THEN
     !    if (atomicPhysics.and.calcSpectrum.and.(.not.monteCarloRT)) then

     IF (cmf) THEN
        CALL setupXarray(grid, xArray, nLambda, atomicDataCube=.TRUE.)
        DO i = 1, 50
           viewVec = VECTOR(SIN(thisInclination), 0.d0, -COS(thisinclination))
           !       gridDistance = 140.d0* pctocm/1.d10
           ang = twoPi * DBLE(i-1)/50.d0
           viewVec =  rotatez(viewVec, ang)
           globalSourceArray(1)%limbDark(1) = 0.d0
           globalSourceArray(1)%limbDark(2) = 0.d0
           CALL calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
                viewVec, DBLE(gridDistance), &
                globalSourceArray, globalnsource, totalflux, " ",occultingDisc=.TRUE.)
        ENDDO
     ELSE
        CALL setupXarray(grid, xArray, nLambda, atomicDataCube=.TRUE.)
        ALLOCATE(flux(1:nLambda))
        viewVec = VECTOR(SIN(thisInclination), 0.d0, -COS(thisinclination))
        !       gridDistance = 140.d0* pctocm/1.d10
        ang = twoPi * DBLE(i-1)/50.d0
        viewVec =  rotatez(viewVec, ang)
        CALL compute_obs_line_flux(lamline, REAL(mHydrogen), DBLE(globalSourceArray(1)%radius), &
             DBLE(TTauriRouter/1.0e10), DBLE(amrGridSize)/2.001d0/SQRT(2.0d0), &
             globalsourcearray(1)%surface, &
             100, 20, 50, 100,  &
             (-1.d0)*viewVec, DBLE(griddistance), grid, 2., .TRUE.,  &
             flux, grid%lamArray, grid%nlambda, "flux.dat", &
             .TRUE., .FALSE., &
             .FALSE., .FALSE., 100, &
             .FALSE., tdisc, .FALSE., .FALSE., &
             .FALSE.)
     ENDIF

     GOTO 666
  ENDIF
#endif


#ifdef PHOTOION
  IF (photoionPhysics.AND.(calcImage.OR.calcMovie)) THEN
     CALL setupXarray(grid, xArray, nLambda, photoion=.TRUE.)
     IF (dustPhysics) CALL setupDust(grid, xArray, nLambda, miePhase, nMumie)
     IF(postsublimate) THEN
        CALL writeInfo("Sublimating dust")
        CALL quickSublimate(grid%octreeRoot)
     END IF

     IF ( splitoverMPI ) THEN
#ifdef MPI
        DO i = 1, nImage
           CALL createImageSplitGrid(grid, globalnSource, globalsourcearray, i)
        ENDDO
#else
        CALL writeFatal("Cannot calculate an image from a domain decomposed grid without MPI")
#endif
     ELSE
        DO i = 1, nImage
           CALL createImagePhotoion(grid, globalnSource, globalsourcearray, i)
        END DO
     END IF
  ENDIF
#endif

#ifdef MOLECULAR
  IF (molecularPhysics.AND.calcDataCube) THEN
     IF (dustPhysics) THEN
        CALL setupXarray(grid, xarray, nLambda, dustRadEq = .TRUE.)
        CALL setupDust(grid, xArray, nLambda, miePhase, nMumie)
        usedust = .TRUE.
        realdust = .TRUE.
     ENDIF
     !        call photoionChemistry(grid, grid%octreeRoot)
     IF ( internalView ) THEN
        CALL make_angular_image(grid)
     ELSE
        CALL calculateMoleculeSpectrum(grid, globalMolecule, dataCubeFilename)
     END IF
     !       viewVec = VECTOR(-0.5d0, 0.5d0, 1.d0/sqrt(2.d0))
     !       nAng = 240
     !       axis = VECTOR(0.d0, 0.d0, 1.d0)
     !!       axis = rotateX(axis, -rotateViewAboutX * degtorad)
     !!       axis = rotateY(axis, -rotateViewAboutY * degtorad)
     !       axis = rotateZ(axis, -rotateViewAboutZ * degtorad)
     !
     !       do i = 187, nAng
     !          ang = fourPi * dble(i-1)/dble(nAng-1)
     !          thisVec = arbitraryRotate(viewVec, ang, axis)
     !          write(dataCubeFilename,'(a,i4.4,a)') "CLU",i+2000,".fits"
     !          call calculateMoleculeSpectrum(grid, globalMolecule, dataCubeFilename, inputViewVec=thisVec)
     !       enddo

  ENDIF

  IF (h21cm .AND. calcDataCube) THEN
     IF ( internalView ) THEN
        CALL make_angular_image(grid)
#ifdef SPH
        CALL map_dI_to_particles(grid)
#endif
     ELSE
        CALL make_h21cm_image(grid)
     END IF
  END IF
#endif

  IF (dustPhysics.AND.(calcspectrum.OR.calcimage.OR.calcMovie).AND.(.NOT.photoionPhysics)) THEN
     mie = .TRUE.
     IF ( calcspectrum ) THEN
        CALL setupXarray(grid, xarray, nLambda, lamMin=SEDlamMin, lamMax=SEDlamMax, &
             wavLin=SEDwavLin, numLam=SEDnumLam, dustRadEq=.TRUE.)
        CALL setupDust(grid, xArray, nLambda, miePhase, nMumie, filestart="sed")
        IF (gasOpacityPhysics) THEN
           ALLOCATE(xArrayDouble(1:nLambda))
           xArrayDouble = DBLE(xArray)
           CALL createAllMolecularTables(nLambda, xArrayDouble, reset=.TRUE.)
           DEALLOCATE(xArrayDouble)
           IF (writeoutput) THEN
              DO j = 1, nLambda
                 CALL returnGasKappaValue(grid, 2000., 1.d0, kappaAbs=kabs, kappaSca=ksca, lambda=xarray(j), ilambda=j)
                 WRITE(*,'(1p,5e13.4)') xArray(j),kabs+ksca,kabs,ksca, ksca/(ksca+kabs)
              ENDDO
           ENDIF
        ENDIF

#ifdef PHOTOION
        IF(postsublimate) THEN
           CALL writeInfo("Sublimating dust")
           CALL quickSublimate(grid%octreeRoot)
        END IF
#endif
!          call getSublimationRadius(grid, rSub)
!          write(message, '(a, f7.3,a )') "Final inner radius is: ",(1.d10*rSub/rSol), " solar radii"
!          call writeInfo(message, FORINFO)
!          write(message, '(a, f7.3,a )') "Final inner radius is: ",(rSub/globalSourceArray(1)%radius), " core radii"
!          call writeInfo(message, FORINFO)


          fastIntegrate=.true.
          nPhotons = nPhotSpec
          call do_phaseloop(grid, .false., 50000, miePhase, globalnsource, globalsourcearray, nmumie)
       end if

       if ((calcImage.or.calcMovie).and.(.not.calcDustCube)) then
          do i = 1, nImage
             nlambda = 1
             lambdaImage = getImageWavelength(i)
             lambdatau = lambdaimage
             call setupXarray(grid, xarray, nlambda, lamMin=lambdaImage, lamMax=lambdaImage, &
                  wavLin=.true., numLam=1, dustRadEq=.true.)

             call setupDust(grid, xArray, nLambda, miePhase, nMumie)
             polarWavelength = lambdaimage
             write(polarFilename,'(a,i6.6,a)') "polar",nint(polarwavelength/1d4),".dat"
             call dumpPolarizability(miePhase, nMuMie, xarray, nLambda)
             if (gasOpacityPhysics) then
                allocate(xArrayDouble(1:nLambda))
                xArrayDouble = dble(xArray)
                call createAllMolecularTables(nLambda, xArrayDouble)
                deallocate(xArrayDouble)
                if (writeoutput) then
                   do j = 1, nLambda
                      call returnGasKappaValue(grid, 2000., 1.d0, kappaAbs=kabs, kappaSca=ksca, lambda=xarray(j), ilambda=j)
                      write(*,'(1p,5e13.4)') xArray(j),kabs+ksca,kabs,ksca, ksca/(ksca+kabs)
                   enddo
                endif
             endif


#ifdef PHOTOION
           IF(postsublimate) THEN
              CALL writeInfo("Sublimating dust")
              CALL quickSublimate(grid%octreeRoot)
           END IF
#endif
             nPhotons = nPhotImage
             fastIntegrate=.true.


!             call setDustInsideRadiusTozero(grid%octreeRoot, dble(rgapinner1))
!             write(*,*) "setting inner dust to zero"


                call do_phaseloop(grid, .false., 50000, &
                     miePhase, globalnsource, globalsourcearray, nmumie, imNum=i)
             enddo
          endif

       if (calcImage.and.calcDustCube) then
          call readLambdaFile(lambdaFilename, lambdaArray, nCubeLambda)

          do i = 1, nImage

             npixels = getImagenPixelsY(i)
             call setCubeParams(getImagenPixelsX(i), 1., -1.)
             call initCube(thisCube, nv=nCubeLambda)
             call addWavelengthAxis(thiscube, lambdaArray(1:nCubeLambda))
             call addSpatialAxes(thisCube, dble(-getImageSizeY(i)/2.e10), dble(getImageSizeY(i)/2.e10), &
                  griddistance, 1.d-30)
             do j = 1, nCubeLambda
                imageSlice = initImage(i)
                nlambda = 1
                lambdaImage = real(lambdaArray(j))
                lambdatau = lambdaimage
                call setupXarray(grid, xarray, nlambda, lamMin=lambdaImage, lamMax=lambdaImage, &
                     wavLin=.true., numLam=1, dustRadEq=.true.)

                call setupDust(grid, xArray, nLambda, miePhase, nMumie)
                fastIntegrate=.true.
                call do_phaseloop(grid, .false., 50000, &
                     miePhase, globalnsource, globalsourcearray, nmumie, imNum=i, returnImage=imageSlice)
                dx = imageSlice%xAxisCentre(2) - imageSlice%xAxisCentre(1)
                allocate(tarray(1:thiscube%nx,1:thiscube%ny))
                tarray = real(imageSlice%pixel(:,:)%i)
                call ConvertArrayToMJanskiesPerStr(tarray, lambdaImage, dx, dble(gridDistance))
                imageSlice%pixel(:,:)%i = tarray(:,:)
                deallocate(tarray)
                call addImageSliceToCube(thisCube, imageSlice, j)
                call freeImage(imageSlice)
             enddo
#ifdef USECFITSIO
           !             call convertSpatialAxes(thisCube, "au")
           IF (myrankGlobal == 0) CALL writeDataCube(thisCube, getImageFilename(i))
           CALL freeDataCube(thisCube)
#else
           CALL writeWarning("Torus was build without FITS support. No data cube written.")
#endif
          enddo
       endif

    endif

!#ifdef CMFATOM
    if (atomicPhysics.and.(calcspectrum.or.calcimage.or.calcMovie)) then
       mie = .false.
       lineEmission = .true.
       grid%lineEmission = .true.
       if ( calcspectrum ) then
          write(*,*) "calling setupxarray ",nv
          call setupXarray(grid, xarray, nv, phaseloop=.true.)
          write(*,*) "nlambda after setupxarray",nlambda,nv
          nlambda = nv
          call do_phaseloop(grid, .true., 50000, miePhase, globalnsource, globalsourcearray, nmumie)
       end if

       if (calcImage.or.calcMovie) then
          if (dustPhysics) then
             do i = 1, nImage
                nlambda = 1
                lambdaImage = getImageWavelength(i)
                call setupXarray(grid, xarray, nlambda, lamMin=lambdaImage, lamMax=lambdaImage, &
                     wavLin=.true., numLam=1, dustRadEq=.true.)

                call do_phaseloop(grid, .false., 50000, &
                     miePhase, globalnsource, globalsourcearray, nmumie, imNum=i)
             enddo
          else
             nv = 200
             write(*,*) "calling setupxarray for image ",nv
             call setupXarray(grid, xarray, nv, phaseloop=.true.)
             write(*,*) "nlambda after setupxarray",nlambda,nv
             nlambda = nv
             do i = 1, nImage
                call do_phaseloop(grid, .true., 50000, miePhase, globalnsource, globalsourcearray, nmumie,imnum=i)
             enddo
          endif


       endif

    endif
!#endif
    if (sourceHistory) then
       call writeSourceHistory(sourceHistoryfilename,globalSourceArray,globalnSource, oldMass, oldAge)
    endif




    if (dowriteRadialFile) then
       if (splitOverMpi) then
#ifdef MPI
        CALL writeRadialFile(radialfilename, grid)
#endif
     ELSE
#ifdef MOLECULAR
        IF (molecularPhysics) THEN
           CALL writeRadialMolecular(grid, radialFilename, globalMolecule, itrans)
        ELSE
#endif
          call writeValuesNoDomain(radialFilename, grid)
!          call writeRadial(grid, radialFilename)
#ifdef MOLECULAR
        ENDIF
#endif
     ENDIF
  ENDIF


666 CONTINUE
END SUBROUTINE doOutputs


  subroutine writeRadial(grid, filename)
    use inputs_mod, only : spherical
    type(GRIDTYPE) :: grid
    integer :: i
    character(len=*) :: filename
    integer :: nr, nOctals
    real(double), allocatable :: rArray(:), rhoArray(:), tArray(:)

    call countVoxels(grid%octreeRoot, nOctals, nr)
    allocate(rArray(1:nr), tArray(1:nr), rhoArray(1:nr))
    nr = 0
    call getRadial(grid%octreeRoot,  nr, rArray, rhoArray, tArray)

    if (writeoutput) then
       open(33, file=filename, status="unknown", form="formatted")
       if (spherical) then
          write(33,'(a)') "# radius (AU), dust temperature (K), density N(H_2)"
       else
          write(33,'(a)') "# x (AU), dust temperature (K), density N(H_2)"
       endif
       do i = 1, nr
          write(33,'(1p,e13.3,e13.3,e13.3)') rArray(i)*1.d10/autocm,tArray(i),rhoArray(i)/(2.33d0*mHydrogen)
       enddo
       close(33)
    endif
    deallocate(rArray, tArray, rhoArray)
  end subroutine writeRadial



  recursive subroutine getRadial(thisOctal, nr, rArray, rhoArray, tArray)
  real(double) :: rArray(:), rhoArray(:), tArray(:)
  integer :: nr
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child
  integer :: subcell, i

  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call  getRadial(child, nr, rArray, rhoArray, tArray)
                exit
             end if
          end do

       else
          nr = nr + 1
          rArray(nr) = modulus(subcellCentre(thisOctal, subcell))
          tArray(nr) = dble(thisOctal%temperature(subcell))
          rhoArray(nr) = thisOctal%rho(subcell)
       endif
    enddo
  end subroutine getRadial

  recursive subroutine setDustInsideRadiusTozero(thisOctal, radius)
  real(double) :: radius
  type(VECTOR) :: centre
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child
  integer :: subcell, i

  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call  setDustInsideRadiusToZero(child, radius)
                exit
             end if
          end do

       else
          centre = subcellCentre(thisOctal, subcell)
          if ((centre%x**2 + centre%y**2) < radius**2) then
             thisOctal%dustTypeFraction(subcell,:) = 1.d-30
          endif
       endif
    enddo
  end subroutine setDustInsideRadiusTozero

  subroutine writeEduard(grid)
    type(GRIDTYPE) :: grid
    integer :: i,j
    type(VECTOR) :: rVec
    integer, parameter :: nPoints=250
    real(double) :: rArray(nPoints),zArray(nPoints), zDash, rDash
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    open(33, file="envelope.txt", status="old", form="formatted")
    do i = 1, nPoints
       read(33,*) rArray(i), zArray(i)
    enddo
    close(33)

    open(33, file="rzt.dat", status="unknown", form="formatted")
    write(33,'(a)') "% R (au), Z (au), T(K)"
    thisOctal => grid%octreeRoot
    do i = 1, nPoints
       do j = 1 , 10
          zDash = (zArray(i) * dble(j-1)/9.d0)*autocm/1.d10
          rDash = rArray(i) * autocm/1.d10

          rVec = VECTOR(rDash, 0.d0, zDash)

          call findSubcellLocal(rVec, thisOctal, subcell)

          write(33, '(1p,3e15.3)') rArray(i), zDash*1.d10/autocm, thisOctal%temperature(subcell)
       enddo
    enddo
    close(33)
  end subroutine writeEduard

  subroutine writeValuesNoDomain(outputFilename, grid)
    type(GRIDTYPE) :: grid
    real(double) ::  tVal
    character(len=*) :: outputFilename
    type(VECTOR) :: rVec, uHat, cVec
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    open(21, file=outputFilename, status="unknown", form="formatted")
    rVec = VECTOR(0.d0, 0.d0, 0.d0)
    uHat = VECTOR(1.d0, 0.d0, 0.d0)
    thisOctal => grid%octreeRoot
    call findSubcellLocal(rVec, thisOctal, subcell)
    do while(inOctal(grid%octreeRoot, rVec))
       call findSubcellLocal(rVec, thisOctal, subcell)
       cVec = subcellCentre(thisOctal, subcell)
       write(21, *) modulus(cVec)*1.d10/autocm, thisOctal%temperature(subcell)
       call distanceToCellBoundary(grid, rVec, uHat, tVal, sOctal=thisOctal)
       rVec = rVec + (tVal + 1.d-3*grid%halfSmallestSubcell) * uHat
    enddo
    close (21)
  end subroutine writeValuesNoDomain

end module outputs_mod
