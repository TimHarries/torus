module outputs_mod


  use gridio_mod
  implicit none

contains

  subroutine doOutputs(grid)
    use gas_opacity_mod
#ifdef CMFATOM
    use cmf_mod, only : calculateAtomSpectrum
    use modelatom_mod, only : globalAtomArray
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
    use inputs_mod, only : sourcelimbaB, sourcelimbbB ,sourcelimbaV, sourcelimbbV
    use sed_mod, only : SEDlamMin, SEDlamMax, SEDwavLin, SEDnumLam
    use image_mod
    use image_utils_mod, only: getImageWavelength, getnImage, getimageFilename, getImagenPixelsX, &
         getImageNPixelsY,getImageSizeX, &
         getImageSizeY, getFluxUnits
    use vtk_mod, only : writeVtkFilenBody
#ifdef MPI
#ifdef HYDRO
    use hydrodynamics_mod, only : checkMaclaurinBenchmark
#endif
#endif
!    use inputs_mod, only : rotateViewAboutX, rotateViewAboutY, rotateViewAboutZ
    use inputs_mod, only : cmf, lamline, ttauriRouter,amrgridsize
    use h21cm_mod, only: h21cm
    use physics_mod, only : setupXarray, setupDust
    use gridanalysis_mod
#ifdef MOLECULAR
    use molecular_mod
    use angularImage
    use angularImage_utils, only: internalView
    use inputs_mod, only : molecularPhysics, useDust, realdust
    use datacube_mod, only: dataCubeFilename
#endif 
    use phasematrix_mod
    use phaseloop_mod, only : do_phaseloop
    use surface_mod, only : surfacetype
    use disc_class, only : alpha_disc
    use blob_mod, only : blobtype
    use dust_mod, only : readLambdaFile, dumpPolarizability
    use setupamr_mod, only : writegridkengo, writeFogel
    use lucy_mod, only : getSublimationRadius
    use inputs_mod, only : fastIntegrate, geometry, intextfilename, outtextfilename, sourceHistoryFilename, lambdatau, itrans
    use inputs_mod, only : lambdaFilename, polarWavelength, polarFilename, nPhotSpec, nPhotImage, nPhotons
    use inputs_mod, only : nDataCubeInclinations, datacubeInclinations, nLamLine, lamLineArray !, rgapinner1
    use formal_solutions, only :compute_obs_line_flux
#ifdef PHOTOION
    use photoion_utils_mod, only: quickSublimate
    use photoion_mod, only: createImagePhotoion
#ifdef MPI
    use inputs_mod, only : columnimagedirection, imodel, columnImageFilename
    use photoionAMR_mod, only : createImageSplitGrid
    use mpi_global_mod, only : loadBalancingThreadGlobal
#endif
#endif

!    real(double) :: rSub
    type(GRIDTYPE) :: grid
    type(DATACUBE) :: thisCube
    real, pointer :: xArray(:)=>null()
    character(len=80) :: tempChar, tempFilename
    type(PHASEMATRIX), pointer :: miePhase(:,:,:) => null()
    integer, parameter :: nMuMie = 180
    integer :: i, j
    character(len=80) :: message
    real(double) :: lambdaArray(2000), dx
    integer :: nimage, nCubeLambda
    type(IMAGETYPE) :: imageSlice

#ifdef MPI
    real(double), pointer :: image(:,:)
    character(len=80) :: thisFile
#endif
    real :: lambdaImage
    real, allocatable :: tarray(:,:)
    real(double), allocatable :: xArrayDouble(:)
    real(double) :: kabs, ksca, oldMass, oldAge
    integer :: iInc
    logical, save :: firstSource = .true.
#ifdef MOLECULAR
!    integer :: nAng
!    type(VECTOR) :: thisVec,  axis
!    real(double) :: ang
#endif
#ifdef CMFATOM
    real(double) :: totalFlux, ang, vFlux, bFlux
    type(ALPHA_DISC) :: tdisc
    type(VECTOR) :: viewVec
    real(double), allocatable :: flux(:)    
#endif
    if (FirstSource) then
       oldMAss = 0.d0
       oldAge = 0.d0
       firstSource = .false.
    endif

    call writeBanner("Creating outputs","-",TRIVIAL)

    if (writegrid.and.(.not.loadBalancingThreadGlobal)) then
          call writeAMRgrid(gridOutputFilename,.false.,grid)

       if (geometry == "kengo") then
          call writegridkengo(grid)
       endif
       if (geometry == "fogel") then
          call writeFogel(grid, intextFilename, outtextFilename)
       endif


    endif

    if (geometry == "envelope") then
       if (writeoutput.and.amr2d) call writeEduard(grid)
       goto 666
    endif

    if (geometry == "mgascii") then
       call writeVtkFilenBody(globalnSource, globalsourceArray, "nbody_rosette.vtk")
       call writeVtkFile(grid, "dump_rosette.vtk", &
         valueTypeString=(/"rho          ","HI           ","temperature  ", &
          "tdust        ", "scatters     ","ioncross     ","tempconv     ","crossings    "/))
    endif


  if (dustPhysics.and.writepolar) then
     call setupXarray(grid, xarray, nLambda, numLam = 1)
     call setupDust(grid, xArray, nLambda, miePhase, nMumie)
     call dumpPolarizability(miePhase, nMuMie, xarray, nLambda)
  endif


    if ( (.not.calcImage).and. &
         (.not.calcSpectrum).and. &
         (.not.calcColumnImage).and. &
         (.not.calcMovie).and. &
         (.not.doAnalysis).and. &
         (.not.doClusterAnalysis).and. &
         (.not.calcDataCube).and. &
         (.not.calcPhotometry).and. &
         (.not.calcBenchmark) .and. &
         (.not. sourceHistory).and. &
         (.not.dowriteRadialfile) &
         ) then
       goto 666
    endif

    nimage = getnImage()

    if (calcColumnImage) then
#ifdef MPI
       if (.not.loadBalancingThreadGlobal) then
          call createColumnDensityImage(grid, columnImageDirection,image)
          call findmultifilename(trim(columnImageFilename), iModel, thisFile)
#ifdef USECFITSIO
          if (writeoutput) call writeFitsColumnDensityImage(image, trim(thisFile))
#else
          call writeWarning("Torus was build without FITS support. No image written.")
#endif
       endif
#endif
    endif

    if (doClusterAnalysis) then
#ifdef PHOTOION
       if (photoionPhysics.or.dustPhysics) then
          call setupXarray(grid, xArray, nLambda, photoion=.true.)
          if (dustPhysics) call setupDust(grid, xArray, nLambda, miePhase, nMumie)
       endif
#endif
#ifdef MPI
       call clusterAnalysis(grid, globalSourceArray, globalnSource, nLambda, xArray, miePhase, nMuMie)
#else
       call clusterAnalysis(grid, globalSourceArray, globalnSource)
#endif
    endif

!    if (calcColumnImage) then
#ifdef MPI
!       if (.not.loadBalancingThreadGlobal) then
!          call createColumnDensityImage(grid, columnImageDirection,image)
!          call findmultifilename(columnImageFilename, iModel, thisFile)
!          if (writeoutput) call writeFitsColumnDensityImage(image, thisFile)
!       endif
#endif
!    endif

    if (doAnalysis) call analysis(grid)
!    if (doAnalysis) call calculateToomreQ(grid%octreeRoot, grid)
!    if (doAnalysis) call writeVtkFile(grid, "toomreq.vtk",  valueTypeString=(/"etacont"/))

    if (calcBenchmark) then
       select case (grid%geometry)
#ifdef MPI
#ifdef HYDRO
          case("maclaurin")
             call checkMaclaurinBenchmark(grid)
#endif
#endif
          case DEFAULT
             write(message,'(a)') "No benchmark calculation for: "//trim(grid%geometry)
       end select
    endif


          
#ifdef CMFATOM
    if (atomicPhysics.and.calcDataCube) then
       call setupXarray(grid, xArray, nLambda, atomicDataCube=.true.)
       if (dustPhysics) call setupDust(grid, xArray, nLambda, miePhase, nMumie)
       do i = 1, ncubes
          do ilambda = 1, nLamLine
             do iInc = 1, nDataCubeInclinations

                viewVec = VECTOR(sin(dataCubeInclinations(iInc)), 0.d0, -cos(datacubeInclinations(iInc)))
                ang = 2.*pi * dble(i-1)/dble(ncubes)
                viewVec =  rotatez(viewVec, ang)

                lamLine = lamLineArray(iLambda)
 
                tempChar = trim(dataCubeFilename)
                if (index(datacubefilename,".fits")/=0) then
                   tempChar = datacubeFilename(1:(index(datacubefilename,".fits")-1))
                endif
                if (ncubes==1) then
                   write(tempFilename,'(a,a,i3.3,a,i5.5)') trim(tempChar),"_", &
                        nint(radtodeg*dataCubeInclinations(iInc)),"_", &
                        nint(lamLine)
                else
                   write(tempFilename,'(a,a,i3.3,a,i5.5,a,i3.3)') trim(tempChar),"_",& 
                        nint(radtodeg*dataCubeInclinations(iInc)),"_", &
                        nint(lamLine),"_",i
                endif
                if (writeoutput) write(*,*) "Calculating spectrum: ",trim(tempFilename)
                call calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
                     viewVec, dble(gridDistance), &
                     globalSourceArray, globalnsource, totalflux, occultingDisc=.false., prefix=tempFilename)
             enddo
          enddo
       enddo
    endif

    if (atomicPhysics.and.calcPhotometry) then
       call setupXarray(grid, xArray, nLambda, atomicDataCube=.true.)
       if (myrankGlobal == 1) then
          open(28, file="phase.dat", status="unknown", form="formatted")
       endif
       do i = 1, 50
          viewVec = VECTOR(sin(thisInclination), 0.d0, -cos(thisinclination))
          !       gridDistance = 140.d0* pctocm/1.d10
          ang = twoPi * dble(i-1)/50.d0
          viewVec =  rotatez(viewVec, ang)

!          globalSourceArray(1)%limbDark(1) = +1.05395E+00
!          globalSourceArray(1)%limbDark(2) = -1.64891E-01
          globalSourceArray(1)%limbDark(1) = sourcelimbaB
          globalSourceArray(1)%limbDark(2) = sourcelimbbB 
         call calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
               viewVec, dble(gridDistance), &
               globalSourceArray, globalnsource,  bflux, " ",forceLambda=4400.d0, occultingDisc=.true.)
!          globalSourceArray(1)%limbDark(1) = +8.29919E-01
!          globalSourceArray(1)%limbDark(2) = +1.62937E-02
         globalSourceArray(1)%limbDark(1) = sourcelimbaV
         globalSourceArray(1)%limbDark(2) = sourcelimbbV


          call calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
               viewVec, dble(gridDistance), &
               globalSourceArray, globalnsource, vflux, " ",forceLambda=5500.d0, occultingDisc=.true.)
          if (myrankGlobal == 1) write(28,'(4f13.4)') ang/twoPi, returnMagnitude(bFlux, "B"), returnMagnitude(vFlux,"V"), &
               returnMagnitude(bFlux, "B") -  returnMagnitude(vFlux,"V")
!          write(*,*)  myrankGlobal,writeoutput, " v, b flux",vflux, bflux
       enddo
       if (myrankGlobal==1) close(28)
    endif

if (.false.) then
!    if (atomicPhysics.and.calcSpectrum.and.(.not.monteCarloRT)) then

       if (cmf) then
          call setupXarray(grid, xArray, nLambda, atomicDataCube=.true.)
          do i = 1, 50
             viewVec = VECTOR(sin(thisInclination), 0.d0, -cos(thisinclination))
             !       gridDistance = 140.d0* pctocm/1.d10
             ang = twoPi * dble(i-1)/50.d0
             viewVec =  rotatez(viewVec, ang)
             globalSourceArray(1)%limbDark(1) = 0.d0
             globalSourceArray(1)%limbDark(2) = 0.d0
             call calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
                  viewVec, dble(gridDistance), &
                  globalSourceArray, globalnsource, totalflux, " ",occultingDisc=.true.)
          enddo
       else
          call setupXarray(grid, xArray, nLambda, atomicDataCube=.true.)
          allocate(flux(1:nLambda))
          viewVec = VECTOR(sin(thisInclination), 0.d0, -cos(thisinclination))
          !       gridDistance = 140.d0* pctocm/1.d10
          ang = twoPi * dble(i-1)/50.d0
          viewVec =  rotatez(viewVec, ang)
          call compute_obs_line_flux(lamline, REAL(mHydrogen), DBLE(globalSourceArray(1)%radius), &
               dble(TTauriRouter/1.0e10), dble(amrGridSize)/2.001d0/SQRT(2.0d0), &
               globalsourcearray(1)%surface, &
               100, 20, 50, 100,  &
               (-1.d0)*viewVec, dble(griddistance), grid, 2., .true.,  &
               flux, grid%lamArray, grid%nlambda, "flux.dat", &
               .true., .false., &
               .false., .false., 100, &
               .false., tdisc, .false., .false., &
               .false.) 
       endif

       goto 666
    endif
#endif


#ifdef PHOTOION
    if (photoionPhysics.and.(calcImage.or.calcMovie)) then
       call setupXarray(grid, xArray, nLambda, photoion=.true.)
       if (dustPhysics) call setupDust(grid, xArray, nLambda, miePhase, nMumie)
       if(postsublimate) then
          call writeInfo("Sublimating dust")
          call quickSublimate(grid%octreeRoot)
       end if
       
       if ( splitoverMPI ) then 
#ifdef MPI
          do i = 1, nImage
             call createImageSplitGrid(grid, globalnSource, globalsourcearray, i)
          enddo
#else
          call writeFatal("Cannot calculate an image from a domain decomposed grid without MPI")
#endif
       else
          do i = 1, nImage
             call createImagePhotoion(grid, globalnSource, globalsourcearray, i)
          end do
       end if
    endif
#endif

#ifdef MOLECULAR
    if (molecularPhysics.and.calcDataCube) then
       if (dustPhysics) then
          call setupXarray(grid, xarray, nLambda, dustRadEq = .true.)
          call setupDust(grid, xArray, nLambda, miePhase, nMumie)
          usedust = .true.
          realdust = .true.
       endif
!        call photoionChemistry(grid, grid%octreeRoot)
       if ( internalView ) then 
          call make_angular_image(grid)
       else
          call calculateMoleculeSpectrum(grid, globalMolecule, dataCubeFilename)
       end if
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

    endif

    if (h21cm .and. calcDataCube) then
       if ( internalView ) then 
          call make_angular_image(grid)
#ifdef SPH
          call map_dI_to_particles(grid)
#endif
       else
          call make_h21cm_image(grid)
       end if
    end if
#endif

    if (dustPhysics.and.(calcspectrum.or.calcimage.or.calcMovie).and.(.not.photoionPhysics)) then
       mie = .true.
       if ( calcspectrum ) then 
          call setupXarray(grid, xarray, nLambda, lamMin=SEDlamMin, lamMax=SEDlamMax, &
               wavLin=SEDwavLin, numLam=SEDnumLam, dustRadEq=.true.)
          call setupDust(grid, xArray, nLambda, miePhase, nMumie, filestart="sed")
          if (gasOpacityPhysics) then
             allocate(xArrayDouble(1:nLambda))
             xArrayDouble = dble(xArray)
             call createAllMolecularTables(nLambda, xArrayDouble, reset=.true.)
             deallocate(xArrayDouble)
             if (writeoutput) then
                do j = 1, nLambda
                   call returnGasKappaValue(grid, 2000., 1.d0, kappaAbs=kabs, kappaSca=ksca, lambda=xarray(j), ilambda=j)
                   write(*,'(1p,5e13.4)') xArray(j),kabs+ksca,kabs,ksca, ksca/(ksca+kabs)
                enddo
             endif
          endif

#ifdef PHOTOION
          if(postsublimate) then
             call writeInfo("Sublimating dust")
             call quickSublimate(grid%octreeRoot)
          end if
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
             if(postsublimate) then
                call writeInfo("Sublimating dust")
                call quickSublimate(grid%octreeRoot)
             end if
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
             if (myrankGlobal == 0) call writeDataCube(thisCube, getImageFilename(i))
             call freeDataCube(thisCube)
#else
             call writeWarning("Torus was build without FITS support. No data cube written.")
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
       call writeRadialFile(radialfilename, grid)
#endif
    else
#ifdef MOLECULAR
       if (molecularPhysics) then
          call writeRadialMolecular(grid, radialFilename, globalMolecule, itrans)
       else
#endif
          call writeRadial(grid, radialFilename)
#ifdef MOLECULAR
       endif
#endif
    endif
    endif


666 continue
  end subroutine doOutputs

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


end module outputs_mod
