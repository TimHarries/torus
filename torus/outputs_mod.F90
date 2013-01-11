module outputs_mod


  use gridio_mod
  implicit none

contains

  subroutine doOutputs(grid)
#ifdef CMFATOM
    use cmf_mod, only : calculateAtomSpectrum
    use modelatom_mod, only : globalAtomArray
#endif
    use source_mod, only : globalNSource, globalSourceArray, writeSourceHistory
    use inputs_mod, only : gridOutputFilename, writegrid, calcPhotometry
    use inputs_mod, only : calcDataCube, atomicPhysics, nAtom, sourceHistory
    use inputs_mod, only : iTransLine, iTransAtom, gridDistance
    use inputs_mod, only : calcImage, calcSpectrum, calcBenchmark
    use inputs_mod, only : photoionPhysics, splitoverMpi, dustPhysics, thisinclination
    use inputs_mod, only : mie, gridDistance, nLambda, nv, ncubes
    use inputs_mod, only : lineEmission, postsublimate
    use inputs_mod, only : monteCarloRT, dowriteradialfile, radialfilename
    use inputs_mod, only : sourcelimbaB, sourcelimbbB ,sourcelimbaV, sourcelimbbV
    use sed_mod, only : SEDlamMin, SEDlamMax, SEDwavLin, SEDnumLam
    use image_utils_mod, only: getImageWavelength, getnImage
#ifdef MPI
#ifdef HYDRO
    use hydrodynamics_mod, only : checkMaclaurinBenchmark
#endif
#endif
!    use inputs_mod, only : rotateViewAboutX, rotateViewAboutY, rotateViewAboutZ
    use inputs_mod, only : cmf, lamline, ttauriRouter,amrgridsize
    use physics_mod, only : setupXarray, setupDust
#ifdef MOLECULAR
    use molecular_mod
    use angularImage
    use inputs_mod, only : molecularPhysics, useDust, realdust, h21cm, internalView
    use datacube_mod, only: dataCubeFilename
#endif 
    use phasematrix_mod
    use phaseloop_mod, only : do_phaseloop
    use surface_mod, only : surfacetype
    use disc_class, only : alpha_disc
    use blob_mod, only : blobtype
    use setupamr_mod, only : writegridkengo, writeFogel
    use lucy_mod, only : getSublimationRadius
    use inputs_mod, only : fastIntegrate, geometry, intextfilename, outtextfilename, sourceHistoryFilename, lambdatau
    use formal_solutions, only :compute_obs_line_flux
#ifdef PHOTOION
    use photoion_utils_mod, only: quickSublimate
    use photoion_mod, only: createImagePhotoion
#ifdef MPI
    use photoionAMR_mod, only : createImageSplitGrid
#endif
#endif

!    real(double) :: rSub
    type(GRIDTYPE) :: grid
    real, pointer :: xArray(:)=>null()
    type(PHASEMATRIX), pointer :: miePhase(:,:,:) => null()
    integer, parameter :: nMuMie = 20
    integer :: i
    character(len=80) :: message
    integer :: nimage
    real :: lambdaImage
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

    call writeBanner("Creating outputs","-",TRIVIAL)

    if (writegrid) then
          call writeAMRgrid(gridOutputFilename,.false.,grid)

       if (geometry == "kengo") then
          call writegridkengo(grid)
       endif
       if (geometry == "fogel") then
          call writeFogel(grid, intextFilename, outtextFilename)
       endif

    endif

    if ( (.not.calcImage).and. &
         (.not.calcSpectrum).and. &
         (.not.calcDataCube).and. &
         (.not.calcPhotometry).and. &
         (.not.calcBenchmark) .and. &
         (.not. sourceHistory).and. &
         (.not.dowriteRadialfile) &
         ) then
       goto 666
    endif

    nimage = getnImage()


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
          viewVec = VECTOR(sin(thisInclination), 0.d0, -cos(thisinclination))
          !       gridDistance = 140.d0* pctocm/1.d10
          ang = 2.*pi * dble(i-1)/dble(ncubes)
          viewVec =  rotatez(viewVec, ang)
          !       gridDistance = 140.d0* pctocm/1.d10
          call calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
               viewVec, dble(gridDistance), &
               globalSourceArray, globalnsource, i, totalflux, occultingDisc=.true.)
       enddo
    endif

    if (atomicPhysics.and.calcPhotometry) then
       call setupXarray(grid, xArray, nLambda, atomicDataCube=.true.)
       if (writeoutput) then
          open(66, file="phase.dat", status="unknown", form="formatted")
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
               globalSourceArray, globalnsource, 1, bflux, forceLambda=4400.d0, occultingDisc=.true.)
!          globalSourceArray(1)%limbDark(1) = +8.29919E-01
!          globalSourceArray(1)%limbDark(2) = +1.62937E-02
         globalSourceArray(1)%limbDark(1) = sourcelimbaV
         globalSourceArray(1)%limbDark(2) = sourcelimbbV


          call calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
               viewVec, dble(gridDistance), &
               globalSourceArray, globalnsource, 1, vflux, forceLambda=5500.d0, occultingDisc=.true.)
          if (writeoutput) write(66,'(4f13.4)') ang/twoPi, returnMagnitude(bFlux, "B"), returnMagnitude(vFlux,"V"), &
               returnMagnitude(bFlux, "B") -  returnMagnitude(vFlux,"V")
!          write(*,*) "v, b flux",vflux, bflux
       enddo
       if (writeoutput) close(66)
    endif

    if (atomicPhysics.and.calcSpectrum.and.(.not.monteCarloRT)) then

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
                  globalSourceArray, globalnsource, i, totalflux, occultingDisc=.true.)
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
    if (photoionPhysics.and.calcImage) then
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

    if (dustPhysics.and.(calcspectrum.or.calcimage).and.(.not.photoionPhysics)) then
       mie = .true.
       if ( calcspectrum ) then 
          call setupXarray(grid, xarray, nLambda, lamMin=SEDlamMin, lamMax=SEDlamMax, &
               wavLin=SEDwavLin, numLam=SEDnumLam, dustRadEq=.true.)
          call setupDust(grid, xArray, nLambda, miePhase, nMumie, filestart="sed")
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
          call do_phaseloop(grid, .false., 100000, miePhase, globalnsource, globalsourcearray, nmumie)
       end if

       if (calcImage) then
          do i = 1, nImage
             nlambda = 1
             lambdaImage = getImageWavelength(i)
             lambdatau = lambdaimage
             call setupXarray(grid, xarray, nlambda, lamMin=lambdaImage, lamMax=lambdaImage, &
                  wavLin=.true., numLam=1, dustRadEq=.true.)

             call setupDust(grid, xArray, nLambda, miePhase, nMumie)
#ifdef PHOTOION
             if(postsublimate) then
                call writeInfo("Sublimating dust")
                call quickSublimate(grid%octreeRoot)
             end if
#endif
             fastIntegrate=.true.
             call do_phaseloop(grid, .false., 100000, &
                  miePhase, globalnsource, globalsourcearray, nmumie, imNum=i)
          enddo
       endif

    endif

#ifdef CMFATOM
    if (atomicPhysics.and.(calcspectrum.or.calcimage)) then
       mie = .false.
       lineEmission = .true.
       grid%lineEmission = .true.
       if ( calcspectrum ) then 
          write(*,*) "calling setupxarray ",nv
          call setupXarray(grid, xarray, nv, atomicDataCube=.true.)
          write(*,*) "nlambda after setupxarray",nlambda,nv
          nlambda = nv
          call do_phaseloop(grid, .true., 100000, miePhase, globalnsource, globalsourcearray, nmumie) 
       end if

       if (calcImage) then
          do i = 1, nImage
             nlambda = 1
             lambdaImage = getImageWavelength(i)
             call setupXarray(grid, xarray, nlambda, lamMin=lambdaImage, lamMax=lambdaImage, &
                  wavLin=.true., numLam=1, dustRadEq=.true.)

             call do_phaseloop(grid, .false., 100000, &
                  miePhase, globalnsource, globalsourcearray, nmumie, imNum=i)
          enddo
       endif

    endif
#endif
    if (sourceHistory) then
       call writeSourceHistory(sourceHistoryfilename,globalSourceArray,globalnSource)
    endif

    if (dowriteRadialFile) then
#ifdef MPI
       call writeRadialFile(radialfilename, grid)
#endif
    endif


666 continue
  end subroutine doOutputs


end module outputs_mod
