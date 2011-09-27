module outputs_mod


  use gridio_mod
  implicit none

contains

  subroutine doOutputs(grid)
    use cmf_mod, only : calculateAtomSpectrum
    use modelatom_mod, only : globalAtomArray
    use source_mod, only : globalNSource, globalSourceArray
    use inputs_mod, only : gridOutputFilename, writegrid, calcPhotometry
    use inputs_mod, only : stokesimage
    use inputs_mod, only : calcDataCube, atomicPhysics, nAtom
    use inputs_mod, only : iTransLine, iTransAtom, gridDistance
    use inputs_mod, only : calcImage, calcSpectrum, calcBenchmark
    use inputs_mod, only : photoionPhysics, splitoverMpi, dustPhysics, thisinclination
    use inputs_mod, only : SEDlamMin, SEDlamMax, SEDwavLin, SEDnumLam
    use inputs_mod, only : lambdaImage, mie, gridDistance, nLambda, nv
    use inputs_mod, only : outfile, npixels, nImage, inclinationArray
    use inputs_mod, only : lamStart, lamEnd, lineEmission
    use inputs_mod, only : monteCarloRT
    use image_utils_mod
#ifdef MPI
    use inputs_mod, only : inclineX, inclineY, inclineZ, singleInclination
    use hydrodynamics_mod, only : checkMaclaurinBenchmark
#endif
!    use inputs_mod, only : rotateViewAboutX, rotateViewAboutY, rotateViewAboutZ
    use inputs_mod, only : cmf, lamline, ttauriRouter,amrgridsize
    use physics_mod, only : setupXarray, setupDust
#ifdef MOLECULAR
    use molecular_mod
    use angularImage, only: make_angular_image, map_dI_to_particles
    use inputs_mod, only : molecularPhysics, useDust, realdust, h21cm, internalView, dataCubeFilename
#endif 
    use phasematrix_mod
    use phaseloop_mod, only : do_phaseloop
    use surface_mod, only : surfacetype
    use disc_class, only : alpha_disc
    use blob_mod, only : blobtype
    use setupamr_mod, only : writegridkengo, writeFogel
    use lucy_mod, only : getSublimationRadius
    use inputs_mod, only : fastIntegrate, geometry, intextfilename, outtextfilename
    use formal_solutions, only :compute_obs_line_flux
#ifdef PHOTOION
    use photoion_mod, only: createImagePhotoion
#ifdef MPI
    use photoionAMR_mod, only : createImageSplitGrid
#endif
#endif

    type(BLOBTYPE) :: tblob(1)
    real(double) :: totalFlux, rSub, ang, vFlux, bFlux
    type(SURFACETYPE) :: tsurface
    type(ALPHA_DISC) :: tdisc
    type(GRIDTYPE) :: grid
    type(VECTOR) :: viewVec
#ifdef PHOTOION
    type(VECTOR) :: observerDirection
#endif
    real, pointer :: xArray(:)=>null()
    type(PHASEMATRIX), pointer :: miePhase(:,:,:) => null()
    integer, parameter :: nMuMie = 20
    integer :: i
!    integer :: nAng
!    type(VECTOR) :: thisVec,  axis
!    real(double) :: ang
    type(VECTOR) :: tvec(1)
    character(len=80) :: message
    real(double), allocatable :: flux(:)

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
         (.not.calcBenchmark) ) then
       goto 666
    endif

    if (calcBenchmark) then
       select case (grid%geometry)
#ifdef MPI
          case("maclaurin")
             call checkMaclaurinBenchmark(grid)
#endif
          case DEFAULT
             write(message,'(a)') "No benchmark calculation for: "//trim(grid%geometry)
       end select
    endif
          

    if (atomicPhysics.and.calcDataCube) then
       call setupXarray(grid, xArray, nLambda, atomicDataCube=.true.)
       if (dustPhysics) call setupDust(grid, xArray, nLambda, miePhase, nMumie)
       do i = 1, 2
          viewVec = VECTOR(sin(thisInclination), 0.d0, -cos(thisinclination))
          !       gridDistance = 140.d0* pctocm/1.d10
          ang = pi * dble(i-1)
          viewVec =  rotatez(viewVec, ang)
          !       gridDistance = 140.d0* pctocm/1.d10
          call calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
               viewVec, dble(gridDistance), &
               globalSourceArray, globalnsource, i, totalflux)
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
          globalSourceArray(1)%limbDark(1) = +1.05395E+00
          globalSourceArray(1)%limbDark(2) = -1.64891E-01
          call calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
               viewVec, dble(gridDistance), &
               globalSourceArray, globalnsource, 1, bflux, forceLambda=4400.d0, occultingDisc=.true.)
          globalSourceArray(1)%limbDark(1) = +8.29919E-01
          globalSourceArray(1)%limbDark(2) = +1.62937E-02
          call calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
               viewVec, dble(gridDistance), &
               globalSourceArray, globalnsource, 1, vflux, forceLambda=5500.d0, occultingDisc=.true.)
          if (writeoutput) write(66,'(4f13.4)') ang/twoPi, returnMagnitude(bFlux, "B"), returnMagnitude(vFlux,"V"), &
               returnMagnitude(bFlux, "B") -  returnMagnitude(vFlux,"V")
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


#ifdef PHOTOION
    if (photoionPhysics.and.calcImage) then
       call setupXarray(grid, xArray, nLambda, photoion=.true.)
       if (dustPhysics) call setupDust(grid, xArray, nLambda, miePhase, nMumie)
       if ( splitoverMPI ) then 
#ifdef MPI
!          observerDirection = VECTOR(0.d0, -1.d0, 0.d0)
!More flexible inclinations
          observerDirection = VECTOR(0.d0, 0.d0, 0.d0)
          if(inclineX) observerDirection%x = sin(singleInclination)
          if(inclineY) observerDirection%y = sin(singleInclination)
          if(inclineZ) observerDirection%z = sin(singleInclination)

          do i = 1, nImage
             call createImageSplitGrid(grid, globalnSource, globalsourcearray, observerDirection, i)
          enddo
#else
          call writeFatal("Cannot calculate an image from a domain decomposed grid without MPI")
#endif
       else
          do i = 1, nImage
             observerDirection = VECTOR(-1.0d0, 0.d0, 0.d0)
             call createImagePhotoion(grid, globalnSource, globalsourcearray, observerDirection, i)
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
          call map_dI_to_particles(grid)
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
          call getSublimationRadius(grid, rSub)
          write(message, '(a, f7.3,a )') "Final inner radius is: ",(1.d10*rSub/rSol), " solar radii"
          call writeInfo(message, FORINFO)
          write(message, '(a, f7.3,a )') "Final inner radius is: ",(rSub/globalSourceArray(1)%radius), " core radii"
          call writeInfo(message, FORINFO)


          fastIntegrate=.true.
          call do_phaseloop(grid, .true., 0., 0., 0.,  &
               VECTOR(0., 0., 0.), 0.d0, 0. , 0., 0., 0.d0, &
               tsurface, 0., 0., tdisc, tvec, 1,       &
               0., 0, .false., 100000, &
               miePhase, globalnsource, globalsourcearray, tblob, nmumie, 0.)
       end if

       if (calcImage) then
          do i = 1, nImage
             nlambda = 1
             stokesImage = .true.
             outfile = getImageFilename(i)
             npixels = getImagenPixels(i)
             lamStart = lambdaImage(i)
             lamEnd = lambdaImage(i)
             call setupXarray(grid, xarray, nlambda, lamMin=lambdaImage(i), lamMax=lambdaImage(i), &
                  wavLin=.true., numLam=1, dustRadEq=.true.)

             call setupDust(grid, xArray, nLambda, miePhase, nMumie)
             fastIntegrate=.true.
             call do_phaseloop(grid, .true., 0., 0., 0.,  &
                  VECTOR(0., 0., 0.), 0.d0, 0. , 0., 0., 0.d0, &
                  tsurface, 0., 0., tdisc, tvec, 1,       &
                  0., 0, .false., 100000, &
                  miePhase, globalnsource, globalsourcearray, tblob, nmumie, 0., &
                  overrideInclinations=inclinationArray(i:i) )
          enddo
       endif

    endif

    if (atomicPhysics.and.(calcspectrum.or.calcimage)) then
       mie = .false.
       lineEmission = .true.
       grid%lineEmission = .true.
       if ( calcspectrum ) then 
          write(*,*) "calling setupxarray ",nv
          call setupXarray(grid, xarray, nv, atomicDataCube=.true.)
          write(*,*) "nlambda after setupxarray",nlambda,nv
          nlambda = nv
          call do_phaseloop(grid, .true., 0., 0., 0.,  &
               VECTOR(0., 0., 0.), 0.d0, 0. , 0., 0., 0.d0, &
               tsurface, 0., 0., tdisc, tvec, 1,       &
               0., 0, .true., 100000, &
               miePhase, globalnsource, globalsourcearray, tblob, nmumie, 0.)
       end if

       if (calcImage) then
          do i = 1, nImage
             nlambda = 1
             stokesImage = .true.
             outfile = getImageFilename(i)
             npixels = getImagenPixels(i)
             lamStart = lambdaImage(i)
             lamEnd = lambdaImage(i)
             call setupXarray(grid, xarray, nlambda, lamMin=lambdaImage(i), lamMax=lambdaImage(i), &
                  wavLin=.true., numLam=1, dustRadEq=.true.)

             call do_phaseloop(grid, .true., 0., 0., 0.,  &
                  VECTOR(0., 0., 0.), 0.d0, 0. , 0., 0., 0.d0, &
                  tsurface, 0., 0., tdisc, tvec, 1,       &
                  0., 0, .false., 100000, &
                  miePhase, globalnsource, globalsourcearray, tblob, nmumie, 0., &
                  overrideInclinations=inclinationArray(i:i))
          enddo
       endif

    endif
666 continue
  end subroutine doOutputs


end module outputs_mod
