module outputs_mod


  use gridio_mod
  implicit none

contains

  subroutine doOutputs(grid)
    use cmf_mod, only : calculateAtomSpectrum
    use modelatom_mod, only : globalAtomArray
    use source_mod, only : globalNSource, globalSourceArray
    use input_variables, only : gridOutputFilename, writegrid
    use input_variables, only : useDust, realdust, stokesimage
    use input_variables, only : calcDataCube, atomicPhysics, nAtom
    use input_variables, only : iTransLine, iTransAtom, gridDistance
    use input_variables, only : imageFilename, calcImage, molecularPhysics, calcSpectrum
    use input_variables, only : photoionPhysics, splitoverMpi, dustPhysics, nImage, thisinclination
    use input_variables, only : SEDlamMin, SEDlamMax, SEDwavLin, SEDnumLam
    use input_variables, only : h21cm, internalView
    use input_variables, only : lambdaImage, npixelsArray, dataCubeFilename, mie, gridDistance, nLambda
    use input_variables, only : outfile, npix, ninclination, nImage, inclinations, inclinationArray
    use input_variables, only : lamStart, lamEnd, lineEmission, nVelocity
!    use input_variables, only : rotateViewAboutX, rotateViewAboutY, rotateViewAboutZ
    use physics_mod, only : setupXarray, setupDust
    use molecular_mod
    use phasematrix_mod
    use phaseloop_mod, only : do_phaseloop
    use surface_mod, only : surfacetype
    use disc_class, only : alpha_disc
    use blob_mod, only : blobtype
    use angularImage, only: make_angular_image, map_dI_to_particles
    use input_variables, only : fastIntegrate
#ifdef MPI
    use input_variables, only : outputImageType
    use photoionAMR_mod, only : createImageSplitGrid
#endif

    type(BLOBTYPE) :: tblob(1)
    real(double) :: totalFlux
    type(SURFACETYPE) :: tsurface
    type(ALPHA_DISC) :: tdisc
    type(GRIDTYPE) :: grid
    type(VECTOR) :: viewVec, observerDirection
    real, pointer :: xArray(:)
    type(PHASEMATRIX), pointer :: miePhase(:,:,:) => null()
    integer, parameter :: nMuMie = 20
    integer :: i
!    integer :: nAng
!    type(VECTOR) :: thisVec,  axis
!    real(double) :: ang
    character(len=80) :: tstring
    type(VECTOR) :: tvec(1)

    if ( (.not.calcImage).and. &
         (.not.calcSpectrum).and. &
         (.not.calcDataCube) ) then
       goto 666
    endif

    call writeBanner("Creating outputs","-",TRIVIAL)

    if (writegrid) then
       call writeAMRgrid(gridOutputFilename,.false.,grid)
    endif

    call setupXarray(grid, xArray, nLambda)



    if (atomicPhysics.and.calcDataCube) then
       if (dustPhysics) call setupDust(grid, xArray, nLambda, miePhase, nMumie)
       viewVec = VECTOR(0.d0, sin(thisInclination), -cos(thisinclination))
       gridDistance = 140.d0* pctocm/1.d10
       call calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
            viewVec, dble(gridDistance), &
            globalSourceArray, globalnsource, 1, totalflux, occultingDisc=.true.)
    endif

    if (photoionPhysics.and.splitoverMPI.and.calcImage) then
       observerDirection = VECTOR(0.d0, -1.d0, 0.d0)
#ifdef MPI
       do i = 1, nImage
          call createImageSplitGrid(grid, globalnSource, globalsourcearray, observerDirection, imageFilename(i), lambdaImage(i), &
               outputImageType(i), nPixelsArray(i))
       enddo
#else
    call writeFatal("Cannot calculate an image from a domain decomposed grid without MPI")
#endif
    endif

    if (molecularPhysics.and.calcDataCube) then
       if (dustPhysics) then
          call setupXarray(grid, xarray, nLambda)
          call setupDust(grid, xArray, nLambda, miePhase, nMumie)
          usedust = .true.
          realdust = .true.
       endif
!        call photoionChemistry(grid, grid%octreeRoot)
       call calculateMoleculeSpectrum(grid, globalMolecule, dataCubeFilename)
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

    if (dustPhysics.and.(calcspectrum.or.calcimage)) then
       mie = .true.
       if ( calcspectrum ) then 
          call setupXarray(grid, xarray, nLambda, lamMin=SEDlamMin, lamMax=SEDlamMax, &
               wavLin=SEDwavLin, numLam=SEDnumLam)
          call setupDust(grid, xArray, nLambda, miePhase, nMumie)
          fastIntegrate=.true.
          call do_phaseloop(grid, .true., 0., 0., 0.,  &
               0., 0., VECTOR(0., 0., 0.), 0.d0, 0. , 0., 0., 0.d0, &
               tsurface,  tstring, 0., 0., tdisc, tvec, 1,       &
               0., 0, .false., 100000, &
               miePhase, globalnsource, globalsourcearray, tblob, nmumie, 0.)
       end if
       allocate(inclinations(1))
       if (calcImage) then
          do i = 1, nImage
             nlambda = 1
             stokesImage = .true.
             outfile = imageFilename(i)
             npix = nPixelsArray(i)
             ninclination = 1
             inclinations(1) = inclinationArray(i)
             lamStart = lambdaImage(i)
             lamEnd = lambdaImage(i)
             call setupXarray(grid, xarray, nlambda, lamMin=lambdaImage(i), lamMax=lambdaImage(i), &
                  wavLin=.true., numLam=1)

             call setupDust(grid, xArray, nLambda, miePhase, nMumie)
             fastIntegrate=.true.
             call do_phaseloop(grid, .true., 0., 0., 0.,  &
                  0., 0., VECTOR(0., 0., 0.), 0.d0, 0. , 0., 0., 0.d0, &
                  tsurface,  tstring, 0., 0., tdisc, tvec, 1,       &
                  0., 0, .false., 100000, &
                  miePhase, globalnsource, globalsourcearray, tblob, nmumie, 0.)
          enddo
       endif

    endif

    if (atomicPhysics.and.(calcspectrum.or.calcimage)) then
       mie = .false.
       lineEmission = .true.
       grid%lineEmission = .true.
       if ( calcspectrum ) then 
          call setupXarray(grid, xarray, nVelocity)
          call do_phaseloop(grid, .true., 0., 0., 0.,  &
               0., 0., VECTOR(0., 0., 0.), 0.d0, 0. , 0., 0., 0.d0, &
               tsurface,  tstring, 0., 0., tdisc, tvec, 1,       &
               0., 0, .true., 100000, &
               miePhase, globalnsource, globalsourcearray, tblob, nmumie, 0.)
       end if
       allocate(inclinations(1))
       if (calcImage) then
          do i = 1, nImage
             nlambda = 1
             stokesImage = .true.
             outfile = imageFilename(i)
             npix = nPixelsArray(i)
             ninclination = 1
             inclinations(1) = inclinationArray(i)
             lamStart = lambdaImage(i)
             lamEnd = lambdaImage(i)
             call setupXarray(grid, xarray, nlambda, lamMin=lambdaImage(i), lamMax=lambdaImage(i), &
                  wavLin=.true., numLam=1)

             call setupDust(grid, xArray, nLambda, miePhase, nMumie)
             call do_phaseloop(grid, .true., 0., 0., 0.,  &
                  0., 0., VECTOR(0., 0., 0.), 0.d0, 0. , 0., 0., 0.d0, &
                  tsurface,  tstring, 0., 0., tdisc, tvec, 1,       &
                  0., 0, .false., 100000, &
                  miePhase, globalnsource, globalsourcearray, tblob, nmumie, 0.)
          enddo
       endif

    endif
666 continue
  end subroutine doOutputs

end module outputs_mod
