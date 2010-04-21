module outputs_mod


  use gridio_mod
  implicit none

contains

  subroutine doOutputs(grid)
    use cmf_mod, only : calculateAtomSpectrum
    use modelatom_mod, only : globalAtomArray
    use source_mod, only : globalNSource, globalSourceArray
    use input_variables, only : gridOutputFilename, writegrid, lineimage
    use input_variables, only : useDust, realdust
    use input_variables, only : calcDataCube, atomicPhysics, nAtom
    use input_variables, only : iTransLine, iTransAtom, gridDistance
    use input_variables, only : imageFilename, calcImage, molecularPhysics, calcSpectrum
    use input_variables, only : photoionPhysics, splitoverMpi, dustPhysics, nImage
    use photoionAMR_mod, only : createImageSplitGrid
    use input_variables, only : lambdaImage, outputimagetype, npixelsArray, dataCubeFilename, mie
    use physics_mod, only : setupXarray, setupDust
    use molecular_mod
    use phasematrix_mod
    use phaseloop_mod, only : do_phaseloop
    use surface_mod, only : surfacetype
    use disc_class, only : alpha_disc
    use blob_mod, only : blobtype
    type(BLOBTYPE) :: tblob(1)
    real(double) :: totalFlux
    type(SURFACETYPE) :: tsurface
    type(ALPHA_DISC) :: tdisc
    type(GRIDTYPE) :: grid
    type(VECTOR) :: viewVec, observerDirection
    real, pointer :: xArray(:)
    integer :: nLambda
    type(PHASEMATRIX), pointer :: miePhase(:,:,:) => null()
    integer, parameter :: nMuMie = 20
    integer :: i
    character(len=80) :: tstring
    type(VECTOR) :: tvec(1)

    call writeBanner("Creating outputs","-",TRIVIAL)

    viewVec = VECTOR(-1.d0, 0.d0, 0.d0)
    if (writegrid) then
       call writeAMRgrid(gridOutputFilename,.false.,grid)
    endif

    call setupXarray(grid, xArray, nLambda)


    if (dustPhysics) call setupDust(grid, xArray, nLambda, miePhase, nMumie)

    if (atomicPhysics.and.calcDataCube) then
       gridDistance = 140.d0* pctocm/1.d10
       call calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
            viewVec, dble(gridDistance), &
            globalSourceArray, globalnsource, 1, totalflux)
    endif

    if (photoionPhysics.and.splitoverMPI.and.calcImage) then
       
       observerDirection = VECTOR(0.d0, -1.d0, 0.d0)
       do i = 1, nImage
          call createImageSplitGrid(grid, globalnSource, globalsourcearray, observerDirection, imageFilename(i), lambdaImage(i), &
               outputImageType(i), nPixelsArray(i))
       enddo
    endif

    if (molecularPhysics.and.calcDataCube) then
       lineimage =  .true.
       if (dustPhysics) then
          call setupXarray(grid, xarray, nLambda)
          call setupDust(grid, xArray, nLambda, miePhase, nMumie)
          usedust = .true.
          realdust = .true.
       endif
!        call photoionChemistry(grid, grid%octreeRoot)
       call calculateMoleculeSpectrum(grid, globalMolecule, dataCubeFilename)
    endif

    if (dustPhysics.and.(calcspectrum.or.calcimage)) then
       mie = .true.
       call setupXarray(grid, xarray, nLambda)
       call setupDust(grid, xArray, nLambda, miePhase, nMumie)
       
       call do_phaseloop(grid, .true., 0., 0., 0.,  &
            0., 0., VECTOR(0., 0., 0.), 0.d0, 0. , 0., 0., 0.d0, &
            tsurface,  tstring, 0., 0., tdisc, tvec, 1,       &
            0., 0, .false., 0., 100000, &
            miePhase, globalnsource, globalsourcearray, tblob, nmumie, 0.)
    endif

  end subroutine doOutputs

end module outputs_mod
