module outputs_mod


  use gridio_mod
  implicit none

contains

  subroutine doOutputs(grid)
    use cmf_mod, only : calculateAtomSpectrum
    use modelatom_mod, only : globalAtomArray
    use source_mod, only : globalNSource, globalSourceArray
    use input_variables, only : gridOutputFilename, writegrid, lineimage
    use input_variables, only : calcDataCube, atomicPhysics, nAtom
    use input_variables, only : iTransLine, iTransAtom, gridDistance
    use input_variables, only : imageFilename, calcImage, molecularPhysics
    use input_variables, only : photoionPhysics, splitoverMpi, dustPhysics
    use photoionAMR_mod, only : createImageSplitGrid
    use physics_mod, only : setupXarray, setupDust
    use molecular_mod
    use phasematrix_mod

    real(double) :: totalFlux
    type(GRIDTYPE) :: grid
    type(VECTOR) :: viewVec, observerDirection
    real, pointer :: xArray(:)
    integer :: nLambda
    type(PHASEMATRIX), pointer :: miePhase(:,:,:) => null()
    integer, parameter :: nMuMie = 20

    
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
       call createImageSplitGrid(grid, globalnSource, globalsourcearray, observerDirection, imageFilename)
    endif

    if (molecularPhysics.and.calcDataCube) then
       lineimage =  .true.
       call calculateMoleculeSpectrum(grid, globalMolecule)
    endif

  end subroutine doOutputs

end module outputs_mod
