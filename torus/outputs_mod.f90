module outputs_mod


  use gridio_mod

contains

  subroutine doOutputs(grid)
    use cmf_mod, only : calculateAtomSpectrum
    use modelatom_mod, only : globalAtomArray
    use source_mod, only : globalNSource, globalSourceArray
    use input_variables, only : gridOutputFilename, writegrid
    use input_variables, only : calcDataCube, atomicPhysics, nAtom
    use input_variables, only : iTransLine, iTransAtom, gridDistance
    real(double) :: totalFlux
    type(GRIDTYPE) :: grid
    type(VECTOR) :: viewVec
    
    viewVec = VECTOR(-1.d0, 0.d0, 0.d0)
    if (writegrid) then
       call writeAMRgrid(gridOutputFilename,.false.,grid)
    endif

    if (atomicPhysics.and.calcDataCube) then
       gridDistance = 140.d0* pctocm/1.d10
       call calculateAtomSpectrum(grid, globalAtomArray, nAtom, iTransAtom, iTransLine, &
            viewVec, dble(gridDistance), &
            globalSourceArray, globalnsource, 1, totalflux)
    endif


  end subroutine doOutputs

end module outputs_mod
