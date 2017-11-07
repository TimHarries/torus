module clump_mod
  ! this is the equivalent to blob_mod, but restricted to the T Tauri
  !   AMR geometry. nhs

  implicit none 
  public
  
  type clumpType
    real :: startTime
    real :: duration
    real :: mDot
    real :: azimuth
    real :: angularSize 
    logical :: activeFlag = .false.
  end type clumpType
  ! if you change this, remember to update the read/write routines
  !   in grid_mod.
  
  type(clumpType), save, dimension(:), allocatable :: clumps

end module clump_mod
