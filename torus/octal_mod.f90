MODULE octal_mod
  ! data type and some routines needed by for adaptive mesh refinement. nhs.
  
  ! these are in a separate module from amr_mod because they need to be
  !   accessed from module(s) that are themselves USEd by amr_mod.

  USE kind_mod
  USE vector_mod
!  USE linked_list_class
  
  
  IMPLICIT NONE

  public :: subcellCentre, within_subcell

!       y                  z
!       |                 /
!       |      __________/______ 
!       |     /        /       /|
!       |    /   7    /   8   / | 
!       |   /________/_______/  |
!       |  /        /       /| 8|   Diagram showing the convention used here for
!       | /   3    /   4   / |  |   numbering the subcells of each octal.
!       |/________/_______/  |  |
!       |        |        | 4| /|
!       |        |        |  |/ |
!       |    3   |   4    |  /  |
!       |        |        | /| 6|
!       |________|________|/ |  /      
!       |        |        |  | /
!       |        |        | 2|/
!       |    1   |   2    |  /
!       |        |        | /
!       |________|________|/________\ x
!                                   /

  TYPE octalWrapper
    TYPE(octal), POINTER  :: content => NULL()
    LOGICAL(KIND=logicKind), DIMENSION(8) :: inUse
  END TYPE octalWrapper
 
  TYPE wrapperArray
    TYPE(octalWrapper), DIMENSION(:), POINTER :: wrappers   ! a number of octal wrappers  
  END TYPE wrapperArray

  TYPE octal

    INTEGER                            :: nDepth       ! depth of octal. root is 1, it's childen are 2...
    INTEGER                            :: nChildren    ! how many pointers to children there are (max 8)
    INTEGER                            :: indexChild(8)! index of child array containing
                                                       !   pointer to each subcell's child (if it exists) 
    TYPE(octal), DIMENSION(:), POINTER :: child
    LOGICAL, DIMENSION(8)              :: hasChild
    TYPE(octal), POINTER               :: parent          
    TYPE(octalVector)                  :: centre

    REAL, DIMENSION(8)                 :: rho            ! density
    TYPE(vector), DIMENSION(8)         :: velocity       ! velocity
    TYPE(vector), DIMENSION(27)        :: cornerVelocity ! velocity at corners of subcells
    REAL, DIMENSION(8)                 :: temperature    ! grid subcell temperatures
    REAL, DIMENSION(8)                 :: distanceGrid   ! distance crossing used by lucy R Eq
    INTEGER, DIMENSION(8)              :: nCrossings     ! no of photon crossings used by lucy R Eq
    REAL, DIMENSION(:,:), POINTER      :: kappaAbs       ! cont absorption opacities
    REAL, DIMENSION(:,:), POINTER      :: kappaSca       ! scattering opacities
    REAL, DIMENSION(8)                 :: chiLine        ! line opacity
    REAL, DIMENSION(8)                 :: etaLine        ! line emissivity
    REAL, DIMENSION(8)                 :: etaCont        ! line emissivity
    REAL, DIMENSION(8)                 :: biasLine3D     ! grid bias distrubtion
    REAL, DIMENSION(8)                 :: biasCont3D     ! grid bias distrubtion
    REAL(KIND=doubleKind), DIMENSION(8) :: probDistLine  ! emissivity probabilty distribution
    REAL(KIND=doubleKind), DIMENSION(8) :: probDistCont  ! emissivity probabilty distribution
    REAL(KIND=doubleKind), DIMENSION(:,:), POINTER ::  N ! stateq level pops
    REAL(KIND=doubleKind), DIMENSION(8) :: Ne            ! electron density
    REAL(KIND=doubleKind), DIMENSION(8) :: nTot          ! total density
    LOGICAL(KIND=logicKind), DIMENSION(8) :: inStar      ! point lies within star
    LOGICAL(KIND=logicKind), DIMENSION(8) :: inFlow      ! inside accretion flow region GET RID OF THIS?    
    INTEGER, DIMENSION(8) :: label                       ! numeric label for each subcell. 
      ! the subcell labels may be useful for debugging the code, but are not needed for
      !   any of the normal AMR routines. They should probably be removed in the future.
    
    REAL(KIND=octalKind)               :: subcellSize    ! the size (length of a vertex) of each subcell

    ! This is used only when we construct the tree from SPH data which
    ! contains the position+density+velocitiy of gas particles.
    ! Should be allocated with # of gas particles in this octal
    INTEGER, POINTER                   :: gas_particle_list(:)  ! SPH index of the particles in this octal
    
  END TYPE octal
 
CONTAINS 
 
  TYPE(octalVector) FUNCTION subcellCentre(thisOctal,nChild)
    ! returns the centre of one of the subcells of the current octal 

    IMPLICIT NONE

    TYPE(octal), INTENT(IN) :: thisOctal 
    INTEGER, INTENT(IN)     :: nChild    ! index (1-8) of the subcell
    
    REAL(KIND=octalKind)    :: d 
    
    d = thisOctal%subcellSize / 2.0_oc

    SELECT CASE (nChild)
      CASE (1)    
        subcellCentre = thisOctal%centre - (d * xHatOctal) - (d * yHatOctal) - (d * zHatOctal)
      CASE (2)    
        subcellCentre = thisOctal%centre + (d * xHatOctal) - (d * yHatOctal) - (d * zHatOctal)
      CASE (3)    
        subcellCentre = thisOctal%centre - (d * xHatOctal) + (d * yHatOctal) - (d * zHatOctal)
      CASE (4)    
        subcellCentre = thisOctal%centre + (d * xHatOctal) + (d * yHatOctal) - (d * zHatOctal)
      CASE (5)    
        subcellCentre = thisOctal%centre - (d * xHatOctal) - (d * yHatOctal) + (d * zHatOctal)
      CASE (6)    
        subcellCentre = thisOctal%centre + (d * xHatOctal) - (d * yHatOctal) + (d * zHatOctal)
      CASE (7)    
        subcellCentre = thisOctal%centre - (d * xHatOctal) + (d * yHatOctal) + (d * zHatOctal)
      CASE (8)    
        subcellCentre = thisOctal%centre + (d * xHatOctal) + (d * yHatOctal) + (d * zHatOctal)
      CASE DEFAULT
        PRINT *, "Invalid nChild passed to subcellCentre"
!        DO ; END DO
    END SELECT   
    
  END FUNCTION subcellCentre


  !
  ! For a given octal object and a x,y,z position,  this
  ! function checks if this position is within this octal.
  !
  function within_subcell(this, subcell, x, y, z) RESULT(out)
    implicit none
    logical :: out
    type(octal), intent(in) :: this
    integer, intent(in) :: subcell   
    double precision, intent(in) :: x, y, z
    !
    TYPE(octalVector)     :: cellCenter
    double precision :: x0, y0, z0  ! cell center
    double precision :: d, dp, dm
    double precision :: eps = 0.0d0
    
    d = (this%subcellSize)/2.0d0
    dp = d+eps
    dm = d-eps
    
    cellCenter = subcellCentre(this,subcell)
    x0=dble(cellCenter%x); y0=dble(cellCenter%y); z0=dble(cellCenter%z)

    ! Fortran check the condidtion from
    ! the top, so this should work, and it is faster...
    if ( x > (x0+dp) ) then
       out = .false.
    else if ( x < (x0-dm)) then
       out = .false.      
    elseif ( y > (y0+dp) ) then
       out = .false.
    elseif ( y < (y0-dm)) then
       out = .false.
    elseif ( z > (z0+dp) ) then
       out = .false.
    elseif ( z < (z0-dm)) then
       out = .false.
    else
       out = .true.
    end if

    !    if ( x <= (x0+dp) .and.  x >= (x0-dm) ) then
    !       if ( y <= (y0+dp) .and.  y >= (y0-dm) )  then
    !	  if ( z <= (z0+dp) .and.  z >= (z0-dm) )  then
    !	     out = .true.
    !	  else
    !	     out = .false.
    !	  end if
    !       else
    !	  out = .false.
    !       end if
    !    else
    !       out = .false.
    !    end if

  end function within_subcell

  
  
END MODULE octal_mod
