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

! OK - now for the nightmare of the 2D case, as implemented by TJH
! started on 25/08/04


!       z
!       |
!       |
!       |
!       |
!       |--------|--------+
!       |        |        |
!       |        |        |
!       |    3   |   4    |
!       |        |        |   Diagram showing the convention used here for
!       |________|________|   numbering the subcells of each octal.
!       |        |        |
!       |        |        |
!       |    1   |   2    |
!       |        |        |
!       |________|________|________\ x
!                                  /



  TYPE octalWrapper
    TYPE(octal), POINTER  :: content => NULL()
    LOGICAL(KIND=logic), DIMENSION(8) :: inUse
  END TYPE octalWrapper
 
  TYPE wrapperArray
    TYPE(octalWrapper), DIMENSION(:), POINTER :: wrappers => NULL()  ! a number of octal wrappers  
  END TYPE wrapperArray

  TYPE octalListElement
    TYPE(octal), POINTER  :: content => NULL()
    LOGICAL(KIND=logic), DIMENSION(8) :: inUse
    TYPE(octalListElement), POINTER :: next => NULL()
  END TYPE octalListElement

  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE copyOctalComponents
  END INTERFACE

  
  TYPE octal

    INTEGER                            :: nDepth       ! depth of octal. root is 1, it's childen are 2...
    INTEGER                            :: nChildren    ! how many pointers to children there are (max 8)
    INTEGER                            :: indexChild(8)! index of child array containing
                                                        ! pointer to each subcell's child (if it exists) 
    LOGICAL                            :: threeD        ! this is a three-dimensional octal
    LOGICAL                            :: twoD          ! this is a two-dimensioanl octal (quartal?!)
    INTEGER                            :: maxChildren   ! this is 8 for three-d and 4 for two-d
    TYPE(octal), DIMENSION(:), POINTER :: child => null()
    LOGICAL, DIMENSION(8)              :: hasChild
    TYPE(octal), POINTER               :: parent => null()         
    TYPE(octalVector)                  :: centre

    REAL(double), DIMENSION(8)         :: rho            ! density
    TYPE(vector), DIMENSION(8)         :: velocity       ! velocity
    TYPE(vector), DIMENSION(27)        :: cornerVelocity ! velocity at corners of subcells
    REAL, DIMENSION(8)                 :: temperature    ! grid subcell temperatures
    REAL, DIMENSION(8)                 :: oldTemperature    ! grid subcell temperatures
    REAL(double), DIMENSION(8)         :: distanceGrid   ! distance crossing used by lucy R Eq
    INTEGER, DIMENSION(8)              :: nCrossings     ! no of photon crossings used by lucy R Eq
    REAL(double), DIMENSION(:,:), POINTER      :: kappaAbs => null() ! cont absorption opacities
    REAL(double), DIMENSION(:,:), POINTER      :: kappaSca => null() ! scattering opacities
    REAL(double), DIMENSION(8)                 :: chiLine        ! line opacity
    REAL(double), DIMENSION(8)                 :: etaLine        ! line emissivity
    REAL(double), DIMENSION(8)                 :: etaCont        ! line emissivity
    REAL(double), DIMENSION(8)                 :: biasLine3D     ! grid bias distrubtion
    REAL(double), DIMENSION(8)                 :: biasCont3D     ! grid bias distrubtion
    real(double), DIMENSION(8) :: probDistLine  ! emissivity probabilty distribution
    real(double), DIMENSION(8) :: probDistCont  ! emissivity probabilty distribution
    real(double), DIMENSION(:,:), POINTER ::  N => null()! stateq level pops
    real(double), DIMENSION(8) :: Ne            ! electron density

    real(double), DIMENSION(8) :: NH            ! total H no density
    real(double), DIMENSION(8) :: NHI            ! neutral H
    real(double), DIMENSION(8) :: NHII            ! HII
    real(double), DIMENSION(8) :: NHeI            ! HeI
    real(double), DIMENSION(8) :: NHeII            ! HeII
    real(double), dimension(8) :: HIheating

    real(double), DIMENSION(8) :: nTot          ! total density
    real, dimension(8) :: oldFrac ! Previous value of dust sublimation fraction
    REAL, DIMENSION(:,:), POINTER      :: departCoeff =>null()! temporary storage for departure coefficients
    LOGICAL(KIND=logic), DIMENSION(8) :: inStar      ! point lies within star
    LOGICAL(KIND=logic), DIMENSION(8) :: inFlow      ! inside accretion flow region GET RID OF THIS?    
    INTEGER, DIMENSION(8) :: label                       ! numeric label for each subcell. 
      ! the subcell labels may be useful for debugging the code, but are not needed for
      !   any of the normal AMR routines. They should probably be removed in the future.
    
    real(oct)               :: subcellSize    ! the size (length of a vertex) of each subcell

    ! This is used only when we construct the tree from SPH data which
    ! contains the position+density+velocitiy of gas particles.
    ! Should be allocated with # of gas particles in this octal
    INTEGER, POINTER                   :: gas_particle_list(:) => null() ! SPH index of the particles in this octal
    LOGICAL(KIND=logic), DIMENSION(8) :: changed     ! octal has changed in some way since previous calculation    
    
    INTEGER, DIMENSION(8)                :: dusttype
    real(double), dimension(:,:), pointer        :: dustTypeFraction => null() ! dust type fraction (sum=1)
    INTEGER :: parentSubcell
    logical :: gasOpacity                            ! use gas rather than dust opacity for this cell
    logical, dimension(8)                 :: diffusionApprox
    logical, dimension(8)                 :: leftHandDiffusionBoundary
    real(double), dimension(8) :: diffusionProb
    logical, dimension(8)                 :: rightHandDiffusionBoundary
    logical, dimension(8) :: undersampled
    real, dimension(8) :: nDiffusion
    real, dimension(8)    :: incidentFlux

  END TYPE octal
 
CONTAINS 
 
  TYPE(octalVector) FUNCTION subcellCentre(thisOctal,nChild)
    ! returns the centre of one of the subcells of the current octal 

    IMPLICIT NONE

    TYPE(octal), INTENT(IN) :: thisOctal 
    INTEGER, INTENT(IN)     :: nChild    ! index (1-8) of the subcell
    
    real(oct)    :: d 
    
    d = thisOctal%subcellSize * 0.5_oc

    if (thisOctal%threeD) then  !do the three-d case as per diagram
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
          PRINT *, "Error:: Invalid nChild passed to subcellCentre threed case"
          PRINT *, "        nChild = ", nChild 
          STOP
       END SELECT
    else
       SELECT CASE (nChild)
       CASE (1)    
          subcellCentre = thisOctal%centre - (d * xHatOctal) - (d * zHatOctal)
       CASE (2)    
          subcellCentre = thisOctal%centre + (d * xHatOctal) - (d * zHatOctal)
       CASE (3)    
          subcellCentre = thisOctal%centre - (d * xHatOctal) + (d * zHatOctal)
       CASE (4)    
          subcellCentre = thisOctal%centre + (d * xHatOctal) + (d * zHatOctal)
       CASE DEFAULT
          PRINT *, "Error:: Invalid nChild passed to subcellCentre twoD case"
          PRINT *, "        nChild = ", nChild 
          do;enddo
       END SELECT
    endif
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
    real(double), intent(in) :: x, y, z
    !
    TYPE(octalVector)     :: cellCenter
    real(double) :: x0, y0, z0  ! cell center
    real(double) :: d, dp, dm
    real(double), parameter :: eps = 0.0d0
    
    d = (this%subcellSize)*0.5d0
    dp = d+eps
    dm = d-eps
    
    cellCenter = subcellCentre(this,subcell)
    x0=dble(cellCenter%x); y0=dble(cellCenter%y); z0=dble(cellCenter%z)

    ! Fortran check the condidtion from
    ! the top, so this should work, and it is faster...

    if (this%threeD) then
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
    else ! two-D case
       if ( x > (x0+dp) ) then
          out = .false.
       else if ( x < (x0-dm)) then
          out = .false.      
       elseif ( z > (z0+dp) ) then
          out = .false.
       elseif ( z < (z0-dm)) then
          out = .false.
       else
          out = .true.
       end if
    endif
  end function within_subcell

  SUBROUTINE copyOctalComponents(dest,source)
    ! copy the components within an octal variable to a new octal variable
    !
    ! WARNING: this does not change the parent and child variables - you
    !   must update those yourself elsewhere.
 
    TYPE(octal), INTENT(INOUT) :: dest
    TYPE(octal), INTENT(IN) :: source

    dest%nDepth           = source%nDepth
    dest%nChildren        = source%nChildren
    dest%indexChild       = source%indexChild
    dest%hasChild         = source%hasChild
    dest%centre           = source%centre
    dest%rho              = source%rho
    dest%velocity         = source%velocity
    dest%cornerVelocity   = source%cornerVelocity
    dest%temperature      = source%temperature
    dest%distanceGrid     = source%distanceGrid
    dest%nCrossings       = source%nCrossings
    dest%chiLine          = source%chiLine
    dest%etaLine          = source%etaLine
    dest%etaCont          = source%etaCont
    dest%biasLine3D       = source%biasLine3D
    dest%biasCont3D       = source%biasCont3D
    dest%probDistLine     = source%probDistLine
    dest%probDistCont     = source%probDistCont
    dest%Ne               = source%Ne
    dest%nTot             = source%nTot
    dest%inStar           = source%inStar
    dest%inFlow           = source%inFlow
    dest%label            = source%label
    dest%subcellSize      = source%subcellSize
    dest%changed          = source%changed
    dest%threeD           = source%threeD 
    dest%twoD             = source%twoD   
    dest%maxChildren      = source%maxChildren
    dest%NH               = source%NH     
    dest%NHI              = source%NHI    
    dest%NHII             = source%NHII             
    dest%NHeI             = source%NHeI             
    dest%NHeII            = source%NHeII            
    dest%HIheating        = source%HIheating        
    dest%oldFrac          = source%oldFrac          
    dest%dustType         = source%dustType         
    dest%parentSubcell    = source%parentSubcell    
    dest%gasOpacity       = source%gasOpacity       
    dest%diffusionApprox  = source%diffusionApprox  
    dest%leftHandDiffusionBoundary = source%leftHandDiffusionBoundary        
    dest%rightHandDiffusionBoundary = source%rightHandDiffusionBoundary       
    dest%diffusionProb    = source%diffusionProb
    dest%underSampled     = source%underSampled
    dest%nDiffusion       = source%nDiffusion
    dest%incidentFlux     = source%incidentFlux

    IF (ASSOCIATED(dest%kappaAbs)) DEALLOCATE (dest%kappaAbs)
    IF (ASSOCIATED(dest%kappaSca)) DEALLOCATE (dest%kappaSca)
    IF (ASSOCIATED(dest%N))        DEALLOCATE (dest%N)
    IF (ASSOCIATED(dest%departCoeff)) DEALLOCATE (dest%departCoeff)
    IF (ASSOCIATED(dest%gas_particle_list)) DEALLOCATE (dest%gas_particle_list)
    
    IF (ASSOCIATED(source%kappaAbs)) THEN                   
      ALLOCATE(dest%kappaAbs( SIZE(source%kappaAbs,1),       &
                              SIZE(source%kappaAbs,2)))
      dest%kappaAbs = source%kappaAbs
    END IF  

    IF (ASSOCIATED(source%kappaSca)) THEN
      ALLOCATE(dest%kappaSca( SIZE(source%kappaSca,1),       &
                              SIZE(source%kappaSca,2)))
      dest%kappaSca = source%kappaSca
    END IF  

    IF (ASSOCIATED(source%N)) THEN 
      ALLOCATE(dest%N( SIZE(source%N,1),              &
                       SIZE(source%N,2)))
      dest%N = source%N
    END IF  

    IF (ASSOCIATED(source%departCoeff)) THEN
      ALLOCATE(dest%departCoeff( SIZE(source%departCoeff,1),    &
                                 SIZE(source%departCoeff,2)))
      dest%departCoeff = source%departCoeff
    END IF  

    IF (ASSOCIATED(source%gas_particle_list)) THEN
      ALLOCATE(dest%gas_particle_list(SIZE(source%gas_particle_list)))
      dest%gas_particle_list = source%gas_particle_list
    END IF  

    IF (ASSOCIATED(source%dustTypeFraction)) THEN                   
      ALLOCATE(dest%kappaAbs( SIZE(source%dustTypeFraction,1),       &
                              SIZE(source%dustTypeFraction,2)))
      dest%kappaAbs = source%kappaAbs
    END IF  

  END SUBROUTINE copyOctalComponents
  
  
END MODULE octal_mod
