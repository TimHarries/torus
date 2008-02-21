MODULE octal_mod
  ! data type and some routines needed by for adaptive mesh refinement. nhs.
  
  ! these are in a separate module from amr_mod because they need to be
  !   accessed from module(s) that are themselves USEd by amr_mod.

  USE kind_mod
  USE vector_mod
  USE messages_mod

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

  TYPE octal

    INTEGER                            :: nDepth       ! depth of octal. root is 1, it's childen are 2...
    integer                            :: mpiThread(8)  ! the thread number associated with this subcell (from 0 to n)
    INTEGER                            :: nChildren    ! how many pointers to children there are (max 8)
    INTEGER                            :: indexChild(8)! index of child array containing
                                                        ! pointer to each subcell's child (if it exists) 
    LOGICAL                            :: threeD        ! this is a three-dimensional octal
    LOGICAL                            :: twoD          ! this is a two-dimensional octal (quartal?!)
    LOGICAL                            :: oneD          ! this is a one-dimensional octal (bital?!)
    LOGICAL                            :: cylindrical   ! A three-d amr grid of x,z,phi
    LOGICAL                            :: splitAzimuthally   ! A three-d amr grid of x,z,phi
    INTEGER                            :: maxChildren   ! this is 8 for three-d and 4 for two-d
    TYPE(octal), DIMENSION(:), POINTER :: child => NULL()
    LOGICAL, DIMENSION(8)              :: hasChild
    TYPE(octal), POINTER               :: parent => null()         
    TYPE(octalVector)                  :: centre
    real(double)                       :: r
    REAL(double), DIMENSION(8)         :: rho            ! density
    TYPE(vector), DIMENSION(8)         :: velocity       ! velocity
    real(double), dimension(8)         :: microturb
    TYPE(vector), DIMENSION(27)        :: cornerVelocity ! velocity at corners of subcells
    REAL, DIMENSION(8)                 :: temperature    ! grid subcell temperatures (gas or dust)
    REAL, DIMENSION(8)                 :: temperaturedust! grid subcell dust temperatures
    REAL, DIMENSION(8)                 :: temperaturegas ! grid subcell gas temperatures
    
    real(double), dimension(8) :: eDens
    REAL, DIMENSION(8)                 :: oldTemperature    ! grid subcell temperatures
    REAL(double), DIMENSION(8)                 :: oldeDens
    REAL(double), DIMENSION(8)                 :: kappaRoss

    REAL(double), DIMENSION(8)         :: distanceGrid   ! distance crossing used by lucy R Eq



    INTEGER, DIMENSION(8)              :: nCrossings     ! no of photon crossings used by lucy R Eq
    REAL(double), DIMENSION(:,:), POINTER      :: kappaAbs => null() ! cont absorption opacities
    REAL(double), DIMENSION(:,:), POINTER      :: kappaSca => null() ! scattering opacities
    REAL(double), DIMENSION(8)                 :: chiLine        ! line opacity
    REAL(double), DIMENSION(8)                 :: etaLine        ! line emissivity
    REAL(double), DIMENSION(8)                 :: etaCont        ! line emissivity
    REAL(double), DIMENSION(8)                 :: biasLine3D     ! grid bias distribution
    REAL(double), DIMENSION(8)                 :: biasCont3D     ! grid bias distribution
    real(double), DIMENSION(8) :: probDistLine  ! emissivity probability distribution
    real(double), DIMENSION(8) :: probDistCont  ! emissivity probability distribution
    real(double), DIMENSION(:,:), POINTER ::  N => null()! stateq level pops
    real(double), DIMENSION(8) :: Ne            ! electron density
    real(double)               :: phi, dphi
    real(double), DIMENSION(8) :: NH            ! total H no density
    real(double), DIMENSION(8) :: NH2            ! total H2 no density
    real(double), pointer :: molecularLevel(:,:) => null() ! molecular level populations
    real(double), pointer :: newmolecularLevel(:,:) => null() ! molecular level populations
    real(double), pointer :: oldmolecularLevel(:,:) => null() ! molecular level populations
    real(double), pointer :: atomLevel(:,:,:) => null() ! atom level populations
    real(double), pointer :: atomAbundance(:,:) ! abundances
    real(double), pointer :: newatomLevel(:,:,:) => null() ! atom level populations
    real(double), pointer :: jnu(:,:) => null() ! Molecular radiation field
    real(double), pointer :: jnuCont(:,:) => null()
    real(double), pointer :: jnuLine(:,:) => null()
    real(double), pointer :: jnugrid(:,:) => null() ! molecular level populations
    real(double), pointer :: bnu(:,:) => null() ! 
    real, pointer :: molAbundance(:) => null() ! molecular abundances ! only 1D because only deal with one molecule at a time
    real(double), DIMENSION(8) :: NHI            ! neutral H
    real(double), DIMENSION(8) :: NHII            ! HII
    real(double), DIMENSION(8) :: NHeI            ! HeI
    real(double), DIMENSION(8) :: NHeII            ! HeII
    real(double), dimension(8) :: Hheating
    real(double), dimension(8) :: Heheating
   
    real(double), dimension(:,:), pointer  :: ionFrac => null()
    real(double), dimension(:,:), pointer  :: photoIonCoeff  => null()
    real(double) :: diffusionCoeff(8)

    real(double), DIMENSION(8) :: nTot          ! total density
    real, dimension(8) :: oldFrac ! Previous value of dust sublimation fraction
    REAL, DIMENSION(:,:), POINTER      :: departCoeff =>null()! temporary storage for departure coefficients
    LOGICAL(KIND=logic), DIMENSION(8) :: inStar      ! point lies within star
    INTEGER, DIMENSION(8) :: nDirectPhotons
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
    real(double), pointer, dimension(:,:,:) :: mpiBoundaryStorage => null()
    INTEGER, DIMENSION(8)                :: dusttype
    real(double), dimension(:,:), pointer        :: dustTypeFraction => null() ! dust type fraction (sum=1)
    INTEGER :: parentSubcell
    logical :: gasOpacity                            ! use gas rather than dust opacity for this cell
    logical, dimension(8)                 :: diffusionApprox
    logical, dimension(8) :: undersampled
    real, dimension(8) :: nDiffusion
    integer :: boundaryCondition(8)

    ! hydrodynamics
    real(double) :: q_i(8), q_i_plus_1(8), q_i_minus_1(8), q_i_minus_2(8)
    real(double) :: x_i(8), x_i_plus_1(8), x_i_minus_1(8)
    real(double) :: u_interface(8), u_i_plus_1(8), u_i_minus_1(8)
    real(double) :: flux_i(8), flux_i_plus_1(8), flux_i_minus_1(8)
    real(double) :: phiLimit(8), rLimit(8)
    logical :: ghostCell(8), feederCell(8)
    logical :: edgeCell(8), refinedLastTime(8)
    real(double) :: rhou(8),  rhov(8), rhow(8), rhoE(8), energy(8)
    real(double) :: pressure_i(8), pressure_i_plus_1(8), pressure_i_minus_1(8)
    real(double), pointer :: tempStorage(:,:) => null()
    type(OCTALVECTOR) :: boundaryPartner(8)

    real(double) :: w(8,4), fluxc(8,3),  a(8,3)
    real(double) :: ac2(8,3), flux(8,3), ac1(8,3)
    real(double) :: qState(8,5), fluxvector(8,5), newFluxVector(8,5), qstate_i_minus_1(8,5)
    real(double) :: rho_i_minus_1(8), rho_i_plus_1(8)
    real(double) :: phi_i(8), phi_i_plus_1(8), phi_i_minus_1(8)
  END TYPE octal
 
CONTAINS 
 
  TYPE(octalVector) FUNCTION subcellCentre(thisOctal,nChild)
    ! returns the centre of one of the subcells of the current octal 

    IMPLICIT NONE

    TYPE(octal), INTENT(IN) :: thisOctal 
    INTEGER, INTENT(IN)     :: nChild    ! index (1-8) of the subcell
    type(OCTALVECTOR) :: rVec
    real(oct)    :: d

    d = thisOctal%subcellSize * 0.5_oc

    if (thisOctal%oneD) then
       select case (nchild)
          case(1)
             subcellCentre = thisOctal%centre - d * xHatOctal
          case(2)
             subcellCentre = thisOctal%centre + d * xHatOctal
          case DEFAULT
             write(*,*) "bug - one-d cell has more than 3 children"
       end select
       goto 666
    endif

    if (thisOctal%threeD) then  !do the three-d case as per diagram
       if (.not.thisOctal%cylindrical) then
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
          rVec = OCTALVECTOR(thisOctal%r,0.d0,thisOctal%centre%z)
          if (thisOctal%splitAzimuthally) then
             SELECT CASE (nChild)
             CASE (1)    
                subcellCentre = rVec - (d * xHatOctal) - (d * zHatOctal)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi-0.25d0*thisOctal%dphi))
             CASE (2)    
                subcellCentre = rVec + (d * xHatOctal) - (d * zHatOctal)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi-0.25d0*thisOctal%dphi))
             CASE (3)    
                subcellCentre = rVec - (d * xHatOctal) + (d * zHatOctal)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi-0.25d0*thisOctal%dphi))
             CASE (4)    
                subcellCentre = rVec + (d * xHatOctal) + (d * zHatOctal)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi-0.25d0*thisOctal%dphi))
             CASE (5)    
                subcellCentre = rVec - (d * xHatOctal) - (d * zHatOctal)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi-0.25d0*thisOctal%dphi))
             CASE (6)    
                subcellCentre = rVec + (d * xHatOctal) - (d * zHatOctal)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi-0.25d0*thisOctal%dphi))
             CASE (7)    
                subcellCentre = rVec - (d * xHatOctal) + (d * zHatOctal)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi-0.25d0*thisOctal%dphi))
             CASE (8)    
                subcellCentre = rVec + (d * xHatOctal) + (d * zHatOctal)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi-0.25d0*thisOctal%dphi))
             CASE DEFAULT
                PRINT *, "Error:: Invalid nChild passed to subcellCentre twoD case 1"
                PRINT *, "        nChild = ", nChild 
                do
                enddo
             END SELECT
          else
             SELECT CASE (nChild)
             CASE (1)    
                subcellCentre = rVec - (d * xHatOctal) - (d * zHatOctal)
                subcellCentre = rotateZ(subcellCentre,-thisOctal%phi)
             CASE (2)    
                subcellCentre = rVec + (d * xHatOctal) - (d * zHatOctal)
                subcellCentre = rotateZ(subcellCentre,-thisOctal%phi)
             CASE (3)    
                subcellCentre = rVec - (d * xHatOctal) + (d * zHatOctal)
                subcellCentre = rotateZ(subcellCentre,-thisOctal%phi)
             CASE (4)    
                subcellCentre = rVec + (d * xHatOctal) + (d * zHatOctal)
                subcellCentre = rotateZ(subcellCentre,-thisOctal%phi)
             CASE DEFAULT
                PRINT *, "Error:: Invalid nChild passed to subcellCentre twoD case 2"
                PRINT *, "        nChild = ", nChild 
                do
                enddo
             END SELECT
          endif
       endif
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
       PRINT *, "Error:: Invalid nChild passed to subcellCentre twoD case 3"
       PRINT *, "        nChild = ", nChild 
       do;enddo
       END SELECT
    endif

666 continue
  END FUNCTION subcellCentre

  real(double) FUNCTION subcellRadius(thisOctal,nChild)
    ! returns the radius of one of the subcells of the current octal 

    IMPLICIT NONE

    TYPE(octal), INTENT(IN) :: thisOctal 
    INTEGER, INTENT(IN)     :: nChild    ! index (1-8) of the subcell
    real(double) :: d

    d = thisOctal%subcellSize * 0.5_oc

    if (thisOctal%oneD) then
       select case(nchild)
          case(1)
             subcellRadius  = thisOctal%centre%x + d
          case(2)
             subcellRadius  = thisOctal%centre%x - d
        end select
        goto 666
     endif



    if (thisOctal%cylindrical) then
       if (thisOctal%splitAzimuthally) then
          SELECT CASE (nChild)
             CASE (1)    
                subcellRadius = thisOctal%r - d
             CASE (2)    
                subcellRadius = thisOctal%r + d
             CASE (3)    
                subcellRadius = thisOctal%r - d
             CASE (4)    
                subcellRadius = thisOctal%r + d
             CASE (5)    
                subcellRadius = thisOctal%r - d
             CASE (6)    
                subcellRadius = thisOctal%r + d
             CASE (7)    
                subcellRadius = thisOctal%r - d
             CASE (8)    
                subcellRadius = thisOctal%r + d
             CASE DEFAULT
                PRINT *, "Error:: Invalid nChild passed to subcellCentre twoD case 4"
                PRINT *, "        nChild = ", nChild 
                do
                enddo
             END SELECT
          else
             SELECT CASE (nChild)
             CASE (1)    
                subcellRadius = thisOctal%r - d
             CASE (2)    
                subcellRadius = thisOctal%r + d
             CASE (3)    
                subcellRadius = thisOctal%r - d
             CASE (4)    
                subcellRadius = thisOctal%r + d
             CASE DEFAULT
                PRINT *, "Error:: Invalid nChild passed to subcellRadius non-split case"
                PRINT *, "        nChild = ", nChild 
                do
                enddo
             END SELECT
          endif
       else
          call writeFatal("subcellRadius called for non-cylindrical grid")
       endif
666 continue
     END FUNCTION subcellRadius


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
    real(double) :: d, dp, dm, r, phi, r0, phi0, dphi
    real(double), parameter :: eps = 0.0d0
    
    d = (this%subcellSize)*0.5d0
    dp = d+eps
    dm = d-eps

    if (this%oneD) then
       write(*,*) "one-d case not implemented in within_subcell"
       stop
    endif
    
    cellCenter = subcellCentre(this,subcell)
    x0=dble(cellCenter%x); y0=dble(cellCenter%y); z0=dble(cellCenter%z)

    ! Fortran check the condidtion from
    ! the top, so this should work, and it is faster...

    if (this%threeD) then
       if (.not.this%cylindrical) then
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
       else
          r = sqrt(x**2+y**2)
          phi = atan2(y,x)
          if (phi < 0.d0) phi = phi + twoPi
          r0 = sqrt(x0**2 + y0**2)
          phi0 = atan2(y0,x0)
          
          if (phi0 < 0.d0) phi0 = phi0 + twoPi
          if (this%splitAzimuthally) then
             dphi = this%dphi/4.d0
          else
             dphi = this%dphi/2.d0
          endif

          if ( z > (z0+dp) ) then
             out = .false.
          elseif ( z < (z0-dm)) then
             out = .false.
          elseif (r > (r0+dp)) then
             out = .false.
          elseif (r < (r0-dm)) then
             out = .false.
          elseif (phi > (phi0+dphi)) then
             out = .false.
          elseif (phi < (phi0-dphi)) then
             out = .false.
          else
             out = .true.
          end if
       endif

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


  real(double) function cellVolume(thisOctal, subcell) result(v)
    type(OCTAL) :: thisOctal
    integer :: subcell
    real(double) :: r1, r2, dphi
    type(OCTALVECTOR) :: rVec
  
    if (thisOctal%oneD) then
       rVec = subcellCentre(thisOctal, subcell)
       r1 = rVec%x - thisOctal%subcellSize/2.d0
       r2 = rVec%x + thisOctal%subcellSize/2.d0
       v = (fourPi / 3.d0) * (r2-r1)**3
       goto 666
    endif

    if (thisOctal%threed) then
       if (.not.thisOctal%cylindrical) then
          v = thisOctal%subcellsize**3
       else
          rVec = subcellCentre(thisOctal,subcell)
          r1 = sqrt(rVec%x**2 + rVec%y**2) - thisOctal%subcellSize/2.d0
          r2 = sqrt(rVec%x**2 + rVec%y**2) + thisOctal%subcellSize/2.d0
!          if (thisOctal%splitAzimuthally) then
!             dPhi = thisOctal%dPhi / 2.d0
!          else
!             dPhi = thisOctal%dPhi 
!          endif
          dPhi = returndPhi(thisOctal)*2.d0
          v = (dphi/dble(twoPi)) * dble(pi) * (r2**2 - r1**2) * thisOctal%subcellSize
       endif
    else
       rVec = subcellCentre(thisOctal,subcell)
       r1 = rVec%x-thisOctal%subcellSize/2.d0
       r2 = rVec%x+thisOctal%subcellSize/2.d0
       v = dble(pi) * (r2**2 - r1**2) * thisOctal%subcellSize
    endif
666 continue
  end function cellVolume


  real(double) function returnDphi(thisOctal)
    type(OCTAL)  :: thisOctal

    if (thisOctal%splitAzimuthally) then
       returndPhi = thisOctal%dPhi / 4.d0
    else
       returndPhi = thisOctal%dPhi / 2.d0
    endif
  end function returnDphi


END MODULE octal_mod
