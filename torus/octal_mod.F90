MODULE octal_mod
  ! data type and some routines needed by for adaptive mesh refinement. nhs.
  
  ! these are in a separate module from amr_mod because they need to be
  !   accessed from module(s) that are themselves USEd by amr_mod.

  USE kind_mod
  use constants_mod
  USE vector_mod
  USE messages_mod
  
  IMPLICIT NONE

  public :: subcellCentre, within_subcell

  logical, public :: cart2d ! 2D rhd in cartesian coords


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

  interface allocateAttribute
     module procedure allocateAttributeDouble
     module procedure allocateAttributeDouble2D
     module procedure allocateAttributeDouble3D
     module procedure allocateAttributeReal
     module procedure allocateAttributeInteger
     module procedure allocateAttributeLogical
     module procedure allocateAttributeVector
  end interface

  interface deallocateAttribute
     module procedure deallocateAttributeDouble
     module procedure deallocateAttributeOctal
     module procedure deallocateAttributeDouble2D
     module procedure deallocateAttributeReal2D
     module procedure deallocateAttributeInteger2D
     module procedure deallocateAttributeDouble3D
     module procedure deallocateAttributeReal
     module procedure deallocateAttributeInteger
     module procedure deallocateAttributeInteger3D
     module procedure deallocateAttributeOctalPointer3D
     module procedure deallocateAttributeVector
     module procedure deallocateAttributeLogical
  end interface

  interface copyAttribute
     module procedure copyAttributeDoublePointer
     module procedure copyAttributeDoublePointer2d
     module procedure copyAttributeDoublePointer3d
     module procedure copyAttributeRealPointer
     module procedure copyAttributeIntegerPointer
     module procedure copyAttributeLogicalPointer
     module procedure copyAttributeVectorPointer
  end interface

  interface readAttributePointer
     module procedure readAttributeDoublePointer
     module procedure readAttributeRealPointer
     module procedure readAttributeIntegerPointer
     module procedure readAttributeVectorPointer
  end interface

  interface writeAttributePointer
     module procedure writeAttributeDoublePointer
     module procedure writeAttributeRealPointer
     module procedure writeAttributeIntegerPointer
     module procedure writeAttributeVectorPointer
  end interface

  interface readAttributeStatic
     module procedure readAttributeStaticDouble
     module procedure readAttributeStaticDoubleSingle
     module procedure readAttributeStaticReal
     module procedure readAttributeStaticInteger
     module procedure readAttributeStaticIntegerSingle
     module procedure readAttributeStaticLogical
     module procedure readAttributeStaticLogicalSingle
     module procedure readAttributeStaticVector
     module procedure readAttributeStaticVectorSingle
  end interface

  interface writeAttributeStatic
     module procedure writeAttributeStaticDouble
     module procedure writeAttributeStaticDoubleSingle
     module procedure writeAttributeStaticReal
     module procedure writeAttributeStaticInteger
     module procedure writeAttributeStaticIntegerSingle
     module procedure writeAttributeStaticLogical
     module procedure writeAttributeStaticLogicalSingle
     module procedure writeAttributeStaticVector
     module procedure writeAttributeStaticVectorSingle
  end interface


  TYPE octalWrapper
    TYPE(octal), POINTER  :: content => NULL()
    LOGICAL, DIMENSION(8) :: inUse
  END TYPE octalWrapper
 
  TYPE wrapperArray
    TYPE(octalWrapper), DIMENSION(:), POINTER :: wrappers => NULL()  ! a number of octal wrappers  
  END TYPE wrapperArray

  TYPE octalListElement
    TYPE(octal), POINTER  :: content => NULL()
    LOGICAL, DIMENSION(8) :: inUse
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
    TYPE(vector)                  :: centre
    real(double)                       :: xMax, xMin, yMax, yMin, zMax, zMin
    real(double)                       :: r
    REAL(double), DIMENSION(8)         :: rho            ! density
    INTEGER, DIMENSION(8) :: label                       ! numeric label for each subcell. 
    integer, pointer :: iEquationOfState(:)  => null()
    real(double), pointer :: gamma(:) => null()
    real(double), pointer :: divV(:) => null()
    REAL, DIMENSION(8)                 :: temperature    ! grid subcell temperatures (gas or dust)
#ifdef PHOTOION
    real, DIMENSION(8)        :: TLastIter     !Temperature at last iteration thaw
    real, dimension(8)        :: TLastLastIter !Temperature at -2 iterations
#endif
    real(oct)               :: subcellSize    ! the size (length of a vertex) of each subcell
#ifdef SPH
    real(double) :: h(8) ! used for cluster geometries where the smoothing length of a cell can be stored
#endif

    LOGICAL, DIMENSION(8)              :: inFlow

    TYPE(OCTALPOINTER), pointer :: neighbourOctal(:,:,:) => null() ! pointer to neighbour
    integer, pointer :: neighbourSubcell(:,:,:)  => null()

    TYPE(vector), DIMENSION(8)         :: velocity       ! velocity
!    TYPE(vector), DIMENSION(8)         :: quadvelocity       ! DAR very temporary velocity expires 28 Feb 2009
!    TYPE(vector), DIMENSION(8)         :: linearvelocity       ! DAR very temporary velocity expires 28 Feb 2009
    TYPE(vector), pointer, DIMENSION(:)    :: cornerVelocity => null() ! velocity at corners of subcells
    real(double), DIMENSION(:), pointer    :: cornerrho => null() ! velocity at corners of subcells
    real(double)               :: phi, dphi, phimin, phimax

    real(double), dimension(:,:,:), pointer :: qViscosity => null()
    
    logical, dimension(:), pointer                 :: diffusionApprox => null()
    logical, dimension(:), pointer                 :: fixedTemperature => null()
    real, dimension(:), pointer :: nDiffusion => null()
    real(double), dimension(:), pointer :: eDens => null()
    real(double), pointer :: diffusionCoeff(:) => null()
    REAL(double), DIMENSION(:), pointer                 :: oldeDens => null()
    INTEGER, DIMENSION(:), pointer :: nDirectPhotons => null()

    logical, dimension(:), pointer :: undersampled => null()
    REAL, DIMENSION(:), pointer                 :: oldTemperature  => null()   ! grid subcell temperatures
    REAL(double), DIMENSION(:), pointer                 :: kappaRoss => null()
    REAL(double), DIMENSION(:), pointer         :: distanceGrid  => null()  ! distance crossing used by lucy R Eq
    real(double), dimension(:,:,:), pointer       :: scatteredIntensity => null()
    real(double), dimension(:), pointer       :: meanIntensity => null()
    INTEGER, DIMENSION(:), pointer              :: nCrossings  => null()    ! no of photon crossings used by lucy R Eq
    real(double), DIMENSION(:), pointer :: nTot => null()          ! total density
    real, dimension(:), pointer :: oldFrac  => null() ! Previous value of dust sublimation fraction

    INTEGER, DIMENSION(:), pointer                :: dusttype => null()
    real(double), dimension(:,:), pointer        :: dustTypeFraction => null() ! dust type fraction (sum=1)


    REAL, DIMENSION(:,:), POINTER      :: departCoeff =>null()! temporary storage for departure coefficients
    LOGICAL, DIMENSION(8) :: inStar      ! point lies within star
    REAL(double), DIMENSION(:,:), POINTER      :: kappaAbs => null() ! cont absorption opacities
    REAL(double), DIMENSION(:,:), POINTER      :: kappaSca => null() ! scattering opacities

    REAL(double), DIMENSION(:), pointer                 :: chiLine  => null()       ! line opacity
    REAL(double), DIMENSION(:), pointer                 :: etaLine  => null()       ! line emissivity
    REAL(double), DIMENSION(:), pointer                 :: etaCont   => null()      ! line emissivity
    REAL(double), DIMENSION(:), pointer                 :: biasLine3D  => null()    ! grid bias distribution
    REAL(double), DIMENSION(:), pointer                 :: biasCont3D  => null()   ! grid bias distribution
    REAL(double), DIMENSION(:), pointer                 :: source  => null()   ! source for poissons equation
    real(double), DIMENSION(:), pointer :: probDistLine  => null() ! emissivity probability distribution
    real(double), DIMENSION(:), pointer :: probDistCont  => null()  ! emissivity probability distribution
    real(double), DIMENSION(:,:), POINTER :: N => null()! stateq level pops
    real(double), DIMENSION(:), pointer :: Ne  => null()           ! electron density
    real(double), DIMENSION(:), pointer :: NH  => null()           ! total H no density
    real(double), pointer :: molecularLevel(:,:) => null() ! molecular level populations
    real(double), pointer :: molcellparam(:,:) => null()
    real(double), pointer :: newmolecularLevel(:,:) => null() ! molecular level populations
    real(double), pointer :: oldmolecularLevel(:,:) => null() ! molecular level populations
    real(double), pointer :: oldestmolecularLevel(:,:) => null() ! molecular level populations

! time dependent RT stuff

    real(double), pointer :: Adot(:) => null()
    real(double), pointer :: uDens(:) => null()
    real(double), pointer :: distanceGridAdot(:) => null()
    real(double), pointer :: distanceGridPhotonFromSource(:) => null()
    real(double), pointer :: distanceGridPhotonFromGas(:) => null()
    real(double), pointer :: photonEnergyDensityFromSource(:) => null()
    real(double), pointer :: photonEnergyDensityFromGas(:) => null()
    real(double), pointer :: photonEnergyDensity(:) => null()
    real(double), pointer :: oldphotonEnergyDensity(:) => null()


    REAL, DIMENSION(:), pointer                 :: temperaturedust=> null() ! grid subcell dust temperatures
    REAL, DIMENSION(:), pointer                 :: temperaturegas => null() ! grid subcell gas temperatures
    real(double), DIMENSION(:), pointer :: NH2 => null()           ! total H2 no density
    real(double), pointer, dimension(:)         :: microturb => null()
    real(double), dimension(:), pointer         :: molmicroturb => null()



    real(double), pointer :: atomLevel(:,:,:) => null() ! atom level populations
    real(double), pointer :: atomAbundance(:,:) => null() ! abundances
    real(double), pointer :: newatomLevel(:,:,:) => null() ! atom level populations
    real(double), pointer :: jnu(:,:) => null() ! Molecular radiation field
    real(double), pointer :: jnuCont(:,:) => null()
    real(double), pointer :: jnuLine(:,:) => null()
    real(double), pointer :: tau(:,:) => null() ! molecular level populations
    real(double), pointer :: bnu(:,:) => null() ! 
    real, pointer :: molAbundance(:) => null() ! molecular abundances ! only 1D because only deal with one molecule at a time
    real, pointer :: convergence(:) => null() ! convergence
    integer, pointer :: levelconvergence(:,:) => null() ! convergence
    integer, pointer :: nsplit(:) => null() ! nsplit
    real(double), pointer, DIMENSION(:) :: NHI => null()           ! neutral H
    real(double),  pointer, DIMENSION(:) :: NHII  => null()           ! HII
    real(double),  pointer, DIMENSION(:) :: NHeI => null()            ! HeI
    real(double),  pointer, DIMENSION(:) :: NHeII => null()            ! HeII
    real(double),  pointer, dimension(:) :: Hheating => null()
    real(double),  pointer, dimension(:) :: Heheating => null()
   
    real(double), dimension(:,:), pointer  :: ionFrac => null()
    real(double), dimension(:,:), pointer  :: photoIonCoeff  => null()
    real(double), dimension(:,:), pointer :: sourceContribution => null()
    real(double), dimension(:,:), pointer :: diffuseContribution => null()
    real(double), dimension(:,:), pointer :: normSourceContribution => null()


      ! the subcell labels may be useful for debugging the code, but are not needed for
      !   any of the normal AMR routines. They should probably be removed in the future.
    

    ! This is used only when we construct the tree from SPH data which
    ! contains the position+density+velocitiy of gas particles.
    ! Should be allocated with # of gas particles in this octal
    INTEGER, POINTER                   :: gas_particle_list(:) => null() ! SPH index of the particles in this octal
    LOGICAL, DIMENSION(:), pointer :: changed => null()    ! octal has changed in some way since previous calculation
    real(double), pointer, dimension(:,:,:) :: mpiBoundaryStorage => null()
    real(double), pointer, dimension(:,:,:) :: mpiCornerStorage => null()
    INTEGER :: parentSubcell
    logical :: gasOpacity                            ! use gas rather than dust opacity for this cell




    ! hydrodynamics
    real(double), pointer :: q_i(:) => null(), q_i_plus_1(:) => null(), q_i_minus_1(:) => null(), q_i_minus_2(:) => null()
    real(double), pointer :: x_i(:) => null(), x_i_plus_1(:) => null(), x_i_minus_1(:) => null(), x_i_minus_2(:) => null()
    real(double), pointer :: u_interface(:) => null(), u_i_plus_1(:) => null(), u_i_minus_1(:) => null(), u_i(:) => null()
    real(double), pointer :: flux_i(:) => null(), flux_i_plus_1(:) => null(), flux_i_minus_1(:) => null()
    real(double), pointer :: phiLimit(:) => null(), rLimit(:) => null()
    logical, pointer :: ghostCell(:) => null(), feederCell(:) => null(), corner(:) => null()
    logical, pointer :: edgeCell(:) => null(), refinedLastTime(:) => null()
    real(double), pointer :: rhou(:) => null(),  rhov(:) => null(), rhow(:) => null(), rhoE(:) => null(), energy(:) => null()
    real(double), pointer :: pressure_i(:) => null(), pressure_i_plus_1(:) => null(), pressure_i_minus_1(:) => null()
    !THAW - temporary
    real(double), pointer :: rhoeLastTime(:) => null()
    real(double), pointer :: tempStorage(:,:) => null()
    type(vECTOR), pointer :: boundaryPartner(:) => null()
    type(vECTOR), pointer :: gravboundaryPartner(:) => null()
    type(VECTOR), pointer :: radiationMomentum(:) => null()
    type(VECTOR), pointer :: kappaTimesFlux(:) => null()
    real(double), pointer :: phi_i(:) => null(), phi_i_plus_1(:) => null(), phi_i_minus_1(:) => null()
    real(double), pointer :: phi_stars(:) => null(), phi_gas(:) => null(), phi_gas_corr(:) => null()
    real(double),pointer :: rho_i_minus_1(:) => null(), rho_i_plus_1(:) => null()
    real(double),pointer :: rhorv_i_minus_1(:) => null(), rhorv_i_plus_1(:) => null()
    type(VECTOR), pointer :: fViscosity(:) => null()
    integer, pointer :: boundaryCondition(:) => null()
    logical, pointer :: boundaryCell(:) => null()

  END TYPE octal
 
  TYPE octalPointer
     TYPE(OCTAL), pointer :: pointer => null()
  end TYPE octalPointer
     
CONTAINS 
 
  TYPE(Vector) FUNCTION subcellCentre(thisOctal,nChild)
    ! returns the centre of one of the subcells of the current octal 

    IMPLICIT NONE

    TYPE(octal), INTENT(IN) :: thisOctal 
    INTEGER, INTENT(IN)     :: nChild    ! index (1-8) of the subcell
    type(VECTOR) :: rVec
    real(oct)    :: d

    d = thisOctal%subcellSize * 0.5_oc

    if (thisOctal%oneD) then
       select case (nchild)
          case(1)
             subcellCentre = thisOctal%centre - d * xHat
          case(2)
             subcellCentre = thisOctal%centre + d * xHat
          case DEFAULT
             write(*,*) "bug - one-d cell has more than 3 children"
       end select
       goto 666
    endif

    if (thisOctal%threeD) then  !do the three-d case as per diagram
       if (.not.thisOctal%cylindrical) then

          SELECT CASE (nChild)
          CASE (1)    
             subcellCentre = thisOctal%centre + d * VECTOR(-1.d0,-1.d0,-1.d0)
          CASE (2)    
             subcellCentre = thisOctal%centre + d * VECTOR(1.d0,-1.d0,-1.d0)
          CASE (3)    
             subcellCentre = thisOctal%centre + d * VECTOR(-1.d0,1.d0,-1.d0)
          CASE (4)    
             subcellCentre = thisOctal%centre + d * VECTOR(1.d0,1.d0,-1.d0)
          CASE (5)    
             subcellCentre = thisOctal%centre + d * VECTOR(-1.d0,-1.d0,1.d0)
          CASE (6)    
             subcellCentre = thisOctal%centre + d * VECTOR(1.d0,-1.d0,1.d0)
          CASE (7)    
             subcellCentre = thisOctal%centre + d * VECTOR(-1.d0,1.d0,1.d0)
          CASE (8)    
             subcellCentre = thisOctal%centre + d * VECTOR(1.d0,1.d0,1.d0)
          CASE DEFAULT
             PRINT *, "Error:: Invalid nChild passed to subcellCentre threed case"
             PRINT *, "        nChild = ", nChild 
             STOP
          END SELECT
       else
          rVec = VECTOR(thisOctal%r,0.d0,thisOctal%centre%z)
          if (thisOctal%splitAzimuthally) then
             SELECT CASE (nChild)
             CASE (1)    
                subcellCentre = rVec - (d * xHat) - (d * zHat)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi-0.25d0*thisOctal%dphi))
             CASE (2)    
                subcellCentre = rVec + (d * xHat) - (d * zHat)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi-0.25d0*thisOctal%dphi))
             CASE (3)    
                subcellCentre = rVec - (d * xHat) - (d * zHat)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi+0.25d0*thisOctal%dphi))
             CASE (4)    
                subcellCentre = rVec + (d * xHat) - (d * zHat)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi+0.25d0*thisOctal%dphi))
             CASE (5)    
                subcellCentre = rVec - (d * xHat) + (d * zHat)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi-0.25d0*thisOctal%dphi))
             CASE (6)    
                subcellCentre = rVec + (d * xHat) + (d * zHat)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi-0.25d0*thisOctal%dphi))
             CASE (7)    
                subcellCentre = rVec - (d * xHat) + (d * zHat)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi+0.25d0*thisOctal%dphi))
             CASE (8)    
                subcellCentre = rVec + (d * xHat) + (d * zHat)
                subcellCentre = rotateZ(subcellCentre,-(thisOctal%phi+0.25d0*thisOctal%dphi))
             CASE DEFAULT
                PRINT *, "Error:: Invalid nChild passed to subcellCentre cylindical 3D case 1"
                PRINT *, "        nChild = ", nChild 
                do
                enddo
             END SELECT
          else
             if(cart2d) then
                SELECT CASE (nChild)
                CASE (1)    
                   subcellCentre = thisOctal%centre + d * VECTOR(-1.d0,0.d0,-1.d0)
                CASE (2)    
                   subcellCentre = thisOctal%centre + d * VECTOR(1.d0,0.d0,-1.d0)
                CASE (3)    
                   subcellCentre = thisOctal%centre + d * VECTOR(-1.d0,0.d0,-1.d0)
                CASE (4)    
                   subcellCentre = thisOctal%centre + d * VECTOR(1.d0,0.d0,-1.d0)
                CASE DEFAULT
                   PRINT *, "Error:: Invalid nChild passed to subcellCentre threed case"
                   PRINT *, "        nChild = ", nChild 
                   STOP
                END SELECT
             else
                SELECT CASE (nChild)
                CASE (1)    
                   subcellCentre = rVec - (d * xHat) - (d * zHat)
                   subcellCentre = rotateZ(subcellCentre,-thisOctal%phi)
                CASE (2)    
                   subcellCentre = rVec + (d * xHat) - (d * zHat)
                   subcellCentre = rotateZ(subcellCentre,-thisOctal%phi)
                CASE (3)    
                   subcellCentre = rVec - (d * xHat) + (d * zHat)
                   subcellCentre = rotateZ(subcellCentre,-thisOctal%phi)
                CASE (4)    
                   subcellCentre = rVec + (d * xHat) + (d * zHat)
                   subcellCentre = rotateZ(subcellCentre,-thisOctal%phi)
                CASE DEFAULT
                   PRINT *, "Error:: Invalid nChild passed to subcellCentre cylindical 3D case 2"
                   PRINT *, "        nChild = ", nChild 
                   do
                   enddo
                END SELECT
             end if
          endif
       endif
    else
    SELECT CASE (nChild)
    CASE (1)    
       subcellCentre = thisOctal%centre - (d * xHat) - (d * zHat)
    CASE (2)    
       subcellCentre = thisOctal%centre + (d * xHat) - (d * zHat)
    CASE (3)    
       subcellCentre = thisOctal%centre - (d * xHat) + (d * zHat)
    CASE (4)    
       subcellCentre = thisOctal%centre + (d * xHat) + (d * zHat)
    CASE DEFAULT
       PRINT *, "Error:: Invalid nChild passed to subcellCentre twoD case 3"
       PRINT *, "        nChild = ", nChild 
       do;enddo
       END SELECT
    endif

666 continue
  END FUNCTION subcellCentre

  subroutine subcellCorners(thisOctal,nChild, corner)
    ! returns the centre of one of the subcells of the current octal 

    IMPLICIT NONE

    TYPE(octal), INTENT(IN) :: thisOctal 
    INTEGER, INTENT(IN)     :: nChild    ! index (1-8) of the subcell
    type(VECTOR) :: rVec, corner(8), cellCentre
    real(oct)    :: d, phi, zp, zm, r1, r2, phistart, phiend, dphi
    corner = subcellCentre(thisOctal, nchild)

    d = thisOctal%subcellSize * 0.5_oc

    if (thisOctal%oneD) then
       select case (nchild)
          case(1)
             cellCentre = thisOctal%centre - d * xHat
          case(2)
             cellCentre = thisOctal%centre + d * xHat
          case DEFAULT
             write(*,*) "bug - one-d cell has more than 3 children"
       end select
       goto 666
    endif

    if (thisOctal%threeD) then  !do the three-d case as per diagram
       if (.not.thisOctal%cylindrical) then

          SELECT CASE (nChild)
          CASE (1)    
             cellCentre = thisOctal%centre + d * VECTOR(-1.d0,-1.d0,-1.d0)
          CASE (2)    
             cellCentre = thisOctal%centre + d * VECTOR(1.d0,-1.d0,-1.d0)
          CASE (3)    
             cellCentre = thisOctal%centre + d * VECTOR(-1.d0,1.d0,-1.d0)
          CASE (4)    
             cellCentre = thisOctal%centre + d * VECTOR(1.d0,1.d0,-1.d0)
          CASE (5)    
             cellCentre = thisOctal%centre + d * VECTOR(-1.d0,-1.d0,1.d0)
          CASE (6)    
             cellCentre = thisOctal%centre + d * VECTOR(1.d0,-1.d0,1.d0)
          CASE (7)    
             cellCentre = thisOctal%centre + d * VECTOR(-1.d0,1.d0,1.d0)
          CASE (8)    
             cellCentre = thisOctal%centre + d * VECTOR(1.d0,1.d0,1.d0)
          CASE DEFAULT
             PRINT *, "Error:: Invalid nChild passed to cellCentre threed case"
             PRINT *, "        nChild = ", nChild 
             STOP
          END SELECT
       else
          rVec = subcellCentre(thisOctal, nChild)
          d = thisOctal%subcellSize/2.d0
          zp = rVec%z + d
          zm = rVec%z - d
          r1 = sqrt(rVec%x**2 + rVec%y**2) - d
          r2 = sqrt(rVec%x**2 + rVec%y**2) + d
          phi = atan2(rVec%y, rVec%x)
          if (phi  < 0.d0) phi = phi + twoPi
          dphi = returndPhi(thisOctal)
          phiStart = phi - dphi
          if (phiStart  < 0.d0) phiStart = phiStart + twoPi
          phiEnd = phi + dphi
          if (phiEnd < 0.d0) phiEnd = phiEnd + twopi
          corner(1) = VECTOR(r1*cos(phiStart), r1*sin(phiStart),zp)
          corner(2) = VECTOR(r2*cos(phiStart), r2*sin(phiStart),zp)
          corner(3) = VECTOR(r1*cos(phiEnd), r1*sin(phiEnd),zp)
          corner(4) = VECTOR(r2*cos(phiEnd), r2*sin(phiEnd),zp)
          corner(5) = VECTOR(r1*cos(phiStart), r1*sin(phiStart),zm)
          corner(6) = VECTOR(r2*cos(phiStart), r2*sin(phiStart),zm)
          corner(7) = VECTOR(r1*cos(phiEnd), r1*sin(phiEnd),zm)
          corner(8) = VECTOR(r2*cos(phiEnd), r2*sin(phiEnd),zm)
       endif
    else
       SELECT CASE (nChild)
       CASE (1)    
          cellCentre = thisOctal%centre - (d * xHat) - (d * zHat)
       CASE (2)    
          cellCentre = thisOctal%centre + (d * xHat) - (d * zHat)
       CASE (3)    
          cellCentre = thisOctal%centre - (d * xHat) + (d * zHat)
       CASE (4)    
          cellCentre = thisOctal%centre + (d * xHat) + (d * zHat)
       CASE DEFAULT
          PRINT *, "Error:: Invalid nChild passed to cellCentre twoD case 3"
          PRINT *, "        nChild = ", nChild 
          do;enddo
       END SELECT
       
       corner(1) = cellCentre + d * xHat + d * zHat
       corner(2) = cellCentre + d * xHat - d * zHat
       corner(3) = cellCentre - d * xHat + d * zHat
       corner(4) = cellCentre - d * xHat - d * zHat

       
    endif
666 continue
  END subroutine subcellCorners
  
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
    use constants_mod, only: twoPi
    implicit none
    logical :: out
    type(octal), intent(in) :: this
    integer, intent(in) :: subcell   
    real(double), intent(in) :: x, y, z
    !
    TYPE(Vector)     :: cellCenter
    real(double) :: x0, y0, z0  ! cell center
    real(double) :: d, dp, dm, r, phi, r0, phi0, dphi
    real(double), parameter :: eps = 0.0d0
    
! Either use a subcell or the whole octal if subcell=0
    if ( subcell > 0 ) then 
       d = (this%subcellSize)*0.5d0
       cellCenter = subcellCentre(this,subcell)
    else
       d = this%subcellSize
       cellCenter = this%centre
    end if

    dp = d+eps
    dm = d-eps

    if (this%oneD) then
       write(*,*) "one-d case not implemented in within_subcell"
       stop
    endif
    
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

          dphi = returndPhi(this)

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
!    use inputs_mod, only: splitOverMPI
    use constants_mod, only: fourPi, twoPi, pi
    type(OCTAL) :: thisOctal
    integer :: subcell
    real(double) :: r1, r2, dphi
    type(VECTOR) :: rVec

    if (thisOctal%oneD) then
!       if(thisOctal%mpiThread(subcell) == 0) then
!          !not domain decomposed
          rVec = subcellCentre(thisOctal, subcell)
          r1 = rVec%x - thisOctal%subcellSize/2.d0
          r2 = rVec%x + thisOctal%subcellSize/2.d0
          v = (fourPi / 3.d0) * (r2**3-r1**3)
!       else
!          v = thisOctal%subcellSize**3
!          rVec = subcellCentre(thisOctal, subcell)
!          r1 = rVec%x - thisOctal%subcellSize/2.d0
!          r2 = rVec%x + thisOctal%subcellSize/2.d0
!          v = (r2**3-r1**3)
!       end if
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

  subroutine allocateAttributeDouble(array, nSize)
    integer :: nSize
    real(double), pointer :: array(:)

    if (.not.associated(array)) then
       allocate(array(1:nSize))
       array = TINY(array)
    endif
  end subroutine allocateAttributeDouble

  subroutine allocateAttributeDouble2d(array, nSize1, nSize2)
    integer :: nSize1, nSize2
    real(double), pointer :: array(:,:)

    if (.not.associated(array)) then
       allocate(array(1:nSize1, 1:nSize2))
       array = TINY(array)
    endif
  end subroutine allocateAttributeDouble2d

  subroutine allocateAttributeDouble3d(array, nSize1, nSize2, nsize3)
    integer :: nSize1, nSize2, nSize3
    real(double), pointer :: array(:,:,:)

    if (.not.associated(array)) then
       allocate(array(1:nSize1, 1:nSize2, 1:nSize3))
       array = TINY(array)
    endif
  end subroutine allocateAttributeDouble3d
  
  subroutine allocateAttributeReal(array, nSize)
    integer :: nSize
    real, pointer :: array(:)

    if (.not.associated(array)) then
       allocate(array(1:nSize))
       array = TINY(array)
    endif
  end subroutine allocateAttributeReal

  subroutine allocateAttributeInteger(array, nSize)
    integer :: nSize
    integer, pointer :: array(:)

    if (.not.associated(array)) then
       allocate(array(1:nSize))
       array = 0
    endif
  end subroutine allocateAttributeInteger

  subroutine allocateAttributeLogical(array, nSize)
    integer :: nSize
    logical, pointer :: array(:)

    if (.not.associated(array)) then
       allocate(array(1:nSize))
       array = .false.
    endif
  end subroutine allocateAttributeLogical

  subroutine allocateAttributeVector(array, nSize)
    integer :: nSize
    type(VECTOR), pointer :: array(:)

    if (.not.associated(array)) then
       allocate(array(1:nSize))
       array = zeroVector
    endif
  end subroutine allocateAttributeVector


  subroutine copyAttributeSingleInteger(dest, source)
    integer :: dest, source
    dest = source
  end subroutine copyAttributeSingleInteger

  subroutine copyAttributeDoublePointer(dest, source)
    real(double), pointer :: dest(:), source(:)

    if (associated(source)) then
       allocate(dest(SIZE(source)))
       dest = source
    endif

  end subroutine copyAttributeDoublePointer

  subroutine copyAttributeDoublePointer2d(dest, source)
    real(double), pointer :: dest(:,:), source(:,:)

    if (associated(source)) then
       allocate(dest(SIZE(source, 1),SIZE(source,2)))
       dest = source
    endif

  end subroutine copyAttributeDoublePointer2d

  subroutine copyAttributeDoublePointer3d(dest, source)
    real(double), pointer :: dest(:,:,:), source(:,:,:)

    if (associated(source)) then
       allocate(dest(SIZE(source, 1),SIZE(source,2),SIZE(source,3)))
       dest = source
    endif

  end subroutine copyAttributeDoublePointer3d

  subroutine copyAttributeRealPointer(dest, source)
    real, pointer :: dest(:), source(:)

    if (associated(source)) then
       allocate(dest(SIZE(source)))
       dest = source
    endif

  end subroutine copyAttributeRealPointer

  subroutine copyAttributeIntegerPointer(dest, source)
    integer, pointer :: dest(:), source(:)

    if (associated(source)) then
       allocate(dest(SIZE(source)))
       dest = source
    endif

  end subroutine copyAttributeIntegerPointer


  subroutine copyAttributeLogicalPointer(dest, source)
    logical, pointer :: dest(:), source(:)

    if (associated(source)) then
       allocate(dest(SIZE(source)))
       dest = source
    endif

  end subroutine copyAttributeLogicalPointer

  subroutine copyAttributeVectorPointer(dest, source)
    type(VECTOR) , pointer :: dest(:), source(:)
    if (associated(source)) then
       allocate(dest(SIZE(source)))
       dest = source
    endif

  end subroutine copyAttributeVectorPointer

  subroutine deallocateAttributeDouble(array)
    real(double), pointer :: array(:)

    if (associated(array)) then
       deallocate(array)
       nullify(array)
    endif
  end subroutine deallocateAttributeDouble

  subroutine deallocateAttributeOctal(thisOctal)
    type(OCTAL), pointer :: thisOctal

    if (associated(thisOctal)) then
       deallocate(thisOctal)
       nullify(thisOctal)
    endif
  end subroutine deallocateAttributeOctal

  subroutine deallocateAttributeDouble2d(array)
    real(double), pointer :: array(:,:)

    if (associated(array)) then
       deallocate(array)
       nullify(array)
    endif
  end subroutine deallocateAttributeDouble2d

  subroutine deallocateAttributeReal2d(array)
    real, pointer :: array(:,:)

    if (associated(array)) then
       deallocate(array)
       nullify(array)
    endif
  end subroutine deallocateAttributeReal2d

  subroutine deallocateAttributeInteger2d(array)
    integer, pointer :: array(:,:)

    if (associated(array)) then
       deallocate(array)
       nullify(array)
    endif
  end subroutine deallocateAttributeInteger2d

  subroutine deallocateAttributeDouble3d(array)
    real(double), pointer :: array(:,:,:)

    if (associated(array)) then
       deallocate(array)
       nullify(array)
    endif
  end subroutine deallocateAttributeDouble3d


  subroutine deallocateAttributeVector(array)
    type(vector), pointer :: array(:)

    if (associated(array)) then
       deallocate(array)
       nullify(array)
    endif
  end subroutine deallocateAttributeVector

  subroutine deallocateAttributeReal(array)
    real, pointer :: array(:)

    if (associated(array)) then
       deallocate(array)
       nullify(array)
    endif
  end subroutine deallocateAttributeReal

  subroutine deallocateAttributeInteger(array)
    integer, pointer :: array(:)

    if (associated(array)) then
       deallocate(array)
       nullify(array)
    endif
  end subroutine deallocateAttributeInteger

  subroutine deallocateAttributeInteger3D(array)
    integer, pointer :: array(:,:,:)

    if (associated(array)) then
       deallocate(array)
       nullify(array)
    endif
  end subroutine deallocateAttributeInteger3D

  subroutine deallocateAttributeOctalPointer3D(array)
    type(OCTALPOINTER), pointer :: array(:,:,:)

    if (associated(array)) then
       deallocate(array)
       nullify(array)
    endif
  end subroutine deallocateAttributeOctalPointer3D

  subroutine deallocateAttributeLogical(array)
    logical, pointer :: array(:)

    if (associated(array)) then
       deallocate(array)
       nullify(array)
    endif
  end subroutine deallocateAttributeLogical



  subroutine readAttributeDoublePointer(lUnit, array, fileFormatted)
    integer :: lUnit
    integer :: nSize
    logical :: valPresent
    real(double), pointer :: array(:)
    logical :: fileFormatted

    if (fileFormatted) then
       read(lUnit,*) valPresent
    else
       read(lUnit) valPresent
    endif

    if (valPresent) then
       if (fileFormatted) then
          read(lunit,*) nSize
          allocate(array(1:nSize))
          read(lUnit,*) array(1:nSize)
       else
          read(lunit) nSize
          allocate(array(1:nSize))
          read(lUnit) array(1:nSize)
       endif
    endif

  end subroutine readAttributeDoublePointer

  subroutine readAttributeIntegerPointer(lUnit, array, fileFormatted)
    integer :: lUnit
    integer :: nSize
    logical :: valPresent
    integer, pointer :: array(:)
    logical :: fileFormatted
    
    if (fileFormatted) then
       read(lUnit,*) valPresent
    else
       read(lUnit) valPresent
    endif

    if (valPresent) then
       if (fileFormatted) then
          read(lunit,*) nSize
          allocate(array(1:nSize))
          read(lUnit,*) array(1:nSize)
       else
          read(lunit) nSize
          allocate(array(1:nSize))
          read(lUnit) array(1:nSize)
       endif
    endif

  end subroutine readAttributeIntegerPointer

  subroutine readAttributeRealPointer(lUnit, array, fileFormatted)
    integer :: lUnit
    integer :: nSize
    logical :: valPresent
    real, pointer :: array(:)
    logical :: fileFormatted
    
    if (fileFormatted) then
       read(lUnit,*) valPresent
    else
       read(lUnit) valPresent
    endif

    if (valPresent) then
       if (fileFormatted) then
          read(lunit,*) nSize
          allocate(array(1:nSize))
          read(lUnit,*) array(1:nSize)
       else
          read(lunit) nSize
          allocate(array(1:nSize))
          read(lUnit) array(1:nSize)
       endif
    endif

  end subroutine readAttributeRealPointer

  subroutine readAttributeVectorPointer(lUnit, array, fileFormatted)
    integer :: lUnit
    integer :: nSize
    logical :: valPresent
    TYPE(VECTOR) , pointer :: array(:)
    logical :: fileFormatted
    
    if (fileFormatted) then
       read(lUnit,*) valPresent
    else
       read(lUnit) valPresent
    endif

    if (valPresent) then
       if (fileFormatted) then
          read(lunit,*) nSize
          allocate(array(1:nSize))
          read(lUnit,*) array(1:nSize)
       else
          read(lunit) nSize
          allocate(array(1:nSize))
          read(lUnit) array(1:nSize)
       endif
    endif

  end subroutine readAttributeVectorPointer



  subroutine readAttributeStaticDouble(lUnit, array, fileFormatted)
    integer :: lUnit
    real(double) :: array(:)
    logical :: fileFormatted

    if (fileFormatted) then
       read(lUnit, *) array
    else
       read(lUnit) array
    endif
  end subroutine readAttributeStaticDouble

  subroutine readAttributeStaticVector(lUnit, array, fileFormatted)
    integer :: lUnit
    type(VECTOR) :: array(:)
    logical :: fileFormatted

    if (fileFormatted) then
       read(lUnit, *) array
    else
       read(lUnit) array
    endif
  end subroutine readAttributeStaticVector

  subroutine readAttributeStaticDoubleSingle(lUnit, array, fileFormatted)
    integer :: lUnit
    real(double) :: array
    logical :: fileFormatted

    if (fileFormatted) then
       read(lUnit, *) array
    else
       read(lUnit) array
    endif
  end subroutine readAttributeStaticDoubleSingle

  subroutine readAttributeStaticVectorSingle(lUnit, array, fileFormatted)
    integer :: lUnit
    type(VECTOR) :: array
    logical :: fileFormatted

    if (fileFormatted) then
       read(lUnit, *) array
    else
       read(lUnit) array
    endif
  end subroutine readAttributeStaticVectorSingle

  subroutine readAttributeStaticIntegerSingle(lUnit, array, fileFormatted)
    integer :: lUnit
    integer :: array
    logical :: fileFormatted

    if (fileFormatted) then
       read(lUnit, *) array
    else
       read(lUnit) array
    endif
  end subroutine readAttributeStaticIntegerSingle

  subroutine readAttributeStaticLogicalSingle(lUnit, array, fileFormatted)
    integer :: lUnit
    logical :: array
    logical :: fileFormatted

    if (fileFormatted) then
       read(lUnit, *) array
    else
       read(lUnit) array
    endif
  end subroutine readAttributeStaticLogicalSingle


  subroutine readAttributeStaticReal(lUnit, array, fileFormatted)
    integer :: lUnit
    real :: array(:)
    logical :: fileFormatted

    if (fileFormatted) then
       read(lUnit, *) array
    else
       read(lUnit) array
    endif
  end subroutine readAttributeStaticReal

  subroutine readAttributeStaticInteger(lUnit, array, fileFormatted)
    integer :: lUnit
    integer :: array(:)
    logical :: fileFormatted

    if (fileFormatted) then
       read(lUnit, *) array
    else
       read(lUnit) array
    endif
  end subroutine readAttributeStaticInteger

  subroutine readAttributeStaticLogical(lUnit, array, fileFormatted)
    integer :: lUnit
    logical :: array(:)
    logical :: fileFormatted

    if (fileFormatted) then
       read(lUnit, *) array
    else
       read(lUnit) array
    endif
  end subroutine readAttributeStaticLogical


  subroutine writeAttributeDoublePointer(lUnit, array, fileFormatted)
    integer :: lUnit
    real(double), pointer :: array(:)
    logical :: fileFormatted

    if (associated(array)) then
       if (fileFormatted) then
          write(lUnit,*) .true.
       else
          write(lUnit) .true.
       endif

       if (fileFormatted) then
          write(lunit,*) SIZE(array)
          write(lUnit,*) array(1:SIZE(array))
       else
          write(lunit) SIZE(array)
          write(lUnit) array(1:SIZE(array))
       endif
    else
       if (fileFormatted) then
          write(lUnit,*) .false.
       else
          write(lUnit) .false.
       endif
    endif
  end subroutine writeAttributeDoublePointer

  subroutine writeAttributeRealPointer(lUnit, array, fileFormatted)
    integer :: lUnit
    real, pointer :: array(:)
    logical :: fileFormatted

    if (associated(array)) then
       if (fileFormatted) then
          write(lUnit,*) .true.
       else
          write(lUnit) .true.
       endif

       if (fileFormatted) then
          write(lunit,*) SIZE(array)
          write(lUnit,*) array(1:SIZE(array))
       else
          write(lunit) SIZE(array)
          write(lUnit) array(1:SIZE(array))
       endif
    else
       if (fileFormatted) then
          write(lUnit,*) .false.
       else
          write(lUnit) .false.
       endif
    endif
  end subroutine writeAttributeRealPointer

  subroutine writeAttributeIntegerPointer(lUnit, array, fileFormatted)
    integer :: lUnit
    integer, pointer :: array(:)
    logical :: fileFormatted

    if (associated(array)) then
       if (fileFormatted) then
          write(lUnit,*) .true.
       else
          write(lUnit) .true.
       endif

       if (fileFormatted) then
          write(lunit,*) SIZE(array)
          write(lUnit,*) array(1:SIZE(array))
       else
          write(lunit) SIZE(array)
          write(lUnit) array(1:SIZE(array))
       endif
    else
       if (fileFormatted) then
          write(lUnit,*) .false.
       else
          write(lUnit) .false.
       endif
    endif
  end subroutine writeAttributeIntegerPointer

  subroutine writeAttributeLogicalPointer(lUnit, array, fileFormatted)
    integer :: lUnit
    logical, pointer :: array(:)
    logical :: fileFormatted

    if (associated(array)) then
       if (fileFormatted) then
          write(lUnit,*) .true.
       else
          write(lUnit) .true.
       endif

       if (fileFormatted) then
          write(lunit,*) SIZE(array)
          write(lUnit,*) array(1:SIZE(array))
       else
          write(lunit) SIZE(array)
          write(lUnit) array(1:SIZE(array))
       endif
    else
       if (fileFormatted) then
          write(lUnit,*) .false.
       else
          write(lUnit) .false.
       endif
    endif
  end subroutine writeAttributeLogicalPointer


  subroutine writeAttributeVectorPointer(lUnit, array, fileFormatted)
    integer :: lUnit
    type(VECTOR), pointer :: array(:)
    logical :: fileFormatted

    if (associated(array)) then
       if (fileFormatted) then
          write(lUnit,*) .true.
       else
          write(lUnit) .true.
       endif

       if (fileFormatted) then
          write(lunit,*) SIZE(array)
          write(lUnit,*) array(1:SIZE(array))
       else
          write(lunit) SIZE(array)
          write(lUnit) array(1:SIZE(array))
       endif
    else
       if (fileFormatted) then
          write(lUnit,*) .false.
       else
          write(lUnit) .false.
       endif
    endif
  end subroutine writeAttributeVectorPointer


  subroutine writeAttributeStaticDouble(lUnit, array, fileFormatted)
    integer :: lUnit
    real(double):: array(:)
    logical :: fileFormatted

    if (fileFormatted) then
       write(lUnit,*) array(1:SIZE(array))
    else
       write(lUnit) array(1:SIZE(array))
    endif
  end subroutine writeAttributeStaticDouble

  subroutine writeAttributeStaticVector(lUnit, array, fileFormatted)
    integer :: lUnit
    type(VECTOR) :: array(:)
    logical :: fileFormatted

    if (fileFormatted) then
       write(lUnit,*) array(1:SIZE(array))
    else
       write(lUnit) array(1:SIZE(array))
    endif
  end subroutine writeAttributeStaticVector

 subroutine writeAttributeStaticReal(lUnit, array, fileFormatted)
    integer :: lUnit
    real :: array(:)
    logical :: fileFormatted

    if (fileFormatted) then
       write(lUnit,*) array(1:SIZE(array))
    else
       write(lUnit) array(1:SIZE(array))
    endif
  end subroutine writeAttributeStaticReal

 subroutine writeAttributeStaticInteger(lUnit, array, fileFormatted)
    integer :: lUnit
    integer :: array(:)
    logical :: fileFormatted

    if (fileFormatted) then
       write(lUnit,*) array(1:SIZE(array))
    else
       write(lUnit) array(1:SIZE(array))
    endif
  end subroutine writeAttributeStaticInteger

 subroutine writeAttributeStaticIntegerSingle(lUnit, array, fileFormatted)
    integer :: lUnit
    integer :: array
    logical :: fileFormatted

    if (fileFormatted) then
       write(lUnit,*) array
    else
       write(lUnit) array
    endif
  end subroutine writeAttributeStaticIntegerSingle

 subroutine writeAttributeStaticLogicalSingle(lUnit, array, fileFormatted)
    integer :: lUnit
    logical :: array
    logical :: fileFormatted

    if (fileFormatted) then
       write(lUnit,*) array
    else
       write(lUnit) array
    endif
  end subroutine writeAttributeStaticLogicalSingle

 subroutine writeAttributeStaticVectorSingle(lUnit, array, fileFormatted)
    integer :: lUnit
    type(VECTOR) :: array
    logical :: fileFormatted

    if (fileFormatted) then
       write(lUnit,*) array
    else
       write(lUnit) array
    endif
  end subroutine writeAttributeStaticVectorSingle


 subroutine writeAttributeStaticDoubleSingle(lUnit, array, fileFormatted)
    integer :: lUnit
    real(double) :: array
    logical :: fileFormatted

    if (fileFormatted) then
       write(lUnit,*) array
    else
       write(lUnit) array
    endif
  end subroutine writeAttributeStaticDoubleSingle


 subroutine writeAttributeStaticLogical(lUnit, array, fileFormatted)
    integer :: lUnit
    logical :: array(:)
    logical :: fileFormatted

    if (fileFormatted) then
       write(lUnit,*) array(1:SIZE(array))
    else
       write(lUnit) array(1:SIZE(array))
    endif
  end subroutine writeAttributeStaticLogical



END MODULE octal_mod
