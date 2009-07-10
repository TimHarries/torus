module disc_class
  
  use kind_mod
  use vector_mod
  use constants_mod
  use octal_mod, only: OCTAL, subcellCentre
  use gridtype_mod, only: GRIDTYPE

  !  
  ! Class definition for a accreation disc (not including the central object)
  !
  ! Created: 05-nov-2003  (Ryuich Kurosawa)
 
  public::  &
       &   new, &
       &   get_parameters, &
       &   component, &
       &   alpha_disc_density, &
       &   ave_alpha_disc_density, &
       &   add_alpha_disc, &
       &   finish_grid, &
       &   turn_on_disc, &
       &   turn_off_disc, &
       &   in_alpha_disc

  private::  &
       &   int_alpha_disc, &
       &   scale_size, &
       &   need_to_split, &
       &   add_new_children, &
       &   assign_values, &
       &   keplerian_velocity
 


  !----------------------------------------------------------
  type alpha_disc
     private 
     ! this is for "alpha-disc" standard model 
     real(double) :: Rh  ! [10^10 cm] inner radius of the disc (hole) 
     real(double) :: Rd  ! [10^10 cm] outer radius of the disc 
     !
     real(double) :: Md  ! mass of the disc in grams
     real(double) :: x0, y0, z0  ! coordinates of the center [10^10cm]
     real(double) :: sx, sy, sz  ! components of spin vector [units?]
     !
     real(double) :: Mcore ! mass of the cnetral object in [Msun]
  end type alpha_disc

  

  ! overload functions to make interface here ----------
  interface new
     module procedure int_alpha_disc_default
     module procedure int_alpha_disc
  end interface


  real, private, parameter :: rho_min = 1.e-19


contains
  ! =====================================================================
  ! constructors
  !======================================================================

  ! initializing the alpha_disc with default parameters
  subroutine int_alpha_disc_default(this)
    implicit none
    type(alpha_disc), intent(inout) :: this 

    this%Rh = 6.0    *6.95      ! 6 R_sun [10^10cm]
    this%Rd = 100.0  *1.5d3     ! 100 AU  [10^10cm]
    this%Md = 0.005  *1.989d33  ! 1/500 M_sun  [g]
    this%x0 = 0.0;  this%y0 = 0.0;  this%z0=1.0
    ! pointing z direction.  
    this%sx = 0.0;  this%sy = 0.0;  this%sz=1.0
    this%Mcore = 1.0d0 ! [Msun]

  end subroutine int_alpha_disc_default
  
  
  ! inilializing with prarameters
  subroutine int_alpha_disc(this, Rh, Rd, Md, x0, y0, z0, sx, sy, sz, Mcore) 
    implicit none 
    
    type(alpha_disc), intent(inout) :: this 
    ! this is for "alpha-disc" standard model 
    real(double), intent(in) :: Rh  ! [10^10 cm] inner radius of the disc (hole) 
    real(double), intent(in) :: Rd  ! [10^10 cm] outer radius of the disc 
    !
    real(double), intent(in) :: Md  ! mass of the disc in grams
    real(double), intent(in) :: x0, y0, z0  ! coordinates of the center [10^10cm]
    real(double), intent(in) :: sx,sy, sz  ! components of spin vector ! [units?]
    real(double), intent(in) :: Mcore ! mass of the cnetral object in [Msun]


    this%Rh = Rh  ! [10^10cm]
    this%Rd = Rd  ! [10^10cm]
    this%Md = Md  ! [g]
    this%x0 = x0;  this%y0 = y0;  this%z0=z0
    ! pointing z direction.  
    this%sx = sx;  this%sy = sy;  this%sz=sz
    this%Mcore = Mcore  ! [Msun]

  end subroutine int_alpha_disc



  !==========================================================================
  ! Accessors
  !==========================================================================

  ! Given a alpha_disc object, it returns all components of the objects.
  subroutine get_parameters(this, Rh, Rd, Md, x0, y0, z0, sx, sy, sz, Mcore) 
    implicit none 
    
    type(alpha_disc), intent(in) :: this 
    ! this is for "alpha-disc" standard model 
    real(double), intent(out) :: Rh  ! [10^10 cm] inner radius of the disc (hole) 
    real(double), intent(out) :: Rd  ! [10^10 cm] outer radius of the disc 
    !
    real(double), intent(out) :: Md  ! mass of the disc in grams
    real(double), intent(out) :: x0, y0, z0  ! coordinates of the center [10^10cm]
    real(double), intent(out) :: sx,sy, sz   ! components of spin vector ! [units?]
    real(double), intent(out) :: Mcore ! mass of the cnetral object in [Msun]

    Rh = this%Rh
    Rd = this%Rd
    Md = this%Md
    x0 = this%x0; y0 = this%y0; z0 = this%z0
    sx = this%sx; sy = this%sy; sz = this%sz
    Mcore = this%Mcore
    
  end subroutine get_parameters


  ! Given a alpha_disc object and the name of its component, this function 
  ! returns the component
  real(double) function component(this, name) 
    implicit none 
    
    type(alpha_disc), intent(in) :: this 
    character(LEN=*), intent(in) :: name 
    
    select case (name)
    case ("Rh")
       component = this%Rh
    case ("Rd")
       component = this%Rd
    case ("Md")
       component = this%Md
    case ("x0")
       component = this%x0
    case ("y0")
       component = this%y0
    case ("z0")
       component = this%z0
    case ("sx")
       component = this%sx
    case ("sy")
       component = this%sy
    case ("sz")
       component = this%sx
    case ("Mcore")
       component = this%Mcore
    case default
       write(*,*) "Error:: Unknown NAME was passed to [disc_class::component]."
       write(*,*) "Exiting the program because of this error!"
       stop
    end select

  end function component


  !=============================================================================
  ! Tools
  !=============================================================================


  ! Given a position (xpos, ypos and zpos) in 10^10 cm, this routine
  ! will return the density [g/cm^3] and the length scale (sizescale).
  !
  ! Note: This routine has been slightly changed from the original 
  !       one (discdensity.f) created by Matthew Bates. 
  !
  ! Matthew's comments:
  !     Adds discs around sink particles
  !     Disc is standard solution Sigma=R^-3/4, H=R^9/8 with H/R=0.05 at Rmax
  !     Returns density at a point in parameter space (0 if not within a disc)
  !     and gives a characteristic sizescale for grid construction.
  
  function alpha_disc_density(this, xpos, ypos, zpos, sizescale) RESULT(density_out)
    implicit none
    real(double) :: density_out  ! in [g/cm^3]
    type(alpha_disc), intent(in) :: this
    real(double), intent(in) :: xpos,ypos,zpos  ! should be in [10^10cm]
    real(double), intent(out):: sizescale       ! should be in [10^10cm]
    !
    real(double), parameter :: pi = 3.14159265d0
    real(double) :: holesize
    real(double) :: discrad, discmass, unitx, unity, unitz
    real(double) :: spinx, spiny, spinz,  x, y, z
    real(double) :: xposi, yposi, zposi, rad, radius2, radius
    real(double) :: unorm, zdist, height, densitymidplane
    ! maximum density allowed
    real(double), parameter :: rho_max = 2.5d50   ! g/cm^3 
!    real(double), parameter :: rho_min = 1.0d-25  ! g/cm^3 
    !
    !
!    real(double), parameter :: rmax = 150.0d00  ! 10^10 cm
    real(double) :: rmax  ! 10^10 cm
    !
    real(double) :: q, alpha, beta, w !, p


    rmax = this%Rd
        
    ! these should be parameterized and put them in this object.
    alpha = 0.75d0   ! used for most of T Tau Halpha paper
    beta = 1.125d0   ! used for most of T Tau Halpha paper

!    alpha = 0.5d0   !  for a flared disc
!    beta = 1.3d0    !  for a flared disc

    !
    !
    holesize = this%Rh  ! [10^10cm]

    !
    density_out = rho_min
    sizescale = 1.0d+10

    ! retrive some disc parameters
    discrad = this%Rd
    discmass = this%Md
    spinx = this%sx
    spiny = this%sy
    spinz = this%sz
    
    ! poistion of the center
    x = this%x0
    y = this%y0
    z = this%z0
       
    xposi = xpos-x
    yposi = ypos-y
    zposi = zpos-z
    w =  xposi**2 + yposi**2
    radius2 = w + zposi**2
    w = SQRT(w)
    radius = SQRT(radius2)

!    if (radius > discrad/2.0d0) beta=4.0d0 ! flared for testing

    if (radius < discrad) then 
       q = SQRT(xposi**2 + yposi**2)
       if (q < 1.1d0*holesize .and. radius < 1.1d0*holesize) then
          sizescale = holesize
          if (q > (holesize/2.0d0)) sizescale=sizescale/2.0d0
       else
          rad = radius/discrad
          !  Distance above disc plane is obtained by doing a dot product with the
          !     unit vector giving the disc orientation
          unitx = spinx
          unity = spiny
          unitz = spinz
          unorm = SQRT(unitx**2 + unity**2 + unitz**2)
          zdist = (xposi*unitx + yposi*unity + zposi*unitz)/unorm
          height = 0.05d0*discrad * (rad**beta)

!          densitymidplane = 5.0d0*discmass/(8.0d0*pi*discrad**2)*  &
!               rad**(-0.75d0)
          densitymidplane = 5.0d0*discmass/(8.0d0*pi*discrad**2)*  &
               rad**(-alpha)

          density_out = densitymidplane*   &
               EXP(-zdist**2/(2.0d0*height**2))/(SQRT(2.0d0*pi)*height)

          
          ! now finding the scale size at (xpos, ypos, zpos)

!          if (density_out>rho_min) then

!             p = abs(zdist/height)
             
!             p = 0.2d0
!             sizescale = (height) * (w/this%Rh)**p

             sizescale = height ! (used until 24-may-05 with small ISM-like grains)
!             sizescale = 5.0d0*height ! testing for flared disc
                
!              ! limit the smallest and largest sizes
              sizescale = MAX(sizescale, 1.0d-1)
!              sizescale = MIN(sizescale, 1.0d02)                

             ! else takes the default density from the top.
!          end if
       endif
       !  given position can only be inside 1 disc at most (discs cannot overlap) 
       goto 10
    endif
10  continue    

    ! converting units from g/(10^10cm)^3 to g/cm^3
    density_out = density_out/1.0d30   ! [g/cm^3]

    ! cutoff density
    ! max and min  density allowed
    if (density_out>rho_max) density_out = rho_max ! [g/cm^3]
    if (density_out<rho_min) density_out = rho_min ! [g/cm^3]

!    if (radius > rmax) density_out = rho_min ! [g/cm^3]
    if (w > rmax) density_out = rho_min ! [g/cm^3]


  end function alpha_disc_density


  

  !
  ! Find the scale size of a disc for given sizes of 
  ! minimum and maximum size scales allowed (ds_min and ds_max), 
  ! and minumum and maximum distance (e.g. radius or height above the disk..), i.e.
  ! s_min and s_max.
  function scale_size(s, s_min, s_max, ds_min, ds_max, logscale) RESULT(out)
    implicit none
    real(double) :: out
    ! distance from center of a star or from the disk plane
    real(double), intent(in) :: s,  s_min, s_max, ds_min, ds_max
    logical, intent(in) :: logscale
    real(double)  :: p,ss

    ss=s
    if (logscale) then
       if (s<s_min) ss =s_min
       if (s>s_max) ss =s_max    
       p = LOG10(ds_max/ds_min) / (s_max - s_min)
       out = ds_min * 10.0d0**(p*(ss-s_min))
    else ! linear scale
       if (s<s_min) ss =s_min
       if (s>s_max) ss =s_max   
       p = (ds_max-ds_min) / (s_max - s_min)
       out = ds_min +  p*(ss-s_min)       
    end if

  end function scale_size





  !
  !  Given a grid object, this routine will add extra grid and assign density (and etc) 
  !  to the grid using the density function defined in this module.
  !  


  RECURSIVE SUBROUTINE add_alpha_disc(thisOctal,grid,this)
    ! uses an external function to decide whether to split a subcell of
    !   the current octal. 

    IMPLICIT NONE

    TYPE(OCTAL), POINTER :: thisOctal
      ! 'limitScalar' is the value the decideSplit function uses to
      !   decide whether or not to split cell.
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid through to the 
    type(alpha_disc), intent(in)  :: this !
    
    !
    TYPE(OCTAL), POINTER :: childPointer  
    INTEGER              :: subcell, i    ! loop counters
    logical :: splitThis
    TYPE(VECTOR)     :: cellCentre 
    integer :: n

    splitThis = .false.
    subcell = 1
    
    if (thisOctal%threeD) then
       n = 8
    else
       n = 4
    endif


    if (thisOctal%hasChild(subcell)) then
       do i = 1, n
          childPointer => thisOctal%child(i)
          CALL add_alpha_disc(childPointer,grid,this)
       end do
    else   
       ! get the size and the position of the centre of the current cell
       cellCentre = subcellCentre(thisOctal, subcell)

       if (modulus(cellCentre) > this%Rh) then 
          do while ((.not.splitThis).and.(subcell <= n))
             IF (need_to_split(thisOctal,subcell,this)) splitThis = .true.           
             subcell = subcell + 1
          end do
          
          if (splitThis) then
             ! using a routine in this module
             CALL add_new_children(thisOctal, grid, this)
             ! find the index of the new child and decend the tree
             do i = 1, n
                childPointer => thisOctal%child(i)
                CALL add_alpha_disc(childPointer,grid,this)
             end do
          end if
       end if
    end IF

    
   END SUBROUTINE add_alpha_disc


  ! This is based decideSplit in amr_mod.f90

  logical function need_to_split(thisOctal,subcell, this)

    IMPLICIT NONE

    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    type(alpha_disc), intent(in) :: this
    
    real(oct)  :: cellSize, x, y, z, w, h
    real(double) :: rho_disc, scale_length
    TYPE(VECTOR)     :: cellCentre 
!    logical:: close_to_disc

    ! get the size and the position of the centre of the current cell
    cellSize = thisOctal%subcellSize*2.0d0
    cellCentre = subcellCentre(thisOctal, subcell)
    x = cellcentre%x
    y = cellcentre%y
    z = cellcentre%z
    w = SQRT(x*x+y*y)
    h = ABS(z)

!    close_to_disc = .false.
!    if (h < 10000.d0) close_to_disc =.true.

!    if (.not. close_to_disc .and. w > this%Rh .and. ABS(x-ThisOctal%subcellSize/2.0d0)>1.0d-3) then
       ! Splits the two cells closer to the disc
!       if (h < cellSize/2.0d0) then
!          need_to_split = .true.
!       end if
!    else
       if (in_alpha_disc(this, cellCentre) )then
          
          rho_disc = alpha_disc_density(this, x, y, z, scale_length)        
          if (cellsize > scale_length) then       
             !       if (cellsize >   scale_size(w,this%Rh, this%Rd, 1.0d0, 1.0d4, .true.)  ) then       
             need_to_split = .true.  ! units in 10^10cm
          else
             need_to_split = .false.
          end if
       else
          need_to_split = .false.
       end if
       
!    end if


  end function need_to_split  
  !
  !
  !
  SUBROUTINE add_new_children(parent, grid, this)
    ! adds all eight new children to an octal

    IMPLICIT NONE
    
    TYPE(octal), POINTER :: parent     ! pointer to the parent octal 
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid to routines that
                                          !   calculate the variables stored in
                                          !   the tree.
    TYPE(alpha_disc), INTENT(IN)  :: this
    INTEGER       :: subcell           ! loop counter
    INTEGER, SAVE :: counter = 9       ! keeps track of the current subcell label
                                       ! - this isn't very clever. might change it. 
    INTEGER       :: newChildIndex     ! the storage location for the new child
    INTEGER       :: error

!    real(double) :: rho
    TYPE(OCTAL), POINTER :: childPointer      
    
    if (any(parent%hasChild(1:8))) write(*,*) "parent already has a child"

    if (parent%threeD) then
       parent%nChildren = 8
       parent%maxChildren = 8
    else
       parent%nChildren = 4
       parent%maxChildren = 4
    endif

    if (parent%threed) then
       ALLOCATE(parent%child(1:8), STAT=error)
    else
       ALLOCATE(parent%child(1:4), STAT=error)
    endif
    IF ( error /= 0 ) THEN
       PRINT *, 'Panic: allocation failed.'
       STOP
    END IF



    ! update the parent octal
    parent%hasChild(1:parent%maxChildren) = .TRUE.
    parent%indexChild(1) = 1
    parent%indexChild(2) = 2
    parent%indexChild(3) = 3
    parent%indexChild(4) = 4
    if (parent%threeD) then
       parent%indexChild(5) = 5
       parent%indexChild(6) = 6
       parent%indexChild(7) = 7
       parent%indexChild(8) = 8
    endif

    do newChildIndex = 1, parent%maxChildren

       ! allocate any variables that need to be  
       ALLOCATE(parent%child(newChildIndex)%kappaAbs(8,grid%nOpacity))
       ALLOCATE(parent%child(newChildIndex)%kappaSca(8,grid%nOpacity))
       parent%child(newChildIndex)%kappaAbs(:,:) = 1.0e-30
       parent%child(newChildIndex)%kappaSca(:,:) = 1.0e-30

       ALLOCATE(parent%child(newChildIndex)%N(8,grid%maxLevels))
       NULLIFY(parent%child(newChildIndex)%child)

       ! set up the new child's variables
       parent%child(newChildIndex)%threed = parent%threed
       parent%child(newChildIndex)%twoD = parent%twod
       parent%child(newChildIndex)%maxChildren = parent%maxChildren
       parent%child(newChildIndex)%parent => parent
       parent%child(newChildIndex)%subcellSize = parent%subcellSize / 2.0_oc
       parent%child(newChildIndex)%hasChild = .false.
       parent%child(newChildIndex)%nChildren = 0
       parent%child(newChildIndex)%indexChild = -999 ! values are undefined
       parent%child(newChildIndex)%nDepth = parent%nDepth + 1
       parent%child(newChildIndex)%centre = subcellCentre(parent,newChildIndex)
       parent%child(newChildIndex)%probDistLine = 0.0
       parent%child(newChildIndex)%probDistCont = 0.0
       parent%child(newChildIndex)%biasLine3D = 1.0 
       parent%child(newChildIndex)%biasCont3D = 1.0 
       parent%child(newChildIndex)%velocity = vector(1.e-30,1.e-30,1.e-30)
       parent%child(newChildIndex)%cornerVelocity = vector(1.e-30,1.e-30,1.e-30)
       parent%child(newChildIndex)%chiLine = 1.e-30
       parent%child(newChildIndex)%etaLine = 1.e-30
       parent%child(newChildIndex)%etaCont = 1.e-30
       parent%child(newChildIndex)%N = 1.e-10
       parent%child(newChildIndex)%Ne = 1.e-10
       parent%child(newChildIndex)%inFlow(:) = .false.
       parent%child(newChildIndex)%temperature = 3.0
       parent%child(newChildIndex)%changed = .false.
       
       ! put some data in the eight subcells of the new child
       DO subcell = 1, parent%maxChildren
          childPointer => parent%child(newChildIndex)
!          rho =ave_alpha_disc_density(childPointer, subcell, this)
!          parent%child(newChildIndex)%rho(subcell) = rho          
          parent%child(newChildIndex)%label(subcell) = counter
          counter = counter + 1
          if (subcell == parent%maxChildren) &
               call fill_velocity_corners(this, childPointer)          

       END DO
 
    enddo

    ! check for a new maximum depth 
    IF (parent%child(1)%nDepth > grid%maxDepth) THEN
       grid%maxDepth = parent%child(1)%nDepth
       grid%halfSmallestSubcell = grid%octreeRoot%subcellSize / &
            2.0_oc**REAL(grid%maxDepth,kind=oct)
       ! we store the value which is half the size of the 
       !   smallest subcell because this is more useful for later
       !   calculations.
    END IF
    
  end SUBROUTINE add_new_children


  !
  ! Recursively assigins some values after the grid is added. 
  !
  RECURSIVE SUBROUTINE finish_grid(thisOctal,grid, this, scaling_tau, &
       sigmaAbs, sigmaSca, meanDustParticleMass)
    ! takes the octree grid that has been created using 'splitGrid'
    !   and calculates all the other variables in the model.
    ! this should be called once the structure of the grid is complete.
    ! the gridConverged variable should be set .TRUE. when the entire
    !   grid has been finished. Until that happens, this subroutine 
    !   will be called repeatedly.
    
    
    IMPLICIT NONE

    TYPE(octal), POINTER   :: thisOctal
    TYPE(gridtype)         :: grid
    TYPE(alpha_disc), INTENT(IN)  :: this
    real, intent(in) :: scaling_tau
    real, intent(in) :: sigmaAbs, sigmaSca
    real, intent(in) :: meanDustParticleMass  ! [g]

    TYPE(octal), POINTER   :: pChild
  
    INTEGER :: subcell, n

    if (thisOctal%threeD) then
       n = 8
    else
       n = 4
    endif
     
    do subcell = 1, n
       if (thisOctal%hasChild(subcell)) then
          ! just decdend the tree branch
          pChild => thisOctal%child(subcell)
          CALL finish_grid(pChild,grid, this, scaling_tau, &
               sigmaAbs, sigmaSca, meanDustParticleMass)
       else
          ! assigin the values to the grid 
          call assign_values(thisOctal,subcell, grid, this, scaling_tau, &
               sigmaAbs, sigmaSca, meanDustParticleMass)
       end if
    end do

  END SUBROUTINE finish_grid
 
  
  !
  !
  !
  subroutine assign_values(thisOctal,subcell, grid, this, scaling_tau, &
       sigmaAbs, sigmaSca, meanDustParticleMass)
    IMPLICIT NONE
    
!    TYPE(octal), intent(inout) :: thisOctal
    TYPE(octal), pointer :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(alpha_disc), INTENT(IN)  :: this
    real, intent(in) :: scaling_tau
    real, intent(in) :: sigmaAbs, sigmaSca    ! [cm^2] photon-dust x-sections
    real, intent(in) :: meanDustParticleMass  ! [g]
    !
    TYPE(VECTOR)     :: cellCentre 
    real :: r
    real, parameter :: rho_bg = rho_min  ! background density
    ! the following should be really taken from a parameter files... 
    ! Need to fix this later. 
    real  :: M
    real(double) :: x, y, z, w, fac
    

    M = 1.989d33*this%Mcore ! [g] mass of the central object


    cellCentre = subcellCentre(thisOctal, subcell)
    x = cellCentre%x; y = cellCentre%y; z = cellCentre%z
    w = x*x + y*y
    r =SQRT(w+z*z)
    w = SQRT(w)

    if (w < this%Rh ) then
       ! then do not totch the grid.
       continue
!    elseif (thisOctal%rho(subcell) > rho_min) then
    elseif (in_alpha_disc(this,cellCentre) ) then 
       thisOctal%rho(subcell) = ave_alpha_disc_density(thisOctal, subcell, this)

       if (thisOctal%rho(subcell) > rho_min) thisOctal%inFlow(subcell) = .true.  
!       ! This is set to false to avoide resampling in the very high optical depth
!       ! in near the mid-plane of the disc
!       if (thisOctal%rho(subcell) > rho_min) thisOctal%inFlow(subcell) = .false.  

       ! assigning dust opacity
       ! This part should be modefied for oneKappa option later.
       ! The followings are devided by 100 since we assume the gas to dust ratio of 100.
       fac = thisOctal%rho(subcell)*1.0e10/100.e0/meanDustParticleMass
       thisOctal%kappaAbs(subcell,:) = sigmaAbs*fac
       thisOctal%kappaSca(subcell,:) = sigmaSca*fac
          
       ! In future, the temperature (and maybe also density) structure(s)
       ! should be calculated by Lucy's algorthem... But for now, we simply 
       ! set it to 1000 K.
       !
       thisOctal%temperature(subcell) = 1000.0
       thisOctal%biasCont3D(subcell) = 1.0e-30  ! should be no emission from disc
       thisOctal%biasLine3D(subcell) = 1.0e-30  ! should be no emission from disc
       thisOctal%etaLine(subcell) = 1.e-30
       thisOctal%etaCont(subcell) = 1.e-30
       thisOctal%chiLine(subcell) = 1.e-30
       thisOctal%velocity = keplerian_velocity(this,cellCentre)

    end if


  end subroutine assign_values
  
  



  real(double) function ave_alpha_disc_density(thisOctal,subcell, this)

    IMPLICIT NONE

    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    type(alpha_disc), intent(in) :: this
    
    real(oct)  :: cellSize, x, y, z
    TYPE(VECTOR)     :: cellCentre 
    real(double) :: c0, c1, c2, c3, c4, c5, c6, c7, c8, d, s

    ! get the size and the position of the centre of the current cell
    cellSize = thisOctal%subcellSize
    cellCentre = subcellCentre(thisOctal, subcell)
    x = cellcentre%x
    y = cellcentre%y
    z = cellcentre%z
    d = cellsize/3.00
    if (.not. thisOctal%threeD)  y = 0.0d0  ! 2D case

    ! assigning density 
    c0 =  alpha_disc_density(this, x, y, z, s)

    if (thisOctal%threeD)  then
       c1 =  alpha_disc_density(this, x+d, y+d, z+d, s)
       c2 =  alpha_disc_density(this, x+d, y+d, z-d, s)
       c3 =  alpha_disc_density(this, x+d, y-d, z+d, s)
       c4 =  alpha_disc_density(this, x+d, y-d, z-d, s)
       c5 =  alpha_disc_density(this, x-d, y+d, z+d, s)
       c6 =  alpha_disc_density(this, x-d, y+d, z-d, s)
       c7 =  alpha_disc_density(this, x-d, y-d, z+d, s)
       c8 =  alpha_disc_density(this, x-d, y-d, z-d, s)

       ave_alpha_disc_density =  &
            (c0 + c1+ c2 + c3 + c4 + c5 + c6 + c7 + c8)/9.0d0
    else
       c1 =  alpha_disc_density(this, x+d, y, z+d, s)
       c2 =  alpha_disc_density(this, x+d, y, z-d, s)
       c3 =  alpha_disc_density(this, x-d, y, z+d, s)
       c4 =  alpha_disc_density(this, x-d, y, z-d, s)

       ave_alpha_disc_density =  &
            (c0 + c1+ c2 + c3 + c4)/5.0d0
    end if

  end function ave_alpha_disc_density


    

  !
  ! Recursively turn off the inflow flag to false if in disc
  !
  RECURSIVE SUBROUTINE turn_off_disc(thisOctal,grid, this)    
    
    IMPLICIT NONE
    
    TYPE(octal), POINTER   :: thisOctal
    TYPE(gridtype)         :: grid
    TYPE(alpha_disc), INTENT(IN)  :: this
    
    TYPE(octal), POINTER   :: pChild
    type(VECTOR) :: rvec
    
    INTEGER :: subcell, n
    
    if (thisOctal%threeD) then
       n = 8
    else
       n =4
    end if
    !
    do subcell = 1, n
       if (thisOctal%hasChild(subcell)) then
          ! just decdend the tree branch
          pChild => thisOctal%child(subcell)
          CALL turn_off_disc(pChild,grid, this)
       else
          ! turnning it off 
!          if (modulus(thisOctal%centre) > this%Rh) thisOctal%inFlow(subcell) = .false.
          rvec = subcellCentre(thisOctal,subcell)
          if (in_alpha_disc(this, rvec)) then
             thisOctal%inFlow(subcell) = .false.
             thisOctal%kappaAbs(subcell,:) = 1.0e-30
             thisOctal%kappaSca(subcell,:) = 1.0e-30          
             thisOctal%temperature(subcell) = 1000.0
             thisOctal%biasCont3D(subcell) = 1.0e-30  ! should be no emission from disc
             thisOctal%biasLine3D(subcell) = 1.0e-30  ! should be no emission from disc
             thisOctal%etaLine(subcell) = 1.e-30
             thisOctal%etaCont(subcell) = 1.e-30
             thisOctal%cornerVelocity = vector(1.e-30,1.e-30,1.e-30)
             thisOctal%velocity = vector(1.e-30,1.e-30,1.e-30)
             thisOctal%rho(subcell) = rho_min
          end if
       end if
    end do
    
  END SUBROUTINE turn_off_disc


  !
  !
  ! Recursively turn off the inflow flag to false if grid > Rh.
  !
  RECURSIVE SUBROUTINE turn_on_disc(thisOctal,grid, this)  
    
    IMPLICIT NONE

    TYPE(octal), POINTER   :: thisOctal
    TYPE(gridtype)         :: grid
    TYPE(alpha_disc), INTENT(IN)  :: this

    TYPE(octal), POINTER   :: pChild
    type(VECTOR) :: rvec
    INTEGER :: subcell, n
    
    if (thisOctal%threeD) then
       n = 8
    else
       n =4
    end if
    !
    do subcell = 1, n
       if (thisOctal%hasChild(subcell)) then
          ! just decdend the tree branch
          pChild => thisOctal%child(subcell)
          CALL turn_on_disc(pChild,grid, this)
       else
          ! turnning it on
!          if (modulus(thisOctal%centre) > this%Rh) thisOctal%inFlow(subcell) = .true.
          rvec = subcellCentre(thisOctal,subcell)
          if (in_alpha_disc(this, rvec)) thisOctal%inFlow(subcell) = .true.
       end if
    end do
    
  END SUBROUTINE turn_on_disc



  !
  !
  !
  !
  SUBROUTINE fill_velocity_corners(this, thisOctal)
    ! store the velocity values at the subcell corners of an octal so
    !   that they can be used for interpolation.

    IMPLICIT NONE
    type(alpha_disc), intent(in) :: this
    TYPE(octal), INTENT(INOUT) :: thisOctal

    real(oct)      :: x1, x2, x3
    real(oct)      :: y1, y2, y3
    real(oct)      :: z1, z2, z3
    
    
    if (thisOctal%threeD) then

       ! we first store the values we use to assemble the position vectors        
       
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       
       y1 = thisOctal%centre%y - thisOctal%subcellSize
       y2 = thisOctal%centre%y
       y3 = thisOctal%centre%y + thisOctal%subcellSize
       
       z1 = thisOctal%centre%z - thisOctal%subcellSize
       z2 = thisOctal%centre%z
       z3 = thisOctal%centre%z + thisOctal%subcellSize
       
       ! now store the 'base level' values    
       thisOctal%cornerVelocity(1) = keplerian_velocity(this, VECTOR(x1,y1,z1))
       thisOctal%cornerVelocity(2) = keplerian_velocity(this, VECTOR(x2,y1,z1))
       thisOctal%cornerVelocity(3) = keplerian_velocity(this, VECTOR(x3,y1,z1))
       thisOctal%cornerVelocity(4) = keplerian_velocity(this, VECTOR(x1,y2,z1))
       thisOctal%cornerVelocity(5) = keplerian_velocity(this, VECTOR(x2,y2,z1))
       thisOctal%cornerVelocity(6) = keplerian_velocity(this, VECTOR(x3,y2,z1))
       thisOctal%cornerVelocity(7) = keplerian_velocity(this, VECTOR(x1,y3,z1))
       thisOctal%cornerVelocity(8) = keplerian_velocity(this, VECTOR(x2,y3,z1))
       thisOctal%cornerVelocity(9) = keplerian_velocity(this, VECTOR(x3,y3,z1))

       ! middle level
  
       thisOctal%cornerVelocity(10) = keplerian_velocity(this, VECTOR(x1,y1,z2))
       thisOctal%cornerVelocity(11) = keplerian_velocity(this, VECTOR(x2,y1,z2))
       thisOctal%cornerVelocity(12) = keplerian_velocity(this, VECTOR(x3,y1,z2))
       thisOctal%cornerVelocity(13) = keplerian_velocity(this, VECTOR(x1,y2,z2))
       thisOctal%cornerVelocity(14) = keplerian_velocity(this, VECTOR(x2,y2,z2))
       thisOctal%cornerVelocity(15) = keplerian_velocity(this, VECTOR(x3,y2,z2))
       thisOctal%cornerVelocity(16) = keplerian_velocity(this, VECTOR(x1,y3,z2))
       thisOctal%cornerVelocity(17) = keplerian_velocity(this, VECTOR(x2,y3,z2))
       thisOctal%cornerVelocity(18) = keplerian_velocity(this, VECTOR(x3,y3,z2))
       
       ! top level
       
       thisOctal%cornerVelocity(19) = keplerian_velocity(this, VECTOR(x1,y1,z3))
       thisOctal%cornerVelocity(20) = keplerian_velocity(this, VECTOR(x2,y1,z3))
       thisOctal%cornerVelocity(21) = keplerian_velocity(this, VECTOR(x3,y1,z3))
       thisOctal%cornerVelocity(22) = keplerian_velocity(this, VECTOR(x1,y2,z3))
       thisOctal%cornerVelocity(23) = keplerian_velocity(this, VECTOR(x2,y2,z3))
       thisOctal%cornerVelocity(24) = keplerian_velocity(this, VECTOR(x3,y2,z3))
       thisOctal%cornerVelocity(25) = keplerian_velocity(this, VECTOR(x1,y3,z3))
       thisOctal%cornerVelocity(26) = keplerian_velocity(this, VECTOR(x2,y3,z3))
       thisOctal%cornerVelocity(27) = keplerian_velocity(this, VECTOR(x3,y3,z3))

    else  ! 2D case

       ! we first store the values we use to assemble the position vectors               
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       
       y1 = 0.0d0
       y2 = 0.0d0
       y3 = 0.0d0
       
       z1 = thisOctal%centre%z - thisOctal%subcellSize
       z2 = thisOctal%centre%z
       z3 = thisOctal%centre%z + thisOctal%subcellSize
       
       ! now store the 'base level' values    
       thisOctal%cornerVelocity(1) = keplerian_velocity(this, VECTOR(x1,y1,z1))
       thisOctal%cornerVelocity(2) = keplerian_velocity(this, VECTOR(x2,y1,z1))
       thisOctal%cornerVelocity(3) = keplerian_velocity(this, VECTOR(x3,y1,z1))
       thisOctal%cornerVelocity(4) = keplerian_velocity(this, VECTOR(x1,y2,z2))
       thisOctal%cornerVelocity(5) = keplerian_velocity(this, VECTOR(x2,y2,z2))
       thisOctal%cornerVelocity(6) = keplerian_velocity(this, VECTOR(x3,y2,z2))
       thisOctal%cornerVelocity(7) = keplerian_velocity(this, VECTOR(x1,y3,z3))
       thisOctal%cornerVelocity(8) = keplerian_velocity(this, VECTOR(x2,y3,z3))
       thisOctal%cornerVelocity(9) = keplerian_velocity(this, VECTOR(x3,y3,z3))

    end if

    
  END SUBROUTINE fill_velocity_corners


  ! Given a position vector rvec(xpos, ypos and zpos) in 10^10 cm, this routine
  ! will return the velocity in [c] that point.
  ! Here we assume rvec to be the position vector originated from the origin (0,0,0) 
  ! of a cartisan coordinate system
  !   
  function keplerian_velocity(this, rvec) RESULT(out)
    implicit none
    type(vector) :: out
    !
    type(alpha_disc), intent(in) :: this
    type(VECTOR), intent(in) :: rvec
    !
    real(double) :: Vphi, w, r, x, y, z
    real(double) :: Vx, Vy, Vz, phi
    real(double), parameter :: c = 2.99792458d10   ! [cm/s] speed of light
    real(double), parameter :: G = 6.67259e-8 ! in cgs
    real(double) :: Mstar, RR
                                                                                                 
    
    x = rvec%x; y =rvec%y; z = rvec%z

    w = x*x + y*y
    r = SQRT(w + z*z)
    w = SQRT(w)

    if (r < this%Rh) then
       out = VECTOR(1.e-20, 1.e-20, 1.e-20)
    elseif (r > this%Rd) then
       out = VECTOR(1.e-20, 1.e-20, 1.e-20)
    else
       Mstar = this%Mcore * 1.9891d33  ! converting from [Msun] to [g]
       RR = w*1.d10                     ! converting from [10^10 cm] to [cm]
       Vphi = SQRT( G*Mstar / RR ) ! [cm/s]
       
       phi = ATAN2(y,x)
       
       Vx = -1.0d0*Vphi*SIN(phi)   ! [cm/s]
       Vy =        Vphi*COS(phi)   ! [cm/s]
       Vz = 0.0d0
       
       Vx = Vx/c  ! [c]
       Vy = Vy/c  ! [c]
       Vz = Vz/c  ! [c]
    
       out = VECTOR(Vx, Vy, Vz)
    end if
    
  end function keplerian_velocity




  !
  ! Function to check if a given point is in alpha disc
  ! 
  function in_alpha_disc(this, rvec) RESULT(out)
    implicit none
    logical :: out
    !
    type(alpha_disc), intent(in) :: this
    type(VECTOR), intent(in) :: rvec
    !
    real(double) ::  x, y, z !, rho, sizescale
!    real(double), parameter :: rho_min = 0.9d-19
    real(double) :: theta_max 
    real(double) :: cos_theta, r, w

!    theta_max = 70.0d0 ! degrees
!    theta_max = theta_max*3.4159d0/180.0d0  ! [radians] ! used by mistake
    theta_max = 80.0d0 ! degrees
    theta_max = theta_max*3.1459d0/180.0d0  ! [radians]

    x = rvec%x; y =rvec%y; z = rvec%z
    w = x*x+y*y
    r = SQRT(w + z*z)
    w =SQRT(w)
    cos_theta = z/r

    if (w > this%Rh .and. r < this%Rd  &
         .and. (ABS(cos_theta) < ABS(COS(theta_max))) ) then
       out =.true.
    else
       out = .false.
    end if

!    rho = alpha_disc_density(this, x, y, z, sizescale)

!    if (w > this%Rh .and. r < this%Rd  &
!         .and. (abs(z) < sizescale) ) then
!       out =.true.
!    else
!       out = .false.
!    end if

    
  end function in_alpha_disc




end module disc_class
