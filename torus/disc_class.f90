module disc_class
  
  use kind_mod
  use octal_mod
  use gridtype_mod
  use grid_mod
  use vector_mod
  use constants_mod

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
       &   turn_off_disc

  private::  &
       &   int_alpha_disc, &
       &   scale_size, &
       &   need_to_split, &
       &   add_new_children, &
       &   assign_values


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
  end type alpha_disc

  

  ! overload functions to make interface here ----------
  interface new
     module procedure int_alpha_disc_default
     module procedure int_alpha_disc
  end interface


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

  end subroutine int_alpha_disc_default
  
  
  ! inilializing with prarameters
  subroutine int_alpha_disc(this, Rh, Rd, Md, x0, y0, z0, sx, sy, sz) 
    implicit none 
    
    type(alpha_disc), intent(inout) :: this 
    ! this is for "alpha-disc" standard model 
    real(double), intent(in) :: Rh  ! [10^10 cm] inner radius of the disc (hole) 
    real(double), intent(in) :: Rd  ! [10^10 cm] outer radius of the disc 
    !
    real(double), intent(in) :: Md  ! mass of the disc in grams
    real(double), intent(in) :: x0, y0, z0  ! coordinates of the center [10^10cm]
    real(double), intent(in) :: sx,sy, sz  ! components of spin vector ! [units?]


    this%Rh = Rh  ! [10^10cm]
    this%Rd = Rd  ! [10^10cm]
    this%Md = Md  ! [g]
    this%x0 = x0;  this%y0 = y0;  this%z0=z0
    ! pointing z direction.  
    this%sx = sx;  this%sy = sy;  this%sz=sz

  end subroutine int_alpha_disc



  !==========================================================================
  ! Accessors
  !==========================================================================

  ! Given a alpha_disc object, it returns all components of the objects.
  subroutine get_parameters(this, Rh, Rd, Md, x0, y0, z0, sx, sy, sz) 
    implicit none 
    
    type(alpha_disc), intent(in) :: this 
    ! this is for "alpha-disc" standard model 
    real(double), intent(out) :: Rh  ! [10^10 cm] inner radius of the disc (hole) 
    real(double), intent(out) :: Rd  ! [10^10 cm] outer radius of the disc 
    !
    real(double), intent(out) :: Md  ! mass of the disc in grams
    real(double), intent(out) :: x0, y0, z0  ! coordinates of the center [10^10cm]
    real(double), intent(out) :: sx,sy, sz   ! components of spin vector ! [units?]
    
    Rh = this%Rh
    Rd = this%Rd
    Md = this%Md
    x0 = this%x0; y0 = this%y0; z0 = this%z0
    sx = this%sx; sy = this%sy; sz = this%sz
    
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
    real(double), parameter :: rho_min = 1.0d-25  ! g/cm^3 
    !
    !
    real(double), parameter :: rmax = 100.0d00  ! 10^10 cm
    !
    logical :: logscale
        
    !
    !
    ! holesize = 12.0d0 * 7.0d00 ! [10^10cm] 12 times the radius of a star
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
    radius2 = xposi**2 + yposi**2 + zposi**2

    radius = SQRT(radius2)
    if (radius < discrad) then          
       if (radius < 1.1d0*holesize) then
          sizescale = holesize
          if (radius > (holesize/2.0d0)) sizescale=sizescale/2.0d0
       else
          rad = radius/discrad
          !  Distance above disc plane is obtained by doing a dot product with the
          !     unit vector giving the disc orientation
          unitx = spinx
          unity = spiny
          unitz = spinz
          unorm = SQRT(unitx**2 + unity**2 + unitz**2)
          zdist = (xposi*unitx + yposi*unity + zposi*unitz)/unorm
          !  Set density and characteristic size scale at position
          height = 0.05d0*discrad*rad**(1.125d0)
          height = height/2.0d0

          densitymidplane = 5.0d0*discmass/(8.0d0*pi*discrad**2)*  &
               rad**(-0.75d0)

          density_out = densitymidplane*   &
               EXP(-zdist**2/(2.0d0*height**2))/(SQRT(2.0d0*pi)*height)

          
          ! now finding the scale size at (xpos, ypos, zpos)

          if (density_out>rho_min) then
                          
             !
             ! Produces the > 500 Mb grid
             !
             
             ! Works with  height = height/2.0d0
             logscale = .false.
             sizescale = scale_size(abs(zdist), 0.0d0,  20.0d0*height, &
                  5.0d0*height,  2.0d0*height, logscale)


                
             ! limit the smallest and largest sizes
             sizescale = MAX(sizescale, 1.0d-1)
             sizescale = MIN(sizescale, 2.0d3)                

             ! else takes the default density from the top.
          end if
       endif
       !  given position can only be inside 1 disc at most (discs cannot overlap) 
       goto 10
    endif
10  continue    

    ! converting units from g/(10^10cm)^3 to g/cm^3
    density_out = density_out/1.0d30   ! [g/cm^3]

    density_out = density_out*10.0   ! [g/cm^3]

!    ! assuming dust density is 100 times smaller than the gas density
!    density_out = density_out/1.0d2   ! [g/cm^3]


    ! cutoff density
    ! max and min  density allowed
    if (density_out>rho_max) density_out = rho_max ! [g/cm^3]
    if (density_out<rho_min) density_out = rho_min ! [g/cm^3]

    if (radius > rmax) density_out = rho_min ! [g/cm^3]



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
       p = LOG(ds_max/ds_min) / (s_max - s_min)
       out = ds_min * EXP(p*(ss-s_min))
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
    TYPE(octalVector)     :: cellCentre 

    splitThis = .false.
    subcell = 1

    if (thisOctal%hasChild(subcell)) then
       ! find the index of the new child and decend the tree
       do i = 1, 8
          childPointer => thisOctal%child(i)
          CALL add_alpha_disc(childPointer,grid,this)
       end do
    else   

       ! get the size and the position of the centre of the current cell
       cellCentre = subcellCentre(thisOctal, subcell)

       if (modulus(cellCentre) > this%Rh) then 
          do while ((.not.splitThis).and.(subcell <= 8))
             IF (need_to_split(thisOctal,subcell,this)) splitThis = .true.           
             subcell = subcell + 1
          end do
          
          if (splitThis) then
             ! using a routine in this module
             CALL add_new_children(thisOctal, grid, this)
             ! find the index of the new child and decend the tree
             do i = 1, 8
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
    
    real(oct)  :: cellSize, x, y, z
    real(double) :: rho_disc, scale_length
    TYPE(octalVector)     :: cellCentre 

    ! get the size and the position of the centre of the current cell
    cellSize = thisOctal%subcellSize*2.0d0
    cellCentre = subcellCentre(thisOctal, subcell)
    x = cellcentre%x
    y = cellcentre%y
    z = cellcentre%z
    

    rho_disc = alpha_disc_density(this, x, y, z, scale_length)
  
    if (cellsize > scale_length) then       
       need_to_split = .true.  ! units in 10^10cm
    else
       need_to_split = .false.
    end if

               

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

    real(double) :: rho
    TYPE(OCTAL), POINTER :: childPointer      
    
    if (any(parent%hasChild(1:8))) write(*,*) "parent already has a child"

    parent%nChildren = 8

    ALLOCATE(parent%child(1:8), STAT=error)
    IF ( error /= 0 ) THEN
       PRINT *, 'Panic: allocation failed.'
       STOP
    END IF

    ! update the parent octal
    
    parent%hasChild(1:8) = .TRUE.
    parent%indexChild(1) = 1
    parent%indexChild(2) = 2
    parent%indexChild(3) = 3
    parent%indexChild(4) = 4
    parent%indexChild(5) = 5
    parent%indexChild(6) = 6
    parent%indexChild(7) = 7
    parent%indexChild(8) = 8

    do newChildIndex = 1, 8

       ! allocate any variables that need to be  
       ALLOCATE(parent%child(newChildIndex)%kappaAbs(8,grid%nOpacity))
       ALLOCATE(parent%child(newChildIndex)%kappaSca(8,grid%nOpacity))
       parent%child(newChildIndex)%kappaAbs(:,:) = 1.0e-30
       parent%child(newChildIndex)%kappaSca(:,:) = 1.0e-30

       ALLOCATE(parent%child(newChildIndex)%N(8,grid%maxLevels))
       NULLIFY(parent%child(newChildIndex)%child)

       ! set up the new child's variables
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
       parent%child(newChildIndex)%inFlow = .false.
       parent%child(newChildIndex)%temperature = 10.0

       
       ! put some data in the eight subcells of the new child
       DO subcell = 1, 8
          childPointer => parent%child(newChildIndex)
          rho =ave_alpha_disc_density(childPointer, subcell, this)
          parent%child(newChildIndex)%rho(subcell) = rho          
          parent%child(newChildIndex)%label(subcell) = counter
          counter = counter + 1

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
  RECURSIVE SUBROUTINE finish_grid(thisOctal,grid, this, scaling_tau)
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

    TYPE(octal), POINTER   :: pChild
  
    INTEGER :: subcell
     
    do subcell = 1, 8
       if (thisOctal%hasChild(subcell)) then
          ! just decdend the tree branch
          pChild => thisOctal%child(subcell)
          CALL finish_grid(pChild,grid, this, scaling_tau)
       else
          ! assigin the values to the grid 
          call assign_values(thisOctal,subcell, grid, this, scaling_tau)
       end if
    end do

  END SUBROUTINE finish_grid
 
  
  !
  !
  !
  subroutine assign_values(thisOctal,subcell, grid, this, scaling_tau)
    IMPLICIT NONE
    
    TYPE(octal), intent(inout) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(alpha_disc), INTENT(IN)  :: this
    real, intent(in) :: scaling_tau
    !
    TYPE(octalVector)     :: cellCentre 
    real :: r
    real, parameter :: rho_min = 1.e-25
    real, parameter :: rho_bg = 1.e-25
    TYPE(Vector)    ::  rvec, vvec
    real :: vel 
    ! the following should be really taken from a parameter files... 
    ! Need to fix this later. 
    TYPE(Vector), parameter    ::  spinAxis = Vector(0.0, 0.0, 1.0)
    real, parameter  :: Mcore= 1.989e33 ! 1 solar mass in gm


    cellCentre = subcellCentre(thisOctal, subcell)

    r = modulus(cellCentre)

    if (r < this%Rh ) then
       ! then do not totch the grid.
       continue
    else       
       if (thisOctal%rho(subcell) > rho_min) then
          thisOctal%inFlow(subcell) = .true.  

          thisOctal%kappaAbs(subcell,:) = 1.e+30/scaling_tau
          thisOctal%kappaSca(subcell,:) = 1.e+30/scaling_tau
          
          thisOctal%temperature(subcell) = 10.0
          thisOctal%velocity = VECTOR(0.,0.,0.)
          thisOctal%biasCont3D(subcell) = 1.0
          thisOctal%biasLine3D(subcell) = 1.0
          thisOctal%etaLine(subcell) = 1.e-30
          thisOctal%etaCont(subcell) = 1.e-30

          ! Assigining the rotational velocity 
          rVec = cellCentre
          rVec = rVec / r  ! r_hat
          vVec = rVec  .cross. spinAxis
          if (modulus(vVec) /= 0.) then
             call normalize(vVec)
             vel = sqrt(bigG * Mcore / (r*1.0e10) ) / cSpeed
             vVec = vel  * vVec
             thisOctal%velocity = vVec
          else
             thisOctal%velocity = VECTOR(0.,0.,0.)
          end if
       else
          
          thisOctal%rho(subcell) =  rho_bg
          
          thisOctal%inFlow(subcell) = .false.
          
          thisOctal%temperature(subcell) = 10.0
          thisOctal%velocity = VECTOR(0.,0.,0.)
          thisOctal%biasCont3D(subcell) = 1.0
          thisOctal%biasLine3D(subcell) = 1.0
          thisOctal%etaLine(subcell) = 1.e-30
          thisOctal%etaCont(subcell) = 1.e-30
          
       end if

    end if


  end subroutine assign_values
  
  



  real(double) function ave_alpha_disc_density(thisOctal,subcell, this)

    IMPLICIT NONE

    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    type(alpha_disc), intent(in) :: this
    
    real(oct)  :: cellSize, x, y, z
    TYPE(octalVector)     :: cellCentre 
    real(double) :: c0, c1, c2, c3, c4, c5, c6, c7, c8, d, s

    ! get the size and the position of the centre of the current cell
    cellSize = thisOctal%subcellSize
    cellCentre = subcellCentre(thisOctal, subcell)
    x = cellcentre%x
    y = cellcentre%y
    z = cellcentre%z
    d = cellsize/3.00

    ! assigning density 
    c0 =  alpha_disc_density(this, x, y, z, s)
    
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

  end function ave_alpha_disc_density


    

  !
  ! Recursively turn off the inflow flag to false if grid > Rh.
  !
  RECURSIVE SUBROUTINE turn_off_disc(thisOctal,grid, this)    
    
    IMPLICIT NONE
    
    TYPE(octal), POINTER   :: thisOctal
    TYPE(gridtype)         :: grid
    TYPE(alpha_disc), INTENT(IN)  :: this
    
    TYPE(octal), POINTER   :: pChild
    
    INTEGER :: subcell
    
    do subcell = 1, 8
       if (thisOctal%hasChild(subcell)) then
          ! just decdend the tree branch
          pChild => thisOctal%child(subcell)
          CALL turn_off_disc(pChild,grid, this)
       else
          ! turnning it off 
          if (modulus(thisOctal%centre) > this%Rh) thisOctal%inFlow(subcell) = .false.
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
  
    INTEGER :: subcell
     
    do subcell = 1, 8
       if (thisOctal%hasChild(subcell)) then
          ! just decdend the tree branch
          pChild => thisOctal%child(subcell)
          CALL turn_on_disc(pChild,grid, this)
       else
          ! turnning it on
          if (modulus(thisOctal%centre) > this%Rh) thisOctal%inFlow(subcell) = .true.        
       end if
    end do
    
  END SUBROUTINE turn_on_disc






end module disc_class
