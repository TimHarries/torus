module jet_class
  
  use kind_mod
  use octal_mod
  use gridtype_mod
  use grid_mod
  use vector_mod

  !  
  ! Class definition for a simple kinematic model of jets from TTauri star
  !
  ! Created: 31-jan-2005  (Ryuich Kurosawa)
 
  public::  &
       &   new, &
       &   get_parameters, &
       &   component, &
       &   jet_density, &
       &   jet_velocity, &
       &   ave_jet_density, &
       &   add_jet, &
       &   turn_on_jet, &
       &   turn_off_jet, &
       &   finish_grid_jet

  private::  &
       &   int_jet, &
       &   scale_size, &
       &   need_to_split, &
       &   need_to_split2, &
       &   add_new_children, &
       &   V_phi, &
       &   V_r, &
       &   in_jet_flow, &
       &   fill_velocity_corners, &
       &   assign_values



  !===============================================================================
  ! This model assumes the following form of the velocity and the density.
  !
  !      ->           ^             ^
  !      V(r) = Vr(r) r  + Vphi(R) phi   -- (1)
  ! 
  !  where Vr(r) = V0 + Vinf (1 - Rin/w)**beta  
  !   and  Vphi(R) = gamma*SQRT(GM_*/w).
  ! Note that r = SQRT(x^2 +y^2 + z^2) and w=SQRT(x^2+z^2)  -- (2)
  !
  ! The density is assume to be in the following form:
  !
  !       rho(r, theta) = rho_r(r)*rho_theta(theta)  --(3)
  ! 
  !Further, we assume the following form of rho_theta
  !
  !       rho_theta(theta) = a * cos(theta)**b  --(4)
  !
  ! Using (1) -- (4), Mass-loss rate in jets, and the mass flux conservation
  ! we find the form of the rho_r(r)
  !
  !
  !                               Mdot ( 1 + b )
  !    rho_r(r) = -----------------------------------------------------
  !                 4Pi r^2 Vr(r) a { 1   - cos(theta_j)**(1+b)  }
  !                                     
  !                                                                          --(5)
  !
  ! where theta_j is the jet openning angle measured from the z-axis.
  !===============================================================================

  type jet
     ! DO NOT MAKE THE COMPONETS PUBLIC! PLEASE ACCESS THE PARAMTERS THROUGH
     ! GET_PARAMETERS OR COMPONENT FUNCTION.
     private  
     real(double) :: Rin     !  [10^10 cm] inner radius of the disc (hole) 
     real(double) :: Rout    !  [10^10 cm] outer radius of jets
     real(double) :: theta_j !  [radian]  jet opening angle
     !
     real(double) :: Mstar   ! [Msun] Mass of the star
     real(double) :: Mdot    ! [Msun/yr] mass loss rate in the jets
     real(double) :: a       ! [-] a parameter in density function 
     real(double) :: b       ! [-] a parameter in density function
     real(double) :: V0      ! [km/s] Base velocity of jets
     real(double) :: Vinf    ! [km/s] Terminal velocity of jets
     real(double) :: beta    ! [-] a parameter in velocity function
     real(double) :: gamma   ! [-] a parameter in velocity function
     !
     real(double) :: limitscalar   ! [-] used for splitting cells
     real(double) :: T       ! [K]  Isothermal temperature of jets

  end type jet

  

  ! overload functions to make interface here ----------
  interface new
     module procedure int_jet_default
     module procedure int_jet
  end interface


  ! the minimum density used in this module
  real(double), private, parameter :: rho_min=1.0d-19  ! background density
  


contains
  ! =====================================================================
  ! constructors
  !======================================================================

  ! initializing the jet with default parameters
  subroutine int_jet_default(this)
    implicit none
    type(jet), intent(inout) :: this 

    this%Rin  = 6.0d0  *6.95d0              ! 6 R_sun [10^10cm]
    this%Rout = 10.0d0 *1.5d3               ! 10 AU  [10^10cm]
    this%theta_j = 2.0d0*ACOS(0.0d0)/3.0d0  ! pi/3 
    this%Mstar = 1.0d0                      ! [Msun]
    this%Mdot = 1.0d-9                      ! [Msun/yr]
    this%a = 0.8d0                          ! [-]
    this%b =2.0d0                           ! [-]
    this%V0 = 20.d0                         ! [km/s]
    this%Vinf = 200.d0                      ! [km/s]
    this%beta=0.5                           ! [-]
    this%gamma=5.0d-2                       ! [-]
    this%limitscalar = 1.e20                ! usually in [g]
    this%T = 1.0d4                          ! [K]
  end subroutine int_jet_default
  
  
  ! inilializing with prarameters
  subroutine int_jet(this, Rin, Rout,theta_j, Mstar, Mdot, a, b,  &
       V0, Vinf, beta, gamma, limitscalar, T)
    implicit none     
    type(jet), intent(inout) :: this 
    !
    real(double), intent(in) :: Rin     !  [10^10 cm] inner radius of the disc (hole) 
    real(double), intent(in) :: Rout    !  [10^10 cm] outer radius of jets
    real(double), intent(in) :: theta_j !  [radian]  jet opening angle
    !
    real(double), intent(in) :: Mstar   ! [Msun] Mass of the star
    real(double), intent(in) :: Mdot    ! [Msun/yr] mass loss rate in the jets
    real(double), intent(in) :: a       ! [-] a parameter in density function 
    real(double), intent(in) :: b       ! [-] a parameter in density function
    real(double), intent(in) :: V0      ! [km/s] Base velocity of jets
    real(double), intent(in) :: Vinf    ! [km/s] Terminal velocity of jets
    real(double), intent(in) :: beta    ! [-] a parameter in velocity function
    real(double), intent(in) :: gamma   ! [-] a parameter in velocity function
    real(double), intent(in) :: limitscalar   ! [-] used for splitting cells
    real(double), intent(in) :: T       ! [K]  Isothermal temperature of jets
    !
    this%Rin  = Rin
    this%Rout = Rout
    this%theta_j = theta_j
    this%Mstar = Mstar
    this%Mdot = Mdot
    this%a = a
    this%b = b
    this%V0 = V0
    this%Vinf = Vinf
    this%beta= beta
    this%gamma = gamma
    this%limitscalar = limitscalar
    this%T = T

  end subroutine int_jet



  !==========================================================================
  ! Accessors
  !==========================================================================

  ! Given a jet object, it returns all components of the objects.
  subroutine get_parameters(this, Rin, Rout,theta_j, Mstar, Mdot, a, b, &
       V0, Vinf, beta, gamma, limitscalar, T)
    implicit none 
    
    type(jet), intent(in) :: this 
    !
    real(double), intent(out) :: Rin     !  [10^10 cm] inner radius of the disc (hole) 
    real(double), intent(out) :: Rout    !  [10^10 cm] outer radius of jets
    real(double), intent(out) :: theta_j !  [radian]  jet opening angle
    !
    real(double), intent(out) :: Mstar   ! [Msun] Mass of the star
    real(double), intent(out) :: Mdot    ! [Msun/yr] mass loss rate in the jets
    real(double), intent(out) :: a       ! [-] a parameter in density function 
    real(double), intent(out) :: b       ! [-] a parameter in density function
    real(double), intent(out) :: V0      ! [km/s] Base velocity of jets
    real(double), intent(out) :: Vinf    ! [km/s] Terminal velocity of jets
    real(double), intent(out) :: beta    ! [-] a parameter in velocity function
    real(double), intent(out) :: gamma   ! [-] a parameter in velocity function
    real(double), intent(out) :: limitscalar   ! [-] used for splitting cell
    real(double), intent(out) :: T       ! [K]  Isothermal temperature of jets

    Rin = this%Rin
    Rout = this%Rout
    theta_j = this%theta_j
    Mstar = this%Mstar
    Mdot = this%Mdot
    a = this%a
    b = this%b
    V0 = this%V0
    beta = this%beta
    gamma = this%gamma
    limitscalar = this%limitscalar
    T = this%T
        
  end subroutine get_parameters


  ! Given a jet object and the name of its component, this function 
  ! returns the component
  real(double) function component(this, name) 
    implicit none 
    
    type(jet), intent(in) :: this 
    character(LEN=*), intent(in) :: name 
    
    select case (name)
    case ("Rin")
       component = this%Rin
    case ("Rout")
       component = this%Rout
    case ("theta_j")
       component = this%theta_j
    case ("Mstar")
       component = this%Mstar
    case ("Mdot")
       component = this%Mdot
    case ("a")
       component = this%a
    case ("b")
       component = this%b
    case ("V0")
       component = this%V0
    case ("Vinf")
       component = this%Vinf
    case ("beta")
       component = this%beta
    case ("gamma")
       component = this%gamma
    case ("limitscalar")
       component = this%limitscalar
    case ("T")
       component = this%T
    case default
       write(*,*) "Error:: Unknown NAME was passed to [jet_class::component]."
       write(*,*) "Exiting the program because of this error!"
       stop
    end select

  end function component


  !=============================================================================
  ! Tools
  !=============================================================================


  ! Given a position vector rvec(xpos, ypos and zpos) in 10^10 cm, this routine
  ! will return the density [g/cm^3].
  ! Here we assume rvec to be the position vector originated from the origin (0,0,0) 
  ! of a cartisan coordinate system
  !   
  function jet_density(this, rvec) RESULT(density_out)
    implicit none
    real(double) :: density_out  ! in [g/cm^3]
    type(jet), intent(in) :: this
    type(octalvector), intent(in) :: rvec  ! componets should be in [10^10cm]
    !
    real(double), parameter :: pi = 3.14159265d0
    !real(double), parameter :: rho_min = 1.0d-22
    real(double), parameter :: rho_max = 1.0d-5
    !
    real(double) :: r, w , x, y, z, r_cm, cos_theta, cos_theta_j
    real(double) :: theta, rho_theta, rho_r, Vr, f, c, Mdot
    
    x = rvec%x; y = rvec%y; z=rvec%z

    w = x*x + y*y
    r = SQRT(w + z*z)  ! radial distance
    w = SQRT(w)        ! distance from z axis

    cos_theta = z/r
    cos_theta_j = COS(this%theta_j)
    theta = ACOS(cos_theta)

    if (.not. in_jet_flow(this, rvec) ) then
       density_out = rho_min ! [g/cm^3]
    else

       Vr = V_r(this, r)  ! [km/s] ...  Using a function in this module  
       Vr = Vr*1.0d5      ! [cm/s] ..   uni conversion

!       rho_theta = (1.0d0 - this%a * cos_theta**this%b)   ! [-]
       rho_theta = this%a * cos_theta**this%b   ! [-]

       c = (1.0d0+this%b)
       f = (1.0d0 - cos_theta_j**c) / c

       Mdot = this%Mdot*1.989d33/3.1536d07   ! Converting from [Msun/yr] to [g/s] 
       r_cm = r*1.0d10                      ! converting from [10^10 cm] to [cm]

       rho_r =  Mdot/ ( 4.0d0*pi*r_cm*r_cm*Vr*f )  ! should be in [g/cm^3]

       density_out = rho_r*rho_theta
    end if


    ! Limits the range of the density for safty
    if (density_out>rho_max) density_out = rho_max ! [g/cm^3]
    if (density_out<rho_min) density_out = rho_min ! [g/cm^3]

  end function jet_density




  !
  ! Given the radial distance (r) from the origin and this object, this function 
  ! will return the polar component of the velocity field in [km/s]
  !
  function V_r(this, r)  RESULT(out)
    implicit none
    real(double) :: out
    type(jet), intent(in) :: this
    real(double), intent(in) :: r  ! [10^10cm] 
    
    out = this%V0 + this%Vinf*(1.0d0-this%Rin/r)**this%beta   ! [km/s]
    
  end function V_r



  !
  ! Given the distance from the z-axis (w) and this object, this function 
  ! will return the polar component of the velocity field in [km/s]
  !
  function V_phi(this, w)  RESULT(out)
    implicit none
    real(double) :: out
    type(jet), intent(in) :: this
    real(double), intent(in) :: w   ! [10^10cm]
    real(double), parameter :: G = 6.67259d-8 ! in cgs
    real(double) :: Mstar, R

    Mstar = this%Mstar * 1.9891d33  ! converting from [Msun] to [g]
    R = w*1.0d10                    ! converting from [10^10 cm] to [cm]

    if (w > this%Rin) then
       out = this%gamma * SQRT( G*Mstar / R )  ! [cm/s]
       out = out*1.0d-5   ! [ km/s]
    else ! gives no Vphi
       out = 1.0d-30
    end if
  end function V_phi
 



  ! Given a position vector rvec(xpos, ypos and zpos) in 10^10 cm, this routine
  ! will return the velocity in [c] that point.
  ! Here we assume rvec to be the position vector originated from the origin (0,0,0) 
  ! of a cartisan coordinate system
  !   
  function jet_velocity(this, rvec) RESULT(out)
    implicit none
    type(vector) :: out
    !
    type(jet), intent(in) :: this
    type(octalvector), intent(in) :: rvec
    !
    real(double) :: Vr, Vphi, w, r, x, y, z
    real(double) :: Vx, Vy, Vz, phi, theta
    real(double), parameter :: c = 2.99792458d5   ! [km/s] speed of light
    
    x = rvec%x; y =rvec%y; z = rvec%z

    w = x*x + y*y
    r = SQRT(w + z*z)
    w = SQRT(w)

    if (.not. in_jet_flow(this, rvec)) then
       out = VECTOR(1.e-20, 1.e-20, 1.e-20)
    else
       ! using the functions in this module
       Vphi = V_phi(this, w)  !  [km/s]
       Vr = V_r(this,r)       !  [km/s]
       
       phi = ATAN2(y,x)
       theta = ACOS(z/r)
       
       Vx = Vr*SIN(theta)*COS(phi) - Vphi*SIN(phi)        ! [km/s]
       Vy = Vr*SIN(theta)*SIN(phi) + Vphi*COS(phi)        ! [km/s]
       Vz = Vr*COS(theta)                                 ! [km/s]

       Vx = Vx/c  ! [c]
       Vy = Vy/c  ! [c]
       Vz = Vz/c  ! [c]
    
       out = VECTOR(Vx, Vy, Vz)
    end if
    
  end function jet_velocity

  !
  !  Given a grid object, this routine will add extra grid and assign density (and etc) 
  !  to the grid using the density function defined in this module.  !  

  RECURSIVE SUBROUTINE add_jet(thisOctal,grid,this)
    ! uses an external function to decide whether to split a subcell of
    !   the current octal. 

    IMPLICIT NONE

    TYPE(OCTAL), POINTER :: thisOctal
      ! 'limitScalar' is the value the decideSplit function uses to
      !   decide whether or not to split cell.
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid through to the 
    type(jet), intent(in)  :: this !
    
    !
    TYPE(OCTAL), POINTER :: childPointer  
    INTEGER              :: subcell, i    ! loop counters
    logical :: splitThis
    TYPE(octalVector)     :: cellCentre 
    integer :: n

    splitThis = .false.
    subcell = 1

    if (thisOctal%threeD) then
       n = 8
    else
       n = 4
    endif


    if (thisOctal%hasChild(subcell)) then
       ! find the index of the new child and decend the tree
       do i = 1, n
          childPointer => thisOctal%child(i)
          CALL add_jet(childPointer,grid,this)
       end do
    else   

       ! get the size and the position of the centre of the current cell
       cellCentre = subcellCentre(thisOctal, subcell)

       if (modulus(cellCentre) > this%Rin .or. modulus(cellCentre) < this%Rout) then 
          do while ((.not.splitThis).and.(subcell <= n))
             IF (need_to_split2(thisOctal,subcell,this)) splitThis = .true.           
             subcell = subcell + 1
          end do
          
          if (splitThis) then
             ! using a routine in this module
             CALL add_new_children(thisOctal, grid, this)
             ! find the index of the new child and decend the tree
             do i = 1, n
                childPointer => thisOctal%child(i)
                CALL add_jet(childPointer,grid,this)
             end do
          end if
       end if
    end IF
    
   END SUBROUTINE add_jet


  ! This is based decideSplit in amr_mod.f90

  logical function need_to_split(thisOctal,subcell, this)

    IMPLICIT NONE

    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    type(jet), intent(in) :: this
    
    real(oct)  :: cellSize
    real(double) :: rho_disc, mass_cell
    TYPE(octalVector)     :: cellCentre

    ! get the size and the position of the centre of the current cell
    cellSize = thisOctal%subcellSize*2.0d0
    cellCentre = subcellCentre(thisOctal, subcell)


    if (.not.in_jet_flow(this, cellCentre) ) then
       need_to_split = .false.
    else
       rho_disc = ave_jet_density(thisOctal, subcell, this)
       if (thisOctal%maxChildren == 8)  then ! 3D case
          mass_cell = rho_disc * (cellSize*1.0d10)**3  ! should be in grams
       else  ! 2-D AMR
          mass_cell = rho_disc * pi *  &
               ( (cellCentre%x+cellsize/2.)**2 - (cellCentre%x-cellsize/2.)**2) &
               * cellsize * 1.d30 ! should be in grams
       end if

       !
       if (mass_cell > this%limitscalar) then       
          need_to_split = .true. 
       else
          need_to_split = .false.
       end if
    end if

  end function need_to_split

  logical function need_to_split2(thisOctal,subcell, this)

    IMPLICIT NONE

    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    type(jet), intent(in) :: this
    
    real(oct)  :: cellSize, d
    real(double) :: rho_disc, mass_cell
    TYPE(octalVector)     :: cellCentre
    integer, parameter :: nr = 100
    real(double) :: r
    integer :: i
    !
    logical, save :: first_time = .true.
    real(double) , save:: rGrid(nr)
    real(double) :: Rmax = 1.5d3  ! [10^10cm]

    need_to_split2 = .false.

    if (first_time) then  ! setup ref. r grid
       do i = 1, nr
          rGrid(i) = log10(this%rIn)+dble(i-1)/dble(nr-1)*(log10(this%rOut)-log10(this%rIn))
!          rGrid(i) = log10(this%rIn)+dble(i-1)/dble(nr-1)*(log10(Rmax)-log10(this%rIn))
       enddo
       do i = 1, nr
          rGrid(i) = 10.d0**rGrid(i)
       end do
       first_time = .false.
    end if

    cellSize = (thisOctal%subcellSize)*2.0d0
    cellCentre = subcellCentre(thisOctal, subcell)


    if (.not.in_jet_flow(this, cellCentre) ) then
       need_to_split2 = .false.
    else
       ! get the size and the position of the centre of the current cell
       r = modulus(cellCentre)
       call locate(rGrid,nr,r,i)
       if (i > (nr-1)) i = nr-1
       d = rGrid(i+1) - rGrid(i)
       if (cellSize > d ) then
          need_to_split2 = .true.
       end if
    endif
  end function need_to_split2


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
    TYPE(jet), INTENT(IN)  :: this
    INTEGER       :: subcell           ! loop counter
    INTEGER, SAVE :: counter = 9       ! keeps track of the current subcell label
                                       ! - this isn't very clever. might change it. 
    INTEGER       :: newChildIndex     ! the storage location for the new child
    INTEGER       :: error

    real(double) :: rho
    TYPE(OCTAL), POINTER :: childPointer      

    if (any(parent%hasChild(1:8))) &
         write(*,*) "parent already has a child  [jet_class::add_new_children]"

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
       PRINT *, 'Error: Allocation failed in [jet_class::add_new_children].'
       STOP
    END IF



    ! update the parent octal
    
    parent%hasChild(1:8) = .TRUE.
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
       parent%child(newChildIndex)%chiLine = 1.e-30
       parent%child(newChildIndex)%etaLine = 1.e-30
       parent%child(newChildIndex)%etaCont = 1.e-30
       parent%child(newChildIndex)%N = 1.e-10
       parent%child(newChildIndex)%Ne = 1.e-10
       parent%child(newChildIndex)%changed = .false.
       ! 
       parent%child(newChildIndex)%temperature(:) = 3000.0
       parent%child(newChildIndex)%biasCont3D(:) = 1.0
       parent%child(newChildIndex)%biasLine3D(:) = 1.0
       parent%child(newChildIndex)%etaLine(:) = 1.e-30
       parent%child(newChildIndex)%etaCont(:) = 1.e-30
       parent%child(newChildIndex)%inFlow(:) =  .false.
       parent%child(newChildIndex)%rho(:) = rho_min   ! [g/cm^3]
       !


       ! put some data in the eight subcells of the new child
       DO subcell = 1, parent%maxChildren
          childPointer => parent%child(newChildIndex)
          ! assigiing velocity
          parent%child(newChildIndex)%velocity(subcell)  &
               = jet_velocity(this, childPointer%Centre)   ! [c]
          if (subcell == parent%maxChildren) &
               call fill_velocity_corners(this, childPointer)          
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
  !
  !
  subroutine assign_values(thisOctal,subcell, this)
    IMPLICIT NONE
    
    TYPE(octal), pointer :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(jet), INTENT(IN)  :: this
    !
    TYPE(octalVector)     :: cellCentre 

    cellCentre = subcellCentre(thisOctal, subcell)
    if (in_jet_flow(this, cellCentre) ) then
       thisOctal%rho(subcell) = ave_jet_density(thisOctal, subcell, this)
       if (thisOctal%rho(subcell) > rho_min) thisOctal%inFlow(subcell) = .true.         
       thisOctal%temperature(subcell)=this%T
       
       thisOctal%velocity = jet_velocity(this,cellCentre)
       if (subcell == thisOctal%maxChildren) &
            call fill_velocity_corners(this, thisOctal)

       thisOctal%biasCont3D(subcell) = 1.0
       thisOctal%biasLine3D(subcell) = 1.0
       thisOctal%etaLine(subcell) = 1.e-30
       thisOctal%etaCont(subcell) = 1.e-30       
    else
       ! don't touch the grid
       continue
    end if


  end subroutine assign_values




  !
  ! Recursively assigins some values after the grid is added. 
  !
  RECURSIVE SUBROUTINE finish_grid_jet(thisOctal, this)
    
    IMPLICIT NONE

    TYPE(octal), POINTER   :: thisOctal
    TYPE(jet), INTENT(IN)  :: this

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
          CALL finish_grid_jet(pChild, this)
       else
          ! assigin the values to the grid 
          call assign_values(thisOctal, subcell, this)
       end if
    end do
    
  END SUBROUTINE finish_grid_jet




  !
  !
  !
  real(double) function ave_jet_density(thisOctal,subcell, this)

    IMPLICIT NONE

    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    type(jet), intent(in) :: this
    
    real(oct)  :: cellSize, x, y, z
    TYPE(octalVector)     :: cellCentre, rvec
    real(double) :: c0, c1, c2, c3, c4, c5, c6, c7, c8, d, s

    ! get the size and the position of the centre of the current cell
    cellSize = thisOctal%subcellSize
    cellCentre = subcellCentre(thisOctal, subcell)
    x = cellcentre%x
    y = cellcentre%y
    z = cellcentre%z
    d = cellsize/3.00
    if (.not.thisOctal%threeD)  y = 0.0d0  ! 2D case

    if (thisOctal%threeD)  then
       ! assigning density 
       rvec = OCTALVECTOR(x, y, z);        c0 =  jet_density(this, rvec)    
       rvec = OCTALVECTOR(x+d, y+d, z+d);  c1 =  jet_density(this, rvec)
       rvec = OCTALVECTOR(x+d, y+d, z-d);  c2 =  jet_density(this, rvec)
       rvec = OCTALVECTOR(x+d, y-d, z+d);  c3 =  jet_density(this, rvec)
       rvec = OCTALVECTOR(x+d, y-d, z-d);  c4 =  jet_density(this, rvec)
       rvec = OCTALVECTOR(x-d, y+d, z+d);  c5 =  jet_density(this, rvec)
       rvec = OCTALVECTOR(x-d, y+d, z-d);  c6 =  jet_density(this, rvec)
       rvec = OCTALVECTOR(x-d, y-d, z+d);  c7 =  jet_density(this, rvec)
       rvec = OCTALVECTOR(x-d, y-d, z-d);  c8 =  jet_density(this, rvec)
       
       !
       ave_jet_density =  &
            (c0 + c1+ c2 + c3 + c4 + c5 + c6 + c7 + c8)/9.0d0
    else ! 2D amr
       ! assigning density 
       rvec = OCTALVECTOR(x, y, z);        c0 =  jet_density(this, rvec)    
       rvec = OCTALVECTOR(x+d, y, z+d);  c1 =  jet_density(this, rvec)
       rvec = OCTALVECTOR(x+d, y, z-d);  c2 =  jet_density(this, rvec)
       rvec = OCTALVECTOR(x-d, y, z+d);  c3 =  jet_density(this, rvec)
       rvec = OCTALVECTOR(x-d, y, z-d);  c4 =  jet_density(this, rvec)       
       !
       ave_jet_density = (c0 + c1+ c2 + c3 + c4)/5.0d0
    end if

  end function ave_jet_density


    

  !
  ! Recursively turn off the inflow flag to false if grid > Rh.
  !
  RECURSIVE SUBROUTINE turn_off_jet(thisOctal,grid, this)    
    
    IMPLICIT NONE
    
    TYPE(octal), POINTER   :: thisOctal
    TYPE(gridtype)         :: grid
    TYPE(jet), INTENT(IN)  :: this
    
    TYPE(octal), POINTER   :: pChild
    
    INTEGER :: subcell, n

    if (thisOctal%threeD) then
       n = 8
    else
       n =4
    end if
    

    do subcell = 1, n
       if (thisOctal%hasChild(subcell)) then
          ! just decdend the tree branch
          pChild => thisOctal%child(subcell)
          CALL turn_off_jet(pChild,grid, this)
       else
          ! turnning it off 
          if (in_jet_flow(this, thisOctal%centre)) thisOctal%inFlow(subcell) = .false.
       end if
    end do
    
  END SUBROUTINE turn_off_jet


  !
  !
  ! Recursively turn off the inflow flag to false if grid > Rh.
  !
  RECURSIVE SUBROUTINE turn_on_jet(thisOctal,grid, this)  
    
    IMPLICIT NONE

    TYPE(octal), POINTER   :: thisOctal
    TYPE(gridtype)         :: grid
    TYPE(jet), INTENT(IN)  :: this

    TYPE(octal), POINTER   :: pChild
  
    INTEGER :: subcell, n

    if (thisOctal%threeD) then
       n = 8
    else
       n =4
    end if

     
    do subcell = 1, n
       if (thisOctal%hasChild(subcell)) then
          ! just decdend the tree branch
          pChild => thisOctal%child(subcell)
          CALL turn_on_jet(pChild,grid, this)
       else
          ! turnning it on
          if (in_jet_flow(this, thisOctal%centre)) thisOctal%inFlow(subcell) = .true.
       end if
    end do
    
  END SUBROUTINE turn_on_jet


  !
  ! Function to check if a given point is in jet flow
  ! 
  function in_jet_flow(this, rvec) RESULT(out)
    implicit none
    logical :: out
    !
    type(jet), intent(in) :: this
    type(octalvector), intent(in) :: rvec
    !
    real(double) ::  w, r, x, y, z
    real(double) ::  cos_theta
    
    x = rvec%x; y =rvec%y; z = rvec%z

    w = x*x + y*y
    r = SQRT(w + z*z)
    w = SQRT(w)

    if (r < this%Rin) then
       out = .false.
    elseif (r > this%Rout) then
       out = .false.
    else
       cos_theta = z/r  ! should be between 0-PI       

       if (ABS(cos_theta) > ABS(COS(this%theta_j))) then
          out = .true.
       else
          out = .false.
       end if
    end if
    
  end function in_jet_flow



  !
  !
  !
  !
  SUBROUTINE fill_velocity_corners(this, thisOctal)
    ! store the velocity values at the subcell corners of an octal so
    !   that they can be used for interpolation.

    IMPLICIT NONE
    type(jet), intent(in) :: this
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
       thisOctal%cornerVelocity(1) = jet_velocity(this, OCTALVECTOR(x1,y1,z1))
       thisOctal%cornerVelocity(2) = jet_velocity(this, OCTALVECTOR(x2,y1,z1))
       thisOctal%cornerVelocity(3) = jet_velocity(this, OCTALVECTOR(x3,y1,z1))
       thisOctal%cornerVelocity(4) = jet_velocity(this, OCTALVECTOR(x1,y2,z1))
       thisOctal%cornerVelocity(5) = jet_velocity(this, OCTALVECTOR(x2,y2,z1))
       thisOctal%cornerVelocity(6) = jet_velocity(this, OCTALVECTOR(x3,y2,z1))
       thisOctal%cornerVelocity(7) = jet_velocity(this, OCTALVECTOR(x1,y3,z1))
       thisOctal%cornerVelocity(8) = jet_velocity(this, OCTALVECTOR(x2,y3,z1))
       thisOctal%cornerVelocity(9) = jet_velocity(this, OCTALVECTOR(x3,y3,z1))

       ! middle level
  
       thisOctal%cornerVelocity(10) = jet_velocity(this, OCTALVECTOR(x1,y1,z2))
       thisOctal%cornerVelocity(11) = jet_velocity(this, OCTALVECTOR(x2,y1,z2))
       thisOctal%cornerVelocity(12) = jet_velocity(this, OCTALVECTOR(x3,y1,z2))
       thisOctal%cornerVelocity(13) = jet_velocity(this, OCTALVECTOR(x1,y2,z2))
       thisOctal%cornerVelocity(14) = jet_velocity(this, OCTALVECTOR(x2,y2,z2))
       thisOctal%cornerVelocity(15) = jet_velocity(this, OCTALVECTOR(x3,y2,z2))
       thisOctal%cornerVelocity(16) = jet_velocity(this, OCTALVECTOR(x1,y3,z2))
       thisOctal%cornerVelocity(17) = jet_velocity(this, OCTALVECTOR(x2,y3,z2))
       thisOctal%cornerVelocity(18) = jet_velocity(this, OCTALVECTOR(x3,y3,z2))
       
       ! top level
       
       thisOctal%cornerVelocity(19) = jet_velocity(this, OCTALVECTOR(x1,y1,z3))
       thisOctal%cornerVelocity(20) = jet_velocity(this, OCTALVECTOR(x2,y1,z3))
       thisOctal%cornerVelocity(21) = jet_velocity(this, OCTALVECTOR(x3,y1,z3))
       thisOctal%cornerVelocity(22) = jet_velocity(this, OCTALVECTOR(x1,y2,z3))
       thisOctal%cornerVelocity(23) = jet_velocity(this, OCTALVECTOR(x2,y2,z3))
       thisOctal%cornerVelocity(24) = jet_velocity(this, OCTALVECTOR(x3,y2,z3))
       thisOctal%cornerVelocity(25) = jet_velocity(this, OCTALVECTOR(x1,y3,z3))
       thisOctal%cornerVelocity(26) = jet_velocity(this, OCTALVECTOR(x2,y3,z3))
       thisOctal%cornerVelocity(27) = jet_velocity(this, OCTALVECTOR(x3,y3,z3))

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
       thisOctal%cornerVelocity(1) = jet_velocity(this, OCTALVECTOR(x1,y1,z1))
       thisOctal%cornerVelocity(2) = jet_velocity(this, OCTALVECTOR(x2,y1,z1))
       thisOctal%cornerVelocity(3) = jet_velocity(this, OCTALVECTOR(x3,y1,z1))
       thisOctal%cornerVelocity(4) = jet_velocity(this, OCTALVECTOR(x1,y2,z2))
       thisOctal%cornerVelocity(5) = jet_velocity(this, OCTALVECTOR(x2,y2,z2))
       thisOctal%cornerVelocity(6) = jet_velocity(this, OCTALVECTOR(x3,y2,z2))
       thisOctal%cornerVelocity(7) = jet_velocity(this, OCTALVECTOR(x1,y3,z3))
       thisOctal%cornerVelocity(8) = jet_velocity(this, OCTALVECTOR(x2,y3,z3))
       thisOctal%cornerVelocity(9) = jet_velocity(this, OCTALVECTOR(x3,y3,z3))

    end if

    
  END SUBROUTINE fill_velocity_corners




end module jet_class
