module cluster_class

  !
  ! Class definition and method for stellar cluster (multiple sources).
  !
  
  use sph_data_class
  use source_mod
  use spectrum_mod
  use octal_mod
  use gridtype_mod
  use isochrone_class
  
  public :: new,  &
       &    kill_all, &
       &    kill_only_stars, &
       &    get_boxsize, &
       &    get_nstar, &
       &    put_a_star, &
       &    get_a_star, &
       &    build_cluster, &
       &    assign_grid_values, &
       &    update_particle_list, &
       &    assign_density, &
       &    find_n_particle_in_subcell

  
  private :: int_cluster, int_cluster_alt
  
  
  type cluster
     private   ! It is better to be private! (The maintenace is easier later...)
     ! size of the cubic box which hold the entire sph particle and stars.
     double precision :: boxsize  
     ! number of stars
     integer :: nstar
     ! stars
     type(sourcetype), pointer :: stars(:)  ! should be allocated with size=nstar
     !
  end type cluster
  
  

  ! overloading functions and creating the interface
  interface new
     module procedure int_cluster
     module procedure int_cluster_alt
  end interface
  

contains
  

  !==============================================================
  ! Constructor
  !==============================================================
  
  ! Initialize the object with default values.  

  subroutine int_cluster(this)
    implicit none
    type(cluster), intent(inout) :: this
    
    this%boxsize = 1
    this%nstar= 1

    NULLIFY(this%stars)
    
  end subroutine int_cluster

  
  !
  ! Initializing with sph data
  !
  subroutine int_cluster_alt(this, sphData, amrGridSize)
    implicit none
    type(cluster), intent(inout) :: this
    type(sph_data), intent(in) :: sphData
    double precision, intent(in) :: amrGridSize

    ! Store important info from the data.    
    this%boxsize = amrGridSize  ! [10^10 cm]
    !
    this%nstar = get_nptmass(sphData)

    ! Allocate star array
    ALLOCATE(this%stars(this%nstar))

  end subroutine int_cluster_alt


  !==============================================================
  ! Destructor
  !==============================================================

  !
  !
  subroutine kill_all(this)
    implicit none
    type(cluster),intent(inout) :: this

    this%nstar=0
    this%boxsize=0.0d0
    DEALLOCATE(this%stars)
    NULLIFY(this%stars)

  end subroutine kill_all

  !
  ! Deallocates the stars array, but keeps other infomations
  subroutine kill_only_stars(this)
    implicit none
    type(cluster),intent(inout) :: this

    DEALLOCATE(this%stars)
    NULLIFY(this%stars)
    
  end subroutine kill_only_stars


  !==============================================================
  ! Accessors
  !==============================================================
  
  !
  ! Return the size of the box which holds the whole particles
  function get_boxsize(this) RESULT(out)
    implicit none
    double precision :: out
    type(cluster), intent(in) :: this
    
    out = this%boxsize
    
  end function get_boxsize
  

  !
  ! Returns the number of stars
  !
  function get_nstar(this) RESULT(out)
    implicit none
    integer :: out
    type(cluster), intent(in) :: this

    out = this%nstar

  end function get_nstar


  !
  ! assigns a star to i-th position of stars array
  ! 
  subroutine put_a_star(this, a_star, i)
    implicit none
    type(cluster), intent(inout) :: this
    type(sourcetype), intent(in) :: a_star
    integer, intent(in) :: i
        
    this%stars(i) = a_star

  end subroutine put_a_star
    

  !
  ! Retrive a star stored in the i-th position of stars array.
  !
  function get_a_star(this, i) RESULT(out)
    implicit none
    type(sourcetype) :: out
    type(cluster), intent(in) :: this
    integer, intent(in) :: i

    ! Needs casting of a pointer???
    out = this%stars(i)

  end function get_a_star
    

  !==============================================================
  ! Tools 
  !==============================================================
  




  !
  ! Creates a cluster fully by assiging T, R, L and so on to each star
  ! 

  subroutine build_cluster(this, sphData, lambda_beg, lambda_end, iso_data)
    implicit none
    type(cluster), intent(inout) :: this
    type(sph_data), intent(in) :: sphData
    ! starting and the ending wavelength in Angstrome!
    double precision, intent(in) :: lambda_beg, lambda_end 
    !
    ! isochrone data to assign R, L, & T os stars
    type(isochrone), intent(in) :: iso_data  
    
    !
    integer :: i, n
    type(sourcetype) :: a_star
    double precision :: mass, age, radius, temperature
    double precision :: luminosity
    double precision :: length_units
    double precision :: x, y, z
    
    n = get_nstar(this)

    ! using a function sph_data_class
    length_units = get_udist(sphData)/1.0d10 ! [10^10 cm]
    
    age = get_time(sphData)*get_utime(sphData) ! in [sec]
    ! converting to years
    age = age/3.1536000d7   ! [years]
    ! assigning stars
    do i = 1, n
       ! preparing the spectrum
       !
       mass = get_pt_mass(sphData, i)*get_umass(sphData) ! [g]
       ! converting it to solar masses
       mass = mass/1.989e33 ! [M_sun]
      

       ! finding the corresponding radius and the surface temperature of
       ! this star.
       ! -- using a routine in isochrone_class
       call  mass_age_to_temp_rad_lum(iso_data, mass, age, temperature, radius, luminosity)

       ! position of this star
       call get_position_pt_mass(sphData, i, x, y, z)
       x = x*length_units; y=y*length_units; z=z*length_units  ! [10^10cm]
       
       ! save it to a star
       a_star%luminosity = luminosity
       a_star%radius = radius
       a_star%teff = temperature
       a_star%position = Vector(x, y, z)

       ! assiging a specturm to this star
       ! ( For now it is just BB specturm)
       ! using the routine in spectrum_mod.f90
       call fillSpectrumBB(a_star%spectrum, temperature, lambda_beg, lambda_end, 1000)
       call normalizedSpectrum(a_star%spectrum)
       
       ! Put this star in a array.
       ! Make sure that the arrays has been already allocated!
       this%stars(i) = a_star

    end do
       

  end subroutine build_cluster



  !
  ! Initialize some parameters in grid
  !
  subroutine initClusterAMR(grid)
    implicit none    
    type(GRIDTYPE), intent(inout)  :: grid
    
    grid%lineEmission = .false.
    
  end subroutine initClusterAMR



  !
  !
  ! 
  !
  !

  
  
  !
  ! Assigns various values which are needed for the run.
  ! This routine should be called only after the grid has been constructed, namely
  ! only in finishGrid routine in amr_mod.f90
  !
  subroutine assign_grid_values(thisOctal,subcell, grid)
    IMPLICIT NONE
    
    TYPE(octal), intent(inout) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    
    thisOctal%temperature(subcell) = 10.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D(subcell) = 1.0
    thisOctal%biasLine3D(subcell) = 1.0
    thisOctal%etaLine(subcell) = 1.e-30
    thisOctal%etaCont(subcell) = 1.e-30
    ! Make sure this is true oterwise it won't be calculated
    ! in Lucy's radiative equiiliburium rouinine!!!    

    if (thisOctal%rho(subcell) > 1.d-20) then
       thisOctal%inFlow(subcell) = .true.  
    else
       thisOctal%inFlow(subcell) = .false.
    endif
       
  end subroutine assign_grid_values
  


  
  !
  ! Adds density values to ths subcell of thisOctal.
  !
  ! This routine can be used for example, in addNewChild and initFirstOctal
  ! in amr_mod.f90
  !
  subroutine assign_density(thisOctal,subcell, sphData, geometry)
    implicit none
    type(octal), intent(inout) :: thisOctal
    integer, intent(in) :: subcell
    !type(gridtype), intent(in) :: grid
    type(sph_data), intent(in) :: sphData
    character(len=*), intent(in) :: geometry
    !
    double precision :: density_ave
    integer :: nparticle    
    !
    !
    ! using the function in amr_mod.f90
    call find_n_particle_in_subcell(nparticle, density_ave, sphData, thisOctal, subcell)
    ! assign density to the subcell
    select case (geometry)
       case ("cluster")
          thisOctal%rho(subcell) = density_ave
          ! This should be in [g/cm^3]
       case("wr104")
          thisOctal%rho(subcell) = max(1.d-30,real(nparticle)/ thisOctal%subcellsize**3)
    end select

  end subroutine assign_density
  
    

  
  !
  ! Updates the sph particle linked list.
  !
  subroutine update_particle_list(thisOctal, subcell, newChildIndex, sphData)
    implicit none
    type(octal), pointer :: thisOctal
    integer, intent(in) :: subcell       ! do loop index of thisOctal
    integer, intent(in) :: newChildIndex ! new child index
    !type(gridtype), intent(in) :: grid
    type(sph_data), intent(in) :: sphData
    !
    integer          :: np, i, j, k
    double precision :: udist  ! units for distance
    double precision :: x, y, z
    double precision :: density_ave
    integer :: nparticle
    
    !
    ! Finding the number of particles in this subcell first
    call find_n_particle_in_subcell(nparticle, density_ave, sphData, thisOctal, subcell)
    ! Allocate the memory for new particle list of the child
    ALLOCATE(thisOctal%child(newChildIndex)%gas_particle_list(nparticle))

    ! # of particles in thisOctal
    np = SIZE(thisOctal%gas_particle_list)

    ! Units of length
    udist = get_udist(sphData) / 1.0d10  ! [10^10cm]
       
    k=0
    do i = 1, np
       ! extract the SPH index of particles in the parent cell.
       j = thisOctal%gas_particle_list(i)
       
       ! copy this value if this particle is in this child.
       ! Using a routine in sph_data_class. 
       call get_position_gas_particle(sphData, j, x, y, z)
       
       ! convert units
       x = x*udist; y = y*udist;  z = z*udist  ! [10^10cm]
       
       ! quick check to see if this gas particle
       ! belongs to this child cell.
       if ( within_subcell(thisOctal, subcell, x, y, z) ) then
          k = k+1
          ! add this particle index to this child using
          ! a routine in linked_list_class.f90
          thisOctal%child(newChildIndex)%gas_particle_list(k) = j
       end if
    end do

    
  end subroutine update_particle_list

  


  
  ! For a given octal object and sph_data_class object, this
  ! function returns the number of gas particles which belongs to a subcell this
  ! of this octal, and the average denisty ( in g/cm^3) of this cell.
  !  
  subroutine find_n_particle_in_subcell(n, rho_ave, sphData, node, subcell)
    implicit none
    integer, intent(out) :: n                ! number of particels in the subcell
    double precision, intent(out) :: rho_ave ! average density of the subcell
    type(sph_data), intent(in)    :: sphData
    type(octal), intent(inout)    :: node
    integer, intent(in)           :: subcell ! index of the subcell
    !
    integer :: npart
    integer :: i, j, counter
    double precision :: x, y, z
    !
    !
    double precision, save :: umass, udist, udent  ! for units conversion   
    double precision, save :: rho_min
    logical, save  :: first_time = .true.
    !
    logical :: restart
    
    ! Carry out the initial calculations
    if (first_time) then       
       ! units of sphData
       umass = get_umass(sphData)  ! [g]
       udist = get_udist(sphData)  ! [cm]
       udent = umass/udist**3

       ! convert units
       udist = udist/1.0d10  ! [10^10cm]

       rho_min = get_rhon_min(sphData)
       first_time = .false.
    end if

    
    ! retriveing the number of the total gas particles
    ! in "node".
    ! -- using the function in linked_list_class.f90
    npart = SIZE(node%gas_particle_list)
    
    counter = 0
    rho_ave = 0.0d0
    
    do i=1, npart
       ! Retriving the sph data index for this paritcle
       j = node%gas_particle_list(i)

       ! retriving this posisiton of the gas particle.
       call get_position_gas_particle(sphData, j, x, y, z)
       
       ! convert units
       x = x*udist; y = y*udist; z = z*udist   ! [10^10cm]
       
       ! quick check to see if this gas particle is
       ! belongs to this cell.
       if ( within_subcell(node, subcell, x, y, z) ) then
          counter = counter + 1
          rho_ave = rho_ave + get_rhon(sphData, j) 
       end if
       
    end do

    n = counter

    
    if (n>0) then
       rho_ave = rho_ave/dble(n)
    else
       rho_ave = rho_min*1.0d-1
    end if
    
    rho_ave = rho_ave*udent  ! [g/cm^3]

    rho_ave = rho_ave*1.0d13 ! just for debug

    
    
  end subroutine find_n_particle_in_subcell


  
 end module cluster_class
