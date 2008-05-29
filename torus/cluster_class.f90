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
  use messages_mod
  use vector_mod
  
  implicit none

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
       &    find_n_particle_in_subcell, &
       &    n_stars_in_octal, &
       &    fill_in_empty_octals, &
       &    remove_too_close_cells, &
       &    disc_density, &
       &    disc_intersects_subcell, &
       &    cell_too_close_to_star, &
       &    average_disc_density, &
       &    average_disc_density_fast, &
       &    max_disc_density_from_sample, &
       &    compute_scale_size, &
       &    include_disc, &
       &    assign_emission_bias, &
       &    restrict, &
       &    reassign_10K_temperature

  
  
  private :: int_cluster, int_cluster_alt,  compute_emission_bias
  
  
  type cluster
     private   ! It is better to be private! (The maintenance is easier later...)
     ! size of the cubic box which hold the entire sph particle and stars.
     real(double) :: boxsize  
     ! number of stars
     integer :: nstar
     ! stars
     type(sourcetype), pointer :: stars(:)  ! should be allocated with size=nstar
     !
     ! Flag to include accreation disc or not...
     logical :: disc_on  
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
    this% disc_on = .false.
    
  end subroutine int_cluster

  
  !
  ! Initializing with sph data
  !
  subroutine int_cluster_alt(this, sphData, amrGridSize, disc_on)
    implicit none
    type(cluster), intent(inout) :: this
    type(sph_data), intent(in) :: sphData
    real(double), intent(in) :: amrGridSize
    logical, intent(in) :: disc_on
    ! Store important info from the data.    
    this%boxsize = amrGridSize  ! [10^10 cm]
    !
    this%nstar = get_nptmass(sphData)

    ! Allocate star array
    ALLOCATE(this%stars(this%nstar))

    !
    this%disc_on = disc_on

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
    real(double) :: out
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
    real(double), intent(in) :: lambda_beg, lambda_end 
    !
    ! isochrone data to assign R, L, & T os stars
    type(isochrone), intent(in) :: iso_data  
    
    !
    integer :: i, n
    type(sourcetype) :: a_star
    real(double) :: mass, age, radius, temperature
    real(double) :: luminosity
    real(double) :: length_units
    real(double) :: x, y, z
    
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
       a_star%luminosity = luminosity     ! erg/sec
!       a_star%luminosity = luminosity * 100.    ! erg/sec
!       write(*,*) "Luminosity fudged"
       a_star%radius = radius             ! 10^10cm
       a_star%teff = temperature          ! Kelvins
       a_star%position = DoubleVector(x, y, z)  ! 10^10cm

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
    
!    thisOctal%temperature(subcell) = 3.0e0
!    thisOctal%temperature(subcell) = 10.0
!    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%etaLine(subcell) = 1.e-30
    thisOctal%etaCont(subcell) = 1.e-30
       
  end subroutine assign_grid_values



  
  !
  ! Adds density values to ths subcell of thisOctal.
  !
  ! This routine can be used for example, in addNewChild and initFirstOctal
  ! in amr_mod.f90
  !
  subroutine assign_density(thisOctal,subcell, sphData, grid, a_cluster)
    implicit none
    type(octal), intent(inout) :: thisOctal
    integer, intent(in) :: subcell
    type(gridtype), intent(in) :: grid
    type(sph_data), intent(in) :: sphData
    character(len=20) :: geometry
    type(cluster), intent(in), optional :: a_cluster
    !
    real(double) :: density_ave, rho_disc_ave, dummy_d
    real(double), parameter :: density_crit = 1d-13
    type(VECTOR) :: velocity_ave
    real :: tem_ave
    integer :: nparticle     
    !
    !
    ! assign density to the subcell
    
    geometry = grid%geometry

    select case (geometry)
       case ("cluster")
          call find_temp_in_subcell(nparticle, tem_ave, sphData, &
               thisOctal, subcell)
          thisOctal%temperature(subcell)  = tem_ave
          call find_n_particle_in_subcell(nparticle, density_ave, sphData, &
               thisOctal, subcell)

          ! if this subcell intersect with any of the stellar disk
          ! in the cluster, then add the contribution from the stellar disc.
          rho_disc_ave = 1.d-30

         if (a_cluster%disc_on) then
             if (stellar_disc_exists(sphData) .and.  &       
                  disc_intersects_subcell(a_cluster, sphData, thisOctal, subcell) ) then

                rho_disc_ave = average_disc_density_fast(sphData, thisOctal, &
                     subcell, a_cluster, dummy_d)

!             rho_disc_ave = average_disc_density(sphData, thisOctal, &
!                  subcell, a_cluster)
                thisOctal%rho(subcell) = density_ave + rho_disc_ave

             ! This should be in [g/cm^3]

!             if (thisOctal%rho(subcell) > 5.d-14) then
                if (thisOctal%rho(subcell) > 1.d-20) then
                   thisOctal%inFlow(subcell) = .true.  
                else
                   thisOctal%inFlow(subcell) = .true.
                endif
             else

                thisOctal%rho(subcell) = density_ave
                ! This should be in [g/cm^3]
                
                if (thisOctal%rho(subcell) > 1.d-27) then
                   thisOctal%inFlow(subcell) = .true.  
                else
                   thisOctal%inFlow(subcell) = .true.
                endif
             end if
          else

             thisOctal%rho(subcell) = density_ave
             ! This should be in [g/cm^3]
             
             if (thisOctal%rho(subcell) > 1.d-27) then
                thisOctal%inFlow(subcell) = .true.  
             else
                thisOctal%inFlow(subcell) = .true.
             endif

          end if


          !!
          !! REMOVE THIS LATER !
          !!
          !!thisOctal%rho(subcell) =  1.0e-28 ! just for testing...
          !!
          !!
          !!

       case ("molcluster")
          call find_temp_in_subcell(nparticle, tem_ave, sphData, &
               thisOctal, subcell)
          thisOctal%temperature(subcell)  = tem_ave
!          call find_n_particle_in_subcell(nparticle, density_ave, sphData, &
!               thisOctal, subcell)
          call find_density(nparticle, density_ave, sphData, &
               thisOctal, subcell) ! DAR routine based on mass in volume not smoothing length averge density

                thisOctal%rho(subcell) = density_ave             
                
                thisOctal%temperature(subcell) = max(10., 10. * ((density_ave / density_crit)**(0.4)))

          call find_velocity(nparticle, velocity_ave, sphData, &
               thisOctal, subcell) ! DAR routine basedon average momentum

          thisOctal%velocity(subcell) = velocity_ave

       case("wr104")
          call find_n_particle_in_subcell(nparticle, density_ave, sphData, &
               thisOctal, subcell)          
          thisOctal%rho(subcell) = max(1.d-30,dble(nparticle)/ cellVolume(thisOctal, subcell))
    end select

  end subroutine assign_density
  
    

  
  !
  ! Updates the sph particle linked list.
  !
  subroutine update_particle_list(thisOctal, subcell, newChildIndex, sphData)
    implicit none
    TYPE(octal), INTENT(INOUT) :: thisOctal ! the parent octal 
!    type(octal), pointer :: thisOctal
    integer, intent(in) :: subcell       ! do loop index of thisOctal
    integer, intent(in) :: newChildIndex ! new child index
    !type(gridtype), intent(in) :: grid
    type(sph_data), intent(in) :: sphData
    !
    integer          :: np, i, j, k
    real(double) :: udist  ! units for distance
    real(double) :: x, y, z
    real(double) :: density_ave
    integer :: nparticle
    
    !
    ! Finding the number of particles in this subcell first
    call find_n_particle_in_subcell(nparticle, density_ave, sphData,&
         thisOctal, subcell)

    if (nParticle > 0) then
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
    endif
    
  end subroutine update_particle_list

  


  
  ! For a given octal object and sph_data_class object, this
  ! function returns the number of gas particles which belongs to a subcell this
  ! of this octal, and the average denisty ( in g/cm^3) of this cell.
  !  
  subroutine find_n_particle_in_subcell(n, rho_ave, sphData, node, subcell)
    implicit none
    integer, intent(out) :: n                ! number of particels in the subcell
    real(double), intent(out) :: rho_ave ! average density of the subcell
    type(sph_data), intent(in)    :: sphData
    type(octal), intent(inout)    :: node
    integer, intent(in)           :: subcell ! index of the subcell
    !
    integer :: npart
    integer :: i, j, counter
    real(double) :: x, y, z
    !
    !
    real(double), save :: umass, udist, udent  ! for units conversion  
    real(double), save :: rho_min
    logical, save  :: first_time = .true.
    !
    !logical :: restart
    !    
    integer, parameter :: nsample =400

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

    if (associated(node%gas_particle_list)) then
       
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
       
    else

       n = 0

    endif
    
    if (n>0) then
       rho_ave = rho_ave/dble(n)
       rho_ave = rho_ave*udent  ! [g/cm^3]
    else
       rho_ave = 1.d-30
    end if
    

  end subroutine find_n_particle_in_subcell

  subroutine find_density(n, dens_ave, sphData, node, subcell)
    implicit none
    integer, intent(out) :: n                ! number of particels in the subcell
    real(double), intent(out) :: dens_ave ! average density of the subcell
    type(sph_data), intent(in)    :: sphData
    type(octal), intent(inout)    :: node
    integer, intent(in)           :: subcell ! index of the subcell
    !
    integer :: npart
    integer :: i, j, counter
    real(double) :: x, y, z
    !
    !
    real(double), save :: umass, udist, udens  ! for units conversion  
    real(double), save :: dens_min 
    logical, save  :: first_time = .true.
    !
    !logical :: restart
    !    
    integer, parameter :: nsample =400

    ! Carry out the initial calculations
    if (first_time) then       
       ! units of sphData
       umass = get_umass(sphData)  ! [g]
       udist = get_udist(sphData)  ! [cm]
       udens = umass/udist**3

       ! convert units
       udist = udist/1.0d10  ! [10^10cm]

       dens_min = get_rhon_min(sphData)
       first_time = .false.
    end if

    if (associated(node%gas_particle_list)) then
       
       ! retriving the number of the total gas particles
       ! in "node".
       ! -- using the function in linked_list_class.f90
       npart = SIZE(node%gas_particle_list)
       
       counter = 0
       dens_ave = 0.0d0
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
             dens_ave = dens_ave + get_mass(sphData, j) 
          end if
          
       end do
       
       n = counter
       
    else

       n = 0

    endif
    
    if (n>0) then
!       dens_ave = dens_ave/dble(n)
       dens_ave = dens_ave*umass/(1e30*cellVolume(node, subcell))  ! [g/cm^3]
    else
       dens_ave = 1.d-30
    end if
    

  end subroutine find_density

  subroutine find_velocity(n, vel_ave, sphData, node, subcell) ! technically finding average momentum?
    implicit none
    integer, intent(out) :: n                ! number of particels in the subcell
    type(VECTOR), intent(out) :: vel_ave ! average density of the subcell
    type(sph_data), intent(in)    :: sphData
    type(octal), intent(inout)    :: node
    integer, intent(in)           :: subcell ! index of the subcell
    !
    integer :: npart
    integer :: i, j, counter
    real(double) :: x, y, z
    real(double) :: mass, mass_sum


    real(double), save :: umass, udist, utime  ! for units conversion  
    real(double), save :: dens_min
    logical, save  :: first_time = .true.

    ! Carry out the initial calculations
    if (first_time) then       
       ! units of sphData
       umass = get_umass(sphData)  ! [g]
       udist = get_udist(sphData)  ! [cm]
       utime = get_utime(sphData)  ! [s]

       ! convert units
       udist = udist/1.0d10  ! [10^10cm]

       dens_min = get_rhon_min(sphData)
       first_time = .false.
    end if

    if (associated(node%gas_particle_list)) then
       
       ! retriving the number of the total gas particles
       ! in "node".
       ! -- using the function in linked_list_class.f90
       npart = SIZE(node%gas_particle_list)
       
       counter = 0
       vel_ave = VECTOR(1d-30,1d-30,1d-30)

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
             mass = get_mass(sphdata, j)
             vel_ave = vel_ave + mass * get_vel(sphData, j)
             mass_sum = mass_sum + mass
          end if
          
       end do
       
       n = counter
       
    else

       n = 0

    endif
    
    if (n>0) then
       vel_ave = (1./mass_sum) * vel_ave
       vel_ave = ((1./(utime*3.1556926e7))*udist*1d10) * vel_ave  ! [cm/s]
       vel_ave = (1./cspeed) * vel_ave  ! fraction of c
    else
       vel_ave = VECTOR(1.d-30,1.d-30,1d-30)
    end if
    
  end subroutine find_velocity

  subroutine find_temp_in_subcell(n, tem_ave, sphData, node, subcell)
    use input_variables, only : TMinGlobal
    implicit none
    integer, intent(out) :: n    ! number of particles in the subcell
    real, intent(out) :: tem_ave ! average temperature of the subcell
    type(sph_data), intent(in)    :: sphData
    type(octal), intent(inout)    :: node
    integer, intent(in)           :: subcell ! index of the subcell
    !
    integer :: npart
    integer :: i, j, counter
    real(double) :: x, y, z
    !
    real(double), save :: udist  ! for units conversion  
    logical, save  :: first_time = .true.

    ! Carry out the initial calculations
    if (first_time) then       
       ! units of sphData
       udist = get_udist(sphData)  ! [cm]

       ! convert units
       udist = udist/1.0d10  ! [10^10cm]

       first_time = .false.
    end if

    if (associated(node%gas_particle_list)) then
       
       ! retriveing the number of the total gas particles
       ! in "node".
       ! -- using the function in linked_list_class.f90
       npart = SIZE(node%gas_particle_list)
       
       counter = 0
       tem_ave = 0.0d0
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
             tem_ave = tem_ave + get_temp(sphData, j) 
          end if
          
       end do
       
       n = counter
       
    else

       n = 0

    endif
    
    if (n>0) then
       tem_ave = tem_ave/real(n)
    else
       tem_ave = TMinGlobal
    end if
    

  end subroutine find_temp_in_subcell

  !
  ! For a given octal object, this function finds and returns the number of stars within in 
  ! this octal.
  !

  function n_stars_in_octal(this, an_octal) result(out)
    implicit none
    integer :: out 
    type(cluster), intent(in) :: this 
    type(octal), intent(in) :: an_octal
    
    integer :: n, i 
    
    n = 0    
    do i = 1, this%nstar
       ! using a fucntion in source_mod
       if ( source_within_octal(this%stars(i), an_octal) ) n = n+1
    end do

    out = n
  end function n_stars_in_octal


  !
  ! Fill in  the empty cells (cells in which there is no sph particle), 
  ! if possible, by using the values of sibling (other childrens who have same parent).
  !
  recursive subroutine fill_in_empty_octals(this, thisOctal, sphData)
    implicit none
    
    type(octal), pointer       :: thisOctal
    type(cluster), intent(in)  :: this
    type(sph_data), intent(in) :: sphData
!    type(cluster), optional, intent(in)  :: stellar_cluster
    type(octal), pointer       :: pChild  
    integer                    :: i           ! loop counter
    integer :: np, nc, subcell, parent_subcell
    real :: rho_ave
    real(double) :: rho_tmp    
    
    nc = thisOctal%nChildren    
    do subcell = 1, 8
       if (thisOctal%hasChild(subcell)) then

          if (this%disc_on) then
             if (.not.disc_intersects_subcell(this, sphData, thisOctal, subcell)) then
                ! checks if this subcell is empty or not
                call find_n_particle_in_subcell(np, rho_tmp, sphData, &
                     thisOctal, subcell)
                if (np<=0) then ! apply correction 
                   parent_subcell = whichSubcell(thisOctal%parent, thisOctal%centre)
                   rho_ave = thisOctal%parent%rho(parent_subcell)
                   thisOctal%rho(subcell)  = rho_ave                    
                   if (thisOctal%parent%nDepth < 3) thisOctal%rho(subcell) = 1.d-19
                else
                   ! do nothing 
                   continue
                end if
             end if
          else
             ! checks if this subcell is empty or not
             call find_n_particle_in_subcell(np, rho_tmp, sphData, &
                  thisOctal, subcell)
             if (np<=0) then ! apply correction 
                ! inherit the value from a paranet.
                parent_subcell = whichSubcell(thisOctal%parent, thisOctal%centre)
                rho_ave = thisOctal%parent%rho(parent_subcell)
                thisOctal%rho(subcell)  = rho_ave
                if (thisOctal%parent%nDepth < 3) thisOctal%rho(subcell) = 1.d-19
             else
                ! do nothing 
                continue
             end if
          end if


          ! find the child
          do i = 1, nc
             if (thisOctal%indexChild(i) == subcell) then
                pChild => thisOctal%child(i)                
                call fill_in_empty_octals(this, pChild, sphData)
                exit
             end if
          end do
       else
          if (this%disc_on) then
             if (.not.disc_intersects_subcell(this, sphData, thisOctal, subcell)) then
                ! checks if this subcell is empty or not
                call find_n_particle_in_subcell(np, rho_tmp, sphData, &
                     thisOctal, subcell)
                if (np<=0) then ! apply correction 
                   parent_subcell = whichSubcell(thisOctal%parent, thisOctal%centre)
                   rho_ave = thisOctal%parent%rho(parent_subcell)
                   thisOctal%rho(subcell)  = rho_ave                   
                   if (thisOctal%parent%nDepth < 3) thisOctal%rho(subcell) = 1.d-19
                else
                   ! do nothing 
                   continue
                end if
             end if
          else
             ! checks if this subcell is empty or not
             call find_n_particle_in_subcell(np, rho_tmp, sphData, &
                  thisOctal, subcell)
             if (np<=0) then ! apply correction 
                ! inherit the value from a paranet.
                parent_subcell = whichSubcell(thisOctal%parent, thisOctal%centre)
                rho_ave = thisOctal%parent%rho(parent_subcell)
                thisOctal%rho(subcell)  = rho_ave
                if (thisOctal%parent%nDepth < 3) thisOctal%rho(subcell) = 1.d-19
             else
                ! do nothing 
                continue
             end if
          end if
       end if
    end do
  
    contains

      FUNCTION whichSubcell(thisOctal,point) RESULT (subcell)
        ! returns the identification number (1-8) of the subcell of the 
        ! current octal which contains a given point
        ! NB this does NOT check that the point lies within the bounds of the octal!
        
        IMPLICIT NONE
        
        TYPE(octal), INTENT(IN)       :: thisOctal
        TYPE(octalVector), INTENT(IN) :: point
        INTEGER                       :: subcell
        
        IF ( point%x < thisOctal%centre%x ) THEN
           IF ( point%y < thisOctal%centre%y ) THEN
              IF ( point%z < thisOctal%centre%z ) THEN
                 subcell = 1
              ELSE 
                 subcell = 5
                 
              ENDIF
           ELSE 
              IF (point%z < thisOctal%centre%z) THEN
                 subcell = 3
              ELSE 
                 subcell = 7
              ENDIF
           END IF
        ELSE
           IF (point%y < thisOctal%centre%y) THEN
              IF (point%z < thisOctal%centre%z) THEN
                 subcell = 2
              ELSE 
                 subcell = 6
              ENDIF
           ELSE 
              IF (point%z < thisOctal%centre%z) THEN
                 subcell = 4
              ELSE 
                 subcell = 8
              ENDIF
           END IF
        ENDIF
        
      END FUNCTION whichSubcell
      
    
  end subroutine fill_in_empty_octals





  !
  ! Removes cells (just setting the inflow flag to .false. if the cells are closer than a 
  ! given radius (in 10^10cm) from star positions specified in cluster_class object.
  recursive subroutine remove_too_close_cells(thisCluster, thisOctal, R_max)
    implicit none

    type(cluster), intent(in)  :: thisCluster    
    type(octal), pointer       :: thisOctal
    real(double), intent(in) :: R_max   ! [10^10 cm]
    !
    type(octal), pointer       :: pChild
    integer                    :: i           ! loop counter
    integer :: nc, subcell
    
    nc = thisOctal%nChildren    
    do subcell = 1, 8
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, nc
             if (thisOctal%indexChild(i) == subcell) then
                pChild => thisOctal%child(i)                
                call remove_too_close_cells(thisCluster, pChild, R_max)
                exit
             end if
          end do
       else
          if (thisCluster%disc_on) then
             continue
          else
             ! checks if it's too close 
             if (cell_too_close_to_star(thisCluster,thisOctal, subcell, R_max) ) then
                thisOctal%inFlow(subcell) = .false.
                thisOctal%rho(subcell) = 1.e-28  ! set it to a very low density
             else
                ! no correction done
                continue             
             end if
          end if
       end if
    end do
    
    
  end subroutine remove_too_close_cells




  !
  ! This routine write the catalog of the stars in a cluster in 
  ! "catalog.dat".
  !
  subroutine write_catalog(this, sphData, filename)
    implicit none
    type(cluster),  intent(in) :: this
    type(sph_data), intent(in) :: sphData
    character(len=*), optional :: filename
    
    integer, parameter :: LUOF=23
    integer :: i
    type(sourcetype) :: a_star
    real(double) :: mass
    real(double), parameter :: M_sun = 1.989e33         ! grams
    real(double), parameter :: L_sun = 3.85d33          !erg/s
    real(double), parameter :: R_sun = 6.96d0           ! 10^10cm

    ! 


    if ( present(filename) ) then
       call writeInfo("Stellar Cluster Catalog written in "//filename//"...",TRIVIAL)
       open(unit=LUOF, file=filename, status="replace")
    else
       call writeInfo("Stellar Cluster Catalog written in catalog.dat...",TRIVIAL)
       open(unit=LUOF, file="catalog.dat", status="replace")
    end if

    ! writing the header
    write(LUOF,'(a)')       "#  This file was written by cluster::write_catalog."
    write(LUOF,'(a)')       "#" 
    write(LUOF,'(a,1x,i15)') "# nstar =", this%nstar
    write(LUOF,'(a,1x,1PE12.3,1x,a10)') &
         &                  "# age   =", get_time(sphData)*get_utime(sphData)/3.1536000d7, " years old"
    write(LUOF,'(a)')       "#" 
    write(LUOF,'(a)')       "#" 
    write(LUOF,11)          "# ", "index", "x[10^cm]", "y[10^cm]", "z[10^cm]", &
         "mass [Msun]", "L [Lsun]", "radius  [Rsun]", "temperature [K]"
11  format(a,1x,a5, 7(2x,a15))
12  format(2x,i5, 7(1x,1PE16.7))
    write(LUOF, '(a)')      "#-----------------------------------------------"//&
         &"---------------------------------------------------------------------------------"
    do i = 1, this%nstar
       a_star = get_a_star(this,i)
       mass = get_pt_mass(sphData, i)*get_umass(sphData) ! [g]
       write(LUOF,12) i, a_star%position%x, a_star%position%y, a_star%position%z, &
            mass/M_sun, a_star%luminosity/L_sun, a_star%radius/R_sun, a_star%teff
    end do

    close(LUOF)

  end subroutine write_catalog



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
  
  function disc_density(xpos,ypos,zpos,sphData, young_cluster ,sizescale) RESULT(density_out)
    implicit none
    real(double) :: density_out  ! in [g/cm^3]
    real(double), intent(in) :: xpos,ypos,zpos  ! should be in [10^10cm]
    type(sph_data), intent(in)   :: sphData
    type(cluster), intent(in)    :: young_cluster
    real(double), intent(out):: sizescale       ! should be in [10^10cm]
    !
    logical,save :: first_time = .true.
    real(double), parameter :: pi = 3.14159265d0
    integer :: nstar, i 
    real(double) :: holesize
    real(double) :: discrad, discmass, unitx, unity, unitz
    real(double) :: spinx, spiny, spinz, starrad, x, y, z
    real(double) :: xposi, yposi, zposi, rad, radius2, radius
    real(double) :: unorm, zdist, height, densitymidplane
    type(sourcetype) :: a_star
    ! maximum density allowed
!    real(double), parameter :: rho_max = 1.0d-12   ! g/cm^3
    real(double), parameter :: rho_max = 2.5d50   ! g/cm^3 (equivalent to no threshold)
    !
    real(double) :: s2 !, s1  
    logical :: logscale
!    real(double) :: Rgap1 = 1000.0 ! 10^10 cm
!    real(double) :: Rgap2 = 2000.0 ! 10^10 cm
    

    ! quick check
    if (first_time) then
       if (stellar_disc_exists(sphData)) then ! OK
          first_time = .false.
       else
          write(*,*) "Error:: You did not read in the stellar disc data! &
               &[cluster::disc_density]."
          stop
       end if
    end if
    
    !
    !  Inner disc radius is 3 times stellar radius (starrad(i) in unit of R_sun)
    !  required to be in units of 10^10 cm (i.e. times 7x10^10)
    !
    nstar = get_nstar(young_cluster)
!    holesize = 3.0d0 * 7.0d+00  !this is the original choice.
    holesize = 6.0d0 * 7.0d+00  !this is the original choice.
!    holesize = 12.0d0 * 7.0d00   ! 12 times the radius of a star



    !
    density_out = 1.0d-27
    sizescale = 1.0d+10
    DO i = 1, nstar   
       a_star = get_a_star(young_cluster,i)
       
       ! retrive some disc parameters from sphData
       call get_stellar_disc_parameters(sphData, i, discrad, discmass, &
            spinx, spiny, spinz)

       ! radius of this star
       starrad = a_star%radius  ! [10^10cm]
       starrad = starrad/7.0d0  ! [R_sun]
       x = a_star%position%x
       y = a_star%position%y
       z = a_star%position%z
       
       xposi = xpos-x
       yposi = ypos-y
       zposi = zpos-z
       radius2 = xposi**2 + yposi**2 + zposi**2

       radius = SQRT(radius2)
       IF (radius.LT.discrad) THEN          
          IF (radius.LT.1.1d0*(holesize*starrad)) THEN
!          IF (radius.LT.1.1d0*(holesize*starrad) .or. &
!               ( radius < Rgap2 .and. radius >Rgap1 ) ) THEN
             sizescale = holesize*starrad
             IF (radius > (holesize*starrad/2.0d0)) sizescale=sizescale/2.0d0
          ELSE
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

             densitymidplane = 5.0d0*discmass/(8.0d0*pi*discrad**2)*  &
                  rad**(-0.75d0)
             density_out = densitymidplane*   &
                  EXP(-zdist**2/(2.0d0*height**2))/(SQRT(2.0d0*pi)*height)
             if (density_out>1.0d-40) then
!                ! original
!                sizescale=0.1d0*height*(densitymidplane/density_out)**(1.0d0/3.0d0)
!                sizescale=1.0d-2*sizescale

!                ! finding the scale size at this location using 
!                ! the function in this module.

!                s1 = compute_scale_size(radius, holesize*starrad, discrad, &
!                     5.0d0*a_star%radius, 1.0d6, logscale)
!                s2 = s1




!!
!! Produces the 90 Mb grid around star 4
!!
                logscale = .false.
                s2 = compute_scale_size(ABS(zdist),0.0d0,  20.0d0*height, &
                    1.0d0*height,  5.0d0*height, logscale)
                sizescale = s2

!                logscale = .false.
!                s2 = compute_scale_size(ABS(zdist),0.0d0,  20.0d0*height, &
!                    1.0d0*height,  4.0d0*height, logscale)
!                sizescale = s2


! 
! This scaling creates 330Mb grid around star 4 and 660M around  star 46
!  
!                logscale = .false.
!                s2 = compute_scale_size(ABS(zdist),0.0d0,  8.0d0*height, &
!                    2.0d0*height,  1.0d0*height, logscale)
!                sizescale = s2


!!
!! For Tim Naylor
!!                logscale = .false.
!                s2 = compute_scale_size(ABS(zdist),0.0d0,  16.0d0*height, &
!                    4.0d0*height,  2.0d0*height, logscale)
!                sizescale = s2

!                logscale = .false.
!                s2 = compute_scale_size(ABS(zdist),0.0d0,  2.0d0*height, &
!                    1.0d0*height,  0.5d0*height, logscale)
!                sizescale = s2




!
!                ! limit the smallest and largest sizes
                sizescale = MAX(sizescale, a_star%radius*1.0d0)
                sizescale = MIN(sizescale, 2.0d3)                

                 ! else takes the default density from the top.
             end if
          ENDIF
          !  Given position can only be inside 1 disc at most (discs cannot overlap) 
          GOTO 10
       ENDIF
    END DO
10  CONTINUE    

    ! converting units from g/(10^10cm)^3 to g/cm^3

    density_out = density_out/1.0d30   ! [g/cm^3]

    ! assuming dust density is 100 times smaller than the gas density
    density_out = density_out/1.0d2   ! [g/cm^3]


    ! cutoff density
    ! maximum density allowed
    if (density_out>rho_max) density_out = rho_max ! [g/cm^3]
    

!    !
!    ! Emulating the larger grain size....
    !   
!    density_out = density_out/1.0d3   ! [g/cm^3]
!
!    ! limiting the minimum density
!    ! should be removed later
!    if (density_out<1.0d-20) density_out = 1.0d-20   ! [g/cm^3]


  end function disc_density




  !
  ! For a given octal, this routine checks if any of the stellar disk 
  ! in a cluster crosses a subcell of a octal.
  ! As a quick check, we consider the disk as a sphere (with r=Rdisc) 
  ! and the cell as a sphere with r=SQRT(3)*(d/2) where d is the subcell size
  function disc_intersects_subcell(this, sphData, octal_in, subcell) RESULT(out)
    implicit none
    logical :: out 
    
    type(cluster), intent(in) :: this
    type(sph_data), intent(in) :: sphData
    type(octal), intent(in) :: octal_in
    integer, intent(in) :: subcell

    real(double) :: x, y, z, xc, yc, zc
    type(sourcetype) :: a_star
    integer :: i, nstar
    real(double) :: dummy_d, R_disc, R2_disc, R2_cell, d2
    TYPE(octalVector)     :: cellCenter

    nstar = get_nstar(this)
    
    R2_cell = 3.0d0*octal_in%subcellSize*octal_in%subcellSize ! square of the cell radius    

    ! position if the cell
    cellCenter = subcellCentre(octal_in,subcell)
    xc=dble(cellCenter%x); yc=dble(cellCenter%y); zc=dble(cellCenter%z)
   
    out = .false.
    do i = 1, nstar
       a_star = get_a_star(this,i)
       ! retrive some disc parameters from sphData
       call get_stellar_disc_parameters(sphData, i, R_disc, dummy_d, &
            dummy_d, dummy_d, dummy_d)       
       
       ! position of the star
       x = a_star%position%x
       y = a_star%position%y
       z = a_star%position%z
       
       ! The square of the distance between the star and the cell center.
       d2 = (x-xc)**2 + (y-yc)**2 + (z-zc)**2
       
       R2_disc = R_disc**2

       if (d2 < R2_disc+R2_cell) then 
          ! the two sphere intersect
          out = .true.
          exit  ! exit this loop.
       end if          
    end do
     
  end function disc_intersects_subcell




  !
  ! For a given octal, this routine checks if the center of the cell is with
  ! in a given radius from stars.
  ! As a quick check, we consider the disk as a sphere (with r=Rdisc) 
  function cell_too_close_to_star(this,octal_in, subcell, R_max) RESULT(out)
    implicit none
    logical :: out 
    
    type(cluster), intent(in) :: this
    type(octal), intent(in) :: octal_in
    integer, intent(in) :: subcell
    real(double), intent(in) :: R_max  ! limiting distance in [10^10cm]

    real(double) :: x, y, z, xc, yc, zc
    type(sourcetype) :: a_star
    integer :: i, nstar
    real(double) :: R_sq, R_max_sq
    TYPE(octalVector)     :: cellCenter

    nstar = get_nstar(this)
    
    ! position if the cell
    cellCenter = subcellCentre(octal_in,subcell)
    xc=dble(cellCenter%x); yc=dble(cellCenter%y); zc=dble(cellCenter%z)

    R_max_sq = R_max**2

    out = .false.
    do i = 1, nstar
       a_star = get_a_star(this,i)
       
       ! position of the star
       x = a_star%position%x
       y = a_star%position%y
       z = a_star%position%z
       
       ! The square of the distance between the star and the cell center.
       R_sq = (x-xc)**2 + (y-yc)**2 + (z-zc)**2       

       if (R_sq < R_max_sq) then 
          ! the two sphere intersect
          out = .true.
          exit  ! exit this loop.
       end if          
    end do
     
  end function cell_too_close_to_star





  !
  ! Find the average density of the subcell by randomly picking points 
  ! in the subcell. Density is eveualated according to "disc_density" function 
  ! in this module, but it does not include the density of from the sph particles.
  !
  function average_disc_density(sphData, node, subcell, stellar_cluster, &
       scale_length) RESULT(out)
    implicit none
    real(double) :: out     ! out put density in [g/cm^3]
    type(sph_data), intent(in)    :: sphData
    type(octal), intent(inout)    :: node
    integer, intent(in)           :: subcell ! index of the subcell
    type(cluster), intent(in)     :: stellar_cluster
    real(double), optional, intent(inout) :: scale_length  ! [10^10cm]

    !
    integer, parameter :: nsample =300
    real(double) :: d, x, y, z, xc, yc, zc, rho_disc_ave, r, dummy_d, rho_sample
    type(octalVector)     :: cellCenter
    logical, save :: first_time = .true.
    integer :: i, n
    real(double) , parameter:: rho_min = 1.0d-23
    
    ! quick check for the first time
    if (first_time) then
       if (stellar_disc_exists(sphData)) then ! OK
          first_time = .false.
       else
          write(*,*) "Error:: You did not read in the stellar disc data! &
               &[cluster::average_disc_density]."
          stop
       end if
    end if
    
    !
    rho_disc_ave =1.0d-30
    if (stellar_disc_exists(sphData)) then
       d = node%subcellsize/2.0d0
       cellCenter = subcellCentre(node, subcell)
       xc=dble(cellCenter%x); yc=dble(cellCenter%y); zc=dble(cellCenter%z)
       
       n = 0       
       do i = 1, nsample
          ! picking a random position in the subcell
          call random_number(r)
          x =  xc - d + r*2.0d0*d
          call random_number(r)
          y =  yc - d + r*2.0d0*d
          call random_number(r)
          z =  zc - d + r*2.0d0*d               
          ! evaluate the disc density density
          ! Taking the max value samples
          rho_sample = disc_density(x, y, z, sphData, stellar_cluster, dummy_d)
          if (rho_sample > rho_disc_ave) rho_disc_ave  = rho_sample

       end do
       

!       n = 0
!       do i = 1, nsample
!          ! picking a random position in the subcell
!          call random_number(r)
!          x =  xc - d + r*2.0d0*d
!          call random_number(r)
!          y =  yc - d + r*2.0d0*d
!          call random_number(r)
!          z =  zc - d + r*2.0d0*d               
!          ! evaluate the disc density density
!          ! exclude the very low density area!
!          rho_sample = disc_density(x, y, z, sphData, stellar_cluster, dummy_d)
!          if (rho_sample > rho_min) then 
!             n = n + 1
!             rho_disc_ave = rho_disc_ave + rho_sample
!          end if
!       end do
!       
!       if (n > 0) then
!          rho_disc_ave = rho_disc_ave/dble(n)
!       else
!          rho_disc_ave = rho_min
!       end if

    end if
    
    ! Finding the scale_length
    if (present(scale_length)) then
       dummy_d = disc_density(xc, yc, zc, sphData, stellar_cluster, scale_length)
    end if


    out = rho_disc_ave  ! this should be in [g/cm^3] as in disc_density.

  end function average_disc_density




  !
  ! Find the average density of the subcell using the 8 corners points and 
  ! the center point of the subcell.
  ! Density is eveualated according to "disc_density" function 
  ! in this module, but it does not include the density of from the sph particles.
  !
  function average_disc_density_fast(sphData, node, subcell, stellar_cluster, &
       scale_length) RESULT(out)
    implicit none
    real(double) :: out     ! out put density in [g/cm^3]
    type(sph_data), intent(in)    :: sphData
    type(octal), intent(inout)    :: node
    integer, intent(in)           :: subcell ! index of the subcell
    type(cluster), intent(in)     :: stellar_cluster
    real(double),  intent(inout) :: scale_length

    !
    real(double) :: d, d2, xc, yc, zc, rho_disc_ave, dummy_d
    type(octalVector)     :: cellCenter
    logical, save :: first_time = .true.
    real(double) :: c0, c1, c2, c3, c4, c5, c6, c7, c8
    real(double) :: c9, c10, c11, c12, c13, c14, c15, c16
    real(double) :: c17, c18, c19, c20, c21, c22
    
    ! quick check for the first time
    if (first_time) then
       if (stellar_disc_exists(sphData)) then ! OK
          first_time = .false.
       else
          write(*,*) "Error:: You did not read in the stellar disc data! &
               &[cluster::daverage_disc_density]."
          stop
       end if
    end if
    
    !
    rho_disc_ave =0.0d0
    if (stellar_disc_exists(sphData)) then
       d = node%subcellsize/3.0d0
       d2 = d*2.0d0
       cellCenter = subcellCentre(node, subcell)
       xc=dble(cellCenter%x); yc=dble(cellCenter%y); zc=dble(cellCenter%z)
       
       c0 =  disc_density(xc, yc, zc, sphData, stellar_cluster, scale_length)

       c1 =  disc_density(xc+d, yc+d, zc+d, sphData, stellar_cluster, dummy_d)
       c2 =  disc_density(xc+d, yc+d, zc-d, sphData, stellar_cluster, dummy_d)
       c3 =  disc_density(xc+d, yc-d, zc-d, sphData, stellar_cluster, dummy_d)
       c4 =  disc_density(xc+d, yc-d, zc+d, sphData, stellar_cluster, dummy_d)
       c5 =  disc_density(xc-d, yc+d, zc+d, sphData, stellar_cluster, dummy_d)
       c6 =  disc_density(xc-d, yc+d, zc-d, sphData, stellar_cluster, dummy_d)
       c7 =  disc_density(xc-d, yc-d, zc-d, sphData, stellar_cluster, dummy_d)
       c8 =  disc_density(xc-d, yc-d, zc+d, sphData, stellar_cluster, dummy_d)

       c9  =  disc_density(xc+d2, yc+d2, zc+d2, sphData, stellar_cluster, dummy_d)
       c10 =  disc_density(xc+d2, yc+d2, zc-d2, sphData, stellar_cluster, dummy_d)
       c11 =  disc_density(xc+d2, yc-d2, zc-d2, sphData, stellar_cluster, dummy_d)
       c12 =  disc_density(xc+d2, yc-d2, zc+d2, sphData, stellar_cluster, dummy_d)
       c13 =  disc_density(xc-d2, yc+d2, zc+d2, sphData, stellar_cluster, dummy_d)
       c14 =  disc_density(xc-d2, yc+d2, zc-d2, sphData, stellar_cluster, dummy_d)
       c15 =  disc_density(xc-d2, yc-d2, zc-d2, sphData, stellar_cluster, dummy_d)
       c16 =  disc_density(xc-d2, yc-d2, zc+d2, sphData, stellar_cluster, dummy_d)

       c17 =  disc_density(xc+d2, yc, zc, sphData, stellar_cluster, dummy_d)
       c18 =  disc_density(xc-d2, yc, zc, sphData, stellar_cluster, dummy_d)
       c19 =  disc_density(xc, yc+d2, zc, sphData, stellar_cluster, dummy_d)
       c20 =  disc_density(xc, yc-d2, zc, sphData, stellar_cluster, dummy_d)
       c21 =  disc_density(xc, yc, zc+d2, sphData, stellar_cluster, dummy_d)
       c22 =  disc_density(xc, yc, zc-d2, sphData, stellar_cluster, dummy_d)

       
       rho_disc_ave = rho_disc_ave + max(c0, c1, c2, c3, c4, c5, c6, c7, c8, &
            c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, c21, c22)

!       rho_disc_ave = rho_disc_ave + (c0+c1+c2+c3+c4+c5+c6+c7+c8 &
!            +c9+c10+c11+c12+c13+c14+c15+c16+c17+c18+c19+c20+c21+c22)/23.0d0
       
    end if

    out = rho_disc_ave  ! this should be in [g/cm^3] as in disc_density.

  end function average_disc_density_fast



  !
  ! Find the average density of the subcell by randomly picking points 
  ! in the subcell. Density is eveualated according to "disc_density" function 
  ! in this module, but it does not include the density of from the sph particles.
  !
  function max_disc_density_from_sample(sphData, node, subcell, stellar_cluster, &
       scale_length) RESULT(out)
    implicit none
    real(double) :: out     ! out put density in [g/cm^3]
    type(sph_data), intent(in)    :: sphData
    type(octal), intent(inout)    :: node
    integer, intent(in)           :: subcell ! index of the subcell
    type(cluster), intent(in)     :: stellar_cluster
    real(double), optional, intent(inout) :: scale_length  ! [10^10cm]

    !
    integer, parameter :: nsample =300
    real(double) :: d, x, y, z, xc, yc, zc, rho_disc_max, r, dummy_d
    real(double) :: rho_tmp
    type(octalVector)     :: cellCenter
    logical, save :: first_time = .true.
    integer :: i
    
    ! quick check for the first time
    if (first_time) then
       if (stellar_disc_exists(sphData)) then ! OK
          first_time = .false.
       else
          write(*,*) "Error:: You did not read in the stellar disc data! &
               &[cluster::average_disc_density]."
          stop
       end if
    end if
    
    !
    rho_disc_max =0.0d0
    if (stellar_disc_exists(sphData)) then
       d = node%subcellsize/2.0d0
       cellCenter = subcellCentre(node, subcell)
       xc=dble(cellCenter%x); yc=dble(cellCenter%y); zc=dble(cellCenter%z)
       
       do i = 1, nsample
          ! picking a random position in the subcell
          call random_number(r)
          x =  xc - d + r*2.0d0*d
          call random_number(r)
          y =  yc - d + r*2.0d0*d
          call random_number(r)
          z =  zc - d + r*2.0d0*d
          ! evaluate the disc density density
          rho_tmp = disc_density(x, y, z, sphData, stellar_cluster, dummy_d)
          rho_disc_max = MAX(rho_disc_max, rho_tmp)
       end do
       
    end if
    
    ! Finding the scale_length
    if (present(scale_length)) then
       dummy_d = disc_density(xc, yc, zc, sphData, stellar_cluster, scale_length)
    end if


    out = rho_disc_max  ! this should be in [g/cm^3] as in disc_density.

  end function max_disc_density_from_sample


  !
  ! Find the scale size of a disc for given sizes of 
  ! minimum and maximum size scales allowed (ds_min and ds_max), 
  ! and minumum and maximum distance (e.g. radius or height above the disk..), i.e.
  ! s_min and s_max.
  function compute_scale_size(s, s_min, s_max, ds_min, ds_max, logscale) RESULT(out)
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

  end function compute_scale_size


  ! 
  ! Just returns if this cluster includes accreation disc or not
  function include_disc(this) RESULT(out)
    implicit none 
    logical :: out
    type(cluster), intent(in) :: this ! cluster

    out = this%disc_on
  end function include_disc

       



  !
  ! Given a octal and subcell index and  a grid, and a wavelength in 
  ! the local frame, this routine computes
  ! the emission bias for contiuum and line photons.
  !
  subroutine compute_emission_bias(thisOctal,subcell, grid, lambda, lambda_array, n_lambda)
    IMPLICIT NONE
    
    TYPE(octal), intent(inout) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(INOUT) :: grid
    real, intent(in) :: lambda ! local wavelength in Angstrom
    integer, intent(in) :: n_lambda
    real, intent(in) :: lambda_array(n_lambda)

    !
    real(double) :: d    ! size of the octal
    real(double) :: tau  ! optical depth across the cell
    integer :: ilambda
    real(double) :: kappaScaReal, KappaAbsReal


    d = thisOctal%subcellSize

    call hunt(lambda_array, n_lambda, lambda, ilambda)

    ! abs and scattering coeff at this wavelength
    if (.not.grid%oneKappa) then
       kappaScaReal = thisOctal%kappaSca(subcell,iLambda)
       kappaAbsReal = thisOctal%kappaAbs(subcell,iLambda)
    else
       kappaScaReal = grid%oneKappaSca(thisOctal%dustType(subcell),iLambda)*thisOctal%rho(subcell)
       kappaAbsReal = grid%oneKappaAbs(thisOctal%dustType(subcell),iLambda)*thisOctal%rho(subcell)
    endif


    if (thisOctal%inFlow(subcell)) then
       tau = d * dble(kappaAbsReal + kappaScaReal)
       if (tau>100.0d0)  tau=100.d0
       if (tau<0.001d0)  tau=0.001d0
    else
       tau = 1.0d100
    end if 

    if (tau <= 0.0d0) then
       write(*,*) "Error:: tau <=0 in [cluster_class::compute_emission_bias]."
       stop
    end if
    
    thisOctal%biasCont3D(subcell) = 1.0d0/tau
    thisOctal%biasLine3D(subcell) = 1.0d0  ! for now this is just 1
       
  end subroutine compute_emission_bias


  !
  ! Recursively assigns the emission bias to all octals.
  !
  recursive subroutine assign_emission_bias(thisOctal,grid, lambda, lambda_array, n_lambda)
    IMPLICIT NONE
    
    TYPE(octal), pointer :: thisOctal
    TYPE(gridtype), INTENT(INOUT) :: grid
    real, intent(in) :: lambda ! local wavelength in Angstrom
    integer, intent(in) :: n_lambda
    real, intent(in) :: lambda_array(n_lambda)
    
    type(octal), pointer   :: child
    integer :: subcell, ichild
     
    do subcell = 1, 8, 1   
       call compute_emission_bias(thisOctal,subcell, grid, lambda, lambda_array, n_lambda)
    end do

   
    do ichild = 1, thisoctal%nchildren, 1
      child => thisoctal%child(ichild)
      call assign_emission_bias(child, grid, lambda, lambda_array, n_lambda)
    end do


  end subroutine assign_emission_bias




  !
  ! Recursively all cells which are outside of a cydinder with radius Rmax and its center 
  ! passing through a source (star) with index i and an observer.  
  !
  ! The sources outside of the cylinder is also removed by setting their luminosity to zero.
  !
  ! NB: If the index is zero, it does not do anyting (just exits the program) with the sources
  !     and the grid/octal untouched. 
  recursive subroutine restrict(in_octal, i, nsource, sources, n_obs, Rmax)
    implicit none
    
    type(octal), pointer :: in_octal                    ! this root node of the octal in grid
    integer, intent(in)  :: i                           ! index in the source arrays (stars) 
    integer, intent(in)  :: nsource                     ! The number of sources.
    type(sourcetype), intent(inout) :: sources(nsource) ! sources (stars)
   type(octalvector), intent(in)    :: n_obs            ! direcion to the observer (should be normalized)
    real(double), intent(in)    :: Rmax        ! The radius of the cylinder (in 10^10cm)
    !
    logical, save :: first_time=.true.
    integer :: j, subcell

    type(octalvector) :: r0, r, rp
    real(double) :: p
    type(octal), pointer  :: child 

    !
    if (i==0)  then 
       return
    end if

    if (first_time) then
       write(*,*) " "
       write(*,*) "Restrcting the SED calculations to ", i, "-th star in cluster"
       write(*,*) " "
       ! assign the luminosity of the sources to be zero
       ! except for the i-th star.
       do j = 1, nsource
          if (j /= i ) sources(j)%luminosity = 1.0d-30
       end do
       !
       ! To avoid a problem in randomsource routine ..
       if (i==1) then
          ! move the star to the second position...
          sources(2) = sources(1)
          sources(1)%luminosity = 0.0d0
       end if

       first_time=.false.       
    end if

    !    
    do subcell = 1, 8
       if (in_octal%hasChild(subcell)) then
          ! find the child
          do j = 1, in_octal%nChildren
             if (in_octal%indexChild(j) == subcell) then
                child => in_octal%child(j)
                call restrict(child, i, nsource, sources, n_obs, Rmax)
                exit
             end if
          end do
       else
          ! this must be a leaf node, so check if it is 
          ! in the cylinder. 
          ! If so,  set the opacity 
          ! and the emissivity to zero. Set the inflow flag to 
          ! .false. as well. 
          r0 = sources(i)%position
          r = subcellCentre(in_octal,subcell)
          rp = r-r0

          p = ABS(rp .dot. n_obs) 

          if ( p > Rmax) then 
! UNCOMMENT THE NEXT TWO LINES FOR UNREDDEN SIMULATION
!             in_octal%inFlow(subcell) = .false.  
!             in_octal%rho(subcell) = 1.e-30
!
             in_octal%temperature(subcell) = 3.0
             in_octal%velocity = VECTOR(0.,0.,0.)
             in_octal%biasCont3D(subcell) = 1.0
             in_octal%biasLine3D(subcell) = 1.0
             in_octal%etaLine(subcell) = 1.e-30
             in_octal%etaCont(subcell) = 1.e-30
          else
             continue
          end if

!          !  For special case...
!          ! Turning off the dust emission beyond Rmax
!          if (modulus(rp) > Rmax) then
!             in_octal%temperature(subcell) = 3.0
!             in_octal%velocity = VECTOR(0.,0.,0.)
!             in_octal%biasCont3D(subcell) = 1.0
!             in_octal%biasLine3D(subcell) = 1.0
!             in_octal%etaLine(subcell) = 1.e-30
!             in_octal%etaCont(subcell) = 1.e-30
!          end if

       end if 

    end do 

  end subroutine restrict




  !
  ! Recursively reassign the artificial 10 K temperature assigned in lucy's routine
  recursive subroutine reassign_10K_temperature(in_octal)
    implicit none
    
    type(octal), pointer :: in_octal                    ! this root node of the octal in grid
    !
    integer :: j, subcell

    type(octal), pointer  :: child 

    !    
    do subcell = 1, 8
       if (in_octal%hasChild(subcell)) then
          ! find the child
          do j = 1, in_octal%nChildren
             if (in_octal%indexChild(j) == subcell) then
                child => in_octal%child(j)
                call reassign_10K_temperature(child)
                exit
             end if
          end do
       else
          ! this must be a leaf node,
          ! Correcting 10 K dust temperaure. (This is set in
          ! Lucy's routine, but it is contributing to the disc emission

          if (in_octal%temperature(subcell) == 10.0) in_octal%temperature(subcell)=3.0
          
       end if 

    end do

  end subroutine reassign_10K_temperature




end module cluster_class
