module formal_solutions

  ! ----------------------------------------------------------------------------------
  ! This is not a class difenition, but a collection of formal
  ! solution solvers.
  ! 
  !----------------------------------------------------------------------------------
  ! created  : aug-08-2003  (R. Kurosawa)
  ! Modefied : 14-may-2005  Added line flux calculation routines (R. Kurosawa)
  !----------------------------------------------------------------------------------


  use constants_mod
  use messages_mod
  use vector_mod             ! vector maths 
  use gridtype_mod, only: GRIDTYPE ! opacity grid
  use disc_class, only: alpha_disc, in_alpha_disc
  use utils_mod, only: locate, linearresample, linearresample_dble

  
  implicit none

  public :: &
       formal_sol_dust_AMR, &
       compute_obs_line_flux
  
  private :: &
       integrate_formal, &
       setup_grid, &
       create_obs_flux_map, &
       find_position_displacement,&
       in_this_area, &
       refine_ray_grid_by_tau, &
       refine_ray_constant_vel
  



contains

  !-----------------------------------------------------------------------------------------
  ! For a given initial intensity (I0) along a ray, this function will returns a intensity 
  ! value at an optical depth (tau) assuming that the source function is known. I.e., we solve 
  ! 
  !              dI
  !             ---- = S - I 
  !             dTau
  !
  !                         -(tau1-tau0)                        -(tau1-t)
  !    I(tau1) = I(tau0) * e               +   Integrate[(S(t)*e        , {t, tau0, tau1}]
  !
  !  
  !    0  I0                          1  I1
  !    *--> ------------------------- * -->
  !   
  !    tau=tau0                       tau=tau1
  !
  !  Normally, tau0=0.
  !
  ! ------------------------------------------------------------------------------------------
  !
  !  NB. This routine is based on integratePathAMR routine. 
  !       For now it only works for the case of continuum rad transfer of dust emission. 
  !   
  !   ===> To avoid the underflow, the caluculations and the results are offset by "offset"!
  !
  function formal_sol_dust_AMR(I0, tau0, wavelength, aVec, uHat, Grid, &
       contPhoton, offset, tau_inf)  RESULT (I1)
    use path_integral, only: intersectCubeAMR
    use amr_mod, only: inOctal, amrGridValues
    use octal_mod, only: OCTAL 

    implicit none

    real(double) ::  I1    ! output intensity.
    !
    real(double), intent(in)  ::  I0    ! input intensity
    real(double), intent(in)  ::  tau0  ! optical depth at the initial point 

    real(double), intent(in)          :: wavelength        ! the wavelength [A]
    type(VECTOR), intent(in)  :: aVec          ! starting position vector
    type(VECTOR), intent(in)  :: uHat          ! direction
    type(GRIDTYPE), intent(in)     :: grid          ! the opacity grid
    logical, intent(in)            :: contPhoton    ! is this a continuum photon?
    real(double), intent(in) :: offset     ! offset scale factor
    real(double), optional, intent(inout) :: tau_inf ! offset scale factor
    
    
    type(VECTOR) :: octVec
    integer :: subcell
    real(double) :: tval
    integer               :: nTau                   ! size of optical depth arrays
    type(VECTOR)     :: rVec                   ! position vector
!    real :: kappaScaReal, kappaAbsReal
    real(double) :: kappaSca, kappaAbs
    real(double) :: etaCont
    real(oct) :: chi, eta, tau
    real(double)  :: integral
    integer :: ilambda
    logical :: escaped
    type(OCTAL), pointer :: thisOctal
    type(OCTAL),pointer :: oldOctal
    real(double) :: dtau, delta
    real(double)::  I00                         ! input intensity

    kappaAbs = 0.d0; kappaSca = 0.d0
    I00 = I0*offset

    !
    ! Just checking ... 
    ! 
    
    if (.not. grid%adaptive) then
       write(*,*) "Error:: Non AMR option has not implemented yet"//&
            &"in [formal_solutions::formal_sol_dust_AMR]."
       write(*,*) "  Exiting program ..."
       stop
    end if
    
    if (.not. contPhoton) then
       write(*,*) "Error:: linePhoton option has not implemented yet"//&
            &" in [formal_solutions::formal_sol_dust_AMR]."
       write(*,*) "  Exiting program ..."
       stop     
    end if


    ! initialize variables
    
    ! locate this wavelength in the grid of wavelengths
    if (grid%flatspec.or.(grid%doRaman)) then
       iLambda = 1
    else
       call locate(grid%lamArray, grid%nLambda, real(wavelength), iLambda)
    endif
    
    
    !
    !
    ! First we find the optical depth to the observer (tau_inf).
    rvec=avec
    nTau = 1
    escaped = .false.
    oldOctal => grid%octreeRoot
    tau_inf = 0.0
    integral = 0.0
    
    do while (.not.escaped) 
       octVec = rVec
       if (.not.inOctal(grid%octreeRoot, octVec)) then
          escaped = .true.
       endif
       if (.not.escaped) then
          call intersectCubeAMR(grid, rVec, uHat, tVal)
          call amrGridValues(grid%octreeRoot, octVec, startOctal=oldOctal,iLambda=iLambda, &
               foundOctal=thisOctal, foundSubcell=subcell, &
               kappaAbs=kappaAbs,kappaSca=kappaSca, etaCont=etaCont, grid=grid)

          nTau = nTau + 1
          
          eta = etaCont  ! [erg/s/sr/cm/cm^3]
          chi = ( kappaSca + kappaAbs )  ! [1/10^10cm]

          tau_inf = tau_inf + tval*chi  ! should be dimensionless
                    
          rVec = rVec + tval * uHat
          oldOctal => thisOctal
          
       endif
    end do


    !
    !  Now performs the integration
    !
    
    rvec=avec
    nTau = 1
    escaped = .false.
    oldOctal => grid%octreeRoot
    tau = 0.0
    integral = 0.0
    
    do while (.not.escaped) 
       octVec = rVec
       if (.not.inOctal(grid%octreeRoot, octVec)) then
          escaped = .true.
       endif
       if (.not.escaped) then
          call intersectCubeAMR(grid, rVec, uHat, tVal)
          call amrGridValues(grid%octreeRoot, octVec, startOctal=oldOctal,iLambda=iLambda, &
               foundOctal=thisOctal, foundSubcell=subcell, &
               kappaAbs=kappaAbs,kappaSca=kappaSca, etaCont=etaCont, grid=grid)

          nTau = nTau + 1
          
!          eta = 0.0
          eta = etaCont  ! [erg/s/sr/cm/cm^3]
          chi = (kappaSca + kappaAbs)  !in [1/10^cm] I hope

          dtau = tval*chi
          tau = tau + dtau  ! now should be in dimensionless units
          
          if (chi > 0.0) then 
             delta = (tau_inf - tau)
             integral = integral + eta/(chi*1.0d-10) * EXP(-delta)*dtau * offset  
             ! ----- [erg/s/sr/cm^2/cm * offset]
          else
             write(*,*) "Error:: chi <= 0 in [formal_solutions::formal_sol_dust_AMR]."
             write(*,*) "  Exiting program ..."  
          end if
          
          rVec = rVec + tval * uHat
          oldOctal => thisOctal
          
       endif
    end do
    
    integral = integral*1.e-8         ! per cm to converting to per Angstrome
    delta = tau_inf-tau0
    I1 = I00*EXP(-delta) + integral   ! should be in [erg/s/sr/cm^2/A * offset]

    if (I1 < 1.0d-250) I1 = 1.0d-200 

    
  end function formal_sol_dust_AMR





  !
  ! Computes the observed flux from a star for a given direction of observer
  ! and the distance to observer as a function of inclination.
  ! Flux units in  [erg cm^-2 s^-2 cm^-1]

  subroutine compute_obs_line_flux(lambda0, mass_ion, R_star, R_max_acc, R_max_wind, surface, &
       nr_wind, nr_acc, nr_core, nphi, dir_obs, dist_obs, grid, sampleFreq, opaqueCore,  &
       flux, lambda, nlam, filename, thin_disc_on, ttau_disc_on, ttau_jet_on, ttau_discwind_on, npix, &
       do_alphadisc_check, alphadisc_par, thinLine, emissOff, pos_disp) 

    use surface_mod, only: SURFACETYPE, whichElement, isHot
    use utils_mod, only: blackBody
    use messages_mod, only : myRankIsZero
    use timing, only: tune
    use parallel_mod

    implicit none
#ifdef MPI
    include 'mpif.h'        
#endif
    real,         intent(in)  :: lambda0       ! [A] Wavelength of the line center
    real,         intent(in)  :: mass_ion      ! mass of atom in [g]
    real(double), intent(in)  :: R_star        ! [10^10 cm] radius of a star
    real(double), intent(in)  :: R_max_acc     ! [10^10 cm] the max radius of accretion flow for flux integration
    real(double), intent(in)  :: R_max_wind    ! [10^10 cm] the max radius for flux integration
    integer, intent(in)                :: nr_wind       ! number of radial points for integration (wind part)
    integer, intent(in)                :: nr_acc        ! number of radial points for integration (accretion)
    integer, intent(in)                :: nr_core       ! number of radial points for integration (core rays)
    integer, intent(in)                :: nphi          ! number of angle  points for integration 

    type(SURFACETYPE),intent(in)       :: surface       ! surface elemet of the core
    real(double), intent(in)           :: dist_obs      ! [cm] distance to the observer
    type(VECTOR), intent(in)      :: dir_obs       ! direction to the observer
    type(GRIDTYPE), intent(in)         :: grid          ! the opacity grid
    real, intent(in)                   :: sampleFreq    ! max. samples per grid cell
    logical, intent(in)                :: opaqueCore    ! T if the core is opaque
    integer, intent(in)                :: nlam          ! number of wavelength/freqeuncy points
    real,    intent(in)                :: lambda(:)  ! wavelength array [A]
    real(double), intent(inout)        :: flux(:)    ! Flux array dimension:
    character(LEN=*), intent(in)       :: filename      ! name of file for the output flux
    logical, intent(in)                :: thin_disc_on  ! T to include thin disc
    logical, intent(in)                :: ttau_disc_on  !
    logical, intent(in)                :: ttau_jet_on   !
    logical, intent(in)                :: ttau_discwind_on   !
    integer, intent(in)                :: npix          ! image size = 2*npix +1 should be an even number
    ! The followings are needed for resampling integration rays
    logical, intent(in)                :: do_alphadisc_check  ! if T the disc check will be done
    type(alpha_disc), intent(in)       :: alphadisc_par       ! parameters for alpha disc
    logical, intent(in)                :: thinLine            ! if T surpress, the abs of cont photon by line
    logical, intent(in)                :: emissOff            ! if T surpress, the line emission
    logical, intent(in)                :: pos_disp            ! if T do position diplacement calulation


    !
    !
    !
    real(double)                       :: dphi            !increment of the angle points
    type(VECTOR), allocatable     :: p_core(:)       ! integration grids
    type(VECTOR), allocatable     :: p_core_orig(:)  ! integration grids
    real(double), allocatable          :: dA_core(:)      ! surface elements at p
    real(double), allocatable          :: dr_core(:)      ! increment for radial grid
    type(VECTOR), allocatable     :: p_acc(:)        ! integration grids
    type(VECTOR), allocatable     :: p_acc_orig(:)   ! integration grids
    real(double), allocatable          :: dA_acc(:)       ! surface elements at p
    real(double), allocatable          :: dr_acc(:)       ! increment for radial grid
    type(VECTOR), allocatable     :: p_wind(:)       ! integration grids
    type(VECTOR), allocatable     :: p_wind_orig(:)  ! integration grids
    real(double), allocatable          :: dA_wind(:)      ! surface elements at p
    real(double), allocatable          :: dr_wind(:)      ! increment for radial grid
    !----------------------------------------------------------------------------
    type(VECTOR), allocatable     :: p_orig(:)       ! integration grids
    real(double), allocatable          :: dA(:)           ! surface elements at p
    real(double), allocatable          :: dr(:)           ! increment for radial grid
    real(double)                       :: R_max           ! max radius for the plot.
    !------------------------------------------------------------------------------
    real(double)                       :: I0  ! [erg cm^-2 s^-2 A^-1 sr^-1] input intensity
    real(double)  :: F_acc, F_core, F_wind, Fi, I1, tau0, freqc
    !
    integer :: i, j, k, nc, na, nw, nray, j_beg, j_end
    logical :: hitcore                ! has the photon hit the core    
    integer :: error                  ! error code returned
    real(double) :: IC_hot, IC_normal, taccretion !, fillfactor, TTauriMdotSurface
    real :: lambda_j
    logical :: continuum
    real(double), allocatable :: flux_map(:,:)      ! Flux map  array dimension= nlam x #of integration points
    real(double), allocatable :: integrated_flux(:) ! freq. integrate Flux map  array (dim=nray)
    character(LEN=80) :: filename_map
    real(double) :: mu  ! dot product between the normal of the stellar surface and the direction of observer
    logical :: logscale
    integer :: i_beg, i_end, ielement
    logical :: coreray, magray
    real(double) :: halfboxsize  ! [10^10cm] a half the size of simulation box
    ! the following should be put in input parameter later
    logical :: limb_darkening 
!    logical, parameter :: limb_darkening = .true.  
    real(double), parameter :: a1= 0.93d0       ! coefficiets in the linear limb dakening law 
    real(double), parameter :: a2= -0.23d0      ! coefficiets in the linear limb dakening law 
    real(double), parameter :: a0= 1.0d0-a1-a2  ! coefficiets in the linear limb dakening law 

#ifdef MPI
    ! For MPI implementations
    integer             ::   my_rank            ! my processor rank
    integer             ::   ncpu               ! The number of processes
    integer             ::   ierr               ! error flag
    integer,allocatable ::   BelongRank(:)  
    integer             ::   tag = 0
    logical             ::   rankComplete
    real(double), allocatable :: buffer1(:), buffer2(:)
    real(double)        :: buffer3(3), buffer4(3)
#endif
    limb_darkening = .false.  
    halfboxsize = grid%octreeRoot%subcellSize 

    ! allocating the grid points
    nc = nr_core*nphi
    ALLOCATE(p_core(nc), p_core_orig(nc), dA_core(nc), dr_core(nc))
    na = nr_acc*nphi
    ALLOCATE(p_acc(na), p_acc_orig(na), dA_acc(na), dr_acc(na))
    nray=nc+na   ! total number of rays


    ! continuum frequency (take the first element of the lamda array)
    freqc = cSpeed_dbl/(lambda(1)*1.d-8)  ! Hz
    call locate(surface%nuArray,SIZE(surface%nuArray),REAL(freqc), k)

!     ! estimating the temperature of the hot ring
!     fillfactor = SUM(surface%element%area,MASK=surface%element%hot) / &
!          SUM(surface%element%area)
!     tAccretion = (real(SUM(surface%totalAccretion),kind=db) / &
!          ((fourPi*(surface%radius*1.e10)**2.*stefanBoltz)*fillFactor) )**0.25

!     !================CHECK UNITS HERE!! ===========================
!     IC_hot = blackbody(REAL(tAccretion), 1.e8*REAL(cSpeed)/surface%nuArray(k)) ! [B_nu]


    if (limb_darkening) then
       ! For now we assume no lim darkening law; 
       !           I(mu) = I0*mu for mu>0. 
       ! Then using the definition of Fnu, one finds:
       !           I(mu) = (3/2pi)*mu*Fnu
       !------------------------------------------------
       ! angular factor will be multiplied later.
       ! BE CAREFULE nNu here is not H_nu, but rather F_nu!
       IC_normal = (1.5d0/Pi)*surface%hNuArray(k)  ! combination of this and angular factor works

       ! THIS SHOULD BE CHANGED ACCORDING TO THE VALUES OF THE LIMB-DARKENING
       ! COEFICIENTS LATER.
       
    else
       ! IF YOU WANT TO BE CONSISTENT WITH MONTE CALRO ROUTINE, THIS SHOULD
       !  BE YOUR CHOICE INSTEAD OF THE ABOVE.
       ! without anguar factor (isotropic for mu>0), then I_nu= F_nu/Pi
       IC_normal = surface%hNuArray(k)/Pi   ! Works without the angular factor
    end if
    
    !
    ! setting up the grid for integration
    !
    ! -- For core rays
    logscale = .false.
    coreray = .true.
    magray = .false.
    call setup_grid(p_core, dA_core, dr_core, dphi, nc, nr_core, nphi, &
         1.0d-5*R_star,  R_star, dir_obs, logscale, coreray, &
         magray, halfboxsize ,p_core_orig)    

    ! -- For rays crossing magnetosphere
    logscale = .false.
    coreray = .false.
    magray = .true.
    if (ttau_disc_on .or. ttau_jet_on .or. ttau_discwind_on) then       
       call setup_grid(p_acc, dA_acc, dr_acc, dphi, na, nr_acc, nphi, &
            R_star,  R_max_acc,  dir_obs, logscale, coreray, &
            magray, halfboxsize,  p_acc_orig) 
    else 
       ! use the small box size ..
       call setup_grid(p_acc, dA_acc, dr_acc, dphi, na, nr_acc, nphi, &
            R_star,  R_max_acc,  dir_obs, logscale, coreray, &
            magray, R_max_acc,  p_acc_orig) 
    end if
    
    ! -- For rays crossing wind/jet/disc
    if (ttau_disc_on .or. ttau_jet_on .or. ttau_discwind_on) then       
       nw = nr_wind*nphi
       ALLOCATE(p_wind(nw), p_wind_orig(nw), dA_wind(nw), dr_wind(nw))
       logscale = .true.
       coreray = .false.
       magray = .false.
       call setup_grid(p_wind, dA_wind, dr_wind, dphi, nw, nr_wind, nphi, &
            R_max_acc, R_max_wind, dir_obs, logscale, coreray,  &
            magray, halfboxsize, p_wind_orig) 
       nray=nc+na+nw ! total number of rays
    end if
       
    ! Allocating memory for flux_map array
    ALLOCATE(flux_map(nlam, nray))

    !
    ! start integration    
    !    
    j_beg=1; j_end=nlam   ! index range of freq. loop for single precessor job
    flux(:) = 0.0d0
    flux_map(:,:) = 0.0d0
    F_core=0.0d0; F_acc=0.0d0; F_wind=0.0d0  ! ONLY rank =z has correct flux!

#ifdef MPI
    ! FOR MPI IMPLEMENTATION=======================================================
    !  Get my process rank # 
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  
    ! Find the total # of precessor being used in this run
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpu, ierr)
    
 
    print *, ' '
    print *, 'Compute_obs_line_flux routine computed by ', ncpu-1, ' processors.'
    print *, ' '
    print *, '         ===> RAY LOOP parallelized <=== '
    print *, ' '

    ! ============================================================================
#endif

    ! Note: The first wavelength bin is used for contiuum calculation!!
    do j = j_beg, j_end
       if (doTuning) call tune(6, "One wavelength  Loop") ! start a stopwatch  
       if (j==1) then
          continuum = .true.
          lambda_j = lambda0  ! use the line center wavelength for contiuum calculation...
       else 
          continuum=.false.  ! line
          lambda_j = lambda(j)
       end if
       !
       ! For core rays
       !
       F_core = 0.0d0

#ifdef MPI
    ! FOR MPI IMPLEMENTATION=======================================================
    if (my_rank == 0) then
       ! we will use an array to store the rank of the process
       !   which will calculate each octal's variables
       allocate(BelongRank(nc))
       call mpiBlockHandout(ncpu,BelongRank,blockDivFactor=1,tag=tag,&
                            maxBlockSize=1,setDebug=.false.)
       deallocate(BelongRank)
    endif
    ! ============================================================================
#endif
       i_beg=1; i_end=nc

#ifdef MPI
 if (my_rank /= 0) then
  block01: do     
  call mpiGetBlock(my_rank,i_beg,i_end,rankComplete,tag,setDebug=.false.)
   if (rankComplete) exit block01 
#endif  

       do i = i_beg, i_end
          ! first check if the beginning of the ray is at hot surface.
          ! -- using the routines in surface_mod.f90
          !
          ielement = whichElement(surface, p_core(i))           

!          !need to get intensity from the surface element of the star.
!          TTauriMdotSurface = TTauriVariableMdot(p_core(i),grid) ! g s^-1
!          if (TTauriMdotSurface > 1.e15) then ! in hot ring     
          if (isHot(surface,ielement)) then
             tAccretion = surface%element(ielement)%temperature
             !================CHECK UNITS HERE!! ===========================
             IC_hot = blackbody(REAL(tAccretion), 1.e8*REAL(cSpeed)/surface%nuArray(k)) ! [B_nu]

             I0  = IC_hot + IC_normal
             !                     [erg cm^-2 s^-1 Hz^-1 sr^-1] input intensity
          else              
             I0  = IC_normal
             !                     [erg cm^-2 s^-1 Hz^-1 sr^-1] input intensity
          end if

          if (limb_darkening) then
             ! Adjust this with angular factor
             mu = ABS(p_core(i) .dot. dir_obs)/MODULUS(p_core(i))
             ! I0 = I0*mu ! initial trial...
             I0 = I0*(a0 + a1*mu + a2*mu*mu)  ! initial trial...
          end if

          tau0 = 0.0d0
          call integrate_formal(I0, I1,  lambda_j,  lambda0,  mass_ion, &
               p_core(i), dir_obs,  continuum, grid, &
               hitCore, thin_disc_on, opaqueCore,  sampleFreq, error, &
               do_alphadisc_check, alphadisc_par, thinLine, emissOff)
          if (error == -10 .or. hitcore) then
             Fi = 1.0d-100
          else
             Fi = I1*dA_core(i)*1.0d20  ! [erg s^-1 Hz^-1 sr^-1] ??
          end if
          F_core = Fi + F_core
          flux_map(j,i) = Fi
       end do

#ifdef MPI
 end do block01        
 end if ! (my_rank /= 0)


 write(*,*) myrankGlobal, " waiting at barrier"
 call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

#endif


       !-------------------------------------------------------------
       !Now accretion flow regions
       !-------------------------------------------------------------

       F_acc = 0.0d0  
#ifdef MPI
    ! ============================================================================
    if (my_rank == 0) then
       ! we will use an array to store the rank of the process
       !   which will calculate each octal's variables
       allocate(BelongRank(na))
       call mpiBlockHandout(ncpu,BelongRank,blockDivFactor=1,tag=tag,&
                            maxBlockSize=1,setDebug=.false.)
       deallocate(BelongRank)
    endif
    ! ============================================================================
#endif       
       i_beg=1; i_end=na

#ifdef MPI
 if (my_rank /= 0) then
  block02: do     
 call mpiGetBlock(my_rank,i_beg,i_end,rankComplete,tag,setDebug=.false.)
   if (rankComplete) exit block02 
#endif
       do i = i_beg, i_end
          !need to get intensity from the surface element of the star.
          I0  = 1.0d-100  ! [erg cm^-2 s^-2 cm^-1 sr^-1] should be zero at the outer boundary
          tau0 = 0.0d0
          call integrate_formal(I0, I1,  lambda_j,  lambda0, mass_ion,  &
               p_acc(i), dir_obs, continuum, grid, &
               hitCore, thin_disc_on, opaqueCore,  sampleFreq, error, &
               do_alphadisc_check, alphadisc_par,thinLine, emissOff)       
          if (error == -10 .or. hitcore) then
             Fi = 1.0d-100
          else
             Fi = I1*dA_acc(i)*1.0d20  ! [erg s^-1 Hz^-1 sr^-1] ??
          end if
          F_acc  = Fi + F_acc
          flux_map(j,nc+i) =  Fi
       end do
#ifdef MPI
 end do block02        
 end if ! (my_rank /= 0)

   call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
#endif


       !---------------------------------------------------
       ! For wind/jet region
       !---------------------------------------------------       

       ! For outer parts (wind and disc regions)
       F_wind = 0.0d0

       if (ttau_disc_on .or. ttau_jet_on .or. ttau_discwind_on) then

#ifdef MPI
    ! ============================================================================
    if (my_rank == 0) then
       ! we will use an array to store the rank of the process
       !   which will calculate each octal's variables
       allocate(BelongRank(nw))
       call mpiBlockHandout(ncpu,BelongRank,blockDivFactor=1,tag=tag,&
                            maxBlockSize=1,setDebug=.false.)
       deallocate(BelongRank)
    endif
    ! ============================================================================
#endif

          i_beg=1; i_end=nw

#ifdef MPI
 if (my_rank /= 0) then
  block03: do     
 call mpiGetBlock(my_rank,i_beg,i_end,rankComplete,tag,setDebug=.false.)
   if (rankComplete) exit block03 
#endif
          do i = i_beg, i_end
             !need to get intensity from the surface element of the star.
             I0  = 1.0d-100  ! [erg cm^-2 s^-2 cm^-1 sr^-1] should be zero at the outer boundary
             tau0 = 0.0d0
             call integrate_formal(I0, I1,  lambda_j,  lambda0, mass_ion,  &
                  p_wind(i), dir_obs, continuum, grid, &
                  hitCore, thin_disc_on, opaqueCore,  sampleFreq, error, &
                  do_alphadisc_check, alphadisc_par,thinLine, emissOff)
             if (error == -10 .or. hitcore) then
                Fi = 1.0d-100
             else
                Fi = I1*dA_wind(i)*1.0d20  ! [erg s^-1 Hz^-1 sr^-1] ??
             end if
             F_wind  = Fi + F_wind
             flux_map(j,nc+na+i) =  Fi
          end do
#ifdef MPI
  end do block03        
  end if ! (my_rank /= 0)
 
   call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
#endif 

       end if



#ifdef MPI  
  ! UPDATING flux for j-th wavelength bin 
  if (j==1) allocate(buffer1(nray), buffer2(nray))
  buffer1(1:nray) = 0.0d0; buffer2(1:nray) = 0.0d0
  buffer1(1:nray) = flux_map(j, 1:nray)
   
  call MPI_REDUCE(buffer1, buffer2, nray ,MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  if (my_rank == 0) then
      flux_map(j, 1:nray) = buffer2(1:nray)  ! Only rank=0 has correct flux!!!
  endif ! 
  if (j==nlam) deallocate(buffer1, buffer2)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
   
   

  buffer3(:) = 0.0d0; buffer4(:) = 0.0d0
  buffer3(1) = F_core; buffer3(2) = F_acc; buffer3(3) = F_wind
   
  call MPI_REDUCE(buffer3, buffer4, 3 ,MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  if (my_rank == 0) then
      F_core = buffer4(1)  
      F_acc  = buffer4(2)  
      F_wind = buffer4(3)  
  endif ! 
  call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
#endif  


   if (myRankIsZero) then
       ! Total Flux 
       ! scaling the flux at the observer's distance.
       flux(j) = (F_core + F_acc + F_wind ) / (r_Star**2*1.d20)!*(R_star*1.0d10/dist_obs)**2  ! [erg/s cm^-2 Hz^-1]
       flux_map(j,1:nray) = flux_map(j,1:nray)*(R_star*1.0d10/dist_obs)**2         ! [erg/s cm^-2 Hz^-1]
       ! Convet it some other units if you want below

   end if
   call torus_mpi_barrier


   if (writeoutput) write(*,*) j, "/", nlam, " of flux integration done,"

       call torus_mpi_barrier

       if (doTuning) call tune(6, "One wavelength  Loop") ! stop a stopwatch  
    end do  ! wavelength loop


    !
    ! Writing the results to files
    !

 if (myRankIsZero) then
    open(unit=89, file=TRIM(ADJUSTL(filename))//".dat", status="replace")
    open(unit=90, file=TRIM(ADJUSTL(filename))//"_v.dat", status="replace")
    write(89,"(a)") "#  Written by [formal_solutions::compute_obs_line_flux]."
    write(89,"(a20, 1PE12.4, a25)") "#  Continuum Flux = ", flux(1),   " [erg/s cm^-2 Hz^-1]" 
    write(89,"(a20, 1PE12.4, a25)") "#   at wavelength = ", lambda0, " [A]" 
    write(90,"(a)") "#  Written by [formal_solutions::compute_obs_line_flux]."
    write(90,"(a20, 1PE12.4, a25)") "#  Continuum Flux = ", flux(1),   " [erg/s cm^-2 Hz^-1]" 
    write(90,"(a20, 1PE12.4, a25)") "#   at wavelength = ", lambda0, " [A]" 
    do j = 2, nlam
       write(89,*) lambda(j), flux(j)
       write(90,*) cspeed*1.0d-5*(lambda(j)-lambda0)/lambda0, flux(j)/flux(1)
    end do
    close(89)
    close(90)

    ! plot a image map here
    ! -- using a routine in this module
    if (ttau_disc_on .or. ttau_jet_on .or. ttau_discwind_on) then
       r_max = r_max_wind
       ALLOCATE(p_orig(nray), dA(nray), dr(nray))
       ! pack the ray info for plotting
       p_orig(1:nc)  =  p_core_orig(1:nc)
       dA(1:nc) = dA_core(1:nc)
       p_orig(nc+1:nc+na)  =  p_acc_orig(1:na)
       dA(nc+1:nc+na) = dA_acc(1:na)
       p_orig(nc+na+1:nray)  =  p_wind_orig(1:nw)
       dA(nc+na+1:nray) = dA_wind(1:nw)
       !
       dr(1:nc) = dr_core(1:nc)
       dr(nc+1:nc+na) = dr_acc(1:na)
       dr(nc+na+1:nray) = dr_wind(1:nw)

    else
       r_max = r_max_acc
       ALLOCATE(p_orig(nray),  dA(nray), dr(nray))
       ! pack the ray info for plotting
       p_orig(1:nc)  =  p_core_orig(1:nc)
       dA(1:nc) = dA_core(1:nc)
       p_orig(nc+1:nc+na)  =  p_acc_orig(1:na)
       dA(nc+1:nc+na) = dA_acc(1:na)
       !
       dr(1:nc) = dr_core(1:nc)
       dr(nc+1:nc+na) = dr_acc(1:na)

    end if
    
    !
    !integrate the flux over the wavelength range
    ALLOCATE(integrated_flux(nray))
    do i = 1, nray
       Fi = 0.0d0
       do j = 1, nlam
          Fi = Fi + flux_map(j,i)
       end do
       integrated_flux(i) = Fi/dA(i)  ! Check units here!
    end do


    filename_map = TRIM(ADJUSTL(filename))//"_xyz_flux.dat"
    open(unit=299, file=filename_map, status="replace")
    do i = 1, nray
       write(299,*) p_orig(i), integrated_flux(i)
    end do
    close(299)



    !
    ! Creates NDF image
    write(*,*) " "
    write(*,*) "Creating and writing flux map (and spectro-astrometry) data... "
    filename_map = "image_"//TRIM(ADJUSTL(filename))
    call create_obs_flux_map(integrated_flux, p_orig,  dr, dphi, &
       nray, R_star, R_max, TRIM(filename_map), npix)
    if (pos_disp) &
         call find_position_displacement(flux_map, p_orig, dA, dr, dphi, &
         nray, nphi, R_star, R_max,  TRIM(filename_map), npix, &
         dist_obs, lambda, lambda0)

    if (ttau_disc_on .or. ttau_jet_on .or. ttau_discwind_on) then
       ! save one for zoomed in image (1/200 of Rmax)
       filename_map = "image_"//TRIM(ADJUSTL(filename))
       call create_obs_flux_map(integrated_flux, p_orig,  dr, dphi, &
            nray,  R_star, MAX(0.5d-2*R_max, 10.d0*R_star), &
            TRIM(filename_map), npix)    
       if (pos_disp) &
            call find_position_displacement(flux_map, p_orig, dA, dr, dphi, &
            nray, nphi, R_star, MAX(1.0d-2*R_max, 10.d0*R_star), &
            TRIM(filename_map), npix, &
            dist_obs, lambda, lambda0)
!       ! save one for zppmed  in image (1/10) of Rmax
!       filename_map = TRIM(ADJUSTL(filename))//"_natural_zoom2_image001"
!       call create_obs_flux_map(integrated_flux, p_orig, dA, dr, dphi, &
!            nlam, nray, nphi, R_star, MAX(1.0d-1*R_max, 10.d0*R_star), &
!            dir_obs, TRIM(filename_map), npix)    
!       if (pos_disp) &
!            call find_position_displacement(flux_map, p_orig, dA, dr, dphi, &
!            nlam, nray, nphi, R_star, MAX(1.0d-1*R_max, 10.d0*R_star), &
!            dir_obs, TRIM(filename_map), npix, &
!            dist_obs, lambda, lambda0)
    end if

    write(*,*) " .. finished creating an writing flux map .. "

    ! clean up
    DEALLOCATE(p_orig)
    DEALLOCATE(dA)
    DEALLOCATE(dr)
    DEALLOCATE(integrated_flux)

 endif ! (my_rank == 0)

 if (myRankIsZero) then
    write(*,*) " "
    write(*,*) "Wating for a root node to finish creating flux_map and "
    write(*,*) "    spectro-astrometry data..."
 endif

 call torus_mpi_barrier
 write(*,*) "... done."

    !
    ! clean up
    DEALLOCATE(p_core)
    DEALLOCATE(p_core_orig)
    DEALLOCATE(dA_core)
    DEALLOCATE(dr_core)    
    DEALLOCATE(p_acc)
    DEALLOCATE(p_acc_orig)
    DEALLOCATE(dA_acc)
    DEALLOCATE(dr_acc)

    if (ttau_disc_on .or. ttau_jet_on .or. ttau_discwind_on) then
       DEALLOCATE(p_wind)
       DEALLOCATE(p_wind_orig)
       DEALLOCATE(dA_wind)
       DEALLOCATE(dr_wind)
    end if

    DEALLOCATE(flux_map)

  end subroutine compute_obs_line_flux



  !
  ! Routine to setup the grid used in the observed flux integration (observed_flux)
  ! For the rays which intersect
  subroutine setup_grid(p, dA, dr, dphi, nray, nr, nphi, R_min, R_max,  dir_obs, &
    logscale, coreray, magray, halfboxsize, p_orig)
    implicit none
    integer, intent(in)                :: nray          ! number of rays
    integer, intent(in)                :: nr            ! number of radial points for integration
    integer, intent(in)                :: nphi          ! number of angle  points for integration
    type(VECTOR), intent(inout)   :: p(nray)       ! integration grids
    real(double), intent(inout)        :: dA(nray)      ! [10^20 cm^2] surface elements at p
    real(double), intent(inout)        :: dr(nray)      ! [10^10 cm] radial increments
    real(double), intent(out)          :: dphi          ! [rad] angle incrment 
    real(double), intent(in)  :: R_min                  ! [10^10cm] minimum radius 
    real(double), intent(in)  :: R_max                  ! [10^10cm] the wavelength 
    
    type(VECTOR), intent(in)      :: dir_obs       ! direction
    logical, intent(in)                :: logscale      ! if T, r is spaced in log scale
    logical, intent(in)                :: coreray       ! if T,  core rays are to be created
    logical, intent(in)                :: magray        ! if T,  this is for rays in magneosphere zone
    real(double), intent(in)           :: halfboxsize   ! [10^10cm] a half of a whole simulation space.
    ! integration grids before rotated to the direction of the observer
    type(VECTOR), intent(inout)   :: p_orig(nray)  
                                                  
    
    integer :: i , j , k
    type(VECTOR)                  :: q             ! 
    real(double)              :: r, phi, drad, rp
    real(double)              :: log_r, log_R_max, log_R_min
    real(oct)              :: x, y, z
    real(double) :: mu, psi
    !
    type(VECTOR) :: refaxis, xhat, zhat, nobs
    
    nobs = dir_obs
    call normalize(nobs)  ! just in case..
    xhat = VECTOR(1.0d0, 0.0d0, 0.0d0)
    zhat = VECTOR(0.0d0, 0.0d0, 1.0d0)

!    open(unit=301, file="p.dat", status="replace")
    ! ===============================================================
    ! Initially setup the grids for an observer is assumed to be on x axis , 
    ! the later, we rotate the vectors to the direction of observers.
    log_R_min = LOG(R_min)
    log_R_max = LOG(R_max)
    dphi = 2.0d0*piDouble/dble(nphi)
    i = 0
    do k = 1, nr
       
       if (logscale) then
          log_r = (log_R_max - log_R_min)*dble(k-1)/dble(nr)  + log_R_min
          r = EXP(log_r)
          log_r = (log_R_max - log_R_min)*dble(k)/dble(nr)  + log_R_min
          rp = EXP(log_r)
          drad = rp-r  ! 10^10cm
          r = r + drad/2.0d0
       else
          drad = (R_max-R_min)/dble(nr) 
          r = R_min + drad*dble(k-1) + drad/2.0d0
       end if


       phi = 0
       do j = 1, nphi
          i = i + 1
          phi = dble(j-1)*dphi + dphi/2.0d0

          y = r*COS(phi);  z = r*SIN(phi)  ! [10^10cm]

          if (coreray) then
             ! the rays intersecting the core
             x = SQRT( max(0.0d0, R_max*R_max - (y*y +z*z))  )
          else if (magray) then
             ! the rays intersecting the magnetosphere
             x = -halfboxsize
          else 
             ! the rays intersecting the wind
             x = -R_max
          end if

          q = VECTOR(1.001*x, y, z)
!          q = VECTOR(x, y, z)
          
          p(i) = q  ! (vector addition) [10^10cm]
          dA(i) = r*dphi*drad   ! [10^20 cm^2]          
          dr(i) = drad


          p_orig(i) = p(i)  ! save the original for later use.
          
          ! now rotate the position vectors.
          ! Method: 1. Find the polar (pi/2-mu) and azimuth (psi) angle of dir_obs
          !         2. Rotate p(i) around z-axis by psi ==> p'(i)
          !         2. Find the reference axis (dir_obs x zhat)
          !         3. Rotate p'(i) around the reference axis by mu


          mu = ASIN(nobs%z)  ! assuming that dir_obs is normalized
          psi = ATAN2(nobs%y, nobs%x)

          ! now rotate around z by psi
          ! function does passive rotation hence minus
          p(i) = rotateZ(p(i), -psi)  

          ! find the new rotational axis
          refaxis = nobs .cross. zhat
          call normalize(refaxis)  ! just in case...
          
          ! Another rotation here.
          p(i) = arbitraryRotate(p(i), mu, refaxis) 


       end do
    end do
!    close(301)
          
  end subroutine setup_grid






  !  Given starting point (aVec) and the direction (uHat), this routine computes: 
  !  the optical depth, emissivity, line segments along the line.
  !  ==> Same as integrate_formal, but use tau space instead of length for integration.
  !-----------------------------------------------------------------------------------------
  ! For a given initial intensity (I0) along a ray, this function will returns a intensity 
  ! value at an optical depth (tau) assuming that the source function is known. I.e., we solve 
  ! 
  !              dI
  !             ---- = S - I 
  !             dTau
  !
  !                         -(tau1-tau0)                        -(tau1-t)
  !    I(tau1) = I(tau0) * e               +   Integrate[(S(t)*e        , {t, tau0, tau1}]
  !
  !  
  !    0  I0                          1  I1
  !    *--> ------------------------- * -->
  !   
  !    tau=tau0                       tau=tau1
  !
  !  Normally, tau0=0.
  !
  ! Same as integral_formal, but this routine resample rays along ray 
  ! using dtau_max values which is the max increment of optical depth along a
  ! integrating path.
  !
  subroutine integrate_formal(I0, I1,  wavelength,  lambda0, mass_ion, aVec, uHat,  &
        continuum, grid, hitCore, thin_disc_on, opaqueCore,  sampleFreq, error, &
        do_alphadisc_check, alphadisc_par, thinLine, emissOff)      
    use amr_mod, only: amr_values_along_ray
    use utils_mod, only: bigGamma, voigtn, linearresample, linearresample_dble

    implicit none   
    real(double), intent(in)  :: I0                     ! Initial intensity
    real(double), intent(out) :: I1                     ! Final intensity
    real, intent(in)          :: wavelength             ! the wavelength  [A] 
    real, intent(in)          :: lambda0                ! rest wavelength of line [A]
    real, intent(in)          :: mass_ion               ! mass of atom in [g]
    type(VECTOR), intent(in)  :: aVec                   ! starting position vector
    type(VECTOR), intent(in)  :: uHat                   ! direction
    logical     , intent(in)  :: continuum              ! if T, just do contiuum (no line)
    type(GRIDTYPE), intent(in):: grid                   ! the opacity grid
    logical, intent(in)       :: thin_disc_on           ! T to include thin disc
    logical, intent(in)       :: opaqueCore             ! is the core opaque
    logical, intent(out)      :: hitcore                ! has the photon hit the core
    real, intent(in)          :: sampleFreq             ! max. samples per grid cell
    integer, intent(out)      :: error                  ! error code returned
    logical, intent(in)          :: do_alphadisc_check  ! if T the disc check will be done
    type(alpha_disc), intent(in) :: alphadisc_par       !      parameters for alpha disc
    logical, intent(in)       :: thinLine               ! if T surpress, the abs of cont photon by line
    logical, intent(in)       :: emissOff               ! if T surpress, the line emission

    !    
    integer                   :: nTau                   ! size of optical depth arrays
    type(VECTOR)         :: aVecOctal              ! VECTOR version of 'aVec'
    type(VECTOR)         :: uHatOctal              ! VECTOR version of 'uHat'
    type(VECTOR)              :: rVec                   ! position vector
    real(double) :: nu, nu_p, nu0, nu0_p           ! frequencies   
!    real :: thisVel
    integer :: ilambda                                ! wavelength index    
    real(double) :: chil                                       ! line opacity
    integer :: i, j                   ! counters
    
    
    real(double) :: deltaNu, DopplerWidth
    real(double) :: a
    real(double) :: Hay, dv,  Vrel
    real(double) :: T_mid, Ne_mid, N_HI_mid, chiline_mid, projVel_mid
    real(double) :: etaCont_mid, etaline_mid, etal
    real(double) :: chiAbs_mid, chiSca_mid, source, source_l, source_c
    integer :: newNtau
!    real :: T, lam, tmp
    real :: sqrt_pi
    real(double) :: dtau, taul, tau_tot, dtau_c, dtau_l
    real(double) :: int_source
    real ::  dtau_max
    logical :: fromDisc
    real(double) :: V_th ! thermal velocity
    !-------------------------------------------------------------------------
    ! WORK ARRAYS
    logical, save :: first_time = .true.
    integer, parameter       :: maxTau = 5000000  ! Use this for a normal opearation
!    integer, parameter       :: maxTau = 7000000  
    real(double), allocatable, save :: rho(:)         ! density (size=maxTau)
    real, allocatable, save :: temperature(:) ! temperature (size=maxTau)
    real(double), allocatable, save :: Ne(:)  ! electron density (size=maxTau)
    type(vector), allocatable, save :: velocity(:) ! size=maxTau       ! 
    real, allocatable, save :: chiLine(:)          ! line optical depth (size=maxTau)
    real, allocatable, save  :: tauAbsLine(:)      ! (size=maxTau)
    real, allocatable, save  :: velocityDeriv(:)   ! directional derivative (size=maxTau)
    real, allocatable, save :: dL(:)     ! distance increment array (size=maxTau)
    real(double), allocatable, save :: projVel(:)   ! (size=maxTau)
    real(double), allocatable, save :: ksca(:)      ! (size=maxTau)
    real(double), allocatable, save :: kabs(:)      ! (size=maxTau)
    real, allocatable, save :: newL(:)              ! (size=maxTau)
    real, allocatable, save :: N_HI(:)              ! size = maxTau
    real, allocatable, save :: tauAbs(:)            ! thermal optical depth
    real, allocatable, save :: tauSca(:)            ! e-scattering optical depth
    real(double), allocatable, save :: etaCont(:)   ! contiuum emissiovity  [check units]
    real(double), allocatable, save :: etaLine(:)   ! line emissivity       [check units]
    real, allocatable, save :: L(:)                 ! path distance array
    real, allocatable, save :: tau(:)               ! total optical depth
    logical, allocatable,save :: inflow(:)
    logical, allocatable,save :: newinflow(:)
    

    ! initialize variables
    sqrt_pi = SQRT(piDouble)
    hitcore = .false.
    rVec = aVec
    error = 0
   

     ! Allocates memory for work arrays for the first time
     ! This should be faster than using automatic arrays which allocates
     ! and deallocates memory every single time.  These array are SAVED
     ! so they should be available next time as well.
     if (first_time) then
        first_time = .false.
        ALLOCATE(rho(maxTau))         
        ALLOCATE(temperature(maxTau)) 
        ALLOCATE(Ne(maxTau))  
        ALLOCATE(velocity(maxTau)) 
        ALLOCATE(chiLine(maxTau)) 
        ALLOCATE(tauAbsLine(maxTau)) 
        ALLOCATE(velocityDeriv(maxTau)) 
        ALLOCATE(dL(maxTau))
        ALLOCATE(projVel(maxTau))
        ALLOCATE(kabs(maxTau))       
        ALLOCATE(ksca(maxTau))
        ALLOCATE(newL(maxTau))   
        ALLOCATE(N_HI(maxTau))
        ALLOCATE(tauAbs(maxTau))       
        ALLOCATE(tauSca(maxTau))
        ALLOCATE(tau(maxTau))       
        ALLOCATE(L(maxTau))
        ALLOCATE(etaCont(maxTau))
        ALLOCATE(etaLine(maxTau))
        ALLOCATE(inflow(maxTau))
        ALLOCATE(newinflow(maxTau))
     end if



    ! locate this wavelength in the grid of wavelengths
    
    if (grid%flatspec.or.(grid%doRaman)) then
       iLambda = 1
    else
       call locate(grid%lamArray, grid%nLambda, wavelength, iLambda)
    endif
    
    if (.not. grid%adaptive) then
       print *, 'integratePathAMR called on a non-adaptive grid!'
       stop
    end if
    
    nTau = 0
    
    ! convert some variables to the appropriate kind for startReturnSamples
    aVecOctal = aVec
    uHatOctal = uHat
    
    
    ! find the projected velocity [c]
!    Vn1 = uHat .dot. vVec         ! projected velocity of the local gas at emission location
!    thisVel = Vn1 + (lambda0-wavelength)/lambda0
!    thisVel = (wavelength-lambda0)/lambda0

    ! Note: temperature and Ne won't be needed here, but they have to be passed
    ! as arguments because NAG compiler do not like it!  (RK)
    ! If the ray intersect with thindisc, aVecOctal (initial position) is shifted
    ! to the intersection of the ray and the thindisc!
    call amr_values_along_ray(aVecOctal,uHatOctal,grid,sampleFreq,nTau,       &
         maxTau,thin_disc_on,opaqueCore,hitcore,fromDisc, .false., iLambda,error,&
         L,kappaAbs=kAbs,kappaSca=kSca,velocity=velocity,&
         velocityDeriv=velocityDeriv,chiLine=chiLine,    &
         rho=rho, temperature=temperature, Ne=Ne, inflow=inflow,&
         etaCont=etaCont, etaLine=etaLine)
    

    if (nTau <= 2) then
       !       print *, 'The code does not yet include routines for handling ',&
       !            'the case where nTau <= 2 (in integratePathAMR)'
       I1 = 1.0d-100
       error = -10
       return
    end if

    
    ! projected velocities 
    forall (i = 1:nTau)
       ! (+ve when moving toward in the direction of photon.)
       projVel(i) = dble(velocity(i) .dot. uHat)  
    end forall
    
    dL(1:nTau-1) = L(2:nTau) - L(1:nTau-1)
    ! Number density of HI N_H = N_HI + Np, but Np=Ne
    N_HI(1:ntau) = rho(1:ntau)/mass_ion - Ne(1:nTau)
    
    !
    !
    ! Find the initial optical depth scale
    !
    !    
    nu0 = cSpeed / (lambda0*angstromtocm)    ! line center frequency
    nu = cSpeed / (wavelength*angstromtocm)  ! freq of this photon
    nu_p = nu  ! freq in the rest frame of local gas
    !
    tauAbs(1) = 1.0e-25
    tauSca(1) = 1.0e-25
    tau(1) = 1.0e-25
    tauAbsLine(1) = 1.0e-25
    int_source = 0.0d0

    if (.not. hitcore) then
       do i = 2, nTau
          if (inflow(i-1)) then
             if (inflow(i)) then
                ! Evaluating the values in the mid point
                T_mid = 0.5d0*(temperature(i-1)+temperature(i))
                Ne_mid = 0.5d0*(Ne(i-1)+Ne(i))
                N_HI_mid = 0.5d0*(N_HI(i-1)+N_HI(i))
                chiline_mid = 0.5d0*(chiline(i-1)+chiline(i))
                etacont_mid = 0.5d0*(etacont(i-1)+etacont(i))
                etaline_mid = 0.5d0*(etaline(i-1)+etaline(i))
                projVel_mid = 0.5d0*(projVel(i-1)+projVel(i))
             else
                T_mid = temperature(i-1)
                Ne_mid = Ne(i-1)
                N_HI_mid = N_HI(i-1)
                chiline_mid = chiline(i-1)
                etacont_mid = etacont(i-1)
                etaline_mid = etaline(i-1)
                projVel_mid = projVel(i-1)
             end if
             
             T_mid = MAX(T_mid, 1000.0d0) ! [K]  To avoid a tiny Doppler width

             V_th = sqrt(2.*kErg*T_mid/mass_ion)  ! [cm/s] theram speed
             
             ! relative velocity wrt the observer
             Vrel = projVel_mid

             
             ! The line centre of absorption profile shifted by Doppler.
!             nu0_p = nu0/(1.0d0-Vrel)  ! [Hz] 
             nu0_p = nu0*(1.0d0+Vrel)  ! [Hz]  binominal expansion

             DopplerWidth = nu0_p/cSpeed * V_th !eq 7  [Hz]
             
             a = bigGamma(N_HI_mid, T_mid, Ne_mid, nu0_p) / (fourPi * DopplerWidth) ! [-]
             deltaNu = nu_p - nu0_p     !  [Hz]
             dv = deltaNu/DopplerWidth  ! [-]
             Hay = voigtn(a,dv)
             chil = chiLine_mid / (sqrt_pi*DopplerWidth) * Hay ! equation 5
             dtau = chil*dL(i-1)
             
             tauAbsLine(i) = tauAbsLine(i-1) +  abs(dtau)
             
             if (inflow(i)) then
                tauSca(i) = tauSca(i-1) + dL(i-1)*(0.5*(ksca(i)+ksca(i-1)))
                tauAbs(i) = tauAbs(i-1) + dL(i-1)*(0.5*(kabs(i)+kabs(i-1)))
             else
                tauSca(i) = tauSca(i-1) + dL(i-1)*ksca(i-1)
                tauAbs(i) = tauAbs(i-1) + dL(i-1)*kabs(i-1)
             end if
             
          else ! not inflow
             tauAbsLine(i) = tauAbsLine(i-1)
             tauSca(i) = tauSca(i-1)
             tauAbs(i) = tauAbs(i-1)             
          end if 

          if (continuum) then
             tau(i) = tauSca(i) + tauAbs(i)                  ! total optical depth
          else  ! line
             tau(i) = tauSca(i) + tauAbs(i) + tauAbsLine(i)  ! total optical depth
          end if

       enddo
       taul =  tauAbsLine(nTau)
       tau_tot =  tau(nTau)

    else  ! ray hits core for some reason (this should not happen though).
       write(*,*) "Error:: A ray hits the stellar core! [forma_solutions::formal_integral]."
       stop
    end if

    !------------------------------------------------------------------------
    ! Now resample rays using tau values
    !------------------------------------------------------------------------
    dtau_max = 0.05  ! Pick this one for a normal operation...
    call refine_ray_grid_by_tau(L, nTau, tau, dtau_max, maxTau, newL, newNTau, &
         inflow, newInFlow, do_alphadisc_check, alphadisc_par, aVecOctal,uHatOctal)
    
    
    ! Now interpolate on to newly sampled ray    
    call linearResample_dble(L, projVel, nTAu, maxTau, newL, newNtau)
    call linearResample_dble(L, kSca, nTAu, maxTau, newL, newNtau)
    call linearResample_dble(L, kAbs, nTAu, maxTau, newL, newNtau)
    call linearResample_dble(L, etaCont, nTAu, maxTau, newL, newNtau)
    call linearResample_dble(L, etaLine, nTAu, maxTau, newL, newNtau)
    call linearResample(L, chiline, nTAu, maxTau, newL, newNtau)
    call linearResample(L, temperature, nTAu, maxTau, newL, newNtau)
    call linearResample(L, N_HI, nTAu, maxTau, newL, newNtau)
    call linearResample_dble(L, Ne, nTAu, maxTau, newL, newNtau)


    ! updating the ray
    nTau = newNtau
    L(1:nTau) = newL(1:nTau)
    L(1) = 1.0e-25
    dL(1:nTau-1) = L(2:nTau) - L(1:nTau-1)
    inFlow(1:nTau) = newInFlow(1:nTau)  

!    Ne(1:ntau) = ABS(Ne(1:ntau))  ! just for safty
    ! Need to change it to below (because of g95 bug).
    do i = 1, ntau
       Ne(i) = ABS(Ne(i))
    end do


     
!     !---------------------------------------------------------------------------
!     ! Additional resampling for the place where the velocity is slowly changing
!     !---------------------------------------------------------------------------
!     call refine_ray_constant_vel(L, nTau, projVel, maxTau, newL, newNTau, &
!          inflow, newInFlow, do_alphadisc_check, alphadisc_par, aVecOctal,uHatOctal)
    
    
!     ! Now interpolate on to newly sampled ray    
!     call linearResample_dble(L, projVel, nTAu, maxTau, newL, newNtau)
!     call linearResample(L, kSca, nTAu, maxTau, newL, newNtau)
!     call linearResample(L, kAbs, nTAu, maxTau, newL, newNtau)
!     call linearResample_dble(L, etaCont, nTAu, maxTau, newL, newNtau)
!     call linearResample_dble(L, etaLine, nTAu, maxTau, newL, newNtau)
!     call linearResample(L, chiline, nTAu, maxTau, newL, newNtau)
!     call linearResample(L, temperature, nTAu, maxTau, newL, newNtau)
!     call linearResample(L, N_HI, nTAu, maxTau, newL, newNtau)
!     call linearResample_dble(L, Ne, nTAu, maxTau, newL, newNtau)


!      ! updating the ray
!      nTau = newNtau
!      L(1:nTau) = newL(1:nTau)
!      L(1) = 1.0e-25
!      dL(1:nTau-1) = L(2:nTau) - L(1:nTau-1)
!      inFlow(1:nTau) = newInFlow(1:nTau)  

!      Ne(1:ntau) = ABS(Ne(1:ntau))  ! just for safty

     !--------------------------------------------------------------------

    !
    ! Now perform the integration with improved ray samples
    !

    tauAbs(1) = 1.0e-25
    tauSca(1) = 1.0e-25
    tau(1) = 1.0e-25
    tauAbsLine(1) = 1.0e-25
    int_source = 0.0d0

!    do i = 2, nTau
    do i = nTau, 2, -1  ! integraing from the observer's side!!!
       j = nTau + 2 -i ! reversed order
       if (inflow(i)) then  
          if (inflow(i-1)) then  
             ! Evaluating the values in the mid point
             T_mid = 0.5d0*(temperature(i-1)+temperature(i))
             Ne_mid = 0.5d0*(Ne(i-1)+Ne(i))
             N_HI_mid = 0.5d0*(N_HI(i-1)+N_HI(i))
             chiline_mid = 0.5d0*(chiline(i-1)+chiline(i))
             etacont_mid = 0.5d0*(etacont(i-1)+etacont(i))
             etaline_mid = 0.5d0*(etaline(i-1)+etaline(i))
             projVel_mid = 0.5d0*(projVel(i-1)+projVel(i))
          else
             T_mid = temperature(i)
             Ne_mid = Ne(i)
             N_HI_mid = N_HI(i)
             chiline_mid = chiline(i)
             etacont_mid = etacont(i)
             etaline_mid = etaline(i)
             projVel_mid = projVel(i)
          end if
             
          T_mid = MAX(T_mid, 1000.0d0) ! [K]  To avoid a tiny Doppler width
          
          V_th = sqrt(2.*kErg*T_mid/mass_ion)  ! [cm/s] theram speed
          
          ! relative velocity wrt the observer
          Vrel = projVel_mid
          

          
          ! The line centre of absorption profile shifted by Doppler.
!          nu0_p = nu0/(1.0d0-Vrel)  ! [Hz] 
          nu0_p = nu0*(1.0d0+Vrel)  ! [Hz]  binominal expansion 
         
          DopplerWidth = nu0_p/cSpeed * V_th !eq 7  [Hz]
          
          a = bigGamma(N_HI_mid, T_mid, Ne_mid, nu0_p) / (fourPi * DopplerWidth) ! [-]
          deltaNu = nu_p - nu0_p     !  [Hz]
          dv = deltaNu/DopplerWidth  ! [-]
          Hay = voigtn(a,dv)
          chil = chiLine_mid / (sqrt_pi*DopplerWidth) * Hay ! equation 5
          dtau_l = chil*dL(i-1)
          
!          tauAbsLine(i) = tauAbsLine(i-1) +  abs(dtau_l)
          ! Now the tau is zero on the observer's side
          tauAbsLine(j) = tauAbsLine(j-1) +  abs(dtau_l)

          if (inflow(i-1)) then  
             chiSca_mid = 0.5*(ksca(i)+ksca(i-1))
             chiAbs_mid = 0.5*(kabs(i)+kabs(i-1))
          else
             chiSca_mid = ksca(i)
             chiAbs_mid = kabs(i)
          end if
          tauSca(j) = tauSca(j-1) + dL(i-1)*chiSca_mid
          tauAbs(j) = tauAbs(j-1) + dL(i-1)*chiAbs_mid
          dtau_c = (chiSca_mid + chiAbs_mid)*dL(i-1)
          dtau   = dtau_c + dtau_l

          
       else  ! not inflow
          tauAbsLine(j) = tauAbsLine(j-1)
          tauSca(j) = tauSca(j-1)
          tauAbs(j) = tauAbs(j-1)
          etaline_mid = 0.0d0
          etacont_mid = 0.0d0
          chiSca_mid =1.0d-25 
          chiAbs_mid = 1.0d-25 
          chil = 1.0d-25
          dtau_l = 1.0d-25
          dtau_c = 1.0d-25
          dtau   = 2.0d-25
       end if


       ! This option is mainly used for debugging
       if (emissOff) then 
          etaline_mid = 0.0d0
       end if

       ! For safty
       if (chiSca_mid<=0.0d0) chiSca_mid = 1.0e-30
       if (chiAbs_mid<=0.0d0) chiAbs_mid = 1.0e-30
       if (chil<=0.0d0) chil = 1.0e-30

       !
       ! Note: eta and chi have been multiplied by 10^10 in stateq_mod.f90, 
       !       so source = eta/chi should have correct units.

       if (continuum) then
          tau(j) = tauSca(j) + tauAbs(j)                  ! total optical depth
          source_c = etacont_mid/(chiSca_mid + chiAbs_mid)
          int_source = int_source +  source_c*EXP( -tau(j) )*dtau_c
       else         
          tau(j) = tauSca(j) + tauAbs(j) + tauAbsLine(j)  ! total optical depth
          ! integrating the source function
          etal = etaline_mid/(sqrt_pi*DopplerWidth) * Hay ! equation 5

          if (thinLine) then
             ! This option is mainly used for debugging
             ! surpresses the absorpion of contiuum by line
             source_c = etacont_mid/(chiSca_mid + chiAbs_mid)
             source_l = etal/chil
             int_source = int_source +  source_c*EXP( -tauSca(j)-tauAbs(j) )*dtau_c &
                  +  source_l*EXP( -tau(j) )*dtau_l
          else
             source = (etacont_mid + etal)/(chiSca_mid + chiAbs_mid + chil)
             int_source = int_source +  source*EXP( -tau(j) )*dtau
          end if
       end if
          
    enddo
    tau_tot =  tau(nTau)
    
    ! In case the ray originates from the stellar core and hits the disc.
    ! The new starting position was chosen at the intersection of the disc 
    ! (in amr_values_along_ray), but the initial intensity has not been adjusted.
    ! It must be set to zero when this occurs.
    if (fromDisc) then 
       I1 = int_source    
    else
       ! final intensity
       if (thinLine) then
          I1 = I0*EXP(-tauAbs(nTau)-tauSca(nTau)) + int_source    
       else
          I1 = I0*EXP(-tau_tot) + int_source    
       end if
    end if
    

  end subroutine integrate_formal



  !
  !
  ! Creates NDF image of flux map
  subroutine create_obs_flux_map(flux_ray, ray, dr, dphi, &
        nray, R_star, R_max, filename, npix)
    use image_mod, only: IMAGETYPE, freeImage, initImage
#ifdef USECFITSIO
    use image_mod, only : writeFitsImage
#endif
    implicit none
    integer, intent(in) :: nray  ! total number of rays
    real(double), intent(in) :: flux_ray(nray) ! Flux map  array dimension= nlam x #of integration points  
    ! integration grids (original ray in setup_grid and setup_core_grid program before roatation.)
    type(VECTOR), intent(in) :: ray(nray)      
    real(double), intent(in)   ::  dr(nray)    ! [10^10 cm]  radial increments
    real(double), intent(in)   ::  dphi      ! [rad]       anglular increment
    real(double), intent(in)   ::  R_star    ! [10^10 cm]  radius of sta
    real(double), intent(in)   ::  R_max     ! [10^10 cm]  max radius for plot
    character(LEN=*),intent(in)   :: filename  ! output filename
    integer, intent(in)           :: npix ! image size = 2*npix +1 should be an even number
    !
    !
    type(imagetype) ::  F_map  ! flux map   []
    integer :: i, j, ix, iy , k, nsize
    real :: x, y, flux, del, flux_min, offset
    character(len=80) :: message
    message = filename

    offset = 1.0e-4*R_star
    ! initializes the map
    F_map = initImage(npix, npix, 2.*REAL(R_max), 2.*REAL(R_max),-1.0, 1.0)

    ! now add flux values to this map
    nsize = 2*npix+1
    del = REAL(2.0d0*R_max)/REAL(nsize)

    flux_min = MAX(1.0d-23, MINVAL(flux_ray) )! limit the minimum value

    ! Very naive search algorithm here! (Slow)
    ! This section should be improved later.    
    do j =1, nsize
       y = -R_max + REAL(j-1)*del + del/2.0 + offset ! mid point
       iy = -npix + (j-1)
       do i =1, nsize
          x = -R_max + REAL(i-1)*del + del/2.0 + offset  ! mid point
          ix = -npix + (i-1)      
          flux = 0.0
          do k = 1, nray
!!!!!             if (in_this_area(dble(x), dble(y), ray(k), dA(k), dr(k), dphi, dir_obs) )then
             ! Flips the value of x to be consistent with Monte Calro's notation....
             ! which is the projected from the center of star not from the observer!
             ! Uncomment above to create the map in observe's perspective.
             if (in_this_area(dble(-x), dble(y), ray(k),  dr(k), dphi) )then
                flux = MAX(flux_min, REAL(flux_ray(k))) + flux
!                exit  ! this loop
             end if
          end do
          F_map%pixel(ix,iy)%i = MAX(flux, flux_min/2.0)
          F_map%pixel(ix,iy)%q = 1.0e-20
          F_map%pixel(ix,iy)%u = 1.0e-20
          F_map%pixel(ix,iy)%v = 1.0e-20
       end do
    end do

    ! now save the image to a flie

#if USECFITSIO
    call writeFitsImage(F_map, filename, 1.d0, "intensity")
#endif

    ! clean up
    call freeImage(F_map)




  end subroutine create_obs_flux_map


  !
  !  
  !  
  !
  ! Use this to simulate the spectro-astrometry data
  !  For a better modelling, implement the followings later:
  !  1. Slit width (now set to infinte size).
  !  2. Slit direction with respect to the stellar coordinates.
  !  3. Combine this with create_obs_flux_map for speed up! (a naive search use
  !     here is very slow!
  subroutine find_position_displacement(flux_map, ray, dA, dr, dphi, &
       nlam, nray,  R_star, R_max, filename, npix, &
       dist_obj, lam, lam0)
    implicit none
    integer, intent(in) :: nlam  ! number of wavelength points
    integer, intent(in) :: nray  ! total number of rays
    real(double), intent(in) :: flux_map(:, :) ! Flux map  array dimension= nlam x #of integration points  
    ! integration grids (original ray in setup_grid and setup_core_grid program before roatation.)
    type(VECTOR), intent(in) :: ray(nray)      
    real(double), intent(in)   ::  dA(:)  ! [10^20 cm]  surface element used for flux integration
    real(double), intent(in)   ::  dr(:)    ! [10^10 cm]  radial increments
    real(double), intent(in)   ::  dphi      ! [rad]       anglular increment
    real(double), intent(in)   ::  R_star    ! [10^10 cm]  radius of sta
    real(double), intent(in)   ::  R_max     ! [10^10 cm]  max radius for plot
    character(LEN=*),intent(in)   :: filename  ! output filename
    integer, intent(in)           :: npix    ! image size = 2*npix +1 should be an even number
    real(double), intent(in)   :: dist_obj   ! physical distance to the object in [cm]
    real(single), intent(in)   :: lam(:)  ! wavelength array [A]
    real(single), intent(in)   :: lam0       ! line wavelength [A]
    !
    !
    integer :: i, j, ix, k, nsize, m
    real :: x, y, offset
    real(double), allocatable ::  FL(:)  ! allocate it with nlam line flux (F-Fc)/Fc
    real(double)  :: lam_pos                ! flux weighted mean wavelenth [A]
    real(double), allocatable:: vel_pos(:)  ! size = nsize flux weighted mean velocity [km/s]
    real(double), allocatable:: pos(:)      ! size = nsize  displacement position [10^10 cm]
    real(double) :: dlam, sum_FL_lam_dlam, sum_FL_dlam
    real(double) :: pixsize           ! pixsize in size in [10^10 cm]
    real(double) :: convert_to_mas    ! conversion factor from 10^10cm to miliarcsec
    real(double) :: convert_to_AU     ! conversion factor from 10^10cm to AU
    real(double), parameter :: c=2.99792458d10 ! speed of light [cm/s]
!    character(len=80) :: tempChar, pos_disp_file
!    real(double) :: vel, Fn
    !
    ! now add flux values to this map
    nsize = 2*npix+1
    pixsize = REAL(2.0d0*R_max)/REAL(nsize)

    !
    ! creating the work arrays
    ALLOCATE(FL(nlam))
    ALLOCATE(vel_pos(nsize))
    ALLOCATE(pos(nsize))

    offset = 1.0e-4*R_star
    

    ! finding the conversion factor.
    convert_to_AU = 1.0d0/1.495979d3  
    convert_to_mas = (1.0d10 / dist_obj) * radtoDeg * 3600.d0 * 1000.d0   

    ! wavelength bin width
    dlam = lam(3)-lam(2) ! [A]
    

    ! Very naive search algorithm here! (Slow)
    ! This section should be improved later.    
    do j =1, nsize
       y = -R_max + REAL(j-1)*pixsize + pixsize/2.0 + offset ! mid point
       pos(j) = y
       FL(:) = 0.0d0
       do i =1, nsize
          x = -R_max + REAL(i-1)*pixsize + pixsize/2.0 + offset  ! mid point
          ix = -npix + (i-1)      
          do k = 1, nray
             if (in_this_area(dble(x), dble(y), ray(k), dr(k), dphi) )then
                ! integrating the profile acroess the x-axis
                FL(1:nlam) = FL(1:nlam) + flux_map(1:nlam,k)/dA(k)
             end if
          end do
       end do
       
       ! set the min values to avoid devided by zero.
       FL(1) = MAX(FL(1), 1.0d-23)


!        ! write the normalized spectrum
!        write(tempChar,'(i3.3)') j
!        pos_disp_file = TRIM(filename)//"_fn_pos_disp_"//TRIM(tempChar)//".dat"
! 5      format(2(1PE15.4, 2x))
!        open(unit=122, file=TRIM(pos_disp_file), status="replace")
!        write(122,'(a)')  "# Written by [formal_solutions::find_position_displacement]."
!        write(122,'(a10,2x,1PE15.4,2x, a2)')  "# OBJECT DISTANCE = ", dist_obj/3.0856776d18, "pc"
!        write(122,'(a)')  "# Format ::  velocity [km/s] --  Normalized flux "
!        do m = 2, nlam
!           Fn =  FL(m) / FL(1)
!           vel = (lam(m)-lam0)/lam0 * c*1.d-5   ! [km/s]
!           write(122, 5)  vel,  Fn
!        END do
!        close(122)


       ! removing the continuum contribution
       FL(2:nlam) = ( FL(2:nlam) - FL(1) )  /   FL(1)

       ! finding the flux weighted mean wavelength                
       sum_FL_lam_dlam=0.0d0; sum_FL_dlam=0.0d0
       do m = 2, nlam
          sum_FL_lam_dlam = sum_FL_lam_dlam + FL(m)*lam(m)*dlam
          sum_FL_dlam = sum_FL_dlam + FL(m)*dlam
       end do
       lam_pos = sum_FL_lam_dlam / sum_FL_dlam           ! [A]
       vel_pos(j) = (lam_pos-lam0)/lam0 * c*1.d-5        ! [km/s]

    end do

    ! writing a result in a file
    
    open(unit=124, file=filename//"_pos_disp.dat", status="replace")    
    write(124,'(a)')  "# Written by [formal_solutions::find_position_displacement]."
    write(124,'(a10,2x,1PE15.4,2x, a2)')  "# OBJECT DISTANCE = ", dist_obj/3.0856776d18, "pc"
    write(124,'(a)')  "# Format ::  velocity [km/s] --  dispacement [AU] -- diplacement [mas] "
10  format(3(1PE15.4, 2x))
    do j = 1, nsize
       write(124,10) vel_pos(j), pos(j)*convert_to_AU, pos(j)*convert_to_mas
    end do
    close(124)

    ! rotated positions (-45 deg) in AU
    open(unit=124, file=filename//"_pos_disp_au_45deg.dat", status="replace")
    write(124,'(a)')  "# Written by [formal_solutions::find_position_displacement]."
    write(124,'(a10,2x,1PE15.4,2x, a2)')  "# OBJECT DISTANCE = ", dist_obj/3.0856776d18, "pc"
    write(124,'(a)')  "# Format ::  dispacement x [AU] -- diplacement y [AU] "
14  format(2(1PE15.4, 2x))
    do j = 1, nsize
       write(124,14) pos(j)*convert_to_AU/sqrt(2.0d0), -pos(j)*convert_to_AU/sqrt(2.0d0)
    end do
    close(124)

    ! rotated positions (-45 deg) in mas
    open(unit=124, file=filename//"_pos_disp_mas_45deg.dat", status="replace")
    write(124,'(a)')  "# Written by [formal_solutions::find_position_displacement]."
    write(124,'(a10,2x,1PE15.4,2x, a2)')  "# OBJECT DISTANCE = ", dist_obj/3.0856776d18, "pc"
    write(124,'(a)')  "# Format ::  dispacement x [mas] -- diplacement y [mas] "
    do j = 1, nsize
       write(124,14) pos(j)*convert_to_mas/sqrt(2.0d0), -pos(j)*convert_to_mas/sqrt(2.0d0)
    end do
    close(124)

    ! clean up
    DEALLOCATE(FL)
    DEALLOCATE(VEL_POS)
    DEALLOCATE(POS)



  end subroutine find_position_displacement




  !
  ! Give a test point, ray, radial, angular step sizes, this routine
  ! will find whether the test point belongs to this flux integration area.
  ! Projected plane is (px-py)  (note this is different from the stellar coordinate)
  ! The axis ade defined as: px = zhat x dir_obs and py = dir_obs x px
  !
  logical function in_this_area(test_px, test_py, p, dr, dphi)
    implicit none
    real(double), intent(in) :: test_px, test_py  ! [10^10] test point coordinates
    type(VECTOR), intent(in)   ::  p        ! integration grids
    real(double),      intent(in)   ::  dr       ! [10^10 cm]  radial increments
    real(double),      intent(in)   ::  dphi     ! [rad]       anglular increment
    !
!    type(VECTOR) ::  px, py        ! unit vector on the project plane coordinates
    real(double) :: p_phi, p_r, pxx, pyy, r_min, r_max, phi_min, phi_max
    real(double) :: test_r, test_phi, pi
    
    ! initialize the value
    in_this_area = .false.
    pi = 2.0d0*ACOS(0.0d0)

!     ! Defining the reference axis (Everthing will be projected on this plane.)
!     px = (VECTOR(0.0d0, 0.0d0, 1.0d0) .cross. dir_obs)
! !    px = dir_obs .cross. VECTOR(1.0d0, 0.0d0, 0.0d0) 
!     call normalize(px)  ! just in case...
!     py = dir_obs .cross. px
!     call normalize(py)  ! just in case...


    ! transforming to polar coordinates
    test_r = SQRT(test_px*test_px + test_py*test_py)
    test_phi = ATAN2(test_py, test_px) 
    if (test_phi < 0.0d0) test_phi=test_phi+2.0d0*pi ! range=> (0, 2pi)

    ! now for the ray point
!    pxx = px .dot. p
!    pyy = py .dot. p 
    ! In the original setup of rays in setup_grid and setup_core_grid, all the rays are in 
    ! +x direction. 
    pxx = p%y  ! 
    pyy = p%z  !
    p_r = SQRT(pxx*pxx + pyy*pyy)

    ! Check if in radial range
    r_min = p_r - dr/2.0d0; r_max = p_r + dr/2.0d0
    if ( test_r >= r_min .and. test_r < r_max) then
       continue  ! 
    else
       in_this_area = .false.
       return
    end if

    ! Now checks the angular range
    !
    p_phi = ATAN2(pyy, pxx)
    if (p_phi<0.0d0) p_phi = p_phi + 2.0d0*pi ! range=> (0, 2pi)

    ! now see if this point is within the area segment
    phi_min = p_phi - dphi/2.0d0; phi_max = p_phi + dphi/2.0d0
    if (phi_min < 0.d0) phi_min = phi_min + 2.0d0*pi
    if (phi_max < 0.d0) phi_max = phi_max + 2.0d0*pi
    if (phi_min > 2.0d0*pi) phi_min = phi_min - 2.0d0*pi
    if (phi_max > 2.0d0*pi) phi_max = phi_max - 2.0d0*pi
    
    if (phi_min > phi_max) then 
       ! the area must contain phi=0 line
       if (test_phi > phi_min .and. test_phi <= 2.0d0*pi) then
          in_this_area = .true.
       elseif (test_phi < phi_max .and. test_phi >=0.0d0) then
          in_this_area = .true.
       else
          in_this_area = .false.
       end if
          
    else
       if (test_phi >= phi_min .and. test_phi < phi_max) then
          in_this_area = .true.
       else
          in_this_area = .false.
       end if
    end if
    return

  end function in_this_area



  !
  !
  ! Given a value of dtau_max, this rotutine splits
  ! the ray segments until optical depth scale in each segment is 
  ! less than dtau_max. 
  ! If the ray segment is in the alpha disc, the ray will not split...
  !
  subroutine refine_ray_grid_by_tau(L, nTau, tau, dtau_max, maxtau, newL, newNTau,&
         inFlow, newInFlow, do_alphadisc_check, alphadisc_par, pos_start, dir_ray)
      integer, intent(in) :: nTau, maxtau
      real, intent(in) :: L(nTau)
      real, intent(in) :: tau(nTau)
      ! maxmum line optical depth increment allowed
      real, intent(in) :: dtau_max
      integer, intent(inout) :: newNtau
      real, intent(inout)  :: newL(maxtau)
      logical, intent(inout)  :: InFlow(nTau)
      logical, intent(inout)  :: newInFlow(maxtau)
      logical, intent(in) :: do_alphadisc_check      ! if T the disc check will be done
      type(alpha_disc), intent(in) :: alphadisc_par  ! parameters for alpha disc
      type(VECTOR), intent(in) :: pos_start     ! startingpoint
      type(VECTOR), intent(in) :: dir_ray       ! direction of ray
      !
      type(VECTOR) :: pos     
      logical :: unsafe_position, ray_starts_inside_disc
      real :: dtau  ! 
      integer :: nAdd 
      integer :: i, j
!      real, parameter :: dvel  = 10.e5/cSpeed
      real ::  dL
 


      !
      newNTau = 0
      unsafe_position = .false.
      ray_starts_inside_disc = .false.
      do i = 1, nTau-1
         ! Initial optical depth in each line sequemnt
         dtau = tau(i+1) - tau(i)
         if (do_alphadisc_check) then
            pos = pos_start + dble(L(i))*dir_ray ! position of ray segment.
            unsafe_position = in_alpha_disc(alphadisc_par, pos)
            ! Checks if the ray start from within the disc
            if (i==1 .and. unsafe_position) ray_starts_inside_disc = .true.            
         end if
         if (dtau > dtau_max .and. (.not. unsafe_position)) then
            nAdd = nint(dtau/dtau_max)  ! this should be at least 1
            dL = (L(i+1)-L(i))/real(nAdd+1)
            do j = 1, nAdd+1
               newNtau = newNtau + 1
               if (newntau > maxtau) then
                  print *, "Error: newntau > maxtau in refine_ray_grid_by_tau. (1)"
                  print *, "         maxtau = ", maxtau
                  print *, "         newtau = ", newNtau
                  stop
               end if
               newL(newNTau) = real(j-1)*dL + L(i)
               newInFlow(newNTau) = inFlow(i)
            enddo
         else
            newNtau = newNtau + 1
            if (newntau > maxtau) then
               print *, "Error: newntau > maxtau in refine_ray_grid_by_tau. (2)"
               stop
            end if
            newL(newNTau) = L(i)
!            if (ray_starts_inside_disc) then
!               ! ========== WARNING =========================================
!               ! Here we assume that a ray does not cross the disc twice and
!               ! there is no emissvity from the disc.
!               ! If so and if the ray begins from within the disc, we ignore
!               ! the grids in the disc to avoid the underflow by EXP(-tau) term!
!               ! ========== WARNING =========================================               
!               newInFlow(newNTau) = .false.
!            else
               newInFlow(newNTau) = inFlow(i)
!            end if
               
         endif
      enddo
      newNtau = newNtau + 1
      newL(newNTau) = L(nTau)
      newInFlow(newNTau) = inFlow(nTau)
    end subroutine refine_ray_grid_by_tau


    !
    !
    ! Add extra points where velocity is changing slowly
    subroutine refine_ray_constant_vel(L, nTau, projVel,  maxtau, newL, newNTau,&
         inFlow, newInFlow, do_alphadisc_check, alphadisc_par, pos_start, dir_ray)

      implicit none

      integer, intent(in) :: nTau, maxtau
      real, intent(in) :: L(nTau)
      real(double), intent(in) :: projVel(nTau)
      integer, intent(inout) :: newNtau
      real, intent(inout)  :: newL(maxtau)
      logical, intent(inout)  :: InFlow(nTau)
      logical, intent(inout)  :: newInFlow(maxtau)
      logical, intent(in) :: do_alphadisc_check      ! if T the disc check will be done
      type(alpha_disc), intent(in) :: alphadisc_par  ! parameters for alpha disc
      type(VECTOR), intent(in) :: pos_start     ! startingpoint
      type(VECTOR), intent(in) :: dir_ray       ! direction of ray
      !
      real(double) :: dProjVel  ! automatic array
      real(double) :: dlam
      integer :: nAdd =5
      integer :: i, j
!      real, parameter :: dvel  = 10.e5/cSpeed
      real :: dvel
      logical :: unsafe_position
      type(VECTOR) :: pos     


      ! the projected speed increment is smaller than dvel
      ! we add nAdd points in between.      
      dvel = 10.e5/cspeed ! 10 km/s

      newNTau = 0
      unsafe_position = .false.
      do i = 2, nTau
         dProjVel = projVel(i) - projVel(i-1)
         if (do_alphadisc_check) then
            pos = pos_start + dble(L(i))*dir_ray ! position of ray segment.
            unsafe_position = in_alpha_disc(alphadisc_par, pos)
         end if

         if (abs(dProjVel) < dVel .and. (.not. unsafe_position)) then
!            nAdd = nint(abs(dProjVel)/dVel)
            dlam = (L(i)-L(i-1))/real(nAdd+1)
            do j = 1, nAdd+1
               newNtau = newNtau + 1
               newL(newNTau) = real(j-1)*dlam + L(i-1)
               newInFlow(newNTau) = inFlow(i-1)
            enddo
         else
            newNtau = newNtau + 1
            newL(newNTau) = L(i-1)
            newInFlow(newNTau) = inFlow(i-1)
         endif
      enddo
      newNtau = newNtau + 1
      newL(newNTau) = L(nTau)
      newInFlow(newNTau) = inFlow(nTau)
    end subroutine refine_ray_constant_vel






!   !
!   !  Based on integrate_formal
!   !  Given starting point (aVec) and the direction (uHat), this routine computes: 
!   !  souce function along a ray.
!   !
!   subroutine source_function_along_ray(wavelength,  lambda0, mass_ion, aVec, uHat,  &
!         continuum, grid, hitCore, thin_disc_on, opaqueCore,  sampleFreq, error, &
!         r_norm, filename)        
!     implicit none   
!     real, intent(in)          :: wavelength             ! the wavelength  [A] 
!     real, intent(in)          :: lambda0                ! rest wavelength of line [A]
!     real, intent(in)          :: mass_ion               ! mass of atom in [g]
!     type(VECTOR), intent(in)  :: aVec                   ! starting position vector
!     type(VECTOR), intent(in)  :: uHat                   ! direction
!     logical     , intent(in)  :: continuum              ! if T, just do contiuum (no line)
!     type(GRIDTYPE), intent(in):: grid                   ! the opacity grid
!     logical, intent(in)       :: thin_disc_on           ! T to include thin disc
!     logical, intent(in)       :: opaqueCore             ! is the core opaque
!     logical, intent(out)      :: hitcore                ! has the photon hit the core
!     real, intent(in)          :: sampleFreq             ! max. samples per grid cell
!     integer, intent(out)      :: error                  ! error code returned
!     ! radius/size scale in [10^10cm] for normalizing length and density
!     real(double), intent(in)    :: r_norm               
!     character(LEN=*),intent(in) :: filename               ! output filename
!     !    
!     integer                   :: nTau                   ! size of optical depth arrays
!     type(VECTOR)         :: aVecOctal              ! VECTOR version of 'aVec'
!     type(VECTOR)         :: uHatOctal              ! VECTOR version of 'uHat'
!     type(VECTOR)              :: rVec                   ! position vector
!     real(double) :: nu, nu_p, nu0, nu0_p           ! frequencies   
! !    real :: thisVel
!     integer :: ilambda                                ! wavelength index    
!     real(double) :: chil                                       ! line opacity
!     integer :: i, j                   ! counters
    
    
!     real(double) :: deltaNu, DopplerWidth
!     real(double) :: a
!     real(double) :: Hay, dv,  Vrel
!     real(double) :: T_mid, Ne_mid, N_HI_mid, chiline_mid, projVel_mid
!     real(double) :: etaCont_mid, etaline_mid, etal
!     real :: T, lam, tmp
!     real :: sqrt_pi
!     logical :: fromDisc
!     real(double) :: V_th ! thermal velocity
!     !-------------------------------------------------------------------------
!     ! WORK ARRAYS
!     integer, parameter       :: maxTau = 500000  ! Use this for a normal opearation
!     real, allocatable, save :: rho(:)         ! density (size=maxTau)
!     real, allocatable, save :: temperature(:) ! temperature (size=maxTau)
!     real(double), allocatable, save :: Ne(:)  ! electron density (size=maxTau)
!     type(vector), allocatable, save :: velocity(:) ! size=maxTau       ! 
!     real, allocatable, save :: chiLine(:)          ! line optical depth (size=maxTau)
!     real, allocatable, save  :: velocityDeriv(:)   ! directional derivative (size=maxTau)
!     real(double), allocatable, save :: projVel(:)   ! (size=maxTau)
!     real, allocatable, save :: ksca(:)              ! (size=maxTau)
!     real, allocatable, save :: kabs(:)              ! (size=maxTau)
!     real, allocatable, save :: newL(:)              ! (size=maxTau)
!     real, allocatable, save :: N_HI(:)              ! size = maxTau
!     real(double), allocatable, save :: etaCont(:)   ! contiuum emissiovity  [check units]
!     real(double), allocatable, save :: etaLine(:)   ! line emissivity       [check units]
!     real, allocatable, save :: L(:)                 ! path distance array
!     logical, allocatable,save :: inflow(:)
    

!     ! initialize variables
!     sqrt_pi = SQRT(piDouble)
!     hitcore = .false.
!     rVec = aVec
!     error = 0
   

!      ! Allocates memory for work arrays 
!     ALLOCATE(rho(maxTau))         
!     ALLOCATE(temperature(maxTau)) 
!     ALLOCATE(Ne(maxTau))  
!     ALLOCATE(velocity(maxTau)) 
!     ALLOCATE(chiLine(maxTau)) 
!     ALLOCATE(velocityDeriv(maxTau)) 
!     ALLOCATE(projVel(maxTau))
!     ALLOCATE(kabs(maxTau))       
!     ALLOCATE(ksca(maxTau))
!     ALLOCATE(newL(maxTau))   
!     ALLOCATE(N_HI(maxTau))
!     ALLOCATE(L(maxTau))
!     ALLOCATE(etaCont(maxTau))
!     ALLOCATE(etaLine(maxTau))
!     ALLOCATE(inflow(maxTau))
!     ALLOCATE(newinflow(maxTau))



!     ! open file and write headers
!     open(unit =55, file=TRIM(filename), status="replace")
!     write(55, "(a)") "# format::  distance (cm)  --- Normalized source function (s(r), s(r=)"

!     ! locate this wavelength in the grid of wavelengths
    
!     if (grid%flatspec.or.(grid%doRaman)) then
!        iLambda = 1
!     else
!        call locate(grid%lamArray, grid%nLambda, wavelength, iLambda)
!     endif
    
!     if (.not. grid%adaptive) then
!        print *, 'integratePathAMR called on a non-adaptive grid!'
!        stop
!     end if
    
!     nTau = 0
    
!     ! convert some variables to the appropriate kind for startReturnSamples
!     aVecOctal = aVec
!     uHatOctal = uHat
    
    
!     ! find the projected velocity [c]
! !    Vn1 = uHat .dot. vVec         ! projected velocity of the local gas at emission location
! !    thisVel = Vn1 + (lambda0-wavelength)/lambda0
! !    thisVel = (wavelength-lambda0)/lambda0

    
!     call amr_values_along_ray(aVecOctal,uHatOctal,grid,sampleFreq,nTau,       &
!          maxTau,thin_disc_on,opaqueCore,hitcore,fromDisc, .false., iLambda,error,&
!          L,kappaAbs=kAbs,kappaSca=kSca,velocity=velocity,&
!          velocityDeriv=velocityDeriv,chiLine=chiLine,    &
!          rho=rho, temperature=temperature, Ne=Ne, inflow=inflow,&
!          etaCont=etaCont, etaLine=etaLine)
!         ! Note: temperature and Ne won't be needed here, but they have to be passed
!         ! as arguments because NAG compiler do not like it!  (RK)
    

!     if (nTau <= 2) then
!        !       print *, 'The code does not yet include routines for handling ',&
!        !            'the case where nTau <= 2 (in integratePathAMR)'
!        I1 = 1.0d-100
!        error = -10
!        return
!     end if

    
!     ! projected velocities 
!     forall (i = 1:nTau)
!        ! (+ve when moving toward in the direction of photon.)
!        projVel(i) = dble(velocity(i) .dot. uHat)  
!     end forall
    
!     ! Number density of HI N_H = N_HI + Np, but Np=Ne
!     N_HI(1:ntau) = rho(1:ntau)/mass_ion - Ne(1:nTau)
    
!     !
!     !
!     ! Find the initial optical depth scale
!     !
!     !    
!     nu0 = cSpeed / (lambda0*angstromtocm)    ! line center frequency
!     nu = cSpeed / (wavelength*angstromtocm)  ! freq of this photon
!     nu_p = nu  ! freq in the rest frame of local gas
!     !
!     if (.not. hitcore) then
!        do i = 2, nTau
!           if (inflow(i-1)) then
!              ! Evaluating the values in the mid point
!              T_mid = 0.5d0*(temperature(i-1)+temperature(i))
!              Ne_mid = 0.5d0*(Ne(i-1)+Ne(i))
!              N_HI_mid = 0.5d0*(N_HI(i-1)+N_HI(i))
!              chiline_mid = 0.5d0*(chiline(i-1)+chiline(i))
!              etacont_mid = 0.5d0*(etacont(i-1)+etacont(i))
!              etaline_mid = 0.5d0*(etaline(i-1)+etaline(i))
!              projVel_mid = 0.5d0*(projVel(i-1)+projVel(i))

!              if (continuum) then
             
!              else
!                 T_mid = MAX(T_mid, 1000.0d0) ! [K]  To avoid a tiny Doppler width
                
!                 V_th = sqrt(2.*kErg*T_mid/mass_ion)  ! [cm/s] theram speed
                
!                 ! relative velocity wrt the observer
!                 Vrel = projVel_mid
                
                
!                 ! The line centre of absorption profile shifted by Doppler.
!                 nu0_p = nu0/(1.0d0-Vrel)  ! [Hz] 
                
!                 DopplerWidth = nu0_p/cSpeed * V_th !eq 7  [Hz]
                
!                 a = bigGamma(N_HI_mid, T_mid, Ne_mid, nu0_p) / (fourPi * DopplerWidth) ! [-]
!                 deltaNu = nu_p - nu0_p     !  [Hz]
!                 dv = deltaNu/DopplerWidth  ! [-]
!                 Hay = voigtn(a,dv)
!                 chil = chiLine_mid / (sqrt_pi*DopplerWidth) * Hay ! equation 5
!                 etal = etaline_mid/(sqrt_pi*DopplerWidth) * Hay ! equation 5
             
!              end if


!        enddo

!     else  ! ray hits core for some reason (this should not happen though).
!        write(*,*) "Error:: A ray hits the stellar core! [formal_soltion::source_function_along_ray]."
!        stop
!     end if

    
!     ! CLEANING UP
!     DEALLOCATE(rho)
!     DEALLOCATE(temperature)
!     DEALLOCATE(Ne)
!     DEALLOCATE(velocity)
!     DEALLOCATE(chiLine)
!     DEALLOCATE(tauAbsLine)
!     DEALLOCATE(velocityDeriv)
!     DEALLOCATE(projVel)
!     DEALLOCATE(kabs)
!     DEALLOCATE(ksca)
!     DEALLOCATE(N_HI)
!     DEALLOCATE(L)
!     DEALLOCATE(etaCont)
!     DEALLOCATE(etaLine)
!     DEALLOCATE(inflow)
    

!     close(55)


!   end subroutine source_function_along_ray




end module formal_solutions

