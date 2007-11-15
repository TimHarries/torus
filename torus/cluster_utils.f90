module cluster_utils

  ! This is a collection of functions an subroutines related to an cluster_class object. 
  ! To avoid some circular dependency, these routines where separated from the main 
  ! cluster_class.f90 



  use formal_solutions
  use cluster_class
  use vector_mod
  use grid_mod
  use atom_mod
  use source_mod
  use filter_set_class
  use kind_mod  
  use path_integral

  implicit none

  public ::  &
       &    analyze_cluster, &
       &    observed_flux
  

  private:: &
       UKIRT_color_magnitudes, &
       compute_colors, &
       setup_grid, &
       compute_Av
  



contains


  !
  ! Given the direction to an observer, a cluster object, a grid, and 
  ! the distance to the cluster, this routine compute the JHKL'M'(UKIRT)
  ! magnitudes the stars and the flux at the pass band of SIRTF
  ! (instruments:MIPS & IRAC).  
  subroutine analyze_cluster(a_cluster, dir_obs, distance, grid)
    implicit none
    type(cluster), intent(in) :: a_cluster 
    type(octalvector), intent(in) :: dir_obs
    real(double), intent(in) :: distance  ! in [pc]
    type(gridtype), intent(in) :: grid    
    !
    type(filter_set) :: SIRTF_MIPS
    type(filter_set) :: SIRTF_IRAC 
    type(filter_set) :: UKIRT_JHKLM
    
    ! making filters
    call make_filter_set(SIRTF_MIPS, 'mips')
    call make_filter_set(SIRTF_IRAC, 'irac')
    call make_filter_set(UKIRT_JHKLM, 'ukirt')
    
    ! compute UKIRT magnitudes and colors
    call compute_Av(a_cluster, dir_obs, distance, grid)
!    call ukirt_color_magnitudes(a_cluster, dir_obs, distance, grid)
!    call compute_colors(a_cluster, dir_obs, distance, grid, UKIRT_JHKLM)
    
    ! compute SIRTF flux and colors
!    call compute_colors(a_cluster, dir_obs, distance, grid, SIRTF_MIPS)
!    call compute_colors(a_cluster, dir_obs, distance, grid, SIRTF_IRAC)

  end subroutine analyze_cluster




  !
  ! Given the direction to an observer, a cluster object, a grid, and 
  ! the distance to the cluster, this routine compute the JHKL'M'(UKIRT)
  ! magnitudes the stars in a cluster. 
  ! It also computes Av (extinction)
  ! The results will be written in a output file.
  subroutine ukirt_color_magnitudes(a_cluster, dir_obs, distance, grid)
    implicit none

    type(cluster), intent(in) :: a_cluster 
    type(octalvector), intent(in) :: dir_obs
    real(double), intent(in) :: distance  ! in [pc]
    type(gridtype), intent(in) :: grid    
    !
    integer :: i
    real(double) :: I0_J, I0_H, I0_K, I0_L, I0_M, I0_V  ! initial intensity
!    real(double) :: I_J, I_H, I_K, I_L, I_M, I_V       ! final intensity
    real(quad) :: I_J, I_H, I_K, I_L, I_M, I_V       ! final intensity
    integer :: nstar
    type(sourcetype) :: a_star
    real(double) :: T  !   temperature in K.    
    type(octalvector) :: position   ! position of stars [10^10cm]
    real(double) :: R            ! radius of star    [10^10cm]
    !
    ! The mean wavelength of the filters
    real(double), parameter :: lam_J =  1.25d4 ! [A]
    real(double), parameter :: lam_H =  1.65d4 ! [A]
    real(double), parameter :: lam_K =  2.2d4  ! [A]
    real(double), parameter :: lam_L =  3.6d4  ! [A]
    real(double), parameter :: lam_M =  10.2d4 ! [A]
    real(double), parameter :: lam_V =  0.55d4 ! [A]
    
    !
    ! zero point flux for each filter.
    real(double), parameter :: offset = 1.0d100
    !
    real(double), parameter :: F0_J =  3.18d-10*offset ! [erg cm^-2 s^-1 A^-1 *offset]
    real(double), parameter :: F0_H =  1.18d-10*offset ! [erg cm^-2 s^-1 A^-1 *offset]
    real(double), parameter :: F0_K =  4.17d-11*offset ! [erg cm^-2 s^-1 A^-1 *offset]
    real(double), parameter :: F0_L =  6.23d-12*offset ! [erg cm^-2 s^-1 A^-1 *offset]
    real(double), parameter :: F0_M =  2.07d-12*offset ! [erg cm^-2 s^-1 A^-1 *offset]
    real(double), parameter :: F0_V =  3.64d-09*offset ! [erg cm^-2 s^-1 A^-1 *offset]

    ! 
    ! Observed flux and magnitudes
!    real(double) :: F_J, F_H, F_K, F_L, F_M, F_V  ! [erg cm^-2 s^-1 A^-1]
    real(quad) :: F_J, F_H, F_K, F_L, F_M, F_V  ! [erg cm^-2 s^-1 A^-1 *offset]
    real(double) :: m_J, m_H, m_K, m_L, m_M, m_V  !  magnitudes of real stars
    real(double) :: m_J0, m_H0, m_K0, m_L0, m_M0  ! mag. of naked stars
    real(double), allocatable :: m_V0(:)
    !
    real(double) :: Av ! extinction

    !
    real(double) :: scale,  Rmax, tau_inf, tmp


    integer, parameter :: nr    = 30
    integer, parameter :: nphi  = 10
    !
    integer, parameter :: LU_MAG = 32
    integer, parameter :: LU_COL = 33    

    real(double), parameter :: pi = 3.141592654d0


    
    ! finding the number of stars in this cluster 
    nstar = get_nstar(a_cluster)


    allocate(m_V0(nstar))


    !
    !  Finding the magnitudes for naked stars
    !
    ! Loop over all stars in the cluster
    do i = 1, nstar

       ! retriving a star from the cluster 
       a_star = get_a_star(a_cluster, i)

       ! temperature of stars.
       T = a_star%Teff  ! [K]

       ! position of the star
       position = a_star%position  ! vector [10^10cm]
       
       ! radius of star 
       R = a_star%radius  ! [10^10cm]

    
       ! Finding the intensity.  (Assuming BB radiation for now..)
       ! Using the function in atom_mod.f90
       I_J = bLambda(lam_J, T) * 1.d-8 * offset  ! [erg cm^-2 s^-2 A^-1 sr^-1]
       I_H = bLambda(lam_H, T) * 1.d-8 * offset  ! [erg cm^-2 s^-2 A^-1 sr^-1]
       I_K = bLambda(lam_K, T) * 1.d-8 * offset  ! [erg cm^-2 s^-2 A^-1 sr^-1]
       I_L = bLambda(lam_L, T) * 1.d-8 * offset  ! [erg cm^-2 s^-2 A^-1 sr^-1]
       I_M = bLambda(lam_M, T) * 1.d-8 * offset  ! [erg cm^-2 s^-2 A^-1 sr^-1]
       I_V = bLambda(lam_V, T) * 1.d-8 * offset  ! [erg cm^-2 s^-2 A^-1 sr^-1]

       scale = ( R / distance/3.08568e8)**2  

       
       ! Note F(R_*) = Pi*B
       F_J = I_J * scale * Pi ! [erg cm^-2 s^-2 A^-1]
       F_H = I_H * scale * Pi ! [erg cm^-2 s^-2 A^-1]
       F_K = I_K * scale * Pi ! [erg cm^-2 s^-2 A^-1]
       F_L = I_L * scale * Pi ! [erg cm^-2 s^-2 A^-1]
       F_M = I_M * scale * Pi ! [erg cm^-2 s^-2 A^-1]
       F_V = I_V * scale * Pi ! [erg cm^-2 s^-2 A^-1]
              
       ! magnitudes
       
       m_J0 = 2.51d0 * ( LOG10(F0_J) - LOG10(F_J) )
       m_H0 = 2.51d0 * ( LOG10(F0_H) - LOG10(F_H) )
       m_K0 = 2.51d0 * ( LOG10(F0_K) - LOG10(F_K) )
       m_L0 = 2.51d0 * ( LOG10(F0_L) - LOG10(F_L) )
       m_M0 = 2.51d0 * ( LOG10(F0_M) - LOG10(F_M) )

       m_V0(i) = 2.51d0 * ( LOG10(F0_V) - LOG10(F_V) )
                     
       ! writing the results in files.

11     format(a14, 1PD13.3, 1x, a3)
12     format(a,  7(a8))
13     format(1x, i8, 6(f8.2))
14     format(a,  5(a8))
15     format(1x, i8, 4(f8.2))
       
16     format(a,  8(a8))
17     format(1x, i8, 7(f8.2))

       if (i==1) then ! open files and write headers
          open(unit=LU_MAG, file ="ukirt_app_mag_naked.dat", status="replace")
          open(unit=LU_COL, file ="ukirt_color_naked.dat", status="replace")

          ! magnitudes file
          write(LU_MAG, '(a)') "# Created by cluster_class::color_magnitudes."
          write(LU_MAG, '(a)') "# ======= NAKED STARS ======================="
          write(LU_MAG,11) "# distance  = ", distance, " pc"
          write(LU_MAG,12) "#", "star ID", "J", "H", "K", "L", "M", "V"

          ! color-color 
          write(LU_COL, '(a)') "# Created by cluster_class::color_magnitudes."
          write(LU_COL,11) "# distance  = ", distance, " pc"
          write(LU_COL,14) "#", "star ID", "H-K", "J-H", "J-K", "K-L"    
       end if
          

       write(LU_MAG, 13) i, m_J0, m_H0, m_K0, m_L0, m_M0, m_V0(i)
       write(LU_COL, 15) i, (m_H0-m_K0), (m_J0-m_H0), (m_J0-m_K0), (m_K0-m_L0)

    end do

    close(LU_MAG)
    close(LU_COL)



    ! Loop over all stars in the cluster
    do i = 1, nstar

       ! retriving a star from the cluster 
       a_star = get_a_star(a_cluster, i)

       ! temperature of stars.
       T = a_star%Teff  ! [K]      

       ! position of the star
       position = a_star%position  ! vector [10^10cm]
       
       ! radius of star 
       R = a_star%radius  ! [10^10cm]

    
       ! Finding the initial intensity.  (Assuming BB radiation for now..)
       ! Using the function in atom_mod.f90
       I0_J = bLambda(lam_J, T) * 1.d-8  ! [erg cm^-2 s^-2 A^-1 sr^-1]
       I0_H = bLambda(lam_H, T) * 1.d-8  ! [erg cm^-2 s^-2 A^-1 sr^-1]
       I0_K = bLambda(lam_K, T) * 1.d-8  ! [erg cm^-2 s^-2 A^-1 sr^-1]
       I0_L = bLambda(lam_L, T) * 1.d-8  ! [erg cm^-2 s^-2 A^-1 sr^-1]
       I0_M = bLambda(lam_M, T) * 1.d-8  ! [erg cm^-2 s^-2 A^-1 sr^-1]
       I0_V = bLambda(lam_V, T) * 1.d-8  ! [erg cm^-2 s^-2 A^-1 sr^-1]



       ! Flux units here are [erg cm^-2 s^-2 A^-1 * offset]
       Rmax = 1.5d13 ! [cm] ... = 1 AU
       F_J = observed_flux(I0_J, lam_J, R*1.d10, Rmax, &
            position,  nr, nphi, dir_obs, distance*3.08568d18, grid, .true., offset)  
       F_H = observed_flux(I0_H, lam_H, R*1.d10, Rmax, &
            position,  nr, nphi, dir_obs, distance*3.08568d18, grid, .true., offset)  
       F_K = observed_flux(I0_K, lam_K, R*1.d10, Rmax, &
            position,  nr, nphi, dir_obs, distance*3.08568d18, grid, .true., offset)  
       F_L = observed_flux(I0_L, lam_L, R*1.d10, Rmax, &
            position,  nr, nphi, dir_obs, distance*3.08568d18, grid, .true., offset)  
       F_M = observed_flux(I0_M, lam_M, R*1.d10, Rmax, &
            position,  nr, nphi, dir_obs, distance*3.08568d18, grid, .true., offset)  
       F_V = observed_flux(I0_V, lam_V, R*1.d10, Rmax, &
            position,  nr, nphi, dir_obs, distance*3.08568d18, grid, .true., offset)  


       
              
       ! magnitudes
       ! SEE page 100 of "Handbook of Space Astronomy and Astrophysics" by Martin V. Zombeck (1990).
       
       m_J = 2.51d0 * ( LOG10(F0_J) - LOG10(F_J) )
       m_H = 2.51d0 * ( LOG10(F0_H) - LOG10(F_H) )
       m_K = 2.51d0 * ( LOG10(F0_K) - LOG10(F_K) )
       m_L = 2.51d0 * ( LOG10(F0_L) - LOG10(F_L) )
       m_M = 2.51d0 * ( LOG10(F0_M) - LOG10(F_M) )
       m_V = 2.51d0 * ( LOG10(F0_V) - LOG10(F_V) )
       
       ! extinction
       !Av = m_V - m_V0(i)
       tau_inf=0.0d0
       tmp = formal_sol_dust_AMR(I0_V, 0.0d0, lam_V, position, dir_obs, &
            grid, .true., offset, tau_inf)
       Av = 1.1d0*tau_inf

       ! writing the results in files.
       
       if (i==1) then ! open files and write headers
          open(unit=LU_MAG, file ="ukirt_app_mag.dat", status="replace")
          open(unit=LU_COL, file ="ukirt_color.dat", status="replace")
          open(unit=39, file ="Av.dat", status="replace")

          ! magnitudes file
          write(LU_MAG, '(a)') "# Created by cluster_class::color_magnitudes."
          write(LU_MAG, '(a)') "# ===== REDDEN STARS ========================"
          write(LU_MAG,11) "# distance  = ", distance, " pc"
          write(LU_MAG,16) "#", "star ID", "J", "H", "K", "L", "M", "V", "Av"

          ! color-color 
          write(LU_COL, '(a)') "# Created by cluster_class::color_magnitudes."
          write(LU_COL,11) "# distance  = ", distance, " pc"
          write(LU_COL,14) "#", "star ID", "H-K", "J-H", "J-K", "K-L"
         
       end if
          

       write(LU_MAG, 17) i, m_J, m_H, m_K, m_L, m_M, m_V, Av
       write(LU_COL, 15) i, (m_H-m_K), (m_J-m_H), (m_J-m_K), (m_K-m_L)
       write(39,*) i, Av

    end do


    close(LU_MAG)
    close(LU_COL)


    deallocate(m_V0)  ! this should be automatic ...


  end subroutine ukirt_color_magnitudes




  !
  ! Given the direction to an observer, a cluster object, a grid, filter_set, and 
  ! the distance to the cluster,this routine computes fluxs at the pass bands of 
  ! the filter. 
  ! It also computes Av (extinction)
  ! The results will be written in a output file.
  subroutine compute_colors(a_cluster, dir_obs, distance, grid, filters)
    implicit none

    type(cluster), intent(in) :: a_cluster 
    type(octalvector), intent(in) :: dir_obs
    real(double), intent(in) :: distance  ! in [pc]
    type(gridtype), intent(in) :: grid    
    type(filter_set), intent(in) :: filters 
    !
    integer :: j, k
    real(double), allocatable :: I0(:) ! initial intensity
    real(double), allocatable :: I(:)  ! final intensity
    integer :: nstar
    type(sourcetype) :: a_star
    real(double) :: T  !   temperature in K.    
    type(octalvector) :: position   ! position of stars [10^10cm]
    real(double) :: R            ! radius of star    [10^10cm]
    ! 
    ! number of filters in a set
    integer :: nfilter
    
    ! The mean wavelength and the band width of the filters
    real(double), allocatable :: lambda(:)      ! [A]
    real(double), allocatable :: dlambda(:)     ! [A]
    !
    ! Observed flux and magnitudes
    real(double), allocatable :: F_lambda(:)     ! [erg cm^-2 s^-1 A^-1]
    real(double), allocatable :: F_Jy(:)         ! [erg cm^-2 s^-1 A^-1]
    real(double), parameter :: offset = 1.0d100
    !
    real(double) :: scale

    integer, parameter :: LU_COL = 32
    integer, parameter :: LU_FLX = 33    

    real(double), parameter :: pi = 3.141592654d0
    real(double), parameter :: c_in_m_per_sec = 2.998d8 ! [m/s]

    real(double) :: Rmax
    real(double) :: convert
    !
    character(LEN=30) :: name_filter_set
    character(LEN=30), allocatable :: name_filter(:)
    character(LEN=30) :: filename_flx, filename_col
    character(LEN=30) :: fmt13, fmt14, fmt16, fmt17
    !
    integer, parameter :: nr    = 50
    integer, parameter :: nphi  = 50

    
    ! Finding the number of filters in the given set
    nfilter = get_nfilter(filters)  
   
    ! finding the number of stars in this cluster 
    nstar = get_nstar(a_cluster)

    ! allocating memory for som arrays
    allocate(I(nfilter), I0(nfilter))
    allocate(F_lambda(nfilter), F_Jy(nfilter))
    allocate(lambda(nfilter), dlambda(nfilter))
    allocate(name_filter(nfilter))

    !
    ! Finding the effective wavelength, and the band width of each filter in the set.
    ! -- using the routines in filter_set_class.f90 
    do j = 1, nfilter
       lambda(j) = lambda_eff_filters(filters,j)   ! [A]
       dlambda(j) = 0.5d0*FWHM_filters(filters,j)  ! [A]
    end do
    
    ! Extracting the name of filters.
    name_filter_set = get_filter_set_name(filters)
    do j = 1, nfilter
       name_filter(j) = get_filter_name(filters, j) 
    end do

    
    ! Some formats for outout  ==========================
11  format(a14, 1PD13.3, 1x, a3)
    
    write(fmt13, *) nfilter-1
    fmt13 = '(a, a10, a18, '//TRIM(ADJUSTL(fmt13))//'(a18)'//')' 
    
    write(fmt14, *) nfilter-2
    fmt14 = '(a, a10, a18, '//TRIM(ADJUSTL(fmt14))//'(a18)'//')'
       
    write(fmt16, *) nfilter-1
    fmt16 = '(1x, i10, '//TRIM(ADJUSTL(fmt16))//'(1pe18.3))'//')'
    
    write(fmt17, *) nfilter-2
    fmt17 = '(1x, i10, '//TRIM(ADJUSTL(fmt17))//'(1pe18.3))'//')'
    !=================================================

    !
    !  Finding the magnitudes for naked stars
    !
    ! Loop over all stars in the cluster
    do k = 1, nstar

       ! retriving a star from the cluster 
       a_star = get_a_star(a_cluster, k)

       ! temperature of stars.
       T = a_star%Teff  ! [K]

       ! position of the star
       position = a_star%position  ! vector [10^10cm]
       
       ! radius of star 
       R = a_star%radius  ! [10^10cm]
    
       ! Finding the intensity.  (Assuming BB radiation for now..)
       ! Using the function in atom_mod.f90 
       scale = ( R / distance/3.08568e8)**2  

       do j = 1, nfilter
          I(j) = bLambda(lambda(j), T) * 1.d-8   ! [erg cm^-2 s^-1 A^-1 sr^-1] 

          ! Note F(R_*) = Pi*B
          F_lambda(j) = I(j) * scale * Pi ! [erg cm^-2 s^-1 A^-1]

          ! Conversion factor (1 Jy = 9.57e-13 erg/s/cm^2/A at 5600 A)
          ! and the fact that 
          !             F_lambda[erg/s/cm^2/A] = lambda^2/c*F_nu [erg/s/cm^2/Hz] 
          !                                    = lambda^2/c*1.0e-23 F_nu[Jy]
          convert =  c_in_m_per_sec/((lambda(j)*1.d-10)**2) * 1.0d-23 * 1.0d-10 

          ! Flux in [Jy]
!          F_Jy(j) = F_lambda(j) / convert / dlambda(j)
          F_Jy(j) = F_lambda(j) / convert  ! it is already per Angstrom!

       end do
              

       if (k==1) then ! open files and write headers
          filename_flx = TRIM(name_filter_set)//"_flux_naked.dat"
          filename_col = TRIM(name_filter_set)//"_color_color_naked.dat"

          open(unit=LU_FLX, file =TRIM(filename_flx), status="replace")
          open(unit=LU_COL, file =TRIM(filename_col), status="replace")

          ! flux file
          write(LU_FLX, '(a)') "#=== Created by cluster_class::compute_colors =============="
          write(LU_FLX, '(a)') "# ======= FOR NAKED STARS === FLUX IN [erg cm^-2 s^-1 A^-1]  ===== ========"
          write(LU_FLX,11) "# distance  = ", distance, " pc"          
          write(LU_FLX,TRIM(fmt13)) "#", "star ID", (TRIM(name_filter(j)), j=2,nfilter)

          ! color-color 
          write(LU_COL, '(a)') "# Created by cluster_class::color_magnitudes."
          write(LU_COL,11)     "# distance  = ", distance, " pc"
          write(LU_COL,TRIM(fmt14)) "#", "star ID",  &
               ("log10("//TRIM(name_filter(j))//"/"//TRIM(name_filter(j+1))//")", &
               j=2,nfilter-1)
       end if
          
       write(LU_FLX, TRIM(fmt16)) k, (F_lambda(j), j=2,nfilter)
       write(LU_COL, TRIM(fmt17)) k, (log10(F_lambda(j))-log10(F_lambda(j+1)), j=2,nfilter-1)

    end do

    close(LU_FLX)
    close(LU_COL)


    !
    ! Now solves the formal solutions.....
    !  
    do k = 1, nstar

       ! retriving a star from the cluster 
       a_star = get_a_star(a_cluster, k)

       ! temperature of stars.
       T = a_star%Teff  ! [K]

       ! position of the star
       position = a_star%position  ! vector [10^10cm]
       
       ! radius of star 
       R = a_star%radius  ! [10^10cm]
    
       ! Finding the intensity.  (Assuming BB radiation for now..)
       ! Using the function in atom_mod.f90 

       scale = ( R / distance/3.08568e8)**2  

       do j = 1, nfilter
          ! initial intensity
          I0(j) = bLambda(lambda(j), T) * 1.d-8   ! [erg cm^-2 s^-2 A^-1 sr^-1] 

          !
          ! Finding the final intensity
          ! Using the routine in formal_solution.f90
          ! usints are in  [erg cm^-2 s^-2 cm^-1 sr^-1 * offset]
          !
!!$          I(j) = formal_sol_dust_AMR(I0(j), 0.0d0, lambda(j), position, dir_obs, grid, .true., offset)
!!$          
!!$
!!$          ! Note F(R_*) = Pi*B
!!$          F_lambda(j) = I(j) * scale * Pi     ! [erg cm^-2 s^-2 A^-1 * offset]


          Rmax = 1.5d13 ! [cm] ... = 1 AU

          F_lambda(j) = observed_flux(I0(j), lambda(j), R*1.d10,  Rmax, &
               position,  nr, nphi, dir_obs, distance*3.08568d18, grid, .true., offset)
          ! [erg cm^-2 s^-2 A^-1 * offset]

          F_lambda(j) =  F_lambda(j)/offset   ! [erg cm^-2 s^-2 A^-1]

          ! Conversion factor (1 Jy = 9.57e-13 erg/s/cm^2/A at 5600 A)
          ! and the fact that 
          !             F_lambda[erg/s/cm^2/A] = lambda^2/c*F_nu [erg/s/cm^2/Hz] 
          !                                    = lambda^2/c*1.0e-23 F_nu[Jy]
          convert =  c_in_m_per_sec/((lambda(j)*1.d-10)**2) * 1.0d-23 * 1.0d-10 

          ! Flux in [Jy]
!          F_Jy(j) = F_lambda(j) / convert / dlambda(j)
          F_Jy(j) = F_lambda(j) / convert   ! It is already in per Angstrom

       end do
              

       if (k==1) then ! open files and write headers
          filename_flx = TRIM(name_filter_set)//"_flux.dat"
          filename_col = TRIM(name_filter_set)//"_color_color.dat"

          open(unit=LU_FLX, file =TRIM(filename_flx), status="replace")
          open(unit=LU_COL, file =TRIM(filename_col), status="replace")

          ! flux file
          write(LU_FLX, '(a)') "#=== Created by cluster_class::compute_colors.==========="
          write(LU_FLX, '(a)') "# ======= REDDEN STARS === FLUX IN [erg cm^-2 s^-1 A^-1] ==========="
          write(LU_FLX,11) "# distance  = ", distance, " pc"          
          write(LU_FLX,TRIM(fmt13)) "#", "star ID", (TRIM(name_filter(j)), j=2,nfilter)

          ! color-color 
          write(LU_COL, '(a)') "# Created by cluster_class::color_magnitudes."
          write(LU_COL,11)     "# distance  = ", distance, " pc"
          write(LU_COL,TRIM(fmt14)) "#", "star ID",  &
               ("log10("//TRIM(name_filter(j))//"/"//TRIM(name_filter(j+1))//")", &
               j=2,nfilter-1)
       end if
          
       write(LU_FLX, TRIM(fmt16)) k, (F_lambda(j), j=2,nfilter) 
       write(LU_COL, TRIM(fmt17)) k, (log10(F_lambda(j))-log10(F_lambda(j+1)), j=2,nfilter-1) 

    end do

    close(LU_FLX)
    close(LU_COL)


  end subroutine compute_colors


  !
  ! Computes the observed flux from a star for a given direction of observer
  ! and the distance to observer.
  ! Flux units in  [erg cm^-2 s^-2 cm^-1 * offset]

  function observed_flux(I0, wavelength, R_star, R_max, pos_star, &
       nr, nphi, dir_obs, dist_obs, amrgrid, contPhoton, offset) RESULT(F_obs)

    implicit none
    
    real(quad) ::  F_obs    ! output intensity.
    !
    real(double), intent(in)  ::  I0  ! [erg cm^-2 s^-2 A^-1 sr^-1] input intensity
    real(double), intent(in)  :: wavelength    ! [A] the wavelength 
    real(double), intent(in)  :: R_star        ! [cm] radius of a star
    real(double), intent(in)  :: R_max         ! [cm] the wavelength (usually ~10 AU)
    integer, intent(in)                :: nr            ! number of radial points for integration
    integer, intent(in)                :: nphi          ! number of angle  points for integration
    type(OCTALVECTOR), intent(in)      :: pos_star      ! position of the star
    real(double), intent(in)  :: dist_obs      ! [cm] distance to the observer
    type(OCTALVECTOR), intent(in)      :: dir_obs       ! direction
    type(GRIDTYPE), intent(in)         :: amrgrid       ! the opacity grid
    logical, intent(in)                :: contPhoton    ! is this a continuum photon?
    real(double), intent(in)  :: offset        ! offset scale factor
    !
    type(OCTALVECTOR), allocatable     :: p(:)          ! integration grids
    real(double), allocatable :: dA(:)         ! surface elements at p
    type(OCTALVECTOR)                  :: q             
    real(quad)                :: F, F_core, Fi, I1, F_sub
    !
    real(double)              :: pi
    integer :: i, np, nc
    integer, parameter    :: nr_core   = 50     ! number of radial points for integration
    integer, parameter    :: nphi_core = 20     ! number of angle  points for integration
    

    pi = 2.0*ACOS(0.0)

    ! allocating the grid points
    nc = nr_core*nphi_core
    ALLOCATE(p(nc), dA(nc))
       
    !
    ! start integration    
    !

    !    ! core flux
    !    I1 = formal_sol_dust_AMR(I0, 0.0d0, wavelength, pos_star, dir_obs, amrgrid, .true., offset)
    !    F_core = I1*pi*R_star**2


    ! setting up the grid for integration
    call setup_grid(p,dA, nr_core, nphi_core, 0.01d0*R_star, 0.99d0*R_star, pos_star, dir_obs)    
    
!$OMP PARALLEL DEFAULT(NONE) & 
!$OMP PRIVATE(i, I1, Fi, F_sub, q)&
!$OMP SHARED(nc, F_core, wavelength, p, dA, dir_obs, amrgrid, offset, I0, R_star, pos_star) 
    F_sub = 0.0
    F_core = 0.0d0  ! core flux is computed in the following loops.
!!$OMP DO SCHEDULE(STATIC)
!$OMP DO SCHEDULE(DYNAMIC)
    do i = 1, nc
       q = VECTOR(p(i)%x, p(i)%y, pos_star%z+R_star/1.0d10)
       I1 = formal_sol_dust_AMR(I0, 0.0d0, wavelength, q, dir_obs, amrgrid, .true., offset)
       Fi = I1*dA(i)
       F_sub = Fi + F_sub
    end do
!$OMP END DO  
!$OMP CRITICAL (update)
    F_core = F_sub + F_core  ! sum from each thread for OpenMP
!$OMP END CRITICAL (update)
!$OMP END PARALLEL

    DEALLOCATE(p)
    DEALLOCATE(dA)
    

    !
    ! Flux in outer parts
    ! 

    np = nr*nphi
    ALLOCATE(p(np), dA(np))

    ! setting up the grid for integration
    call setup_grid(p,dA, nr, nphi, R_star, R_max, pos_star, dir_obs)


!$OMP PARALLEL DEFAULT(NONE) & 
!$OMP PRIVATE(i, I1, Fi, F_sub)&
!$OMP SHARED(np, F, wavelength, p, dA, dir_obs, amrgrid, offset,  pos_star) 
    F_sub = 0.0
    F = 0.0
!!$OMP DO SCHEDULE(STATIC)
!$OMP DO SCHEDULE(DYNAMIC)
    do i = 1, np
       I1 = formal_sol_dust_AMR(0.0d0, 0.0d0, wavelength, p(i), dir_obs, amrgrid, .true., offset)
       Fi = I1*dA(i)
       F_sub = Fi + F_sub
    end do
!$OMP END DO  
!$OMP CRITICAL (update)
    F = F_sub + F  ! sum from each thread for OpenMP
!$OMP END CRITICAL (update)
!$OMP END PARALLEL


    !
    ! Total Flux 
    !

    ! scaling the flux at the observer's distance.
    F_obs = (F+F_core) / (dist_obs*dist_obs)  ! [erg cm^-2 s^-2 cm^-1 * offset]
    

    ! Cleaning up.
    DEALLOCATE(p)
    DEALLOCATE(dA)

  end function observed_flux



  !
  ! Routine to setup the grid used in the observed flux integration (observed_flux)
  !
  ! THIS IS STILL DEVELOPMENTAL STAGE!!
  subroutine setup_grid(p, dA,nr, nphi, R_min, R_max, pos_star, dir_obs)
    implicit none
    type(OCTALVECTOR), intent(inout)     :: p(:)          ! integration grids
    real(double), intent(inout) :: dA(:)         ! surface elements at p
    integer, intent(in)                :: nr            ! number of radial points for integration
    integer, intent(in)                :: nphi          ! number of angle  points for integration
    real(double), intent(in)  :: R_min         ! [cm] minimum radius (usually 5*R_star)
    real(double), intent(in)  :: R_max         ! [cm] the wavelength (usually ~10 AU)
    
    type(OCTALVECTOR), intent(in)      :: pos_star      ! position of the sta
    type(OCTALVECTOR), intent(in)      :: dir_obs       ! direction
    
    integer :: i , j , k
    type(OCTALVECTOR)                  :: q             ! 
    real(double)              :: pi
    real(double)              :: r, phi, dphi, dr, rp, rm
    real(double)              :: log_r, log_R_max, log_R_min
    real(oct)              :: x, y, z
    
    pi = 2.0*ACOS(0.0)

    ! ===============================================================
    ! FOR NOW the observer is assumed to be on Z axis 
    ! THIS SHOULD BE GENERAALIZED LATER!!
    log_R_min = LOG(R_min)
    log_R_max = LOG(R_max)
    dphi = 2.0*pi/dble(nphi-1)
    i = 0
    do k = 1, nr
       log_r = (log_R_max - log_R_min)*dble(k-1)/dble(nr-1)  + log_R_min
       r = EXP(log_r)/1.d10  ! 10^10cm
       log_r = (log_R_max - log_R_min)*dble(k-1+0.5)/dble(nr-1)  + log_R_min
       rp = EXP(log_r)/1.d10  ! 10^10cm
       log_r = (log_R_max - log_R_min)*dble(k-1-0.5)/dble(nr-1)  + log_R_min
       rm = EXP(log_r)/1.d10  ! 10^10cm
       dr = rp-rm             ! 10^10cm

       phi = 0
       do j = 1, nphi
          i = i + 1
          phi = phi + dble(j-1)*dphi
          
          x = r*COS(phi);  y = r*SIN(phi);   z = - R_max/1.0d10 ! [10^10cm]
          q = VECTOR(x, y, z)
          
          p(i) = pos_star + q  ! (vector addition) [10^10cm]
          dA(i) = r*dphi*dr *1.0d20    ! [cm^2]
          
       end do
    end do
          
  end subroutine setup_grid


  !
  ! Computing extinction
  !

  subroutine compute_Av(a_cluster, dir_obs, distance, grid)
    implicit none

    type(cluster), intent(in) :: a_cluster 
    type(octalvector), intent(in) :: dir_obs
    real(double), intent(in) :: distance  ! in [pc]
    type(gridtype), intent(in) :: grid    
    !
    integer :: i
    real(double) :: I0_V  ! initial intensity
    integer :: nstar
    type(sourcetype) :: a_star
    real(double) :: T  !   temperature in K.    
    type(octalvector) :: position   ! position of stars [10^10cm]
    real(double) :: R            ! radius of star    [10^10cm]
    !
    ! The mean wavelength of the filters
    real(double), parameter :: lam_V =  0.55d4 ! [A]
    
    real(double), parameter :: offset = 1.0d100
    !
    real(double) :: Av ! extinction
    !
    real(double) ::  tau_inf, tmp

    !
    !
    type(VECTOR)  :: coolStarPosition
    
    integer, parameter :: maxTau  = 1000
    integer, parameter :: nLambda = 1000
    ! local variables
    real :: lambda(maxTau)    ! path distance array
    real :: tauExt(maxTau)    ! optical depth
    real :: tauAbs(maxTau)    ! optical depth
    real :: tauSca(maxTau)    ! optical depth
    real :: tauCont(maxTau,nLambda)
    !
    
    integer  :: nTau        ! size of optical depth arrays
    real     :: escProb     ! the escape probability
    logical  :: hitcore     ! has the photon hit the core
    integer  :: error       ! error code returned
    logical  :: useinterp =.false.
    logical      :: opaqueCore  =.true.       ! is the core opaque
    real         :: lamStart, lamEnd
    !
    logical      :: thinLine=.false.           ! ignore line absorption of continuum
    integer      :: nUpper = 3
    integer      :: nLower =2
    real         :: sampleFreq=2.0        ! max. samples per grid cell
    real :: junk        


    ! finding the number of stars in this cluster 
    nstar = get_nstar(a_cluster)


    ! Loop over all stars in the cluster
    do i = 1, nstar

       ! retriving a star from the cluster 
       a_star = get_a_star(a_cluster, i)

       ! temperature of stars.
       T = a_star%Teff  ! [K]      

       ! position of the star
       position = a_star%position  ! vector [10^10cm]
       
       ! radius of star 
       R = a_star%radius  ! [10^10cm]

    
       ! Finding the initial intensity.  (Assuming BB radiation for now..)
       ! Using the function in atom_mod.f90
       I0_V = bLambda(lam_V, T) * 1.d-8  ! [erg cm^-2 s^-2 A^-1 sr^-1]
       
       ! extinction
       !Av = m_V - m_V0(i)
       tau_inf=0.0d0
!       call IntegratePathAMR(real(lam_V),  real(lam_V), OCTALVECTOR(1.,1.,1.), &
!            position, &
!            dir_obs, grid, lambda, tauExt, tauAbs, &
!            tauSca, maxTau, nTau, opaqueCore, escProb, .true. , &
!            lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., &
!            .false., nUpper, nLower, 0., 0., 0., junk,&
!            sampleFreq,error)
!
!       tau_inf = tauExt(ntau)


       tau_inf=0.0d0
       tmp = formal_sol_dust_AMR(I0_V, 0.0d0, lam_V, position, dir_obs, &
            grid, .true., offset, tau_inf)
       Av = 1.1d0*tau_inf

       ! writing the results in files.
       
       if (i==1) then ! open files and write headers
          open(unit=39, file ="Av.dat", status="replace")         
       end if
          
       write(39,*) i, Av

    end do

    close(39)


  end subroutine compute_Av








end module cluster_utils
