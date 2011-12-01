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

