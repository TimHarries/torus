module formal_solutions

  ! ----------------------------------------------------------
  ! This is not a class difenition, but a collection of formal
  ! solution solvers.
  ! 
  !-----------------------------------------------------------
  ! created : aug-08-2003  (R. Kurosawa)
  !-----------------------------------------------------------


  use vector_mod             ! vector maths 
  use gridtype_mod           ! opacity grid
  use path_integral          ! 
  use octal_mod
  use kind_mod

  
  public :: formal_sol_dust_AMR


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

    implicit none
    
    real(kind=quadkind) ::  I1    ! output intensity.
    !
    real(kind=doublekind), intent(in)  ::  I0    ! input intensity
    real(kind=doublekind), intent(in)  ::  tau0  ! optical depth at the initial point 

    real(kind=doublekind), intent(in)          :: wavelength        ! the wavelength [A]
    type(OCTALVECTOR), intent(in)  :: aVec          ! starting position vector
    type(OCTALVECTOR), intent(in)  :: uHat          ! direction
    type(GRIDTYPE), intent(in)     :: grid          ! the opacity grid
    logical, intent(in)            :: contPhoton    ! is this a continuum photon?
    real(kind=doublekind), intent(in) :: offset     ! offset scale factor
    real(kind=doublekind), optional, intent(inout) :: tau_inf ! offset scale factor
    
    
    type(OCTALVECTOR) :: octVec
    integer :: subcell
    real(kind=doublekind) :: tval
    integer               :: nTau                   ! size of optical depth arrays
    type(octalvector)     :: rVec                   ! position vector
    real :: kappaScaReal, kappaAbsReal, etaCont
    real(kind=octalKind) :: chi, eta, tau
    real(kind=quadkind)  :: integral
    integer :: ilambda
    logical :: escaped
    type(OCTAL), pointer :: thisOctal
    type(OCTAL),pointer :: oldOctal
    real(kind=quadkind) :: dtau, delta
    real(kind=doublekind)::  I00                         ! input intensity
    real(kind=doublekind), parameter :: c=2.99792458d10  ! speed of light in cm/s
    
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
               kappaAbs=kappaAbsReal,kappaSca=kappaScaReal, etaCont=etaCont, grid=grid)

          nTau = nTau + 1
          
          eta = etaCont  ! [erg/s/sr/cm/cm^3]
          chi = ( kappaScaReal + kappaAbsReal )  ! [1/10^10cm]

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
               kappaAbs=kappaAbsReal,kappaSca=kappaScaReal, etaCont=etaCont, grid=grid)

          nTau = nTau + 1
          
!          eta = 0.0
          eta = etaCont  ! [erg/s/sr/cm/cm^3]
          chi = (kappaScaReal + kappaAbsReal)  !in [1/10^cm] I hope

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




end module formal_solutions

