module density_mod

  !
  ! Generic function to return a density for a given point
  ! and geometry type.
  !
  ! If you add a new geometry, please don't forget to add the name
  ! in the list in 
  !
  
  use utils_mod
  use vector_mod
  use gridtype_mod
  use jets_mod
  use wr104_mod
  use ostar_mod

  public :: density, TTauriDensity,spiralWindDensity
  ! the specific definition of density functions should really be private, 
  ! but for now they are public... Better yet, they should be in their own module.
  ! See jets_mod for example. 
  
  private :: print_geometry_list
  
contains

  !
  !  For a given point (as a vector) with each component in
  !  10^10 cm and gridtype object, it will return the density in g/cm^3
  function density(r_vec, grid) RESULT(out)
    implicit none
    double precision :: out
    type(octalVector), intent(in) :: r_vec
    type(gridtype), intent(in) :: grid

    select case (grid%geometry)
       
    case ("jets")
       ! using a routine in jets_mod.f90
       out = JetsDensity(r_vec, grid)/1000.d0 
       !                              ^^^^^^^
       !                          Converting kg/m^3 to g/cm^3

!       ! test function....
!       out = density_inverse_sq(r_vec)

    case ("ttauri")
       !  using a routine this module
       out = TTauriDensity(r_vec, grid) ! [g/cm^3]
       ! Double check the units

    case("testamr")
       out = testDensity(r_vec,grid)

    case("proto")
       out = protoDensity(r_vec,grid)

    case("spiralwind")
       out = spiralWindDensity(r_vec, grid)

    case("benchmark")
       out = benchmarkDensity(r_vec, grid)

    case("melvin")
       out = melvinDensity(r_vec, grid)

    case("shakara","aksco")
       out = shakaraSunyaevDisc(r_vec, grid)

    case("clumpydisc")
       out = clumpydisc(r_vec, grid)
       
    case default
       print *, 'Error:: Geometry option passed to [density_mod::density] '
       print *, '       (through grid) has not been implemented yet. Sorry.'
       print *, 'Your geometry = ', grid%geometry
       print *, 'Avilable gemetries are: '
       call print_geometry_list()
       print *, 'Terminating the program ...'
       stop 
    end select

  end function density
    


  !
  ! List the available geometry
  ! 
  subroutine print_geometry_list()    
    implicit none
    ! # of option available currently.
    integer, parameter :: n = 4
    character(LEN=30) :: name(n)
    integer :: i

    ! List here
    name(1) =  'jets'
    name(2) =  'wr104'
    name(3) =  'ttauri'
    name(4) =  'testamr'
    
    do i = 1, n
       write(*, *)  '   ', i, '. ', name(i)
    end do
       
  end subroutine print_geometry_list


  !
  !
  !
    FUNCTION TTauriDensity(point,grid) RESULT(rho)
    ! calculates the density at a given point for a model of a T Tauri 
    !   star with magnetospheric accretion
    !   see Hartman, Hewett & Calvet 1994ApJ...426..669H 

    use parameters_mod

    IMPLICIT NONE

    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid
    REAL                          :: rho

    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec

    REAL :: r, rM, theta, y

    starPosn = grid%starPos1
    pointVec = (point - starPosn) * 1.e10_oc
    r = modulus( pointVec ) 

    ! test if the point lies within the star
    IF ( r < TTauriRstar ) THEN
      rho = 1.e-25
      RETURN
    END IF
    
    ! test if the point lies too close to the disk
    IF ( ABS(pointVec%z) < 4.0e-2 * TTauriRstar) THEN
      ! we will fake the coordinates so that the density
      !  remains constant within this region 
      pointVec%z = SIGN( 4.0e-2 * TTauriRstar, REAL(pointVec%z) )
      r = modulus( pointVec ) 
    END IF
   
    theta = ACOS( pointVec%z  / r )
    IF (ABS(MODULO(theta,pi)) > 1.e-10 ) THEN 
      rM  = r / SIN(theta)**2
    ELSE
      rM = HUGE(rM)
    END IF
     
    ! test if the point lies outside the accretion stream
    IF  ((rM > TTauriRinner) .AND. (rM < TTauriRouter )) THEN
      y = SIN(theta)**2 
  
      rho = (TTauriMdot * TTauriRstar) / (4.0 * pi * &
              (TTauriRStar/TTauriRinner - TTauriRstar/TTauriRouter))
      rho = rho * r**(-5.0/2.0) / SQRT( 2.0 * bigG * TTauriMstar ) 
      rho = rho * SQRT( 4.0 - 3.0*y) / SQRT( 1.0 - y) 
    ELSE
      rho = 1.e-25
    END IF
    
  END FUNCTION TTauriDensity


  function testDensity(point, grid)
    use input_variables
    real :: testDensity
    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid
    real :: r
    r = modulus(point)
    testDensity = 1.e-30
    if ((r > grid%rInner).and.(r < grid%rOuter)) then
       testDensity = rho * (grid%rInner / r)**2
    endif
  end function testDensity

  function protoDensity(point, grid)
    use input_variables
    real :: testDensity
    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid
    real :: r
    r = modulus(point)
    testDensity = 1.e-30
    if ((r > grid%rInner).and.(r < grid%rOuter)) then
       testDensity = grid%densityScaleFac*rho * (grid%rInner / r)**rPower
    endif
  end function protoDensity

  function benchmarkDensity(point, grid)

    use input_variables
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(octalVector), INTENT(IN) :: point
    real :: r, hr, rd
    
    rd = rOuter / 2.
    benchmarkDensity = 1.e-30
    r = sqrt(point%x**2 + point%y**2)
    if ((r > rInner).and.(r < rOuter)) then
       hr = height * (r/rd)**1.125
       benchmarkDensity = rho * ((r / rd)**(-1.))*exp(-pi/4.*(point%z/hr)**2)
    endif
!    if ((r < rInner).and.(abs(point%z) < rInner)) then
!       r = rInner
!       hr = height * (r/rd)**1.125
!       benchmarkDensity = rho * ((r / rd)**(-1.))*exp(-pi/4.*(point%z/hr)**2)*100.
!    endif
  end function benchmarkDensity

  function shakaraSunyaevDisc(point, grid) result (rhoOut)
    use input_variables
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(octalVector), INTENT(IN) :: point
    real :: r, h, rhoOut, warpHeight,alphaDisc,betaDisc
    real :: kspiral
    real :: xpoint,ypoint,rscale,r1,fac
    integer :: nspiral1
    real :: phase(10)

    nSpiral1 = 3
    do i = 1, nspiral1
       phase(i)=twoPi*real(i-1)/real(nSpiral1)
    enddo

    kspiral = grid%rOuter/Pi

    alphaDisc = 2.25
    betaDisc = 1.25

    rhoOut = 1.e-30
    r = sqrt(point%x**2 + point%y**2)
    phi = atan2(point%y,point%x)
    warpHeight = 0. !cos(phi) * rInner * sin(30.*degtorad) * sqrt(rinner / r)
    if ((r > rinner).and.(r < rOuter)) then
       h = height * (r / (100.*autocm/1.e10))**betaDisc
       rhoOut = rho0 * (rInner/r)**alphaDisc * exp(-0.5 * ((point%z-warpheight)/h)**2)
    endif
    rhoOut = max(rhoOut,1.e-30)
    xpoint = point%x
    ypoint = point%y
    fac = 1.
    do i = 1, nSpiral1
       r1 = spiraldist(xpoint,ypoint,kspiral, phase(i))
       rScale = r/10.
       fac = max(fac,1.+99.*exp(-r1/rscale)) 
    enddo
    rhoOut = rhoOut*fac

  end function shakaraSunyaevDisc

  real(kind=doubleKind) function melvinDensity(point, grid) 

    use input_variables
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(octalVector), INTENT(IN) :: point
    real :: r,  mStar, rCav, r_c, r_d
    real :: openingAngle, mu, mu_0, bigR
    real :: z, H_R, H_0, alpha
    real :: rhoEnv, rhoDisc
    real :: cavZ

! Envelope with cavity plus alpha disc model, as presented
! by Alvarez, Hoare and Lucas (2004).


    
    rStar = 10. * rSol
    H_0 = 10. * autocm
    mStar = 10. * mSol
    mDot = 1.11e-4 * mSol / (365.25*24.*3600.) * dusttogas
    rCav = 50.*auTocm
    r_c = 50. * auTocm
    r_d = 1000. * auTocm
    openingAngle = 20. * degtoRad
    alpha = 1.875
    beta  = 1.125
    rho0 = 2.e-7


    melvinDensity = 1.e-30
    r = modulus(point) * 1.e10
    mu = point%z * 1.e10 / r
    z = point%z * 1.e10

    bigR = sqrt(point%x**2 + point%y**2)*1.e10
    H_R = H_0 * (bigR / (100.*auToCm))**beta

    rhoEnv = 1.d-30
    if (r > r_c) then
       mu_0 = rtnewt(-0.2 , 1.5 , 1.e-4, r/r_c, abs(mu))
       rhoEnv = (Mdot/(8. * pi * r_c * sqrt(bigG * Mstar)))  ! Equation 1
       rhoEnv = rhoEnv * 1./sqrt(1.+mu_0) * 1./sqrt(r)
    endif
    rhoDisc = 1.d-30
    if ((bigR < r_d).and.(bigR > 10.*rStar)) then
       rhoDisc = rho0 * ((bigR/Rstar)**(-alpha))*exp(-(z**2 / (2.*H_R**2))) ! Eq 3
    endif

    melvinDensity = max(rhoDisc, rhoEnv)

    if (bigR < Rcav) melvinDensity = 1.e-30

    if (bigR > Rcav) then
       cavZ = tan(pi/2.-0.5*openingAngle)*(bigR - Rcav) ! Eq 4
       if (abs(z) > cavZ) then
          melvinDensity = 1.e-30
       endif
    endif


  end function melvinDensity

  real function rtnewt(x1,x2,xacc, p1, p2) result(junk)

    real :: x1, x2, xacc, p1, p2
    integer :: jmax, j
    real ::  dx, f, df
    parameter (jmax=20)
    junk = 0.5 * (x1+x2)
    do j=1,jmax
       call equation2(junk,f,df,p1,p2)
       dx=f/df
       junk=junk-dx
       if((x1-junk)*(junk-x2).lt.0.) then
          write(*,*) 'RTNEWT: jumped out of brackets',p1,p2,junk
          stop
       endif
       if(abs(dx).lt.xacc) return
    enddo
    write(*,*) 'rtnewt exceeding maximum iterations'
  end function rtnewt
  

  real function Equation2(mu0, eq2, deq2, r, mu)
    real :: r, mu, mu0
    real :: eq2, deq2

    eq2 = mu0**3 + (r-1.)*mu0 -r*mu
    deq2 = 3.*mu0**2 + r - 1.

  end function Equation2

  function clumpyDisc(rVec, grid) result (rhoOut)
    type(OCTALVECTOR) :: rVec
    type(GRIDTYPE) :: grid
    real :: rhoOut
    
    call findFactor(rhoOut, rVec, grid%gArray, grid%ng)
    rhoOut = rhoOut *   grid%densityScalefac
    if (modulus(rVec) < grid%rInner) rhoOut = 1.e-33

  end function clumpyDisc



end module density_mod 
