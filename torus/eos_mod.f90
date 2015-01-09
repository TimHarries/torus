module eos_mod



!-----------------------------------------------------------------------
! Data module for saving equation of state values
! PJC 22/05/2008
! DHF 09/09/2009
!-----------------------------------------------------------------------
      use kind_mod
      use messages_mod
      implicit none
      private
      public :: get_eos_info
! Units
      integer :: nrhopoints,nUpoints
      real(double),parameter :: Boltzmannk = 1.3807d-16
      real(double),parameter :: mH = 1.6726d-24
      real(double),parameter :: G = 6.672041d-8
      real(double),parameter :: udist = 1.50d13
      real(double),parameter :: umass = 1.99d33
      real(double),parameter :: utime = sqrt((udist**3)/(G*umass))
      real(double),parameter :: uergg = udist*udist/(utime*utime)

! Arrays
      real(double),allocatable,dimension(:) :: gammamuT
      real(double),allocatable,dimension(:,:,:) :: eostable
      real(double),allocatable,dimension(:,:) :: cstab


contains
  subroutine eosread
    !-----------------------------------------------------------------------
    ! Reads the equation of state tables
    ! Updated to fortran 90 PJC 22/05/2008
    ! Extended to opacities DHF 01/09/2009
    !
    ! Note the columns (== k values) are as follows:
    !   1) density
    !   2) temperature
    !   3) internal energy
    !   4) mean molecular weight
    !   5) Mass Weighted Opacity
    !   6) Mass-weighted Opacity
    !
    !-----------------------------------------------------------------------

    use unix_mod
    implicit none
    integer :: i,j,k,check
    real(double) :: cv, gamma
    character(len=80) :: fname, dataDirectory
    ! Open eos data file and read in values

    call unixGetenv("TORUS_DATA", dataDirectory, i)
    fname = trim(dataDirectory)//"/"//"myeos.dat"

    open(50,file=fname, &
         status='old',iostat=check,action='read')
    if (check /= 0) then
       print*, 'Input file myeos.dat not found'
       stop
    endif
    if (Writeoutput) print*, "Reading equation of state tables"

    read(50,*) nrhopoints, nUpoints

    ! Allocate input data array
    allocate(eostable(nrhopoints,nUpoints,6))
    allocate(cstab(nrhopoints,nUpoints))
    do i = 1,nrhopoints
       do j = 1,nUpoints
          read(50,*) (eostable(i,j,k),k=1,6)

          !     Now calculate sound speed values for interpolation
          cv = eostable(i,j,3)/eostable(i,j,2)
          gamma = 1.0d0 + (Boltzmannk / &
               (mH*eostable(i,j,4)*cv))
          cstab(i,j) = sqrt(gamma*(gamma-1.0)*eostable(i,j,3))
       enddo
    enddo
    close(50)

    if (writeoutput) then
       print*, "Equation of state tables read in successfully"
       print*, "-----------------------------------------------"
    endif

  end subroutine eosread

! get_eos_info
! Returns midplane density, omega and temperature for a given accretion rate
! and radius
! Cass Hall Nov 2014
!=========================================================================
SUBROUTINE get_eos_info(Mstar,Mdot,metallicity,x,y,z,rho,Temp,Omega)
implicit none
logical,save :: firstTime = .true.
real(double), parameter  :: tolerance = 1.0e-5 !Tolerance limit
real(double), parameter  :: Q  = 1.9           !Toomre paramter, fixed
real(double), parameter  :: pi = 3.14159265285 !Pi
real(double), parameter  :: sigma_SB = 5.67e-5 !Stefan Boltzmann
real(double),intent(in)  :: Mstar              !Central star mass 
real(double),intent(in)  :: Mdot               !Accretion rate 
real(double),intent(in)  :: x,y,z              !Position
real(double),intent(in)  :: metallicity        !Relative to solar
real(double),intent(out) :: Omega              !Angular frequency
real(double),intent(out) :: rho                !Density 
real(double),intent(out) :: Temp               !Temperature

real(double) :: sigma_old                      !Old surface density
real(double) :: sigma                          !Surface density
real(double) :: dT,fine,ntries,dsig            !For tolerance
real(double) :: Mdot_try                       !Compare to imposed Mdot
real(double) :: cs                             !Sound speed
real(double) :: H                              !Scale height
real(double) :: r                              !Radius
real(double) :: rhomid                         !Midplane density
real(double) :: arg                            !Junk argument
real(double) :: gamma,mu,T,kappa,tau           !Info from equation of state
real(double) :: T_irr                          !Irradiated temp
real(double) :: betac                          !Beta for cooling
real(double) :: alpha                          !Gravitational alpha
real(double) :: Mstar1,Mdot1,z1                !Local copies of (IN)
real(double) :: Omega1                         !Local copies of (OUT)
real(double) :: spiralA,spiralB                !a and b in spiral equation
real(double) :: dSigma, dSigma_max             !Variation in surface density
real(double) :: theta_s, theta_loc, phi1       !Spiral angle, current location,
                                               !Difference between them
real(double) :: x1,y1                          !Local copies of x and y
integer :: iter
!------------------------------------------------------------

  !Convert to cgs
  Mstar1 = Mstar*umass
  Mdot1  = Mdot*umass/3.15e7 
  r      = sqrt(x**2 + y**2)
  r      = r*udist
  z1     = abs(z)*udist
  x1     = x*udist
  y1     = y*udist
  Omega1 = sqrt(G*Mstar1/r**3)

  ! gammamuT columns:
  ! 1) gamma
  ! 2) mu
  ! 3) T
  ! 4) Opacity
  ! 5) Mass Weighted Opacity (Stamatellos et al 2007)

  !Read in equation of state tables.
  if (firstTime) then
     allocate(gammamuT(5))
     CALL eosread
     firstTime = .false.
  endif

  !Guess Sigma 

  sigma_old = 50000.0
  sigma     = 2.0*sigma_old
  dT        = 1.0e30
  ntries    = 0.0
  fine      = 0.01
  iter = 0


  DO WHILE(ABS(dT)> tolerance)
           if (iter > 20000) then
              write(*,*) "exiting after 20000 iterations ",abs(dt)
              exit
           endif
           iter = iter + 1
           !Calculate sound speed assuming fixed Q


           cs = Q*pi*G*sigma/Omega1

           !Calculate scale height
           H = cs/Omega1

           !Midplane density
           rhomid = sigma/(2.0*H)
           arg    = z1/H
           rho    = rhomid/(cosh(arg)*cosh(arg))

           !Use EoS to calculate tau, gamma, betac
           CALL eos_cs(rhomid, cs)

           gamma = gammamuT(1)
           mu    = gammamuT(2)
           T     = gammamuT(3)
           kappa = (gammamuT(4)*metallicity)
           tau   = sigma*kappa

           !Assume no irraditation
           T_irr = 0.0

           ! Calculate cooling timescale for these parameters
           betac = (tau+1.0/tau)*(cs*cs)*Omega1/&
                (sigma_SB*(T**4.0-T_irr**4.0)*gamma*(gamma-1.0))

           !Calculate alpha from this value --> accretion rate
           alpha = 4.0/(9.0*gamma*(gamma-1.0)*betac)

           !Compare calculated mdot with imposed mdot
           Mdot_try = 3.0*pi*alpha*cs*cs*sigma/Omega1
           dT = (Mdot1-Mdot_try)/Mdot1

           IF(ntries> 1000) THEN 
              fine = fine/10.0
              ntries = 0.0
           ENDIF

           dsig = (sigma-sigma_old)/sigma_old   
      
           !Exit if percentage change in sigma very small
           IF(ABS(dsig)<1.0e-5.and.ntries>500) exit  

           sigma_old = sigma
           sigma     = sigma*(1.0 +dT/(abs(dT))*fine)
           ntries    = ntries+1

        ENDDO
      
      Omega = Omega1
      Temp  = T
	
  !write(564,*)x,y,z
  !xcart=x*cos(y)
  !ycart=x*sin(y)	
  !write(564,*)xcart,ycart


  !Theta for x and y
  if(x1 .GE. 0.0)then
     theta_loc = asin(y1/r)
  else
     theta_loc = -1.*asin(y1/r) + 3.14159265359
  end if

  !Spiral stuff
  spiralA    = 20.0
  spiralB    = 1.0
  r          = r/udist
  theta_s    = (1./spiralB)*log(r/spiralA)
  r          = r*udist
  phi1        = theta_s - theta_loc

  dSigma_max = sqrt(sqrt(alpha)**2.)*sigma
  dSigma     = -1.*dSigma_max*cos(2.*(phi1))
  sigma      = sigma + dSigma
  rhomid     = sigma/(2.0*H)
  arg        = z1/H
  rho        = rhomid/(cosh(arg)*cosh(arg))
  cs         = Q*pi*G*sigma/Omega1
  !Use EoS to calculate tau, gamma, betac
  CALL eos_cs(rho, cs)
  Temp          = gammamuT(3)


END SUBROUTINE get_eos_info

      subroutine eos_cs(rho,cs)
!-----------------------------------------------------------------------
! Reads and interpolates the equation of state tables to give
! temperatures, opacities etc
! Updated to fortran 90 PJC 22/05/2008
! Extended to opacities DHF 01/09/2009
!-----------------------------------------------------------------------
      

      implicit none

      integer :: i,j
      real(double) :: mT,mT1,mT2,cT,cT1,cT2,T1,T2
      real(double) :: mmu,mmu1,mmu2,cmu,cmu1,cmu2,mu1,mu2
      real(double) :: mkap,mkap1,mkap2,ckap,ckap1,ckap2,kap1,kap2 
      real(double) :: mkbar,mkbar1,mkbar2,ckbar,ckbar1,ckbar2,kbar1,kbar2
      real(double) :: rho, cs
                              
! gammamuT columns:
! 1) gamma
! 2) mu
! 3) T
! 4) Opacity
! 5) Mass Weighted Opacity (Stamatellos et al 2007)
      gammamuT = 0.0
!print*,"in eos"
! Find the relevant records in the table...
! ... for rho
         if (rho < 1.0e-24) rho = 1.0e-24
         i = 1
         do 
            if ((eostable(i,1,1) >= rho).or.(i == nrhopoints)) exit
            i = i + 1
         enddo
                        
! ... and for internal energy
         if (cs < cstab(1,1)) cs = cstab(1,1)
         j = 1
         do 
            if ((cstab(i-1,j) >= cs).or.(j==nUpoints)) exit
            j = j + 1
         enddo

        IF(j==1) j = j+1
                        
! Interpolate over the j value at i-1
         mT1 = (eostable(i-1,j-1,2) - eostable(i-1,j,2))/ &
              (cstab(i-1,j-1) - cstab(i-1,j))
         cT1 = eostable(i-1,j,2) - mT1*cstab(i-1,j)
         T1 = mT1*cs + cT1

         mmu1 = (eostable(i-1,j-1,4) - eostable(i-1,j,4))/ &
              (cstab(i-1,j-1) - cstab(i-1,j))
         cmu1 = eostable(i-1,j,4) - mmu1*cstab(i-1,j)
         mu1 = mmu1*cs + cmu1
                 
         mkap1 = (eostable(i-1,j-1,6) - eostable(i-1,j,6))/ &
              (cstab(i-1,j-1) - cstab(i-1,j))
         ckap1 = eostable(i-1,j,6) - mkap1*cstab(i-1,j)
         kap1 = mkap1*cs + ckap1
                 
         mkbar1 = (eostable(i-1,j-1,5) - eostable(i-1,j,5))/ &
              (cstab(i-1,j-1) - cstab(i-1,j))
        ckbar1 = eostable(i-1,j,5) - mkbar1*cstab(i-1,j)
        kbar1 = mkbar1*cs + ckbar1
                                 
! Then interpolate over the j value at i
! Update j value as necessary
         j = 1
         do 
            if ((cstab(i,j) >= cs).or.(j==nUpoints)) exit
            j = j + 1
         enddo

        IF(j==1) j = j+1
                        
         mT2 = (eostable(i,j-1,2) - eostable(i,j,2))/ &
              (cstab(i,j-1) - cstab(i,j))
         cT2 = eostable(i,j,2) - mT2*cstab(i,j)
         T2 = mT2*cs + cT2 
         
         mmu2 = (eostable(i,j-1,4) - eostable(i,j,4))/ &
              (cstab(i,j-1) - cstab(i,j))
         cmu2 = eostable(i,j,4) - mmu2*cstab(i,j)
         mu2 = mmu2*cs + cmu2
                 
         mkap2 = (eostable(i,j-1,6) - eostable(i,j,6))/ &
              (cstab(i,j-1) - cstab(i,j))
         ckap2 = eostable(i,j,6) - mkap2*cstab(i,j)
         kap2 = mkap2*cs + ckap2
                 
         mkbar2 = (eostable(i,j-1,5) - eostable(i,j,5))/ &
              (cstab(i,j-1) - cstab(i,j))
         ckbar2 = eostable(i,j,5) - mkbar2*cstab(i,j)
         kbar2 = mkbar2*cs + ckbar2

! Finally interpolate over i at the fractional j value
         mT = (T2 - T1) / (eostable(i,1,1)-eostable(i-1,1,1))
         cT = T2 - mT*eostable(i,1,1)
         gammamuT(3) = mT*rho + cT

         mmu = (mu2 - mu1) / (eostable(i,1,1)-eostable(i-1,1,1))
         cmu = mu2 - mmu*eostable(i,1,1)
         gammamuT(2) = mmu*rho + cmu
                 
        mkap = (kap2 - kap1) / (eostable(i,1,1)-eostable(i-1,1,1))
         ckap = kap2 - mkap*eostable(i,1,1)
         gammamuT(4) = mkap*rho + ckap
                 
        mkbar = (kbar2 - kbar1) / (eostable(i,1,1)-eostable(i-1,1,1))
         ckbar = kbar2 - mkbar*eostable(i,1,1)
         gammamuT(5) = mkbar*rho + ckbar

! Evaluate the gamma value from gamma = 1 + k/(mH*mu*cv), with cv = U/T
         
         gammamuT(1) = cs*cs*gammamuT(2)*mH/(Boltzmannk*gammamuT(3))
     

      end subroutine eos_cs

end module eos_mod

