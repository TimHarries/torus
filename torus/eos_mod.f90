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
SUBROUTINE get_eos_info(Mstar,Mdot,metallicity,x,y,z,rho,Temp,Omega,irrchoice,Q_irr,T_irr,fixatQcrit,fixatalphasat,nspiralarms)
use utils_mod, only:gasdev
use mpi
implicit none
logical,save :: firstTime = .true.
logical,save :: first1 = .true.
logical,save :: first2 = .true.
logical,save :: first3 = .true.
logical,save :: first4 = .true.
logical,save :: first5 = .true.
logical      :: firstcheck = .true.
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
real(double) :: betac                          !Beta for cooling
real(double) :: alpha                          !Gravitational alpha
real(double) :: Mstar1,Mdot1,z1                !Local copies of (IN)
real(double) :: Omega1                         !Local copies of (OUT)
real(double) :: spiralA1,spiralB1              !a and b in spiral equation
real(double) :: dSigma, dSigma_max             !Variation in surface density
real(double) :: theta_s, theta_loc, phi1       !Spiral angle, current location,
                                               !Difference between them
real(double) :: x1,y1                          !Local copies of x and y
integer :: iter

!Change
real(double),parameter :: gamma_sigma=1.0e30 !Accretion timescale(orbital periods)
real(double),parameter :: gamma_omega=1.0e30 !Omega variation timescale(in orb p)
real(double),parameter :: Qfrag = 2.0        !Fragmentation limit
real(double),parameter :: Q_irrcrit = 1.99   !Critical irradiation
real(double),parameter :: alpha_sat = 0.1   !Torque saturates
real(double),parameter :: gamma_Jcrit = -5   !Min for fragmentation
real(double),parameter :: gamma_Qcrit = 10   !abs(gamma_Q)> this = fragment
real(double),parameter :: mjup = 1.8986e30   !Jupiter mass
real,intent(in)        :: Q_irr              !Fix irradiation to Q_irr
real,intent(in)        :: T_irr              !Fix irradition to T_irr
integer,intent(in)     :: irrchoice          !Radiation choice
integer,intent(in)     :: nspiralarms        !Number of spiral arms
real                   :: nspiralarms1       !Local copy of nspiralarms
real                   :: Q_irr1,T_irr1      !Local copies
real(double)           :: cs_irr             !Irradiated sound speed
real(double)           :: Msol               !Mass in solar masses
real(double)           :: rAU                !Radius in AU
real(double)           :: TLin               !Ida Lin temperature
real(double)           :: deltasigma         !Rho perturbations
real(double)           :: mjeans,ljeans      !Jeans mass and length
real(double)           :: rhill              !Hill radius
real(double)           :: gamma_J,gamma_Q    !To determine fragmentation
real(double)           :: csterm             !For calculating
real(double)           :: gamma_sigma_max    !Max allowed for steady disc
real(double)           :: alpha_grav         !Grav alpha
real(double)           :: Q_var              !Q that varies in the disc.
logical,intent(in)     :: fixatQcrit         !Fix Q to critical value
logical,intent(in)     :: fixatalphasat      !Fix alpha to saturation
integer                :: frag,selfgrav      !Flags for self-grav and fragment
integer                :: taskid1,ierr1
real(double)           :: alphasave          !alpha check for frag check
real(double)           :: sigmasave          !sigma check for frag check

!irrchoice options
! 1 = no irradiation
! 2 = irradiation at fixed Q_irr
! 3 = T_irr (irradiation temp)
! 4 = IdaLin prescription


  !Convert to cgs, and cast things
  Mstar1       = Mstar*umass
  Mdot1        = Mdot*umass/3.15e7 
  r            = sqrt(x**2 + y**2)
  r            = r*udist
  z1           = abs(z)*udist
  x1           = x*udist
  y1           = y*udist
  Omega1       = sqrt(G*Mstar1/r**3)
  nspiralarms1 = real(nspiralarms)
!  spiralA1     = 2400.0
!  spiralB1     = -0.1103178

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
     open(564,file="ralpha.dat",status="replace")
     open(594,file="xyzsigma.dat",status="replace")
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

           cs = Qfrag*pi*G*sigma/Omega1

           !Set Q_var = Qfrag as a first guess
           Q_var = Qfrag

           !Calculate scale height
           !H = cs/Omega1

           !Scale height of self-grav disc is different.
           H = (cs*cs)/(pi*G*sigma)

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

           ! If no irradiation, then set all irradiation parameters to zero
           IF(irrchoice==1) THEN
              T_irr1 = 0.0
              Q_irr1 = 0.0
              cs_irr = 0.0
           ENDIF

           ! If irradiation at fixed Q_irr, calculate other irradiation parameters
           IF(irrchoice==2) THEN
              !Assign local copy a value
              Q_irr1 = Q_irr
              cs_irr = Q_irr1*pi*G*sigma/Omega1
              T_irr1 = 0.0
              IF(cs_irr/=0.0) THEN
                 CALL eos_cs(rhomid,cs_irr)                 
                 T_irr1 = gammamuT(3)
              ENDIF
 
              ! If irradiation uses T_irr, calculate this and other parameters
           ELSE IF(irrchoice>2) THEN
              T_irr1 = T_irr  
            
              ! Ida and Lin prescription modifies Tirr
              IF(irrchoice==4) THEN
                 rAU = r/udist
                 Msol = Mstar/umass
                 TLin = 280.0*Msol/sqrt(rAU)                 
                 ! Must account for optical depth
                 IF(tau/=0.0)TLin = TLin/(tau+1.0/tau)**0.25
                 T_irr1 = max(TLin, T_irr1)
              ENDIF

              cs_irr = SQRT(gamma*Boltzmannk*T_irr1/(mu*mH))
              Q_irr1 = cs_irr*omega1/(pi*G*sigma)

              ! If Q_irr > crit then fix it to crit
              ! If logical flag on to do so
                 IF(Q_irr1>Q_irrcrit) THEN
                    !Let Q be fixed
                    !Change to not fix Q now 23rd april
                   ! Q_irr1 = Q_irrcrit

                    cs_irr = Q_irr1*pi*G*sigma/Omega1
                    !Don's set irradiation temp to zero,
                    !Leave it at the input.
                    
                    T_irr1 = T_irr
                    !T_irr1 = 0.0
                    IF(cs_irr/=0.0) THEN
                       CALL eos_cs(rhomid,cs_irr)                 
                    !   T_irr1 = gammamuT(3)
                    ENDIF
                 ENDIF

           ENDIF

           !if((T**4.0-T_irr1**4.0).LT. 0.0)then
            !if((T-T_irr1).LT. 0.5)then 
           !Change this for small radii
           !Basically the temperature goes below but it's ok we just 
           !Do this sanity check

           !This should only happen in irradiated case
           !Looks like also only happens for low accretion rates.
            if((T-T_irr1).LT. 0.0)then   
                 !Set alpha to backgroud vale 10^-2
                 !Cass change
                 !change 17 August

                  alpha  = 1.0e-2 
!-----------------Change for fixed Beta
                  !betac = 6.0
                  !alpha    = 4.0/(9.0*gamma*(gamma-1.0)*betac) +1.0e-2
!-----------------End fixed beta
                  !This T_irr is set in the input file, it's the lowest it can go
                  !Maybe I want T alone here.
                 ! cs_irr = SQRT(gamma*Boltzmannk*T/(mu*mH))

!-----------------Changed back to T_irr 1oct 2015
                  !This way we set sound speed by whatever is highest
                  cs_irr = SQRT(gamma*Boltzmannk*T_irr1/(mu*mH))

                  !Set cs = cs_irr for further down
                 
                  !cs = cs_irr!Commented out 1oct 2015

                  sigma = Mdot1*omega1/(3.*pi*cs_irr*cs_irr*alpha)

                  !Reset scale height
                  H = (cs_irr*cs_irr)/(pi*G*sigma)                  

                  Q_var = 3.*alpha*(cs_irr*cs_irr*cs_irr)/(G*mdot1)
                  !T = T_irr
     
                  EXIT

            else 

                 betac = (tau+1.0/tau)*(cs*cs)*omega1/&
                      &(sigma_SB*(T**4.0-T_irr1**4.0)*gamma*(gamma-1.0))

!----------------Change for fixed beta
                ! betac = 6.0
!----------------End fixed beta

                 !Calculate alpha from this value --> accretion rate
                 !alpha   = alpha_grav + alpha_visc
                 alpha    = 4.0/(9.0*gamma*(gamma-1.0)*betac) +1.0e-2
                 if(alpha.LT. 1.0e-2)alpha=1.0e-2

                 mdot_try = 3.0*pi*alpha*cs*cs*sigma/omega1

                 Q_var = 3.*alpha*(cs*cs*cs)/(G*mdot_try)
           end if           

           dT = (Mdot1-Mdot_try)/Mdot1
           IF(ntries> 1000) THEN 
              fine = fine/10.0
              ntries = 0.0
           ENDIF

           dsig = (sigma-sigma_old)/sigma_old   
      
           sigma_old = sigma
           sigma     = sigma*(1.0 +dT/(abs(dT))*fine)
           ntries    = ntries+1

!----------Fixed a bug that meant couldn't correctly guess sigma at large radii
!----------If alpha large, double check correct
!----------Because this is based on a code that worked on sigma_old being
!           from a smaller radii, this isn't great at large radii.
!           in this case, want to double check that it is actually behaving

           if(alpha .GT. alpha_sat)then
              if(firstcheck)then
                 alphasave=alpha
                 sigmasave=sigma
                 sigma_old=sigma/2.0
                 sigma = sigma*(1.0 +dT/(abs(dT))*fine)
                 firstcheck = .false.

                 !Now go to end of loop and check
                 CYCLE
              else !Not first check

                 !If percentage change is small, we are good
                 if((sigmasave-sigma)/sigma <0.00001)then

                    !Check as normal
                    IF(ABS(dsig)<1.0e-5.and.ntries>500) exit 

                    !If not satisfied,try again
                    firstcheck = .true.

                    !To next END DO
                    CYCLE 

                 !else percentage change not small   
                 else

                    !Set flag to true to check so it divide again
                    firstcheck=.TRUE.
                    
                    !Go to next END DO
                    CYCLE
                 end if
              !End first check   
              end if
            !End frag check  
            end if

!-----------------------------------------------------------------------------

           !Exit if percentage change in sigma very small
           IF(ABS(dsig)<1.0e-5.and.ntries>500) exit 

        ENDDO

      
      Omega = Omega1
      Temp  = T

        IF(T>1000.0) THEN
           alpha = 0.01
           sigma = Mdot_try*omega1/(3.0*pi*alpha*cs*cs)
              print*,'MRI active ', r/udist,T, alpha, sigma
 
        ENDIF
           !If alpha exceed saturation,set sigma=0
           !As such a disc would fragment
           
           if(alpha> alpha_sat) then
              
              !Only fix it at alpha sat if flag is on
              if(fixatalphasat) then 
                 rho   = 0.0
                 Omega = 0.0
                 Temp  = 0.0
                 sigma = 0.0
                 dsigma=0.0                
               end if

           else
                !Calc gravitational alpha
                !alpha_grav = alpha_code - background
                !So dsigma becomes zero if alpha_grav 0
                alpha_grav = alpha - 1.0e-2

  		!Theta for x and y
  		if(x1 .GE. 0.0)then
     			theta_loc = asin(y1/r)
  		else
	    	 	theta_loc = -1.*asin(y1/r) + 3.14159265359
 		end if

  		!Spiral stuff
                !The constants now set in params filie
		!For loose spirals A=20,B=1. For tight,A=2400,B=-0.1103178
  		spiralA1   = 13.5
  		spiralB1   = 0.3832
  		r          = r/udist
 		theta_s    = (1./spiralB1)*log(r/spiralA1)
  		r          = r*udist
  		phi1       = theta_s - theta_loc

  		dSigma_max = sqrt(sqrt(alpha_grav)**2.)*sigma
  		dSigma     = -1.*dSigma_max*cos(nspiralarms1*(phi1))

                !Now add some noise to this dSigma.
                !Gasdev returns number with mean=0, variance(and SD)=1
                !So multiply by the desired SD, and add the desired mean.
                !This gives us a gaussian centered on dsigma, with SD=15%
                !Commented out noise, not sure need it.
                !dsigma     = dsigma+(0.15*dsigma)*gasdev()
  		sigma      = sigma + dSigma
  		rhomid     = sigma/(2.0*H)
  		arg        = z1/H
  		rho        = rhomid/(cosh(arg)*cosh(arg))
  		cs         = Qfrag*pi*G*sigma/Omega1
                !Changed to depends on Q_var 19th August
                !cs         = Q_var*pi*G*sigma/Omega1
  		!Use EoS to calculate tau, gamma, betac
  		CALL eos_cs(rho, cs)
                !Watch this, I don't know if this need to be here.
                !changed on 19th August
                !Leave this in for T_irr as it catches local variations
  		Temp          = gammamuT(3)

                !If temp is below irradiated temp, make it equal Tirr
                if(Temp - T_irr1 .LT. 0.0) Temp=T_irr1

                 !print*,"alpha_grav",alpha_grav,r/udist,sigma,dsigma
        end if

        call MPI_COMM_RANK(MPI_COMM_WORLD, taskid1, ierr1)
        !Write to ralpha.dat
        if(taskid1 .EQ. 0) write(564,*)r/udist,mdot1*3.15e7/umass,alpha,sigma,dsigma,betac,Temp,cs,T,rho,tau,Q_var
        !Write to xyzsigma.dat 
	if(taskid1 .EQ.0) write(594,*)x1,y1,z1,sigma,rho,dSigma,Temp,tau      

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
