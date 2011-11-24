  !
  ! Compares the Lorentz profile and Voigt Profile
  !
  subroutine test_profiles()
    use atom_mod, only: bigGamma
    implicit none
    
    integer, parameter :: nbin = 100
    real(double) :: Lorentz(nbin), lambda(nbin), freq(nbin), Voigt_prof(nbin)
    integer :: i, j  
    integer :: nsample = 10000000
    real(double) :: nu0, nu_rand, lam0, lam_min, lam_max, Gamma
    real(double) :: lam_rand
    real(double) :: c=2.99792458e10  ! cm/s
    real(double) :: tmp, dlam, a, doppler, dnu, dv, T, N_HI, Ne
    
    lam_min = 6550.0d-8  ! cm
    lam0    = 6562.8d-8  ! cm
    lam_max = 6580.0d-8  ! cm
    
    nu0 = c/lam0
    
    ! sets up wavelength array
    tmp = (lam_max-lam_min)/dble(nbin-1)
    do i = 1, nbin
       lambda(i) = lam_min + tmp*dble(i-1)
       freq(i) = c/lambda(i)
    end do
    
     
    
    ! take random samples
!    Gamma = 3.48d11
    T = 6000.0; N_HI=1.e-19; Ne=1.e-9
    Gamma = bigGamma(N_HI, T, Ne, nu0)
    !Gamma = 0.0d0

    
    dlam = lambda(2) - lambda(1)
    
    Lorentz(1:nbin) = 0.0d0
    do i = 1, nsample
       
       nu_rand = random_Lorentzian_frequency(nu0, Gamma)
       
       lam_rand = c/nu_rand + dlam/2.0
       
       call locate(lambda, nbin, lam_rand, j)
       
       Lorentz(j) = Lorentz(j) + 1.0d0
       
    end do
    
    
    ! No computes the Voigt profile
    doppler = nu0/cSpeed * sqrt(2.*kErg*T/mHydrogen)
    a = Gamma/4.0d0/3.141593d0/doppler
    do i = 1, nbin
       dnu = nu0 - freq(i)
       dv = dnu/doppler
       voigt_prof(i) = voigtn(a, dv)
    end do
    
  
    ! Renomarlize the profiles so that the max is 1
    tmp = MAXVAL(Lorentz(1:nbin)) 
    if (tmp/=0) Lorentz(:) = Lorentz(:)/tmp

    tmp = MAXVAL(voigt_prof(1:nbin)) 
    if (tmp/=0) voigt_prof(:) = voigt_prof(:)/tmp
  
    
    ! writing the results in a file
    
    open(unit=66, file="lorentz_voigt.dat", status="replace")
    do i = 1, nbin
       write(66,*) lambda(i)*1.0d8, lorentz(i), voigt_prof(i)
       !             [a]              
    end do
    close(66)
    
  end subroutine test_profiles
