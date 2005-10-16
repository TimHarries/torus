module spectrum_mod

  use constants_mod
  use atom_mod
  use utils_mod

  implicit none

  public

  type SPECTRUMTYPE
     real(double), pointer :: flux(:)
     real(double), pointer :: normflux(:)
     real(double), pointer :: normflux2(:)
     real(double), pointer :: lambda(:)
     real(double), pointer :: prob(:)
     real(double), pointer :: dlambda(:)
     integer :: nLambda
  end type SPECTRUMTYPE
  

  contains

    subroutine getWavelength(spectrum, wavelength)

      type(SPECTRUMTYPE) :: spectrum
      real(double) :: wavelength
      real(double) :: r, t
      integer :: i

      call random_number(r)
      call locate(spectrum%prob, spectrum%nLambda, r, i)
      t = (r - spectrum%prob(i))/(spectrum%prob(i+1)-spectrum%prob(i))
      wavelength = spectrum%lambda(i) + t*(spectrum%lambda(i+1)-spectrum%lambda(i))

    end subroutine getWavelength


    subroutine fillSpectrumBB(spectrum, teff, lamStart, lamEnd, nLambda, biasToLyman)

      type(SPECTRUMTYPE) :: spectrum
      integer :: nLambda
      real(double) :: lamStart, lamEnd, teff
      logical, optional :: biasToLyman
      real(double) :: logLamStart, logLamEnd
      integer :: i

      allocate(spectrum%flux(1:nLambda))
      allocate(spectrum%lambda(1:nLambda))
      allocate(spectrum%dlambda(1:nLambda))
      allocate(spectrum%prob(1:nLambda))

      logLamStart = log10(lamStart)
      logLamEnd = log10(lamEnd)
      
      do i = 1, nLambda
         spectrum%lambda(i) = logLamStart + real(i-1)/real(nLambda-1)*(logLamEnd - logLamStart)
         spectrum%lambda(i) = 10.**spectrum%lambda(i)
      enddo
      do i = 2, nLambda-1
         spectrum%dlambda(i) = 0.5*((spectrum%lambda(i+1)+spectrum%lambda(i))-(spectrum%lambda(i)+spectrum%lambda(i-1)))
      enddo
      spectrum%dlambda(1) = spectrum%lambda(2)-spectrum%lambda(1)
      spectrum%dlambda(nLambda) = spectrum%lambda(nlambda)-spectrum%lambda(nLambda-1)
      
      do i = 1, nLambda
         spectrum%flux(i) = bLambda(spectrum%lambda(i), dble(teff))
      enddo
      spectrum%nLambda = nLambda
      where(spectrum%flux(1:spectrum%nLambda) == 0.d0) spectrum%flux = 1.e-30

      call probSpectrum(spectrum, biasToLyman)
    end subroutine fillSpectrumBB

    subroutine readSpectrum(spectrum, filename)
      type(SPECTRUMTYPE) :: spectrum
      character(len=*) :: filename
      real :: fTemp(1200000),xTemp(1200000)
      integer :: nLambda, i

      open(20,file=filename,form="formatted",status="old")
      nLambda = 1
10    continue
      read(20,*,end=20) xtemp(nLambda),fTemp(nLambda)
      nLambda = nLambda + 1
      goto 10
20    continue
      close(20)
      nLambda = nLambda - 1
      allocate(spectrum%flux(1:nLambda))
      allocate(spectrum%lambda(1:nLambda))
      allocate(spectrum%dlambda(1:nLambda))
      allocate(spectrum%prob(1:nLambda))
      spectrum%nLambda = nLambda
      spectrum%flux(1:nLambda) = fTemp(1:nLambda)
      spectrum%lambda(1:nLambda) = xTemp(1:nLambda)
      do i = 2, nLambda-1
         spectrum%dlambda(i) = 0.5*((xtemp(i+1)+xtemp(i))-(xtemp(i)+xtemp(i-1)))
      enddo
      spectrum%dlambda(1) = xtemp(2)-xtemp(1)
      spectrum%dlambda(nLambda) = xtemp(nlambda)-xtemp(nLambda-1)
      call probSpectrum(spectrum)
    end subroutine readSpectrum

    subroutine probSpectrum(spectrum, biasToLyman)
      type(SPECTRUMTYPE) :: spectrum
      logical, optional :: biasToLyman
      real(double) :: fac
      integer :: i
      spectrum%prob = 0.d0
      do i = 2, spectrum%nLambda
         fac = 1.d0
         if (present(biasToLyman)) then
            if (biasToLyman) then
               if (spectrum%lambda(i) < 912.) then
                  fac = 10.d0
               else
                  fac = 1.d0
               endif
            endif
         endif
         spectrum%prob(i) = spectrum%prob(i-1) + spectrum%flux(i) * spectrum%dLambda(i) * fac
      enddo
      spectrum%prob(1:spectrum%nLambda) = spectrum%prob(1:spectrum%nLambda) / spectrum%prob(spectrum%nLambda)
    end subroutine probSpectrum

    subroutine normalizedSpectrum(spectrum)
      type(SPECTRUMTYPE) :: spectrum
      allocate(spectrum%normflux(1:spectrum%nLambda))
      spectrum%normflux(1:spectrum%nLambda) = spectrum%flux(1:spectrum%nLambda) &
           / SUM(spectrum%flux(1:spectrum%nLambda)*spectrum%dlambda(1:spectrum%nLambda))
    end subroutine normalizedSpectrum

    real(double) function returnNormValue(spectrum, lambda)
      type(SPECTRUMTYPE) :: spectrum
      real(double) :: lambda, t
      logical :: ok
      integer :: i

      if ((lambda < spectrum%lambda(1)).or.(lambda > spectrum%lambda(spectrum%nlambda))) then
         returnNormValue = 0. 
      else
!         i = findIlambda(real(lambda), real(spectrum%lambda), spectrum%Nlambda, ok)
         call locate(spectrum%lambda,spectrum%nLambda, lambda, i)
         t = (lambda - spectrum%lambda(i))/(spectrum%lambda(i+1) - spectrum%lambda(i))
!         t = 0.
         returnNormValue = spectrum%normFlux(i) + t * (spectrum%normFlux(i+1) - spectrum%normFlux(i))
!         returnNormValue = interpLogLinearDouble(spectrum%lambda, spectrum%normFlux, spectrum%nLambda, lambda)
      endif
    end function returnNormValue
         
  end module spectrum_mod
