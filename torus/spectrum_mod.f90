module spectrum_mod

  use constants_mod
  use messages_mod
  use mpi_global_mod
  use utils_mod, only: locate
  use unix_mod, only: unixGetenv
  implicit none
  public

  type SPECTRUMTYPE
     real(double), pointer :: flux(:) => null()
     real(double), pointer :: normflux(:) => null()
     real(double), pointer :: normflux2(:) => null()
     real(double), pointer :: lambda(:) => null()
     real(double), pointer :: prob(:) => null()
     real(double), pointer :: dlambda(:) => null()
     integer :: nLambda
     real(double), pointer :: ppw(:) => null() !Thaw - photon packet weight
  end type SPECTRUMTYPE
  

  contains

    subroutine freeSpectrum(spectrum)
      type(SPECTRUMTYPE) :: spectrum
      if (associated(spectrum%flux)) deallocate(spectrum%flux)
      if (associated(spectrum%lambda)) deallocate(spectrum%lambda)
      if (associated(spectrum%dlambda)) deallocate(spectrum%dlambda)
      if (associated(spectrum%prob)) deallocate(spectrum%prob)
      if (associated(spectrum%normflux)) deallocate(spectrum%normflux)
      if (associated(spectrum%normflux2)) deallocate(spectrum%normflux2)
      if (associated(spectrum%ppw)) deallocate(spectrum%ppw)
      spectrum%nlambda = 0
    end subroutine freeSpectrum

    subroutine getWavelength(spectrum, wavelength, photonPacketWeight)
      use random_mod, only: randomNumberGenerator
      type(SPECTRUMTYPE) :: spectrum
      real(double), intent(out) :: wavelength
      real(double) :: r, t
      real(double), intent(out) :: photonPacketWeight
      integer :: i

      call randomNumberGenerator(getDouble=r)
      call locate(spectrum%prob, spectrum%nLambda, r, i)
!      photonPacketWeight = spectrum%ppw(i)
      t = (r - spectrum%prob(i))/(spectrum%prob(i+1)-spectrum%prob(i))
      wavelength = spectrum%lambda(i) + t*(spectrum%lambda(i+1)-spectrum%lambda(i))
      photonPacketWeight = spectrum%ppw(i)! + t*(spectrum%ppw(i+1)-spectrum%ppw(i))
    end subroutine getWavelength

    function integrateNormSpectrumOverBand(spectrum, lam1 , lam2) result(tot)
      real(double) :: tot
      type(SPECTRUMTYPE) :: spectrum
      real(double) :: lam1, lam2, tlam1, tlam2
      integer :: i1, i2, i

      tlam1 = lam1
      tlam2 = lam2

      tlam1 = max(lam1, spectrum%lambda(1))
      tlam2 = min(lam2, spectrum%lambda(spectrum%nLambda))

      call locate(spectrum%lambda, spectrum%nLambda, tlam1, i1)
      call locate(spectrum%lambda, spectrum%nLambda, tlam2, i2)
      

      write(*,*) "spectrum%nlambda",spectrum%nlambda
      write(*,*) lam1,tlam1,lam2,tlam2
      write(*,*) i1,i2
      if (i1 == i2) then 
         tot = spectrum%normFlux(i1)*(tlam2 - tlam1)
      else
         tot = 0.d0
         tot = tot + spectrum%normFlux(i1)*(spectrum%lambda(i1+1)-tlam1)
         tot = tot + spectrum%normFlux(i2)*(tlam2-spectrum%lambda(i2))
         do i = i1, i2-1
            tot = tot + 0.5d0*(spectrum%normFlux(i+1)+spectrum%normFlux(i)) * &
                 (spectrum%lambda(i+1)-spectrum%lambda(i))
         enddo
      endif
    end function integrateNormSpectrumOverBand

    function sumPhotonsOverBand(spectrum, lam1 , lam2) result(tot)
      real(double) :: tot
      type(SPECTRUMTYPE) :: spectrum
      real(double) :: lam1, lam2, tlam1, tlam2, ePhoton
      integer :: i1, i2, i

      tlam1 = lam1
      tlam2 = lam2

      tlam1 = max(lam1, spectrum%lambda(1))
      tlam2 = min(lam2, spectrum%lambda(spectrum%nLambda))

      call locate(spectrum%lambda, spectrum%nLambda, tlam1, i1)
      call locate(spectrum%lambda, spectrum%nLambda, tlam2, i2)
      

      tot = 0.d0
      ePhoton = hCgs*cspeed/(tlam1*angstromtocm)
      tot = tot + spectrum%flux(i1)*(spectrum%lambda(i1+1)-tlam1)/ephoton
      ePhoton = hCgs*cspeed/(tlam2*angstromtocm)
      tot = tot + spectrum%flux(i2)*(tlam2-spectrum%lambda(i2))/ephoton
      do i = i1, i2-1
         ePhoton = hCgs*cspeed/(0.5*(spectrum%lambda(i+1)+spectrum%lambda(i))*angstromtocm)
         tot = tot + (0.5d0*(spectrum%flux(i+1)+spectrum%flux(i)) * &
              (spectrum%lambda(i+1)-spectrum%lambda(i))/ePhoton)
      enddo
    end function sumPhotonsOverBand

    function integrateSpectrumOverBand(spectrum, lam1 , lam2) result(tot)
      real(double) :: tot
      type(SPECTRUMTYPE) :: spectrum
      real(double) :: lam1, lam2, tlam1, tlam2
      integer :: i1, i2, i

      tlam1 = lam1
      tlam2 = lam2

      tlam1 = max(lam1, spectrum%lambda(1))
      tlam2 = min(lam2, spectrum%lambda(spectrum%nLambda))

      call locate(spectrum%lambda, spectrum%nLambda, tlam1, i1)
      call locate(spectrum%lambda, spectrum%nLambda, tlam2, i2)
      

      if (i1 == i2) then 
         tot = spectrum%flux(i1)*(tlam2 - tlam1)
      else
         tot = 0.d0
         tot = tot + spectrum%flux(i1)*(spectrum%lambda(i1+1)-tlam1)
         tot = tot + spectrum%flux(i2)*(tlam2-spectrum%lambda(i2))
         do i = i1, i2-1
            tot = tot + 0.5d0*(spectrum%flux(i+1)+spectrum%flux(i)) * &
                 (spectrum%lambda(i+1)-spectrum%lambda(i))
         enddo
      endif
    end function integrateSpectrumOverBand

    subroutine getWavelengthOverBand(spectrum, lam1 , lam2, wavelength)
      use random_mod, only: randomNumberGenerator
      real(double) :: wavelength
      type(SPECTRUMTYPE) :: spectrum
      real(double) :: lam1, lam2, tlam1, tlam2
      integer :: i1, i2, i
      real(double), allocatable :: prob(:), lam(:)
      real(double) :: t, r
      integer :: nProb

      tlam1 = lam1
      tlam2 = lam2

      tlam1 = max(lam1, spectrum%lambda(1))
      tlam2 = min(lam2, spectrum%lambda(spectrum%nLambda))

      call locate(spectrum%lambda, spectrum%nLambda, tlam1, i1)
      call locate(spectrum%lambda, spectrum%nLambda, tlam2, i2)
      
      if (i1 == i2) then 
         call randomNumberGenerator(getDouble=r)
         r = r - 0.5d0
         wavelength = spectrum%lambda(i1) + r * spectrum%dlambda(i1)
      else
         allocate(prob(spectrum%nLambda), lam(spectrum%nLambda))
         nProb = 1
         prob(nProb) = spectrum%normFlux(i1)*(spectrum%lambda(i1+1)-tlam1)
         lam(nProb) = 0.5d0*(tlam1+spectrum%lambda(i1+1))
         do i = i1, i2-1
            nProb = nProb + 1
            prob(nProb) = 0.5d0*(spectrum%normFlux(i+1)+spectrum%normFlux(i)) * &
                 (spectrum%lambda(i+1)-spectrum%lambda(i))
            lam(nProb) = 0.5d0*(spectrum%lambda(i+1)+spectrum%lambda(i))
         enddo
         nProb = nProb + 1
         prob(nProb) = spectrum%normFlux(i2)*(tlam2-spectrum%lambda(i2))
         lam(nProb) = 0.5d0*(tlam2+spectrum%lambda(i2))
         prob(1:nProb) = prob(1:nProb) - prob(1)
         prob(1:nProb) = prob(1:nProb)/prob(nProb)
         call randomNumberGenerator(getDouble=r)
         call locate(prob, nProb, r, i1)
         t = (r - prob(i1))/(prob(i1+1)-prob(i1))
         wavelength = lam(i1) + t * (lam(i1+1)-lam(i1))
         deallocate(prob, lam)
      endif
    end subroutine getWavelengthOverBand


    subroutine fillSpectrumBB(spectrum, teff, lamStart, lamEnd, nLambda, lamArray)
      use atom_mod, only: bLambda
      type(SPECTRUMTYPE) :: spectrum
      integer :: nLambda
      real(double) :: lamStart, lamEnd, teff
      real, optional :: lamArray(:)
      real(double) :: logLamStart, logLamEnd
      integer :: i
      
      if (associated(spectrum%flux)) deallocate(spectrum%flux)
      if (associated(spectrum%lambda)) deallocate(spectrum%lambda)
      if (associated(spectrum%dlambda)) deallocate(spectrum%dlambda)
      if (associated(spectrum%prob)) deallocate(spectrum%prob)
      if (associated(spectrum%ppw)) deallocate(spectrum%ppw)

      allocate(spectrum%flux(1:nLambda))
      allocate(spectrum%lambda(1:nLambda))
      allocate(spectrum%dlambda(1:nLambda))
      allocate(spectrum%prob(1:nLambda))
      allocate(spectrum%ppw(1:nLambda))

      logLamStart = log10(lamStart)
      logLamEnd = log10(lamEnd)
      
      if (.not.PRESENT(lamArray)) then
         do i = 1, nLambda
            spectrum%lambda(i) = logLamStart + real(i-1)/real(nLambda-1)*(logLamEnd - logLamStart)
            spectrum%lambda(i) = 10.**spectrum%lambda(i)
         enddo
      else
         spectrum%lambda(1:nLambda) = dble(lamArray(1:nLambda))
      endif
      do i = 2, nLambda-1
         spectrum%dlambda(i) = 0.5*((spectrum%lambda(i+1)+spectrum%lambda(i))-(spectrum%lambda(i)+spectrum%lambda(i-1)))
      enddo
      spectrum%dlambda(1) = spectrum%lambda(2)-spectrum%lambda(1)
      spectrum%dlambda(nLambda) = spectrum%lambda(nlambda)-spectrum%lambda(nLambda-1)
      
      do i = 1, nLambda
         spectrum%flux(i) = max(1.d-30,pi*bLambda(spectrum%lambda(i), dble(teff)) * 1.d-8) ! (Per cm to per angstrom)
      enddo
      spectrum%nLambda = nLambda
      where(spectrum%flux(1:spectrum%nLambda) == 0.d0) spectrum%flux = 1.d-100

      call probSpectrum(spectrum)
    end subroutine fillSpectrumBB

    subroutine addToSpectrumBB(spectrum, tBB, frac)
      use atom_mod, only: bLambda
      type(SPECTRUMTYPE) :: spectrum
      real(double) :: tBB, frac
      integer :: i

      do i = 1, spectrum%nLambda
         spectrum%flux(i) = spectrum%flux(i) + frac * pi*bLambda(spectrum%lambda(i), dble(tBB)) * 1.d-8 ! (per cm to per angstrom)
      enddo
      where(spectrum%flux(1:spectrum%nLambda) == 0.d0) spectrum%flux = 1.d-100

      call probSpectrum(spectrum)
    end subroutine addToSpectrumBB

    subroutine addXray(spectrum, frac)
      type(SPECTRUMTYPE) :: spectrum
      real(double) :: frac, totalFlux, xRayFlux, xRayFluxPerAngs
      real(double) :: lamStart, lamEnd, bolFlux, lambda
      integer :: i1, i2, i

      totalflux = integrateSpectrumOverBand(spectrum, 0.d0, 1.d30)

      xrayFlux = (frac / (1.d0 - frac)) * totalflux
      lamStart = (cSpeed / (10.d0 * 1000.d0 * evtoerg/hcgs))*1.d8
      lamEnd = (cSpeed / (0.09d0 * 1000.d0 * evtoerg/hcgs))*1.d8
      
      do i = 10, 1, -1
         lambda = 10.d0**(log10(lamStart) + (log10(lamend)-log10(lamstart))*dble(i-1)/9.d0)
         call insertWavelength(spectrum, lambda)
      enddo

      lamStart = (cSpeed / (10.d0 * 1000.d0 * evtoerg/hcgs))*1.d8
      lamEnd = (cSpeed / (0.1d0 * 1000.d0  * evtoerg/hcgs))*1.d8
      
      xRayFluxPerAngs = xRayFlux / (lamEnd - lamStart)


      call locate(spectrum%lambda, spectrum%nLambda, lamStart, i1)
      call locate(spectrum%lambda, spectrum%nLambda, lamEnd, i2)
      do i = i1, i2
         spectrum%flux(i) = spectrum%flux(i) + xRAyFluxPerAngs
      enddo
      bolFlux =  integrateSpectrumOverBand(spectrum, 0.d0, 1.d30)
      if (writeoutput) write(*,*) "X-ray flux added LX/LBol",xRayFlux/bolFlux
      

      call probSpectrum(spectrum)

    end subroutine addXray

      

    subroutine readSpectrum(spectrum, filename, ok)
      type(SPECTRUMTYPE) :: spectrum
      logical :: ok
      character(len=*) :: filename
      real :: fTemp(120000),xTemp(120000), x, f
      character(len=80) :: cLine
      integer :: nLambda, i

      ok = .true.
      open(20,file=filename,form="formatted",status="old", err=666)
      nLambda = 1
10    continue
      read(20,'(a)',end=20) cline
      if (len(trim(cline)).gt.0) then
         if (cLine(1:1) /="#") then
            read(cLine,*,err=10) x, f
         else 
            goto 10
         endif
      else
         goto 10
      endif
      xtemp(nLambda) = x
      fTemp(nLambda) = max(f,1.e-30)
      nLambda = nLambda + 1
      goto 10
20    continue
      close(20)
      nLambda = nLambda - 1
      if (nLambda == 0) then
         call writeFatal("Error reading continuum flux file: "//trim(filename))
      endif
      allocate(spectrum%flux(1:nLambda))
      allocate(spectrum%lambda(1:nLambda))
      allocate(spectrum%dlambda(1:nLambda))
      allocate(spectrum%prob(1:nLambda))
      allocate(spectrum%ppw(1:nLambda))
      spectrum%nLambda = nLambda
      spectrum%flux(1:nLambda) = fTemp(1:nLambda)
      spectrum%lambda(1:nLambda) = xTemp(1:nLambda)
      do i = 2, nLambda-1
         spectrum%dlambda(i) = 0.5*((xtemp(i+1)+xtemp(i))-(xtemp(i)+xtemp(i-1)))
      enddo
      spectrum%dlambda(1) = xtemp(2)-xtemp(1)
      spectrum%dlambda(nLambda) = xtemp(nlambda)-xtemp(nLambda-1)
      call probSpectrum(spectrum)

      goto 999

666   continue
      ok = .false.
      call writeFatal("Error in opening file: "//trim(filename))
999   continue

    end subroutine readSpectrum

    subroutine readSpectrumFromDump(spectrum, lunit)
      type(SPECTRUMTYPE) :: spectrum
      integer ::lunit, nLambda
      
      read(lunit) nLambda
!      write(*,*) myrankWorldGlobal, " reading nlambda ",nlambda
      allocate(spectrum%flux(1:nLambda))
      allocate(spectrum%lambda(1:nLambda))
      allocate(spectrum%dlambda(1:nLambda))
      allocate(spectrum%prob(1:nLambda))
      allocate(spectrum%ppw(1:nLambda))
      spectrum%nLambda = nLambda
      read(lunit) spectrum%flux(1:nLambda)
!      write(*,*) myrankWorldGlobal, " flux " , spectrum%flux(1:nlambda)
      read(lunit) spectrum%lambda(1:nLambda)
      read(lunit) spectrum%dlambda(1:nLambda)
      read(lunit) spectrum%prob(1:nLambda)
      read(lunit) spectrum%ppw(1:nLambda)
    end subroutine readSpectrumFromDump

    subroutine writeSpectrumToDump(spectrum, lunit)
      type(SPECTRUMTYPE) :: spectrum
      integer ::lunit
      
      write(lunit) spectrum%nLambda
      write(lunit) spectrum%flux(1:spectrum%nLambda)
      write(lunit) spectrum%lambda(1:spectrum%nLambda)
      write(lunit) spectrum%dlambda(1:spectrum%nLambda)
      write(lunit) spectrum%prob(1:spectrum%nLambda)
      write(lunit) spectrum%ppw(1:spectrum%nLambda)
    end subroutine writeSpectrumToDump

    subroutine insertWavelength(spectrum, lambda)
      use utils_mod, only: loginterp_dble
      type(SPECTRUMTYPE) :: spectrum, tmpSpectrum
      logical :: ok
      integer :: i, newnLambda
      real(double) :: lambda, flux

      ok = .true.
      do i = 1, spectrum%nlambda
         if (lambda == spectrum%lambda(i)) ok = .false.
      enddo
      if (.not. ok) goto 666

      if (lambda < spectrum%lambda(1)) then 
         flux = 1.d-30
      else if (lambda > spectrum%lambda(spectrum%nLambda)) then
         flux = 1.d-30
      else
         flux = loginterp_dble(spectrum%flux, spectrum%nlambda, spectrum%lambda, lambda)
      endif
      
      newnLambda = spectrum%nlambda + 1
      
      allocate(tmpspectrum%flux(1:newnLambda))
      allocate(tmpspectrum%lambda(1:newnLambda))
      allocate(tmpspectrum%dlambda(1:newnLambda))
      allocate(tmpspectrum%prob(1:newnLambda))
      allocate(tmpspectrum%ppw(1:newnLambda))
      tmpSpectrum%nLambda = newNLambda

      if (lambda < spectrum%lambda(1)) then 
         tmpSpectrum%lambda(1) = lambda
         tmpSpectrum%flux(1) = flux 
         tmpSpectrum%lambda(2:newnLambda) = spectrum%lambda(1:spectrum%nLambda)
         tmpSpectrum%flux(2:newnLambda) = spectrum%flux(1:spectrum%nLambda)

     else if (lambda > spectrum%lambda(spectrum%nLambda)) then

         tmpSpectrum%lambda(newnLambda) = lambda
         tmpSpectrum%flux(newnLambda) = flux 
         tmpSpectrum%lambda(1:newnLambda-1) = spectrum%lambda(1:spectrum%nLambda)
         tmpSpectrum%flux(1:newnLambda-1) = spectrum%flux(1:spectrum%nLambda)


      else

       call locate(spectrum%lambda, spectrum%nlambda, lambda, i)

         tmpSpectrum%lambda(1:i) = spectrum%lambda(1:i)
         tmpSpectrum%flux(1:i) = spectrum%flux(1:i)

         tmpSpectrum%lambda(i+1) = lambda
         tmpSpectrum%flux(i+1) = flux 

         tmpSpectrum%lambda(i+2:newnLambda) = spectrum%lambda(i+1:spectrum%nLambda)
         tmpSpectrum%flux(i+2:newnLambda) = spectrum%flux(i+1:spectrum%nLambda)
       
      endif


      call freeSpectrum(spectrum)
      call copySpectrum(Spectrum, tmpspectrum)
      do i = 2, spectrum%nLambda-1
         spectrum%dlambda(i) = 0.5*((spectrum%lambda(i+1)+spectrum%lambda(i)) - &
              (spectrum%lambda(i)+spectrum%lambda(i-1)))
      enddo
      spectrum%dlambda(1) = spectrum%lambda(2)-spectrum%lambda(1)
      spectrum%dlambda(spectrum%nLambda) = spectrum%lambda(spectrum%nlambda) - &
           spectrum%lambda(spectrum%nLambda-1)
      
      
      666 continue
     end subroutine insertWavelength



    subroutine copySpectrum(a, b)
      type(SPECTRUMTYPE) :: a, b
      integer :: nLambda
      nLambda = b%nlambda
      allocate(a%flux(1:nLambda))
      allocate(a%lambda(1:nLambda))
      allocate(a%dlambda(1:nLambda))
      allocate(a%prob(1:nLambda))
      allocate(a%ppw(1:nLambda))
      a%nLambda = nLambda
      a%flux = b%flux
      a%lambda = b%lambda
      a%dlambda = b%dlambda
      a%prob = b%prob
      a%ppw = b%ppw
    end subroutine copySpectrum

    subroutine probSpectrum(spectrum)
      use inputs_mod, only : biasToLyman, biasMagnitude
      type(SPECTRUMTYPE) :: spectrum
      real(double) :: fac
      integer :: i
      real(double) :: unbiasedProb(spectrum%nLambda)
      spectrum%prob = 0.d0
      unbiasedProb = 0.d0
      spectrum%ppw = 0.d0

      do i = 2, spectrum%nLambda
         fac = 1.d0
         if (biasToLyman) then
            if (spectrum%lambda(i) < 912.) then
               fac = biasMagnitude
            else
               fac = 1.d0
            endif
            !record both the biased and unbiased probabilities
            spectrum%prob(i) = spectrum%prob(i-1) + spectrum%flux(i) * spectrum%dLambda(i) * fac
            unbiasedProb(i) = unbiasedprob(i-1) + spectrum%flux(i) * spectrum%dLambda(i)
         else
            spectrum%prob(i) = spectrum%prob(i-1) + spectrum%flux(i) * spectrum%dLambda(i) * fac
            spectrum%ppw(i) = 1.d0
         endif
      enddo
      
      spectrum%prob(1:spectrum%nLambda) = spectrum%prob(1:spectrum%nLambda) / spectrum%prob(spectrum%nLambda)

      if(biasToLyman) then
         !Normalize biased probability distribution if required
         unbiasedProb(1:spectrum%nLambda) = unbiasedProb(1:spectrum%nLambda) / unbiasedProb(spectrum%nLambda)

         !Packet weight is the ratio of normalised unbiased to biased probabilities
         !These will be interpolated between later !Thaw
         do i = 2, spectrum%nLambda
            spectrum%ppw(i) = abs(unbiasedProb(i+1) - unbiasedProb(i)) / abs(spectrum%prob(i+1) - spectrum%prob(i))
         end do
      end if
     
    end subroutine probSpectrum

    subroutine normalizedSpectrum(spectrum)
      type(SPECTRUMTYPE) :: spectrum
      allocate(spectrum%normflux(1:spectrum%nLambda))
      spectrum%normflux(1:spectrum%nLambda) = spectrum%flux(1:spectrum%nLambda) &
           / SUM(spectrum%flux(1:spectrum%nLambda)*spectrum%dlambda(1:spectrum%nLambda))
    end subroutine

    real(double) function returnNormValue2(spectrum, lambda, lam1, lam2)
      type(SPECTRUMTYPE) :: spectrum
      real(double) :: lambda, t, lam1, lam2, tot
      integer :: i, i1, i2

      if ((lambda < lam1).or.(lambda > lam2)) then
         returnNormValue2 = 0.d0 
      else

         call locate(spectrum%lambda, spectrum%nLambda, lam1, i1)
         call locate(spectrum%lambda, spectrum%nLambda, lam2, i2)
      
         if (i1 == i2) then 
            returnNormValue2 = spectrum%flux(i1)/(spectrum%flux(i1)*(lam2-lam1))
         else
            tot = 0.d0
            tot = tot + spectrum%flux(i1)*(spectrum%lambda(i1+1)-lam1)
            tot = tot + spectrum%flux(i2)*(lam2-spectrum%lambda(i2))
            do i = i1, i2-1
               tot = tot + 0.5d0*(spectrum%flux(i+1)+spectrum%flux(i)) * &
                    (spectrum%lambda(i+1)-spectrum%lambda(i))
            enddo
            call locate(spectrum%lambda, spectrum%nLambda, lambda, i)
            t = (lambda - spectrum%lambda(i))/(spectrum%lambda(i+1)-spectrum%lambda(i))
            t = (spectrum%flux(i)+t*(spectrum%flux(i+1)-spectrum%flux(i)))
            returnNormValue2 = t/tot
         endif
      endif

    end function returnNormValue2
         
    subroutine fillSpectrumKurucz(spectrum, teff, mass, radius, freeup)
      type(SPECTRUMTYPE) :: spectrum, spec1, spec2
      real(double) :: teff, mass, radius
      real(double) :: logg
      logical, optional :: freeUp
      logical :: ok1, ok2, ok
      integer, parameter :: nFiles = 60
      real,save :: teffArray(nFiles)
      integer :: i, j
      real(double) :: t
      real :: loggArray(11)
      logical, save :: firstTime = .true.
      character(len=200) :: thisFile1, thisFile2, dataDirectory
      character(len=80) :: label1, label2, message
      integer :: i1, i2
      integer, parameter :: nKurucz = 410
      type(SPECTRUMTYPE),save :: kSpectrum(nKurucz)
      character(len=80),save :: klabel(nKurucz)
      if (firstTime) call  readKuruczGrid(klabel, kspectrum, nKurucz)
         


      logg = log10(bigG * mass / radius**2)

      ok = .true.
      call unixGetenv("TORUS_DATA", dataDirectory, i)


      loggArray = (/ 000., 050., 100., 150., 200., 250., 300., 350., 400., 450., 500. /)
      if (firsttime) then
         open(31, file=trim(dataDirectory)//"/Kurucz/filelist.dat", form="formatted", status="old")
         do i = 1, nFiles
            read(31, *) teffArray(i)
         end do
         close(31)
         firstTime = .false.
      endif




      if ((teff > teffArray(1)).and.(teff < teffArray(nFiles))) then
         call locate(teffArray, nFiles, real(teff), i)
         call locate(loggArray, 11, real(logg*100.), j)

         t = ((logg*100.) - loggArray(j))/(loggArray(j+1) - loggArray(j))
         if (t  > 0.5) j = j + 1
         t = (teff - teffArray(i))/(teffArray(i+1)-teffArray(i))
         i1 = i
         i2 = i + 1
         call createKuruczFilename(teffArray(i1), loggArray(j), thisFile1, label1)
         call createKuruczFilename(teffArray(i2), loggArray(j), thisFile2, label2)
         
         write(message, '(a,a)') "Interpolating Kurucz atmospheres between: ",trim(label1)
         call writeInfo(message, TRIVIAL)
         write(message, '(a,a)') "                                        : ",trim(label2)
         call writeInfo(message, TRIVIAL)
         
         call readKuruczSpectrum(spec1, label1, klabel, kspectrum, nKurucz, ok1)
         if (.not.ok1) then
            if (writeoutput) then
               write(*,*) "Can't find kurucz spectrum: ",trim(thisfile1)," ",trim(label1)
               !            do i = 1, nKurucz
               !               write(*,*) trim(label1), " and ", trim(klabel(i))
               !            enddo
!            stop
            endif
         endif
         call readKuruczSpectrum(spec2, label2, klabel, kspectrum, nKurucz, ok2)
         if (.not.ok2) then
            if (writeoutput) write(*,*) "Can't find kurucz spectrum: ",trim(thisfile2), " ", trim(label2)
         endif
         
         if (ok1.and.ok2) then
            call createInterpolatedSpectrum(spectrum, spec1, spec2, t)
         else
            call fillSpectrumBB(spectrum, teff, 10.d0, 1.d7, 200)
         endif
      else
            call fillSpectrumBB(spectrum, teff, 10.d0, 1.d7, 200)
         endif


         if (PRESENT(freeUp)) then
         if (freeup) then
            do i = 1, nKurucz
               call freeSpectrum(kspectrum(i))
            enddo
            firstTime = .true.
         endif
      endif
      call freeSpectrum(spec1)
      call freeSpectrum(spec2)
      call probSpectrum(spectrum)

    end subroutine fillSpectrumKurucz

     subroutine createInterpolatedSpectrum(spectrum, spec1, spec2, t)
       use utils_mod, only: loginterp_dble
       type(SPECTRUMTYPE) :: spectrum, spec1, spec2
       real(double) :: t
       integer :: i 
       real(double) :: y1, y2
       call copySpectrum(spectrum, spec1)

       do i = 1, spectrum%nLambda
          y1 = spec1%flux(i)

          y2 = loginterp_dble(spec2%flux, spec2%nlambda, spec2%lambda, spec1%lambda(i))

          spectrum%flux(i) = 10.d0**(log10(y1) + t * (log10(y2) - log10(y1)))
       enddo
     end subroutine createInterpolatedSpectrum


     subroutine createKuruczFileName(teff, logg, thisfile, thislabel)
       real :: teff, logg
       integer :: i
       character(len=*) thisfile, thisLabel
       character(len=80) :: fluxfile, dataDirectory

       call unixGetenv("TORUS_DATA", dataDirectory, i)


       if (teff < 10000.) then
          write(fluxfile,'(a,i4,a,i3.3,a)') "f",int(teff),"_",int(logg),".dat"
       else
          write(fluxfile,'(a,i5,a,i3.3,a)') "f",int(teff),"_",int(logg),".dat"          
       endif
       thisLabel = fluxfile
       thisFile = trim(dataDirectory)//"/Kurucz/"//trim(fluxfile)
     end subroutine createKuruczFileName

     subroutine readKuruczGrid(label, spectrum, nFiles)
       character(len=*) :: label(:)
       type(SPECTRUMTYPE) :: spectrum(:)
       integer :: nFiles
       character(len=200) :: tfile,fluxfile,dataDirectory 
       logical :: ok
       integer :: i
       ok = .true.

       call unixGetenv("TORUS_DATA", dataDirectory, i)
       
       call writeInfo("Reading Kurucz grid...",TRIVIAL)

       tfile = trim(dataDirectory)//"/Kurucz/files.dat"
       open(31, file = tfile, status = "old", form="formatted")
       do i = 1, nFiles
          read(31,*) fluxfile
          label(i) = trim(fluxfile)
          tfile = trim(dataDirectory)//"/Kurucz/"//trim(fluxfile)
          call readSpectrum(spectrum(i), tfile, ok)
          spectrum(i)%flux = spectrum(i)%flux * pi ! astrophysical to real fluxes
       enddo
       close(31)
       call writeInfo("Done.",TRIVIAL)

     end subroutine readKuruczGrid

     subroutine readKuruczSpectrum(thisSpectrum,thisLabel, label, spectrum, nFiles, ok)
       character(len=*) :: label(:), thisLabel
       type(SPECTRUMTYPE) :: spectrum(:), thisSpectrum
       integer :: nFiles
       logical :: ok
       integer :: i
       
       ok = .false.
       do i = 1, nFiles
          if (trim(label(i)).eq.trim(thisLabel)) then
             call copySpectrum(thisSpectrum, spectrum(i))
             ok = .true.
             exit
          endif
       enddo
     end subroutine readKuruczSpectrum


  end module spectrum_mod
