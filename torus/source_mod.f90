 module source_mod
  
  use spectrum_mod
  use constants_mod
  use vector_mod
  use utils_mod
  use octal_mod
  use gridtype_mod
  use surface_mod, only : surfacetype, ishot, getelement, writesurface, readsurface, &
       getphotovec, buildspherenbody, emptysurface, buildsphere

  implicit none

  public 


  type SOURCETYPE
     type(VECTOR)    :: position   ! [10^10cm]
     type(VECTOR)    :: force
     type(VECTOR)    :: velocity
     real(double) :: radius     ! [10^10cm]
     real(double) :: luminosity ! [erg/s]
     real(double) :: teff       ! [K]
     real(double) :: initialMass ! [solar]
     real(double) :: mass       ! [solar]
     real(double) :: age        ! [years]
     real(double) :: mdot        ! [mDot g/s]
     real(double) :: mdotWind    ! [mDot g/s]
     type(SPECTRUMTYPE)    :: spectrum   ! [???]
     type(SURFACETYPE) :: surface
     real(double) :: limbDark(2)
     logical :: outsideGrid
     logical :: onEdge
     logical :: onCorner
     logical :: pointSource
     real(double) :: distance
     real(double) :: prob ! probability of packet from this source
     logical :: stellar
     logical :: viscosity
     logical :: diffuse ! isrf, cmb
     type(VECTOR) :: angMomentum
     real(double) :: accretionRadius
     real(double) :: time
  end type SOURCETYPE

  type(SOURCETYPE), pointer :: globalSourceArray(:) => null()
  integer :: globalnSource

  contains

    subroutine writeSourceHistory(rootfilename, source, nSource)
      use inputs_mod, only : iModel
      type(SOURCETYPE) :: source(:)
      character(len=*) :: rootFilename
      character(len=80) :: filename
      integer :: nSource
      integer :: i
      logical, save :: firstTime = .true.
      if (writeoutput) then
         do i = 1, nSource
            write(filename,'(a,i3.3,a)') trim(rootFilename),i,".dat"
            write(*,*) "filename ",trim(filename)
            if (firstTime) then
               open(22, file=filename, status="unknown", form="formatted")
               write(22,'(a4,a12,a9,a9,a9,a12,a12,a12,a12,a12,a12,a12)') &
                    "#  N","Age (yr)","Mass","Radius","Teff","Lum","Mdot","x","y","z"
            else
               open(22, file=filename, status="old", form="formatted", position="append")
            endif
            write(22,'(i4.4, 1p, e12.3, 0p, f9.4, f9.4, f9.1, 1p, e12.3, e12.3, e12.3, e12.3, e12.3,e12.3,e12.3)') &
                 iModel,source(i)%time*secstoYears, &
                 source(i)%mass/msol, &
                 source(i)%radius*1.d10/rsol, source(i)%teff, source(i)%luminosity/lsol, source(i)%mdot/msol/secstoyears, &
                 source(i)%position%x, &
                 source(i)%position%y,source(i)%position%z
            close(22)
         enddo
         if (firstTime) firstTime = .false.
      endif
    end subroutine writeSourceHistory

    subroutine writeSourceList(source, nSource)
      type(SOURCETYPE) :: source(:)
      integer :: nSource
      real(double) :: fromSpec
      real(double) :: tot
      integer :: i, j


      if ( nsource<1 ) return

      if (writeoutput) &
           write(*,'(a4,a2,a8,a8,a12,a12,a12,a12,a12,a12,a12,a12,a12,a12)') &
           "Ran ","N","M","R","L","L(spec)","x","y","z","vx","vy","vz"

      
      do i = 1, nSource
         fromSpec = sumSourceLuminosity(source(i:i), 1, 1.0, 1.e30)
         if (writeoutput) &
              write(*, &
         '(i3.3,1x,i2.2, f8.3, f8.3, 1p, e12.3, e12.3, e12.3, e12.3,e12.3,e12.3,e12.3,e12.3,e12.3)') &
              myrankWorldGlobal, i,&
              source(i)%mass/msol, &
              source(i)%radius*1.d10/rsol, source(i)%luminosity/lsol, fromspec/lsol, source(i)%position%x, &
              source(i)%position%y,source(i)%position%z, source(i)%velocity%x/1.d5, source(i)%velocity%y/1.d5, &
              source(i)%velocity%z/1.d5
         tot = 0.d0
         do j = 1, source(i)%spectrum%nLambda
            if (j < source(i)%spectrum%nlambda) then
               tot = tot + source(i)%spectrum%flux(j) * source(i)%spectrum%dlambda(j)
            endif
         enddo
         tot = tot * fourPi*source(i)%radius**2 * 1.d20
!         write(*,*) myrankWorldGlobal, " source lum ",tot/lsol
      enddo
    end subroutine writeSourceList


    subroutine writeIonizingFluxes(source, nSource)
      type(SOURCETYPE) :: source(:)
      integer :: nSource
 !     real(double) :: fromSpec
      real(double) :: tot, max
      integer :: i!, j
      



      if ( nsource<1 ) return

!      if (writeoutput) &!
!           write(*,'(a4,a2,a8,a8,a12,a12,a12,a12,a12,a12,a12,a12,a12,a12)') &
!           "Ran ","N","M","R","L","L(spec)","x","y","z","vx","vy","vz"

      print *, "ionizing fluxes are"
      tot = 0.d0
      max = 0.d0
      do i = 1, nSource
!         fromSpec = sumSourceLuminosity(source(i:i), 1, 1.0, 1.e30)
!         print *, "mass ", source(i)%mass/msol
         if(source(i)%mass/msol >= 20.d0) then

!            if (writeoutput) &
                 write(*,*)  myrankWorldGlobal, i, ionizingflux(source(i))
                 tot = tot + ionizingflux(source(i))
                 if(ionizingflux(source(i)) > max) then
                    max = ionizingflux(source(i))
                 endif
         endif
!         tot = 0.d0
!         do j = 1, source(i)%spectrum%nLambda
         !            if (j < source(i)%spectrum%nlambda) then
!               tot = tot + source(i)%spectrum%flux(j) * source(i)%spectrum%dlambda(j)
!            endif
!         enddo
 !        tot = tot * fourPi*source(i)%radius**2 * 1.d20
!         write(*,*) myrankWorldGlobal, " source lum ",tot/lsol
      enddo
      print *, "TOTAL IONIZING FLUX FROM SINKS IS " ,tot
      print *, "HIGHEST IONIZING FLUX IS ", max
    end subroutine writeIonizingFluxes

    subroutine writeSourceArray(rootfilename)
      use inputs_mod, only : iModel
      use utils_mod, only : findMultiFilename
      character(len=*) :: rootFilename
      character(len=80) :: filename
      integer :: iSource
      integer, parameter :: lunit = 65

      call findMultiFilename(rootFilename, iModel, filename)



      open(lunit, file=filename, form="unformatted", status="unknown")
      write(lunit) globalnSource
      do isource = 1, globalnSource
         call writesource(globalsourceArray(isource), lunit)
      enddo
      close(lunit)
    end subroutine writeSourceArray

    subroutine readSourceArray(nSource, source, rootfilename)
      use inputs_mod, only : iModel
      use utils_mod, only : findMultiFilename
      type(SOURCETYPE) :: source(:)
      integer :: nSource
      character(len=*) :: rootfilename
      character(len=120) :: filename, message
      integer :: iSource
      integer, parameter :: lunit = 65


      call findMultiFilename(rootFilename, iModel, filename)
      write(message,*) "Reading sources from file: ",trim(filename)
!      write(*,*) myrankWorldGlobal,trim(message)
      call writeInfo(message,TRIVIAL)

      if (SIZE(source) > 0) then
         do iSource = 1, SIZE(source)
            call freeSource(source(isource))
         enddo
      endif

      open(lunit, file=filename, form="unformatted", status="old")
      read(lunit) nSource

      do iSource = 1, nSource
         call freeSource(source(isource))
         call readSource(source(isource), lunit)
      enddo
      close(lunit)
    end subroutine readSourceArray

    subroutine writeSource(source, lunit)
      type(SOURCETYPE) :: source
      integer :: lunit

      write(lunit) source%position
      write(lunit) source%velocity
      write(lunit) source%force
      write(lunit) source%radius
      write(lunit) source%luminosity
      write(lunit) source%teff
      write(lunit) source%initialMass
      write(lunit) source%mass
      write(lunit) source%age
      write(lunit) source%mdot
      write(lunit) source%outsideGrid
      write(lunit) source%onCorner
      write(lunit) source%distance
      write(lunit) source%stellar
      write(lunit) source%viscosity
      write(lunit) source%diffuse
      write(lunit) source%angMomentum
      write(lunit) source%accretionRadius
      write(lunit) source%time
      write(lunit) source%mdotwind
      call writeSpectrumToDump(source%spectrum,lunit)
      call writeSurface(source%surface, lunit)
    end subroutine writeSource

    subroutine readSource(source, lunit)
      type(SOURCETYPE) :: source
      integer :: lunit

      read(lunit) source%position
      read(lunit) source%velocity
      read(lunit) source%force
      read(lunit) source%radius
      read(lunit) source%luminosity
      read(lunit) source%teff
      read(lunit) source%initialMass
      read(lunit) source%mass
      read(lunit) source%age
      read(lunit) source%mdot
      read(lunit) source%outsideGrid
      read(lunit) source%onCorner
      read(lunit) source%distance
      read(lunit) source%stellar
      read(lunit) source%viscosity
      read(lunit) source%diffuse
      read(lunit) source%angMomentum
      read(lunit) source%accretionRadius
      read(lunit) source%time
      read(lunit) source%mdotwind
      call freeSpectrum(source%spectrum)
      call readSpectrumFromDump(source%spectrum,lunit)
      call readSurface(source%surface, lunit)
    end subroutine readSource


    subroutine freeGlobalSourceArray()
      integer :: i
      if (globalnSource > 0) then
         do i =  1, globalNSource
            call freeSource(globalSourceArray(i))
         enddo
         deallocate(globalSourcearray)
      endif
    end subroutine freeGlobalSourceArray

    subroutine freeSource(source)
      type(SOURCETYPE) :: source
      call freeSpectrum(source%spectrum)
      call emptySurface(source%surface)
    end subroutine freeSource

    function ionizingFlux(source) result(flux)
      type(sourcetype) :: source
      real(double) :: flux
      flux = sumPhotonsOverBand(source%spectrum, 10.d0, 912.d0)
      if (source%outsidegrid) then
         flux = flux * (source%radius*1.d10)**2/(source%distance**2)
      else
         flux = flux * fourPi*(source%radius*1.d10)**2
      endif
    end function ionizingFlux


!    subroutine randomSource(source, nSource, iSource)
!      integer :: nSource
!      type(SOURCETYPE) :: source(1:nSource)
!      integer :: iSource
!      real, save, allocatable :: prob(:)
!      real :: r, t
!      real :: lRatio
!      integer :: i
!      logical, save :: first_time = .true.
!
!      if (nSource == 1) then
!         iSource = 1
!      else if (nSource == 2) then
!         lRatio = source(1)%luminosity / (source(1)%luminosity + source(2)%luminosity)
!         call randomNumberGenerator(getDouble=r)
!         if ( r  < lRatio ) then
!            iSource = 1
!         else
!            iSource = 2
!         endif
!      else if (nSource > 2) then
!	 if (first_time) then
!	    ! allocate array
!	    ALLOCATE(prob(1:nSource))
!	    ! Create the prob. dist. function.
!	    prob(1:nSource) = source(1:nSource)%luminosity/lsol
!	    do i = 2, nSource
!	       prob(i) = prob(i) + prob(i-1)
!	    enddo
!	    prob(1:nSource) = prob(1:nSource) - prob(1)
!	    prob(1:nSource) = prob(1:nSource) / prob(nSource)
!
!	    first_time = .false.
!	    
!	 end if
!	 
!         call randomNumberGenerator(getDouble=r)
!         call locate(prob, nSource, r, iSource)
!         if (iSource < nSource) then
!            t = (r - prob(iSource))/(prob(iSource+1) - prob(iSource))
!            if (t > 0.5) iSource = iSource + 1
!         endif
!      endif
!
!    end subroutine randomSource


    subroutine randomSource(source, nSource, iSource, weight,lamArray, nLambda, initialize, lambdaMono)

      integer :: nSource
      type(SOURCETYPE) :: source(:)
      integer, intent(out) :: iSource
      real(double), save, allocatable :: prob(:), weightArray(:)
      real, optional :: lamArray(:)
      real(double) :: weight
      integer,optional :: nlambda
      real(double) :: r
      integer :: i
      logical, optional :: initialize
! Wavelength used when initialize is present and nLambda=1
      real, optional    :: lambdaMono

!$OMP THREADPRIVATE (prob, weightArray)

      if (nSource == 1) then
         iSource = 1
         weight = 1.d0
      else
         if (PRESENT(initialize)) then
	    ! allocate array
            if (allocated(prob)) then
               deallocate(prob)
            endif
            if (allocated(weightArray)) then
               deallocate(weightArray)
            endif
            ALLOCATE(prob(1:nSource))
            prob = 0.d0
	    ! Create the prob. dist. function.
            if (nLambda > 1) then
               do i = 1, nSource
                  prob(i) = integrateSpectrumOverBand(source(i)%spectrum, dble(lamArray(1)) , &
                       dble(lamArray(nLambda))) * (fourPi * (source(i)%radius*1.d10)**2) /lsol
               enddo
            else
               if (present(lambdaMono)) then
                  do i = 1, nSource
                     prob(i) = sourceLuminosityMonochromatic(source(i), dble(lambdaMono)) * &
                          (fourPi * (source(i)%radius*1.d10)**2)
                  enddo
               else
                  do i = 1, nSource
                     prob(i) = sourceLuminosityMonochromatic(source(i), dble(lamArray(1))) * &
                          (fourPi * (source(i)%radius*1.d10)**2)
                  enddo
               end if
            endif
            prob(1:nSource) = prob(1:nSource)/SUM(prob(1:nSource))
            allocate(weightArray(1:nSource))
            if (SUM(source(1:nSource)%prob) == 0.d0) then
               weightArray(1:nSource) = 1.d0
            else
               weightArray(1:nSource) = prob(1:nSource) / source(1:nSource)%prob
               prob(1:nSource) = source(1:nSource)%prob
            endif
            if (nSource > 2) then
               do i = 2, nSource
                  prob(i) = prob(i) + prob(i-1)
               enddo
               prob(1:nSource) = prob(1:nSource) / prob(nSource)
            endif
!            if (writeoutput) write(*,*) "prob ",prob(1:nSource)
!            if (writeoutput) write(*,*) "source%prob ",source(1:nSource)%prob
!            if (writeoutput) write(*,*) "weight ",weightArray(1:nSource)
         end if
	 
         if (nSource > 2) then
!            do i = 1, nSource
!               if (writeoutput) write(33,*) i,prob(i)
!            enddo
!            if (writeoutput) stop
            call randomNumberGenerator(getDouble=r)
            if (r < prob(1)) then
               iSource = 1
            else
               call locate(prob, nSource, r, iSource)
               iSource = iSource + 1
            endif
            weight = weightArray(isource)
         else
            call randomNumberGenerator(getDouble=r)
            if (r < prob(1)) then
               iSource = 1
            else 
               isource = 2
            endif
            weight = weightArray(isource)
         endif
      endif
    end subroutine randomSource




    real(double) function sumSourceLuminosity(source, nsource, lam1, lam2) result (tot)
      integer :: nSource
      type(SOURCETYPE) :: source(:)
      real :: lam1, lam2
      integer :: i
      character(len=80) :: message

      tot = 0
      do i = 1, nSource
         if (associated(source(i)%spectrum%lambda)) then
            tot = tot + integrateSpectrumOverBand(source(i)%spectrum, dble(lam1) , &
                 dble(lam2)) * (fourPi*(1.d10*source(i)%radius)**2)
         else
            write(message,*) "source_mod: spectrum%lambda not allocated for source", i
            call writeFatal(message)
         endif
      enddo
    end function sumSourceLuminosity

    subroutine testRandomSource(source, nsource)
    integer :: nSource, isource
    type(SOURCETYPE) :: Source(:)
    real(double) :: totlum, tot
    integer :: i, n
    real(double) :: eps, packetweight
    real :: lamArray(2)
    integer :: nLambda

    nlambda = 2
    lamArray(1) = 100.d0
    lamArray(2) = 1.d7
      call randomSource(source, nSource, isource, weight=packetWeight, lamArray=lamArray, nLambda=nLambda, &
          initialize = .true.)
    n = 1000000
    totlum = sum(source%luminosity)
    eps = totLum / dble(n)
    
    tot = 0.d0
    do i = 1, n
    call randomSource(source, nSource, iSource, weight=packetWeight)
       tot = tot + eps * packetWeight
    enddo
    if (writeoutput) then
       write(*,*) "Test random source: expect luminosity ",totlum
       write(*,*) "Test random source: found luminosity ",tot
    endif
    if (abs(totlum-tot)/tot > 0.01d0) then
       call writeFatal("Error in randomSource")
    endif
  end subroutine testRandomSource

    subroutine testSourceSpectrum(source, nsource)
    integer :: nSource, isource
    type(SOURCETYPE) :: Source(:)
    integer :: i
    real(double) :: packetweight, wavelength
    do isource = 1, nSource
       do i = 1 , 10
           call getWavelength(source(isource)%spectrum, wavelength, PacketWeight)
           write(*,*) i, wavelength,packetweight
        enddo
     end do
   end subroutine testSourceSpectrum

    real(double) function sumSourceLuminosityMonochromatic(grid, source, nsource, lam) result (tot)
      integer :: nSource
      type(SOURCETYPE) :: source(:)
      type(GRIDTYPE) :: grid
      real(double) :: lam, flux
      integer :: i !, j

      tot = 0.d0
      do i = 1, nSource
         if (lam > source(i)%spectrum%lambda(source(i)%spectrum%nLambda)) then
           tot = tot + 1.d-200
         else
!           call locate(source(i)%spectrum%lambda, source(i)%spectrum%nLambda, lam, j)

           flux =  loginterp_dble(source(i)%spectrum%flux(1:source(i)%spectrum%nLambda), &
                source(i)%spectrum%nLambda, source(i)%spectrum%lambda(1:source(i)%spectrum%nLambda), lam)
!THAW -old
!              tot = tot + flux * fourPi * (source(i)%radius*1.d10)**2 

!THaw - new
           if (.not.source(i)%outsideGrid) then                                                    
              tot = tot + flux * fourPi * (source(i)%radius*1.d10)**2
           else                          
              if(cart2d) then
                 if(grid%octreeroot%twod) then
                    tot = tot + flux *  (2.d0*grid%octreeRoot%subcellSize*1.d20*2.d0* &
                         grid%halfsmallestsubcell) * &                  
                         (source(i)%radius*1.d10)**2 / (source(i)%distance**2)                 
                 else
                    tot = tot + flux *  (1.d20*(2.d0* &
                         grid%halfsmallestsubcell)**2) * &  
                         (source(i)%radius*1.d10)**2 / (source(i)%distance**2)    
                 end if
              else
                 tot = tot + flux *  (2.d0*grid%octreeRoot%subcellSize*1.d10)**2 * &                  
                      (source(i)%radius*1.d10)**2 / (source(i)%distance**2)
              end if
           endif                 

        endif
     enddo
    end function sumSourceLuminosityMonochromatic

    real(double) function sourceLuminosityMonochromatic(source, lam) result (tot)
      type(SOURCETYPE) :: source
      real(double) :: lam, flux

      tot = 0.d0
      if (lam > source%spectrum%lambda(source%spectrum%nLambda)) then
         tot = tot + 1.d-200
      else
         flux =  loginterp_dble(source%spectrum%flux(1:source%spectrum%nLambda), &
              source%spectrum%nLambda, source%spectrum%lambda(1:source%spectrum%nLambda), lam)
         tot = tot + flux * fourPi * (source%radius*1.d10)**2 
      endif
    end function sourceLuminosityMonochromatic


    subroutine distanceToSource(source, nSource, rVec, uHat, hitSource, distance, sourcenumber)
      integer :: nSource
      type(SOURCETYPE) :: source(:)
      type(VECTOR) :: rVec, uHat
      integer, optional, intent(out) :: sourceNumber
      real(double), intent(out) :: distance
      real(double) :: cosTheta, sintheta
      logical, intent(out) :: hitSource
      integer :: i

      hitSource = .false.
      if (present(sourceNumber)) sourceNumber = 0 
      distance = 1.d30
      mainloop: do i = 1, nSource
         if (source(i)%pointsource) cycle
         cosTheta = (uHat .dot. (source(i)%position - rVec))/modulus(source(i)%position - rVec)
         if (cosTheta < 0.) then ! source is behind you
            cycle mainloop
         endif
         sinTheta = sqrt(max(0.d0,1.d0 - cosTheta*cosTheta))
         if (source(i)%radius > (modulus(source(i)%position - rVec)*sinTheta)) then
            distance = min(distanceToSphere(rVec, uHat, source(i)%position, source(i)%radius),distance)
            hitSource = .true. ! bang
            if (present(sourcenumber)) sourceNumber = i
         endif
      enddo mainloop
    end subroutine distanceToSource

    !thap
    subroutine getMelvinPositionDirection(source, position, direction, grid, photonPacketWeight)
      type(GRIDTYPE) :: grid
      real(double) :: r, photonPacketWeight
      type(SOURCETYPE) :: source
      type(VECTOR),intent(out) :: position, direction

      if (.not.source%outsideGrid) then
            position = source%position
            call melvinUnitVector(direction, photonPacketWeight)  
      else
         position%x = -grid%octreeRoot%subcellSize+1.d-10*grid%octreeRoot%subcellSize
         call randomNumberGenerator(getDouble=r)
         r = 2.d0 * r - 1.d0
         position%y = r * grid%octreeRoot%subcellSize
         call randomNumberGenerator(getDouble=r)
         r = 2.d0 * r - 1.d0
         position%z = r * grid%octreeRoot%subcellSize
         direction = VECTOR(1.d0, 0.d0, 0.d0)
      endif
    end subroutine getMelvinPositionDirection

    subroutine createSourceISRF(source, grid)
      use amr_utils_mod, only : gridArea
      type(SOURCETYPE) :: source
      type(GRIDTYPE) :: grid
      character(len=120) :: ifilename,dataDirectory
      integer :: i
      logical :: ok

      source%stellar = .false.
      source%viscosity = .false.
      source%diffuse = .true.

      call unixGetenv("TORUS_DATA", dataDirectory, i)
      ifilename = trim(dataDirectory)//"/"//"isrf.dat"

      call readSpectrum(source%spectrum, ifilename, ok)
      if (.not.ok) then
         call writeFatal("Cannot read isrf.dat")
         stop
      endif

      source%spectrum%lambda = source%spectrum%lambda * micronsToAngs
      source%spectrum%dlambda = source%spectrum%dlambda * micronsToAngs
      source%spectrum%flux = source%spectrum%flux / micronsToAngs

      source%spectrum%flux = source%spectrum%flux / 4.d0  ! from fourpi jnu to pi jnu
      
      source%luminosity = integrateSpectrumOverBand(source%spectrum,0.d0,1.d30)

      source%luminosity = source%luminosity * gridArea(grid)
    end subroutine createSourceISRF

    subroutine createSourceCMB(source, grid)
      use amr_utils_mod, only : gridArea
      type(GRIDTYPE) :: grid
      type(SOURCETYPE) :: source

      source%stellar = .false.
      source%viscosity = .false.
      source%diffuse = .true.

      call fillSpectrumBB(source%spectrum, 2.8d0, 1000.d0, 1.d9, 1000)
      source%luminosity = integrateSpectrumOverBand(source%spectrum,0.d0,1.d30)

      source%luminosity = source%luminosity * gridArea(grid)
    end subroutine createSourceCMB
      
      


    subroutine getIsotropicPositionDirection(position, direction, grid)
      use amr_utils_mod, only : distanceToGridFromOutside
      type(GRIDTYPE) :: grid
      real(double) :: z, ang, distToGrid
      type(VECTOR) :: entryPoint, norm, v, n, zAxis, thisPoint
      type(VECTOR),intent(out) :: position, direction
      logical :: ok

      zAxis = VECTOR(0.d0, 0.d0, 1.d0)
      entryPoint = randomUnitVector()
      norm = (-1.d0)*entryPoint
      ok = .false.

      do while (.not.ok)

         call randomNumberGenerator(getDouble=z)
         z = sqrt(z)
         ang = acos(z)
         n = norm .cross. zAxis
         call normalize(n)
         v = norm
         v = arbitraryRotate(v, ang, n)
         call randomNumberGenerator(getDouble=ang)
         ang = ang * twoPi
         v = arbitraryRotate(v, ang, norm)
         
         direction = v
         thisPoint = entryPoint*4d0*grid%octreeRoot%subcellSize
         
         distToGrid = distanceToGridFromOutside(grid, thisPoint, direction, hitGrid=ok)

      enddo
      position = thisPoint + (distToGrid+1.d-3*grid%halfSmallestSubcell)*direction
    end subroutine getIsotropicPositionDirection

    subroutine getViscosityPositionDirection(position, direction, grid)
      use amr_utils_mod, only : distanceToGridFromOutside
      type(GRIDTYPE) :: grid
      real(double) :: z, ang, distToGrid
      type(VECTOR) :: entryPoint, norm, v, n, zAxis, thisPoint
      type(VECTOR),intent(out) :: position, direction
      logical :: ok

      zAxis = VECTOR(0.d0, 0.d0, 1.d0)
      entryPoint = randomUnitVector()
      norm = (-1.d0)*entryPoint
      ok = .false.

      do while (.not.ok)

         call randomNumberGenerator(getDouble=z)
         z = sqrt(z)
         ang = acos(z)
         n = norm .cross. zAxis
         call normalize(n)
         v = norm
         v = arbitraryRotate(v, ang, n)
         call randomNumberGenerator(getDouble=ang)
         ang = ang * twoPi
         v = arbitraryRotate(v, ang, norm)
         
         direction = v
         thisPoint = entryPoint*4d0*grid%octreeRoot%subcellSize
         
         distToGrid = distanceToGridFromOutside(grid, thisPoint, direction, hitGrid=ok)

      enddo
      position = thisPoint + (distToGrid+1.d-3*grid%halfSmallestSubcell)*direction
    end subroutine getViscosityPositionDirection


    subroutine getPhotonPositionDirection(source, position, direction, rHat, grid, weight)
      use inputs_mod, only : biasPhiDirection, biasPhiProb, biasPhiInterval 
      use inputs_mod, only : amrgridcentrey, amrgridcentrex, hotSpot, amrGridSize
      type(GRIDTYPE) :: grid
      real(double) :: r, t, u, v, w, ang, planetBiasWeight
      real(double), optional :: weight
      type(SOURCETYPE) :: source
      type(VECTOR),intent(out) :: position, direction, rHat
      type(VECTOR) :: cornerDir
      type(VECTOR), parameter :: zAxis=VECTOR(0.d0, 0.d0, 1.d0)
!      real(double) :: prob, chance
      logical :: ok

      if (PRESENT(weight)) weight = 1.d0

      if (source%stellar) then
         if (.not.source%outsideGrid) then
            if (source%pointSource) then
               !      ! simply treating as a point source
               position = source%position
               direction = randomUnitVector()
               
               if (biasPhiDirection > 0.d0) then
                  if (.not.PRESENT(weight)) then
                     call writeFatal("getPhotonPositionDirection called without weight when weighted direction required")
                  endif
                  call randomNumberGenerator(getDouble=r)
                  w = 2.d0*r-1.d0
                  call randomNumberGenerator(getDouble=r)
                  if (r < biasPhiProb) then
                     call randomNumberGenerator(getDouble=r)
                     ang = biasPhiDirection + (r-0.5d0) * biasPhiInterval
                     weight = (biasPhiInterval/twoPi) / biasPhiProb
                  else
                     ang = biasPhiDirection
                     do while (abs(ang-biasPhiDirection) < biasPhiInterval/2.d0)
                        call randomNumberGenerator(getDouble=r)
                        ang = r * twoPi
                     enddo
                     weight = (1.d0-(biasPhiInterval/twoPi)) / (1.d0-biasPhiProb)
                  endif
                  t = sqrt(1.d0-w*w)
                  u = t*cos(ang)
                  v = t*sin(ang)
                  direction = VECTOR(u, v, w)
                  
               endif
               

!               if (present(weight)) then
!                  call randomNumberGenerator(getDouble=r)
!                  chance = 1.d0 - twoPi * (cos(80.d0*degtorad)-cos(100.d0*degtoRad))/fourPi
!                  prob = 0.01d0
!                  if (r < prob) then
!                     direction = randomUnitVector()
!                     ang = acos(direction.dot.zAxis)*radtodeg
!                     do while ((ang >= 80.d0).and.(ang <= 100.d0))
!                        direction = randomUnitVector()
!                        ang = acos(direction.dot.zAxis)*radtodeg
!                     enddo
!                     weight = chance/prob
!                  else
!                     direction = randomUnitVector()
!                     ang = acos(direction.dot.zAxis)*radtodeg
!                     do while ((ang <= 80.d0).or.(ang >= 100.d0))
!                        direction = randomUnitVector()
!                        ang = acos(direction.dot.zAxis)*radtodeg
!                     enddo
!                     weight = (1.d0-chance)/(1.d0-prob)
!                  endif
!                  ang = acos(direction.dot.zAxis)*radtodeg
!               endif

               
            else
               rHat = randomUnitVector()
               position = source%position + source%radius*rHat
               ! A limb darkening law should be applied here for 
               ! for general case here.....
               direction = fromPhotoSphereVector(rHat)
               if (hotSpot) then
                  call getPhotoVec(source%surface, position, direction)
               endif
               
            endif
         !Thaw- if source is on a corner, direct all photons into computational domain
         !These are re-weighted to account for bias in photoionAMR_mod
         if (source%onCorner) then
            cornerDir = source%position - grid%octreeRoot%centre
            if(cornerDir%x == 0.0 .and. cornerDir%y == 0.0 .and. cornerDir%z == 0.0 ) then
               print *, "Corner Error! "
               stop
            end if
            !2D    print *, "cornerDir = ", cornerDir
            ! print *, "source%position = ", source%position
            ! print *, "grid%octreeRoot%centre", grid%octreeRoot%centre
            ! print *, "preDirection = ", direction
            if (cornerDir%x > 0 .and. direction%x > 0) then
               direction%x = -direction%x 
               ! print *, "called A"
            else if (cornerDir%x < 0.d0 .and. direction%x < 0.d0) then
               direction%x = -direction%x
               ! print *, "called B"
            end if
            if (cornerDir%z > 0.d0 .and. direction%z > 0.d0) then
               direction%z = -direction%z
               !print *, "called C"
            else if (cornerDir%z < 0.d0 .and. direction%z < 0.d0) then
               direction%z = -direction%z
               !  print *, "called D"
            end if
            if (cornerDir%y > 0.d0 .and. direction%y > 0.d0) then
               direction%y = -direction%y
               ! print *, "called E"
            else if (cornerDir%y < 0.d0 .and. direction%y < 0.d0) then
               direction%y = -direction%y
               !! print *, "called F"
            end if
            !Currently dont handle edges
            
            !            print *, "postDirection : ", direction
            
            
            !direction%z = abs(direction%z)
            !direction%y = abs(direction%y)
            !direction%x = abs(direction%x)
         end if
      else
         position%x = amrgridcentrex-grid%octreeRoot%subcellSize-(1.d-10*grid%halfsmallestsubcell)
         call randomNumberGenerator(getDouble=r)
         r = 2.d0 * r - 1.d0
         position%y = r * grid%octreeRoot%subcellSize*2.d0
         call randomNumberGenerator(getDouble=r)
         r = 2.d0 * r - 1.d0
         position%z = r * grid%octreeRoot%subcellSize*2.d0
!         direction = VECTOR(1.d0, 0.d0, 0.d0)
         direction = VECTOR(1.d0, 0.d0, 0.d0)
         if(grid%octreeroot%twod) then
            position%y = amrgridcentrey
            if(grid%geometry == "SB_gasmix") then
               ok = .true.
!               if(position%z < 2.6d-10*pctocm .and. position%z > 1.4d-10*pctocm) then                  
!               if(position%z < 2.625d-10*pctocm .and. position%z > 1.375d-10*pctocm) then              
               if(position%z < (6.d8 + (10.d0*6.d8/32.d0)) .and. position%z > &
                    (6.d8 - (10.d0*6.d8/32.d0))) then
!               if(position%z < 2.75d-10*pctocm .and. position%z > 4.312d8) then              
                  !In the shadow
                  ok = .false.
               end if

               if(.not. ok) then
                  do while (.not. ok)
                     call randomNumberGenerator(getDouble=r)
                     r = 2.d0 * r - 1.d0
                     position%z = r * grid%octreeRoot%subcellSize*2.d0
!                     if(position%z > 2.6d-10*pctocm .or. position%z < 1.4d-10*pctocm) then                  
!                     if(position%z > 2.75d-10*pctocm .or. position%z < 1.25d-10*pctocm) then                  
                     if(position%z > (6.d8 + (10.d0*6.d8/32.d0)) .or. position%z < &
                          (6.d8 - (10.d0*6.d8/32.d0))) then
                        !In the shadow
                        ok = .true.
                     end if
                  end do
               end if
            end if
         else if(grid%octreeroot%oned) then

               rHat = randomUnitVector()
               position = source%position + source%radius*rHat
               ! A limb darkening law should be applied here for 
               ! for general case here.....
               direction = fromPhotoSphereVector(rHat)
               if (hotSpot) then
                  call getPhotoVec(source%surface, position, direction)
               endif
               call biasTowardsPlanet(position, rHat, VECTOR(0.d0, 0.d0, 0.d0), dble(amrGridSize), &
                    planetBiasweight, direction)
               if (PRESENT(weight)) weight = weight * planetBiasWeight

         end if
      endif
      !Thaw - photons from corner sources should be directed into the grid

   else
      call getIsotropicPositionDirection(position, direction, grid)
   endif
        
    end subroutine getPhotonPositionDirection

  function distanceToSphere(rVec, uHat, centre, radius) result (distance)
    real(double) :: distance, radius
    real(double) :: x1, x2,a,b,c,d,costheta
    logical :: ok
    type(VECTOR) :: rVec, uHat, centre

    d = modulus(centre - rVec)
    cosTheta = (uHat .dot. (centre-rVec)) / d

    a = 1.d0
    b = -2.d0*d*cosTheta
    c = d**2 - radius**2
    call solveQuadDble(a, b, c, x1, x2, ok)
    if (.not.ok) then 
       write(*,*) "! Distance to sphere called with no intersection"
    endif
    if ((x1 > 0.).and.(x2 > 0.)) then
       distance = min(x1, x2)
    else
       if (x1 > 0.) then
          distance = x1
       else
          distance = x2
       endif
    endif
    if (distance < 0.) then
       write(*,*) "! screw up in distance to sphere"
    endif
  end function distanceToSphere




    !
    ! For a given octal and sourcetype objects, this routine checks 
    ! if the source is with in the range this octal.
    
    function source_within_octal(this, an_octal) RESULT(out)
      implicit none
      logical :: out 
      type(sourcetype), intent(in) :: this
      type(octal), intent(in) :: an_octal 
      
      real(oct) :: x, y, z    ! position of the source.
      
      x = this%position%x
      y = this%position%y
      z = this%position%z
 
      out = within_subcell(an_octal, 0, x, y, z)
 
  end function source_within_octal




  !
  ! Pick a random direction (with length 1) from the surface of a sphere. 
  ! If the direction points inward of the sphere, it will be rejected, and a new one will 
  ! be calculated. 
  !
  ! note: rvec is a vector which connects a point on the surface of the sphere and its origin.
  !       The length of the rvec does not matther when passed to this routine as long as 
  !        it has a correct direction, this routine should work.
  function random_direction_from_sphere(rvec) RESULT(out)
    implicit none
    type(VECTOR) :: out 
    !
    type(VECTOR), intent(in) :: rvec
    !
    type(VECTOR) :: dir 
    real(oct) :: inner_product
    
    inner_product = -1.0d0
    do while (inner_product<0.0d0) 
       dir = randomUnitVECTOR()
       inner_product = dir .dot. rvec
    end do
    
    out = dir

  end function random_direction_from_sphere

  real(double) function I_nu(source, nu, iElement, mu) 
    type(SOURCETYPE) :: source
    real(double) :: nu, fnu, mu !, flambda, lam
    integer :: i, iElement
    real(double) :: tAccretion, ic_hot, fac

!    lam = 1.d8 * cSpeed/ nu ! angs
!    if (lam < source%spectrum%lambda(1)) then
!       I_nu = tiny(i_nu)
!    endif
!    if (lam > source%spectrum%lambda(source%spectrum%nlambda)) then
!       I_nu = tiny(i_nu)
!    endif


!    call locate(source%spectrum%lambda, source%spectrum%nLambda, lam, i)
!    fLambda = source%spectrum%flux(i)
!    fnu = flambda * (cSpeed*1.d8) /nu**2 ! 1.d8 to go from cm/s to angs/c
!    i_nu = 1.5d0 * fnu / pi

!    write(*,*) "i_nu ",nu, source%surface%nuArray(1), &
!         source%surface%nuArray(source%surface%nNuHotFlux)

    if (nu < source%surface%nuArray(1)) then
       i_nu = tiny(i_nu)
!       write(*,*) nu,source%surface%nuArray(1)
    else if (nu > source%surface%nuArray(source%surface%nNuHotFlux)) then
       i_nu = tiny(i_nu)
!       write(*,*) nu,source%surface%nuArray(source%surface%nNuHotFlux)
    else
       call locate(source%surface%nuArray, source%surface%nNuHotFlux, nu, i)
       fnu = logint(nu, source%surface%nuArray(i), source%surface%nuArray(i+1), &
            source%surface%hnuArray(i), source%surface%hnuArray(i+1))
       i_nu = fnu / pi
    endif

    fac = (2.d0/12.d0)*(6.d0 - 2.d0*source%limbDark(1) - source%limbDark(2))

    i_nu = i_nu * (1.d0 - source%limbDark(1)*(1.d0 - mu) - source%limbDark(2) * (1.d0-mu)**2)/fac

    if (isHot(source%surface,ielement)) then
       tAccretion = source%surface%element(ielement)%temperature
       !================CHECK UNITS HERE!! ===========================
       IC_hot = blackbody(REAL(tAccretion), 1.e8*REAL(cSpeed/nu))! [B_nu]
!       IC_hot = bNu(nu, dble(tAccretion))
!       write(*,*) "hit hotspot",ic_hot,i_nu
       i_nu  = I_nu + ic_hot
    endif


  end function I_nu

  logical function insideSource(thisOctal, subcell, nsource, source)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    type(VECTOR) :: rVec, corner
    integer :: nsource
    type(SOURCETYPE) :: source(:)
    integer :: i
    real(double) :: r

    insideSource = .false.
    rvec = subcellCentre(thisOctal, subcell)

    do i = 1, nSource
       
       if (thisOctal%twoD) then
          corner = rVec + thisOctal%subcellSize/rootTwo * VECTOR(oneOverrootTwo,0.d0,oneOverrootTwo)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/rootTwo * VECTOR(oneOverrootTwo,0.d0,-oneOverrootTwo)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/rootTwo * VECTOR(-oneOverrootTwo,0.d0,oneOverrootTwo)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/rootTwo * VECTOR(-oneOverrootTwo, 0.d0, -oneOverrootTwo)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
       else if (thisOctal%oneD) then
          corner = rVec - VECTOR(thisOctal%subcellSize/2.d0,0.d0,0.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
       else
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(1.d0, 1.d0, 1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(1.d0, -1.d0, 1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(1.d0, 1.d0, -1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(1.d0, -1.d0, -1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(-1.d0, -1.d0, -1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(-1.d0, 1.d0, 1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(-1.d0, 1.d0, -1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(-1.d0, -1.d0, 1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif

       endif
    enddo
  end function insideSource



  end module source_mod
