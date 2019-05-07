module starburst_mod

! A module to create sources, and their distribution, for
! starburst dust/photoionization calculations

! started by tjh on 23/8/2006



  use kind_mod
  use citations_mod
  use constants_mod
  use messages_mod
  use utils_mod
  use source_mod
  use unix_mod, only: unixGetenv

  implicit none

  type TRACKTABLE
     character(len=80) :: label
     integer :: nMass
     real(double), pointer :: initialMass(:) ! Msol
     integer, pointer :: nAges(:)
     real(double), pointer :: age(:,:) ! yr
     real(double), pointer :: massAtAge(:,:) ! Msol
     real(double), pointer :: logL(:,:) ! log(Lsol)
     real(double), pointer :: logTeff(:,:) ! log(K)
     real(double), pointer :: mDot(:,:) ! log(Msol/yr)
     real(double), pointer :: logRadius(:,:) ! log(Rsol)
  end type TRACKTABLE

contains

  function randomMassFromIMF(imfType, minMass, maxMass, alpha) result (mass)

    ! returns a stellar mass (in solar masses) drawn from an IMF

    character(len=*) :: imfType
    character(len=80) :: message
    real(double) :: minMass, maxMass, mass, alpha
    integer, parameter :: nMass = 100
    integer :: i
    real(double) :: massArray(nMass), prob(nMass), dMass(nMass), r, t
    
    do i = 1, nMass
       massArray(i) = log10(minMass) + (dble(i-1)/dble(nMass-1))*(log10(maxMass)-log10(minMass))
       massArray(i) = 10.d0**massArray(i)
    enddo

    do i = 2, nMass-1
       dMass(i) = massArray(i+1)-massArray(i)
    enddo
    dMass(1) = 2.d0*(massArray(2)-massArray(1))
    dMass(nMass) = 2.d0*(massArray(nMass)-massArray(nMass-1))
       

    select case (imfType)
    
      case("salpeter")
         do i = 1, nMass
            prob(i) = (massArray(i)**alpha) * dMass(i)
         enddo

      ! todo
!     case("chabrier")      
!     case("kroupa")      

      case DEFAULT
         write(message,'(a,a)') "IMF not recognised: ", trim(imfType)
         call writefatal(message)

    end select


    do i = 2, nMass
       prob(i) = prob(i) + prob(i-1)
    enddo
    prob(1:nMass) = prob(1:nMass) - prob(1)
    prob(1:nMass) = prob(1:nMass) / prob(nMass)
    call randomNumberGenerator(getDouble=r)
    call locate(prob, nMass, r, i)
    t = (r - prob(i))/(prob(i+1)-prob(i))
    mass = log10(massArray(i)) + t * (log10(massArray(i+1))-log10(massArray(i)))
    mass = 10.d0**mass
  end function randomMassFromIMF

  subroutine populateClusters(clusters, nClusters, age) 
    use inputs_mod, only : criticalMass ! [g]
    type(SOURCETYPE), pointer :: clusters(:)
    real(double) :: totalCreatedMass ! [msol]
    real(double) :: reservoir  ! [g]
    real(double) :: age ! [yr]
    integer :: nClusters, i

    do i = 1, nClusters
       reservoir = clusterReservoir(clusters(i))
       if (reservoir >= criticalMass) then
          if (writeoutput) write(*,*) "Creating subsources for cluster ", i
          call createSources(clusters(i)%nSubsource, clusters(i)%subsourceArray, "instantaneous", age, reservoir/msol, 0.d0, &
                totalCreatedMass, zeroNsource=.false.)
          if (writeoutput) write(*,*) " ... created ", totalCreatedMass, " Msol" 
!          if (totalCreatedMass > 0.d0) write(*,*) myrankglobal, " CREATED ", clusters(i)%subsourceArray(:)%mass/msol

          clusters(i)%subsourceArray(1:clusters(i)%nSubsource)%position = clusters(i)%position
          clusters(i)%subsourceArray(1:clusters(i)%nSubsource)%velocity = clusters(i)%velocity

          call writeClusterIMF(clusters, nClusters)
       endif
    enddo


  end subroutine populateClusters

  subroutine writeClusterIMF(clusters, nClusters) 
    type(SOURCETYPE), pointer :: clusters(:)
    integer :: nClusters, i, j
    character(len=80) :: fn

    if (writeoutput) then
       do i = 1, nClusters
          if (clusters(i)%nSubsource > 0) then
             write(fn, '(a,i4.4,a)') "imf_", i, ".dat"
             open(67,file=trim(fn),status="replace",form="formatted")
             write(67, '(a4,a9,a12)') "#i", "M(msol)", "age(yr)"
             do j = 1, clusters(i)%nSubsource
                write(67, '(i4.4,f9.3,es12.5)') j, clusters(i)%subsourceArray(j)%mass/msol, clusters(i)%subsourceArray(j)%age
             enddo
             close(67)
          endif
       enddo
    endif
  end subroutine writeClusterIMF
  
  ! calculate cluster spectrum from already calculated subsource spectra
  subroutine setClusterSpectra(clusters, nClusters)
    type(SOURCETYPE), pointer :: clusters(:)
    integer :: i, j, nClusters

    do i = 1, nClusters
       if (clusters(i)%nSubsource > 0) then
          ! update subsource spectra
          call setSourceSpectra(clusters(i)%subsourceArray, clusters(i)%nSubsource)
          ! now do cluster spectrum (weighted sum of subsource fluxes)
          call freeSpectrum(clusters(i)%spectrum)
          call copySpectrum(clusters(i)%spectrum, clusters(i)%subsourceArray(1)%spectrum)
          clusters(i)%spectrum%flux = clusters(i)%spectrum%flux * clusters(i)%subsourceArray(1)%luminosity
          if (clusters(i)%nSubsource > 1) then
             do j = 2, clusters(i)%nSubsource
                call addSpectrum(clusters(i)%spectrum, clusters(i)%subsourceArray(j)%spectrum,&
                   clusters(i)%subsourceArray(j)%luminosity)
             enddo
          endif
          clusters(i)%luminosity = sum(clusters(i)%subsourceArray(1:clusters(i)%nsubsource)%luminosity)
          clusters(i)%spectrum%flux = clusters(i)%spectrum%flux / clusters(i)%luminosity ! divide by sum of weights

          call probSpectrum(clusters(i)%spectrum)
          call normalizedSpectrum(clusters(i)%spectrum)

          ! effective radius, Teff given this luminosity and spectrum
          clusters(i)%radius = sqrt(clusters(i)%luminosity/(fourPi * integrateSpectrumOverBand(clusters(i)%spectrum, 10.d0, 1.d7)))
          clusters(i)%radius = clusters(i)%radius/1.d10
          clusters(i)%teff = (clusters(i)%luminosity/(fourPi * (clusters(i)%radius*1.d10)**2 * stefanBoltz))**0.25

       else
          ! zero flux
          call freeSpectrum(clusters(i)%spectrum)
          call newSpectrum(clusters(i)%spectrum, 100.d0, 1.d7, 1000)
          clusters(i)%luminosity = 0.d0
       endif
    enddo
  end subroutine setClusterSpectra

  subroutine setSourceSpectra(sources, nSource)
    use inputs_mod, only : clusterSinks
    type(SOURCETYPE) :: sources(:) 
    integer :: i, nSource

    do i = 1, nSource
      if (clustersinks) then
         ! TODO use tlusty/kurucz spectra (but currently easier to add BB spectra together)
         call fillSpectrumBB(sources(i)%spectrum, sources(i)%teff, 100.d0, 1.d7, 1000)
      else
         ! update spectrum. If tlusty spectrum is not found for a source, kurucz spectrum is used instead
         call fillSpectrumTlusty(sources(i)%spectrum, sources(i)%teff, sources(i)%mass, sources(i)%radius*1.d10)
      endif
    enddo
  end subroutine setSourceSpectra

  real(double) function clusterReservoir(cluster) ! [g]
    type(SOURCETYPE) :: cluster
    if (cluster%nsubsource > 0) then
       clusterReservoir = cluster%mass - sum(cluster%subsourceArray(1:cluster%nSubsource)%mass) 
    else
       clusterReservoir = cluster%mass
    endif
  end function clusterReservoir

  subroutine createSources(nSource, source, burstType, burstAge, burstMass, sfRate, totMass,zeroNsource)
    use inputs_mod, only : imfType, clusterSinks
    integer :: nSource, thisNsource, initialNsource
    type(SOURCETYPE), pointer :: source(:)
    type(SOURCETYPE), allocatable :: tempSourceArray(:)
    type(TRACKTABLE) :: thisTable
    character(len=80) :: message
    character(len=*) :: burstType
    real(double) :: burstAge  ! yr
    real(double) :: burstMass ! Msol
    real(double) :: thisMass, totMass    ! Msol
    real(double) :: sfRate    ! Msol/yr
    real(double), allocatable :: initialMasses(:), temp(:) ! Msol
    logical, optional :: zeroNsource
    integer :: i, j
    integer :: nDead, nSupernova, nOB
!    integer, parameter :: nKurucz=69 !nKurucz = 410
    logical :: thirtyFound, converged
!    type(SPECTRUMTYPE) :: kSpectrum(nKurucz)
!    character(len=80) :: klabel(nKurucz)
    character(len=80) :: filename

   
!    call  readKuruczGrid(klabel, kspectrum, nKurucz)
!    call  readTlustyGrid(klabel, kspectrum, nKurucz)

    if (present(zeroNsource)) then
       ! may want to keep nsource as it is (e.g. if creating more subsources in a cluster)
       if (zeroNsource) then
          nSource = 0
       endif
    else
       nSource = 0
    endif

    thisNsource = 0
    initialNsource = nSource


    ! set up the initial number of stars and their masses and ages

    select case(burstType)

       case("continuous")
          burstMass = sfRate * burstAge
          totMass = 0.d0
          do while (totMass < burstMass)
             nSource = nSource + 1
             call randomNumberGenerator(getDouble=source(nSource)%age)
             source(nSource)%age = source(nSource)%age * burstAge
             source(nSource)%initialmass = randomMassFromIMF(imfType, 0.8d0, 120.d0, -2.35d0)
             totMass = totMass + source(nSource)%initialmass
          enddo

       case("instantaneous")
          allocate(initialMasses(1:1000))
          thirtyFound = .false.
          do while(.not.thirtyFound)
             initialMasses = 0.d0
             totMass = 0.d0
             thisNsource = 0
             converged = .false.
             do while (.not. converged)
                thisMass = randomMassFromIMF(imfType, 0.8d0, 120.d0, -2.35d0)
                if ((thisMass + totMass) > burstMass) then
!                   if (writeoutput) write(*,*) "CONVERGED ", thisMass, totMass, burstMass, thirtyFound, initialMasses(1:thisNsource)
                   converged = .true.
                else
                   converged = .false.
                   thisNsource = thisNsource + 1
                   initialMasses(thisNsource) = thisMass 
                   if (thisMass >= 30.d0) thirtyFound = .true.
                   totMass = totMass + thisMass 
                endif
             enddo
          enddo
          allocate(temp(1:thisNsource))
          temp = initialMasses(1:thisNsource) 
          call sort(thisNsource, temp) ! sort in ascending order
          temp = temp(thisNsource:1:-1) ! reverse, i.e. sort in descending order

          if (writeoutput) write(*,*) "TEMP SOURCE ", temp
          if (.not.associated(source)) then
             allocate(source(1:thisNsource))
             if (writeoutput) write(*,*) "ALLOCATING ", size(source)
          endif
          ! extend source array if necessary
          if ((initialnSource+thisNsource) > size(source)) then
             if (writeoutput) write(*,*) "RESIZING (A) ", initialnsource, thisNsource, size(source)
             allocate(tempSourceArray(1:initialnSource))
             tempSourceArray(1:initialnSource) = source(1:initialnSource)
             deallocate(source)
             allocate(source(1:initialnSource+thisNsource))
             source(1:initialnSource) = tempSourceArray(1:initialnSource)
             deallocate(tempSourceArray)
             if (writeoutput) write(*,*) "RESIZED ", size(source)
          endif

          ! update nSource
          nSource = nSource + thisNsource
          if (writeoutput) write(*,*) "NEW NSOURCE", nSource, thisNsource
          
          ! add to actual source array
          source(initialNsource+1:nSource)%initialMass = temp(1:thisNsource) 
          source(initialNsource+1:nSource)%age = burstAge

          if (writeoutput) write(*,*) "NEW SOURCE ARRAY ", source(:)%initialMass
          deallocate(initialMasses)
          deallocate(temp)
          
       case("supernovatest")
             nSource = 1
             source(1)%initialMass = 40.d0
             totMass = 40.d0
             source(1:nSource)%age = burstAge
       
       case("singlestartest")
             nSource = 1
             source(1)%initialMass = burstMass
             totmass = burstmass
             source(1)%age = burstAge

       case DEFAULT
         write(message,'(a,a)') "Burst type not recognised: ", trim(burstType)
         call writefatal(message)
             
      end select
      
      source(initialnSource+1:nSource)%nSubsource = 0

      if (Writeoutput) then
         write(*,*) "number of sources in this burst ", nSource-initialNsource
         write(*,*) "burst mass ",totMass
         write(*,*) "using ", imfType, " imf"
      endif

      ! now get actual masses, temps, and luminosities, and radii for age from evolution tracks

      call readinTracks("mist", thisTable)
      call writeInfo("MIST tracks successfully read", FORINFO)
      nDead = 0
      nSupernova = 0
      i = initialnSource+1
      do while (i <= nSource)
         if (source(i)%nSubsource == 0) then 
            if (.not.isSourceDead(source(i), thisTable)) then
!               write(message, '(a,i4)') "Setting properties of source ", i
!               call writeInfo(message, TRIVIAL)
               call setSourceProperties(source(i))
               i = i + 1
            else
               write(message, '(a,i4)') "Removing source ", i
               call writeInfo(message, TRIVIAL)
               if (clusterSinks) then
!                  i = i + 1 ! ignore
                  ! TODO make this work for clustersinks
                  call torus_stop("todo")
               else
                  nDead = nDead + 1
                  if (source(i)%initialMass > 8.d0) then
                     nSupernova = nSupernova + 1
                  endif
                  call removeSource(source, nSource, i)
               endif
            endif
         endif
      enddo
      if (Writeoutput) then
         write(*,*) "Number of dead sources",  nDead
         write(*,*) "Number of supernova",  nSupernova
         write(*,*) "Burst luminosity",SUM(source(initialnSource+1:nSource)%luminosity)/lsol
      endif

      nOB = 0
      do i = initialnSource+1, nSource
         if (source(i)%mass > 15.*msol) then
            nOB = nOB + 1
         endif
      enddo
      if (writeoutput) then
         write(*,*) "Number of OB stars (>15 msol): ",nOB
      endif
!      do i = 1, nSource
!         ! use tlusty spectrum - if it's not found for a source, kurucz spectrum is used instead
!         call fillSpectrumTlusty(source(i)%spectrum, source(i)%teff, source(i)%mass, source(i)%radius*1.d10)
!      enddo

      if (.not. clusterSinks) then
         if (writeoutput) then
            do i = initialnSource+1, nSource
               write(filename,'(a, i3.3, a)') "spectrum_source", i,".dat"
               open(67,file=filename,status="replace",form="formatted")
               do j = 1, source(i)%spectrum%nlambda
                  write(67,'(2es12.5)') source(i)%spectrum%lambda(j), source(i)%spectrum%flux(j)
               enddo
               close(67)
            enddo
         endif
      endif

!      do i = initialnSource+1, nSource
!         if (source(i)%mass/msol > 15.d0) then
!            source(i)%mDotWind = 1.d-6 * msol / (365.25d0 * 24.d0 * 3600.d0)
!         else
!            source(i)%mdotWind = 0.d0
!         endif
!      enddo

      source(initialnSource+1:nSource)%stellar = .true.
      source(initialnSource+1:nSource)%viscosity = .false.
      source(initialnSource+1:nSource)%diffuse = .false.
      source(initialnSource+1:nSource)%outsideGrid = .false.
      source(initialnSource+1:nSource)%prob = 0.d0 ! 1.d0/dble(nsource)
      call writeInfo("Photons will be sampled according to source luminosity", TRIVIAL)
    end subroutine createSources

    subroutine dumpSources(source, nsource, label)
      type(SOURCETYPE) :: source(:)
      integer :: nSource
      integer :: i
      integer, optional :: label
      character(len=80) :: filename

      if (writeoutput) then
         if (present(label)) then
            write(filename,'(a, i4.4, a)') "starburst_", label, ".dat"
            open(32, file=filename, form="formatted", status="unknown")
         else
            open(32, file="starburst.dat", form="formatted", status="unknown")
         endif
         write(32,'(a)') "    #      mass    teff  radius  luminosity    position "
         write(32,'(a)') "    #    (Msol)     (K)  (Rsol)      (Lsol)   (10^10cm)"
         do i = 1, nSource
            write(32, '(i5,  f10.3, i8, f8.1, 1pe12.2, 1p, 3e12.2)') i, source(i)%mass/msol, nint(source(i)%teff), &
                 source(i)%radius*1.d10/rsol, source(i)%luminosity/lsol, source(i)%position
         enddo
         close(32)
      endif
    end subroutine dumpSources
         


    logical function isSourceDead(source, thisTable) result (dead)
      type(SOURCETYPE) :: source
      type(TRACKTABLE) :: thisTable
      real(double) :: t, t1, t2, deadAge
      integer :: i
      call locate(thisTable%initialMass, thisTable%nMass, source%initialMass, i)
      t1 = thisTable%age(i,thisTable%nAges(i))
      t2 = thisTable%age(i+1,thisTable%nAges(i+1))

      t = (source%initialMass - thisTable%initialMass(i))/(thisTable%initialMass(i+1) - thisTable%initialMass(i))
      deadAge = t1 + t * (t2 - t1)
      if (source%age > deadAge) then
         dead = .true.
      else
         dead = .false.
      endif
    end function isSourceDead


    subroutine checkSourceSupernova(nSource, source, nSupernova, supernovaIndex, ejectaMass, ke)

      ! checks if any sources have gone supernova, looping through all sources                                                            
      ! requires sources to be dead and initially above 8 solar masses                                                                    
      ! counts number gone supernova, tabulates each supernova source index                                                               

      type(SOURCETYPE) :: source(:)
      type(TRACKTABLE),save :: thisTable
      logical,save :: firstTime = .true.
      real(double) :: t, t1, t2, deadAge, remnantMin1, remnantMin2, absMetallicity
      integer :: i, j
      integer :: nSource, nSupernova, supernovaIndex(:)
      real(double) :: ejectaMass(:), ke(:)


      if (firstTime) then
         call readinTracks("mist", thisTable)
         firstTime = .false.
      endif


      nSupernova = 0

      do i=1, nSource

         ! locates sources lower mass bound evolutionary track, reads
         ! and assigns death ages of upper/lower tracks

         call locate(thisTable%initialMass, thisTable%nMass, source(i)%initialMass, j)
         t1 = thisTable%age(j,thisTable%nAges(j))
         t2 = thisTable%age(j+1,thisTable%nAges(j+1))

         ! interpolates death age for mass inbetween 2 consecutive evolutionary tracks                                                    

         t = (source(i)%initialMass - thisTable%initialMass(j))/(thisTable%initialMass(j+1) - thisTable%initialMass(j))
         deadAge = t1 + t * (t2 - t1)
         if (source(i)%initialMass > 8.d0 .and. writeoutput) then
           write(*,'(a, i4, a, 1pe12.5)') "Source ", i, " Time until supernova (yr): ", deadAge-source(i)%age
         endif
         ! checks if each source is dead and initially > 8 solar mass, adds to SN count, tabulates index                                  


         if (source(i)%age > deadAge .and. source(i)%initialMass > 8.d0) then
         if (writeoutput) write(*,*) "Source " ,i, " explodes as a supernova"
            nSupernova=nSupernova+1

         ! ejecta mass algorithm for delayed supernova model, calculates ejecta mass based on remnant mass using various
         ! equations dependant on ZAMS (initial stellar mass) range. Model includes metallicity dependance, assumed to be
         ! = 1 (ratio to solar for algorithm) within torus currently (14/11/2014).

         absMetallicity = 1     ! ratio to solar metallicity

            if (source(i)%initialMass .lt. 11.d0) then     ! mstar < 11
                ejectaMass(nSupernova) = source(i)%mass - (1.28 + ((source(i)%initialMass - 8.d0)/3.d0) * 0.08) * msol
            else if (source(i)%initialMass .ge. 11.d0 .and. source(i)%initialMass .lt. 30.d0) then     ! 11 <= mstar < 30
                ejectaMass(nSupernova) = source(i)%mass - (1.1 + 0.2 * exp((source(i)%initialMass - 11.d0)/4.d0) - &
                   (2.0 + absMetallicity) * exp(0.4 * (source(i)%initialMass-26.d0))) * msol
            else if (source(i)%initialMass .ge. 30.d0 .and. source(i)%initialMass .lt. 50.d0) then     ! 30 <= mstar < 50
                remnantmin1 = 33.35 + (4.75 + 1.25 * absMetallicity) * (source(i)%initialMass-34d0)
                remnantmin2 = source(i)%initialMass - ((absMetallicity)**0.5) * (1.3 * source(i)%initialMass - 18.35)
                ejectaMass(nSupernova) = source(i)%mass - min(remnantmin1, remnantmin2) * msol
            else if (source(i)%initialMass .ge. 50.d0 .and. source(i)%initialMass .lt. 90.d0) then     ! 50 <= mstar < 90
                ejectaMass(nSupernova) = source(i)%mass - (1.8 + 0.04 * (90.d0 - source(i)%initialMass)) * msol
            else if (source(i)%initialMass .ge. 90.d0) then                                            ! mstar >= 90
                ejectaMass(nSupernova) = source(i)%mass - (1.8 + log10(source(i)%initialMass - 89.d0)) * msol
            end if

            supernovaIndex(nSupernova)=i

         ! Ejecta energy algorithm for delayed supernova model, calculates kinetic energy based on each SN sources
         ! ZAMS (initial stellar mass) range. Based on computational model ("Progenitor explosion connection.."
         ! Ugliano et al. 2012) with similar evolution properties to Schaller. Model assumed to be unitary in
         ! SN probability and continuous in energy distribution, energy extrapolated within range (16/01/2015)

            if (source(i)%initialMass .lt. 14.d0) then                                ! mstar < 14
                ke(nSupernova) = (1.5 + ((source(i)%initialMass - 8.d0)/6.d0)*0.25)*1.d51
            else if (source(i)%initialMass .lt. 16.2 .and. source(i)%initialMass .ge. 14.d0) then     ! 14 <= mstar < 16.2
                ke(nSupernova) = (1.75 - ((source(i)%initialMass - 14.d0)/2.2)*1.05)*1.d51
            else if (source(i)%initialMass .lt. 17.d0 .and. source(i)%initialMass .ge. 16.2) then     ! 16.2 <= mstar < 17
                ke(nSupernova) = (0.7 + ((source(i)%initialMass - 16.2)/0.8)*0.8)*1.d51
            else if (source(i)%initialMass .lt. 18.5 .and. source(i)%initialMass .ge. 17.d0) then     ! 17 <= mstar < 18.5
                ke(nSupernova) = (1.5 - ((source(i)%initialMass - 17.d0)/1.5)*0.6)*1.d51
            else if (source(i)%initialMass .lt. 19.5 .and. source(i)%initialMass .ge. 18.5) then      ! 18.5 <= mstar < 19.5
                ke(nSupernova) = (0.9 + ((source(i)%initialMass - 18.5)/1)*0.6)*1.d51
            else if (source(i)%initialMass .lt. 25 .and. source(i)%initialMass .ge. 19.5) then        ! 19.5 <= mstar < 25
                ke(nSupernova) = (1.5 - ((source(i)%initialMass - 19.5)/5.5)*0.6)*1.d51
            else if (source(i)%initialMass .lt. 26 .and. source(i)%initialMass .ge. 25.d0) then       ! 25 <= mstar < 26
                ke(nSupernova) = (0.9 + ((source(i)%initialMass - 25.d0)/1)*0.6)*1.d51
            else if (source(i)%initialMass .lt. 31.d0 .and. source(i)%initialMass .ge. 26.d0) then    ! 26 <= mstar < 31
                ke(nSupernova) = (1.5 - ((source(i)%initialMass - 26.d0)/5)*0.3)*1.d51
            else if (source(i)%initialMass .lt. 35.5 .and. source(i)%initialMass .ge. 31.d0) then     ! 31 <= mstar < 35.5
                ke(nSupernova) = (1.2 + ((source(i)%initialMass - 31.d0)/4.5)*0.3)*1.d51
            else if (source(i)%initialMass .lt. 40.d0 .and. source(i)%initialMass .ge. 35.5) then     ! 35.5 <= mstar <= 40
                ke(nSupernova) = (1.5 - ((source(i)%initialMass - 35.5)/4.5)*0.3)*1.d51
            else if (source(i)%initialMass .lt. 90.d0 .and. source(i)%initialMass .ge. 40.d0) then    ! 40 <= mstar < 90
                ke(nSupernova) = (1.2 - ((source(i)%initialMass - 40.d0)/50)*0.7)*1.d51
            else if (source(i)%initialMass .ge. 90.d0) then                                          ! 90 <= mstar
                ke(nSupernova) = 0.5d51
            endif

         endif
      end do
    end subroutine checkSourceSupernova



    subroutine setSourceProperties(source)
      use inputs_mod, only : mStarburst, clusterRadius, accretionRadius, smallestCellSize, burstType, burstPosition
      type(SOURCETYPE) :: source
      real(double) :: r
      type(VECTOR) :: vVec
      real(double) :: sigmaVel

      ! set source mass, age, luminosity, Teff, radius, mass-loss rate, spectrum
      call updateSourceProperties(source)

      source%position = VECTOR(0.d0, 0.d0, 0.d0)
      source%velocity = VECTOR(0.d0, 0.d0, 0.d0)
      if (clusterRadius > 0.d0) then
         source%position = randomUnitVector()
         call randomNumberGenerator(getDouble=r)
         r = r**2
         source%position = source%position * (clusterRadius / 1.d10) * r
         sigmaVel = sqrt(bigG * ((Mstarburst+1000.d0)*mSol)/(2.d0*clusterRadius))
   !      if (writeoutput) write(*,*) "Sigma velocity ",sigmaVel/1.e5
         vVec = randomUnitVector()
         r = gasdev()
         source%velocity = r * sigmaVel * vVec
      endif
      source%accretionRadius = accretionRadius*smallestCellsize*1.d10

      select case(burstType)
         case("supernovatest")
            source%position = VECTOR(0.d0, 0.d0, 0.25d0*clusterRadius/1.d10)
            source%velocity = VECTOR(0.d0, 0.d0, 0.d0)
         case("singlestartest")
            source%position = burstPosition 
            source%velocity = VECTOR(0.d0, 0.d0, 0.d0)
         case DEFAULT
      end select

    end subroutine setSourceProperties

    subroutine updateSourceProperties(source)
      type(SOURCETYPE) :: source
      type(TRACKTABLE),save :: thisTable
      logical,save :: firstTime = .true.
      integer :: i, j
      real(double) :: t, u, mass1, logL1, logT1, logmdot1
      real(double) :: mass2, logL2, logT2, logmdot2


      if (firstTime) then
         call readinTracks("mist", thisTable)
         firstTime = .false.
         call writeInfo("MIST tracks successfully read", FORINFO)
      endif


      ! FIXME - mist models start at pre-MS - need to start as ZAMS

      ! find relevant Schaller track file given source's initial mass 
      call locate(thisTable%initialMass, thisTable%nMass, source%initialmass, i)
      
      ! interpolate between rows j,j+1 in lower mass boundary file (i)
      call locate(thisTable%age(i,1:thisTable%nAges(i)), thisTable%nAges(i), source%age, j)
      t = (source%age - thisTable%age(i, j))/(thisTable%age(i,j+1)-thisTable%age(i,j))
      mass1 = thisTable%massAtAge(i,j) + t * (thisTable%massAtAge(i,j+1)-thisTable%massAtAge(i,j))
      logL1 = thisTable%logL(i,j) + t * (thisTable%logL(i,j+1)-thisTable%logL(i,j))
      logT1 = thisTable%logTeff(i,j) + t * (thisTable%logTeff(i,j+1)-thisTable%logTeff(i,j))
      logMdot1 = thisTable%mdot(i,j) + t * (thisTable%mdot(i,j+1)-thisTable%mdot(i,j))

      ! interpolate between rows j,j+1 in upper mass boundary file (i+1)
      call locate(thisTable%age(i+1,1:thisTable%nAges(i+1)), thisTable%nAges(i+1), source%age, j)
      if (source%age < thisTable%age(i+1, thisTable%nAges(i+1))) then
         ! if source age < final age of reference star, interpolate as normal
         t = (source%age - thisTable%age(i+1, j))/(thisTable%age(i+1,j+1)-thisTable%age(i+1,j))
         mass2 = thisTable%massAtAge(i+1,j) + t * (thisTable%massAtAge(i+1,j+1)-thisTable%massAtAge(i+1,j))
         logL2 = thisTable%logL(i+1,j) + t * (thisTable%logL(i+1,j+1)-thisTable%logL(i+1,j))
         logT2 = thisTable%logTeff(i+1,j) + t * (thisTable%logTeff(i+1,j+1)-thisTable%logTeff(i+1,j))
         logMdot2 = thisTable%mdot(i+1,j) + t * (thisTable%mdot(i+1,j+1)-thisTable%mdot(i+1,j))
      else
         ! otherwise just set to final values (EOF) 
         mass2 = thisTable%massAtAge(i+1, thisTable%nAges(i+1)) 
         logL2 = thisTable%logL(i+1, thisTable%nAges(i+1)) 
         logT2 = thisTable%logTeff(i+1, thisTable%nAges(i+1)) 
         logMdot2 = thisTable%mdot(i+1, thisTable%nAges(i+1)) 
      endif

      ! interpolate between the two files and set source properties
      u = (source%initialmass - thisTable%initialMass(i))/(thisTable%initialMass(i+1) - thisTable%initialMass(i))
      
      source%mass = (mass1 + (mass2 - mass1) * u) * mSol
      source%luminosity = logL1 + (logL2 - logL1) * u
      source%luminosity = (10.d0**source%luminosity) * lSol
      source%teff = logT1 + (logT2  - logT1) * u
      source%teff = 10.d0**source%teff
      source%radius = sqrt(source%luminosity / (fourPi * stefanBoltz * source%teff**4))/1.d10
      source%mdotWind = 0.d0
      if (.not.ANY(thisTable%mdot(i:i+1,j:j+1) == 0.d0)) then
         source%mdotWind = 10.d0**(logmdot1 + (logmdot2  - logmdot1) * u)
         source%mDotWind = source%mDotWind * msol/(365.25*24.d0*3600.d0)
      endif

!     ! update spectrum. If tlusty spectrum is not found for a source, kurucz spectrum is used instead
!     call fillSpectrumTlusty(source%spectrum, source%teff, source%mass, source%radius*1.d10)
      
      call emptySurface(source%surface)
      call buildSphereNBody(source%position, source%accretionRadius/1.d10, source%surface, 20)

    end subroutine updateSourceProperties

    subroutine setSourceArrayProperties(source, nSource, fractionOfAccretionLum)
      type(SOURCETYPE) :: source(:)
      integer :: nSource, i
      real(double) :: lumAcc, tAcc, fractionOfAccretionLum


      do i = 1, nSource
         source(i)%age = 1.e5
         source(i)%initialmass = source(i)%mass/msol
         call getHosokawaProperties(source(i))
         lumAcc = bigG *source(i)%mass * source(i)%mdot / (source(i)%radius*1.d10) * fractionOfAccretionLum
         source(i)%luminosity = source(i)%luminosity + lumAcc
         tAcc = (lumAcc / (fourPi*stefanBoltz*source(i)%radius**2*1.d20))**0.25d0
         if (Writeoutput .and. nSource < 1000) then
            write(*,*) "Information for source: ",i
            write(*,*) "Mass: ",source(i)%mass/msol
            write(*,*) "radius: ",source(i)%radius*1.d10/rsol
            write(*,*) "lum: ", source(i)%luminosity/lSol
            write(*,*) "teff: ",source(i)%teff
            write(*,*) "accretion temp: ",tacc
            write(*,*) "accretion lum: ",lumAcc/lsol
            write(*,*) "fraction ",fractionOfAccretionLum
         endif
         call fillSpectrumkurucz(source(i)%spectrum, source(i)%teff, source(i)%mass, source(i)%radius*1.d10)
         if (tAcc > 0.d0) call addToSpectrumBB(source(i)%spectrum, tAcc, 1.d0)
         call normalizedSpectrum(source(i)%spectrum)
      enddo
      if (Writeoutput) call writeSourceArray("tempsource.dat")

         if (writeoutput) then
            open(77,file="spec.dat",status="unknown",form="formatted")
            do i = 1, source(1)%spectrum%nLambda
               write(77,*) source(1)%spectrum%lambda(i),source(1)%spectrum%flux(i)
            enddo
            close(77)
         endif
    end subroutine setSourceArrayProperties

    subroutine setSourceLumTemp(source, thisTable)

      type(SOURCETYPE) :: source
      type(TRACKTABLE) :: thisTable
      integer :: i, j
      real(double) :: t, u, mass1, logL1, logT1
      real(double) :: mass2, logL2, logT2

      call locate(thisTable%initialMass, thisTable%nMass, source%initialmass, i)
      
      call locate(thisTable%age(i,1:thisTable%nAges(i)), thisTable%nAges(i), source%age, j)
      t = (source%age - thisTable%age(i, j))/(thisTable%age(i,j+1)-thisTable%age(i,j))
      mass1 = thisTable%massAtAge(i,j) + t * (thisTable%massAtAge(i,j+1)-thisTable%massAtAge(i,j))
      logL1 = thisTable%logL(i,j) + t * (thisTable%logL(i,j+1)-thisTable%logL(i,j))
      logT1 = thisTable%logTeff(i,j) + t * (thisTable%logTeff(i,j+1)-thisTable%logTeff(i,j))

      call locate(thisTable%age(i+1,1:thisTable%nAges(i+1)), thisTable%nAges(i+1), source%age, j)
      t = (source%age - thisTable%age(i+1, j))/(thisTable%age(i+1,j+1)-thisTable%age(i+1,j))
      mass2 = thisTable%massAtAge(i+1,j) + t * (thisTable%massAtAge(i+1,j+1)-thisTable%massAtAge(i+1,j))
      logL2 = thisTable%logL(i+1,j) + t * (thisTable%logL(i+1,j+1)-thisTable%logL(i+1,j))
      logT2 = thisTable%logTeff(i+1,j) + t * (thisTable%logTeff(i+1,j+1)-thisTable%logTeff(i+1,j))

      u = (source%initialmass - thisTable%initialMass(i))/(thisTable%initialMass(i+1) - thisTable%initialMass(i))
      
      source%mass = (mass1 + ((mass2 - mass1) * u))*msol
      source%luminosity = logL1 + (logL2 - logL1) * u
      source%luminosity = (10.d0**source%luminosity) * lSol
      source%teff = logT1 + (logT2  - logT1) * u
      source%teff = 10.d0**source%teff
      source%radius = sqrt(source%luminosity / (fourPi * stefanBoltz * source%teff**4))/1.d10
      

    end subroutine setSourceLumTemp

      


    subroutine readinTracks(tracks, thisTable)
      character(len=*) :: tracks
      type(TRACKTABLE) :: thisTable
      character(len=200) :: tfile, thisFile, dataDirectory
      integer :: i

      select case(tracks)
         case("schaller")
            thisTable%nMass = 21
            allocate(thisTable%nAges(21))
            allocate(thisTable%initialMass(21))
            allocate(thisTable%age(21, 52))
            allocate(thisTable%massAtAge(21, 52))
            allocate(thisTable%logL(21, 52))
            allocate(thisTable%logteff(21, 52))
            allocate(thisTable%mDot(21, 52))
            thisTable%label = "Schaller (1992) evolutionary tracks"
            call readSchallerModel(thisTable, 21, 120.d0, "table1")
            call readSchallerModel(thisTable, 20, 85.d0, "table2")
            call readSchallerModel(thisTable, 19, 60.d0, "table3")
            call readSchallerModel(thisTable, 18, 40.d0, "table4")
            call readSchallerModel(thisTable, 17, 25.d0, "table5")
            call readSchallerModel(thisTable, 16, 20.d0, "table6")
            call readSchallerModel(thisTable, 15, 15.d0, "table7")
            call readSchallerModel(thisTable, 14, 12.d0, "table8")
            call readSchallerModel(thisTable, 13,  9.d0, "table9")
            call readSchallerModel(thisTable, 12, 7.d0, "table10")
            call readSchallerModel(thisTable, 11, 5.d0, "table11")
            call readSchallerModel(thisTable, 10, 4.d0, "table12")
            call readSchallerModel(thisTable, 9, 3.d0, "table13")
            call readSchallerModel(thisTable, 8, 2.5d0, "table14")
            call readSchallerModel(thisTable, 7, 2.d0, "table15")
            call readSchallerModel(thisTable, 6, 1.7d0, "table16")
            call readSchallerModel(thisTable, 5, 1.5d0, "table17")
            call readSchallerModel(thisTable, 4, 1.25d0, "table18")
            call readSchallerModel(thisTable, 3, 1.d0, "table20")
            call readSchallerModel(thisTable, 2, 0.9d0, "table21")
            call readSchallerModel(thisTable, 1, 0.8d0, "table22")

         case("mist")
            call addBibcode("2016ApJ...823..102C", "MESA MIST evolutionary tracks")
            thisTable%nMass = 196
            allocate(thisTable%nAges(thisTable%nMass))
            allocate(thisTable%initialMass(thisTable%nMass))
            allocate(thisTable%age(thisTable%nMass, 1710))
            allocate(thisTable%massAtAge(thisTable%nMass, 1710))
            allocate(thisTable%logL(thisTable%nMass, 1710))
            allocate(thisTable%logteff(thisTable%nMass, 1710))
            allocate(thisTable%logradius(thisTable%nMass, 1710))
            allocate(thisTable%mDot(thisTable%nMass, 1710))
            thisTable%label = "MIST"
            ! read filename from list, read the track from each file 
            call unixGetenv("TORUS_DATA", dataDirectory)
            write(tfile, '(a,a)') trim(dataDirectory), "/mist/filelist.txt"
            open(30, file=tfile, status="old", form="formatted")
            i = 0
10 continue
            read(30,*,end=55) thisFile 
            i = i + 1
            call readMistModel(thisTable, i, trim(thisFile))
            goto 10
55 continue
            close(30)
      end select

     end subroutine readinTracks

     subroutine getHosokawaProperties(source)
       type(SOURCETYPE) :: source
       integer, parameter :: nSteps = 176
       real(double),save :: mStar(nSteps), rStar(nSteps), rPhot(nSteps), lstar(nSteps), ltot(nSteps), tstep(nSteps)
       real(double) :: t 
       integer :: i, j
       character(len=80) :: message
       character(len=200) :: dataDirectory, tFile
       logical, save :: firstTime = .true.

       if (firstTime) then
          call unixGetenv("TORUS_DATA", dataDirectory, i)
          tfile = trim(dataDirectory)//"/md3.dat"
          open(43, file=tfile, status="old", form="formatted")
          read(43,'(a)') message
          do i = 1, nSteps
             read(43,*) j, mStar(i), rStar(i), rPhot(i), lStar(i), lTot(i), tStep(i)
          enddo
          close(43)
          mstar = mStar * mSol
          rstar = rStar * rSol
          rPhot = rPhot * rSol
          lStar = lStar * lSol
          lTot = lTot * lSol
          firstTime = .false.
       endif
       if (source%mass < mstar(1)) then
          i = 1
          t = 0.d0
       else
          call locate(mStar, nSteps, source%mass, i)
          t = (source%mass - mStar(i))/(mStar(i+1)-mstar(i))
       endif
       source%radius = (rStar(i) + t * (rStar(i+1) - rStar(i)))/1.d10
       source%luminosity = max(1.d0, lStar(i) + t * (lStar(i+1) - lStar(i)))
       source%teff = (source%luminosity / (fourpi * source%radius**2 * 1.d20 * stefanBoltz))**0.25d0

     end subroutine getHosokawaProperties
       


     subroutine readSchallerModel(thisTable, nMass, initMass, thisfile)
       type(TRACKTABLE) :: thisTable
       integer :: nMass
       real(double) :: initMass
       character(len=*) :: thisfile
       character(len=200) :: tFile, datadirectory
       integer :: nt, i
       character(len=254) :: cLine

       thisTable%initialMass(nMass) = initMass
       

       call unixGetenv("TORUS_DATA", dataDirectory, i)

       tfile = trim(dataDirectory)//"/Schaller/"//trim(thisFile)
       nt = 1
       open(31, file=tfile, status="old", form="formatted")
10 continue
       read(31,'(a)',end=55) cline
       read(cline, *) i, thisTable%age(nMass, nt),  thisTable%massAtAge(nMass, nt), &
            thisTable%logL(nMass, nt), thisTable%logteff(nMass, nt)
       read(cline,'(145X,F7.3)') thisTable%mdot(nMass,nt)
       nt = nt + 1
       goto 10
55     continue
       nt = nt - 1
       close(31)
       thisTable%nAges(nMass) = nt
     end subroutine readSchallerModel

     subroutine readMistModel(thisTable, nMass, thisfile)
       type(TRACKTABLE) :: thisTable
       integer :: nMass
       character(len=*) :: thisfile
       character(len=200) :: tFile, datadirectory
       integer :: nt, i
       character(len=254) :: cLine, junk

       call unixGetenv("TORUS_DATA", dataDirectory, i)
       ! nb: rotation is available with vvcrit0.4
       write(tfile, *) trim(dataDirectory), "/mist/MIST_v1.2_feh_p0.00_afe_p0.0_vvcrit0.0_EEPS/", trim(thisFile)

       nt = 1
       open(31, file=tfile, status="old", form="formatted")
       read(31,'(a)',end=55) cline
       read(cline, *) junk, thisTable%initialMass(nMass) 
10 continue
       read(31,'(a)',end=55) cline
       read(cline, *) thisTable%age(nMass, nt), thisTable%massAtAge(nMass, nt), thisTable%mdot(nMass,nt), &
         thisTable%logL(nMass, nt), thisTable%logteff(nMass, nt), thisTable%logRadius(nMass, nt)
       ! table gives mdot as a negative rate in Msol/yr - we want log(mdot)
       thisTable%mdot(nMass, nt) = log10(-thisTable%mdot(nMass, nt))
       
       nt = nt + 1
       goto 10
55     continue
       nt = nt - 1
       close(31)
       thisTable%nAges(nMass) = nt
     end subroutine readMistModel


     subroutine removeSource(source, nSource, n)
       use inputs_mod, only : smallestCellSize
       type(SOURCETYPE) :: source(:)
       integer :: nSource, n, i

       do i = 1, nSource
          call freeSpectrum(source(i)%spectrum)
          call emptySurface(source(i)%surface)
       enddo

       if ( n /= nSource) then
          do i = n, nSource - 1
             source(i) = source(i+1)
          enddo
       endif

       nSource = nSource - 1
       ! todo clustersinks
       do i = 1, nSource
          call buildSphereNBody(source(i)%position, 2.5d0*smallestCellSize, source(i)%surface, 20)
!          call fillSpectrumkurucz(source(i)%spectrum, source(i)%teff, source(i)%mass, source(i)%radius*1.d10)
          call fillSpectrumTlusty(source(i)%spectrum, source(i)%teff, source(i)%mass, source(i)%radius*1.d10)
       enddo
     end subroutine removeSource


!     subroutine fillSpectrum(source, nKurucz, kLabel, kSpectrum)
!       type(SOURCETYPE) :: source
!       real(double) :: logg
!       integer, parameter :: nFiles = 60
!       real,save :: teff(nFiles)
!       integer :: i, j
!       real(double) :: t
!       real :: loggArray(11)
!       logical, save :: firstTime = .true.
!       character(len=200) :: thisFile = " ", dataDirectory = " "
!       logical :: ok 
!       integer :: nKurucz
!       character(len=*) :: kLabel(:)
!       logical,save :: firstWarning = .true.
!       type(SPECTRUMTYPE) :: kSpectrum(:)
!
!       ok = .true.
!       call unixGetenv("TORUS_DATA", dataDirectory, i)
!
!
!       loggArray = (/ 000., 050., 100., 150., 200., 250., 300., 350., 400., 450., 500. /)
!       if (firsttime) then
!          open(31, file=trim(dataDirectory)//"/Kurucz/filelist.dat", form="formatted", status="old")
!          do i = 1, nFiles
!             read(31, *) teff(i)
!          end do
!          close(31)
!          firstTime = .false.
!       endif
!
!
!
!
!       logg = (bigG*source%mass*mSol)/((source%radius*1.d10)**2)
!       logg = log10(logg)
!       call locate(teff, nFiles, real(source%teff), i)
!       call locate(loggArray, 11, real(logg*100.), j)
!
!       t = (source%teff - teff(i))/(teff(i+1)-teff(i))
!       if (t > 0.5) i = i + 1
!       t = ((logg*100.) - loggArray(j))/(loggArray(j+1) - loggArray(j))
!       if (t > 0.5) j = j + 1
!       call createKuruczFilename(teff(i), loggArray(j), thisFile)
!!       call readSpectrum(source%spectrum, thisFile, ok)
!
!       call readKuruczSpectrum(source%spectrum, thisFile, klabel, kspectrum, nKurucz, ok)
!       if (ok) then
!          call normalizedSpectrum(source%spectrum)
!       else
!
!          ! try a higher gravity
!
!          if (j < 11) then
!             j = j + 1
!             call createKuruczFilename(teff(i), loggArray(j), thisFile)
!             call readKuruczSpectrum(source%spectrum, thisFile, klabel, kspectrum, nKurucz, ok)
!             if (ok) then
!                call normalizedSpectrum(source%spectrum)
!             endif
!          endif
!          if (.not. ok) then
!             do j = 1, 11
!                call createKuruczFilename(teff(i), loggArray(j), thisFile)
!                call readKuruczSpectrum(source%spectrum, thisFile, klabel, kspectrum, nKurucz, ok)
!                if (ok) then
!                   call normalizedSpectrum(source%spectrum)
!                   exit
!                endif
!             enddo
!          endif
!          if (.not.ok) then
!             if (firstWarning) then
!                call writeInfo("Cannot find appropriate model atmosphere for source: "//trim(thisFile), IMPORTANT)
!                firstWarning  = .false.
!             endif
!             call fillSpectrumBB(source%spectrum,source%teff, 100.d0, 1.d7, 1000)
!             call normalizedSpectrum(source%spectrum)
!          endif
!       endif
!
!     end subroutine fillSpectrum
!
!
!     subroutine createKuruczFileName(teff, logg, thisfile)
!       real :: teff, logg
!       integer :: i
!       character(len=*) thisfile
!       character(len=80) :: fluxfile, dataDirectory
!
!       call unixGetenv("TORUS_DATA", dataDirectory, i)
!
!
!       if (teff < 10000.) then
!          write(fluxfile,'(a,i4,a,i3.3,a)') "f",int(teff),"_",int(logg),".dat"
!       else
!          write(fluxfile,'(a,i5,a,i3.3,a)') "f",int(teff),"_",int(logg),".dat"          
!       endif
!       thisFile = trim(dataDirectory)//"/Kurucz/"//trim(fluxfile)
!     end subroutine createKuruczFileName
!
!     subroutine readKuruczGrid(label, spectrum, nFiles)
!       character(len=*) :: label(:)
!       type(SPECTRUMTYPE) :: spectrum(:)
!       integer :: nFiles
!       character(len=200) :: tfile,fluxfile,dataDirectory = " "
!       logical :: ok
!       integer :: i
!       ok = .true.
!
!       call unixGetenv("TORUS_DATA", dataDirectory, i)
!       
!       call writeInfo("Reading Kurucz grid...",TRIVIAL)
!
!       tfile = trim(dataDirectory)//"/Kurucz/files.dat"
!       open(31, file = tfile, status = "old", form="formatted")
!       do i = 1, nFiles
!          read(31,*) fluxfile
!          label(i) = trim(fluxfile)
!          tfile = trim(dataDirectory)//"/Kurucz/"//trim(fluxfile)
!          call readSpectrum(spectrum(i), tfile, ok)
!       enddo
!       close(31)
!       call writeInfo("Done.",TRIVIAL)
!
!     end subroutine readKuruczGrid
!
!     subroutine readKuruczSpectrum(thisSpectrum,thisLabel, label, spectrum, nFiles, ok)
!       character(len=*) :: label(:), thisLabel
!       type(SPECTRUMTYPE) :: spectrum(:), thisSpectrum
!       integer :: nFiles
!       logical :: ok
!       integer :: i
!       
!       ok = .false.
!       do i = 1, nFiles
!          if (trim(label(i)).eq.trim(thisLabel)) then
!             call copySpectrum(thisSpectrum, spectrum(i))
!             write(*,*) "spectrum copied"
!             ok = .true.
!             exit
!          endif
!       enddo
!     end subroutine readKuruczSpectrum

   subroutine testTracks
      real(double) :: t, dt, totalCreatedMass
      integer :: i, j
      character(len=100) :: fn
      type(SOURCETYPE), pointer :: src=>null()

      ! create sources
      if (associated(globalSourceArray)) then
         deallocate(globalSourceArray)
         globalSourceArray => null()
      endif
      globalnsource = 2
      allocate(globalsourceArray(1:globalnsource))
      globalSourcearray(1:globalnsource)%mass = 500.d0 * msol
      globalSourcearray(1:globalnsource)%age = 0.d0

      do i = 1, globalnsource
         call createSources(globalsourceArray(i)%nSubsource, globalsourceArray(i)%subsourceArray, "instantaneous", & 
            0.d0, 120.d0, 0.d0, totalCreatedMass, zeroNsource=.false.)
         if (writeoutput) write(*,*) i, " ... created ", totalCreatedMass, " Msol" 
      enddo

      ! write headers
      if (writeoutput) then
         do i = 1, globalnSource
            do j = 1, globalsourceArray(i)%nSubsource
               write(fn, '(a,i4.4,a,i4.4,a)') "track_", i, "_", j, ".dat"
               open(400, file=trim(fn), status="replace", form="formatted")
               write(400, '(6(a12,1x))') "# age", "mass", "mdot", "logL", "logTeff", "logR"
               close(400)
             enddo
         enddo
      endif

      ! update ages
      t = 0.d0
      dt = 1.d4
      do while (t <= 3.d6)
         ! update from track
         do i = 1, globalnSource
            do j = 1, globalsourceArray(i)%nSubsource
               call updateSourceProperties(globalsourcearray(i)%subsourceArray(j))
            enddo
         enddo
         

         ! write out 
         if (writeoutput) then
            do i = 1, globalnSource
               do j = 1, globalsourceArray(i)%nSubsource
                  write(fn, '(a,i4.4,a,i4.4,a)') "track_", i, "_", j, ".dat"
                  open(400, file=trim(fn), status="old", position="append", form="formatted")
                  src => globalsourcearray(i)%subsourcearray(j)
                  write(400, '(6(es12.5,1x))') src%age, src%mass/msol, src%mdotwind/msol/secstoyears,&
                      log10(src%luminosity/lsol), log10(src%teff), log10(src%radius*1.d10/rsol)
                  close(400)
               enddo
            enddo
         endif

         ! evolve
         globalSourceArray(1:globalnSource)%age = globalSourceArray(1:globalnSource)%age + dt
         do i = 1, globalnSource
            if (globalSourceArray(i)%nSubsource > 0) then 
               globalSourceArray(i)%subsourceArray(1:globalSourceArray(i)%nSubsource)%age = & 
               globalSourceArray(i)%subsourceArray(1:globalSourceArray(i)%nSubsource)%age + dt 
            endif
         enddo
         t = t + dt
      enddo
      stop
   end subroutine testTracks


end module starburst_mod
