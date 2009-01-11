module starburst_mod

! A module to create sources, and their distribution, for
! starburst dust/photoionization calculations

! started by tjh on 23/8/2006



  use kind_mod
  use constants_mod
  use messages_mod
  use utils_mod
  use source_mod
  use unix_mod, only: unixGetenv

  implicit none

  type TRACKTABLE
     character(len=80) :: label
     integer :: nMass
     real(double), pointer :: initialMass(:)
     integer, pointer :: nAges(:)
     real(double), pointer :: age(:,:)
     real(double), pointer :: massAtAge(:,:)
     real(double), pointer :: logL(:,:)
     real(double), pointer :: logTeff(:,:)
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

      case DEFAULT
         write(message,'(a,a)') "IMF not recognised: ", trim(imfType)
         call writefatal(message)

    end select


    do i = 2, nMass
       prob(i) = prob(i) + prob(i-1)
    enddo
    prob(1:nMass) = prob(1:nMass) - prob(1)
    prob(1:nMass) = prob(1:nMass) / prob(nMass)
    call random_number(r)
    call locate(prob, nMass, r, i)
    t = (r - prob(i))/(prob(i+1)-prob(i))
    mass = log10(massArray(i)) + t * (log10(massArray(i+1))-log10(massArray(i)))
    mass = 10.d0**mass
  end function randomMassFromIMF

  subroutine createSources(nSource, source, burstType, burstAge, burstMass, sfRate)
    integer :: nSource
    type(SOURCETYPE) :: source(:)
    type(TRACKTABLE) :: thisTable
    character(len=80) :: message
    character(len=*) :: burstType
    real(double) :: burstAge
    real(double) :: burstMass, totMass
    real(double) :: sfRate
    integer :: i
    integer :: nDead, nSupernova, nOB
    integer, parameter :: nKurucz = 410
    type(SPECTRUMTYPE) :: kSpectrum(nKurucz)
    character(len=80) :: klabel(nKurucz)

    kLabel = " "
    thisTable%label = " "
    call  readKuruczGrid(klabel, kspectrum, nKurucz)


    nSource = 0

    ! set up the initial number of stars and their masses and ages

    select case(burstType)

       case("continuous")
          burstMass = sfRate * burstAge
          totMass = 0.d0
          do while (totMass < burstMass)
             nSource = nSource + 1
             call random_number(source(nSource)%age)
             source(nSource)%age = source(nSource)%age * burstAge
             source(nSource)%initialmass = randomMassFromIMF("salpeter", 0.8d0, 120.d0, -2.35d0)
             totMass = totMass + source(nSource)%initialmass
          enddo

       case("instantaneous")
          totMass = 0.d0
          do while (totMass < burstMass)
             nSource = nSource + 1
             source(nSource)%initialMass = randomMassFromIMF("salpeter", 0.8d0, 120.d0, -2.35d0)
             totMass = totMass + source(nSource)%initialMass
          enddo
          source(1:nSource)%age = burstAge

       case DEFAULT
         write(message,'(a,a)') "Burst type not recognised: ", trim(burstType)
         call writefatal(message)
             
      end select

      write(*,*) "number of sources in burst", nSource
      write(*,*) "burst mass",totMass

      ! now get actual masses, temps, and luminosities, and radii for age from evolution tracks

      call readinTracks("schaller", thisTable)
      call writeInfo("Schaller tracks successfully read", FORINFO)
      i = 1
      nDead = 0
      nSupernova = 0
      do while (i <= nSource)
         if (.not.isSourceDead(source(i), thisTable)) then
            call setSourceProperties(source(i), thisTable)
            i = i + 1
         else
            nDead = nDead + 1
            if (source(i)%initialMass > 8.d0) then
               nSupernova = nSupernova + 1
            endif
            call removeSource(source, nSource, i)
         endif
      enddo
      write(*,*) "Number of dead sources",  nDead
      write(*,*) "Number of supernova",  nSupernova
      write(*,*) "Burst luminosity",SUM(source(1:nSource)%luminosity)/lsol

      nOB = 0
      do i = 1, nSource
         if (source(i)%mass > 15.) then
            nOB = nOB + 1
         endif
      enddo
      write(*,*) "Number of OB stars (>15 msol): ",nOB

      do i = 1, nSource
         call fillSpectrum(source(i), nKurucz, kLabel, kSpectrum)
      enddo
      

    end subroutine createSources


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

    subroutine setSourceProperties(source, thisTable)
      type(SOURCETYPE) :: source
      type(TRACKTABLE) :: thisTable
      integer :: i, j
      real(double) :: r, t, u, mass1, logL1, logT1
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
      
      source%mass = mass1 + (mass2 - mass1) * u
      source%luminosity = logL1 + (logL2 - logL1) * u
      source%luminosity = (10.d0**source%luminosity) * lSol
      source%teff = logT1 + (logT2  - logT1) * u
      source%teff = 10.d0**source%teff
      source%radius = sqrt(source%luminosity / (fourPi * stefanBoltz * source%teff**4))/1.d10
      
      source%position = randomUnitVector()
      call random_number(r)
      r = sqrt(r)
      source%position = source%position * (5.d0 * pctocm / 1.d10) * r

    end subroutine setSourceProperties

      

    subroutine readinTracks(tracks, thisTable)
      character(len=*) :: tracks
      type(TRACKTABLE) :: thisTable

      select case(tracks)
         case("schaller")
            thisTable%nMass = 21
            allocate(thisTable%nAges(21))
            allocate(thisTable%initialMass(21))
            allocate(thisTable%age(21, 52))
            allocate(thisTable%massAtAge(21, 52))
            allocate(thisTable%logL(21, 52))
            allocate(thisTable%logteff(21, 52))
            thisTable%label = "Scharer (1992) evolutionary tracks"
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
       end select

     end subroutine readinTracks



     subroutine readSchallerModel(thisTable, nMass, initMass, thisfile)
       type(TRACKTABLE) :: thisTable
       integer :: nMass
       real(double) :: initMass
       character(len=*) :: thisfile
       character(len=200) :: tFile, datadirectory = " "
       integer :: nt, i
       thisTable%initialMass(nMass) = initMass
       

       call unixGetenv("TORUS_DATA", dataDirectory, i)

       tfile = trim(dataDirectory)//"/Schaller/"//trim(thisFile)
       nt = 1
       open(31, file=tfile, status="old", form="formatted")
10 continue
       read(31, *, end=55) i, thisTable%age(nMass, nt),  thisTable%massAtAge(nMass, nt), &
            thisTable%logL(nMass, nt), thisTable%logteff(nMass, nt)
       nt = nt + 1
       goto 10
55     continue
       nt = nt - 1
       close(31)
       thisTable%nAges(nMass) = nt
     end subroutine readSchallerModel


     subroutine removeSource(source, nSource, n)
       type(SOURCETYPE) :: source(:)
       integer :: nSource, n, i

       if ( n /= nSource) then
          do i = n, nSource - 1
             source(i) = source(i+1)
          enddo
       endif
       nSource = nSource - 1
     end subroutine removeSource


     subroutine fillSpectrum(source, nKurucz, kLabel, kSpectrum)
       type(SOURCETYPE) :: source
       real(double) :: logg
       integer, parameter :: nFiles = 60
       real,save :: teff(nFiles)
       integer :: i, j
       real(double) :: t
       real :: loggArray(11)
       logical, save :: firstTime = .true.
       character(len=200) :: thisFile = " ", dataDirectory = " "
       logical :: ok 
       integer :: nKurucz
       character(len=*) :: kLabel(:)
       type(SPECTRUMTYPE) :: kSpectrum(:)

       ok = .true.
       call unixGetenv("TORUS_DATA", dataDirectory, i)


       loggArray = (/ 000., 050., 100., 150., 200., 250., 300., 350., 400., 450., 500. /)
       if (firsttime) then
          open(31, file=trim(dataDirectory)//"/Kurucz/filelist.dat", form="formatted", status="old")
          do i = 1, nFiles
             read(31, *) teff(i)
          end do
          close(31)
          firstTime = .false.
       endif




       logg = (bigG*source%mass*mSol)/((source%radius*1.d10)**2)
       logg = log10(logg)
       call locate(teff, nFiles, real(source%teff), i)
       call locate(loggArray, 11, real(logg*100.), j)

       t = (source%teff - teff(i))/(teff(i+1)-teff(i))
       if (t > 0.5) i = i + 1
       t = ((logg*100.) - loggArray(j))/(loggArray(j+1) - loggArray(j))
       if (t > 0.5) j = j + 1
       call createKuruczFilename(teff(i), loggArray(j), thisFile)
       call readSpectrum(source%spectrum, thisFile, ok)

       call readKuruczSpectrum(source%spectrum, thisFile, klabel, kspectrum, nKurucz, ok)
       if (ok) then
          call normalizedSpectrum(source%spectrum)
       else

          ! try a higher gravity

          if (j < 11) then
             j = j + 1
             call createKuruczFilename(teff(i), loggArray(j), thisFile)
             call readKuruczSpectrum(source%spectrum, thisFile, klabel, kspectrum, nKurucz, ok)
             if (ok) then
                call normalizedSpectrum(source%spectrum)
             endif
          endif
          if (.not. ok) then
             do j = 1, 11
                call createKuruczFilename(teff(i), loggArray(j), thisFile)
                call readKuruczSpectrum(source%spectrum, thisFile, klabel, kspectrum, nKurucz, ok)
                if (ok) then
                   call normalizedSpectrum(source%spectrum)
                   exit
                endif
             enddo
          endif
          if (.not.ok) then
             call writeFatal("Cannot find appropriate model atmosphere for source")
          endif
       endif

     end subroutine fillSpectrum


     subroutine createKuruczFileName(teff, logg, thisfile)
       real :: teff, logg
       character(len=*) thisfile
       character(len=80) :: fluxfile
       if (teff < 10000.) then
          write(fluxfile,'(a,i4,a,i3.3,a)') "f",int(teff),"_",int(logg),".dat"
       else
          write(fluxfile,'(a,i5,a,i3.3,a)') "f",int(teff),"_",int(logg),".dat"          
       endif
       thisFile = trim(fluxfile)
     end subroutine createKuruczFileName

     subroutine readKuruczGrid(label, spectrum, nFiles)
       character(len=*) :: label(:)
       type(SPECTRUMTYPE) :: spectrum(:)
       integer :: nFiles
       character(len=200) :: tfile,fluxfile,dataDirectory = " "
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
             thisSpectrum = spectrum(i)
             ok = .true.
             exit
          endif
       enddo
     end subroutine readKuruczSpectrum


end module starburst_mod
