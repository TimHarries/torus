#ifdef SPH
module sph_data_class

  use vector_mod
  use constants_mod

  ! 
  ! Class definition for Mathew's SPH data.
  ! 
  ! LOGS: 
  ! Created on Jan-10-2003 (R. Kurosawa)
  ! 

  implicit none

  public:: &
       kill, &
       info_sph, &
       get_udist, &
       get_umass, &
       get_utime, &
       get_npart, &
       get_nptmass, &
       get_time, &
       get_vel, &
       get_position_gas_particle, &
       Put_position_gas_particle, &
       get_rhon, &
       put_rhon, &
       get_position_pt_mass, &
       get_pt_position, &
       get_pt_mass, &
       get_pt_velocity, &
       get_rhon_min, &
       get_rhon_max, &
       get_spins, &
       get_stellar_disc_parameters, &
       stellar_disc_exists, &
       find_inclinations, &
       ClusterParameter, &
       isAlive, &
       sphData, &
       sphVelocityPresent, &
       read_sph_data_wrapper, & ! only the wrapper routine is public
       sph_mass_within_grid, &
       npart, init_sph_data, & ! for distributeSphDataOverMPI
       deallocate_sph

#ifdef WITHSPHNG
  public :: init_sphtorus
#endif


! Default is private  
  private

  ! At a given time (time)
  type sph_data
     logical      :: inUse=.false.          ! Flag to indicate if this object is in use.
     real(double) :: udist, umass, utime, uvel, utemp    ! Units of distance, mass, time in cgs
     real(double) :: codeVelocitytoTORUS    ! Conversion from SPH code velocity units to Torus units
     real(double) :: codeEnergytoTemperature    ! Conversion from SPH code velocity units to Torus units
     !                                          ! (umass is M_sol, udist=0.1 pc)
     integer      :: npart                  ! Total number of gas particles (field+disc)
     real(double) :: time                   ! Time of sph data dump (in units of utime)
     integer          :: nptmass                ! Number of stars/brown dwarfs
     real(double), pointer, dimension(:) :: gasmass            ! Mass of each gas particle ! DAR changed to allow variable mass
     real(double)                        :: totalgasmass       ! Total gas mass summed over all SPH particles.
     real(double)                        :: totalHImass        ! Total HI mass summed over all SPH particles. 
     real(double)                        :: totalMolmass       ! Total molecular gas mass summed over all SPH particles.
     ! Positions of gas particles
     real(double), pointer, dimension(:) :: xn,yn,zn
     real(double), pointer, dimension(:) :: vxn,vyn,vzn
     real(double), pointer, dimension(:) :: hn                 ! Smoothing length
     ! Density of the gas particles
     real(double), pointer, dimension(:) :: rhon
     ! Dust fraction
     real(double), pointer, dimension(:) :: dustfrac
     real(double), pointer, dimension(:,:) :: dustfracs
     ! Density of H2
     real(double), pointer, dimension(:) :: rhoH2 => null()
     ! Density of CO
     real(double), pointer, dimension(:) :: rhoCO => null()
     ! Temperature of the gas particles
     real(double), pointer, dimension(:) :: temperature
     ! Smoothing lengths of the sink particles
     real(double), pointer, dimension(:) :: hpt
     ! Positions of stars
     real(double), pointer, dimension(:) :: x,y,z
     real(double), pointer, dimension(:) :: vx,vy,vz
     !
     real(double), pointer, dimension(:) :: ptmass ! Masses of stars
     ! 
     !
     ! Some extra data for the stellar disk.
     ! Note: the units used here are different from the ones used above!     
     logical :: have_stellar_disc            ! T if the following data are assigned
     real(double), pointer, dimension(:) :: discrad ! in [10^10cm]
     real(double), pointer, dimension(:) :: discmass   ! in [g]
     real(double), pointer, dimension(:) :: spinx, spiny, spinz
          
  end type sph_data

! Mean molecular weight, used for calculating temperature from internal energy.
! This assumes a 10:1 H:He ratio by number
    real(double), parameter, private :: gmw = 14.0/11.0 

! The following variables will be set by molecular_mod
    character(len=10), public, save :: sph_molecule = "NOT SET   "
    real(double), public, save      :: sph_molecularWeight = -1.0e30
    
  real(double), allocatable :: PositionArray(:,:), OneOverHsquared(:), &
                               RhoArray(:), TemArray(:), VelocityArray(:,:), &
                               RhoH2Array(:), rhoCOarray(:), dustfrac(:)

  type(sph_data), save :: sphdata
  integer, save :: npart

! These are only read by clusterparameter so don't need to be threadprivate
  real(double), save :: rcrit, rmax
  real(double), allocatable :: xarray(:), etaarray(:)
  real(double) :: maxx, maxy, maxz, maxr2 ! maximum extents of particles

! This is the maximum line length when reading an ASCII file, specified as a module parameter
! as it appears in multiple subroutines. If the line is longer than this length itype is likely 
! to be lost leading to out of bounds array access. 
  integer, parameter, private :: MaxAsciiLineLength=1024
! Maximum number of words per line when reading an ASCII file
  integer, parameter, private :: MaxWords=50
! Flag to specify when we are reading a Gadget ASCII dump 
  logical :: gadget = .false.

  interface kill
     module procedure kill_sph_data
  end interface
  
  
contains  


  ! 
  ! Initializes an object with parameters (if possible).
  ! 
  subroutine init_sph_data(udist, umass, utime,  time, nptmass, uvel, utemp)
    use inputs_mod, only: discardSinks
    implicit none

    real(double), intent(in)  :: udist, umass, utime    ! Units of distance, mass, time in cgs
    real(double), optional  :: uvel, utemp    ! Units of distance, mass, time in cgs
    !                                                       ! (umass is M_sol, udist=0.1 pc)
    real(double), intent(in)  :: time                   ! Time of sph data dump (in units of utime)
    integer, intent(in)           :: nptmass                ! Number of stars/brown dwarfs

    ! Indicate that this object is in use
    sphdata%inUse = .true.

    ! save these values in this object
    sphdata%udist = udist
    sphdata%umass = umass
    sphdata%utime = utime
    if(present(uvel)) then
       sphdata%uvel = uvel
    endif
    if(present(utemp)) then
       sphdata%utemp = utemp
    endif

    sphdata%npart = npart
    sphdata%time = time


    ! allocate arrays
    ! -- for gas particles
    ALLOCATE(sphdata%xn(sphdata%npart))
    ALLOCATE(sphdata%yn(sphdata%npart))
    ALLOCATE(sphdata%zn(sphdata%npart))


    ALLOCATE(sphdata%vxn(sphdata%npart))
    ALLOCATE(sphdata%vyn(sphdata%npart))
    ALLOCATE(sphdata%vzn(sphdata%npart))

    Allocate(sphdata%rhon(sphdata%npart))
    ALLOCATE(sphdata%temperature(sphdata%npart))
    ALLOCATE(sphdata%gasmass(sphdata%npart))
    ALLOCATE(sphdata%hn(sphdata%npart))

! If sinks are discarded then the point mass arrays will be set to 
! zero size rather than left unallocated
    if (discardSinks) then
       call writeInfo("Sink particles will be discarded", FORINFO)
       sphdata%nptmass = 0
    else
       if (nptmass > 0) call writeInfo("Sink particles will be stored", FORINFO)
       sphdata%nptmass = nptmass
    endif

    ! -- for star positions
    ALLOCATE(sphdata%x(nptmass))
    ALLOCATE(sphdata%y(nptmass))
    ALLOCATE(sphdata%z(nptmass))
       
    ALLOCATE(sphdata%vx(nptmass))
    ALLOCATE(sphdata%vy(nptmass))
    ALLOCATE(sphdata%vz(nptmass))
    ALLOCATE(sphdata%hpt(nptmass))

       
    ! -- for mass of stars
    ALLOCATE(sphdata%ptmass(nptmass))

    !
    sphdata%have_stellar_disc = .false. 

  end subroutine init_sph_data

#ifdef WITHSPHNG
  ! 
  ! Initializes an object with parameters when torus is called as a subroutine from sphNG.
  ! 
  subroutine init_sphtorus(udist, umass, utime, b_num_gas, time, nptmass, &
        b_npart, b_idim, b_iphase, b_xyzmh, b_rho, b_temp, b_totalgasmass)
    implicit none

! Arguments --------------------------------------------------------------------
    real(double), intent(in)  :: udist, umass, utime    ! Units of distance, mass, time in cgs
    !                                                   ! (umass is M_sol, udist=0.1 pc)
    integer, intent(in)           :: b_num_gas          ! Number of gas particles (field+disc)
    real(double), intent(in)  :: time                   ! Time of sph data dump (in units of utime)
    integer, intent(in)           :: nptmass            ! Number of stars/brown dwarfs

    integer, intent(in)   :: b_npart, b_idim
    integer(kind=1), intent(in) :: b_iphase(b_idim)
    real(double), intent(in)    :: b_xyzmh(5,b_idim)
    real, intent(in)    :: b_rho(b_idim)
    real(double), intent(in)    :: b_temp(b_idim)
    real(double), intent(in)    :: b_totalgasmass
! Local variables --------------------------------------------------------------
    integer :: iii, iiipart, iiigas

! Begin executable statements --------------------------------------------------

    npart = b_num_gas

    ! Indicate that this object is in use
    sphData%inUse = .true.

! no conversion required
    sphdata%codeEnergytoTemperature = 1.0 

    ! save these values in this object
    sphData%udist = udist
    sphData%umass = umass
    sphData%utime = utime
    sphData%npart = npart
    sphData%time = time
    sphData%nptmass = nptmass
    sphData%totalgasmass = b_totalgasmass


    ! allocate arrays
    ! -- for gas particles
    ALLOCATE(sphData%xn(npart))
    ALLOCATE(sphData%yn(npart))
    ALLOCATE(sphData%zn(npart))
    ALLOCATE(sphData%hn(npart))
    ALLOCATE(sphData%rhon(npart))
    ALLOCATE(sphData%temperature(npart))
    ALLOCATE(sphData%gasmass(npart))


    ! -- for star positions
    ALLOCATE(sphData%x(nptmass))
    ALLOCATE(sphData%y(nptmass))
    ALLOCATE(sphData%z(nptmass))
    
    ! -- for mass of stars
    ALLOCATE(sphData%ptmass(nptmass))

    !
    sphData%have_stellar_disc = .false. 
    
    iiipart=0
    iiigas=0

    ! Set up density and position of gas and point mass particles
    do iii=1,b_npart
       if (b_iphase(iii) == 0) then
          iiigas=iiigas+1
          if (iiigas > npart) then
             write (*,*) 'Error: iiigas>npart',iiigas, npart
             cycle
          endif
          sphData%rhon(iiigas)        = b_rho(iii)
          sphData%temperature(iiigas) = b_temp(iii)
          sphData%xn(iiigas)          = b_xyzmh(1,iii)
          sphData%yn(iiigas)          = b_xyzmh(2,iii)
          sphData%zn(iiigas)          = b_xyzmh(3,iii)
          sphData%gasmass(iiigas)     = b_xyzmh(4,iii)
          sphData%hn(iiigas)          = b_xyzmh(5,iii)

       elseif (b_iphase(iii) > 0) then
          iiipart=iiipart+1
          if (iiipart > nptmass) then
             write (*,*) 'Error: iiipart>nptmass',iiipart,nptmass
             cycle
          endif
          sphData%x(iiipart)      = b_xyzmh(1,iii)
          sphData%y(iiipart)      = b_xyzmh(2,iii)
          sphData%z(iiipart)      = b_xyzmh(3,iii)
          sphData%ptmass(iiipart) = b_xyzmh(4,iii)
       endif
    end do

    ! Check the number of gas particles and point masses are as expected
    IF ( iiigas /= npart ) THEN
       print *,'Error: expecting to find', npart, 'gas particles but found', iiigas
       do; end do
    ENDIF
    IF ( iiipart /= nptmass ) THEN
       print *,'Error: expecting to find', nptmass, 'point masses but found', iiipart
       do; end do 
    ENDIF

  end subroutine init_sphtorus
#endif

! Wrapper subroutine which calls the appropriate read subroutine 
! according to the string inputFileFormat
! D. Acreman, June 2012
  subroutine read_sph_data_wrapper
    use inputs_mod, only: inputFileFormat, sphdatafilename, limitScalar

    real(double) :: minMass, maxMass
    character(len=150) :: message

    select case (inputFileFormat)

    case("mpi","binary")
       call read_sph_data_mpi(sphdatafilename)

    case("clump")
       call read_sph_data_clumpfind(sphdatafilename)

    case("ascii")
       call new_read_sph_data(sphdatafilename)

    case("ascii-gadget")
       gadget=.true.
       call new_read_sph_data(sphdatafilename)

    case("gadget2")
       call read_gadget2_data(sphdatafilename)

       

    case default
       call writeFatal("Unrecognised file format "//inputFileFormat)

    end select

! Report particle masses and check mass splitting condition against the particle masses
! Do this in the wrapper routine as it is the same for all types of input file
    minMass = minVal(sphdata%gasmass) * sphdata%umass
    maxMass = maxVal(sphdata%gasmass) * sphdata%umass
    call writeInfo("",forInfo)
    if (minMass == maxMass) then
       write(message,'(a,es11.4,a,es11.4,a)') "Equal mass particles of ", minMass, "g = ", minMass/mSol, " solar masses" 
       call writeInfo(message,forInfo)
    else
       write(message,'(a,es11.4,a,es11.4,a)') "Minimum particle mass  ", minMass, "g = ", minMass/mSol, " solar masses"
       call writeInfo(message,forInfo)
       write(message,'(a,es11.4,a,es11.4,a)') "Maximum particle mass  ", maxMass, "g = ", maxMass/mSol, " solar masses"
       call writeInfo(message,forInfo)
    endif

    write(message,'(a,es11.4,a,es11.4,a)') "Mass splitting condition  ", limitScalar, "g = ", limitScalar/mSol, " solar masses"
    call writeInfo(message,forInfo)
    if (limitscalar >= minMass) then
       write(message,'(a,f15.2,a)') "Mass splitting condition is ", limitscalar/minMass, " particle masses"
       call writeInfo(message,forInfo)
    else
       call writeWarning("Mass splitting condition is less than minimum particle mass")
    endif
    call writeInfo("",forInfo)

! Initialise arrays to be used by clusterParameter function
    if (get_npart() > 0) call init_clusterparameter_arrays

  end subroutine read_sph_data_wrapper

! Read SPH data from a splash ASCII dump.
  subroutine new_read_sph_data(rootfilename)
    use inputs_mod, only: convertRhoToHI, sphh2col, sphmolcol, sphwithchem, iModel, discardSinks
    use inputs_mod, only: dragon, sphdensitycut, sphVelOffset, sphPosOffset
    use angularImage_utils, only:  internalView, galaxyPositionAngle, galaxyInclination
    use utils_mod, only : findMultiFilename
    use parallel_mod, only: torus_abort
    
    implicit none

    character(LEN=*), intent(in)  :: rootfilename
    character(LEN=80) :: filename
    character(len=24) :: word(MaxWords), unit(MaxWords), pType(MaxWords), unitName(maxWords)
    integer :: nword, nunit, npType, status
    !   
    integer, parameter  :: LUIN = 10 ! logical unit # of the data file
    real(double) :: udist, umass, utime,  time, uvel, utemp, ufac
    real(double) :: xn, yn, zn, vx, vy, vz, gaspartmass, rhon, u, h, h2ratio, gmw
    integer :: itype ! splash particle type, different convention to SPHNG 
    integer :: ipart, icount, iptmass, igas, istar, idead, i
    integer :: nptmass, nstar, nother, nlines
    integer :: irequired, irejected ! Record number of gas particles required and not required
    integer, parameter :: maxNumVals=50 ! Maximum number of values per line (first dim of junkArray)
    real(double), allocatable :: junkArray(:,:)
    character(LEN=1)  :: junkchar
    character(LEN=150) :: message
    character(len=MaxAsciiLineLength) :: namestring, unitString, pTypeString, unitNameString
    integer :: ix, iy, iz, ivx, ivy, ivz, irho, iu, iitype, ih, imass, iUoverT, ipType, iDustTemperature
    logical :: haveUandUoverT, haveDustTemperature
    integer :: iDustfrac
    real(double) :: dustfrac
    integer, allocatable :: pNumArray(:) ! Array of particle numbers read from header

!
! For SPH simulations with chemistry
!
! CO fraction. 
    real(double)       :: COfrac

! Account for time in Galactic plane surveys
    logical, parameter :: multiTime=.false.
    real(double) :: extraPA
    

!
! 1. Read in data from an ASCII SPH dump
!

    if (gadget) then
       call writeInfo("Reading ASCII dump from Gadget", FORINFO)
    else if (dragon) then 
       call writeInfo("Reading ASCII dump from Dragon", FORINFO)
    else 
       call writeInfo("Reading ASCII dump", FORINFO)
    end if

! 1.1 Open the file

    call findMultiFilename(rootfilename, iModel, filename)
    open(unit=LUIN, file=TRIM(filename), form="formatted",status="old")

! 1.2 Read the header

    read(LUIN,*) 
    read(LUIN,*)
    read(LUIN,*)
    read(LUIN,*) junkchar, time, utime
    read(LUIN,*)

!
! Read in particle types and number of each type
!

! Read in particle types
    read(LUIN,'(a)') pTypeString
! Discard the first 8 characters ("# npart:"). 
    pTypeString = pTypeString(9:)
    write(*,*) "*",trim(pTypeString)
    call splitintoWords(pTypeString, pType, npType, wordLen=13, adjL=.false.)
    do i = 1, npType
       ptype(i) = adjustl(ptype(i))
    enddo
    
! Read in number of particles
    allocate(pNumArray(npType))
    read(LUIN,*) junkchar, pNumArray(:)

    write(*,*) "Particle types: "
    do i = 1, npType
       write(*,*) trim(pType(i)),": ",pNumArray(i)
    enddo

! Get the number of each type of particle. The names depend on which code the dump came from (see the read_data_*.f90 files
! in the splash source code). However there are only three cases we care about here:
! 1. Gas particles
! 2. Sinks which will be treated as sources
! 3. Stars which are discarded by default but this could be changed if required
! 4. Everything else (ghost, dead/unknown). These will be discarded
    
    ! itype: 1 - gas particles
    ipType = indexWord("gas",pType,npType)
    if (ipType/=0) then
       npart = pNumArray(ipType)
    else
       npart = 0
    endif

    ! itype: 3 - sinks, treat as sources
    if (gadget) then
       ipType = indexWord("sink / black",pType,npType)
       if (ipType/=0) then
          nptmass = pNumArray(ipType)
       else
          nptmass = 0
       end if
    else
       ipType = indexWord("sink",pType,npType)
       if (ipType/=0) then
          nptmass = pNumArray(ipType)
       else
          nptmass = 0
       end if
    endif

    ! itype: 4 - stars
    ipType = indexWord("star",pType,npType)
    if (ipType/=0) then
       nstar = pNumArray(ipType)
    else
       nstar = 0
    end if

    ! Everything else is discarded
    nother=SUM(pNumArray(:))-npart-nptmass-nstar

! Read in units
! Initialise the unit string so that printing the itype unit doesn't give junk
    unit(:) = "No unit            "
    unitname(:) = "No unit            "
    read(LUIN,*)
    read(LUIN,'(a)') unitString
    unitString = unitstring(2:)
    read(LUIN,'(a)') unitnameString
    unitnameString = unitnamestring(2:)
    read(LUIN,*)
    read(LUIN,'(a)') nameString
    nameString = nameString(2:)

    call splitIntoWords(unitnameString, unitname, nunit)
    call splitintoWords(unitString, unit, nunit)
    call splitIntoWords(nameString, word, nWord)
    write(*,*) nword," Words: "
    do i = 1, nWord
       write(*,*) trim(word(i)),": ",trim(unit(i)), " ",trim(unitname(i))
    enddo

! Gadget doesn't have itype so just set iitype to zero explicitly to avoid 
! a warning when calling indexWord
    if (gadget) then 
       iitype = 0
    else
       iitype = indexWord("itype",word,nWord)
    endif

    if ( ((iiType==0).and.(nUnit /= nWord)).or. &
         ((iiType/=0).and.(nUnit /= (nWord-1))) )then
       call writeFatal("Error reading splash file")
       write(*,*) nUnit, nWord
       write(*,*) "Words: "
       do i = 1, nWord
          write(*,*) trim(word(i)),": ",trim(unit(i))
       enddo
       stop
    endif
    
    if (internalView.and.multiTime) then
       call rotateForTime
    else
       extraPA = 0.0
    endif

    ix = indexWord("x",word,nWord)
    iy = indexWord("y",word,nWord)
    iz = indexWord("z",word,nWord)

    read(unit(ix),*) uDist
    uDist =  toTorusUnits(unitName(ix)) * uDist
    write(*,*) "!!!!udist ",udist
! Check how the velocity columns are labelled, there are two possibilities so handle both
    if ( wordIsPresent("v\dx",word,nWord) ) then
       ivx = indexWord("v\dx",word,nWord)
       ivy = indexWord("v\dy",word,nWord)
       ivz = indexWord("v\dz",word,nWord)
    else if ( wordIsPresent("v_x",word,nWord) ) then
       ivx = indexWord("v_x",word,nWord)
       ivy = indexWord("v_y",word,nWord)
       ivz = indexWord("v_z",word,nWord)
    else
       call writeFatal("Did not find velocity columns in SPH file.")
    end if

    idustfrac = 0
    if (wordIsPresent("dustfrac",word,nWord)) then
       idustfrac = indexWord("dustfrac",word,nWord)
    endif


    read(unit(ivx),*) uvel
    uVel = toTorusUnits(unitName(ivx))

    imass = indexWord("particle mass",word,nWord)

    read(unit(imass),*) umass
    umass = toTorusUnits(unitName(imass)) * umass

    if(dragon) then
       iu = indexWord("temperature",word,nWord)
    else
       iu = indexWord("u",word,nWord)
    endif
    utemp = 1.d0
!    if (iu /=0) read(unit(iu),*) utemp
    
    irho = indexWord("density",word,nWord)
    ih = indexWord("h",word,nWord)

! A Gadget ASCII dump doesn't contain the physical unit information so hardwire here
    if (gadget) then
       call writeInfo("Resetting physical unit conversion for Gadget dump", FORINFO)
       udist = 3.085678e18_db
       umass = 1.989e33_db
       uvel  = 1.0e5_db
       utime = udist/uvel
       write(message,*) "udist/umass/uvel=", udist, umass, uvel
       call writeInfo(message,FORINFO)
    end if
    
! Check whether this SPH dump has u and u/T  available (e.g. SPH with radiative transfer). 
! If they are both present then set the temperature from these columns. 
    if ( wordIsPresent("u",word,nWord) .and. wordIsPresent("u/T",word,nWord) ) then
       call writeinfo("Found u and u/T in the SPH file",TRIVIAL)
       iUoverT = indexWord("u/T",word,nWord)
       haveUandUoverT = .true.
       if ( wordIsPresent("DustTemperature",word,nWord) ) then
          call writeinfo("Found DustTemperature in the SPH file",TRIVIAL)
          iDustTemperature = indexWord("DustTemperature",word,nWord)
          haveDustTemperature = .true.
       else
          haveDustTemperature = .false.
       endif
    else
       haveUandUoverT = .false.
       haveDustTemperature = .false.
    end if

!
! 1.3 Read in particle data
!
    nlines = npart + nptmass + nstar + nother ! nlines now equal to no. lines - 12 = sum of particles dead or alive

    ! Allocate an array for storing the values we are about to read in from the SPH file
    allocate(junkArray(maxNumVals,nlines))

    ! Zero counters and accumulators
    iptmass = 0
    icount = 12 ! There are 12 header lines
    igas = 0
    istar = 0
    idead = 0
    status=0

    sphdata%totalgasmass = 0.d0
    sphdata%totalHImass  = 0.d0
    sphdata%totalMolmass = 0.d0
    
    ! Read in SPH data and store it in junkArray
    call writeinfo("Reading particle data", TRIVIAL)
    
    do ipart=1, nlines
       icount = icount + 1
       read(LUIN,*,iostat=status) junkArray(1:nWord,ipart)
       if (status /= 0) then
          write(message,*) "Error reading from line ", icount
          call writeFatal(message)
          STOP
       endif
    end do

    ! Finished reading SPH dump so close the file
    write(message,*) "Read ",icount, " lines from file"
    call writeinfo(message, TRIVIAL)
    call writeinfo("Closing ASCII file", TRIVIAL)
    close(LUIN)

    ! 1.4 Work out which particles need to be stored
    if ( sphdensitycut > 0 ) then 
       irequired=0
       irejected=0
       do ipart=1, nlines
          if ( iitype == 0 ) then
             itype=get_gadget_itype(ipart, npart, nptmass) 
          else
             itype = int(junkArray(iitype,ipart))
          endif
          ! Only consider gas particles
          if ( itype /= 1 ) cycle
          rhon = junkArray(irho,ipart)
          if ( rhon > sphdensitycut ) then
             irequired=irequired+1
          else
             irejected=irejected+1
          end if
       end do

       write(message,*) "Applied density cut. ", irequired, "particles retained and ", irejected, " rejected"
       call writeinfo(message, TRIVIAL)
       if ( irequired==0 ) then
          call writeFatal("No particles remaining after density cut")
          call torus_abort
       endif
       
    ! Set npart to the number of particles required so that arrays are allocated to correct size
       npart = int(irequired)
       
    else
       ! No density cut applied
       irequired=npart
       irejected=0
    endif

! Say what is going to be stored. In init_sph_data sphdata%npart is set to npart and gas particle arrays are 
! allocated to be sphdata%npart in size. The point mass arrays are allocated as nptmass in size. 
    if (discardSinks) then 
       write(message,*) "Allocating ", npart, " gas particles"
    else
       write(message,*) "Allocating ", npart, " gas particles and ", nptmass, " sink particles"
    endif
    call writeinfo(message, TRIVIAL)

! Allocate storage for SPH data     
    call init_sph_data(udist, umass, utime, time, nptmass, uvel, utemp)
    ! velocity unit is derived from distance and time unit (converted to seconds from years)
    write(*,*) "uvel ",uvel
    sphdata%codeVelocitytoTORUS = uvel / cspeed
    write(*,*) "code velocity to torus ",uvel/cspeed

!
! Chemistry section: allocates memory for chemistry data if required
!
    
! Decide how to convert between internal energy and temperature. If there is information about
! the H2 fraction then use this to calculate the molecular weight later on.
    if (convertRhoToHI) then
       write (message,'(a,i2)') "Will convert total density to HI density using H2 fraction from column ", sphh2col
       call writeInfo(message,FORINFO)
       write (message,'(a)') "Will calculate molecular weight based on H2 fraction"
       call writeInfo(message,FORINFO)
       sphdata%codeEnergytoTemperature = 1.0
    else if (sphWithChem) then
       write(message,'(a)') "Will read chemistry data but will not convert total density to HI density"
       call writeInfo(message,FORINFO)
       sphdata%codeEnergytoTemperature = 1.0
    else if (haveUandUoverT) then
       call writeInfo("Calculating temperature from u and u/T",FORINFO)
       sphdata%codeEnergytoTemperature = 1.0
    else
       if(dragon) then
          sphdata%codeEnergytoTemperature = utemp
       else
          sphdata%codeEnergytoTemperature = utemp * 1.9725e-8 ! temperature from molcluster! 2. * 2.46 * (u * 1d-7) / (3. * 8.314472)
       endif
       write(message,*) "Conversion factor between u and temperature (assumes molecular weight of 2.46): ", &
            sphdata%codeEnergytoTemperature
       call writeInfo(message, FORINFO)
    end if

    if (sphWithChem) then
       write (message,'(a,i2)') "Reading molecular abundance from column ", sphmolcol
       call writeInfo(message,FORINFO)
       if (sphmolcol > nWord) then 
          write(message,*) "SPH file only has", nWord, " columns"
          call writeFatal(message)
       endif
       allocate(sphData%rhoCO(npart))
    end if

    if (convertRhoToHI.or.sphwithChem) then 
       write (message,'(a,i2)') "Will store particle H2 density. H2 fraction is from column ", sphh2col
       call writeInfo(message,FORINFO)
       allocate(sphdata%rhoH2(npart))
    end if

    if (idustfrac/=0) then
       ALLOCATE(sphdata%dustfrac(sphdata%npart))
    endif

    
! Populate Torus's SPH data structure
part_loop: do ipart=1, nlines
       
       xn = junkArray(ix,ipart)
       yn = junkArray(iy,ipart)
       zn = junkArray(iz,ipart)
       vx = junkArray(ivx,ipart)
       vy = junkArray(ivy,ipart)
       vz = junkArray(ivz,ipart)

       if (internalView) call rotate_particles(galaxyPositionAngle+extraPA, galaxyInclination)

       u = 0.d0
       if (iu /= 0) then
          u = junkArray(iu,ipart)
          read(unit(iu),*) ufac
          u = u / ufac
       endif
       rhon = junkArray(irho,ipart)
       h = junkArray(ih,ipart)


       dustfrac = 0.1
       if (idustFrac /= 0) then
          dustfrac = junkArray(idustfrac,ipart)
       endif

       if ( iitype == 0 ) then
          itype=get_gadget_itype(ipart, irequired+irejected, nptmass)
       else
          itype = int(junkArray(iitype,ipart))
       endif
       
       gaspartmass = junkArray(imass,ipart)

! 1=gas particle
       !FOR MHD DUMPS MAY WISH TO EXCLUDE LOW DENSITY ENVELOPE
       !if ((itype == 1) .AND. (rhon > 5.e-5))  then
       if (itype == 1) then

          ! Exclude particles with density below mininum threshold
          if ( rhon < sphdensitycut ) cycle
          
          igas = igas + 1

          sphdata%xn(igas) = xn + sphPosOffset%x
          sphdata%yn(igas) = yn + sphPosOffset%y
          sphdata%zn(igas) = zn + sphPosOffset%z


          sphdata%gasmass(igas) = gaspartmass

          sphdata%vxn(igas) = vx + sphVelOffset%x
          sphdata%vyn(igas) = vy + sphVelOffset%y
          sphdata%vzn(igas) = vz + sphVelOffset%z

          if (idustfrac/=0) sphdata%dustfrac = dustfrac

! For SPH simulations with chemistry we need to set up H2
          if ( convertRhoToHI.or.sphwithChem ) then
             h2ratio = junkArray(sphh2col,ipart)
          endif

! Set up temperature or internal energy
          if (convertRhoToHI.or.sphwithChem ) then 
             gmw = (2.*h2ratio+(1.-2.*h2ratio)+0.4) / (0.1+h2ratio+(1.-2.*h2ratio))
             if(dragon) then
                sphdata%temperature(igas) = utemp
             else
                sphdata%temperature(igas) = (2.0/3.0) * u * ( gmw / Rgas) * utemp
             endif
          else if (haveDustTemperature) then
             sphdata%temperature(igas) = junkArray(iDustTemperature,ipart)
          else if (haveUandUoverT) then 
             sphdata%temperature(igas) = u / junkArray(iUoverT,ipart)
          else
             sphdata%temperature(igas) = u
          end if

          sphdata%hn(igas) = h 

          if(rhon .gt. 0.d0) then
             sphdata%rhon(igas) = rhon
          else
             sphdata%rhon(igas) = gaspartmass / h**3
          endif

          if (sphwithChem) then 
             COfrac = junkArray(sphmolcol,ipart)
             sphdata%rhoCO(igas) = COfrac * rhon
          endif

          if (convertRhoToHI.or.sphwithChem) then
             sphdata%rhoH2(igas) = h2ratio*2.*rhon*5./7.
          end if

! Calculate total gas mass summed over all particles and (if required) HI and CO masses. 
          sphdata%totalgasmass = sphdata%totalgasmass + gaspartmass
          if (convertRhoToHI.or.sphwithChem) then
             sphdata%totalHImass  = sphdata%totalHImass + (1.0-2.0*h2ratio)*gaspartmass*5.0/7.0
          endif
          if (sphwithChem) then
             sphdata%totalMolmass = sphdata%totalMolmass + COfrac*gaspartmass
          endif
          
! 3=sink
       else if(itype .eq. 3) then

! Update counter before cycling loop so that we count the sinks correctly even if we don't use them
          iptmass = iptmass + 1
          if (discardSinks) cycle part_loop

          sphdata%x(iptmass) = xn
          sphdata%y(iptmass) = yn
          sphdata%z(iptmass) = zn

          sphdata%vx(iptmass) = vx
          sphdata%vy(iptmass) = vy
          sphdata%vz(iptmass) = vz

          sphdata%ptmass(iptmass) = gaspartmass 
          sphdata%hpt(iptmass) = h

!          sphdata%totalgasmass = sphdata%totalgasmass + gaspartmass

! This write statement is checked by the test suite so update checkSphToGrid.pl if it changes.
          write(message,*) "Sink Particle number", iptmass," - mass", gaspartmass, " Msol - Index", iptmass + igas
          call writeinfo(message, TRIVIAL)
          write(98,*) iptmass, xn*udist*1e-10, yn*udist*1e-10, zn*udist*1e-10, gaspartmass

! 4=stars
       else if(itype .eq. 4) then

          istar = istar + 1

! 2=ghost, 5=dead/unknown
       else

          idead = idead + 1
       
       endif

    enddo part_loop

 call writeInfo("Summary of SPH particles read: ", TRIVIAL)
 write(message,*) igas," are gas particles"
 call writeinfo(message, TRIVIAL)
 write(message,*)  iptmass," are sink particles"
 call writeinfo(message, TRIVIAL)
  write(message,*) istar," are star particles"
 call writeinfo(message, TRIVIAL)
 write(message,*) idead, " are unknown/dead"
 call writeinfo(message, TRIVIAL)
 write(message,*) "Total Mass in gas particles: ", sphdata%totalgasmass * umass/mSol, " Msol"
 call writeinfo(message, TRIVIAL)

! Warn if we didn't find the expected number of particles or point masses
 if (igas /= npart ) then
    write(message,*) "Expected ", npart, " gas particles but found ", igas
    call writeWarning(message)
 endif
 
 if (iptmass /= nptmass) then
    write(message,*) "Expected ", nptmass+nstar, " sinks but found ", iptmass
    call writeWarning(message)
 endif

 if (istar /= nstar) then
    write(message,*) "Expected ", nstar, " stars but found ", istar
    call writeWarning(message)
 endif

 if (idead /= nother ) then
    write(message,*) "Expected ", nother, " unknown/dead particles but found ", idead
    call writeWarning(message)
 endif

 deallocate(junkArray)
   
 contains

   subroutine rotate_particles(thisPositionAngle, thisInclination)
     real(double), intent(in) ::  thisPositionAngle, thisInclination
     TYPE(VECTOR) :: orig_sph, rot_sph

     orig_sph = VECTOR(xn,  yn,  zn)
     rot_sph  = rotateZ( orig_sph,  thisPositionAngle*degToRad )
     rot_sph  = rotateY( rot_sph, thisInclination*degToRad )
     xn = rot_sph%x; yn=rot_sph%y; zn=rot_sph%z

     orig_sph = VECTOR(vx,  vy,  vz)
     rot_sph  = rotateZ( orig_sph,  thisPositionAngle*degToRad )
     rot_sph  = rotateY( rot_sph, thisInclination*degToRad )
     vx = rot_sph%x; vy=rot_sph%y; vz=rot_sph%z

   end subroutine rotate_particles


   subroutine rotateForTime
     implicit none
     real(double), parameter :: refTime = 5.3010001E+00 * 4.7165268E+07
     real(double), parameter :: omega   = 1.0e-15 ! (rad s^-1)
     real(double) :: deltaT

     write(message,*) "Time of this dump: ", time*utime/1.0e6, " Myr"
     call writeInfo(message,FORINFO)
     write(message,*) "Reference time: ", refTime/1.0e6, " Myr"
     call writeInfo(message,FORINFO)
     deltaT = refTime - time*utime
     write(message,*) "Time offset: ", deltaT/1.0e6, " Myr"
     call writeInfo(message,FORINFO)
     write(message,*) "Omega: ", omega, " rad s^-1"
     call writeInfo(message,FORINFO)
     extraPA = omega * (deltaT /secsToYears)  * radToDeg
     write(message,*) "Angle offset: ", extraPA, " degrees"
     call writeInfo(message,FORINFO)

   end subroutine rotateForTime
   
  end subroutine new_read_sph_data

   ! Return an integer which specifies the particle type for Gadget.
   ! Gadget output doesn't include this information so we work it out based
   ! on the order of particles in the file. Assume splash has written out the
   ! particles in the order gas, sinks, others. 
   ! Moved to function (D. Acreman, March 2017)
  integer function get_gadget_itype(this_index, this_npart, this_nptmass)

    integer, intent(in) :: this_index   ! Index of this particle
    integer, intent(in) :: this_npart   ! number of gas particles (for gadget)
    integer, intent(in) :: this_nptmass ! number of point masses  (for gadget)
     
    if (this_index<=this_npart) then
       get_gadget_itype = 1
    else if (this_index <= this_npart+this_nptmass) then 
       get_gadget_itype = 3
    else
       get_gadget_itype = 5
    end if

  end function get_gadget_itype
  
!---------------------------------------------------------------------------------------------------------------------------------

! Read data from a Gadget2 snapshot file.
! D. Acreman September 2014

  subroutine read_gadget2_data(filename)
    use inputs_mod, only: convertRhoToHI, sphwithChem, discardSinks
    use parallel_mod, only: torus_abort

    character(len=*), intent(in) :: filename
    integer, parameter  :: LUIN = 10 ! logical unit # of the data file
    integer :: i, igas, iSink
    integer :: nptmass, nTotal, nDead, nGasTotal
    integer :: foundnDead, iStart, iEnd
    character(len=80) :: message
    real(double) :: udist, umass, utime, uvel, utemp 
    real(double) :: numDen, rhoHI, h2ratio
    integer, parameter :: sinkType=6 ! Which particle type to treat as sinks?
! Gadget data
    integer, parameter :: nPartType=6  ! Number of particle types in Gadget
    integer :: g_npart(nPartType) ! Number of particles of each type in this file
    integer :: g_nall(nPartType)  ! Total number of particles of each type in the simulation
    integer :: g_nFiles           ! Number of files in each snapshot
    integer :: g_nVarMass         ! Number of particles with variable mass 
    real(double) :: g_massTable(nPartType), time
    real(single), allocatable :: g_pos(:,:), g_vel(:,:), g_m(:), g_u(:), g_rho(:), g_h(:)
    real(single), allocatable :: aH(:), aH2(:), aCO(:)
    integer, allocatable      :: g_idNum(:)
    character(len=4) :: blockLabel

    logical :: isBlockLabelled ! Does this dump use block labelled format?
    logical :: foundH, foundH2
    integer :: myiostat

! Dummy variables for reading the header
    real(double) :: dummyDouble
    integer      :: dummyInt, dummyInt6(6)

!
! 1. Read the file header and perform some checks
!
    call writeInfo("",FORINFO)
    write(message,*) "Reading Gadget dump from "//trim(filename)
    call writeInfo(message,FORINFO)
    open(unit=LUIN,status="old",file=filename,form="unformatted")

! Check if this dump uses block labelled format or standard format
    read(LUIN) blockLabel
    if ( blockLabel=="HEAD" ) then
       isBlockLabelled=.true.
       write(message,*) "This is a block labelled dump"
       call writeInfo(message,TRIVIAL)
    else
       isBlockLabelled=.false.
       rewind(LUIN)
       write(message,*) "This is not a block labelled dump"
       call writeInfo(message,TRIVIAL)
    end if

! Read the header
    read (LUIN) g_npart, g_massTable, time, dummyDouble, dummyInt, dummyInt, g_nall, &
         dummyInt, g_nFiles, dummyDouble, dummyDouble, dummyDouble, dummyDouble, &
         dummyInt, dummyInt, dummyInt6, utime, umass, udist

    call writeInfo("",TRIVIAL)
    write(message,*) "Gadget particle list: "
    call writeInfo(message,TRIVIAL)
    write(message,*) "      Gas      Halo      Disk     Bulge     Stars  Boundary"
    call writeInfo(message,TRIVIAL)
    write(message,"(6(1x,i9))") g_npart(:)
    call writeInfo(message,TRIVIAL)
    call writeInfo("",TRIVIAL)

! Check that there is only one file in this snapshot ...
    if ( g_nFiles /= 1 ) then
       write(message,*) "Multiple files per snapshot is not supported ", g_nFiles
       call writeFatal(message)
    endif

! ... which means that the number of particles in this file is the same as the total number of particles in the simulation.
    do i=1, nPartType
       if (g_npart(i) /= g_nall(i)) then
          write(message,*) "This file does not contain all the particles ", i, g_npart(i), g_nall(i)
          call writeWarning(message)
          exit
       end if
    end do

!
! 2. Set up units and report the values used
!    Storing units in the header is not in the public Gadget2 but any problems will shown up here
!
    uvel  = udist/utime
    utemp = (udist/utime)**2 ! Units of internal energy
    sphdata%codeEnergytoTemperature = 1.0 ! Do the conversion in this routine
    sphdata%codeVelocitytoTORUS = uvel / cspeed

    call writeInfo("Units: ",TRIVIAL)
    write(message,'(a,es11.4,a,es11.4,a)') "Time units (from header):     ", utime, " s  = ", utime*secsToYears, " years"
    call writeInfo(message,TRIVIAL)
    write(message,'(a,es11.4,a,es11.4,a)') "Mass units (from header):     ", umass, " g  = ", umass/mSol, " mSol"
    call writeInfo(message,TRIVIAL)
    write(message,'(a,es11.4,a,es11.4,a)') "Distance units (from header): ", udist, " cm = ", udist/pctocm, " pc"
    call writeInfo(message,TRIVIAL)
    call writeInfo("",TRIVIAL)

    if (utime==0.0 .or. umass==0.0 .or. udist==0.0 ) then
       call writeFatal("Found zero in code units")
    endif
!
! 3. Read in particle information from the dump file
!
    ntotal     = sum(g_npart) ! Total number of particles 
    nGasTotal  = g_npart(1)   ! Total number of gas particles dead or alive

! Count the number of particles with variable masses
    g_nVarMass = 0
    do i=1, nPartType
       if (g_massTable(i) == 0.0 ) g_nVarMass = g_nVarMass + g_npart(i)
    end do
    write(message,'(a,i9)') "Number of particles with variable masses: ", g_nVarMass
    call writeInfo(message,TRIVIAL)

! Particle positions: total number of particles in file
    if (isBlockLabelled) then
       read(LUIN) blockLabel
       write(message,*) "Reading positions block: ", blockLabel
       call writeInfo(message,TRIVIAL)
    endif
    allocate(g_pos(3,nTotal))
    read(LUIN) g_pos

! Particle velocities: total number of particles in file
    if (isBlockLabelled) then
       read(LUIN) blockLabel
       write(message,*) "Reading velocities block: ", blockLabel
       call writeInfo(message,TRIVIAL)
    endif
    allocate(g_vel(3,nTotal))
    read(LUIN) g_vel

! ID numbers: total number of particles in file
    if (isBlockLabelled) then
       read(LUIN) blockLabel
       write(message,*) "Reading IDs block: ", blockLabel
       call writeInfo(message,TRIVIAL)
    endif
    allocate(g_idNum(nTotal))
    read(LUIN) g_idNum

! Particle masses: number of particles with variable mass
    if (isBlockLabelled) then
       read(LUIN) blockLabel
       write(message,*) "Reading mass block: ", blockLabel
       call writeInfo(message,TRIVIAL)
    endif
    allocate(g_m(g_nVarMass))
    read(LUIN) g_m

! Internal energy: number of gas particles
    if (isBlockLabelled) then
       read(LUIN) blockLabel
       write(message,*) "Reading energy block: ", blockLabel
       call writeInfo(message,TRIVIAL)
    endif
    allocate(g_u(nGasTotal))
    read(LUIN) g_u

! Density: number of gas particles
    if (isBlockLabelled) then
       read(LUIN) blockLabel
       write(message,*) "Reading density block: ", blockLabel
       call writeInfo(message,TRIVIAL)
    endif
    allocate(g_rho(nGasTotal))
    read(LUIN) g_rho

! Smoothing length: number of gas particles
    if (isBlockLabelled) then
       read(LUIN) blockLabel
       write(message,*) "Reading smoothing length block: ", blockLabel
       call writeInfo(message,TRIVIAL)
    endif
    allocate(g_h(nGasTotal))
    read(LUIN) g_h

! Optionally read the chemistry components
    if (convertRhoToHI.or.sphwithChem ) then
       write (message,'(a,i2)') "Will store particle H2 density"
       call writeInfo(message,FORINFO)
       allocate(sphdata%rhoH2(nGasTotal))
       allocate(aH(nGasTotal))
       allocate(aH2(nGasTotal))

       if (isBlockLabelled) then
          foundH=.false.
          foundH2=.false.
          do
             read(LUIN,iostat=myiostat) blockLabel
             if (myiostat<0) then
                call torus_abort("Reached end of file and didn't find H and/or H2 data")
             else if (myiostat>0) then
                call torus_abort("Error reading block label")
             end if

             select case (blockLabel)

             case('nH')
                write(message,*) "Reading H block: ", blockLabel
                call writeInfo(message,TRIVIAL)
                read(LUIN) aH
                foundH=.true.

             case('nH2')
                write(message,*) "Reading H2 block: ", blockLabel
                call writeInfo(message,TRIVIAL)
                read(LUIN) aH2
                foundH2=.true.

             case default
                write(message,*) "Ignoring block: ", blockLabel
                call writeInfo(message,TRIVIAL)
                read(LUIN)

             end select

             if (foundH.and.foundH2) exit
          end do

       else
! If the blocks are not labelled then assume the next two blocks H, H2 respectively
          read(LUIN) aH
          read(LUIN) aH2
       end if

    end if

    if (sphWithChem) then
       write (message,'(a,i2)') "Will store CO fraction"
       call writeInfo(message,FORINFO)

       if (isBlockLabelled) then
          read(LUIN) blockLabel
          write(message,*) "Reading CO block: ", blockLabel
          call writeInfo(message,TRIVIAL)
       endif
       allocate(sphData%rhoCO(nGasTotal))
       allocate(aCO(nGasTotal))
       read(LUIN) aCO
    end if

! At this point we're done reading the Gadget file
    close(LUIN)

!
! 4. Determine some more information about particle numbers
!

! Count the number of dead gas particles ...
    nDead = 0
    do i=1,nGasTotal
       if ( g_idNum(i) < 0 ) nDead=nDead+1
    end do

! ... and set npart to the number of live gas particles
    npart = g_npart(1) - nDead

! Component number 'sinkType' will be treated as point masses
    nptmass = g_npart(sinkType)

    write(message,'(a,i9,a,i9,a,i9,a)') "Expecting ", npart, " active gas particles and ", nDead, " dead particles"
    call writeInfo(message,TRIVIAL)

!
! 5. Set up the SPH data strucutre and populate it
!

! Initialiase SPH data structure 
    call init_sph_data(udist, umass, utime, time, nptmass, uvel, utemp)

! Zero the particle type counters
    igas = 0; foundnDead=0
    sphData%totalGasMass=0.0; sphData%totalGasMass=0.0; sphData%totalMolMass=0.0

! Gas particles are stored first
gaspart: do i=1, nGasTotal

! Dead particles have negative id numbers 
       if ( g_idNum(i) < 0 ) then 

          foundnDead = foundnDead + 1

       else

          igas = igas + 1

          if (igas > npart ) then
             write(message,*) "Error igas > npart"
             call writeFatal(message)
             write(message,*) "foundnDead, igas, nGasTotal=", foundnDead, igas, nGasTotal
             call writeFatal(message)
          end if

          sphdata%xn(igas) = g_pos(1,i)
          sphdata%yn(igas) = g_pos(2,i)
          sphdata%zn(igas) = g_pos(3,i)

          sphdata%vxn(igas) = g_vel(1,i)
          sphdata%vyn(igas) = g_vel(2,i)
          sphdata%vzn(igas) = g_vel(3,i)

! Particle masses are either from a data block (variable masses) or read from the header (fixed value for all particles)
          if (g_massTable(1) == 0.0 ) then
             sphdata%gasmass(igas) = g_m(i)
          else
             sphdata%gasmass(igas) = g_massTable(1)
          end if
          sphData%totalGasMass = sphData%totalGasMass + sphdata%gasmass(igas)

          sphdata%temperature(igas) = g_u(i) * utemp * (2.0/3.0) * ( gmw / Rgas)
          sphdata%rhon(igas)        = g_rho(i)
! Halve the smoothing lengths as the gadget definition is different to other SPH codes
          sphdata%hn(igas)          = 0.5 * g_h(i)
          
! Calculate molecular densities. aH2/aCO are fractions by number
          if (convertRhoToHI.or.sphwithChem ) then
             numDen = sphdata%rhon(igas) / (gmw*amu) ! number density of H,H2+He
! aH2 is the H2 fraction by number relative to total gas number density
             sphdata%rhoH2(igas) = aH2(i)*numDen * 2.0*amu
! Next two lines follow the calculation of HI density in clusterParameter
             h2ratio       = 0.5 * (7.0/5.0) *  sphdata%rhoH2(igas) / sphdata%rhon(igas)
             rhoHI         = (1.0-2.0*h2ratio)*sphdata%rhon(igas)*5.0/7.0
! Accumulate total HI mass
             sphdata%totalHImass = sphdata%totalHImass + &
!                 Total gas mass        * (fraction of the mass which is HI)
                  sphdata%gasmass(igas) * (rhoHI / sphdata%rhon(igas))
! CO
             if (sphwithChem) then
                sphdata%rhoCO(igas) = aCO(i)*numDen * 28.0*amu
                sphData%totalMolMass = sphData%totalMolMass + &
!                    Total gas mass        * (fraction of the mass which is CO)
                     sphdata%gasmass(igas) * (sphdata%rhoCO(igas)/sphdata%rhon(igas))
             endif
          endif


       endif

    end do gaspart

    if (igas == npart .and. foundnDead == nDead ) then
       write(message,'(a,i9,a,i9,a)') "Found     ", igas, " active gas particles and ", foundnDead, " dead particles"
       call writeinfo(message, TRIVIAL)
    else
       call writeFatal("Did not find the expected number of live and dead gas particles")
       write(message,'(a,i9,a,i9,a)') "Found     ", igas, " active gas particles and ", foundnDead, " dead particles"
       call writeFatal(message)
    endif
    write(message,*) "Total gas mass in active particles: ", sphData%totalGasMass*umass/mSol, " solar masses"
    call writeInfo(message,TRIVIAL)
    call writeInfo("",TRIVIAL)

!
! Sink particles
!
    if (.not.discardSinks.and.nptmass /= 0) then

       write(message,*) "Reading ", nptmass, " point masses from particle type ", sinkType 
       call writeInfo(message, FORINFO)
       iStart = sum(g_npart(1:sinkType-1))+1
       iEnd   = sum(g_npart(1:sinkType))

       iSink=0
       sinks: do i=iStart, iEnd
          iSink=iSink+1

          sphdata%x(iSink) = g_pos(1,i)
          sphdata%y(iSink) = g_pos(2,i)
          sphdata%z(iSink) = g_pos(3,i)
          
          sphdata%vx(iSink) = g_pos(1,i)
          sphdata%vy(iSink) = g_pos(2,i)
          sphdata%vz(iSink) = g_pos(3,i)

          if (g_massTable(sinkType) == 0.0 ) then
             sphdata%ptmass(iSink) = g_m(i)
          else
             sphdata%ptmass(iSink) = g_massTable(sinkType)
          end if
       end do sinks

       if ( iSink /= nptmass ) then 
          write(message,*) "Expected ", nptmass, "point masses but found ", iSink
          call writeWarning(message)
       endif

    endif

! 6. Wrap up and exit
    deallocate(g_pos, g_vel, g_m, g_u, g_rho, g_h)
    if (allocated(aH))  deallocate (aH)
    if (allocated(aH2)) deallocate (aH2)
    if (allocated(aCO)) deallocate (aCO)

  end subroutine read_gadget2_data

!Read in clumpfind SPH style
! Cass Hall May 2015
!--------------------------------------------------
  subroutine read_sph_data_clumpfind(filename)
    use inputs_mod, only: amrgridcentrex, amrgridcentrey, amrgridcentrez, amrgridsize, splitovermpi, discardSinks,&
         convertrhotohi, sphWithChem
    use angularImage_utils, only:  internalView, galaxyPositionAngle, galaxyInclination
#ifdef MPI
    use mpi
#endif
    implicit none

    
    character(len=*), intent(in) :: filename
    character(LEN=150) :: message

    integer(kind=8) :: iiigas, irejected, irequired, i, j, iiisink
    integer :: iblock ! MPI block number
    integer :: nums(8), numssink(8), numsrt(8), numsmhd(8)

    real(double) :: udist, umass, utime
    real(double) :: time
    integer, parameter :: nptmass=0

    INTEGER(kind=4)  :: int1, int2, i1
    integer(kind=4)  :: number,n1,n2,nreassign,naccrete,nkilltot,nblocks,nkill
    integer(kind=4)  :: nblocktypes
    REAL(kind=8)     :: r1
    integer(kind=8)  :: blocknpart, blocknptmass, blocknradtrans, blocknmhd, blocksum_npart
    integer  :: blocksum_nptmass, i_pt_mass
    CHARACTER(len=100) ::  fileident

    integer(kind=1), parameter  :: LUIN = 10 ! logical unit # of the data file

    integer(kind=1), allocatable :: iphase(:)
    integer, allocatable         :: isteps(:)
    real(kind=8), allocatable    :: xyzmh(:,:)
    real(kind=8), allocatable    :: vxyzu(:,:)
    real(kind=8), allocatable    :: uoverTarray(:)
    real(kind=4), allocatable    :: rho(:)
    real(kind=8), allocatable    :: h2ratio(:)
    real(kind=8), allocatable    :: COfrac(:)

    integer,allocatable :: listpm(:)
    real(kind=8),allocatable    :: spinx(:)
    type(VECTOR) :: centre
    real(double) :: halfSize

    integer :: iThread
    integer :: loopIndex
! The MPI communication handles multiple hydro threads so I've keyed it out for non-hydro builds
! For Galactic plane survey runs build with hydro=no
#ifdef MPI
#ifdef HYDRO
    integer :: ierr
    integer, parameter :: tag = 33
    integer :: status(MPI_STATUS_SIZE)
#endif
#endif
    real(double) :: gmwWithH2

#ifdef MPI
#ifdef HYDRO
    if (myrankWorldGlobal == 0 .or. (.not.splitOverMPI)) then
#endif
#endif

! Part 1: Read in the dump
    write(message,*) "Reading SPH data from "//trim(filename)
    call writeinfo(message, FORINFO)

!---Open file
    open(unit=LUIN, file=filename, form='unformatted', status='old') ! removed little endian statement
    
    read(LUIN) int1, r1, int2, i1, int1
    read(LUIN)  fileident
    write(message,*) "fileident=", fileident
    call writeinfo(message, TRIVIAL)
    
read(LUIN) number
    IF (number==6) THEN
       READ (LUIN) npart,n1,n2,nreassign,naccrete,nkill
       nblocks = 1
    ELSE
       read(LUIN) npart, n1, n2, nreassign, naccrete, nkilltot, nblocks
    ENDIF 
! Read in integers of size int*1, int*2, int*4, int*8
      do i = 1, 8
         read (10) number
      end do

    read(LUIN) udist, umass, utime
    read(LUIN) number

    nblocktypes = number/nblocks
    if (number.ne.nblocktypes*nblocks) then
       write(message,*) "Error in number of block types", number, nblocktypes, nblocks
       call writeFatal(message)
       stop
    endif

! Allocate storage for data read from dump file
    allocate( isteps(npart) )
    allocate( iphase(npart) )
    allocate( xyzmh(5,npart))
    allocate( vxyzu(4,npart))
!---Allocate uoverT for rad trans part - may not be present
    allocate( uoverTarray(npart))
    allocate( rho(npart)    )
    if (ConvertRhoToHI .or. sphWithChem) then
       write (message,'(a)') "Will read H2 fraction"
       call writeInfo(message,TRIVIAL)
       allocate ( h2ratio(npart) )
    endif
    if (sphWithChem) then
       write (message,'(a)') "Will read CO abundance"
       call writeInfo(message,TRIVIAL)
       allocate ( COfrac(npart) )
    endif

    iiigas = 0
    blocksum_npart = 0
    blocksum_nptmass = 0 
    mpi_blocks: do iblock=1, nblocks

       write(message,*) "Reading MPI block ", iblock
       call writeInfo(message,TRIVIAL)

!------Read block containing particles
       read(LUIN) blocknpart, (nums(i), i=1,8)
       write(message,*) "This block has npart=", blocknpart
       call writeInfo(message,TRIVIAL)
       blocksum_npart = blocksum_npart + blocknpart

!------Read block containing pointmasses
       read(LUIN) blocknptmass, (numssink(i), i=1,8)
       write(message,*) "This block has nptmass=", blocknptmass
       call writeInfo(message,TRIVIAL)
       blocksum_nptmass = blocksum_nptmass + int(blocknptmass)

!------Read radiative transfer blocks if present
       if (nblocktypes.GE.3) then
          read(LUIN) blocknradtrans, (numsrt(i), i=1,8)
          write(message,*) "Found RT blocks ", blocknradtrans
          if (blocknradtrans.ne.blocknpart) then
             write (message,*) "blocknradtrans.ne.blocknpart ",blocknradtrans,blocknpart
             call writeFatal(message)
             stop
          endif
          call writeInfo(message,TRIVIAL)
       endif

!------Read MHD blocks if present
       if (nblocktypes.GE.4) then
          read(LUIN) blocknmhd, (numsmhd(i), i=1,8)
          write(message,*) "Found MHD blocks ", blocknmhd
          call writeInfo(message,TRIVIAL)
       endif

!------Read default integers
       READ (LUIN) (isteps(i), i=iiigas+1, iiigas+blocknpart)

!------Read integer*1s
       READ (LUIN) (iphase(i), i=iiigas+1, iiigas+blocknpart)

!------Read default reals
       do j=1,5
          read(LUIN) ( xyzmh(j,i), i=iiigas+1,iiigas+blocknpart)
       end do
       do j=1,4
          read(LUIN) ( vxyzu(j,i), i=iiigas+1,iiigas+blocknpart)
       end do

!-----Read density
      read(LUIN) (rho(i), i=iiigas+1,iiigas+blocknpart)
!-----Read dgrav - added by duncan
      read(LUIN) !(dgrav(i), i=1,npart)

!-----Read in some junk
       allocate(listpm(blocknptmass), spinx(blocknptmass))
       READ (LUIN) (listpm(i),i=1,blocknptmass)
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       deallocate(listpm, spinx)

       iiigas = iiigas+blocknpart

    end do mpi_blocks

    close(LUIN)
    if (blocksum_npart == npart ) then 
       write(*,*) "Read ", npart, " particles"
!       call writeInfo(message, TRIVIAL)
    else
       write(message,*) "Read ",blocksum_npart, "particles but expected ", npart
       call writeWarning(message)
    endif

! If the grid is split over MPI processes then only the rank zero process will get here. 
! It then needs to loop over all the hydro threads to hand out the particles. 
! Otherwise all MPI process get here and execute this loop once. 
    if (splitOverMPI) then
       loopIndex = nHydroThreadsGlobal
    else
       loopIndex = 1
    endif

#ifdef MPI
#ifdef HYDRO
hydroThreads: do iThread = 1, loopIndex
#endif
#endif


!
! Part 2: Set up SPH data structure

    centre = VECTOR(amrgridcentrex, amrgridcentrey, amrgridcentrez)
    halfSize = amrgridSize / 2.d0
    if (splitOverMpi) then
       call  domainCentreAndSize(iThread, centre, halfsize)
    endif
    
!

! Work out how many particles we actually need
    irequired=0
    irejected=0
    i_pt_mass =0 
    do i=1, blocksum_npart
! Rotate/tilt for generating Galactic plane surveys
       if (internalView) call rotate_particles_mpi(galaxyPositionAngle, galaxyInclination)

       if (iphase(i) > 0) then
          if ( particleInBox(xyzmh(:,i),uDist, centre, halfSize) .and. (.not.discardSinks)) then 
             i_pt_mass = i_pt_mass + 1
          endif
       endif
       if (iphase(i) == 0) then
          if ( particleRequired(xyzmh(:,i),uDist, centre, halfSize)) then 
             irequired = irequired + 1
          else
             irejected = irejected + 1
          endif
       end if
    end do
! irequired is the number of particles we are going to store in sphData
! int is explicit convertion from kind=8 to default integer to avoid compiler warning
    npart = int(irequired)

    write(*,*) myrankGlobal, "nptmass ", i_pt_mass, npart
    call init_sph_data(udist, umass, utime, time, i_pt_mass)
    sphdata%codeVelocitytoTORUS = (udist / utime) / cspeed
! Allocate chemistry related components
    if (sphwithChem) then
       allocate ( sphData%rhoCO(npart) )
    endif
    if (sphwithChem.or.convertRhoToHI) then
       allocate ( sphData%rhoH2(npart) )
    endif

    if (nblocktypes.GE.3) then
       sphdata%codeEnergytoTemperature = 1.0
    else if (convertRhoToHI) then
       write(message,*) "Will convert density to HI density and set molecular weight accordingly"
       call writeInfo(message,FORINFO)
       sphdata%codeEnergytoTemperature = 1.0
    else if (sphWithChem) then
       write(message,'(a)') "Will set molecular weight using H2 fraction but not convert total density to HI density"
       call writeInfo(message,FORINFO)
       sphdata%codeEnergytoTemperature = 1.0
    else
       sphdata%codeEnergytoTemperature = (2.0/3.0) * (gmw / Rgas) * ( (udist**2)/(utime**2) )
       write(message,*) "Conversion factor between u and temperature (assumes molecular weight of", gmw,"): ", &
            sphdata%codeEnergytoTemperature
       call writeInfo(message,FORINFO)
    endif
    iiigas=0
    iiisink = 0 
    sphdata%totalgasmass = 0.0
    sphdata%totalHImass  = 0.0
    sphdata%totalMolmass = 0.0

! This loop needs to be over all the particles read in
    do i=1, blocksum_npart
       if (iphase(i) == 0) then
          if ( particleRequired(xyzmh(:,i),uDist, centre, halfSize)) then 

             iiigas=iiigas+1

             sphData%rhon(iiigas)  = rho(i)

! Calculate chemistry related components
             if (convertRhoToHI.or.sphWithChem) then 
                sphData%rhoH2(iiigas) = h2ratio(i)*2.0*rho(i)*5./7.
             endif
             if (sphWithChem) then
                sphData%rhoCO(iiigas) = rho(i) * COfrac(i)
             endif

             sphdata%vxn(iiigas)         = vxyzu(1,i)
             sphdata%vyn(iiigas)         = vxyzu(2,i)
             sphdata%vzn(iiigas)         = vxyzu(3,i)
             sphData%temperature(iiigas) = vxyzu(4,i)
! If the radiative transfer block exists, set temperatures as u / (u/T)
             if (nblocktypes.GE.3) then
                if (uoverTarray(i).le.0.0) then
                   call writeFatal("u over T is <=0")
                   write(*,*) "u over t ",i,uOverTarray(i)
                   stop
                endif
                sphData%temperature(iiigas) = sphData%temperature(iiigas) / uoverTarray(i)
             endif
! Convert internal energy to temperature using the mean molecular weight derived from the H2 fraction
! sphdata%codeEnergytoTemperature = 1.0 so this needs to be temperature
             if (convertRhoToHI.or.sphWithChem) then
                gmwWithH2 = (2.*h2ratio(i)+(1.-2.*h2ratio(i))+0.4) / (0.1+h2ratio(i)+(1.-2.*h2ratio(i)))
                sphdata%temperature(iiigas) = (2.0/3.0) * sphData%temperature(iiigas) * ( gmwWithH2 / Rgas) &
                     * ( (udist**2)/(utime**2) )
             endif

             sphData%xn(iiigas)          = xyzmh(1,i)
             sphData%yn(iiigas)          = xyzmh(2,i)
             sphData%zn(iiigas)          = xyzmh(3,i)
             sphData%gasmass(iiigas)     = xyzmh(4,i)
             sphData%hn(iiigas)          = xyzmh(5,i)

             sphdata%totalgasmass = sphdata%totalgasmass + sphData%gasmass(iiigas)
             if (convertRhoToHI.or.sphwithChem) then
                sphdata%totalHImass  = sphdata%totalHImass + (1.0-2.0*h2ratio(i))*sphData%gasmass(iiigas)*5.0/7.0
             endif
             if (sphwithChem) then
                sphdata%totalMolmass = sphdata%totalMolmass + COfrac(i)*sphData%gasmass(iiigas)
             endif
  endif
       else if(iphase(i) > 0) then
          if ( particleInBox(xyzmh(:,i),uDist, centre, halfSize).and. (.not.discardSinks)) then 
       
             iiisink=iiisink+1


             sphdata%x(iiisink)         = xyzmh(1,i)
             sphdata%y(iiisink)         = xyzmh(2,i)
             sphdata%z(iiisink)         = xyzmh(3,i)
             sphData%ptmass(iiisink)    = xyzmh(4,i)
             sphdata%vx(iiisink)         = vxyzu(1,i)
             sphdata%vy(iiisink)         = vxyzu(2,i)
             sphdata%vz(iiisink)         = vxyzu(3,i)
          endif
       endif
    end do
    sphdata%nptmass = int(iiisink)

    write(message,*) "Stored ", irequired, " gas particles which are required"
    call writeinfo(message, TRIVIAL)    

    write(message,*) "Rejected ", irejected, " gas particles which do not influence the grid"
    call writeinfo(message, TRIVIAL) 
    write(message,*) "Total number of gas particles is ", irequired+irejected
    call writeinfo(message, TRIVIAL)    

    write(message,*) "Total number of required sink particles is ", iiisink
    call writeinfo(message, TRIVIAL)    

    write(message,*) "Minimum x=", minval(xyzmh(1,:)) * udist
    call writeinfo(message, TRIVIAL)
    write(message,*) "Maximum x=", maxval(xyzmh(1,:)) * udist
    call writeinfo(message, TRIVIAL)

    write(message,*) "Minimum y=", minval(xyzmh(2,:)) * udist
    call writeinfo(message, TRIVIAL)
    write(message,*) "Maximum y=", maxval(xyzmh(2,:)) * udist
    call writeinfo(message, TRIVIAL)

    write(message,*) "Minimum z=", minval(xyzmh(3,:)) * udist
    call writeinfo(message, TRIVIAL)
    write(message,*) "Maximum z=", maxval(xyzmh(3,:)) * udist
    call writeinfo(message, TRIVIAL)

#ifdef MPI
#ifdef HYDRO

    if (splitOverMPI) then

    call mpi_send(uDist, 1, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(uMass, 1, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(uTime, 1, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphData%codeVelocityToTorus, 1, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)


    call mpi_send(sphdata%nPart, 1, MPI_INTEGER, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%rhon, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)

    call mpi_send(sphdata%vxn, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%vyn, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%vzn, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)

call mpi_send(sphdata%temperature, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)

    call mpi_send(sphdata%xn, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%yn, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%zn, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)

    call mpi_send(sphdata%gasmass, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%hn, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)

    call mpi_send(sphdata%nptmass, 1, MPI_INTEGER, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%ptmass, sphData%nptmass, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)

    call mpi_send(sphdata%x, sphData%nptmass, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%y, sphData%nptmass, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%z, sphData%nptmass, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)

    call mpi_send(sphdata%vx, sphData%nptmass, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%vy, sphData%nptmass, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%vz, sphData%nptmass, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)


    end if

 enddo hydroThreads
    deallocate(rho)
    deallocate(vxyzu)
    deallocate(xyzmh)
    deallocate(uoverTarray)
    deallocate(iphase)
    deallocate(isteps)
    if (ConvertRhoToHI .or. sphWithChem) then
       deallocate (h2ratio)
    endif
    if (sphWithChem) then
       deallocate ( COfrac)
    endif


    else
       call mpi_recv(sphData%uDist, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       call mpi_recv(sphData%uMass, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       call mpi_recv(sphData%uTime, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       call mpi_recv(sphData%codeVelocitytoTorus, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

       call mpi_recv(sphData%nPart, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
       npart = sphData%nPart
       allocate(sphdata%rhon(1:sphData%nPart))
       call mpi_recv(sphData%rhon, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

       allocate(sphdata%vxn(1:sphData%nPart))
       allocate(sphdata%vyn(1:sphData%nPart))
       allocate(sphdata%vzn(1:sphData%nPart))
call mpi_recv(sphData%vxn, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       call mpi_recv(sphData%vyn, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       call mpi_recv(sphData%vzn, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

       allocate(sphdata%temperature(1:sphData%nPart))
       call mpi_recv(sphData%temperature, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

       allocate(sphdata%xn(1:sphData%nPart))
       allocate(sphdata%yn(1:sphData%nPart))
       allocate(sphdata%zn(1:sphData%nPart))
       call mpi_recv(sphData%xn, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       call mpi_recv(sphData%yn, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       call mpi_recv(sphData%zn, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

    allocate(sphdata%gasmass(1:sphData%nPart))
       call mpi_recv(sphData%gasmass, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

       allocate(sphdata%hn(1:sphData%nPart))
       call mpi_recv(sphData%hn, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       write(*,*) myrankglobal, " received ", sphData%npart

       call mpi_recv(sphData%nptmass, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)

       allocate(sphData%ptmass(1:sphData%nptmass))
       call mpi_recv(sphData%ptmass, sphData%nptmass, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

       allocate(sphData%x(1:sphData%nptmass))
       call mpi_recv(sphData%x, sphData%nptmass, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       allocate(sphData%y(1:sphData%nptmass))
       call mpi_recv(sphData%y, sphData%nptmass, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

       allocate(sphData%z(1:sphData%nptmass))
       call mpi_recv(sphData%z, sphData%nptmass, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

       allocate(sphData%vx(1:sphData%nptmass))
       call mpi_recv(sphData%vx, sphData%nptmass, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       allocate(sphData%vy(1:sphData%nptmass))
       call mpi_recv(sphData%vy, sphData%nptmass, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       allocate(sphData%vz(1:sphData%nptmass))
       call mpi_recv(sphData%vz, sphData%nptmass, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)


    endif
    call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
#endif

contains

   subroutine rotate_particles_mpi(thisPositionAngle, thisInclination)
     real(double), intent(in) ::  thisPositionAngle, thisInclination
     TYPE(VECTOR) :: orig_sph, rot_sph

     orig_sph%x = xyzmh(1,i)
     orig_sph%y = xyzmh(2,i)
     orig_sph%z = xyzmh(3,i)
     rot_sph  = rotateZ( orig_sph,  thisPositionAngle*degToRad )
     rot_sph  = rotateY( rot_sph, thisInclination*degToRad )
     xyzmh(1,i) = rot_sph%x
     xyzmh(2,i) = rot_sph%y
     xyzmh(3,i) = rot_sph%z

     orig_sph%x = vxyzu(1,i)
     orig_sph%y = vxyzu(2,i)
     orig_sph%z = vxyzu(3,i)
     rot_sph  = rotateZ( orig_sph,  thisPositionAngle*degToRad )
     rot_sph  = rotateY( rot_sph, thisInclination*degToRad )
     vxyzu(1,i) = rot_sph%x
     vxyzu(2,i) = rot_sph%y
     vxyzu(3,i) = rot_sph%z

   end subroutine rotate_particles_mpi





   
end subroutine read_sph_data_clumpfind
!------------------------------------------

!---------------------------------------------------------------------------------------------------------------------------------

! Read in SPH data from an SPH-NG MPI dump file
! D. Acreman, June 2012
! M.R. Bate, November 2025
!
  subroutine read_sph_data_mpi(filename)
    use inputs_mod, only: amrgridcentrex, amrgridcentrey, amrgridcentrez, amrgridsize, splitovermpi, discardSinks,&
         convertrhotohi, sphWithChem
    use angularImage_utils, only:  internalView, galaxyPositionAngle, galaxyInclination
#ifdef MPI
    use mpi
#endif
    implicit none

    
    character(len=*), intent(in) :: filename
    character(LEN=150) :: message

    integer(kind=8) :: iiigas, irejected, irequired, i, j, iiisink
    integer :: iblock ! MPI block number
    integer :: nums(8), numssink(8), numsrt(8), numsmhd(8)

    real(double) :: udist, umass, utime
    real(double) :: time
    integer, parameter :: nptmass=0

    INTEGER(kind=4)  :: int1, int2, i1, int1o
    integer(kind=4)  :: number,n1,n2,nreassign,naccrete,nkilltot,nblocks,nkill
    integer(kind=4)  :: nblocktypes
    integer(kind=4)  :: iyr,idum,iplanetesimals,irotpot,idragscheme,HY09_ndust_bins,idustFluid,ndusttypes
    integer(kind=8)  :: iuniquemax
    REAL(kind=8)     :: r1
    REAL(kind=4)     :: r4
    integer(kind=8)  :: blocknpart, blocknptmass, blocknradtrans, blocknmhd, blocksum_npart
    integer  :: blocksum_nptmass, i_pt_mass
    CHARACTER(len=100) ::  fileident
    integer(kind=4), parameter  :: nmaxtags = 128
    CHARACTER(len=16)  ::  tagsreal(nmaxtags),tagi
    integer(kind=4)  :: rheader(nmaxtags)
    logical(kind=4)  :: tagged

    integer(kind=1), parameter  :: LUIN = 10 ! logical unit # of the data file
    integer(kind=4), parameter  :: NTAB = 32
    integer(kind=4)  :: iv(NTAB)
    
    integer(kind=1), allocatable :: iphase(:)
    integer, allocatable         :: isteps(:)
    real(kind=8), allocatable    :: xyzmh(:,:)
    real(kind=8), allocatable    :: vxyzu(:,:)
    real(kind=8), allocatable    :: uoverTarray(:)
    real(kind=4), allocatable    :: rho(:)
    real(kind=8), allocatable    :: h2ratio(:)
    real(kind=8), allocatable    :: COfrac(:)
    real(kind=8), allocatable    :: dustfracs(:,:)

    integer,allocatable :: listpm(:)
    real(kind=8),allocatable    :: spinx(:)
    type(VECTOR) :: centre
    real(double) :: halfSize

    integer :: thisNumGas
    integer :: nlistinactive, listinactive(1)
    integer :: iThread
    integer :: loopIndex
! The MPI communication handles multiple hydro threads so I've keyed it out for non-hydro builds
! For Galactic plane survey runs build with hydro=no
#ifdef MPI
#ifdef HYDRO
    integer :: ierr
    integer, parameter :: tag = 33
    integer :: status(MPI_STATUS_SIZE)
#endif
#endif
    real(double) :: gmwWithH2

#ifdef MPI
#ifdef HYDRO
    if (myrankWorldGlobal == 0 .or. (.not.splitOverMPI)) then
#endif
#endif

!
! Part 1: Read in the dump
!

    write(message,*) "Reading SPH data from "//trim(filename)
    call writeinfo(message, FORINFO)

! ifort and gfortran can specify endian conversion in the open statment e.g.
!    open(unit=LUIN, file=filename, form='unformatted', status='old', convert="BIG_ENDIAN")
! but this is not standard and not supported by some compilers so use the appropriate 
! environement variable instead (F_UFMTENDIAN for ifort)
    open(unit=LUIN, file=filename, form='unformatted', status='old')! removed this,convert="LITTLE_ENDIAN")

    read(LUIN) int1, r1, int2, i1, int1o

! Check that the endian is correct. int1o is the same in version 1 and version 2 SPH dumps. 
    if (int1o/=690706) then
       call writeFatal("Problem opening file. Check endian")
    endif

! Check file format: version 2 is incompatible with version 1 and we can't read it yet.
    if (i1==690706) then
       write(message,*) "Version 1 dump opened OK"
       call writeinfo(message, TRIVIAL)
    else if (i1==1) then
       write(message,*) "Version 2 dump: ",int1,r1,int2,i1,int1o
!       call writeFatal("Version 2 dumps not supported yet. Convert to ASCII.")
       call writeinfo(message, TRIVIAL)
    else
       write(message,*) "Unrecognised version number", i1
       call writeWarning(message)
    endif

    read(LUIN) fileident
    write(message,*) "fileident=", fileident
    call writeinfo(message, TRIVIAL)
    if (fileident(2:2).EQ.'T') then
         tagged = .true.
    endif

    read(LUIN) number
    if (number==15+NTAB) then
       read(LUIN)  ! Skip tags
       read(LUIN) npart,n1,n2,nreassign,naccrete,nkill,nblocks,iyr,idum,&
            (iv(i),i=1,NTAB),iplanetesimals,irotpot,idragscheme,&
            HY09_ndust_bins,idustFluid,ndusttypes
    else
       write(message,*) "Unexpected length of header info ",number
       call writeFatal("Fatal: Unexpected length of header info")
    endif
    write(message,*) "Total number of particles= ", npart
    print*,"Total number of particles= ", npart
    call writeInfo(message,TRIVIAL)
    write(message,*) "Number of MPI blocks= ", nblocks
    print*,"Number of MPI blocks= ", nblocks
    call writeInfo(message,TRIVIAL)

! Allocate storage for data read from dump file
    allocate( isteps(npart) )
    allocate( iphase(npart) )
    allocate( xyzmh(5,npart))
    allocate( vxyzu(4,npart))
    allocate( uoverTarray(npart))
    allocate( rho(npart)    )
    allocate( dustfracs(ndusttypes,npart) )
    if (ConvertRhoToHI .or. sphWithChem) then
       write (message,'(a)') "Will read H2 fraction"
       call writeInfo(message,TRIVIAL)
       allocate ( h2ratio(npart) )
    endif
    if (sphWithChem) then
       write (message,'(a)') "Will read CO abundance"
       call writeInfo(message,TRIVIAL)
       allocate ( COfrac(npart) )
    endif

    do i=1,3
       read(LUIN)
    end do
    read(LUIN) number
    write(message,*) "number int8 ", number
    call writeInfo(message,TRIVIAL)
    if (number.eq.1) then
       read(LUIN)
       read(LUIN) iuniquemax
       write(message,*) "iuniquemax = ", iuniquemax
       call writeInfo(message,TRIVIAL)
    endif
    read(LUIN) number
    write(message,*) "number default real = ", number
    call writeInfo(message,TRIVIAL)

    read(LUIN) tagsreal(1:MIN(number,nmaxtags))
    read(LUIN) (rheader(i),i=1,min(number,nmaxtags))
    time = rheader(1)
    write(message,*) "Dump time=", time
    call writeInfo(message,TRIVIAL)

    write(message,*) tagsreal(40)," ",rheader(40) ! HY09_size_min

    read(LUIN) ! Skip number of REAL*4
    read(LUIN) ! Skip number of REAL*8
    read(LUIN) ! Skip tags for units
    read(LUIN) udist, umass, utime
    read(LUIN) number ! This is the number of blocks*nblocktypes

    nblocktypes = number/nblocks

    if (number.ne.nblocktypes*nblocks) then
       write(message,*) "Error in number of block types", number, nblocktypes, nblocks
       call writeFatal(message)
       stop
    endif

    write(message,*) "SPH code units (dist, mass, time): ", udist, umass, utime
    call writeInfo(message,TRIVIAL)

    iiigas = 0
    blocksum_npart = 0
    blocksum_nptmass = 0 
    mpi_blocks: do iblock=1, nblocks

       write(message,*) "Reading MPI block ", iblock
       call writeInfo(message,TRIVIAL)

       read(LUIN) blocknpart, (nums(i), i=1,8)
       write(message,*) "This block has npart=", blocknpart
       call writeInfo(message,TRIVIAL)
       blocksum_npart = blocksum_npart + blocknpart

       read(LUIN) blocknptmass, (numssink(i), i=1,8)
       write(message,*) "This block has nptmass=", blocknptmass
       call writeInfo(message,TRIVIAL)
       blocksum_nptmass = blocksum_nptmass + int(blocknptmass)

       if (nblocktypes.GE.3) then
          read(LUIN) blocknradtrans, (numsrt(i), i=1,8)
          write(message,*) "Found RT blocks ", blocknradtrans
          if (blocknradtrans.ne.blocknpart) then
             write (message,*) "blocknradtrans.ne.blocknpart ",blocknradtrans,blocknpart
             call writeFatal(message)
             stop
          endif
          call writeInfo(message,TRIVIAL)
       endif

       if (nblocktypes.GE.4) then
          read(LUIN) blocknmhd, (numsmhd(i), i=1,8)
          write(message,*) "Found MHD blocks ", blocknmhd
          call writeInfo(message,TRIVIAL)
       endif

       READ (LUIN) ! Skip tags
       READ (LUIN) (isteps(i), i=iiigas+1, iiigas+blocknpart)

       do j=1,nums(1)-1
          READ (LUIN) ! Skip tags
          READ (LUIN) (nlistinactive, listinactive(i), i=1,1)
       END DO
       
       READ (LUIN) tagi ! Skip tag: iphase
       IF (TRIM(tagi).NE.'iphase') THEN
          write (message,*) "Expected tag: iphase ",tagi
          call writeFatal(message)
       ENDIF
       READ (LUIN) (iphase(i), i=iiigas+1, iiigas+blocknpart)

       thisNumGas = 0
       do i=iiigas+1, iiigas+blocknpart
          if (iphase(i)==0) thisNumGas = thisNumGas+1
       end do
       write (message,*) "There are ", thisNumGas, " active gas particles"
       call writeInfo(message,TRIVIAL)

       read(LUIN) ! Skip tag: iunique
       read(LUIN) 

       do j=1,5
          read(LUIN) ! Skip tags: x,y,z,m,h
          read(LUIN) ( xyzmh(j,i), i=iiigas+1,iiigas+blocknpart)
       end do

       do j=1,4
          read(LUIN) ! Skip tags: vx, vy, vz, u
          read(LUIN) ( vxyzu(j,i), i=iiigas+1,iiigas+blocknpart)
       end do

! Read H2 fraction and CO abundance if required
       if ( ConvertRhoToHI .or. sphWithChem ) then
          read(LUIN) ! Skip tags: h2ratio
          READ (LUIN) ( h2ratio(i), i=iiigas+1,iiigas+blocknpart)
          do j=1,3
             read(LUIN) ! Skip tags
             READ (LUIN)
          end do
          if (sphWithChem) then
             read(LUIN) ! Skip tags
             READ (LUIN) ( COfrac(i), i=iiigas+1,iiigas+blocknpart)
          else
             read(LUIN) ! Skip tags
             READ (LUIN)
          end if
          do j=1,nums(6)-15
             read(LUIN) ! Skip tags
             READ (LUIN)
          end do
       endif
       do j=1,nums(6)-10
          read(LUIN) tagi ! Read tags
          if (TRIM(tagi(1:8)).EQ.'dustfrac') GOTO 200
          if (TRIM(tagi(1:5)).EQ.'Dust:') GOTO 200
          READ (LUIN)
       end do
       GOTO 500

! Read one-fluid dust, including multigran
200    write(message,*) "Found dustfrac, ndusttypes=",ndusttypes,tagi(1:10)
       call writeInfo(message,TRIVIAL)
       if (TRIM(tagi(1:5)).NE.'Dust:') read (LUIN)

       do j = 1, ndusttypes
          if (TRIM(tagi(1:5)).NE.'Dust:') then
             read (LUIN) tagi
             if ('dustfrac'.NE.TRIM(tagi)) then
                write(message,*) "ERROR - expected dustfracs ",tagi
                call writeFatal(message)
             endif
          endif
          read (LUIN) (dustfracs(j,i), i=iiigas+1,iiigas+blocknpart)
       end do

500    do j=1,nums(6)
          read(LUIN) tagi ! Read tags
          if (TRIM(tagi).EQ.'rho') GOTO 600
          READ (LUIN)
       end do

       write (message,*) "Could not find tag: rho ",tagi
       call writeFatal(message)

600    read(LUIN) ( rho(i), i=iiigas+1,iiigas+blocknpart) 

       do j=1,nums(7)-1
          READ (LUIN) tagi
          write(message,*) "Skipping ", tagi
          call writeInfo(message,TRIVIAL)
          READ (LUIN)
       end do
!       read (LUIN) (r4, i=iiigas+1,iiigas+blocknpart)

       allocate(listpm(blocknptmass), spinx(blocknptmass))
       READ (LUIN) tagi
       IF (TRIM(tagi).NE.'listpm') THEN
          write(message,*) "ERROR - expected listpm ",tagi
          call writeFatal(message)
       ENDIF
       READ (LUIN) (listpm(i),i=1,blocknptmass)
       READ (LUIN) ! Skip tags
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) ! Skip tags
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) ! Skip tags
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) ! Skip tags
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) ! Skip tags
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) ! Skip tags
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) ! Skip tags
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) ! Skip tags
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       READ (LUIN) ! Skip tags
       READ (LUIN) (spinx(i),i=1,blocknptmass)
       deallocate(listpm, spinx)

       DO i = 1, numssink(6)-10
          READ (LUIN) ! Skip tags
          READ (LUIN)
       END DO
       DO i = 1, numssink(8)
          READ (LUIN) ! Skip tags
          READ (LUIN)
       END DO

       if (nblocktypes.GE.3) then
          call writeInfo ("Reading RT data")
          do j=1,2
             READ (LUIN) tagi
             write(message,*) "Skipping ", tagi
             call writeInfo(message,TRIVIAL)
             READ (LUIN)
          end do

          READ (LUIN) tagi
          write(message,*) "Reading ", tagi
          call writeInfo(message,TRIVIAL)
          READ (LUIN) (uoverTarray(i),i=iiigas+1,iiigas+blocknpart)

          do j = 1, numsrt(6)-3
             READ (LUIN) tagi
             write(message,*) "Skipping ", tagi
             call writeInfo(message,TRIVIAL)
             READ (LUIN)
          end do

          do j = 1, numsrt(7)
             READ (LUIN) tagi
             write(message,*) "Skipping ", tagi
             call writeInfo(message,TRIVIAL)
             READ (LUIN)
          end do

          call writeInfo ("Read RT data")
       endif

       if (nblocktypes.GE.4) then
          call writeInfo ("Reading MHD data")
          do i=1,8
             do j=1,nums(i)
                READ (LUIN) ! Skip tags
                READ (LUIN)
             end do
          end do
          call writeInfo ("Read MHD data")
       endif

       iiigas = iiigas+blocknpart

    end do mpi_blocks

    close(LUIN)

    if (blocksum_npart == npart ) then 
       write(*,*) "Read ", npart, " particles"
!       call writeInfo(message, TRIVIAL)
    else
       write(message,*) "Read ",blocksum_npart, "particles but expected ", npart
       call writeWarning(message)
    endif

! If the grid is split over MPI processes then only the rank zero process will get here. 
! It then needs to loop over all the hydro threads to hand out the particles. 
! Otherwise all MPI process get here and execute this loop once. 
    if (splitOverMPI) then
       loopIndex = nHydroThreadsGlobal
    else
       loopIndex = 1
    endif

#ifdef MPI
#ifdef HYDRO
hydroThreads: do iThread = 1, loopIndex
#endif
#endif


!
! Part 2: Set up SPH data structure

    centre = VECTOR(amrgridcentrex, amrgridcentrey, amrgridcentrez)
    halfSize = amrgridSize / 2.d0
    if (splitOverMpi) then
       call  domainCentreAndSize(iThread, centre, halfsize)
    endif
    
!

! Work out how many particles we actually need
    irequired=0
    irejected=0
    i_pt_mass =0 
    do i=1, blocksum_npart

! Rotate/tilt for generating Galactic plane surveys
       if (internalView) call rotate_particles_mpi(galaxyPositionAngle, galaxyInclination)

       if (iphase(i) > 0) then
          if ( particleInBox(xyzmh(:,i),uDist, centre, halfSize) .and. (.not.discardSinks)) then 
             i_pt_mass = i_pt_mass + 1
          endif
       endif


       if (iphase(i) == 0) then
          if ( particleRequired(xyzmh(:,i),uDist, centre, halfSize)) then 
             irequired = irequired + 1
          else
             irejected = irejected + 1
          endif
       end if
    end do

! irequired is the number of particles we are going to store in sphData
! int is explicit convertion from kind=8 to default integer to avoid compiler warning
    npart = int(irequired)

    write(*,*) myrankGlobal, "nptmass ", i_pt_mass, npart
    call init_sph_data(udist, umass, utime, time, i_pt_mass)
    sphdata%codeVelocitytoTORUS = (udist / utime) / cspeed

! Allocate chemistry related components
    if (sphwithChem) then
       allocate ( sphData%rhoCO(npart) )
    endif
    if (sphwithChem.or.convertRhoToHI) then
       allocate ( sphData%rhoH2(npart) )
    endif

    if (nblocktypes.GE.3) then
       sphdata%codeEnergytoTemperature = 1.0
    else if (convertRhoToHI) then
       write(message,*) "Will convert density to HI density and set molecular weight accordingly"
       call writeInfo(message,FORINFO)
       sphdata%codeEnergytoTemperature = 1.0
    else if (sphWithChem) then
       write(message,'(a)') "Will set molecular weight using H2 fraction but not convert total density to HI density"
       call writeInfo(message,FORINFO)
       sphdata%codeEnergytoTemperature = 1.0
    else
       sphdata%codeEnergytoTemperature = (2.0/3.0) * (gmw / Rgas) * ( (udist**2)/(utime**2) )
       write(message,*) "Conversion factor between u and temperature (assumes molecular weight of", gmw,"): ", &
            sphdata%codeEnergytoTemperature
       call writeInfo(message,FORINFO)
    endif

    iiigas=0
    iiisink = 0 
    sphdata%totalgasmass = 0.0
    sphdata%totalHImass  = 0.0
    sphdata%totalMolmass = 0.0

! This loop needs to be over all the particles read in
    do i=1, blocksum_npart
       if (iphase(i) == 0) then
          if ( particleRequired(xyzmh(:,i),uDist, centre, halfSize)) then 

             iiigas=iiigas+1

             sphData%rhon(iiigas)  = rho(i)

! Calculate chemistry related components
             if (convertRhoToHI.or.sphWithChem) then 
                sphData%rhoH2(iiigas) = h2ratio(i)*2.0*rho(i)*5./7.
             endif
             if (sphWithChem) then
                sphData%rhoCO(iiigas) = rho(i) * COfrac(i)
             endif

             sphdata%vxn(iiigas)         = vxyzu(1,i)
             sphdata%vyn(iiigas)         = vxyzu(2,i)
             sphdata%vzn(iiigas)         = vxyzu(3,i)
             sphData%temperature(iiigas) = vxyzu(4,i)

! If the radiative transfer block exists, set temperatures as u / (u/T)
             if (nblocktypes.GE.3) then
                if (uoverTarray(i).le.0.0) then
                   call writeFatal("u over T is <=0")
                   write(*,*) "u over t ",i,uOverTarray(i)
                   stop
                endif
                sphData%temperature(iiigas) = sphData%temperature(iiigas) / uoverTarray(i)
             endif

! Convert internal energy to temperature using the mean molecular weight derived from the H2 fraction
! sphdata%codeEnergytoTemperature = 1.0 so this needs to be temperature
             if (convertRhoToHI.or.sphWithChem) then
                gmwWithH2 = (2.*h2ratio(i)+(1.-2.*h2ratio(i))+0.4) / (0.1+h2ratio(i)+(1.-2.*h2ratio(i)))
                sphdata%temperature(iiigas) = (2.0/3.0) * sphData%temperature(iiigas) * ( gmwWithH2 / Rgas) &
                     * ( (udist**2)/(utime**2) )
             endif

             sphData%xn(iiigas)          = xyzmh(1,i)
             sphData%yn(iiigas)          = xyzmh(2,i)
             sphData%zn(iiigas)          = xyzmh(3,i)
             sphData%gasmass(iiigas)     = xyzmh(4,i)
             sphData%hn(iiigas)          = xyzmh(5,i)

             sphdata%totalgasmass = sphdata%totalgasmass + sphData%gasmass(iiigas)
             if (convertRhoToHI.or.sphwithChem) then
                sphdata%totalHImass  = sphdata%totalHImass + (1.0-2.0*h2ratio(i))*sphData%gasmass(iiigas)*5.0/7.0
             endif
             if (sphwithChem) then
                sphdata%totalMolmass = sphdata%totalMolmass + COfrac(i)*sphData%gasmass(iiigas)
             endif

          endif
       else if(iphase(i) > 0) then
          if ( particleInBox(xyzmh(:,i),uDist, centre, halfSize).and. (.not.discardSinks)) then 
       
             iiisink=iiisink+1


             sphdata%x(iiisink)         = xyzmh(1,i)
             sphdata%y(iiisink)         = xyzmh(2,i)
             sphdata%z(iiisink)         = xyzmh(3,i)
             sphData%ptmass(iiisink)    = xyzmh(4,i)
             sphdata%vx(iiisink)         = vxyzu(1,i)
             sphdata%vy(iiisink)         = vxyzu(2,i)
             sphdata%vz(iiisink)         = vxyzu(3,i)
          endif
       endif
    end do
    sphdata%nptmass = int(iiisink)

    write(message,*) "Stored ", irequired, " gas particles which are required"
    call writeinfo(message, TRIVIAL)    

    write(message,*) "Rejected ", irejected, " gas particles which do not influence the grid"
    call writeinfo(message, TRIVIAL)    

    write(message,*) "Total number of gas particles is ", irequired+irejected
    call writeinfo(message, TRIVIAL)    

    write(message,*) "Total number of required sink particles is ", iiisink
    call writeinfo(message, TRIVIAL)    

    write(message,*) "Minimum x=", minval(xyzmh(1,:)) * udist
    call writeinfo(message, TRIVIAL)
    write(message,*) "Maximum x=", maxval(xyzmh(1,:)) * udist
    call writeinfo(message, TRIVIAL)

    write(message,*) "Minimum y=", minval(xyzmh(2,:)) * udist
    call writeinfo(message, TRIVIAL)
    write(message,*) "Maximum y=", maxval(xyzmh(2,:)) * udist
    call writeinfo(message, TRIVIAL)

    write(message,*) "Minimum z=", minval(xyzmh(3,:)) * udist
    call writeinfo(message, TRIVIAL)
    write(message,*) "Maximum z=", maxval(xyzmh(3,:)) * udist
    call writeinfo(message, TRIVIAL)

#ifdef MPI
#ifdef HYDRO

    if (splitOverMPI) then

    call mpi_send(uDist, 1, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(uMass, 1, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(uTime, 1, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphData%codeVelocityToTorus, 1, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)


    call mpi_send(sphdata%nPart, 1, MPI_INTEGER, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%rhon, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)

    call mpi_send(sphdata%vxn, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%vyn, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%vzn, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)

    call mpi_send(sphdata%temperature, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)

    call mpi_send(sphdata%xn, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%yn, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%zn, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)

    call mpi_send(sphdata%gasmass, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%hn, sphData%nPart, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)

    call mpi_send(sphdata%nptmass, 1, MPI_INTEGER, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%ptmass, sphData%nptmass, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)

    call mpi_send(sphdata%x, sphData%nptmass, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%y, sphData%nptmass, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%z, sphData%nptmass, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)

    call mpi_send(sphdata%vx, sphData%nptmass, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%vy, sphData%nptmass, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)
    call mpi_send(sphdata%vz, sphData%nptmass, MPI_DOUBLE_PRECISION, ithread, tag, localWorldCommunicator,ierr)


    end if

 enddo hydroThreads

    deallocate(rho)
    deallocate(vxyzu)
    deallocate(xyzmh)
    deallocate(uoverTarray)
    deallocate(iphase)
    deallocate(isteps)
    deallocate(COfrac,h2ratio)

    else

       call mpi_recv(sphData%uDist, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       call mpi_recv(sphData%uMass, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       call mpi_recv(sphData%uTime, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       call mpi_recv(sphData%codeVelocitytoTorus, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

       call mpi_recv(sphData%nPart, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
       npart = sphData%nPart
       allocate(sphdata%rhon(1:sphData%nPart))
       call mpi_recv(sphData%rhon, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

       allocate(sphdata%vxn(1:sphData%nPart))
       allocate(sphdata%vyn(1:sphData%nPart))
       allocate(sphdata%vzn(1:sphData%nPart))
       call mpi_recv(sphData%vxn, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       call mpi_recv(sphData%vyn, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       call mpi_recv(sphData%vzn, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

       allocate(sphdata%temperature(1:sphData%nPart))
       call mpi_recv(sphData%temperature, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

       allocate(sphdata%xn(1:sphData%nPart))
       allocate(sphdata%yn(1:sphData%nPart))
       allocate(sphdata%zn(1:sphData%nPart))
       call mpi_recv(sphData%xn, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       call mpi_recv(sphData%yn, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       call mpi_recv(sphData%zn, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)


       allocate(sphdata%gasmass(1:sphData%nPart))
       call mpi_recv(sphData%gasmass, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

       allocate(sphdata%hn(1:sphData%nPart))
       call mpi_recv(sphData%hn, sphData%nPart, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       write(*,*) myrankglobal, " received ", sphData%npart

       call mpi_recv(sphData%nptmass, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)

       allocate(sphData%ptmass(1:sphData%nptmass))
       call mpi_recv(sphData%ptmass, sphData%nptmass, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

       allocate(sphData%x(1:sphData%nptmass))
       call mpi_recv(sphData%x, sphData%nptmass, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       allocate(sphData%y(1:sphData%nptmass))
       call mpi_recv(sphData%y, sphData%nptmass, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       allocate(sphData%z(1:sphData%nptmass))
       call mpi_recv(sphData%z, sphData%nptmass, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

       allocate(sphData%vx(1:sphData%nptmass))
       call mpi_recv(sphData%vx, sphData%nptmass, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       allocate(sphData%vy(1:sphData%nptmass))
       call mpi_recv(sphData%vy, sphData%nptmass, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
       allocate(sphData%vz(1:sphData%nptmass))
       call mpi_recv(sphData%vz, sphData%nptmass, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)


    endif
    call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
#endif

contains

   subroutine rotate_particles_mpi(thisPositionAngle, thisInclination)
     real(double), intent(in) ::  thisPositionAngle, thisInclination
     TYPE(VECTOR) :: orig_sph, rot_sph

     orig_sph%x = xyzmh(1,i)
     orig_sph%y = xyzmh(2,i)
     orig_sph%z = xyzmh(3,i)
     rot_sph  = rotateZ( orig_sph,  thisPositionAngle*degToRad )
     rot_sph  = rotateY( rot_sph, thisInclination*degToRad )
     xyzmh(1,i) = rot_sph%x
     xyzmh(2,i) = rot_sph%y
     xyzmh(3,i) = rot_sph%z

     orig_sph%x = vxyzu(1,i)
     orig_sph%y = vxyzu(2,i)
     orig_sph%z = vxyzu(3,i)
     rot_sph  = rotateZ( orig_sph,  thisPositionAngle*degToRad )
     rot_sph  = rotateY( rot_sph, thisInclination*degToRad )
     vxyzu(1,i) = rot_sph%x
     vxyzu(2,i) = rot_sph%y
     vxyzu(3,i) = rot_sph%z

   end subroutine rotate_particles_mpi


  end subroutine read_sph_data_mpi

  logical function particleRequired(positionArray, uDist, cen, halfsize)
    use inputs_mod, only: sphboxcut, sphboxxmin, sphboxxmax, sphboxymin, sphboxymax, &
         sphboxzmin, sphboxzmax, sphspherecut, sphspherex, sphspherey, sphspherez, sphsphereradius
    type(VECTOR) :: cen
    real(double) :: halfSize
    real(kind=8) :: positionArray(5)
    real(db) :: x, y, z, h
    real(db) :: distFromGridCentre, maxDist, distFromSphereCentre
    real(db), intent(in) :: uDist

    ! Particle position in Torus units
    x = positionArray(1) * uDist / 1.0e10
    y = positionArray(2) * uDist / 1.0e10
    z = positionArray(3) * uDist / 1.0e10
    h = positionArray(5) * uDist / 1.0e10

! If we are limiting the particles to a user specified box and this particle
! is outside the box then it is not required. Otherwise continue to see whether
! the particle influences the grid.
    if (sphboxcut) then 
       if ( x<sphboxxmin .or. x>sphboxxmax .or. &
            y<sphboxymin .or. y>sphboxymax .or. &
            z<sphboxzmin .or. z>sphboxzmax) then
          particleRequired = .false.
          return
       endif
    endif

! Ditto for limiting particles to a sphere
    if (sphspherecut) then
       distFromSphereCentre = sqrt ( (x-sphspherex)**2 + (y-sphspherey)**2 + (z-sphspherez)**2 )
       if (distFromSphereCentre > sphsphereradius) then
          particleRequired = .false.
          return
       endif
    endif

    distFromGridCentre = sqrt ( (x-cen%x)**2 + (y-cen%y)**2 + (z-cen%z)**2)

! Maximum distance at which a particle can influence the grid. This is taken to be the distance from the 
! grid centre to a grid corner (in 3D) plus 3 particle smoothing lengths
    maxDist = sqrt(3.0)*halfSize + 3.0 * h

    if (distFromGridCentre < maxDist) then 
       particleRequired = .true.
    else
       particleRequired = .false.
    endif

  end function particleRequired

  logical function particleInBox(positionArray, uDist, cen, halfsize)
    type(VECTOR) :: cen
    real(double) :: halfSize
    real(kind=8) :: positionArray(5)
    real(db) :: x, y, z
    real(db), intent(in) :: uDist

    ! Particle position in Torus units
    x = positionArray(1) * uDist / 1.0e10
    y = positionArray(2) * uDist / 1.0e10
    z = positionArray(3) * uDist / 1.0e10

    particleInBox = .true.
    if (x > (cen%x+halfSize)) then
       particleInBox = .false.
    else if (y > (cen%y+halfSize)) then
       particleInBox = .false.
    else if (z > (cen%z+halfSize)) then
       particleInBox = .false.
    else if (x < (cen%x-halfSize)) then
       particleInBox = .false.
    else if (y < (cen%y-halfSize)) then
       particleInBox = .false.
    else if (z < (cen%z-halfSize)) then
       particleInBox = .false.
    endif

  end function particleInBox

  !
  !
  ! accessors
  !
  
  ! returns program units of distance in cm 
  Function get_udist() RESULT(out)
    implicit none
    real(double) :: out
    out = sphdata%udist
  end function get_udist

  ! returns program units of mass in g
  function get_umass() RESULT(out)
    implicit none
    real(double) :: out
    out = sphdata%umass
  end function get_umass

  ! returns program units of time in s
  function get_utime() RESULT(out)
    implicit none
    real(double) :: out
    out = sphdata%utime
  end function get_utime
    
  ! returns program units of velocity in cm/s
  function get_uvel() RESULT(out)
    implicit none
    real(double) :: out
    out = sphdata%uvel
  end function get_uvel

  ! returns the number of gas particles
  function get_npart() RESULT(out)
    implicit none
    integer ::out 
    out = sphdata%npart
  end function get_npart
    

  ! returns the number of point masses
  function get_nptmass() RESULT(out)
    implicit none
    integer :: out
    out = sphdata%nptmass
  end function get_nptmass
    

  ! returns the time of dump time in the unit of [utime]
 function get_time() RESULT(out)
    implicit none
    real(double) :: out
    out = sphdata%time
  end function get_time
    

  !
  ! Returns the position of the i_th gas particle. 
  !
  ! x, y, z are in [udist] ... See the type definition section.
  subroutine get_position_gas_particle(i, x, y, z)
    implicit none    
    integer, intent(in) :: i
    real(double), intent(out) :: x, y, z
    
    X = sphdata%xn(i)
    y = sphdata%yn(i)
    z = sphdata%zn(i)
    
  end subroutine get_position_gas_particle

  !
  ! saves the position of the i_th gas particle in this object. 
  !
  ! "name" must be one of the following: "x", "y", "z"
  ! x, y, z are in [udist] ... See the type definition section.
  subroutine put_position_gas_particle(i, name, value)
    implicit none    
    integer, intent(in) :: i
    character(LEN=*), intent(in) :: name
    real(double), intent(in) :: value
    
    select case(name)
    case ("x", "X")
       sphdata%xn(i) = value
    case ("y", "Y") 
       sphdata%yn(i) = value
    case ("z", "Z") 
       sphdata%zn(i) = value
    case default
       write(*,*) "Error: Unknown name passed to sph_data_class::put_position_gas_particle."
       stop
    end select
    
  end subroutine put_position_gas_particle


  ! Returns the density of gas particle at the postion of
  ! i-th particle.
  
  function get_rhon(i) RESULT(out)
    implicit none
    real(double) :: out 
    integer, intent(in) :: i

    out  = sphdata%rhon(i)
    
  end function get_rhon

!!$  ! Returns the temperature of gas particle at the postion of
!!$  ! i-th particle.
!!$  
!!$  function get_temp(i) RESULT(out)
!!$    implicit none
!!$    real(double) :: out 
!!$    Integer, intent(in) :: i
!!$
!!$    out  = sphdata%temperature(i)
!!$    
!!$  end function get_temp

  function sphVelocityPresent () RESULT(out)
    implicit none
    logical :: out

    if (associated(sphdata%vxn) .and. associated(sphdata%vyn) .and. associated(sphdata%vzn) ) then 
       out = .true.
    else
       out=.false.
    end if

  end function sphVelocityPresent

  function get_vel(i) RESULT(out)
    implicit none
    type(VECTOR) :: out 
    integer, intent(in) :: i

    out  = VECTOR(sphdata%vxn(i),sphdata%vyn(i),sphdata%vzn(i))
    
  end function get_vel

  ! Assigns the density of gas particle at the postion of
  ! i-th particle.
  
  subroutine put_rhon(i, value)
    implicit none
    integer, intent(in) :: i
    real(double), intent(in) :: value

    sphdata%rhon(i) = value
    
  end subroutine put_rhon
     
  !
  ! Rerurns the postions of the i-th point mass.
  !
  ! x, y, z are in [udist] ... See the type definition section.
  subroutine get_position_pt_mass(i, x, y, z)
    implicit none    
    integer, intent(in) :: i
    real(double), intent(out) :: x, y, z
    
    x = sphdata%x(i) 
    y = sphdata%y(i)
    z = sphdata%z(i)
    
  end subroutine get_position_pt_mass
 

  !
  ! Returns the mass of the point mass (star) in
  ! the i-th position of data.
  !
  function get_pt_mass(i) RESULT(out)
    implicit none
    real(double)  :: out 
    integer, intent(in):: i

    out = sphdata%ptmass(i)  !in program units [umass]. See above.

  end function get_pt_mass

  !
  ! Returns the position of the point mass (star) in
  ! the i-th position of data.
  !
  function get_pt_position(i) RESULT(out)
    implicit none
    type(vector)  :: out 
    integer, intent(in):: i

    out = VECTOR(sphdata%x(i),sphdata%y(i),sphdata%z(i))  !in program units [umass]. See above.

  end function get_pt_position

  !
  ! Returns the velocity of the point mass (star) in
  ! the i-th position of data.
  !
  function get_pt_velocity(i) RESULT(out)
    implicit none
    type(vector)  :: out 
    integer, intent(in):: i

    out = VECTOR(sphdata%vx(i),sphdata%vy(i),sphdata%vz(i))  !in program units [umass]. See above.

  end function get_pt_velocity
  
  
  

  !
  !  Destructor
  ! 
  !  Deallocates the array memories

  subroutine kill_sph_data()
    implicit none
    
    DEALLOCATE(sphdata%xn, sphdata%yn, sphdata%zn)
    DEALLOCATE(sphdata%rhon, sphdata%temperature)
    DEALLOCATE(sphdata%x, sphdata%y, sphdata%z)
    DEALLOCATE(sphdata%ptmass)

    NULLIFY(sphdata%xn, sphdata%yn, sphdata%zn)
    NULLIFY(sphdata%rhon, sphdata%temperature)
    NULLIFY(sphdata%x, sphdata%y, sphdata%z)
    NULLIFY(sphdata%ptmass)

    if ( ASSOCIATED(sphdata%gasmass) ) then 
       DEALLOCATE(sphdata%gasmass)
       NULLIFY(sphdata%gasmass)
    end if

    if ( ASSOCIATED(sphdata%hn) ) then 
       DEALLOCATE(sphdata%hn)
       NULLIFY(sphdata%hn)
    end if

    if ( ASSOCIATED(sphdata%vxn) ) then 
       DEALLOCATE(sphdata%vxn)
       NULLIFY(sphdata%vxn)
    end if

    if ( ASSOCIATED(sphdata%vyn) ) then 
       DEALLOCATE(sphdata%vyn)
       NULLIFY(sphdata%vyn)
    end if

    if ( ASSOCIATED(sphdata%vzn) ) then 
       DEALLOCATE(sphdata%vzn)
       NULLIFY(sphdata%vzn)
    end if

    if ( ASSOCIATED(sphdata%hpt) ) then 
       DEALLOCATE(sphdata%hpt)
       NULLIFY(sphdata%hpt)
    end if

    if ( ASSOCIATED(sphdata%vx) ) then 
       DEALLOCATE(sphdata%vx)
       NULLIFY(sphdata%vx)
    end if

    if ( ASSOCIATED(sphdata%vy) ) then 
       DEALLOCATE(sphdata%vy)
       NULLIFY(sphdata%vy)
    end if

    if ( ASSOCIATED(sphdata%vz) ) then 
       DEALLOCATE(sphdata%vz)
       NULLIFY(sphdata%vz)
    end if

    if ( ASSOCIATED(sphdata%rhoH2) ) then 
       DEALLOCATE(sphdata%rhoH2)
       NULLIFY(sphdata%rhoH2)
    end if

    if ( associated (sphData%rhoCO) ) then 
       deallocate(sphData%rhoCO)
       nullify(sphData%rhoCO)
    end if

    sphdata%inUse = .false.

  end subroutine kill_sph_data


  subroutine deallocate_sph

    if (.not.sphdata%inUse) return

    call writeInfo("Deallocating SPH data", TRIVIAL)
    call kill_sph_data

    call writeInfo("Deallocating clusterParameter storage", TRIVIAL)
    call deallocate_clusterparameter_arrays

  end subroutine deallocate_sph

  !
  ! Retuns the minimum value of rhon
  !
  function get_rhon_min() RESULT(out)
    implicit none
    real(double) :: out

    out = MINVAL(sphdata%rhon)
    
  end function get_rhon_min

  
  !
  ! retuns the maximum value of rhon
  !
  function get_rhon_max() RESULT(out)
    implicit none
    real(double) :: out

    out = MAXVAL(sphdata%rhon)
    
  end function get_rhon_max



  !
  ! retuns the spin direction of i_th star
  !
  subroutine get_spins(i, sx, sy, sz) 
    implicit none
    integer, intent(in) :: i 
    real(double), intent(out) :: sx, sy, sz 
    
    sx = sphdata%spinx(i)
    sy = sphdata%spiny(i)
    sz = sphdata%spinz(i)
    
  end subroutine get_spins



  !
  ! Prints basic infomation
  !
  ! if filename is '*' then it prints on screen.
  subroutine info_sph(filename)
    use inputs_mod, only: convertRhoToHI, sphwithchem
    implicit none
    character(LEN=*), intent(in) :: filename
    integer :: UN
    real(double) :: tmp
    real(double) :: myTotalMass, myHImass, myMolmass

    if (filename(1:1) == '*') then
       UN = 6   ! prints on screen
    else
       UN = 69
       open(unit=UN, file = TRIM(filename), status = 'replace')
    end if

    ! time of the data dump
    tmp = get_time()*get_utime()/(60.0d0*60.0d0*24.0d0*365.0d0*1.0d6)
    
    ! Total gas mass, HI mass and molecular mass in solar masses
    myTotalMass = sphdata%totalGasMass * get_umass() / mSol
    myHImass    = sphdata%totalHIMass  * get_umass() / mSol
    myMolMass   = sphdata%totalMolMass * get_umass() / mSol

    Write(UN,'(a)') ' '
    write(UN,'(a)') '################################################################'
    write(UN,'(a)') 'SPH data info :'
    write(UN,'(a)') ' '    
    write(UN,*)     'Units of length            : ', get_udist(), ' [cm]'
    write(UN,*)     'Units of mass              : ', get_umass(), ' [g]'
    write(UN,*)     'Units of time              : ', get_utime(),  ' [s]' 
    write(UN,*)     'Units of velocity          : ', get_uvel(),   '[cm/s]'
    write(UN,'(a)') ' '    
    write(UN,*)     '# of stars                 : ',  get_nptmass()
    write(UN,*)     '# of gas particles (total) : ',  get_npart()   
    write(UN,*)     'Time of data dump          : ',  tmp, ' [Myr]'    
    write(UN,*)     'Total gas mass             : ',  myTotalMass, ' [M_Sun]'
    if (convertRhoToHI.or.sphwithchem) &
          write(UN,*)     'Total HI mass              : ',  myHIMass, ' [M_Sun]'
    if (sphwithchem) &
          write(UN,*)     'Total molecular mass       : ',  myMolMass, ' [M_Sun]'
    write(UN,'(a)') '############################################################'
    write(UN,'(a)') ' '
    write(Un,'(a)') ' '
    
    if (filename(1:1) /= '*')  close(UN)

  end subroutine info_sph



  !
  ! Returns the parameters for stellar disc of i-th star
  ! in the data.
  ! x, y, z are in [udist] ... See the type definition section.
  subroutine get_stellar_disc_parameters(i, discrad, discmass, &
       spinx, spiny, spinz)
    implicit none    
    integer, intent(in) :: i
    real(double), intent(out) :: discrad    ! in [10^10cm] 
    real(double), intent(out) :: discmass   ! in [grams] 
    real(double), intent(out) :: spinx, spiny, spinz
    logical, save :: first_time = .true.
    
    ! quick check
    if (first_time) then
       if (sphdata%have_stellar_disc) then
          first_time = .false.
       else
          write(*,*) "Error:: You did not read in the stellar disc data! &
               &[get_stellar_disc_parameters]."
          stop
       end if
    end if
    
    discrad = sphdata%discrad(i)    ! in [10^10cm] 
    discmass = sphdata%discmass(i)  ! in [grams] 
    spinx = sphdata%spinx(i)
    spiny = sphdata%spiny(i)
    spinz = sphdata%spinz(i)
        
  end subroutine get_stellar_disc_parameters


  !
  !
  !
  function stellar_disc_exists() RESULT(out)
    implicit none
    logical :: out
    out = sphdata%have_stellar_disc
  end function stellar_disc_exists




  !
  ! Compute the inclination angles  angle between the spin ]
  ! axis and the observer. The results will be written in 
  ! a file. 
  !
  subroutine find_inclinations(obs_x, obs_y, obs_z, outfilename)
    implicit none
    real(double), intent(in) :: obs_x,  obs_y, obs_z  ! directions cosines of of observer.
    character(LEN=*), intent(in) :: outfilename 
    real(double) :: r1, r2, inc, dp, pi
    real(double) :: sx, sy, sz ! spins
    integer :: i

    open(unit=43, file=TRIM(ADJUSTL(outfilename)), status="replace")

    pi = 2.0d0*acos(0.0d0)

    write(43, '(a)') "#   star ID    ------   inclination [deg]"

    ! just in case the vectors are not normalized.
    r1 = SQRT(obs_x*obs_x + obs_y*obs_y + obs_z*obs_z)
    do i = 1, sphdata%nptmass       
       call get_spins(i , sx, sy, sz)
       ! just in case the vectors are not normalized.
       r2 = SQRT(sx*sx+sy*sy+sz*sz)
    
       ! inclinations
       dp = obs_x*sx + obs_y*sy + obs_z*sz  ! dot product

       inc = ACOS(dp/r1/r2)*(180.d0/Pi)   ! degrees
       if (inc>180.d0) inc = inc-180.d0
       if (inc>90.d0) inc = 90.0 - (inc-90.d0)
       write(43, '(1x, i5, 2x, f6.1)') i, inc
       
    end do
    

    close(43)

  end subroutine find_inclinations
 

  !
  ! Interface function to get the status of SPH data object
  !
  pure function isAlive() RESULT(out)
    implicit none
    logical :: out
    out = sphdata%inUse
  end function isAlive

  subroutine sortbyx(xarray,ind)
    use utils_mod, only: sortdouble2index

    real(double) :: xarray(:)
    integer :: ind(:)
    integer :: i
      
    do i=1,npart
       ind(i) = i
    Enddo
    
    call sortdouble2index(xarray,ind)

  end subroutine sortbyx

  subroutine FindCriticalValue(array,critval,percentile,output)
!    use timing
    use utils_mod, only: dquicksort
    real(double) :: array(:), critval
    real(double), intent(in) :: percentile
    integer :: i
    logical, optional :: output 

    character(len=100) :: message

!    npart = size(array)

    if(output .and. writeoutput) then
    
!       call writeinfo("Start sort", TRIVIAL)
!       call tune(6, "Sort")  ! start a stopwatch
       call dquicksort(array)
!       call tune(6, "Sort")  ! stop a stopwatch
!       call writeinfo("Start sort", TRIVIAL)

!       open(unit=50, file='sortedarray.dat', status="replace")
       open(unit=49, file='percentiles.dat', status="unknown")
    
!       do i=1,npart
!          write(50,*) array(i)
!       enddo

       do i=100,1,-1
! Using an explcit conversion to real prevents calculation of 
! npart * i which can overflow for large particle numbers
          write(49,*) i,"%",array(nint (npart*real(i/100.)) )
       enddo
    
       write(message,*) "Maximum Value",array(npart)
       call writeinfo(message, FORINFO)
       write(message,*) "Minimum Value",array(1)
       call writeinfo(message, FORINFO)
    else
       call dquicksort(array)
    endif

    critval = array(max(1,nint(percentile*npart)))

  end subroutine FindCriticalValue

  TYPE(vector)  function Clusterparameter(point, thisoctal, subcell, rho_out, rhoH2_out, rhoCO_out, temp_out, dustfrac_out)
    USE inputs_mod, only: sph_norm_limit, convertRhoToHI, sphToGridSimple
    USE constants_mod, only: tcbr
    use octal_mod, only: OCTAL

    type(vector), intent(in) :: point
    type(vector),save :: oldpoint = VECTOR(0.d0,0.d0,0.d0)
    type(vector),save :: oldvel = VECTOR(2.d0,2.d0,2.d0)
#ifdef SPH2GRID_PARALLEL
!$OMP threadprivate(oldpoint, oldvel)
#endif
    real(double), optional, intent(out) :: rho_out
    real(double), optional, intent(out) :: rhoH2_out
    real(double), optional, intent(out) :: rhoCO_out
    real(double), optional, intent(out) :: temp_out
    real(double), optional, intent(out) :: dustfrac_out
    real(double) :: rhoH2_local

    real(double) :: r
    real(double) :: d
    real(double) :: fac
    real(double) :: sumWeight
    real(double) :: vx, vy, vz
    real(double) :: h2ratio

    integer :: i

    integer :: nparticles
    integer :: indexArray(npart)
    real(double) :: partArray(npart), q2array(npart)

    type(octal), pointer :: thisOctal
    integer :: subcell

! Initialise optional parameters to missing data values
    if (present(rho_out))      rho_out=1d-30
    if (present(rhoH2_out))    rhoH2_out=1d-37
    if (present(rhoCO_out))    rhoCO_out=1d-99
    if (present(temp_out))     temp_out=tcbr
    if (present(dustfrac_out)) dustfrac_out=1d-99
    rhoH2_local=1d-30

! If there are no particles then there is nothing to do so just return. 
! Return values for optional parameters have been set to missing values above
    if (get_npart() == 0) then
          Clusterparameter = VECTOR(0.d0,0.d0,0.d0)
       return
    endif
              
    d = min(thisoctal%h(subcell), rmax)! using the placeholder h from splitgrid

!   d = thisoctal%subcellsize ! the splitgrid routine effectively picks the grid size based on smoothing length (mass condition)

    ! If point is outside particles then there will be no information so just return.
    ! Hope the compiler optimizes this (for a disc, first condition is simpler to test)
    if(abs(point%z) .gt. maxz + rmax .or. &
         (point .dot. point) .gt. maxr2) then 
       if(modulus(point - oldpoint) .lt. 2. * d) then! cells nearby?
          Clusterparameter = oldvel ! Velocity of nearby cells
       else
 ! 2c is not a real velocity and everything should be trying to handle this error
          Clusterparameter = VECTOR(2.d0,2.d0,2.d0)
       endif
       return
    endif

    if (npart == 0) then
       write(*,*) "Error npart is zero"
       write(*,*) "rank ",myrankGlobal
       write(*,*) "sph ",sphdata%npart
    endif

    if (sphToGridSimple) then

       ! A simple algorithm for particle selection parallelised with OpenMP
       call findParticles_simple(point, nparticles, indexArray, partarray, sumweight)

    else

    ! --------------------------------------------------------------------------------------------
    ! Build a list of particles to use for calculating kernel smoothed values. We will make 
    ! three attempts to find a list of particles with a sum of smoothing weights above 1.0e-3. 
    ! This method trades off time for accuracy in low density regions. See Dave Rundle's thesis
    ! for full details. 
    ! --------------------------------------------------------------------------------------------

    ! Attempt 1: 4d is far enough away to have particles with their smoothing lengths captured, 
    ! rcrit puts an upper limit on time
       nparticles = 0; sumweight = 0.d0
       r = min(4.d0 * d, rcrit) 
       call findNearestParticles(point, nparticles, r, indexArray, q2array)
       call doweights(sumweight, nparticles, q2array, partarray, indexArray)

    ! Attempt 2: 
    ! 2 rcrit is essential for the mass to be correctly done... (in my case it had to be hcrit = 99%)
       if (sumweight .le. 1d-3) then
          r = min(max(r * 4.d0, 2.d0 * rcrit), rmax * 0.1)
          call findNearestParticles(point, nparticles, r, indexArray, q2array)
          call doweights(sumweight, nparticles, q2array, partarray, indexArray)
       end if

    ! Attempt 3: Search out to 0.5*rmax. If the sum of weights is still too small then this is 
    ! considered an emnpty cell
       if (sumweight .le. 1d-3) then
          r = rmax * 0.5d0
          call findNearestParticles(point, nparticles, r, indexArray, q2array)
          call doweights(sumweight, nparticles, q2array, partarray, indexArray)
       end if

    endif
    

    ! --------------------------------------------------------------------------------------------
    ! Calculate kernel smoothed values for velocity and any other required parameters 
    ! --------------------------------------------------------------------------------------------

    if(sumweight .gt. 0.d0) then

! Velocity: always normalise by sum of weights
       fac = 1.d0 / sumweight
       vx=0; vy=0; vz=0
       do i = 1, nparticles
          vx = vx + partArray(i) * VelocityArray(1, indexArray(i))
          vy = vy + partArray(i) * VelocityArray(2, indexArray(i))
          vz = vz + partArray(i) * VelocityArray(3, indexArray(i))
       enddo
          
       Clusterparameter = VECTOR(vx * fac, vy * fac, vz * fac) ! Velocity
       oldvel = clusterparameter ! store this point in case it's nearby an empty cell (for velocity)
       oldpoint = point

! Scalars: decide whether to normalise by sum of weights
       if(sumweight .gt. sph_norm_limit) then
          fac = 1.d0 / sumWeight
       else
          fac = 1.0
       end if

       ! Density
       if (present(rho_out)) then
          rho_out = 0.0
          do i = 1, nparticles
             rho_out = rho_out + partArray(i) * rhoArray(indexArray(i)) 
          enddo
          rho_out=rho_out * fac
       end if

       ! H2 density
       if (present(rhoH2_out).or.(convertRhoToHI.and.present(rho_out))) then
          if (allocated(rhoH2Array)) then
             rhoH2_local = 0.0
             do i = 1, nparticles
                rhoH2_local = rhoH2_local + partArray(i) * rhoH2Array(indexArray(i)) 
             enddo
             rhoH2_local=rhoH2_local * fac
          end if
       end if

       if (present(rhoH2_out)) rhoH2_out = rhoH2_local

       ! Convert density to HI density if required. Done here so that particle rho stores total density
       ! even if the grid needs rho_HI. This way the grid refinement is done using total density. 
       if (convertRhoToHI.and.present(rho_out)) then 
          h2ratio = 0.5 * (7.0/5.0) * rhoH2_local / rho_out
          rho_out = (1.0-2.0*h2ratio) * rho_out * 5.0/7.0
          rho_out = max(rho_out,1.0e-60_db)
       endif

       ! CO density
       if (present(rhoCO_out).and.allocated(rhoCOArray)) then
          rhoCO_out = 0.0
          do i = 1, nparticles
             rhoCO_out = rhoCO_out + partArray(i) * rhoCOArray(indexArray(i)) 
          enddo
          rhoCO_out=rhoCO_out * fac
       end if

       ! Temperature
       if (present(temp_out)) then
          temp_out = 0.0
          do i = 1, nparticles
             temp_out = temp_out + partArray(i) * TemArray(indexArray(i)) 
          enddo
          temp_out=temp_out * fac
       end if

       ! dust fraction
       if (present(dustfrac_out).and.allocated(dustfrac)) then
          dustfrac_out = 0.0
          do i = 1, nparticles
             dustfrac_out = dustfrac_out + partArray(i) * dustfrac(indexArray(i)) 
          enddo
          dustfrac_out=dustfrac_out * fac
       end if

    else

       if(modulus(point - oldpoint) .lt. 2. * d) then! cells nearby?
          Clusterparameter = oldvel ! Velocity of nearby cells
       else
          ! 2c is not a real velocity and everything should be trying to handle this error
          Clusterparameter = VECTOR(2.d0,2.d0,2.d0)
       endif

    endif

  end function Clusterparameter

  subroutine init_clusterparameter_arrays
    USE inputs_mod, only: hcritPercentile, hmaxPercentile, variableEta

    implicit none

    real(double) :: codeVelocitytoTORUS, codeLengthtoTORUS, codeDensitytoTORUS, udist, umass, utime
    real(double) :: hcrit, hmax
    real(double) :: Harray(npart)
    integer      :: ind(npart)
    character(len=100) :: message

    udist = get_udist()
    utime = get_utime()
    umass = get_umass()
    codeLengthtoTORUS =  udist/1.d10
    codeVelocitytoTORUS = sphdata%codeVelocitytoTORUS
!    codeDensitytoTORUS = umass / ((udist/1.d10) ** 3)
    codeDensitytoTORUS = umass / ((udist) ** 3)
    write(*,*) "code length to torus ",codeLengthToTorus
    write(*,*) "code vel to torus ",codeVelocityToTorus
    write(*,*) "code density to torus ",codeDensityToTorus
       
    allocate(PositionArray(3,npart)) ! allocate memory
    allocate(xArray(npart))
    allocate(RhoArray(npart))
    allocate(TemArray(npart))
    allocate(RhoH2Array(npart))
    allocate(OneOverHsquared(npart))
    if (associated(sphData%rhoCO)) allocate (rhoCOarray(npart))
    if (associated(sphData%dustfrac)) allocate (dustfrac(npart))
          
    PositionArray = 0.d0; hArray = 0.d0; ind = 0
    if (.not.associated(sphdata%xn)) then
       write(*,*) "sphdata not associated"
       write(*,*) "rank ",myrankGlobal
       write(*,*) "npart ",sphData%npart
       write(*,*) "ptmass ",sphData%nptmass
    endif
          
    xArray(:) = sphdata%xn(:) * codeLengthtoTORUS! fill with x's to be sorted       
    call sortbyx(xarray(:),ind(:)) ! sort the x's and recall their indices

    Positionarray(1,:) = sphdata%xn(ind(:)) * codeLengthtoTORUS ! fill with x's in sorted order
    PositionArray(2,:) = sphdata%yn(ind(:)) * codeLengthtoTORUS ! y's go with their x's
    PositionArray(3,:) = sphdata%zn(ind(:)) * codeLengthtoTORUS ! z's go with their x's

    
    ! Decide if we need to set velocities for this configuration
    if ( associated(sphData%vxn) ) then 
          
       allocate(VelocityArray(3,npart))
       VelocityArray(:,:) = 0.d0
       VelocityArray(1,:) = sphdata%vxn(ind(:)) * codeVelocitytoTORUS ! velocities
       
       if ( associated(sphData%vyn) ) then 
          VelocityArray(2,:) = sphdata%vyn(ind(:)) * codeVelocitytoTORUS 
       else
          STOP
          !             call torus_abort("Error in function Clusterparameter: vxn is associated but vyn is not")
       end if
             
       if ( associated(sphData%vzn) ) then 
          VelocityArray(3,:) = sphdata%vzn(ind(:)) * codeVelocitytoTORUS
       else
          STOP
          !             call torus_abort("Error in function Clusterparameter: vxn is associated but vzn is not")
       end if
             
       write(message, *) "Max/min Vx", maxval(VelocityArray(1,:)), minval(VelocityArray(1,:))          
       call writeinfo(message, TRIVIAL)
       write(message, *) "Max/min Vy", maxval(VelocityArray(2,:)), minval(VelocityArray(2,:))
       call writeinfo(message, TRIVIAL)
       write(message, *) "Max/min Vz", maxval(VelocityArray(3,:)), minval(VelocityArray(3,:))
       call writeinfo(message, TRIVIAL)
          
    end if
          
    maxx = max(maxval(PositionArray(1,:)), abs(minval(PositionArray(1,:))))
    maxy = max(maxval(PositionArray(2,:)), abs(minval(PositionArray(2,:))))
       
    maxz = max(maxval(PositionArray(3,:)), abs(minval(PositionArray(3,:))))
    maxr2 = maxx**2 + maxy**2 + maxz**2
       
    write(message, *) "Max/min x", maxval(PositionArray(1,:)), minval(PositionArray(1,:)) 
    call writeinfo(message, TRIVIAL)
    write(message, *) "Max/min y", maxval(PositionArray(2,:)), minval(PositionArray(2,:)) 
    call writeinfo(message, TRIVIAL)
    write(message, *) "Max/min z", maxval(PositionArray(3,:)), minval(PositionArray(3,:))
    call writeinfo(message, TRIVIAL)
    
    Harray(:) = sphdata%hn(ind(:)) ! fill h array
       
    call FindCriticalValue(harray, hcrit, real(hcritPercentile,kind=db), output = .true.) ! find hcrit as percentile of total h
    call FindCriticalValue(harray, hmax,  real(hmaxPercentile, kind=db), output = .false.) ! find hmax as percentile of total h
    Harray(:) = sphdata%hn(ind(:)) ! repeat second time to fix messed up array by sort
    RhoArray(:) = sphdata%rhon(ind(:))
       
    if (variableEta) then
       call writeInfo("Will calculate variable eta to allow for eta !=1.2")
       ! added to fix smoothing length discrepancy - h != 1.2 (rho/m)^(1/3)
       ! etaarray stores 1/(eta**3)
       allocate(etaArray(npart)) 
       etaarray(:) = sphdata%gasmass(ind(:)) / (rhoarray(:) * harray(:)**3)
       write(message, *) "Max/min eta", (1.d0/maxval(etaArray(:)))**(1.0/3.0), &
            (1.d0/minval(etaArray(:)))**(1.0/3.0)
       call writeinfo(message, TRIVIAL)
    else
       call writeInfo("Assuming eta=1.2 in smoothing length calculation")
    endif

    write(message, *) "Critical smoothing Length in code units", hcrit
    call writeinfo(message, TRIVIAL)
    write(message, *) "Maximum smoothing Length in code units", hmax
    Call writeinfo(message, TRIVIAL)
          
    RhoArray(:) = Rhoarray(:) * codeDensitytoTORUS
    TemArray(:) = sphdata%temperature(ind(:)) * sphdata%codeEnergytoTemperature
    
    write(message, *) "Max/min rho", maxval(RhoArray(:)), minval(RhoArray(:))
    call writeinfo(message, TRIVIAL)
    write(message, *) "Max/min Temp", maxval(TemArray(:)), minval(TemArray(:))
    call writeinfo(message, TRIVIAL)
    
    ! Set to H2 density from particles if available or flag as missing data otherwise. 
    if ( associated(sphData%rhoH2) ) then 
       RhoH2Array(:) = sphdata%rhoH2(ind(:)) * codeDensitytoTORUS
    else
       RhoH2Array(:) = -1.0e-30_db
    end if
          
    if (allocated (rhoCOarray) .and. associated(sphdata%rhoCO) ) then 
       rhoCOarray(:) = sphdata%rhoCO(ind(:)) * codeDensityToTorus
    end if
    if (allocated (dustfrac) .and. associated(sphdata%dustfrac) ) then 
       dustfrac(:) = sphdata%dustfrac(ind(:)) 
    end if
          
    hcrit = hcrit * codeLengthtoTORUS          
    hmax = hmax * codeLengthtoTORUS

    Harray(:) = Harray(:) * codeLengthtoTORUS  ! fill h array
    
    OneOverHsquared(:) = 1.d0 / (Harray(:)**2)

    write(message,*) "Critical smoothing Length in 10^10cm", hcrit
    call writeinfo(message, TRIVIAL)
    write(message,*) "Maximum smoothing Length in 10^10cm", hmax
    call writeinfo(message, TRIVIAL)
          
    rcrit = 2.d0 * hcrit ! edge of smoothing sphere
    rmax = 2.d0 * hmax ! edge of smoothing sphere

  end subroutine init_clusterparameter_arrays

  subroutine deallocate_clusterparameter_arrays

    if (allocated(PositionArray))   deallocate(PositionArray)
    if (allocated(xArray))          deallocate(xArray)
    if (allocated(RhoArray))        deallocate(RhoArray)
    if (allocated(Temarray))        deallocate(Temarray)
    if (allocated(RhoH2Array))      deallocate(RhoH2Array)
    if (allocated(OneOverHsquared)) deallocate(OneOverHsquared)
    if (allocated(etaarray))        deallocate(etaarray)
    if (allocated(rhoCOarray))      deallocate(rhoCOarray)
    if (allocated(dustfrac))        deallocate(dustfrac)
    if (allocated(VelocityArray))   deallocate (VelocityArray)

  end subroutine deallocate_clusterparameter_arrays

! A simpler version of the particle selection algorithm parallelised 
! with OpenMP. This was added to allow a comparison with the Rundle algorithm. 
! This alogithm is hard-wired to use a cut-off at q=3. 
! D. Acreman, May 2018
  subroutine findparticles_simple(pos, partcount, indexarray, partarray, sumweight)

    use inputs_mod, only : kerneltype, variableEta

    type(VECTOR), intent(in)  :: pos
    integer,      intent(out) :: partcount
    integer,      intent(out) :: indexarray(npart)
    real(double), intent(out) :: partarray(npart)
    real(double), intent(out) :: sumweight

    real(double) :: q2
    real(double), parameter :: qmax =3.0
    real(double), parameter :: qmax2=qmax**2
    real(double) :: thisWeight
    integer      :: i

    ! This is eta^-nu with eta=1.2 and nu=3
    real(double), parameter :: num = 0.578703703d0 ! (5/6)^3 = (1/1.2^3)

    ! This constant is: eta^-nu * 1/(pi^(nu/2))
    ! eta = 1.2 = 6/5 comes from the expression for the variable smoothing length
    ! nu = 3 is the number of spatial dimensions
    real(double), parameter :: scalar = 0.103927732d0 ! 5 / 6 * one over sqrtpicubed

    partcount = 0
    sumweight = 0.0
!$OMP PARALLEL DO default(none) private(i, q2, thisWeight) &
!$OMP shared(pos, partcount, indexArray, partArray, npart, OneOverHsquared, positionArray) &
!$OMP shared(variableEta, kerneltype, etaarray, xarray) &
!$OMP reduction (+: sumweight) 
    do i = 1, npart 

       q2 = (positionArray(1,i) - pos%x)**2 &
          + (positionArray(2,i) - pos%y)**2 &
          + (positionArray(3,i) - pos%z)**2 

       q2 = q2 * OneOverHsquared(i)

       if ( q2 < qmax2 ) then 

          if(kerneltype .eq. 0) then
             if (variableEta) then
                thisWeight = etaarray(i) * (1.0/pi**1.5) * exp(-q2)
             else
                thisWeight = scalar * exp(-q2)
             end if
          elseif( kerneltype .eq. 1) then
             if (variableEta) then
                thisWeight = etaarray(i) * SmoothingKernel3d(sqrt(q2))
             else
                thisWeight = num * SmoothingKernel3d(sqrt(q2))
             endif
          endif
          sumweight  = sumweight + thisWeight

!$OMP CRITICAL
          partcount = partcount + 1
          indexArray(partcount) = i
          partarray(partcount)  = thisWeight
!$OMP END CRITICAL

       endif

    enddo
!$OMP END PARALLEL DO

  end subroutine findparticles_simple

  subroutine findnearestparticles(pos, partcount, r, indexarray, q2array)

    use inputs_mod, only : kerneltype
    use utils_mod, only: locate_double_f90

    type(VECTOR), intent(in)  :: pos
    real(double), intent(in)  :: r
    integer,      intent(out) :: partcount
    integer,      intent(out) :: indexarray(npart)
    real(double), intent(out) :: q2array(npart)

    real(double) :: x,y,z
    integer :: i
    integer :: nupper, nlower
    integer :: closestXindex, testIndex
    real(double) :: r2test, q2test
    real(double) :: rtest
    real(double) :: ydiff, zdiff
    integer :: stepsize, sense
    logical :: test, prevtest, up
    real(double) :: fac, fac2

    if(kerneltype .eq. 1) then
       fac = 2.d0
       fac2 = 4.d0
    else
       fac = 5.d0
       fac2 = 25.d0
    endif

    x = pos%x
    y = pos%y
    z = pos%z


    if (npart < 50) then
       partcount = npart
       do i = 1, partcount
          indexArray(i) = i
       enddo
       nlower = 1
       nupper = partcount
       closestXindex = 0
       goto 1001
    endif

    call locate_double_f90(xarray, x, closestXindex) ! find the nearest particle to your point
  
    nupper = 1
    nlower = -1
    
    closestxIndex = min(max(1, closestXindex), npart)

    if(closestxIndex .lt. npart) then
  
    up = .true.
    nupper = 16
    stepsize = 16 ! changed from 1 to reflect fact that number of neighbours for *MOST* particles is 50. Could use 32?
    sense = 1
    prevtest = .true.

    do while (stepsize .ge. 1)
       if (min(npart,closestXindex+nUpper) < 1) then
          write(*,*) "rank ",myrankGlobal
          write(*,*) "npart ",npart,sphData%nPart
          write(*,*) "size(xarray) ",size(xArray)
          write(*,*) "closestXindex ",closestXindex
          write(*,*) "nUpper ",nUpper
       endif
       test = abs(xArray(min(npart,closestXindex + nupper)) - x) .le. r
 
       if(test .and. (nupper .eq. npart - closestXindex)) exit
       if(closestXindex + nupper .ge. npart) nupper = npart - closestXindex

       if(.not. (test .eqv. prevtest)) then
          sense = -sense ! forwards or backwards
          up = .false.
       endif
       
       if(up) then
          stepsize = 2 * stepsize
       else
          stepsize = stepsize / 2
       endif
    
       if(stepsize .lt. 1) then
          if (test) then
             nupper = nupper
          else
             nupper = nupper - 1 ! always has to be lower
          endif
          exit
       endif

       nupper = min(nupper + sense * stepsize, npart - closestXindex)

       prevtest = test
    enddo

    Else
       nupper = 0
    endif
! repeat for nlower
    if(closestXindex .gt. 1) then

    up = .true.
    nlower = -16
    stepsize = 16
    sense = -1
    prevtest = .true.

    do while (stepsize .ge. 1)
       test = abs(xArray(max(1,closestXindex + nlower)) - x) .le. r

       If(test .and. (nlower .eq. 1 - closestXindex)) exit
       if(closestXindex + nlower .le. 1) nlower = 1 - closestXindex

       if(.not. (test .eqv. prevtest)) then
          sense = -sense ! forwards or backwards
          up = .false.
       endif
       
       if(up) then
          stepsize = 2 * stepsize
       else
          stepsize = stepsize / 2
       endif
    
       if(stepsize .lt. 1) then
          if (test) then
             nlower = nlower
          else
             nlower = nlower + 1 ! always has to be higher
          endif
          exit
       endif

       nlower = max(nlower + sense * stepsize, 1 - closestXindex)

       prevtest = test
    enddo

    else
       nlower = 0
    endif

1001 continue

    rtest = r
    partcount = 0          

    do i = nlower, nupper ! search over all those particles we just found
       
       testIndex = closestXindex + i
             
       ydiff = abs(PositionArray(2,testindex) - y)

       if(ydiff .le. rtest) then ! if it's near in y then
          zdiff = abs(PositionArray(3,testIndex) - z)

          if(zdiff .le. rtest) then !only if it's near in z as well then work out contribution

             r2test = (xarray(testIndex) - x)**2 + & !The kernel will do the rest of the work for those outside the sphere
                  ydiff ** 2 + zdiff ** 2

!             write(*,*) "r^2 h^2 x ", r2test,1.d0/OneOverHsquared(testindex),x
             q2test = r2test * OneOverHsquared(testIndex) ! dimensionless parameter that we're interested in
             
             if(q2test .lt. fac2) then
                partcount = partcount + 1
                indexArray(partcount) = testIndex
                q2array(partcount) = q2test
             endif
          endif
       endif
    enddo

  end subroutine findnearestparticles

  Subroutine doWeights(sumweight, nparticles, q2array, partarray, indexArray)

    use inputs_mod, only : kerneltype, variableEta

    ! This is eta^-nu with eta=1.2 and nu=3
    real(double), parameter :: num = 0.578703703d0 ! (5/6)^3 = (1/1.2^3)

    ! This constant is: eta^-nu * 1/(pi^(nu/2))
    ! eta = 1.2 = 6/5 comes from the expression for the variable smoothing length
    ! nu = 3 is the number of spatial dimensions
    real(double), parameter :: scalar = 0.103927732d0 ! 5 / 6 * one over sqrtpicubed

    integer,      intent(in)  :: nparticles
    real(double), intent(in)  :: q2array(npart)
    integer,      intent(in)  :: indexArray(npart)
    real(double), intent(out) :: partarray(npart)
    real(double), intent(out) :: sumweight

    integer :: i
    real(double) :: sqrtq

    if(nparticles .gt. 0) then
              
       if(kerneltype .eq. 0) then

          if (variableEta) then
             do i = 1, nparticles
                partarray(i) = etaarray(indexArray(i)) * (1.0/pi**1.5) * exp(-q2array(i))
             end do
          else
             partarray(1:nparticles) = scalar * exp(-q2array(1:nparticles))
          end if
          
       elseif( kerneltype .eq. 1) then
          
          do i = 1, nparticles
             sqrtq = sqrt(q2array(i))
             if (variableEta) then
                partarray(i) = etaarray(indexArray(i)) * SmoothingKernel3d(sqrtq)
             else
                partarray(i) = num * SmoothingKernel3d(sqrtq)
             endif
          enddo
          
       endif
       
       Sumweight = sum(partarray(1:nparticles))
       
    else
       sumweight = 0.d0
    endif
    
  end subroutine doWeights

  real(double) pure function SmoothingKernel3d(q) result(Weight)
      
    real(double), intent(in) :: q
    real(double) :: qminus2, qminus1

    qminus2 = q - 2.d0
    qminus1 = q - 1.d0
   
    if(qminus2 .lt. 0.d0) then
       Weight = (-qminus2) ** 3
       if(qminus1 .lt. 0.d0) then
          Weight = Weight + 4.d0 * ((qminus1) ** 3)
       endif
       Weight = OneOnFourPi * Weight
    else
       Weight = 0.d0
    endif
    
  end function SmoothingKernel3d
  
 subroutine splitIntoWords(longString, word, nWord, wordLen, adjL)
   character(len=24) :: word(MaxWords)
   integer :: nWord, thisLen
   character(len=*) :: longString
   character(len=MaxAsciiLineLength) :: tempString
   logical :: stillSplitting, doAdjL
   integer, optional, intent(in) :: wordLen
   logical, optional, intent(in) :: adjL

! Set the word length. The default is 16 characters. 
   if (present(wordLen)) then
      thisLen = wordLen
   else
      thisLen = 24
      thisLen = 16
   end if

! Decide whether to left adjust the string we have been given. Default is to left adjust.
   if (present(adjL)) then
      doAdjL=adjL
   else
      doAdjL=.true.
   endif
   
   if (doAdjL) then
      tempString = ADJUSTL(longString)
   else
      tempString = longString
   endif
   
   nWord = 0
   stillSplitting = .true.
   do while(stillSplitting)
      nWord = nWord + 1
      if (nWord > MaxWords) then 
         call writeFatal("Trying to read too many words from ASCII file")
         call writeFatal("Recompile with larger MaxWords in sph_data_class")
      endif
      word(nWord) = tempString(1:thisLen)
      tempString = tempString(thisLen+1:)
      if (len(trim(tempString)) == 0) stillSplitting = .false.
   end do
 end subroutine splitIntoWords

 integer function indexWord(inputword, wordarray, nword) 
   character(len=*) :: inputWord
   character(len=24) :: wordArray(:)
   integer :: nWord, i, j
   
   indexWord = 0 
   do i = 1, nWord
      j = len(inputWord)
      if (inputWord(1:j) == adjustl(wordArray(i)(1:j))) then
         indexWord = i
         exit
      endif
   enddo

   if ( indexWord == 0 ) then 
      call writeWarning("indexWord: "//inputword//" not found")
   end if

 end function indexWord

! Check if the specified word is present in the list
 logical function wordIsPresent(inputword, wordarray, nword)
   character(len=*) :: inputWord
   character(len=24) :: wordArray(:)
   integer :: nWord, i
   
   wordIsPresent = .false.
   do i = 1, nWord
      if (inputWord(1:len(inputWord)) == wordArray(i)(1:len(inputWord))) then
         wordIsPresent = .true.
         exit
      endif
   enddo

 end function wordIsPresent

! Calculate the total mass of SPH particles found within the AMR grid. Use 
! sphdata%rhon to calculate the mass as this is used to populate  the grid 
! density in assign_grid_values. Assumes h = eta * (m/rho)^(1/3) with eta=1.2.
! D. Acreman, November 2010.
 function sph_mass_within_grid(grid)
   use gridtype_mod, only: gridtype
   use amr_utils_mod, only: inOctal

   implicit none

   type(GRIDTYPE), intent(in) :: grid
   real(db)     :: sph_mass_within_grid
   integer      :: ipart
   real(db)     :: sphDistToTorus
   real(db)     :: thisMass
   real(db), parameter :: eta=1.2
   TYPE(VECTOR) :: thisPos

   sph_mass_within_grid = 0.0
   if ( .not. sphData%inUse ) return

   sphDistToTorus = sphData%udist * 1.0e-10_db

   do ipart=1, sphdata%npart
      thisPos = VECTOR(sphdata%xn(ipart), sphdata%yn(ipart), sphdata%zn(ipart))
      thisPos = thisPos * sphDistToTorus
      if ( inOctal(grid%octreeRoot, thisPos ) ) then
         thisMass = sphdata%rhon(ipart) * (sphData%hn(ipart) / eta)**3 
         sph_mass_within_grid =  sph_mass_within_grid + thisMass
       end if
   end do

   sph_mass_within_grid = sph_mass_within_grid * sphData%umass / mSol 

 end function sph_mass_within_grid

   subroutine domainCentreAndSize(iThread, domainCentre, halfdomainSize)
    use inputs_mod, only: amrgridcentrex, amrgridcentrey, amrgridcentrez, amrgridsize
     type(VECTOR) :: centre, domainCentre
     real(double) :: halfdomainSize, d
     integer :: j
     integer :: iThread

     centre = VECTOR(amrgridcentrex, amrgridcentrey, amrgridcentrez)
     d = 0.5d0 * amrgridsize
     if (nHydroThreadsGlobal == 8) then
        halfdomainSize = d
        select case (iThread)
          CASE (1)    
             domainCentre = centre + d * VECTOR(-1.d0,-1.d0,-1.d0)
          CASE (2)    
             domainCentre = centre + d * VECTOR(1.d0,-1.d0,-1.d0)
          CASE (3)    
             domainCentre = centre + d * VECTOR(-1.d0,1.d0,-1.d0)
          CASE (4)    
             domainCentre = centre + d * VECTOR(1.d0,1.d0,-1.d0)
          CASE (5)    
             domainCentre = centre + d * VECTOR(-1.d0,-1.d0,1.d0)
          CASE (6)    
             domainCentre = centre + d * VECTOR(1.d0,-1.d0,1.d0)
          CASE (7)    
             domainCentre = centre + d * VECTOR(-1.d0,1.d0,1.d0)
          CASE (8)    
             domainCentre = centre + d * VECTOR(1.d0,1.d0,1.d0)
          end select
       endif

     if (nHydroThreadsGlobal == 64) then
        d = 0.5d0 * amrgridsize
        j =  (iThread-1)/8 + 1
        select case (j)
          CASE (1)    
             domainCentre = centre + d * VECTOR(-1.d0,-1.d0,-1.d0)
          CASE (2)    
             domainCentre = centre + d * VECTOR(1.d0,-1.d0,-1.d0)
          CASE (3)    
             domainCentre = centre + d * VECTOR(-1.d0,1.d0,-1.d0)
          CASE (4)    
             domainCentre = centre + d * VECTOR(1.d0,1.d0,-1.d0)
          CASE (5)    
             domainCentre = centre + d * VECTOR(-1.d0,-1.d0,1.d0)
          CASE (6)    
             domainCentre = centre + d * VECTOR(1.d0,-1.d0,1.d0)
          CASE (7)    
             domainCentre = centre + d * VECTOR(-1.d0,1.d0,1.d0)
          CASE (8)    
             domainCentre = centre + d * VECTOR(1.d0,1.d0,1.d0)
          end select
        j =  iThread - (j-1)*8
        d = d / 2.d0
        halfdomainSize = d
        select case (j)
          CASE (1)    
             domainCentre = domaincentre + d * VECTOR(-1.d0,-1.d0,-1.d0)
          CASE (2)    
             domainCentre = domaincentre + d * VECTOR(1.d0,-1.d0,-1.d0)
          CASE (3)    
             domainCentre = domaincentre + d * VECTOR(-1.d0,1.d0,-1.d0)
          CASE (4)    
             domainCentre = domaincentre + d * VECTOR(1.d0,1.d0,-1.d0)
          CASE (5)    
             domainCentre = domaincentre + d * VECTOR(-1.d0,-1.d0,1.d0)
          CASE (6)    
             domainCentre = domaincentre + d * VECTOR(1.d0,-1.d0,1.d0)
          CASE (7)    
             domainCentre = domaincentre + d * VECTOR(-1.d0,1.d0,1.d0)
          CASE (8)    
             domainCentre = domaincentre + d * VECTOR(1.d0,1.d0,1.d0)
          end select
       endif
     end subroutine domainCentreAndSize

     real(double) function toTorusUnits(unitName)
       character(len=*) :: unitName
       select case (unitName)
          case("[au]")
             toTorusUnits = autocm
          case("[cm]")
             toTorusUnits = 1.d-10
          case("[M_{Sun}]")
             toTorusUnits = mSol
          case("[g/cm^3]")
             toTorusUnits = 1.d0
          case("[km/s]")
             toTorusUnits = 1.d5
          case("[erg/g]")
             toTorusUnits = 1.d0
          case("[g]")
             toTorusUnits = 1.d0
          case("[cm/s]")
             toTorusUnits = 1.d0
          case DEFAULT
             write(*,*) "Unit not recognised ",trim(unitName)
             write(*,*) "Setting conversion factor to unity"
             toTorusUnits  = 1.d0
          end select
        end function toTorusUnits
     
end module sph_data_class

#endif
