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
       read_sph_data, &
       get_udist, &
       get_umass, &
       get_utime, &
       get_npart, &
       get_nptmass, &
       get_time, &
       get_position_gas_particle, &
       Put_position_gas_particle, &
       get_rhon, &
       put_rhon, &
       get_position_pt_mass, &
       get_pt_mass, &
       get_rhon_min, &
       get_rhon_max, &
       get_spins, &
       max_distance, &
       info, &
       get_stellar_disc_parameters, &
       stellar_disc_exists, &
       find_inclinations, &
       ClusterParameter, &
       isAlive

  ! At a given time (time)
  type sph_data
!     private  ! Believe me. It's better to be private!    
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
     ! Positions of gas particles
     real(double), pointer, dimension(:) :: xn,yn,zn
     real(double), pointer, dimension(:) :: vxn,vyn,vzn
     real(double), pointer, dimension(:) :: hn                 ! Smoothing length
     ! Density of the gas particles
     real(double), pointer, dimension(:) :: rhon
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
     ! Note: the units used here are diffrent from the ones used above!     
     logical :: have_stellar_disc            ! T if the following data are assigned
     real(double), pointer, dimension(:) :: discrad ! in [10^10cm]
     real(double), pointer, dimension(:) :: discmass   ! in [g]
     real(double), pointer, dimension(:) :: spinx, spiny, spinz
          
     ! Does the temperature of the SPH particles get used?
     logical :: useSphTem 
     

  end type sph_data

  real(double), allocatable :: PositionArray(:,:), OneOverHsquared(:), &
                               RhoArray(:), TemArray(:), VelocityArray(:,:), Harray(:), RhoH2Array(:)
  real(double), pointer :: rhoCOarray(:) => null()

  logical, allocatable :: HullArray(:)
  type(sph_data), save :: sphdata
  integer, save :: npart

  real(double), allocatable :: partArray(:), q2array(:), xarray(:), etaarray(:)
  Integer, allocatable :: indexArray(:)
  integer,save :: nparticles
  real(double), save :: rcrit, rmax
  type(VECTOR), save :: prevpos

  real(double) :: maxx, maxy, maxz, maxr2 ! maximum extents of particles
  
  private:: &
       kill_sph_data


  !
  !
  interface kill
     module procedure kill_sph_data
  end interface
  
  
contains  


  ! 
  ! Initializes an object with parameters (if possible).
  ! 
  subroutine init_sph_data(udist, umass, utime,  time, nptmass, uvel, utemp)
    implicit none

    real(double), intent(in)  :: udist, umass, utime    ! Units of distance, mass, time in cgs
    real(double), optional  :: uvel, utemp    ! Units of distance, mass, time in cgs
    !                                                       ! (umass is M_sol, udist=0.1 pc)
    real(double), intent(in)  :: time                   ! Time of sph data dump (in units of utime)
    integer, intent(in)           :: nptmass                ! Number of stars/brown dwarfs

    ! Indicate that this object is in use
    sphdata%inUse = .true.

! Do not use temperature from SPH particles to initialise temperature grid
    sphData%useSphTem = .false. 

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
    sphdata%nptmass = nptmass


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

! Use temperature from SPH particles to initialise temperature grid
    sphData%useSphTem = .true. 
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

! Wrapper subroutine which calls the appropriate read subroutine 
! according to the string inputFileFormat
! D. Acreman, June 2012
  subroutine read_sph_data_wrapper
    use inputs_mod, only: inputFileFormat, sphdatafilename

    select case (inputFileFormat)
    case("binarywithChem","galaxy")
       call read_sph_data_withChem(sphdatafilename)
       
    case("mpi")
       call read_sph_data_mpi(sphdatafilename)

    case("ascii")
       call new_read_sph_data(sphdatafilename)

    case default
       call writeFatal("Unrecognised file format "//inputFileFormat)

    end select

  end subroutine read_sph_data_wrapper

  !
  !
  ! Read in the data from a file, allocate the array memory, and store the
  ! the number of gas and stars in this object.
  !
  ! Note: This routine only works with the 'old' (unknown how old) dump file
  ! format. Use new_read_sph_data to read in ASCII created by SPLASH.

  subroutine read_sph_data(this, rootfilename)
    use inputs_mod, only : iModel
    use utils_mod, only : findMultiFilename
    implicit none
    type(sph_data), intent(inout) :: this
    character(LEN=*), intent(in)  :: rootfilename
    !   
    character(len=80) :: filename
    integer, parameter  :: LUIN = 10 ! logical unit # of the data file
!    real(double) :: udist, umass, utime,  time,  gaspartmass, discpartmass
!    integer*4 :: npart,  nsph, nptmass
    real(double) :: udist, umass, utime
    real(double) :: gaspartmass, time
    integer :: nptmass, n1, n2
    real(double), allocatable :: dummy(:)     


    call findMultiFilename(rootfilename, iModel, filename)
    open(unit=LUIN, file=TRIM(filename), form='unformatted')


    ! reading in the first line
    READ(LUIN) udist, umass, utime, npart, n1, n2, time, nptmass, gaspartmass

    ! initilaizing the sph_data object (allocating arrays, saving parameters and so on....)
    call init_sph_data(udist, umass, utime, time, nptmass)


    ! reading the positions  of gas particles and stars,
    ALLOCATE(dummy(npart))

    write(*,*) ' '
    write(*,*) 'Reading Matthew''s SPH data....'
    write(*,*) ' '
    READ(LUIN) this%xn
    READ(LUIN) this%yn
    READ(LUIN) this%zn

    READ(LUIN) dummy   ! Vx
    READ(LUIN) dummy   ! Vy
    READ(LUIN) dummy   ! Vz

    READ(LUIN) this%rhon
    
    READ(LUIN) this%x
    READ(LUIN) this%y
    READ(LUIN) this%z

    READ(LUIN) this%ptmass
   


    DEALLOCATE(dummy)

  end subroutine read_sph_data

! Read SPH data from a splash ASCII dump.
  subroutine new_read_sph_data(rootfilename)
    use inputs_mod, only: internalView, convertRhoToHI, ih2frac, sphwithchem
    use inputs_mod, only : iModel, galaxyPositionAngle, galaxyInclination
    use utils_mod, only : findMultiFilename

    implicit none

    character(LEN=*), intent(in)  :: rootfilename
    character(LEN=80) :: filename
    character(len=20) :: word(40), unit(40)
    integer :: nword, nunit
    !   
    integer, parameter  :: LUIN = 10 ! logical unit # of the data file
    real(double) :: udist, umass, utime,  time, uvel, utemp
    real(double) :: xn, yn, zn, vx, vy, vz, gaspartmass, rhon, u, h, h2ratio, gmw
    integer :: itype ! splash particle type, different convention to SPHNG 
    integer :: ipart, icount, iptmass, igas, idead
    integer :: nptmass, nghost, nstar, nunknown, nlines
    real(double) :: junkArray(50) !, junk
    character(LEN=1)  :: junkchar
    character(LEN=150) :: message
    character(len=500) :: namestring, unitString
    integer :: ix, iy, iz, ivx, ivy, ivz, irho, iu, iitype, ih, imass, iUoverT
    logical :: haveUandUoverT

!
! For SPH simulations with chemistry
!
! CO fraction. The column will not be labelled by splash
! so need to specify which column to use.
    logical :: useCO
    integer, parameter :: iCO = 15
    real(double)       :: COfrac
! Do we need particle H2?
    logical :: useH2

! Account for time in Galactic plane surveys
    logical, parameter :: multiTime=.false.
    real(double) :: extraPA

    if (sphwithchem) then 
       useCO = .true.
       useH2 = .true.
    else
       useCO = .false.
       useH2 = .false.
    endif

    call findMultiFilename(rootfilename, iModel, filename)
    open(unit=LUIN, file=TRIM(filename), form="formatted",status="old")
    read(LUIN,*) 
    read(LUIN,*)
    read(LUIN,*)
    read(LUIN,*) junkchar, time, utime
    read(LUIN,*)
    read(LUIN,*)
    read(LUIN,*) junkchar, npart, nghost, nptmass, nstar, nunknown
    read(LUIN,*)
!    read(LUIN,*) junkchar, udist, junk, junk, umass, junk, junk, junk, junk, junk, uvel, junk, junk, utemp
    read(LUIN,'(a)') unitString
    unitString = unitstring(2:)
!    read(LUIN,*) junkchar, udist, junk, junk, umass
    read(LUIN,*)
    read(LUIN,*)
    read(LUIN,'(a)') nameString
    nameString = nameString(2:)
    call splitintoWords(unitString, unit, nunit)
    call splitIntoWords(nameString, word, nWord)

    if (nUnit /= nWord-1) then
       call writeFatal("Error reading splash file")
       write(*,*) nUnit, nWord
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

    ivx = indexWord("v\dx",word,nWord)
    ivy = indexWord("v\dy",word,nWord)
    ivz = indexWord("v\dz",word,nWord)

    read(unit(ivx),*) uvel

    imass = indexWord("particle mass",word,nWord)

    read(unit(imass),*) umass

    iu = indexWord("u",word,nWord)

    read(unit(iu),*) utemp
    
    irho = indexWord("density",word,nWord)
    ih = indexWord("h",word,nWord)
    iitype = indexWord("itype",word,nWord)

    write(message,*) "Allocating ", npart-nptmass, " gas particles and ", nptmass, " sink particles"
    call writeinfo(message, TRIVIAL)
    call init_sph_data(udist, umass, utime, time, nptmass, uvel, utemp)
    ! velocity unit is derived from distance and time unit (converted to seconds from years)
    sphdata%codeVelocitytoTORUS = uvel / cspeed 

! Check whether this SPH dump has u and u/T  available (e.g. SPH with radiative transfer). 
! If they are both present then set the temperature from these columns. 
    if ( wordIsPresent("u",word,nWord) .and. wordIsPresent("u/T",word,nWord) ) then
       call writeinfo("Found u and u/T in the SPH file",TRIVIAL)
       iUoverT = indexWord("u/T",word,nWord)
       haveUandUoverT = .true.
    else
       haveUandUoverT = .false.
    end if

! Decide how to convert between internal energy and temperature
    if (convertRhoToHI) then
       write (message,'(a,i2)') "Converting density to HI density using H2 fraction from column ", ih2frac
       call writeInfo(message,FORINFO)
       ! This affects the conversion of internal energy to temperature via the mean molecular weight
       ! so needs to be handled differently to cases with fixed abundances. 
       sphdata%codeEnergytoTemperature = 1.0
    else if (haveUandUoverT) then
       call writeInfo("Calculating temperature from u and u/T",FORINFO)
       sphdata%codeEnergytoTemperature = 1.0
    else
       sphdata%codeEnergytoTemperature = utemp * 1.9725e-8 ! temperature from molcluster! 2. * 2.46 * (u * 1d-7) / (3. * 8.314472)
       write(message,*) "Conversion factor between u and temperature (assumes molecular weight of 2.46): ", &
            sphdata%codeEnergytoTemperature
       call writeInfo(message, FORINFO)
    end if

    if (useCO) then
       write (message,'(a,i2)') "Reading CO fraction from column ", iCO
       call writeInfo(message,FORINFO)
       allocate(sphData%rhoCO(npart))
    end if

    if (useH2) then 
       write (message,'(a,i2)') "Will store particle H2 density. H2 fraction is from column ", ih2frac
       call writeInfo(message,FORINFO)
       allocate(sphdata%rhoH2(npart))
    end if

    sphData%useSphTem = .true.
    sphdata%totalgasmass = 0.d0

    nlines = npart + nghost + nptmass + nstar + nunknown ! npart now equal to no. lines - 12 = sum of particles dead or alive

    write(message,*) "Reading SPH data from ASCII...."
    call writeinfo(message, TRIVIAL)

    iptmass = 0
    icount = 12 ! header lines
    igas = 0
    idead = 0

    do ipart=1, nlines

!       read(LUIN,*) xn, yn, zn, gaspartmass, h, rhon, junk, junk, junk, vx, vy, vz, u, junk, junk, junk, junk, junk, junk, junk, &
!            junk, junk, junk,junk, junk,junk,junk,junk,itype
       read(LUIN,*) junkArray(1:nWord)

       xn = junkArray(ix)
       yn = junkArray(iy)
       zn = junkArray(iz)
       vx = junkArray(ivx)
       vy = junkArray(ivy)
       vz = junkArray(ivz)

       if (internalView) call rotate_particles(galaxyPositionAngle+extraPA, galaxyInclination)

       u = junkArray(iu)
       rhon = junkArray(irho)
       h = junkArray(ih)
       itype = int(junkArray(iitype))
       gaspartmass = junkArray(imass)

       icount = icount + 1

! 1=gas particle
       if (itype == 1) then 
          igas = igas + 1

          sphdata%xn(igas) = xn
          sphdata%yn(igas) = yn
          sphdata%zn(igas) = zn


          sphdata%gasmass(igas) = gaspartmass

          sphdata%vxn(igas) = vx
          sphdata%vyn(igas) = vy
          sphdata%vzn(igas) = vz

! For SPH simulations with chemistry we need to set up H2
          if ( convertRhoToHI .or. useH2 ) then
             h2ratio = junkArray(ih2frac)
          endif

          if (convertRhoToHI) then 
             gmw = (2.*h2ratio+(1.-2.*h2ratio)+0.4) / (0.1+h2ratio+(1.-2.*h2ratio))
             sphdata%temperature(igas) = (2.0/3.0) * u * ( gmw / Rgas) * utemp
          else if (haveUandUoverT) then 
             sphdata%temperature(igas) = u / junkArray(iUoverT)
          else
             sphdata%temperature(igas) = u
          end if

          sphdata%hn(igas) = h

          if(rhon .gt. 0.d0) then
             sphdata%rhon(igas) = rhon
          else
             sphdata%rhon(igas) = gaspartmass / h**3
          endif

          sphdata%totalgasmass = sphdata%totalgasmass + gaspartmass

          if ( convertRhoToHI ) then
             sphdata%rhon(igas) = (1.0-2.0*h2ratio)*rhon*5.0/7.0
          end if

          if (useCO) then 
             COfrac = junkArray(iCO)
             sphdata%rhoCO(igas) = COfrac * rhon
          endif

          if (useH2) then
             sphdata%rhoH2(igas) = h2ratio*2.*rhon*5./7.
          end if
          
! 3=sink, 4=star
       else if(itype .eq. 3 .or. itype .eq. 4) then

          iptmass = iptmass + 1

          sphdata%x(iptmass) = xn
          sphdata%y(iptmass) = yn
          sphdata%z(iptmass) = zn

          sphdata%vx(iptmass) = vx
          sphdata%vy(iptmass) = vy
          sphdata%vz(iptmass) = vz

          sphdata%ptmass(iptmass) = gaspartmass 
          sphdata%hpt(iptmass) = h

!          sphdata%totalgasmass = sphdata%totalgasmass + gaspartmass

          write(message,*) "Sink Particle number", iptmass," - mass", gaspartmass, " Msol - Index", iptmass + igas
          call writeinfo(message, TRIVIAL)
          write(98,*) iptmass, xn*udist*1e-10, yn*udist*1e-10, zn*udist*1e-10, gaspartmass

! 2=ghost, 5=dead/unknown
       else

          idead = idead + 1
       
       endif

 enddo

 write(message,*) "Read ",icount, " lines"
 call writeinfo(message, TRIVIAL)

 write(message,*)  iptmass," are sink particles and ",igas," are gas particles and ", idead, " are dead"
 call writeinfo(message, TRIVIAL)

 write(message,*) "Total Mass in all particles, ", sphdata%totalgasmass * umass/mSol, " Msol"
 call writeinfo(message, TRIVIAL)
    
 close(LUIN)
   
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

! Read in SPH data from an MPI dump file
! D. Acreman, June 2012
  subroutine read_sph_data_mpi(filename)
    use inputs_mod, only: amrgridcentrex, amrgridcentrey, amrgridcentrez, amrgridsize, splitovermpi
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
    REAL(kind=4)     :: r4
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

    integer,allocatable :: listpm(:)
    real(kind=8),allocatable    :: spinx(:)
    type(VECTOR) :: centre
    real(double) :: halfSize

    integer :: thisNumGas
    integer :: nlistinactive, listinactive(1)
    integer :: iThread
#ifdef MPI
    integer :: ierr
    integer, parameter :: tag = 33
    integer :: status(MPI_STATUS_SIZE)
#endif
! Mean molecular weight, used for calculating temperature from internal energy.
! This assumes a 10:1 H:He ratio by number
    real(double), parameter :: gmw = 14.0/11.0 

#ifdef MPI
    if (myrankWorldGlobal == 0) then
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
    open(unit=LUIN, file=filename, form='unformatted', status='old')

    read(LUIN) int1, r1, int2, i1, int1
    read(LUIN) fileident
    write(message,*) "fileident=", fileident
    call writeinfo(message, TRIVIAL)
    
    read(LUIN) number
    IF (number==6) THEN
       READ (LUIN) npart,n1,n2,nreassign,naccrete,nkill
       nblocks = 1
    ELSE
       read(LUIN) npart, n1, n2, nreassign, naccrete, nkilltot, nblocks
    ENDIF       
    write(message,*) "Total number of particles= ", npart
    call writeInfo(message,TRIVIAL)
    write(message,*) "Number of MPI blocks= ", nblocks
    call writeInfo(message,TRIVIAL)

! Allocate storage for data read from dump file
    allocate( isteps(npart) )
    allocate( iphase(npart) )
    allocate( xyzmh(5,npart))
    allocate( vxyzu(4,npart))
    allocate( uoverTarray(npart))
    allocate( rho(npart)    )

    do i=1,6
       read(LUIN)
    end do

    read(LUIN) time 
    write(message,*) "Dump time=", time
    call writeInfo(message,TRIVIAL)

    read(LUIN)
    read(LUIN)
    read(LUIN) udist, umass, utime
    read(LUIN) number

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

       READ (LUIN) (isteps(i), i=iiigas+1, iiigas+blocknpart)

       do j=1,nums(1)-1
          READ (LUIN) (nlistinactive, listinactive(i), i=1,1)
       END DO
       READ (LUIN) (iphase(i), i=iiigas+1, iiigas+blocknpart)

       thisNumGas = 0
       do i=iiigas+1, iiigas+blocknpart
          if (iphase(i)==0) thisNumGas = thisNumGas+1
       end do
       write (message,*) "There are ", thisNumGas, " active gas particles"
       call writeInfo(message,TRIVIAL)

       read(LUIN) 

       do j=1,5
          read(LUIN) ( xyzmh(j,i), i=iiigas+1,iiigas+blocknpart)
       end do

       do j=1,4
          read(LUIN) ( vxyzu(j,i), i=iiigas+1,iiigas+blocknpart)
       end do

       do j=1,nums(6)-9
          READ (LUIN)
       end do

       read(LUIN) ( rho(i), i=iiigas+1,iiigas+blocknpart) 

       do j=1,nums(7)-2
          READ (LUIN)
       end do
       read (LUIN) (r4, i=iiigas+1,iiigas+blocknpart)

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

       DO i = 1, numssink(6)-9
          READ (LUIN)
       END DO
       DO i = 1, numssink(8)
          READ (LUIN)
       END DO

       if (nblocktypes.GE.3) then
          call writeInfo ("Reading RT data")
          do j=1,2
             READ (LUIN)
          end do

          READ (LUIN) (uoverTarray(i),i=iiigas+1,iiigas+blocknpart)

          do j = 1, numsrt(6)-3
             READ (LUIN)
          end do
          call writeInfo ("Read RT data")
       endif

       if (nblocktypes.GE.4) then
          call writeInfo ("Reading MHD data")
          do i=1,8
             do j=1,nums(i)
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


#ifdef MPI
    do iThread = 1, nHydroThreadsGlobal
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
       if (iphase(i) > 0) then
          if ( particleInBox(xyzmh(:,i),uDist, centre, halfSize)) then 
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
    sphData%useSphTem = .true. 
    sphdata%codeVelocitytoTORUS = (udist / utime) / cspeed

    if (nblocktypes.GE.3) then
       sphdata%codeEnergytoTemperature = 1.0
    else
       sphdata%codeEnergytoTemperature = (2.0/3.0) * (gmw / Rgas) * ( (udist**2)/(utime**2) )
    endif

    iiigas=0
    iiisink = 0 



! This loop needs to be over all the particles read in
    do i=1, blocksum_npart
       if (iphase(i) == 0) then
          if ( particleRequired(xyzmh(:,i),uDist, centre, halfSize)) then 

             iiigas=iiigas+1

             sphData%rhon(iiigas)  = rho(i)

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

             sphData%xn(iiigas)          = xyzmh(1,i)
             sphData%yn(iiigas)          = xyzmh(2,i)
             sphData%zn(iiigas)          = xyzmh(3,i)
             sphData%gasmass(iiigas)     = xyzmh(4,i)
             sphData%hn(iiigas)          = xyzmh(5,i)

          endif
       else if(iphase(i) > 0) then
          if ( particleInBox(xyzmh(:,i),uDist, centre, halfSize)) then 
       
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


    call kill_sph_data()
 enddo
    deallocate(rho)
    deallocate(vxyzu)
    deallocate(xyzmh)
    deallocate(uoverTarray)
    deallocate(iphase)
    deallocate(isteps)


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





  end subroutine read_sph_data_mpi

  logical function particleRequired(positionArray, uDist, cen, halfsize)
    type(VECTOR) :: cen
    real(double) :: halfSize
    real(kind=8) :: positionArray(5)
    real(db) :: x, y, z, h
    real(db) :: distFromGridCentre, maxDist
    real(db), intent(in) :: uDist

    ! Particle position in Torus units
    x = positionArray(1) * uDist / 1.0e10
    y = positionArray(2) * uDist / 1.0e10
    z = positionArray(3) * uDist / 1.0e10
    h = positionArray(5) * uDist / 1.0e10

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

! Read in SPH data from dump file with Chemistry. Errors reading from file could be due to incorrect endian. 
  subroutine read_sph_data_withChem(filename)
    use inputs_mod, only: internalView, galaxyPositionAngle, galaxyInclination

    implicit none
    
    character(len=*), intent(in) :: filename
    character(LEN=150) :: message

    integer :: i, j, iiigas
    real(double) :: hI_mass

    real(double) :: udist, umass, utime,  time
    integer, parameter :: nptmass=0

    INTEGER(kind=4)  :: int1, int2, i1
    integer(kind=4)  :: number,n1,n2,nreassign,naccrete,nkilltot,nblocks
    REAL(kind=8)     :: r1, dummy
    integer :: intarray(8)
    CHARACTER(len=100) ::  fileident

    integer(kind=1), parameter  :: LUIN = 10 ! logical unit # of the data file

    integer(kind=1), allocatable :: iphase(:)
    integer, allocatable   :: isteps(:)
    real(kind=8), allocatable    :: xyzmh(:,:)
    real(kind=8), allocatable    :: vxyzu(:,:)
    real(kind=4), allocatable    :: rho(:)
    real(kind=8), allocatable    :: h2rho(:)
    real(kind=8), allocatable    :: h1rho(:)
    real(kind=8), allocatable    :: COrho(:)

    integer :: status
    type(VECTOR) :: orig_sph, rot_sph

!
! Read in data from file
!

    write(message,*) "Reading SPH data from "//trim(filename)
    call writeinfo(message, TRIVIAL)

    open(unit=LUIN, file=filename, form='unformatted', status='old')

    read(LUIN) int1, r1, int2, i1, int1
    read(LUIN) fileident

    read(LUIN) number
    read(LUIN) npart, n1, n2, nreassign, naccrete, nkilltot, nblocks

    do i=1,4
       read(LUIN) number
    end do

    read(LUIN) number
    
    read(LUIN) dummy, dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,&
         dummy,dummy,dummy,dummy,dummy,dummy

    read(LUIN) number

    read(LUIN) number
    read(LUIN) udist, umass, utime, dummy

    time=0.0

    call init_sph_data(udist, umass, utime, time, nptmass)
    allocate(sphdata%rhoH2(npart)) 

    sphdata%codeVelocitytoTORUS = (udist / utime) / cspeed
    sphdata%codeEnergytoTemperature = 1.0 ! no conversion required

! Arrays

    read(LUIN) number, nblocks

! Array length 1 header
    read(LUIN) dummy, intarray(1:8)

! Array length 2 header
    read(LUIN) dummy, intarray(1:8)

! isteps and iphase
    allocate(isteps(npart))
    allocate(iphase(npart))
! Note: this should not be a do loop
    read(LUIN) ( isteps(i), i=1,npart)
    read(LUIN) ( iphase(i), i=1, npart)

! xyzmh, vxyzu
    allocate( xyzmh(5,npart) )
    allocate( vxyzu(4,npart) )

    do j=1,5
       read(LUIN) ( xyzmh(j,i), i=1, npart)
    end do

    do j=1,4
       read(LUIN) ( vxyzu(j,i), i=1, npart)
    end do

! rho
    allocate(rho(npart))
    read(LUIN) ( rho(i), i=1, npart) 

! h2rho and h1rho
    allocate(h2rho(npart))
    allocate(h1rho(npart))
    read(LUIN) ( h2rho(i), i=1, npart)
    read(LUIN) ( h1rho(i), i=1, npart)

! See if there is CO data in this dump
    allocate(COrho(npart))
    read(LUIN, iostat=status) ( COrho(i), i=1,npart)
    if (status == 0) then 
       call writeInfo ("Found CO data")
       allocate(sphdata%rhoCO(npart))
    else
       call writeInfo("No CO data found in this dump")
    end if

    close(LUIN)


!
! Set up SPH data structure
!

    ! Indicate that this object is in use
    sphData%inUse = .true.

! Use temperature from SPH particles to initialise temperature grid
    sphData%useSphTem = .true. 
    sphData%npart = npart
    sphData%time = time

    iiigas=0
    hI_mass = 0.0
    do i=1, npart
       if (iphase(i) == 0) then

          iiigas=iiigas+1
          hI_mass = hI_mass + xyzmh(4,i) * (h1rho(i) / rho(i))

          sphData%rhon(iiigas)  = h1rho(i)
          sphData%rhoH2(iiigas) = h2rho(i)
          if ( associated (sphData%rhoCO) ) sphData%rhoCO(iiigas) = COrho(i)

          if ( internalView ) then 

! For the internal view rotate the galaxy so that we are not looking along cell boundaries
             orig_sph = VECTOR( vxyzu(1,i),  vxyzu(2,i),  vxyzu(3,i) )
             rot_sph  = rotateZ( orig_sph,  galaxyPositionAngle*degToRad )
             rot_sph  = rotateY( rot_sph, galaxyInclination*degToRad )

             sphdata%vxn(iiigas)         = rot_sph%x
             sphdata%vyn(iiigas)         = rot_sph%y
             sphdata%vzn(iiigas)         = rot_sph%z
             sphData%temperature(iiigas) = vxyzu(4,i)
             
             orig_sph = VECTOR( xyzmh(1,i),  xyzmh(2,i),  xyzmh(3,i) )
             rot_sph  = rotateZ( orig_sph,  galaxyPositionAngle*degToRad )
             rot_sph  = rotateY( rot_sph, galaxyInclination*degToRad )
             
             sphData%xn(iiigas)          = rot_sph%x
             sphData%yn(iiigas)          = rot_sph%y
             sphData%zn(iiigas)          = rot_sph%z
             sphData%gasmass(iiigas)     = xyzmh(4,i)
             sphData%hn(iiigas)          = xyzmh(5,i)

          else

             sphdata%vxn(iiigas)         = vxyzu(1,i)
             sphdata%vyn(iiigas)         = vxyzu(2,i)
             sphdata%vzn(iiigas)         = vxyzu(3,i)
             sphData%temperature(iiigas) = vxyzu(4,i)

             sphData%xn(iiigas)          = xyzmh(1,i)
             sphData%yn(iiigas)          = xyzmh(2,i)
             sphData%zn(iiigas)          = xyzmh(3,i)
             sphData%gasmass(iiigas)     = xyzmh(4,i)
             sphData%hn(iiigas)          = xyzmh(5,i)

          end if

       end if
    end do

    write(message,*) "Read ", iiigas, " gas particles"
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

    write(message,*) "Total HI mass=", hI_mass * umass / mSol
    call writeinfo(message, TRIVIAL)

    deallocate(h1rho)
    deallocate(h2rho)
    deallocate(COrho)
    deallocate(rho)
    deallocate(vxyzu)
    deallocate(xyzmh)
    deallocate(iphase)
    deallocate(isteps)

  end subroutine read_sph_data_withChem


  !
  !
  ! Read in the data from a file, allocate the array memory, and store the
  ! the number of gas and stars in this object.
  !
  subroutine read_stellar_disc_data(this, filename)
    implicit none
    type(sph_data), intent(inout) :: this
    character(LEN=*), intent(in)  :: filename
    ! 
    integer, parameter  :: LUIN = 10 ! logical unit # of the data file
    integer :: nstar
    integer :: i, dum_i
    character(LEN=1) :: dum_a
    real(double), parameter :: M_sun = 1.989e33 ! [grams]

    open(unit=LUIN, file=TRIM(filename), status='old')

    nstar = this%nptmass

    ! reading in the header
    do i = 1, 6
       READ(LUIN, *) dum_a
    end do


    ! allocating arrays
    ! --- for stellar discs
    ALLOCATE(this%discrad(nstar))
    ALLOCATE(this%discmass(nstar))
    ALLOCATE(this%spinx(nstar))
    ALLOCATE(this%spiny(nstar))
    ALLOCATE(this%spinz(nstar))


    write(*,*) ' '
    write(*,*) 'Reading Matthew''s stellar disc data....'
    write(*,*) ' '

    
    do i=1, nstar 
       read(luin, *) dum_i, this%discrad(i), this%discmass(i), &
            this%spinx(i), this%spiny(i), this%spinz(i)
       ! convert the mass into grams
       this%discmass(i) = this%discmass(i)*M_sun  ![g]
    end do


    this%have_stellar_disc = .true.    


  end subroutine read_stellar_disc_data
  

    

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

  ! Returns the temperature of gas particle at the postion of
  ! i-th particle.
  
  function get_temp(i) RESULT(out)
    implicit none
    real(double) :: out 
    Integer, intent(in) :: i

    out  = sphdata%temperature(i)
    
  end function get_temp

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
    




  !
  ! find the maximum distance between the pt masses
  !

  function max_distance() RESULT(out)
    implicit none
    real(double) :: out
    !
    real(double) :: d_max, d, x, y, z
    integer :: i, j, n
    
    d_max= -1.0
    n=get_nptmass() ! function in this moudle
    
    do i = 1, n-2
       do j = i+1, n
          call get_position_pt_mass(j, x, y, z)
          d = (x*x+y*y+z*z)  ! omit SQRT here cus it costs too much.
          d_max = MAX(d_max, d)
       end do
    end do 
    d_max = SQRT(d_max)
    out = d_max
  end function max_distance


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
  subroutine info(filename)
    implicit none
    character(LEN=*), intent(in) :: filename
    integer :: UN
    real(double) :: tmp
    
    if (filename(1:1) == '*') then
       UN = 6   ! prints on screen
    else
       UN = 69
       open(unit=UN, file = TRIM(filename), status = 'replace')
    end if

    ! time of the data dump
    tmp = get_time()*get_utime()/(60.0d0*60.0d0*24.0d0*365.0d0*1.0d6)
    
    Write(UN,'(a)') ' '
    write(UN,'(a)') '######################################################'
    write(UN,'(a)') 'SPH data info :'
    write(UN,'(a)') ' '    
    write(UN,*)     'Units of length            : ', get_udist(), ' [cm]'
    write(UN,*)     'Units of mass              : ', get_umass(), ' [g]'
    write(UN,*)     'Units of time              : ', get_utime(),  ' [s]' 
    write(UN,'(a)') ' '    
    write(UN,*)     '# of stars                 : ',  get_nptmass()
    write(UN,*)     '# of gas particles (total) : ',  get_npart()   
    write(UN,*)     'Time of data dump          : ',  tmp, ' [Myr]'    
    write(UN,'(a)') '#######################################################'
    write(UN,'(a)') ' '
    write(Un,'(a)') ' '
    
    if (filename(1:1) /= '*')  close(UN)

  end subroutine info



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

  TYPE(vector)  function Clusterparameter(point, thisoctal, subcell, theparam, isdone)
    USE inputs_mod, only: hcritPercentile, hmaxPercentile, sph_norm_limit, useHull
    USE constants_mod, only: tcbr
    use octal_mod, only: OCTAL

    type(vector), intent(in) :: point
    type(vector) :: posvec
    type(vector),save :: oldpoint = VECTOR(0.d0,0.d0,0.d0)
    type(vector),save :: oldvel = VECTOR(2.d0,2.d0,2.d0)

    integer :: i

    integer, optional :: theparam
    integer :: param

    logical, optional :: isdone
    logical :: done
    
    integer, allocatable, save :: ind(:)
    Real(double), save :: hcrit, hmax, OneOverHcrit, OneOverhMax
    real(double) :: codeVelocitytoTORUS, codeLengthtoTORUS, codeDensitytoTORUS, udist, umass, utime
    real(double), save :: r
    real(double) :: d
    real(double) :: fac
    real(double), save :: sumWeight
    real(double) :: bigx, bigy, bigz, Hedge
    real(double) :: paramValue(4)
    
    character(len=100) :: message

    logical, save :: notfound
    logical :: reuse

    integer :: rcounter

    type(octal), pointer, save :: previousOctal => null()
    type(octal), pointer :: thisOctal

    logical, save :: firsttime = .true.
    integer :: subcell
    integer, save :: prevsubcell

    if(present(theparam)) then
       param = theparam
    else
       param = 1
    endif

    if(present(isdone)) then
       done = isdone
    else
       done = .false.
    endif

    if(firsttime .and. .not. done) then

       if (get_npart() > 0) then
          udist = get_udist()
          utime = get_utime()
          umass = get_umass()
          codeLengthtoTORUS = udist * 1d-10
          codeVelocitytoTORUS = sphdata%codeVelocitytoTORUS
          codeDensitytoTORUS = umass / ((udist) ** 3)
          
          allocate(PositionArray(3,npart)) ! allocate memory
          allocate(xArray(npart))
          allocate(q2Array(npart))
          allocate(etaArray(npart)) ! added to fix smoothing length discrepancy - h != 1.2 (rho/m)^(1/3)
          allocate(RhoArray(npart))
          allocate(TemArray(npart))
          allocate(RhoH2Array(npart))
          allocate(Harray(npart))
          allocate(ind(npart))
          allocate(OneOverHsquared(npart))
          if (associated(sphData%rhoCO)) allocate (rhoCOarray(npart))
          
          allocate(HullArray(npart))
          
          
          PositionArray = 0.d0; hArray = 0.d0; ind = 0; q2array = 0.d0; etaarray = 0.d0
          if (.not.associated(sphdata%xn)) then
             write(*,*) "sphdata not associated"
             write(*,*) "rank ",myrankGlobal
             write(*,*) "npart ",sphData%npart
             write(*,*) "ptmass ",sphData%nptmass
          endif
          
          Positionarray(1,:) = sphdata%xn(:) * codeLengthtoTORUS! fill with x's to be sorted
          xArray(:) = sphdata%xn(:) * codeLengthtoTORUS! fill with x's to be sorted
          
          call sortbyx(xarray(:),ind(:)) ! sort the x's and recall their indices
          
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
          etaarray(:) = ((50.d0 / npart) / rhoarray(:)) / harray(:)**3
          !       write(*,*) sum(etaarray(:)) / npart, sqrt(sum(etaarray(:)**2)/npart), &
          !            sqrt(sum(etaarray(:)**2)/npart) - (sum(etaarray(:)) / npart)
          
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
          
          if (associated (rhoCOarray) .and. associated(sphdata%rhoCO) ) then 
             rhoCOarray(:) = sphdata%rhoCO(ind(:)) * codeDensityToTorus
          end if
          
          hcrit = hcrit * codeLengthtoTORUS
          OneOverHcrit = 1.d0 / hcrit
          
          hmax = hmax * codeLengthtoTORUS
          OneOverHmax = 1.d0 / hmax
          
          Harray(:) = Harray(:) * codeLengthtoTORUS  ! fill h array
          
          bigx = maxval(PositionArray(1,:)) - minval(PositionArray(1,:))
          bigy = maxval(PositionArray(2,:)) - minval(PositionArray(2,:))
          bigz = maxval(PositionArray(3,:)) - minval(PositionArray(3,:))
          
          Hedge = ((bigX * bigY * bigZ / npart) ** (1.d0/3.d0))
          
          do i=1, npart
             if(Harray(i) .ge. Hedge) then
                HullArray(i) = .true.
             else
                HullArray(i) = .false.
             endif
          enddo
          
          if (useHull) then 
             write(message,*) "Hull smoothing Length in 10^10cm ", Hedge
             call writeinfo(message, TRIVIAL)
             
             write(message,*) "Percentage of Hull particles ", 100. * real(count(HullArray)) / real(npart)
             call writeinfo(message, TRIVIAL)
          endif
          
          OneOverHsquared(:) = 1.d0 / (Harray(:)**2)

          write(message,*) "Critical smoothing Length in 10^10cm", hcrit
          call writeinfo(message, TRIVIAL)
          write(message,*) "Maximum smoothing Length in 10^10cm", hmax
          call writeinfo(message, TRIVIAL)
          
          rcrit = 2.d0 * hcrit ! edge of smoothing sphere
          rmax = 2.d0 * hmax ! edge of smoothing sphere
       
          allocate(partarray(npart), indexarray(npart))
          
          firsttime = .false.
       endif
    endif
       
    if (get_npart() == 0) then
       if(param .eq. 1) then
          Clusterparameter = VECTOR(0.d0,0.d0,0.d0)
       elseif(param .eq. 2) then
          Clusterparameter = VECTOR(1d-37, tcbr, 1d-37)  ! density ! stays as vector for moment
       elseif(param .eq. 3) then 
          Clusterparameter = VECTOR(1d-99, 0.d0, 0.d0)
       endif
       goto 666
    endif


    if(done) then
       if (allocated(PositionArray)) then
          deallocate(xArray, PositionArray, harray, RhoArray, Temarray, ind, &
               q2Array, HullArray, etaarray, RhoH2Array)
          
          deallocate(OneOverHsquared)
          deallocate(partarray, indexarray)

          if (associated(rhoCOarray)) deallocate(rhoCOarray)
          if (allocated(VelocityArray)) deallocate (VelocityArray)
       endif
       nullify(previousOctal)
       firsttime = .true.
       goto 666
    endif
              
    posVec = point

    d = min(thisoctal%h(subcell), rmax)! using the placeholder h from splitgrid

!   d = thisoctal%subcellsize ! the splitgrid routine effectively picks the grid size based on smoothing length (mass condition)

    ! hope the compiler optimizes this (for a disc, first condition is simpler to test)
    if(abs(posvec%z) .gt. maxz + rmax .or. &
         (posvec .dot. posvec) .gt. maxr2) then ! if point is outside particles then there will be no information
       if(param .eq. 1) then
          if(modulus(posvec - oldpoint) .lt. 2. * d) then! cells nearby?
             Clusterparameter = oldvel ! Velocity of nearby cells
          else
 ! 2c is not a real velocity and everything should be trying to handle this error
             Clusterparameter = VECTOR(2.d0,2.d0,2.d0)
          endif
          return
       elseif(param .eq. 2) then
          Clusterparameter = VECTOR(1d-37, tcbr, 1d-37)  ! density ! stays as vector for moment
          return
       elseif(param .eq. 3) then 
          Clusterparameter = VECTOR(1d-99, 0.d0, 0.d0)
          return
       endif
    endif

    if(posvec .eq. prevpos) goto 500

    if(associated(previousoctal)) then
       if(thisoctal%centre .eq. previousoctal%centre .and. nparticles .ne. 0 .and. subcell .eq. prevsubcell) then 
          reuse = .true.
       else
          reuse = .false.
       endif
    else
       reuse = .false.
    endif

    previousoctal => thisoctal
    prevsubcell = subcell

    notfound = .false.
    nparticles = 0
    rcounter = 0
    sumweight = 0.d0

    if (npart == 0) then
       write(*,*) "Error npart is zero"
       write(*,*) "rank ",myrankGlobal
       write(*,*) "sph ",sphdata%npart
    endif

    if(reuse) then
       notfound = .false.
       call findNearestParticles(posvec, nparticles, r, shouldreuse = reuse) ! redo but with exponential kernel
       if(nparticles .gt. 0) call doweights(sumweight)
    else
       ! tradeoff time for accuracy in low density regions
       r = min(4.d0 * d, rcrit) ! 4d is far enough away to have particles with their 
       !smoothing lengths captured, rcrit puts an upper limit on time
       ! 2 rcrit is essential for the mass to be correctly done... (in my case it had to be hcrit = 99%)
 
       do while (sumweight .le. 1d-3)
          call findNearestParticles(posvec, nparticles, r) ! redo but with exponential kernel
          if(nparticles .gt. 0) call doweights(sumweight)
          if(r .ge. rmax * 0.5d0) exit
          rcounter = rcounter + 1
          if(rcounter .eq. 1) then
             r = min(max(r * 4.d0, 2.d0 * rcrit), rmax * 0.1)
          elseif(rcounter .eq. 2) then
             r = rmax * 0.5d0
          endif
       enddo

    endif
    
500 continue

    prevpos = point

    if(sumweight .le. 0.d0) notfound = .true.
    
    paramvalue(:) = 0.d0

    if(.not. notfound) then
!       thisoctal%inflow(subcell) = .false. ! octal is occupied by something at least cf. inflow = .true.
       if(param .eq. 1) then

          fac = 1.d0 / sumweight
!          fac = 1.d0

          do i = 1, nparticles
             Paramvalue(1) = paramValue(1) + partArray(i) * VelocityArray(1, indexArray(i)) ! Vx
             paramValue(2) = paramValue(2) + partArray(i) * VelocityArray(2, indexArray(i)) ! Vy
             paramValue(3) = paramValue(3) + partArray(i) * VelocityArray(3, indexArray(i)) ! Vz
          enddo
          
          Clusterparameter = VECTOR(paramValue(1) * fac, paramValue(2) * fac, paramValue(3) * fac) ! Velocity 
          oldvel = clusterparameter ! store this point in case it's nearby an empty cell (for velocity)
          oldpoint = posvec

       elseif(param .eq. 2) then

          if ( useHull ) then 
             if(any(HullArray(indexarray(1:nparticles)))) then
                fac = 1.d0
             else
                fac = 1.d0 / sumWeight
             endif
          
             if(sumweight .gt. sph_norm_limit) then
                HullArray(indexarray(1:nparticles)) = .false.
             else
                HullArray(indexarray(1:nparticles)) = .true.
             endif
          else
             if(sumweight .gt. sph_norm_limit) then
                fac = 1.d0 / sumWeight
             else
                fac = 1.0
             end if
          end if
 

          do i = 1, nparticles
             paramValue(3) = paramValue(3) + partArray(i) * TemArray(indexArray(i)) ! Temperature
             paramValue(4) = paramValue(4) + partArray(i) * RhoArray(indexArray(i)) ! rho
             paramValue(1) = paramValue(1) + partArray(i) * rhoH2Array(indexArray(i)) ! H2 density
          enddo
          
          Clusterparameter = VECTOR(paramValue(4)*fac, paramValue(3)*fac, paramValue(1)*fac)  ! density ! stays as vector for moment
          

       elseif (param .eq. 3) then 

          if(sumweight .gt. sph_norm_limit) then
             fac = 1.d0 / sumWeight
          else
             fac = 1.0
          end if

          do i = 1, nparticles
             paramValue(1) = paramValue(1) + partArray(i) * rhoCOArray(indexArray(i)) ! CO density
          enddo
          
          Clusterparameter = VECTOR(paramValue(1)*fac, 0.d0, 0.d0)  ! density ! stays as vector for moment

       endif
    else
!       thisoctal%inflow(subcell) = .true. ! what to do about stuff with nothing in it...
       if(param .eq. 1) then
          if(modulus(posvec - oldpoint) .lt. 2. * d) then! cells nearby?
             Clusterparameter = oldvel ! Velocity of nearby cells
          else
 ! 2c is not a real velocity and everything should be trying to handle this error
             Clusterparameter = VECTOR(2.d0,2.d0,2.d0)
          endif
       elseif(param .eq. 2) then
          Clusterparameter = VECTOR(1d-37, tcbr, 1d-37)  ! density, temperature and H2 density 
       elseif(param .eq. 3) then
          Clusterparameter = VECTOR(1d-99, 0.d0, 0.d0) ! CO density
       endif
    endif
666 continue
  End function Clusterparameter

  subroutine findnearestparticles(pos, partcount, r, shouldreuse)

    use inputs_mod, only : kerneltype
    use utils_mod, only: locate_double_f90

    type(VECTOR) :: pos
    real(double) :: x,y,z
    integer :: i
    integer, save :: nupper, nlower
    integer, save :: closestXindex, testIndex
    integer, intent(out) :: partcount
    real(double) :: r2test, q2test, r
    real(double), save :: rtest
    real(double) :: ydiff, zdiff
    integer :: stepsize, sense
    logical :: test, prevtest, up
    logical, optional :: shouldreuse
    logical :: reuse
    real(double) :: fac, fac2

    if(kerneltype .eq. 1) then
       fac = 2.d0
       fac2 = 4.d0
    else
       fac = 5.d0
       fac2 = 25.d0
    endif
    
    if(present(shouldreuse)) then
       reuse = shouldreuse
    else
       reuse = .false.
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

    if(reuse) then
!       rtest = max(2.d0 * maxval(harray(closestXindex+nlower:closestXindex+nupper)),rmax)
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

    rtest = r
1001 continue

    partcount = 0          

    do i = nlower, nupper ! search over all those particles we just found
       
       testIndex = closestXindex + i
             
       ydiff = abs(PositionArray(2,testindex) - y)

       if(ydiff .le. rtest) then ! if it's near in y then
          zdiff = abs(PositionArray(3,testIndex) - z)

          if(zdiff .le. rtest) then !only if it's near in z as well then work out contribution

             r2test = (xarray(testIndex) - x)**2 + & !The kernel will do the rest of the work for those outside the sphere
                  ydiff ** 2 + zdiff ** 2

             q2test = r2test * OneOverHsquared(testIndex) ! dimensionless parameter that we're interested in
             
             if(q2test .lt. fac2) then
                partcount = partcount + 1
                indexArray(partcount) = testIndex
                q2array(partcount) = q2test
             endif
          endif
       endif
    enddo
!666 continue
  end subroutine findnearestparticles

  Subroutine doWeights(sumweight)

    use inputs_mod, only : kerneltype

!    real(double), parameter :: num = 0.578703703d0 ! (5/6)^3 ! (1/1.2^3)
    real(double), parameter :: scalar = 0.103927732d0 ! 5 / 6 * one over sqrtpicubed

    real(double) :: sumweight

    integer :: i
    real(double) :: sqrtq

    if(nparticles .gt. 0) then
              
       if(kerneltype .eq. 0) then
       
          partarray(1:nparticles) = scalar * exp(-q2array(1:nparticles))
       
       elseif( kerneltype .eq. 1) then
          
          do i = 1, nparticles
             sqrtq = sqrt(q2array(i))
!             partarray(i) = num * SmoothingKernel3d(sqrtq)
             partarray(i) = etaarray(i) * SmoothingKernel3d(sqrtq)
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
  
 subroutine splitIntoWords(longString, word, nWord)
   character(len=20) :: word(:)
   integer :: nWord
   character(len=*) :: longString
   character(len=500) :: tempString
   logical :: stillSplitting

   tempString = ADJUSTL(longString)
   nWord = 0
   stillSplitting = .true.
   do while(stillSplitting)
      nWord = nWord + 1
      word(nWord) = tempString(1:16)
      tempString = tempString(17:)
      if (len(trim(tempString)) == 0) stillSplitting = .false.
   end do
 end subroutine splitIntoWords

 integer function indexWord(inputword, wordarray, nword) 
   character(len=*) :: inputWord
   character(len=20) :: wordArray(:)
   integer :: nWord, i
   
   indexWord = 0 
   do i = 1, nWord
      if (inputWord == wordArray(i)) then
         indexWord = i
         exit
      endif
   enddo

   if ( indexWord == 0 ) then 
      call writeFatal("indexWord: "//inputword//" not found")
   end if

 end function indexWord

! Check if the specified word is present in the list
 logical function wordIsPresent(inputword, wordarray, nword)
   character(len=*) :: inputWord
   character(len=20) :: wordArray(:)
   integer :: nWord, i
   
   wordIsPresent = .false.
   do i = 1, nWord
      if (inputWord == wordArray(i)) then
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
     use gridtype_mod, only : gridtype
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

end module sph_data_class

#endif
