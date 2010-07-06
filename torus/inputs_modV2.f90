module inputs_mod

  use vector_mod
  use unix_mod
  use messages_mod
  use kind_mod
  use input_variables
  use constants_mod
  use utils_mod, only: file_line_count

  implicit none

contains
  
  subroutine  inputs()

    implicit none

    integer :: nLines
    logical :: ok
    character(len=80), allocatable :: cLine(:) 
    character(len=80) :: message, thisLine
    logical :: done
    logical, allocatable :: fLine(:)

    ! character vars for unix environ

    character(len=80) :: dataDirectory

    character(len=80) :: paramFile

    integer :: error

    datadirectory = " "
    done = .false.
    ok = .true.
    oneKappa = .false.

    grainSize = 1.
    nDustType = 1
    freefreeImage = .false.
    filter_set_name = "natural"
    probDust = 0.1
    probContPhoton = 1.
    imagescale = 1.
    noscattering = .false.

    nBlobs = 0
    nLines = 0
    inputnMonte = 0

    call unixGetEnv("TORUS_JOB_DIR",absolutePath)
!   call get_environment_variable("TORUS_JOB_DIR",absolutePath)


    if (command_argument_count() == 1) then
       call get_command_argument(1, paramFile)
    else
       paramFile = trim(absolutePath)//"parameters.dat"
    endif
    if (writeoutput) then
       call writeInfo("Parameters file is: "//trim(paramFile),TRIVIAL)
    endif

    nLines = file_line_count(paramfile)
    ! nLines+1 element of cline is used in loop which reads parameter file
    allocate ( cLine(nLines+1) )
    nLines = 0 

    open(unit=32, file=paramfile, status='old', iostat=error)
    if (error /=0) then
       print *, 'Panic: parameter file open error, file: ',trim(paramFile) ; stop
    end if

    do
       nLines = nLines + 1
       read(32,'(a80)',end=10) cLine(nLines)
       thisLine = trim(adjustl(cLine(nLines)))
       if ( thisLine(1:1) == "%" .or. thisLine(1:1) == "!" ) nLines = nLines - 1   ! % or ! is a comment
       if (thisLine(1:1) == " ") nLines = nLines - 1
    end do
10  continue
    nLines = nLines - 1
    write(message,'(a,i4,a)') "Read in ", nLines, " parameters"
    call writeInfo(message,TRIVIAL)

    allocate(Fline(1:nLines))
    fLine = .false. 

    call getInteger("verbosity", verbosityLevel, cLine, fLine, nLines, &
         "Verbosity level: ", "(a,i8,1x,a)", 3, ok, .false.)

! the grid setup. Either we read the grid in or set it up from scratch

    call writeBanner("Grid setup parameters","*",TRIVIAL)

    call getLogical("splitovermpi", splitOverMPI, cLine, fLine, nLines, &
         "Grid is domain decomposed over MPI: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("debug", debug, cLine, fLine, nLines, &
         "Output debug information: ","(a,1l,1x,a)", .false., ok, .false.)


    call getLogical("readgrid", readGrid, cLine, fLine, nLines, &
         "Read grid file: ","(a,1l,1x,a)", .false., ok, .true.)

    if (.not.readgrid) then
       call readGridInitParameters(cLine, fLine, nLines)
       call readGeometrySpecificParameters(cLine, fLine, nLines)
    else
       call getString("inputfile", gridInputFilename, cLine, fLine, nLines, &
                  "Grid input filename: ","(a,a,1x,a)","none", ok, .true.)
    endif


! the physical ingredients of the model

    call writeBanner("Physical ingredients of model","*",TRIVIAL)

    call getLogical("dustphysics", dustPhysics, cLine, fLine, nLines, &
         "Include dust physics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("atomicphysics", atomicPhysics, cLine, fLine, nLines, &
         "Include atomic physics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("molecularphysics", molecularPhysics, cLine, fLine, nLines, &
         "Include molecular physics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("photoionphysics", photoionPhysics, cLine, fLine, nLines, &
         "Include photoionization physics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

!    if (.not.(dustPhysics.or.atomicPhysics.or.molecularPhysics.or.photoionPhysics)) then
!       call writeFatal("Must include one of: dustPhysics, atomicPhysics, molecularPhysics, photoionPhysics")
!    endif


    if (dustPhysics) call readDustPhysicsParameters(cLine, fLine, nLines)
    if (atomicPhysics) call readAtomicPhysicsParameters(cLine, fLine, nLines)
    if (molecularPhysics) call readMolecularPhysicsParameters(cLine, fLine, nLines)
    if (photoionPhysics) call readPhotoionPhysicsParameters(cLine, fLine, nLines)

    if (checkPresent("nsource", cLine, nLines)) then
       call readSourceParameters(cLine, fLine, nLines)
    endif


! the type of calculation

    call writeBanner("Type of model calculation","*",TRIVIAL)

    call getLogical("stateq", statisticalEquilibrium, cLine, fLine, nLines, &
         "Perform a statistical equililbrium calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("radeq", radiativeEquilibrium, cLine, fLine, nLines, &
         "Perform a radiative equilibrium calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("photoioneq", photoionEquilibrium, cLine, fLine, nLines, &
         "Perform a photoionization calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("hydrodynamics", hydrodynamics, cLine, fLine, nLines, &
         "Perform a hydrodynamics calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    if (statisticalEquilibrium.and.(.not.(molecularPhysics.or.atomicPhysics))) then
       call writeFatal("Must include either molecularPhysics or atomicPhysics for statistical equilibrium calculation")
    endif

    if (radiativeEquilibrium.and.(.not.dustPhysics)) then
       call writeFatal("Can only perform radiative equilibrium using dust physics")
    endif

    if (photoionEquilibrium.and.(.not.photoionPhysics)) then
       call writeFatal("Can only perform a photoionization calculation using photoionization physics")
    endif

    if (statisticalEquilibrium.and.molecularPhysics) call readMolecularLoopParameters(cLine, fLine, nLines)
    if (statisticalEquilibrium.and.atomicPhysics) call readAtomicLoopParameters(cLine, fLine, nLines)
    if (radiativeEquilibrium) call readRadiativeEquilibriumParameters(cLine, fLine, nLines)
    if (photoionEquilibrium) call readPhotoionEquilibriumParameters(cLine, fLine, nLines)
    if (hydrodynamics) call readHydrodynamicsParameters(cLine, fLine, nLines)

! now do we dump the output grid

    call writeBanner("Output file details","*",TRIVIAL)

    call getLogical("writegrid", writeGrid, cLine, fLine, nLines, &
         "Write grid file: ","(a,1l,1x,a)", .false., ok, .true.)

    if (writeGrid) then
       call getString("outputfile", gridOutputFilename, cLine, fLine, nLines, &
                  "Grid output filename: ","(a,a,1x,a)","none", ok, .true.)
    endif

! now we select the outputs

    call getLogical("datacube", calcDataCube, cLine, fLine, nLines, &
         "Calculate a data cube: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("image", calcImage, cLine, fLine, nLines, &
         "Calculate an image: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("spectrum", calcSpectrum, cLine, fLine, nLines, &
         "Calculate a spectrum: ","(a,1l,1x,a)", .false., ok, .false.)


    if (calcDataCube) call readDataCubeParameters(cLine, fLine, nLines)
    if (calcImage) call readImageParameters(cLine, fLine, nLines)
    if (calcSpectrum) call readSpectrumParameters(cLine, fLine, nLines)

    if (writeoutput) write(*,*) " "

    call writeUnusedKeywords(cLine, fLine, nLines)

    ! Close input file
    close(32)

    deallocate (cLine)

  end subroutine inputs

  subroutine writeUnusedKeywords(cLine, fline, nLines)
    character(len=*) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    integer :: i,k
    character(len=80) :: message
    
    if (ANY(.not.fLine(1:nLines))) then
       call writeBanner("Unused keywords found in parameters file","!")
       do i = 1, nLines
          if (.not.fLine(i)) then
             k = index(cline(i)," ")-1
             message = "Unused keyword: "//trim(cLine(i)(1:k))
             call writeWarning(message)
          endif
       enddo
    endif
  end subroutine writeUnusedKeywords

  subroutine readGeometrySpecificParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
    logical :: fLine(:)
    logical :: ok
    TYPE(VECTOR) :: gridCentre
    character(len=100) :: message

    select case(geometry)

       case("ttauri")
       call getReal("teff", teff, 1., cLine, fLine, nLines, &
            "Source temperature (K) : ","(a,f8.0,a)", 0.0, ok, .true.)
       call getReal("ttaurirstar", TTauriRstar, real(rsol), cLine, fLine, nLines, &
            "T Tauri stellar radius (in R_sol): ","(a,f7.1,1x,a)", 2.0, ok, .true.)
       rcore = TTauriRstar/1.0e10       ! [10^10cm]
       call getReal("ttaurimstar", TTauriMstar, real(msol), cLine, fLine, nLines, &
            "T Tauri stellar mass (in M_sol): ","(a,f7.1,1x,a)", 0.8, ok, .true.)
       call getReal("ttaurirouter", TTauriRouter, TTaurirStar, cLine, fLine, nLines, &
            "T Tauri outer flow radius (in R_star): ","(a,f7.1,1x,a)", 3.0, ok, .true.)
       call getReal("ttauririnner", TTauriRinner, TTaurirStar, cLine, fLine, nLines, &
            "T Tauri inner flow radius (in R_star): ","(a,f7.1,1x,a)", 2.2, ok, .true.)

       call getLogical("ttauridisc", ttauriDisc, cLine, fLine, nLines, &
            "Dusty disc around magnetosphere: ","(a,1l,1x,a)", .false., ok, .false.)

       if (ttauriDisc) then
          call getReal("rinner", rInner, ttauriRstar/1.e10, cLine, fLine, nLines, &
               "Inner Radius (stellar radii): ","(a,f7.3,a)", 12., ok, .true.)
          call getReal("router", rOuter, real(autocm/1.e10), cLine, fLine, nLines, &
               "Outer Radius (AU): ","(a,f5.1,a)", 20., ok, .true.)
          
          call getReal("rsub", rSublimation, ttauriRstar/1.e10, cLine, fLine, nLines, &
               "Sublimation radius (rstar): ","(a,f5.1,a)", 20., ok, .true.)

          call getReal("height", height, real(autocm/1e10), cLine, fLine, nLines, &
               "Scale height (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

          call getReal("mdisc", mDisc, real(msol), cLine, fLine, nLines, &
               "Disc mass (solar masses): ","(a,f6.4,a)", 1.e-4, ok, .true.)

          call getReal("alphadisc", alphaDisc, 1., cLine, fLine, nLines, &
               "Disc alpha parameter: ","(a,f5.3,a)", 2.25, ok, .true.)

          call getReal("betadisc", betaDisc, 1., cLine, fLine, nLines, &
               "Disc beta parameter: ","(a,f5.3,a)", 1.25, ok, .true.)
       endif

       call getReal("curtainsphi1s", curtainsPhi1s, 1., cLine, fLine, nLines, &
            "Curtains 1: Phi start: (degrees): ","(a,f7.1,1x,a)", 30.0, ok, .false.)
       call getReal("curtainsphi1e", curtainsPhi1e, 1., cLine, fLine, nLines, &
            "Curtains 1: Phi end: (degrees): ","(a,f7.1,1x,a)", 150.0, ok, .false.)
       call getReal("curtainsphi2s", curtainsPhi2s, 1., cLine, fLine, nLines, &
            "Curtains 2: Phi start: (degrees): ","(a,f7.1,1x,a)", 210.0, ok, .false.)
       call getReal("curtainsphi2e", curtainsPhi2e, 1., cLine, fLine, nLines, &
            "Curtains 2: Phi end: (degrees): ","(a,f7.1,1x,a)", 330.0, ok, .false.)

       !  converting the angles in radians  (RK) 
       curtainsPhi1s =    curtainsPhi1s * (pi/180.0) 
       curtainsPhi1e =    curtainsPhi1e * (pi/180.0) 
       curtainsPhi2s =    curtainsPhi2s * (pi/180.0) 
       curtainsPhi2e =    curtainsPhi2e * (pi/180.0) 


       ! The following two are used for "constantcurtans" geometry  (RK)
       call getInteger("curtain_number", curtain_number, cLine, fLine, nLines, &
            "Number of curtains : ","(a,i8,a)", 2, ok, .false.)
       call getReal("curtain_width", curtain_width, 1., cLine, fLine, nLines, &
            "Width of each curtain (degree) : ","(a,f7.1,1x,a)", 120.0, ok, .false.)
       ! converting the curtain width from degrees to radians.
       curtain_width =  curtain_width*Pi/180.0


       call getString("mdottype", mDotType, cLine, fLine, nLines, &
            "T Tauri accretion rate model: ","(a,a,1x,a)","constant", ok, .true.)
       call getReal("mdotpar1", MdotParameter1, 1., cLine, fLine, nLines, &
            "1st parameter for accretion rate: ", "(a,e9.3,1x,a)", 1.0, ok, .true.)
       call getReal("mdotpar2", MdotParameter2, 1., cLine, fLine, nLines, &
            "2nd parameter for accretion rate: ", "(a,e9.3,1x,a)", 1.0, ok, .false.)
       call getReal("mdotpar3", MdotParameter3, 1., cLine, fLine, nLines, &
            "3rd parameter for accretion rate: ", "(a,e9.3,1x,a)", 1.0, ok, .false.)
       call getReal("mdotpar4", MdotParameter4, 1., cLine, fLine, nLines, &
            "4th parameter for accretion rate: ", "(a,e9.3,1x,a)", 1.0, ok, .false.)
       call getReal("mdotpar5", MdotParameter5, 1., cLine, fLine, nLines, &
            "5th parameter for accretion rate: ", "(a,e9.3,1x,a)", 1.0, ok, .false.)
       call getReal("mdotpar6", MdotParameter6, 1., cLine, fLine, nLines, &
            "6th parameter for accretion rate: ", "(a,e9.3,1x,a)", 1.0, ok, .false.)

       call getLogical("ttauriwarp", ttauriwarp, cLine, fLine, nLines, &
            "Include warped disc around magnetosphere: ","(a,1l,1x,a)", .false., ok, .false.)

       call getDouble("hoverr", hoverr, 1.d0, cLine, fLine, nLines, &
            "Warped disc H/R: ","(a,f7.4,1x,a)", 0.3d0, ok, .false.)

       call getLogical("lte", lte, cLine, fLine, nLines, &
            "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("lycontthick", LyContThick, cLine, fLine, nLines, &
            "Optically thick Lyman continuum:","(a,1l,1x,a)", .false., ok, .false.)
       call getString("contflux", contFluxFile, cLine, fLine, nLines, &
            "Continuum flux filename (primary): ","(a,a,1x,a)","none", ok, .true.)
       call getString("popfile", popFilename, cLine, fLine, nLines, &
            "Grid populations filename: ","(a,a,1x,a)","none", ok, .true.)
       call getLogical("writepops", writePops, cLine, fLine, nLines, &
            "Write populations file: ","(a,1l,1x,a)", .true., ok, .true.)
       call getLogical("readpops", readPops, cLine, fLine, nLines, &
            "Read populations file: ","(a,1l,1x,a)", .true., ok, .true.)
       call getLogical("writephasepops", writePhasePops, cLine, fLine, nLines, &
            "Write populations file at each phase: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("readphasepops", readPhasePops, cLine, fLine, nLines, &
            "Read populations file (specific phase): ","(a,1l,1x,a)", .false., ok, .false.)
       !call getLogical("curtains", curtains, cLine, fLine, nLines, &
       !         "Curtains of accretion: ","(a,1l,1x,a)", .false., ok, .false.)
       call getReal("dipoleoffset", dipoleOffset, real(degtorad), cLine, fLine, nLines, &
            "Dipole offset (degrees): ","(a,f7.1,1x,a)", 1.e14, ok, .true.)
       call getLogical("enhance", enhance, cLine, fLine, nLines, &
            "Accretion enhancement: ","(a,1l,1x,a)", .false., ok, .false.)
       call getInteger("nlower", nLower, cLine, fLine, nLines,"Lower level: ","(a,i2,a)",2,ok,.true.)
       call getInteger("nupper", nUpper, cLine, fLine, nLines,"Upper level: ","(a,i2,a)",3,ok,.true.)
       call getLogical("usehartmanntemp", useHartmannTemp, cLine, fLine, nLines, &
            "Use temperatures from Hartmann paper:","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("isotherm", isoTherm, cLine, fLine, nLines, &
            "Use isothermal temperature :","(a,1l,1x,a)", .false., ok, .false.)
       call getReal("isothermtemp", isoThermTemp, 1., cLine, fLine, nLines, &
            "Isothermal temperature (K): ","(a,f7.1,1x,a)", 6500.0, ok, .false.)

       call getLogical("ttauriwind", ttauriWind, cLine, fLine, nLines, &
            "T Tauri disc wind present:","(a,1l,1x,a)", .false., ok, .false.)

       if (ttauriwind) then
          call getDouble("DW_Rmin", DW_Rmin, 1.d0, cLine, fLine, nLines, &
               "Disc wind:: Inner radius of the disc wind [magnetospheric radii]: ", &
               "(a,es9.3,1x,a)", 70.0d0, ok, .true.) 
          call getDouble("DW_Rmax", DW_Rmax, 1.d0, cLine, fLine, nLines, &
               "Disc wind:: Outer radius of the disc [disc wind inner radii]: ", &
               "(a,es9.3,1x,a)", 700.0d0, ok, .true.) 
          call getDouble("DW_Mdot", DW_Mdot, 1.d0,  cLine, fLine, nLines, &
               "Disc wind:: Total mass-loss rate from disc [mass accretion rate]: ", &
               "(a,es9.3,1x,a)", 1.0d-8, ok, .true.) 
          call getDouble("DW_theta", DW_theta, 1.d0, cLine, fLine, nLines, &
               "Disc wind:: Disc wind angle [degrees]: ", &
               "(a,es9.3,1x,a)", 0.5d0, ok, .true.) 
          call getDouble("DW_Twind", DW_temperature, 1.d0, cLine, fLine, nLines, &
               "Disc wind:: Isotherma temperature of disc wind [K]: ", &
               "(a,es9.3,1x,a)", 5000.0d0, ok, .true.) 
          DW_rMin = DW_rmin * ttauriRouter/1.d10
          DW_rMax = DW_rmax * DW_rMin
          DW_theta = 60.d0 * degtoRad
       endif


       if (useHartmannTemp .and. isoTherm) then 
          if (writeoutput)  write(*,'(a)') "WARNING: useHartmannTemp and isoTherm both specified!"
          stop
       end if
       if (useHartmannTemp) &
            call getReal("maxharttemp", maxHartTemp, 1., cLine, fLine, nLines, &
            "Maximum of Hartmann temperature: ","(a,f7.1,1x,a)", 7436., ok, .false.)
       ! sub options for ttauri geometry
       if (ttau_discwind_on) then   ! commnted out here to make ttaur_turn_off_discwind to work
          ! --- parameters for ttauri wind
          call getDouble("DW_d", DW_d, 1.d0, cLine, fLine, nLines, &
               "Disc wind:: Wind soudce displacement [10^10cm]: ", &
               "(a,es9.3,1x,a)", 70.0d0, ok, .true.) 
          call getDouble("DW_Rmin", DW_Rmin,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Inner radius of the disc [10^10cm]: ", &
               "(a,es9.3,1x,a)", 70.0d0, ok, .true.) 
          call getDouble("DW_Rmax", DW_Rmax,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Outer radius of the disc [10^10cm]: ", &
               "(a,es9.3,1x,a)", 700.0d0, ok, .true.) 
          call getDouble("DW_Tmax", DW_Tmax,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Temperature of disc at inner radius [K]: ", &
               "(a,es9.3,1x,a)", 2000.0d0, ok, .true.) 
          call getDouble("DW_gamma", DW_gamma,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Exponent in the disc temperature power law [-]: ", &
               "(a,es9.3,1x,a)", -0.5d0, ok, .true.) 
          call getDouble("DW_Mdot", DW_Mdot,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Total mass-loss rate from disc [Msun/yr]: ", &
               "(a,es9.3,1x,a)", 1.0d-8, ok, .true.) 
          call getDouble("DW_alpha", DW_alpha,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Exponent in the mass-loss rate per unit area [-]: ", &
               "(a,es9.3,1x,a)", 0.5d0, ok, .true.) 
          call getDouble("DW_beta", DW_beta,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Exponent in the modefied beta-velocity law [-]: ", &
               "(a,es9.3,1x,a)", 0.5d0, ok, .true.) 
          call getDouble("DW_Rs", DW_Rs,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Effective accerelation length [10^10cm]: ", &
               "(a,es9.3,1x,a)", 50.0d0*DW_Rmin, ok, .true.) 
          call getDouble("DW_f", DW_f,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Scaling on the terminal velocity [-]: ", &
               "(a,es9.3,1x,a)", 2.0d0, ok, .true.) 
          call getDouble("DW_Twind", DW_Twind,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Isotherma temperature of disc wind [K]: ", &
               "(a,es9.3,1x,a)", 5000.0d0, ok, .true.) 
       endif
       if (ttau_jet_on) then  ! commented out here to make ttaur_turn_off_jet to work
          ! --- parameters for ttauri wind
          call getDouble("JET_Rmin", JET_Rmin,  1.d0, cLine, fLine, nLines, &
               "Minmium radius of Jet [10^10 cm]: ", &
               "(a,es9.3,1x,a)", TTauriRouter/1.0d10, ok, .false.) 
          call getDouble("JET_theta_j", JET_theta_j,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [deg]  jet opening angle: ", &
               "(a,es9.3,1x,a)", 80.0d0, ok, .true.) 
          JET_theta_j = JET_theta_j * (Pi/180.0)  ! converting [deg] to [radians]

          call getDouble("JET_Mdot", JET_Mdot,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [Msun/yr] mass loss rate in the jets: ", &
               "(a,es9.3,1x,a)", 1.0d-9, ok, .true.) 
          call getDouble("JET_a_param", JET_a_param,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [-] a parameter in density function: ", &
               "(a,es9.3,1x,a)", 0.8d0, ok, .true.) 
          call getDouble("JET_b_param", JET_b_param,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [-] b parameter in density function: ", &
               "(a,es9.3,1x,a)", 2.0d0, ok, .true.) 
          call getDouble("JET_Vbase", JET_Vbase,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [km/s] Base velocity of jets: ", &
               "(a,es9.3,1x,a)", 20.0d0, ok, .true.) 
          call getDouble("JET_Vinf", JET_Vinf,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [km/s] Terminal velocity of jets: ", &
               "(a,es9.3,1x,a)", 200.0d0, ok, .true.) 
          call getDouble("JET_beta", JET_beta,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [-] a parameter in velocity function: ", &
               "(a,es9.3,1x,a)", 0.5d0, ok, .true.) 
          call getDouble("JET_gamma", JET_gamma,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [-] a parameter in velocity function: ", &
               "(a,es9.3,1x,a)", 0.05d0, ok, .true.) 
          call getDouble("JET_T", JET_T,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [K]  Isothermal temperature of jets: ", &
               "(a,es9.3,1x,a)", 1.0d4, ok, .true.) 
       endif


       case("fogel")
       call getString("asciifile", textFilename, cLine, fLine, nLines, &
            "Ascii file for abundance data: ","(a,a,a)","none", ok, .true.)
       case("clumpyagb")
          call getReal("rinner", rinner, real(rsol/1.e10), cLine, fLine, nLines, &
               "Inner radius (solar radii): ","(a,1pe8.1,1x,a)", 1000., ok, .true.) 

          call getReal("router", router, rinner, cLine, fLine, nLines, &
               "Outer radius (solar radii): ","(a,1pe8.1,1x,a)", 1000., ok, .true.) 

          call getReal("vterm", vterm, 1.e5, cLine, fLine, nLines, &
               "Terminal velocity (km/s): ","(a,1pe8.1,1x,a)", 1000., ok, .true.) 

          call getReal("mdot", mdot, real(mSol) /( 365.25 * 24. * 3600.),  cLine, fLine, nLines, &
               "Mass-loss rate (solar masses per year): ","(a,1pe8.1,1x,a)", 1000., ok, .true.) 

       case("lexington")
          call getReal("rinner", rInner, 1.e7, cLine, fLine, nLines, &
               "Inner Radius (10^17cm): ","(a,1pe8.1,1x,a)", 30., ok, .true.)

       case("benchmark")
          call getReal("radius1", rCore, real(rsol/1.e10), cLine, fLine, nLines, &
               "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)

          call getReal("rinner", rInner, real(autocm/1.e10), cLine, fLine, nLines, &
               "Inner Radius (AU): ","(a,f5.1,a)", 12., ok, .true.)

          call getReal("router", rOuter, real(autocm/1.e10), cLine, fLine, nLines, &
               "Outer Radius (AU): ","(a,f8.2,a)", 20., ok, .true.)

          call getReal("height", height, real(autocm/1.e10), cLine, fLine, nLines, &
               "Scale height (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

          call getReal("rho", rho, 1., cLine, fLine, nLines, &
               "Density: ","(a,e12.5,a)", 1., ok, .true.)

       case("molebench")
          call getReal("rinner", rInner, 1., cLine, fLine, nLines, &
               "Inner Radius for dumpresults (10^10cm): ","(a,1pe8.2,a)", 1e4, ok, .true.)
          
          call getReal("router", rOuter, 1., cLine, fLine, nLines, &
               "Outer Radius (10^10cm): ","(a,1pe8.2,a)", 1e6, ok, .true.)


       case("molcluster")
          call getString("sphdatafilename", sphdatafilename, cLine, fLine, nLines, &
               "Input sph data file: ","(a,a,1x,a)","sph.dat.ascii", ok, .true.)

          call getReal("hcritPercentile", hcritPercentile, 1., cLine, fLine, nLines, &
               "Percentile for hcrit: ", "(a,f10.4,1x,f10.4)", 0.80, ok, .false.)

          call getReal("hmaxPercentile", hmaxPercentile, 1., cLine, fLine, nLines, &
               "Percentile for hmax: ", "(a,f10.4,1x,f10.4)", 0.99, ok, .false.)

          call getReal("sphNormLimit", sph_norm_limit, 1., cLine, fLine, nLines, &
               "Limit for SPH normalisation: ", "(a,f10.4,1x,f10.4)", 0.3, ok, .false.)

          call getInteger("kerneltype", kerneltype, cLine, fLine, nLines, &
               "Kernel type (0 is exponential/1 is spline): ","(a,i1,a)",0, ok, .false.)

    case("shakara")

       oneKappa = .true.
       fastIntegrate = .true.
       call getLogical("gasopacity", includeGasOpacity, cLine, fLine, nLines, &
            "Include gas opacity: ","(a,1l,a)", .false., ok, .false.)

       call getReal("mdot", mdot,  real(msol * secstoyears), cLine, fLine, nLines, &
            "Mass accretion  rate (msol/yr): ","(a,1p,e12.5,a)", 0.0,  ok, .false.)

       call getLogical("noscat", noScattering, cLine, fLine, nLines, &
            "No scattering opacity in model: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("smoothinneredge", smoothInnerEdge, cLine, fLine, nLines, &
            "Smooth density drop at inner edge: ","(a,1l,1x,a)", .false., ok, .false.)

       call getReal("radius1", rCore, real(rsol/1.e10), cLine, fLine, nLines, &
            "Core radius (solar radii): ","(a,f7.3,a)", 10., ok, .true.)

       call getReal("rinner", rInner, rCore, cLine, fLine, nLines, &
            "Inner Radius (stellar radii): ","(a,f7.3,a)", 12., ok, .true.)

       call getReal("router", rOuter, real(autocm/1.d10), cLine, fLine, nLines, &
            "Outer Radius (AU): ","(a,f5.1,a)", 20., ok, .true.)

       call getReal("height", height, real(autocm/1.d10), cLine, fLine, nLines, &
            "Scale height (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

       call getReal("mass1", mCore, real(msol), cLine, fLine, nLines, &
            "Core mass (solar masses): ","(a,f6.4,a)", 0.5, ok, .true.)

       call getReal("mdisc", mDisc, real(msol), cLine, fLine, nLines, &
            "Disc mass (solar masses): ","(a,f6.4,a)", 1.e-4, ok, .true.)

       call getReal("alphadisc", alphaDisc, 1., cLine, fLine, nLines, &
            "Disc alpha parameter: ","(a,f5.3,a)", 2.25, ok, .true.)

       call getReal("betadisc", betaDisc, 1., cLine, fLine, nLines, &
            "Disc beta parameter: ","(a,f5.3,a)", 1.25, ok, .true.)

       call getLogical("vardustsub", variableDustSublimation, cLine, fLine, nLines, &
            "Variable dust sublimation temperature: ", "(a,1l,1x,a)", .false., ok, .true.)

       call getLogical("hydro", solveVerticalHydro, cLine, fLine, nLines, &
            "Solve vertical hydrostatical equilibrium: ","(a,1l,1x,a)", .false., ok, .false.)

       if (solveVerticalHydro) then
          call getInteger("nhydro", nhydro,  cline, fLine, nLines, &
               "Max number of hydro iterations : ","(a,i4,a)", 5, ok, .true.)
       endif


       call getReal("heightsplitfac", heightSplitFac, 1., cLine, fLine, nLines, &
            "Splitting factor for scale height (local scale heights): ","(a,f5.2,a)", 0.2, ok, .false.)

!       rCore = rCore * rSol / 1.e10
!       rinner = (rinner * (rCore * 1e10)) / autocm
!
!       rho0 = densityfrommass(mdisc, height, rinner, router, 100.0, alphadisc, betadisc)
!
!       rInner = rInner * autocm / 1.e10
!       rOuter = rOuter * autoCm / 1.e10
!       height = height * autoCm / 1.e10
!       mCore = mCore * mSol
!       mDisc = mDisc * mSol

       rho0  = mDisc *(betaDisc-alphaDisc+2.) / ( twoPi**1.5 * (height*1.e10)/(100.d0*autocm)**betaDisc  &
            * (rInner*1.d10)**alphaDisc * &
            (((rOuter*1.e10)**(betaDisc-alphaDisc+2.)-(rInner*1.e10)**(betaDisc-alphaDisc+2.))) )

       case("theGalaxy")
          call getString("sphdatafilename", sphdatafilename, cLine, fLine, nLines, &
               "Input sph data file: ","(a,a,1x,a)","sph.dat.ascii", ok, .true.)
          
          call getString("inputFileFormat", inputFileFormat, cLine, fLine, nLines, &
               "Input file format: ","(a,a,1x,a)","binary", ok, .false.)
          
          call getReal("hcritPercentile", hcritPercentile, 1., cLine, fLine, nLines, &
               "Percentile for hcrit: ", "(a,f10.4,1x,f10.4)", 0.80, ok, .false.)

          call getReal("hmaxPercentile", hmaxPercentile, 1., cLine, fLine, nLines, &
               "Percentile for hmax: ", "(a,f10.4,1x,f10.4)", 0.99, ok, .false.)
          
          call getReal("sphNormLimit", sph_norm_limit, 1., cLine, fLine, nLines, &
               "Limit for SPH normalisation: ", "(a,f10.4,1x,f10.4)", 0.3, ok, .false.)

          call getInteger("kerneltype", kerneltype, cLine, fLine, nLines, &
               "Kernel type (0 is exponential/1 is spline): ","(a,i1,a)",0, ok, .false.)

          call getLogical("internalView", internalView, cLine, fLine, nLines, &
               "View as our Galaxy:", "(a,1l,1x,a)", .false., ok, .true.)

 ! Read parameters used by galactic plane survey 
          if ( internalView ) then 
             call getDouble("intPosX", intPosX,  1.0_db, cLine, fLine, nLines, "Observer x position (x10^10cm)", &
                  "(a,e10.4,1x,a)", 0.d0, ok, .false.)
             call getDouble("intPosY", intPosY,  1.0_db, cLine, fLine, nLines, "Observer y position (x10^10cm)", &
                  "(a,e10.4,1x,a)", 2.2e12_db, ok, .false.)
             call getDouble("intPosZ", intPosZ, 1.0_db,  cLine, fLine, nLines, "Observer z position (x10^10cm)", &
                  "(a,e10.4,1x,a)", 0.d0, ok, .false.)

             call getDouble("intDeltaVx", intDeltaVx, 1.0_db,  cLine, fLine, nLines, "Observer x velocity boost (km/s)", &
                  "(a,f8.2x,a)", 0.d0, ok, .false.)
             call getDouble("intDeltaVy", intDeltaVy, 1.0_db,  cLine, fLine, nLines, "Observer y velocity boost (km/s)", &
                  "(a,f8.2x,a)", 0.d0, ok, .false.)
             call getDouble("intDeltaVz", intDeltaVz, 1.0_db,  cLine, fLine, nLines, "Observer z velocity boost (km/s)", &
                  "(a,f8.2x,a)", 0.d0, ok, .false.)

             ! For the internal case use these parameters to rotate the galaxy so we are not looking along cell boundaries. 
             ! Rotation about y-axis
             call getDouble("galaxyInclination", galaxyInclination, 1.0_db,  cLine, fLine, nLines, &
                  "Galaxy Inclination:", "(a,f4.1,1x,a)", 45.d0, ok, .false.)
             ! Rotation about z-axis
             call getDouble("galaxyPositionAngle", galaxyPositionAngle, 1.0_db, cLine, fLine, nLines, &
                  "Galaxy position angle:", "(a,f4.1,1x,a)", 0.d0, ok, .false.)

             ! Rotate the centre of the grid so it covers the required domain
             gridCentre     = VECTOR(amrGridCentreX, amrGridCentreY, amrGridCentreZ)
             write(message,'(a,3(ES12.3,2x),a)') "Grid centre is ", gridCentre
             call writeInfo(message)
             gridCentre     = rotateZ( gridCentre, galaxyPositionAngle*degToRad )
             gridCentre     = rotateY( gridCentre, galaxyInclination*degToRad   )
             write(message,'(a,3(ES12.3,2x),a)') "Modified grid centre is ", gridCentre
             call writeInfo(message)
             amrGridCentreX = gridCentre%x
             amrGridCentreY = gridCentre%y
             amrGridCentreZ = gridCentre%z

          else
             call getDouble("galaxyInclination", galaxyInclination, 1.0_db, cLine, fLine, nLines, &
                  "Galaxy Inclination:", "(a,f4.1,1x,a)", 50.d0, ok, .false.)
             call getDouble("galaxyPositionAngle", galaxyPositionAngle, 1.0_db, cLine, fLine, nLines, &
                  "Galaxy position angle:", "(a,f4.1,1x,a)", 20.d0, ok, .false.)
             rotateViewAboutX = 90.0 - galaxyPositionAngle 
             rotateViewAboutY = 90.0 + galaxyPositionAngle
             rotateViewAboutZ = 0.0
             call getDouble("dataCubeVelocityOffset", dataCubeVelocityOffset, 1.0_db, cLine, fLine, nLines, &
                  "Data cube velocity offset:", "(a,f8.1,1x,a)", 0.d0, ok, .true.)
          end if

    end select
  end subroutine readGeometrySpecificParameters
         

  subroutine readGridInitParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok


    gridUsesAMR = .true.

    call getLogical("cylindrical", cylindrical, cLine, fLine, nLines, &
         "Grid uses 3D cylindical  coords: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("amr1d", amr1d, cLine, fLine, nLines, &
         "AMR grid is in one-dimensions only: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("amr2d", amr2d, cLine, fLine, nLines, &
         "AMR grid is in two-dimensions only: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("amr3d", amr3d, cLine, fLine, nLines, &
         "AMR grid is in three-dimensions: ","(a,1l,1x,a)", .false., ok, .false.)

    if (.not.(amr1d.or.amr2d.or.amr3d)) then
       call writeFatal("AMR dimensionality must be specified")
       stop
    endif
    if (amr1d.and.(amr2d.or.amr3d)) &
       call writeFatal("Only one of amr1d, amr2d or amr3d should be set to true")

    if (amr2d.and.(amr1d.or.amr3d)) &
       call writeFatal("Only one of amr1d, amr2d or amr3d should be set to true")

    if (amr3d.and.(amr1d.or.amr2d)) &
       call writeFatal("Only one of amr1d, amr2d or amr3d should be set to true")


    call getInteger("mindepthamr", minDepthAMR, cLine, fLine, nLines, "Minimum cell depth of AMR grid: ", &
         & "(a,i3,a)",5,ok,.false.)

    call getInteger("maxdepthamr", maxDepthAMR, cLine, fLine, nLines, "Maximum cell depth of AMR grid: ", &
         & "(a,i3,a)",31,ok,.false.)

    call getReal("amrgridsize", amrGridSize, 1., cLine, fLine, nLines, &
         "Size of adaptive mesh grid: ","(a,1pe8.1,1x,a)", 1000., ok, .true.) 
    call getDouble("amrgridcentrex", amrGridCentreX, 1.d0 , cLine, fLine, nLines, &
         "Grid centre X-coordinate: ","(a,es9.3,1x,a)", 0.0d0, ok, .false.) 
    call getDouble("amrgridcentrey", amrGridCentreY, 1.d0 , cLine, fLine, nLines, &
         "Grid centre Y-coordinate: ","(a,es9.3,1x,a)", 0.0d0, ok, .false.) 
    call getDouble("amrgridcentrez", amrGridCentreZ, 1.d0, cLine, fLine, nLines, &
         "Grid centre Z-coordinate: ","(a,es9.3,1x,a)", 0.0d0, ok, .false.) 

    if (amr2d) amrGridCentrex = amrGridSize/2.d0


    call getDouble("limitscalar", limitScalar, 1.d0, cLine, fLine, nLines, &
         "Scalar limit for subcell division: ","(a,es9.3,1x,a)", 1000._db, ok, .false.) 

    call getDouble("limittwo", limitScalar2, 1.d0, cLine, fLine, nLines, &
         "Second scalar limit for subcell division: ","(a,es9.3,1x,a)", 0._db, ok, .false.) 

    call getLogical("doVelocitySplit", doVelocitySplit, cLine, fLine, nLines, &
         "Use velocity splitting condition for SPH to grid: ","(a,1l,1x,a)",.false., ok, .false.)

    call getString("geometry", geometry, cLine, fLine, nLines, &
         "Geometry: ","(a,a,a)","sphere",ok, .true.)

  end subroutine readGridInitParameters

  subroutine readDustPhysicsParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok
    real :: grainFracTotal
    integer :: i
    character(len=20) :: grainTypeLabel, grainFracLabel, aMinLabel, &
         aMaxLabel, qDistLabel, pdistLabel, a0label

       oneKappa = .true.

       call getReal("dusttogas", dusttoGas, 1., cLine, fLine, nLines, &
            "Dust to gas ratio: ","(a,f5.3,a)",0.01,ok,.false.)

       call getInteger("ndusttype", nDustType, cLine, fLine, nLines,"Number of different dust types: ","(a,i12,a)",1,ok,.false.)
       if (nDustType .gt. maxDustTypes) then
          if (writeoutput) write (*,*) "Max dust types exceeded: ", maxDustTypes
          stop
       end if


       grainFracTotal = 0.
       do i = 1, nDustType
          write(grainTypeLabel, '(a,i1.1)') "graintype",i
          write(grainFracLabel, '(a,i1.1)') "grainfrac",i
          write(aMinLabel, '(a,i1.1)') "amin",i
          write(aMaxLabel, '(a,i1.1)') "amax",i
          write(qDistLabel, '(a,i1.1)') "qdist",i
          write(pDistLabel, '(a,i1.1)') "pdist",i
          write(a0Label, '(a,i1.1)') "a0",i
          !       if (writeoutput) write(*,'(a,i1.1)') "Dust properties for grain ",i
          !       if (writeoutput) write(*,'(a,i1.1)') "-------------------------------"
          !       if (writeoutput) write(*,*)
          call getString(grainTypeLabel, grainType(i), cLine, fLine, nLines, &
               "Grain type: ","(a,a,1x,a)","sil_dl", ok, .true.)

          call getReal(grainFracLabel, grainFrac(i), 1., cLine, fLine, nLines, &
               "Grain fractional abundance: ","(a,f8.5,1x,a)",1. , ok, .false.)
          grainFracTotal = grainFracTotal + grainFrac(i)

          call getReal(aminLabel, aMin(i), 1., cLine, fLine, nLines, &
               "Min grain size (microns): ","(a,f8.5,1x,a)", 0.005, ok,  .true.)

          call getReal(amaxLabel, aMax(i), 1., cLine, fLine, nLines, &
               "Max grain size (microns): ","(a,f10.5,1x,a)", 0.25, ok, .true.)

          call getReal(qDistLabel, qdist(i), 1., cLine, fLine, nLines, &
               "Grain power law: ","(a,f4.1,1x,a)", 3.5, ok, .true. )

          call getReal(a0Label, a0(i), 1., cLine, fLine, nLines, &
               "Scale length of grain size (microns): ","(a,f8.5,1x,a)", 1.0e20, ok, .false.)


          call getReal(pdistLabel, pdist(i), 1., cLine, fLine, nLines, &
               "Exponent for exponential cut off: ","(a,f4.1,1x,a)", 1.0, ok, .false. )
          if (writeoutput) write(*,*)
       enddo

       call getLogical("iso_scatter", isotropicScattering, cLine, fLine, nLines, &
         "Isotropic scattering: ","(a,1l,1x,a)", .false., ok, .false.)


  end subroutine readDustPhysicsParameters

  subroutine readAtomicPhysicsParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    integer :: i 
    character(len=20) :: keyword
    logical :: ok

    call writeBanner("Atomic physics data","#",TRIVIAL)

    call getInteger("natom", nAtom, cLine, fLine, nLines, &
         "Number of model atoms to solve for: ","(a,i1,a)",1,ok,.true.)
    do i = 1, nAtom
       write(keyword, '(a,i1)') "atom",i
       call getString(keyword, atomFileName(i), cLine, fLine, nLines, &
            "Use atom filename: ","(a,a,1x,a)","none", ok, .true.)
    enddo
    call getInteger("itrans", itransLine, cLine, fLine, nLines, &
         "Index of line transition: ","(a,i4,a)",4,ok,.true.)

    call getInteger("iatom", itransAtom, cLine, fLine, nLines, &
         "Atom of line transition: ","(a,i4,a)",4,ok,.true.)



  end subroutine readAtomicPhysicsParameters

  subroutine readSourceParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:), message
    logical :: fLine(:)
    integer :: nLines
    integer :: i 
    character(len=20) :: keyword
    logical :: ok


    call writeBanner("Stellar source  data","#",TRIVIAL)

    call getInteger("nsource", inputNSource, cLine, fLine, nLines, &
         "Number of sources: ","(a,i2,a)",1,ok,.true.)
    do i = 1, inputnSource
       if (writeoutput) write(*,*) " "
       write(message,'(a,i1)') "Source number: ",i
       call writeInfo(message,TRIVIAL)
       write(keyword, '(a,i1)') "radius",i
       call getDouble(keyword, sourceRadius(i), rsol/1.d10, cLine, fLine, nLines, &
            "Source radius (solar radii) : ","(a,f7.2,a)",1.d0, ok, .true.)

       write(keyword, '(a,i1)') "teff",i
       call getDouble(keyword, sourceTeff(i), 1.d0, cLine, fLine, nLines, &
            "Source temperature (K) : ","(a,f8.0,a)",1.d0, ok, .true.)

       write(keyword, '(a,i1)') "mass",i
       call getDouble(keyword, sourceMass(i), mSol, cLine, fLine, nLines, &
            "Source mass (solar masses) : ","(a,f7.2,a)",1.d0, ok, .true.)

       write(keyword, '(a,i1)') "contflux",i
       call getString(keyword, inputcontFluxFile(i), cLine, fLine, nLines, &
            "Continuum flux file: ","(a,a,a)","none", ok, .true.)

       write(keyword, '(a,i1)') "sourcepos",i
       call getVector(keyword, sourcePos(i), 1.d0, cLine, fLine, nLines, &
            "Source position (10^10 cm): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)
    enddo

  end subroutine readSourceParameters

  subroutine readMolecularPhysicsParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

       call getLogical("lte", lte, cLine, fLine, nLines, &
            "Read in LTE grid: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("densitysubsample", densitysubsample, cLine, fLine, nLines, &
               "Use density interpolation: ","(a,1l,a)", .false., ok, .false.)
       call getLogical("suppresswarnings", suppressWarnings, cLine, fLine, nLines, &
            "Suppress Warnings: ","(a,l,1x,a)",.false., ok, .false.) 
       call getLogical("addnewmoldata", addnewmoldata, cLine, fLine, nLines, &
            "Add new molecular data to non-molecular grid: ","(a,l,1x,a)",.false., ok, .false.)
       call getString("moleculefile", moleculefile, cLine, fLine, nLines, &
            "Input molecule filename: ","(a,a,1x,a)","none", ok, .true.)
       call getReal("distance", gridDistance, real(pctocm), cLine, fLine, nLines, &
            "Grid distance (pc): ","(a,f6.1,1x,a)", 1., ok, .true.)
       call getInteger("initnray", initnray, cLine, fLine, nLines, &
               "Number of fixed rays for stage 1: ","(a,i4,a)", 100, ok, .false.)
       call getLogical("dongstep", dongstep, cLine, fLine, nLines, &
               "Use Ng Acceleration: ","(a,1l,a)", .false., ok, .false.)
       call getLogical("quasi", quasi, cLine, fLine, nLines, &
               "Use Quasirandom numbers: ","(a,1l,a)", .false., ok, .false.)
       call getReal("tolerance", tolerance, 1., cLine, fLine, nLines, &
            "Maximum Fractional Change in level populations:","(a,f4.1,1x,a)", 0.01, ok, .false.)
       call getReal("vturb", vturb, 1., cLine, fLine, nLines, &
            "Subsonic turbulent velocity (km/s):","(a,f4.1,1x,a)", 0.3, ok, .false.)
       call getLogical("noturb", noturb, cLine, fLine, nLines, &
            "No microturbulence","(a,1l,a)",.false., ok, .false.)
       call getInteger("setmaxlevel", setmaxlevel, cLine, fLine, nLines, &
            "Maximum molecular level to be considered:","(a,i2,1x,a)", 0, ok, .false.)
       call getReal("molAbundance", molAbundance, 1., cLine, fLine, nLines, &
            "Molecular Abundance:","(a,e12.5,1x,a)", 1e-9, ok, .false.)
       call getLogical("isinlte", isinlte, cLine, fLine, nLines, &
            "Assume LTE: ", "(a,1l,1x,a)", .false., ok, .false.)
       call getReal("dusttogas", dusttoGas, 1., cLine, fLine, nLines, &
            "Dust to gas ratio: ","(a,f5.3,a)",0.01,ok,.false.)
       call getLogical("plotlevels", plotlevels, cLine, fLine, nLines, &
            "Plot Molecular Levels ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("realdust", realdust, cLine, fLine, nLines, &
            "Use realistic dust model: ", "(a,1l,1x,a)", .true., ok, .false.)

       call getLogical("outputconvergence", outputconvergence, cLine, fLine, nLines, &
            "Write out convergence data : ", "(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("dotune", dotune, cLine, fLine, nLines, &
            "Write out convergence data : ", "(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("blockhandout", blockHandout, cLine, fLine, nLines, &
            "Use blockhandout for parallel computations ", "(a,1l,1x,a)", .true., ok, .false.)


       call getLogical("doCOchemistry", doCOchemistry, cLine, fLine, nLines, &
            "Use drop profile to model CO depletion: ", "(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("useDust", useDust, cLine, fLine, nLines, &
            "Calculate continuum emission from dust:", "(a,1l,1x,a)", .false., ok, .false.)

       if(doCOchemistry) then

          call getReal("fracCOdepletion", x_D, 1., cLine, fLine, nLines, &
               "Fraction of CO depletion", "(a,1l,1x,a)", 0.1, ok, .true.)
          
       endif




  end subroutine readMolecularPhysicsParameters


  subroutine readPhotoionPhysicsParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call getLogical("quickthermal", quickThermal, cLine, fLine, nLines, &
         "Compute photoionization equilibrium: ","(a,1l,a)", .false., ok, .false.)

    call getLogical("usemetals", usemetals, cLine, fLine, nLines, &
         "Use metals in photoionization calculation: ","(a,1l,a)", .true., ok, .false.)


    call getReal("h_abund", h_abund, 1., cLine, fLine, nLines, &
         "Hydrogen abdunance: ","(a,1PF8.3,a)", &
         1., ok, .false.)

    call getReal("he_abund", he_abund, 1., cLine, fLine, nLines, &
         "Helium abdunance: ","(a,1PF8.3,a)", &
         0.1, ok, .false.)

    call getReal("c_abund", c_abund, 1., cLine, fLine, nLines, &
         "Carbon abdunance: ","(a,1PF8.3,a)", &
         22.e-5, ok, .false.)

    call getReal("n_abund", n_abund, 1., cLine, fLine, nLines, &
         "Nitrogen abdunance: ","(a,1PF8.3,a)", &
         4.e-5, ok, .false.)

    call getReal("o_abund", o_abund, 1., cLine, fLine, nLines, &
         "Oxygen abdunance: ","(a,1PF8.3,a)", &
         33.e-5, ok, .false.)

    call getReal("ne_abund", ne_abund, 1., cLine, fLine, nLines, &
         "Neon abdunance: ","(a,1PF8.3,a)", &
         5.e-5, ok, .false.)

    call getReal("s_abund", s_abund, 1., cLine, fLine, nLines, &
         "Sulphur abdunance: ","(a,1PF8.3,a)", &
         0.9e-5, ok, .false.)

  end subroutine readPhotoionPhysicsParameters


  subroutine readMolecularLoopParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call getLogical("hydrovelconv", hydroVelocityConv, cLine, fLine, nLines, &
         "Take velocity data from hydrodynamics: ","(a,1l,1x,a)", .false., ok, .false.)

  end subroutine readMolecularLoopParameters

  subroutine readAtomicLoopParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call getLogical("lte", lte, cLine, fLine, nLines, &
         "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
    call getReal("vturb", vturb, real(kmstoc), cLine, fLine, nLines, &
         "Turbulent velocity (km/s):","(a,f4.1,1x,a)", 50., ok, .true.)


  end subroutine readAtomicLoopParameters

  subroutine readRadiativeEquilibriumParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok
    call getInteger("nlucy", nLucy, cLine, fLine, nLines,"Number of photons per lucy iteration: ","(a,i12,a)",0,ok,.false.)

    call getInteger("iterLucy", iterLucy, cline, fLine, nlines, "Minimum number of Lucy iterations: ", "(a,i3,a)",3,ok,.false.)

    call getLogical("forceLucyConv", forceLucyConv, cLine, fLine, nLines, &
         "Force convergence of Lucy algorithm: ","(a,1l,1x,a)", .false., ok, .false.)

    call getReal("lucy_undersampled", lucy_undersampled, 1., cLine, fLine, nLines, &
         "Minimum percentage of undersampled cell in lucy iteration: ", &
         "(a,f4.2,a)",0.0,ok,.false.)

    call getReal("diffdepth", diffDepth, 1., cLine, fLine, nLines, &
         "Depth of diffusion zone (in Rosseland optical depths): ", &
         "(a,f5.1,a)",10.0,ok,.false.)

    call getInteger("mincrossings", minCrossings, cLine, fLine, nLines, &
         "Minimum crossings required for cell to be sampled: ","(a,i12,a)",100,ok,.false.)

    call getReal("tminglobal", TMinGlobal, 1., cLine, fLine, nLines, &
         "Minimum Temperature (K): ","(a,f4.1,1x,a)", 2.8, ok, .false.)

    call getReal("taudiff", tauDiff, 1., cLine, fLine, nLines, &
         "Mininum optical depth of cell to be in diffusion approx : ","(a,f7.1,a)",100., ok, .false.)

    call getReal("tauforce", tauForce, 1., cLine, fLine, nLines, &
         "Forced optical depth of cell to be in diffusion approx : ","(a,f7.1,a)",10., ok, .false.)

    call getInteger("mincrossings", minCrossings, cLine, fLine, nLines, &
         "Minimum crossings required for cell to be sampled: ","(a,i12,a)",100,ok,.false.)

    call getLogical("dosmoothgrid", doSmoothGrid, cLine, fLine, nLines, &
         "Smooth AMR grid: ","(a,1l,1x,a)", .false., ok, .false.)

    if (doSmoothGrid) then
       call getReal("smoothfactor", smoothFactor, 1.0, cLine, fLine, nLines, &
            "Inter-cell maximum ratio before smooth: ","(a,f6.1,1x,a)", 5., ok, .true.)
    else
       smoothFactor = 0.0      
    end if

    call getLogical("smoothgridtau", doSmoothGridtau, cLine, fLine, nLines, &
         "Smooth AMR grid using tau: ","(a,1l,1x,a)", .false., ok, .false.)

    if (dosmoothgridtau) then
       call getReal("lambdasmooth", lambdasmooth, 1.0, cLine, fLine, nLines, &
            "Lambda for tau smoothing: ","(a,1PE10.3,1x,a)", 5500.0, ok, .false.)
       call getReal("taumax", tauSmoothMax, 1.0, cLine, fLine, nLines, &
            "Maximum tau for smoothing: ","(a,f10.1,1x,a)", 0.5, ok, .false.)
       call getReal("taumin", tauSmoothMin, 1.0, cLine, fLine, nLines, &
            "Minimum tau for smoothing: ","(a,f10.1,1x,a)", 0.001, ok, .false.)
    endif

    call getReal("edenstol", eDensTol, 1., cLine, fLine, nLines, &
         "Fractional change in energy density for convergence: ","(a,f7.1,a)",0.001, ok, .false.) ! used for gauss-seidel sweep also

    call getReal("scatteredlightwavelength", scatteredLightWavelength, 1., cLine, fLine, nLines, &
         "Wavelength of scattered light","(a,f7.1,a)",1.e4, ok, .false.)


  end subroutine readRadiativeEquilibriumParameters

  subroutine readPhotoionEquilibriumParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call getInteger("nmonte", inputnMonte, cLine, fLine, nLines, &
         "Number of photons in image","(a,i8,a)", 0, ok, .false.)

    call getReal("taudiff", tauDiff, 1., cLine, fLine, nLines, &
         "Mininum optical depth of cell to be in diffusion approx : ","(a,f7.1,a)",100., ok, .false.)

  end subroutine readPhotoionEquilibriumParameters

  subroutine readHydrodynamicsParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call getReal("cfl", cflNumber, 1., cLine, fLine, nLines, &
         "Courant number:","(a,f4.1,1x,a)", 0.3, ok, .false.)

    call getDouble("griddistancescale", gridDistanceScale, 1.d0, cLine, fLine, nLines, &
         "Distance grid scale:","(a,e12.3,1x,a)", 1.d10, ok, .false.)

  end subroutine readHydrodynamicsParameters

  subroutine readDataCubeParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call getReal("inclination", thisinclination, real(degtorad), cLine, fLine, nLines, &
         "Inclination angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)
    call getReal("positionangle", positionAngle(1), real(degtorad), cLine, fLine, nLines, &
         "Position angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)
    call getString("datacubefile", datacubeFilename, cLine, fLine, nLines, &
         "Output datacube  filename: ","(a,a,1x,a)","none", ok, .true.)
    call getReal("imageside", imageside, 1., cLine, fLine, nLines, &
         "Image size (x10^10cm):","(a,es9.3,1x,a)", 5e7, ok, .true.)
    call getReal("distance", gridDistance, real(pcToCm), cLine, fLine, nLines, &
         "Grid distance (pc): ","(a,f4.1,1x,a)", 100., ok, .false.)
    call getInteger("npixels", npixels, cLine, fLine, nLines, &
         "Number of pixels per row: ","(a,1x,i4,a)", 50, ok, .true.)
    call getInteger("nv", nv, cLine, fLine, nLines, &
         "Number of velocity bins ","(a,i4,a)", 50, ok, .true.)
    call getDouble("maxVel", maxVel, 1.d0, cLine, fLine, nLines, &
         "Maximum Velocity Channel (km/s): ","(a,f5.1,1x,a)", 1.0d0, ok, .true.)
    call getDouble("minVel", minVel, 1.0_db, cLine, fLine, nLines, &
         "Minimum Velocity Channel (km/s): ","(a,f4.1,1x,a)", -1.0d0*maxVel, ok, .false.)
    call getLogical("h21cm", h21cm, cLine, fLine, nLines, &
         "Calculate data cube of 21cm emission: ","(a,1l,a)", .false., ok, .false.)
    call getLogical("splitCubes", splitCubes, cLine, fLine, nLines, &
         "Split intensity into +/- components: ","(a,1l,a)", .false., ok, .false.)

    if ( h21cm ) then 
       call getInteger("nSubpixels", nSubpixels, cLine, fLine, nLines, &
            "Subpixel splitting (0 denotes adaptive)","(a,i4,a)", 1, ok, .false.)
       call getLogical("densitysubsample", densitysubsample, cLine, fLine, nLines, &
            "Use density interpolation: ","(a,1l,a)", .false., ok, .false.)
       call getDouble("centrevecx", centrevecx, 1.d0, cLine, fLine, nLines, &
            "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .true.)
       call getDouble("centrevecy", centrevecy, 1.d0, cLine, fLine, nLines, &
            "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .true.)
       call getLogical("wanttau", wanttau, cLine, fLine, nLines, &
            "Write Tau information to datacube: ","(a,1l,1x,a)", .false., ok, .false.)
    end if
 
   if (molecularPhysics) then
       call getInteger("nSubpixels", nSubpixels, cLine, fLine, nLines, &
            "Subpixel splitting (0 denotes adaptive)","(a,i4,a)", 1, ok, .false.)
       call getLogical("densitysubsample", densitysubsample, cLine, fLine, nLines, &
            "Use density interpolation: ","(a,1l,a)", .false., ok, .false.)
       call getLogical("lineimage", lineImage, cLine, fLine, nLines, &
            "Line emission: ","(a,1l,a)", .true., ok, .false.)
       if(.not. lineimage) then
          call getReal("lamline", lamLine, 1.e4,cLine, fLine, nLines, &
               "Line emission wavelength (um): ","(a,f6.1,1x,a)", 850., ok, .true.)
       endif
       
       call getLogical("rgbcube", rgbCube, cLine, fLine, nLines, &
         "Create an RGB data cube (reverses velocity axis): ","(a,1l,a)", .false., ok, .false.)
       call getInteger("itrans", itrans, cLine, fLine, nLines, &
            "Molecular Line Transition","(a,i4,a)", 1, ok, .true.)
       call getReal("beamsize", beamsize, 1., cLine, fLine, nLines, &
            "Beam size (arcsec): ","(a,f4.1,1x,a)", 1000., ok, .false.)
       call getDouble("rotateviewaboutx", rotateViewAboutX, 1.d0, cLine, fLine, nLines, &
            "Angle to rotate about X (deg): ","(a,f4.1,1x,a)", 0.d0, ok, .false.)
       call getDouble("rotateviewaboutz", rotateViewAboutZ, 1.d0, cLine, fLine, nLines, &
            "Angle to rotate about Z (deg): ","(a,f4.1,1x,a)", 0.d0, ok, .false.)
       call getDouble("centrevecx", centrevecx, 1.d0, cLine, fLine, nLines, &
            "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .true.)
       call getDouble("centrevecy", centrevecy, 1.d0, cLine, fLine, nLines, &
            "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .true.)
       call getDouble("centrevecz", centrevecz, 1.d0, cLine, fLine, nLines, &
            "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .true.)
       call getLogical("wanttau", wanttau, cLine, fLine, nLines, &
            "Write Tau information to datacube: ","(a,1l,1x,a)", .false., ok, .false.)
    endif
       
  end subroutine readDataCubeParameters

  subroutine readImageParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:), message
    logical :: fLine(:)
    integer :: nLines
    logical :: ok
    integer :: i
    character(len=20) :: keyword


    call getBigInteger("nphotons", nphotons, cLine, fLine, nLines, &
         "Number of photons in image: ","(a,i8,a)", 10000, ok, .true.)

    call getReal("distance", gridDistance, 1., cLine, fLine, nLines, &
         "Grid distance (pc): ","(a,f4.1,1x,a)", 100., ok, .false.)

    call getInteger("nimage", nimage, cLine, fLine, nLines, &
         "Number of images to calculate: ","(a,i8,a)", 1, ok, .false.)

    call getLogical("inarcsec", imageinArcsec, cLine, fLine, nLines, &
         "Write image distances in arcseconds: ","(a,1l,1x,a)", .false., ok, .false.)

    call getReal("imagesize", setImageSize, real(autocm), cLine, fLine, nLines, &
         "Image size (AU): ", "(a,1pe10.2,1x,a)", 0., ok, .false.)


    if (nimage == 1) then

       call getString("imagefile", imageFilename(1), cLine, fLine, nLines, &
            "Output image  filename: ","(a,a,1x,a)","none", ok, .true.)

       call getReal("lambdaimage", lambdaImage(1),1., cLine, fLine, nLines, &
         "Wavelength for monochromatic image (A):","(a,f8.1,1x,a)", 6562.8, ok, .false.)
       if (photoionPhysics) then
          call getString("imagetype", outputimageType(1), cLine, fLine, nLines, &
               "Type of output image: ","(a,a,1x,a)","none", ok, .true.)
       endif
       call getInteger("npixels", npixelsArray(1), cLine, fLine, nLines, &
            "Number of pixels per side in image","(a,i8,a)", 200, ok, .false.)

       call getReal("inclination", inclinations(1), real(degtorad), cLine, fLine, nLines, &
            "Inclination angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)

    call getReal("positionangle", positionAngle(1), real(degtorad), cLine, fLine, nLines, &
         "Position angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)

    else
       do i = 1, nImage

          write(message,'(a,i1.1)') "Details for image: ",i
          call writeInfo(message)
          call writeInfo(" ")

          write(keyword,'(a,i1.1)') "imagefile",i
          call getString(keyword, imageFilename(i), cLine, fLine, nLines, &
               "Output image  filename: ","(a,a,1x,a)","none", ok, .true.)
          write(keyword,'(a,i1.1)') "lambdaimage",i
          call getReal(keyword, lambdaImage(i),1., cLine, fLine, nLines, &
               "Wavelength for monochromatic image (A):","(a,f8.1,1x,a)", 6562.8, ok, .false.)
          if (photoionPhysics) then
             write(keyword,'(a,i1.1)') "imagetype",i
             call getString(keyword, outputimageType(i), cLine, fLine, nLines, &
                  "Type of output image: ","(a,a,1x,a)","none", ok, .true.)
          endif
          write(keyword,'(a,i1.1)') "npixels",i
          call getInteger(keyword, npixelsArray(i), cLine, fLine, nLines, &
               "Number of pixels per side in image","(a,i8,a)", 200, ok, .false.)
          write(keyword,'(a,i1.1)') "inclination",i
          call getReal(keyword, inclinationArray(i), real(degtorad), cLine, fLine, nLines, &
               "Inclination of image: ","(a,f4.1,1x,a)",0., ok, .true.)
          write(keyword,'(a,i1.1)') "positionangle",i
          call getReal(keyword, positionAngle(i), real(degtorad), cLine, fLine, nLines, &
               "Position angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)

    enddo
 end if
       
  end subroutine readImageParameters

  subroutine readSpectrumParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok


    call getBigInteger("nphotons", nPhotons, cLine, fLine, nLines, &
         "Number of photons in SED: ", "(a,i15,1x,a)", 100000, ok, .false.)

    call getInteger("ninc", nInclination, cLine, fLine, nLines, &
         "Number of inclination angles: ", "(a,i3,1x,a)", 1, ok, .false.)

    allocate(inclinations(nInclination))
    call findRealArray("inclinations", inclinations, cLine, fLine, nLines, ok)
    if (ok) then
       call getRealArray("inclinations", inclinations, 1.0, cLine, fLine, nLines, &
            "Inclinations (deg): ",90., ok, .false.)
       inclinations(:) = inclinations(:) * degToRad
    else
       deallocate(inclinations)
       call getReal("firstinc", firstInclination, 1.0, cLine, fLine, nLines, &
            "First inclination angle (deg): ","(a,f4.1,1x,a)", 10., ok, .true.)
       firstInclination = firstInclination * degToRad
       if (nInclination > 1) &
            call getReal("lastinc", lastInclination, 1.0, cLine, fLine, nLines, &
            "Last inclination angle (deg): ","(a,f4.1,1x,a)", 80., ok, .true.)
       lastInclination = lastInclination * degToRad
    end if


    call getReal("distance", gridDistance, 1., cLine, fLine, nLines, &
         "Grid distance (pc): ","(a,f4.1,1x,a)", 100., ok, .false.)

    call getInteger("nphase", nPhase, cLine, fLine, nLines, &
         "Number of phases: ", "(a,i3,1x,a)", 1, ok, .false.)

    call getInteger("nstart", nStartPhase, cLine, fLine, nLines, &
         "Start at phase: ", "(a,i3,1x,a)", 1, ok, .false.)

    call getInteger("nend", nEndPhase, cLine, fLine, nLines, &
         "End at phase: ", "(a,i3,1x,a)", nPhase, ok, .false.)

    ! Used for optical depth tests in phaseloop
    call getReal("lambdatau", lambdatau, 1.0, cLine, fLine, nLines, &
         "Lambda for tau test: ","(a,1PE10.3,1x,a)", 5500.0, ok, .false.)

    ! Parameters of output file
    call getString("filename", outFile, cLine, fLine, nLines, &
         "Output spectrum filename: ","(a,a,1x,a)","spectrum.dat", ok, .false.)

    call getLogical("sed", sed, cLine, fLine, nLines, &
         "Write spectrum as lambda vs lambda Flambda: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("sised", sised, cLine, fLine, nLines, &
         "Write spectrum as lambda (microns) vs lambda F_lambda (microns * W/m^2):","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("jansky", jansky, cLine, fLine, nLines, &
         "Write spectrum in janskies: ","(a,1l,1x,a)", .false., ok, .false.)

    call getReal("sedlammin", SEDlamMin, 1.0e4, cLine, fLine, nLines, &
         "Minimum wavelength output to SED (microns)","(a,1PE10.3,1x,a)", 0.1, ok, .false.)

    call getReal("sedlammax", SEDlamMax, 1.0e4, cLine, fLine, nLines, &
         "Maximum wavelength output to SED (microns)","(a,1PE10.3,1x,a)", 2000.0, ok, .false.)

    call getLogical("sedwavlin", SEDwavLin, cLine, fLine, nLines, &
         "Linear wavelength spacing in SED: ","(a,1l,1x,a)", .false., ok, .false.)

    call getInteger("sednumlam", SEDnumLam, cLine, fLine, nLines, &
         "Number of SED points: ", "(a,i3,1x,a)", 200, ok, .false.)

  end subroutine readSpectrumParameters



!-----------------------------------------------------------------------------------------


subroutine findReal(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 real :: value
 character(len=80) :: cLine(*)
 logical :: fLine(:)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  k = index(cline(i)," ")-1
  j = len_trim(name)
  if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
     if (.not.ok) then
        ok = .true.
        read(cLine(i)(j+1:80),*) value
        fLine(i) = .true.
     else
        call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
        stop
     endif
  endif
 end do
 end subroutine findReal

 subroutine findDouble(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 real(double) :: value
 character(len=80) :: cLine(*)
 logical :: fLine(:)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  j = len_trim(name)
  k = index(cline(i)," ")-1
  if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
     if (.not.ok) then
        ok = .true.
        read(cLine(i)(j+1:80),*) value
        fLine(i) = .true.
     else
     call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
     stop
  endif
endif
end do
 end subroutine findDouble

 subroutine findVECTOR(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 type(VECTOR) :: value
 character(len=80) :: cLine(*)
 logical :: fLine(:)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  j = len_trim(name)
  k = index(cline(i)," ")-1
  if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
     if (.not.ok) then
        ok = .true.
        read(cLine(i)(j+1:80),*) value%x, value%y, value%z
        fline(i) = .true.
  else
     call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
     stop
  endif
     endif


 end do
end subroutine findVECTOR

subroutine findInteger(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 integer :: value
 character(len=80) :: cLine(*)
 logical :: fLine(:)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  j = len_trim(name)
  k = index(cline(i)," ")-1
  if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
     if (.not.ok) then
        ok = .true.
        read(cLine(i)(j+1:),*) value
        fLine(i) = .true.
     else
        call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
        stop
     endif
  endif
end do
end subroutine findInteger

 subroutine findBigInteger(name, value, cLine, fLine, nLines, ok)
   implicit none
   character(len=*) :: name
   logical :: fLine(:)
   integer(kind=bigInt) :: value
   character(len=80) :: cLine(*)
   integer :: nLines
   logical :: ok
   integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  j = len_trim(name)
  k = index(cline(i)," ")-1
  if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
     if (.not.ok) then
        ok = .true.
        read(cLine(i)(j+1:),*) value
        fLine(i) = .true.
     else
        call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
        stop
     endif
  endif
   end do
 end subroutine findBigInteger

subroutine findLogical(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 logical :: value
 character(len=80) :: cLine(*)
 logical :: fLine(:)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  k = index(cline(i)," ")-1
  j = len_trim(name)
  if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
     if (.not.ok) then
        ok = .true.
        read(cLine(i)(j+1:),*) value
        fLine(i) = .true.
     else
        call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
        stop
     endif
  endif
end do
 end subroutine findLogical

subroutine findString(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 character(len=*) :: value
 character(len=80) :: cLine(*)
 logical :: fLine(:)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  j = len_trim(name)
  k = index(cline(i)," ")-1
     if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
  if (.not.ok) then
        ok = .true.
        read(cLine(i)(j+1:),*) value
        fLine(i) = .true.
  else
     call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
     stop
  endif
endif
 end do
 value = trim(value)
 end subroutine findString

subroutine findRealArray(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 real :: value(:)
 character(len=80) :: cLine(*)
 logical :: fLine(:)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  j = len_trim(name)
  k = index(cline(i)," ")-1
     if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
  if (.not.ok) then
        ok = .true.
        read(cLine(i)(j+1:),*) value
        fLine(i) = .true.
  else
     call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
     stop
  endif
endif
 end do
 end subroutine findRealArray

 subroutine getInteger(name, ival, cLine, fLine, nLines, message, format, idef, ok, &
                       musthave)
  character(len=*) :: name
  integer :: ival
  character(len=80) :: cLine(*)
  logical :: Fline(:)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default
  logical :: musthave
  integer :: idef
  logical :: ok
  ok = .true.
  default = " "
  call findInteger(name, ival, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput) write(*,'(a,a)') name, " must be defined"
       stop
    endif
    ival = idef
    default = " (default)"
  endif
  if (musthave.or.(ival /= idef)) then
     write(output,format) trim(message),ival,default
     call writeInfo(output, TRIVIAL)
  endif
 end subroutine getInteger

 subroutine getBigInteger(name, ival, cLine, fLine, nLines, message, format, idef, ok, &
                       musthave)
  character(len=*) :: name
  integer(kind=bigInt) :: ival
  character(len=80) :: cLine(*)
  logical :: fLine(:)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default
  logical :: musthave
  integer :: idef
  logical :: ok
  ok = .true.
  default = " "
  call findBigInteger(name, ival, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput) write(*,'(a,a)') name, " must be defined"
       stop
    endif
    ival = idef
    default = " (default)"
  endif
  if (musthave.or.(ival /= idef)) then
     write(output,format) trim(message),ival,default
     call writeInfo(output, TRIVIAL)
  endif
end subroutine getBigInteger

 subroutine getReal(name, rval, unitConversion, cLine, fLine, nLines, message, format, rdef, ok, &
                    musthave)
  character(len=*) :: name
  real :: rval
  real :: unitConversion
  logical :: musthave
  character(len=80) :: cLine(*)
  logical :: fLine(:)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default
  real :: rdef
  logical :: ok
  ok = .true.
  default = " "
  call findReal(name, rval, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput) write(*,'(a,a)') name, " must be defined"
       stop
    endif
    rval = rdef
    default = " (default)"
 endif
 if (musthave) then
   write(output,format) trim(message)//" ",rval,default
    call writeInfo(output, TRIVIAL)
 endif
 rval = rval * unitConversion
 end subroutine getReal

 logical function checkPresent(keyword, cline, nlines)
   character(len=*) :: keyword
   character(len=80) :: cLine(:)
   integer :: nLines
   integer :: i, j 
   checkPresent = .false.

   j = len(trim(keyword))
   do i = 1, nLines
      if (cline(i)(1:j) == trim(keyword)) then
         checkPresent = .true.
      endif
   enddo
 end function checkPresent


 subroutine getDouble(name, dval, unitConversion, cLine, fLine, nLines, message, format, ddef, ok, &
                    musthave)
  character(len=*) :: name
  real(double) :: dval, unitConversion
  logical :: musthave
  character(len=80) :: cLine(*)
  logical :: fLine(:)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default
  real(double) :: ddef
  logical :: ok
  ok = .true.
  default = " "
  call findDouble(name, dval, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput) write(*,'(a,a)') name, " must be defined"
       stop
    endif
    dval = ddef
    default = " (default)"
 endif
 if (musthave) then
    write(output,format) trim(message),dval,default
    call writeInfo(output, TRIVIAL)
 endif
 dval = dval  * unitConversion
 end subroutine getDouble

 subroutine getVector(name, dval, unitConversion, cLine, fLine, nLines, message, format, ddef, ok, &
                    musthave)
  character(len=*) :: name
  real(double) :: unitConversion
  type(VECTOR) :: dval, ddef
  logical :: musthave
  character(len=80) :: cLine(*)
  logical :: fLine(:)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default
  logical :: ok
  ok = .true.
  default = " "
  call findVECTOR(name, dval, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput) write(*,'(a,a)') name, " must be defined"
       stop
    endif
    dval = ddef
    default = " (default)"
 endif
 if (musthave) then
    write(output,format) trim(message),dval,default
    call writeInfo(output, TRIVIAL)
 endif
 dval = dval  * unitConversion
end subroutine getVector


 subroutine getString(name, rval, cLine, fLine, nLines, message, format, rdef, ok, &
                      musthave)
  character(len=*) :: name
  character(len=*) :: rval
  logical :: musthave
  character(len=80) :: cLine(*)
  logical :: Fline(:)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default
  character(len=*) :: rdef
  logical :: ok
  ok = .true.
  default = " "
  call findString(name, rval, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput)  write(*,'(a,a)') name, " must be defined"
       stop
    endif
    rval = rdef
    default = " (default)"
  endif
  write(output,format) trim(message)//" ",trim(rval),default
  call writeInfo(output, TRIVIAL)
 end subroutine getString


 subroutine getLogical(name, rval, cLine, fLine, nLines, message, cformat, rdef, ok, &
                       musthave)
  character(len=*) :: name
  logical :: rval
  logical :: musthave
  character(len=80) :: cLine(*)
  logical :: fLine(:)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, cformat
  character(len=10) :: default
  logical :: rdef
  character(len=6) :: trueOrFalse, cf
  logical :: ok, thisIsDefault

  cf = cformat
  ok = .true.
  default = " "
  thisIsDefault = .false.
  call findLogical(name, rval, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput) write(*,'(a,a)') name, " must be defined"
       stop
    endif
    rval = rdef
    default = " (default)"
    thisIsDefault = .true.
  endif
  if (rVal) then
     trueOrFalse = " True"
     else
     trueOrFalse = " False"
  endif


  if (musthave .or. .not.thisIsDefault) then
     write(output,'(a,a,a)') trim(message),trueOrFalse,default
     call writeInfo(output, TRIVIAL)
  endif
 end subroutine getLogical

 subroutine getRealArray(name, rval, unitConversion, cLine, fLine, nLines, message, rdef, ok, &
                      musthave)
  character(len=*) :: name
  real :: rval(:), unitConversion
  logical :: musthave
  character(len=80) :: cLine(*)
  logical :: fLine(:)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message
  character(len=10) :: default
  real :: rdef
  logical :: ok
  ok = .true.
  default = " "
  call findRealArray(name, rval, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput)  write(*,'(a,a)') name, " must be defined"
       stop
    endif
    rval(:) = rdef
    default = " (default)"
  endif
  write(output,*) trim(message)//" ",rval,default
  call writeInfo(output, TRIVIAL)
  rVal = rVal * unitConversion
 end subroutine getRealArray

         
      

end module inputs_mod

