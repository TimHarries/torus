module inputs_mod

  use vector_mod
  use unix_mod
  use messages_mod
  use kind_mod
  use constants_mod
  use units_mod
  use utils_mod
  use parallel_mod
#ifdef USEZLIB
  use zlib_mod, only : uncompressedDumpFiles, buffer, nbuffer
#endif
  implicit none

  include "input_variables.f90"

! Maximum line length in parameters file
    integer, parameter, private :: lencLine=200
! Corresponding format specifier in read statement
    character(len=*), parameter, private :: clineFormat="(a200)"

contains

  subroutine  inputs()

#ifdef MPI
    use mpi
#endif
    implicit none

    integer :: nLines
#ifdef MPI
    integer :: ierr
#endif
    logical :: ok
#ifdef USEZLIB
    logical :: compresseddumpfiles
#endif
    character(len=lencLine), allocatable :: cLine(:)
    character(len=lencLine) :: thisLine
    character(len=80) :: message
    logical :: done
    logical, allocatable :: fLine(:)

    ! character vars for unix environ

    character(len=80) :: dataDirectory

    character(len=100) :: paramFile, arg1

    integer :: error, i
    logical :: parameterCheckMode
    
    datadirectory = " "
    done = .false.
    ok = .true.
    oneKappa = .false.
    monteCarloRT = .false.
    vtkincludeGhosts = .true.
    biasToLyman = .false.

    nDustType = 1
    nDiscModule = 1
    tMinGlobal = 3.
    restartLucy = .false.
    filter_set_name = "natural"
    noscattering = .false.
    forceFirstScat = .false.
    usebias = .true.
    tthresh = 0.
    intProFilename = "none"
    nBlobs = 0
    nLines = 0
    inputnMonte = 0
    inflowTemp = 10.d0

    call setupUnits()
    
    call unixGetEnv("TORUS_JOB_DIR",absolutePath)
!   call get_environment_variable("TORUS_JOB_DIR",absolutePath)

! Parse command line arguments. Currently these will be either the name of the parameters file
! and/or the 'check' argument to run in parameter checking mode.
    parameterCheckMode = .false.
    if (command_argument_count() == 1) then
       call get_command_argument(1, paramFile)
       if ( paramFile == "check" ) then
          parameterCheckMode = .true.
          paramFile = trim(absolutePath)//"parameters.dat"
       endif
    else if (command_argument_count() == 2) then
       call get_command_argument(1, arg1)
       call get_command_argument(2, paramFile)
       if (trim(arg1) == "check") then
          parameterCheckMode = .true.
       else
          call writeWarning("Unrecognised argument "//trim(arg1))
       endif
    else
       paramFile = trim(absolutePath)//"parameters.dat"
    endif
    if (writeoutput) then
       call writeInfo("Parameters file is: "//trim(paramFile),TRIVIAL)
    endif

    call testFileExists(paramFile)

    nLines = file_line_count(paramfile)
    ! nLines+1 element of cline is used in loop which reads parameter file
    ! Needs to be +2 in case the last line doesn't have a newline at the end
    allocate ( cLine(nLines+2) )
    nLines = 0

    open(unit=32, file=paramfile, status='old', iostat=error)
    if (error /=0) then
       print *, 'Panic: parameter file open error, file: ',trim(paramFile)
       call torus_stop
    end if

    do
       nLines = nLines + 1
       read(32,cLineFormat,end=10) cLine(nLines)
       if (index(cline(nLines), char(9))/=0) then
          call writeFatal("Input file must not contain tab characters")
          call torus_abort
       endif
       thisLine = trim(adjustl(cLine(nLines)))
       if ( thisLine(1:1) == "%" .or. thisLine(1:1) == "!" .or. thisLine(1:1) == "#" ) then
          nLines = nLines - 1   ! % # or ! is a comment
          cycle
       endif
       if (thisLine(1:1) == " ") then
          nLines = nLines - 1
          cycle
       endif
       if (index(thisLine,"%") /= 0) then
          i = index(thisLine,"%")
          thisLine = thisLine(:(i-1))
          cLine(nLines) = thisLine
       endif
       
       if (index(thisLine,"!") /= 0) then
          i = index(thisLine,"!")
          thisLine = thisLine(:(i-1))
          cLine(nLines) = thisLine
       endif
       
       if (index(thisLine,"#") /= 0) then
          i = index(thisLine,"#")
          thisLine = thisLine(:(i-1))
          cLine(nLines) = thisLine
       endif

       
    end do
10  continue
    nLines = nLines - 1
    write(message,'(a,i4,a)') "Read in ", nLines, " parameters"
    call writeInfo(message,TRIVIAL)

    allocate(Fline(1:nLines))
    fLine = .false.

    call getInteger("verbosity", verbosityLevel, cLine, fLine, nLines, &
         "Verbosity level: ", "(a,i8,1x,a)", 3, ok, .false.)

    call getbigInteger("seed", inputSeed, cLine, fLine, nLines, &
         "Random number seed","(a,i12,a)", 0_bigInt, ok, .false.)


    call getLogical("binaryxml", useBinaryXMLVTKfiles, cLine, fLine, nLines, &
         "Use binary XML VTK files: ","(a,1l,1x,a)", .true., ok, .false.)



    call getLogical("novtkgrid", noVtkGrid, cLine, fLine, nLines, &
         "Suppress VTK grid files: ","(a,1l,1x,a)", .false., ok, .false.)

#ifdef USEZLIB
    call getLogical("compressdumps", compressedDumpFiles, cLine, fLine, nLines, &
         "Use compressed dump files: ","(a,1l,1x,a)", .false., ok, .false.)
    uncompressedDumpFiles = .not.compressedDumpFiles
    nbuffer = 0
    buffer = 0
#endif
    call getLogical("parallelvtu", parallelVTUfiles, cLine, fLine, nLines, &
         "Use parallel VTU files: ","(a,1l,1x,a)", .false., ok, .false.)

! the grid setup. Either we read the grid in or set it up from scratch

    call writeBanner("Grid setup parameters","*",TRIVIAL)

    call getLogical("gridshuffle", gridshuffle, cLine, fLine, nLines, &
         "Initial grid takes account of quantities: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("splitovermpi", splitOverMPI, cLine, fLine, nLines, &
         "Grid is domain decomposed over MPI: ","(a,1l,1x,a)", .false., ok, .false.)

    if (splitOverMPI) then
       call getInteger("nhydrothreads", nHydroThreadsInput, cLine, fLine, nLines, &
            "Number of threads for domain decomposition: ","(a,i3,a)", 0, ok, .true.)
       call getLogical("loadbalancing", loadBalancing, cLine, fLine, nLines, &
            "Employ load balancing MPI methods: ","(a,1l,1x,a)", .true., ok, .false.)
       if (loadBalancing) then
          call getString("loadmethod", loadBalancingMethod, cLine, fLine, nLines, &
               "Load-balancing photoionization step by: ","(a,a,1x,a)", "crossings", ok, .false.)
       endif
    endif



    call getLogical("debug", debug, cLine, fLine, nLines, &
         "Output debug information: ","(a,1l,1x,a)", .false., ok, .false.)


    call getLogical("multimodels", multimodels, cLine, fLine, nLines, &
         "Perform calculation on multiple models: ","(a,1l,1x,a)", .false., ok, .false.)

    nModelStart = 1
    nModelEnd = 1
    modelStep = 1
    iModel = 1
    if (multimodels) then
       call getInteger("nstart", nModelStart, cLine, fLine, nLines, &
            "Starting number for model: ","(a,i7,a)", 1, ok, .true.)
       call getInteger("nend", nModelEnd, cLine, fLine, nLines, &
            "Starting number for model: ","(a,i7,a)", 1, ok, .true.)
       call getInteger("modelstep", modelStep, cLine, fLine, nLines, &
            "Step for multimodel loop: ","(a,i7,a)", 1, ok, .false.)
    endif


    call getLogical("readgrid", readGrid, cLine, fLine, nLines, &
         "Read grid file: ","(a,1l,1x,a)", .false., ok, .true.)


!    if (multimodels.and.(.not.readgrid)) then
!       write(message,*) "Multiple models specified by readgrid set to false"
!       call writeFatal(message)
!       stop
!    endif

    call getLogical("rhofromtable", rhoFromTable, cLine, fLine, nLines, &
         "Interpolate the density distribution from a table: ","(a,1l,1x,a)", .false., ok, .false.)

    if (rhoFromTable) call getString("rhofile", rhofile, cLine, fLine, nLines, &
         "Grid input filename: ","(a,a,1x,a)","none", ok, .true.)

    if(rhofromtable)        call getInteger("nrholines", nrhoLines, cLine, fLine, nLines, &
            "Number of lines in rhofile: ","(a,i7,a)", 1, ok, .true.)


    call getLogical("singlemegaphoto", singleMegaPhoto, cLine, fLine, nLines, &
         "Do a long photoionization calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    if(readgrid) call getLogical("justdump", justDump, cLine, fLine, nLines, &
         "Dump a vtk file and abort: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("rhospec", densitySpectrum, cLine, fLine, nLines, &
         "Dump a density spectrum and abort: ","(a,1l,1x,a)", .false., ok, .false.)

    call getDouble("normfac", normFac, 1.d0, cLine, fLine, nLines, &
            "RhoSpec normalization factor: ","(a,e12.3,1x,a)", 1.d0, ok, .false.)

    if(readgrid) call getLogical("dumpBisbas", dumpBisbas, cLine, fLine, nLines, &
         "Dump the grid for use in 3D-PDR: ","(a,1l,1x,a)", .false., ok, .false.)

    if (readgrid) call getString("inputfile", gridInputFilename, cLine, fLine, nLines, &
                  "Grid input filename: ","(a,a,1x,a)","none", ok, .true.)

!    if (multimodels.and.(index(gridInputFilename, "*") == 0)) then
!       write(message,*) "Multiple models but input filename has no asterixes"
!       call writeFatal(message)
!       stop
!    endif

    call readGridInitParameters(cLine, fLine, nLines)

    call getLogical("readsources", readsources, cLine, fLine, nLines, &
         "Read sources from a file: ","(a,1l,1x,a)", .false., ok, .false.)


    if (.not. readsources) then
!       if (checkPresent("nsource", cLine, nLines)) then
          call readSourceParameters(cLine, fLine, nLines)
!       endif
    else
       call getString("sourcefile", sourceFilename, cLine, fLine, nLines, &
                  "Source filename: ","(a,a,1x,a)","none", ok, .true.)
    endif



    call readGeometrySpecificParameters(cLine, fLine, nLines)



! the physical ingredients of the model

    call writeBanner("Physical ingredients of model","*",TRIVIAL)

    call getLogical("biophysics", bioPhysics, cLine, fLine, nLines, &
         "Include biophysics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("chemistry", doChemistry, cLine, fLine, nLines, &
         "Include chemistry in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("dustphysics", dustPhysics, cLine, fLine, nLines, &
         "Include dust physics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("atomicphysics", atomicPhysics, cLine, fLine, nLines, &
         "Include atomic physics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("molecularphysics", molecularPhysics, cLine, fLine, nLines, &
         "Include molecular physics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("photoionphysics", photoionPhysics, cLine, fLine, nLines, &
         "Include photoionization physics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("gasopacityphysics", gasOpacityPhysics, cLine, fLine, nLines, &
         "Include gas opacity (TiO, VO etc) in calculation: ","(a,1l,1x,a)", .false., ok, .false.)
    includeGasOpacity = gasopacityphysics

    call getLogical("nbodyphysics", nBodyPhysics, cLine, fLine, nLines, &
         "Include n-body physics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getDouble("timeunit", timeUnit, 1.d0, cLine, fLine, nLines, &
         "Code unit of time: ","(a,e12.3,1x,a)", 1.d0, ok, .false.)

    call getDouble("massunit", massUnit, 1.d0, cLine, fLine, nLines, &
         "Code unit of mass: ","(a,e12.3,1x,a)", 1.d0, ok, .false.)

    call getDouble("lengthunit", lengthUnit, 1.d0, cLine, fLine, nLines, &
         "Code unit of length: ","(a,e12.3,1x,a)", 1.d0, ok, .false.)

    call getDouble("slabwidth", slabwidth, 1.d0, cLine, fLine, nLines, &
         "Width of slab for bubble calc (pc): ","(a,e12.3,1x,a)", 1.d0, ok, .false.)

    slabwidth = slabwidth*pctocm/1.d10

    call getbigInteger("nmonte", inputnMonte, cLine, fLine, nLines, &
         "Number of photon packets","(a,i12,a)", 0_bigInt, ok, .false.)

    call getInteger("maxiter", maxPhotoionIter, cLine, fLine, nLines, &
         "Maximum number of iterations","(a,i8,a)", 10, ok, .false.)

    call getLogical("starburst", starburst, cLine, fline, nLines, &
         "Generate sources as starburst: ", "(a,1l,1x,a)", .false., ok, .false.)

    if (starburst) then
       call getDouble("feedbackdelay", feedbackDelay, 1.d0, cLine, fLine, nLines, &
            "Delay feedback after starburst (tff): ","(a,f6.1,a)", 0.1d0, ok, .false.)

       call getDouble("mstarburst", mStarburst, 1.d0, cLine, fLine, nLines, &
            "Starburst mass (solar masses): ","(a,f6.1,a)", 1000.d0, ok, .true.)

       call getDouble("burstage", burstAge, 1.d0, cLine, fLine, nLines, &
            "Starburst age (years): ","(a,f10.1,a)", 1000.d0, ok, .true.)

       call getDouble("bursttime", burstTime, yearstoSecs, cLine, fLine, nLines, &
            "Starburst time (years): ","(a,f10.1,a)", 0.d0, ok, .false.)

       call getDouble("clusterradius", clusterRadius, pctocm, cLine, fLine, nLines, &
            "Burst cluster radius (pc): ","(a,f10.2,a)", 0.d0, ok, .false.)

       call getString("bursttype", burstType, cLine, fLine, nLines, &
            "Star burst type: ","(a,a,1x,a)","none", ok, .true.)

       if (burstType == 'singlestartest') then
          call getVector("burstposition", burstPosition, 1.d0, cLine, fLine, nLines, &
               "Star position (10^10 cm): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .false.)
       endif

       call getString("imf", imfType, cLine, fLine, nLines, &
            "Initial mass function: ","(a,a,1x,a)","salpeter", ok, .false.)
       call getDouble("imfmin", imfMin, 1.d0, cLine, fLine, nLines, &
            "IMF minimum mass (msol): ","(a,f6.1,a)", 0.8d0, ok, .false.)
       call getDouble("imfmax", imfMax, 1.d0, cLine, fLine, nLines, &
            "IMF maximum mass (msol): ","(a,f6.1,a)", 120.d0, ok, .false.)

    endif

#ifdef CHEMISTRY
    if (doChemistry) then
       call readChemistryParameters(cline, fline ,nLines)
    endif
#endif

    call getLogical("clustersinks", clusterSinks, cLine, fLine, nLines, &
         "Sinks represent clusters: ", "(a,1l,1x,a)", .false., ok, .false.)
    if (clusterSinks) then
       call getString("population", populationMethod, cLine, fLine, nLines, &
            "Cluster population method: ","(a,a,1x,a)","threshold", ok, .true.)

       call getDouble("criticalmass", criticalMass, mSol, cLine, fLine, nLines, &
            "Critical mass for creating subsources (msol): ","(a,f6.1,a)", 600.d0, ok, .false.)

       call getDouble("sfe", starFormationEfficiency, 1.d0, cLine, fLine, nLines, &
            "Clustersink star formation efficiency: ","(a,f6.1,a)", 1.d0, ok, .true.)

       call getLogical("readimf", readimf, cLine, fLine, nLines, &
            "Read IMF from file: ", "(a,1l,1x,a)", .false., ok, .false.)

       if (readIMF) then
          call getString("imffile", imfFilename, cLine, fLine, nLines, &
               "Filename of IMF to read in: ","(a,a,1x,a)","imf.dat", ok, .true.)
       else
          call getDouble("popmass", populationMass, mSol, cLine, fLine, nLines, &
               "Total mass to generate in IMF list (msol): ","(a,e12.3,a)", 1.d5, ok, .false.)
       endif

       call getString("imf", imfType, cLine, fLine, nLines, &
            "Initial mass function: ","(a,a,1x,a)","chabrier", ok, .false.)

       call getDouble("imfmin", imfMin, 1.d0, cLine, fLine, nLines, &
            "IMF minimum mass (msol): ","(a,f6.1,a)", 0.8d0, ok, .false.)

       call getDouble("imfmax", imfMax, 1.d0, cLine, fLine, nLines, &
            "IMF maximum mass (msol): ","(a,f6.1,a)", 120.d0, ok, .false.)
    endif

    if (nBodyPhysics .or. starburst .or. clusterSinks) then
       call getString("spectrumtype", sourceSpectrumType, cLine, fLine, nLines, &
            "Type of stellar spectrum: ","(a,a,1x,a)","tlusty", ok, .false.)
    endif

    if (nBodyPhysics) then
       call getUnitDouble("tend", tEnd, "time", cLine, fLine, nLines, &
            "End time for calculation: ","(a,e12.3,1x,a)", 1.d10, ok, .false.)

       call getUnitDouble("tdump", tDump, "time", cLine, fLine, nLines, &
            "Time between dump files: ","(a,e12.3,1x,a)", 0.d0, ok, .false.)

       call getDouble("eps", inputEps, 1.d0, cLine, fLine, nLines, &
            "Gravity softening length (cm): ","(a,e12.3,1x,a)", 0.d0, ok, .false.)

       call getDouble("accradius", accretionRadius, 1.d0, cLine, fLine, nLines, &
            "Accretion radius of sinks (smallest cell): ","(a,e12.3,1x,a)", 2.6d0, ok, .false.)

       call getLogical("addsinks", addSinkParticles, cLine, fLine, nLines, &
            "Add sink particles: ", "(a,1l,1x,a)", .true., ok, .false.)

       call getLogical("hosokawa", hosokawaTracks, cLine, fLine, nLines, &
            "Use Hosokawa evolutionary tracks: ", "(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("mergesinks", mergeBoundSinks, cLine, fLine, nLines, &
            "Merge bound sinks: ", "(a,1l,1x,a)", .true., ok, .false.)

       call getLogical("movesinks", moveSources, cLine, fLine, nLines, &
            "Allow sources to move: ", "(a,1l,1x,a)", .true., ok, .false.)
    endif

    call getLogical("evolvesources", evolveSources, cLine, fLine, nLines, &
         "Stars follow evolutionary tracks: ", "(a,1l,1x,a)", .true., ok, .false.)

    call getDouble("rhothresh", rhoThreshold, 1.d0, cLine, fLine, nLines, &
         "Threshold density for sink creation (g/cc): ","(a,e12.3,1x,a)", 1.d-14, ok, .false.)

    call getLogical("blockhandout", blockHandout, cLine, fLine, nLines, &
         "Use blockhandout for parallel computations ", "(a,1l,1x,a)", .false., ok, .false.)


    if(photoionphysics) then
       call getLogical("checkForPhoto", checkforphoto, cLine, fLine, nLines, &
            "Check whether a photoionization loop is necessary:", "(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("dustonly", dustOnly, cLine, fLine, nLines, &
            "Use just dust opacities: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("packetsplitting", usePacketSplitting, cLine, fLine, nLines, &
            "Use photon packet splitting: ", "(a,1l,1x,a)", .true., ok, .false.)

       call getLogical("regularVTU", dumpregularVTUS, cLine, fLine, nLines, &
            "Dump VTU files regularly: ", "(a,1l,1x,a)", .false., ok, .false.)
       call getInteger("nsmallpackets", inputnSmallPackets, cLine, fLine, nLines, &
            "Number of small packets per big packet: ","(a,i8,a)", 100, ok, .false.)
    end if

    if (molecularPhysics.and.atomicPhysics) then
       call writeFatal("Cannot conduct any calculation that simultaneously includes atoms and molecules")
    endif

    if (molecularPhysics.and.photoionPhysics) then
       call writeFatal("Cannot conduct any calculation that simultaneously includes photoionization and molecules")
    endif

!    if (.not.(dustPhysics.or.atomicPhysics.or.molecularPhysics.or.photoionPhysics)) then
!       call writeFatal("Must include one of: dustPhysics, atomicPhysics, molecularPhysics, photoionPhysics")
!    endif


    if (dustPhysics) call readDustPhysicsParameters(cLine, fLine, nLines)
    if (atomicPhysics) call readAtomicPhysicsParameters(cLine, fLine, nLines)
    if (molecularPhysics) call readMolecularPhysicsParameters(cLine, fLine, nLines)
    if (bioPhysics) call readBiophysicsParameters(cLine, fLine, nLines)
    if (bioPhysics) call readDetectorParameters(cLine, fLine, nLines)
#ifdef PHOTOION
    if (photoionPhysics) call readPhotoionPhysicsParameters(cLine, fLine, nLines)
#endif

! the type of calculation

    call writeBanner("Type of model calculation","*",TRIVIAL)

    call getLogical("stateq", statisticalEquilibrium, cLine, fLine, nLines, &
         "Perform a statistical equililbrium calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("radeq", radiativeEquilibrium, cLine, fLine, nLines, &
         "Perform a radiative equilibrium calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("photoioneq", photoionEquilibrium, cLine, fLine, nLines, &
         "Perform a photoionization calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("timedep", timeDependentRT, cLine, fLine, nLines, &
         "Time-dependent RT: ", "(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("xray", xraycalc, cLine, fLine, nLines, &
         "Include x-ray treatment: ","(a,1l,1x,a)", .false., ok, .false.)

!    if(xraycalc) then
    call getLogical("useionparam", useionparam, cLine, fLine, nLines, &
         "Use ionization parameter in x-ray treatment: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("xheat", xheat, cLine, fLine, nLines, &
         "do xray heating?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("xrayonly", xrayonly, cLine, fLine, nLines, &
         "Only use xrays: ","(a,1l,1x,a)", .false., ok, .false.)
!    end if


    call getLogical("pdr", pdrcalc, cLine, fLine, nLines, &
         "Include pdr treatment with healpix: ","(a,1l,1x,a)", .false., ok, .false.)

    if(pdrcalc) then
    call getLogical("uvfromphoto", uvfromphoto, cLine, fLine, nLines, &
         "Use UV field from photoionization calculation: ","(a,1l,1x,a)", .false., ok, .false.)

       call getInteger("hlevel", hlevel, cLine, fLine, nLines, &
            "Level of healpix refinement : ","(a,i8,a)", 2, ok, .false.)

       call getDouble("v_turb", v_turb, 1.d0, cLine, fLine, nLines, &
            "Turbulent velocity (cm/s): ", "(a,f7.1,1x,a)", 1.d0, ok, .false.)
    end if

    call getLogical("hydrodynamics", hydrodynamics, cLine, fLine, nLines, &
         "Perform a hydrodynamics calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("nbody", donBodyOnly, cLine, fLine, nLines, &
         "Perform an n-body (bigG=1) calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("nbodytest", nBodyTest, cLine, fLine, nLines, &
         "Test the nbody forces: ","(a,1l,1x,a)", .false., ok, .false.)

    if(photoionEquilibrium .and. hydrodynamics) radiationhydrodynamics=.true.

    call getLogical("radiationHydrodynamics", radiationHydrodynamics, cLine, fLine, nLines, &
         "Perform a radiation-hydrodynamics calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("mchistories", MChistories, cLine, fLine, nLines, &
         "Update Monte Carlo estimator histories: ","(a,1l,1x,a)", .false., ok, .false.)

    if (MChistories) then
       call getDouble("radiationtimescale", radiationTimeScale, 1.d0, cLine, fLine, nLines, &
            "ratio of radiation to hydro timescales: ", "(a,f7.1,1x,a)", 1.0d0, ok, .true.)

       call getInteger("shotnoiseweight", ShotNoiseWeight, cLine, fLine, nLines, &
            "#crossings for rad history weighted<current: ","(a,i8,a)", 400, ok, .true.)
    endif

    call getLogical("caklineopacity", CAKlineOpacity, cLine, fLine, nLines, &
         "use Abbot82 temp invarient form of line driving: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("radforcemonte", RadForceMonte, cLine, fLine, nLines, &
         "use a path length based estimation for the rad pressure: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("habingflux", habingFlux, cLine, fLine, nLines, &
         "Calculate flux between 912 and 2400 A for sources: ","(a,1l,1x,a)", .false., ok, .false.)

    call getDouble("radforcethresh", RadForceThresh, 1.d0, cLine, fLine, nLines, &
         "use a path length based estimation for the rad pressure: ","(a,f7.1,1x,a)", 1.0d-15, ok, .false.)

    call getLogical("accretionfeedback", AccretionFeedback, cLine, fLine, nLines, &
         "re-inject some of the accreted material into the domain: ","(a,1l,1x,a)", .false., ok, .false.)

    call getInteger("accfeedbackcells", AccFeedbackCells, cLine, fLine, nLines, &
         "number of cells to inject into : ","(a,i8,a)", 4, ok, .false.)

    call getDouble("feedbackopenangle", FeedbackTheta0, 1.d0, cLine, fLine, nLines, &
         "Theta_0 in feedback distribution: ", "(a,f7.1,1x,a)", 1.d-2, ok, .false.)

    call getDouble("feedbackfw", Feedbackfw, 1.d0, cLine, fLine, nLines, &
         "f_w in feedback distribution: ", "(a,f7.1,1x,a)", 0.27d0, ok, .false.)

    call getDouble("feedbackfv", Feedbackfv, 1.d0, cLine, fLine, nLines, &
         "f_v in feedback distribution: ", "(a,f7.1,1x,a)", 0.333d0, ok, .false.)

    call getDouble("feedbackstartmass", FeedbackStartMass, msol, cLine, fLine, nLines, &
         "Mass that stellar feedback starts: ", "(a,f7.1,1x,a)", 1.0d0, ok, .false.)

    call getInteger("nhydroperphoto", nHydroPerPhoto, cLine, fLine, nLines, &
         "Number of hydro steps per photoionisation loop: ","(a,i4,a)", 1, ok, .false.)

    if (clusterSinks) then
       call getInteger("nhydroperspectra", nHydroPerSpectra, cLine, fLine, nLines, &
            "Number of hydro steps per spectrum calculation: ","(a,i4,a)", 1, ok, .false.)
    endif

    call getLogical("doselfgrav", doselfgrav, cLine, fLine, nLines, &
         "Use self gravity: ","(a,1l,1x,a)", .false., ok, .false.)

    call getDouble("gravtol", gravTol, 1.d0, cLine, fLine, nLines, &
            "Tolerance for gravity solver: ", "(a,e12.5,1x,a)", 1.d-6, ok, .false.)

    call getLogical("forcevcycle", forceVcycle, cLine, fLine, nLines, &
         "Always do multigrid V cycle: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("redogravonread", redoGravOnRead, cLine, fLine, nLines, &
         "Re-solve self-gravity on read: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("simplegrav", simplegrav, cLine, fLine, nLines, &
         "Do simple self gravity calculation: ","(a,1l,1x,a)", .false., ok, .false.)

!Otherwise periodic boundaries are automatically used
    call getLogical("dirichlet", dirichlet, cLine, fLine, nLines, &
         "Use dirichlet boundary conditions: ","(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("severedamping", severedamping, cLine, fLine, nLines, &
         "Severe damping: ","(a,1l,1x,a)", .false., ok, .false.)

    call getDouble("hydrospeedlimit", hydroSpeedLimit, 1.d5, cLine, fLine, nLines, &
            "Speed limit for hydrodynamic advection (km/s): ", "(a,f7.1,1x,a)", 0.d0, ok, .false.)

    call getLogical("dumpradial", dumpRadial, cLine, fLine, nLines, &
         "Dump a radial slice: ","(a,1l,1x,a)", .false., ok, .false.)

    call getDouble("zetacutoff", zetacutoff, 1.d0, cLine, fLine, nLines, &
            "Dimensionless cutoff radius for BES: ", "(a,es9.3,1x,a)", 3.0d0, ok, .false.)


    if (.not. readsources) then
!       if (checkPresent("nsource", cLine, nLines)) then
          call readSourceParameters(cLine, fLine, nLines)
!       endif
    else
       call getString("sourcefile", sourceFilename, cLine, fLine, nLines, &
                  "Source filename: ","(a,a,1x,a)","none", ok, .true.)
    endif

    call getReal("metallicity", stellarMetallicity, 1., cLine, fLine, nLines, &
         "Metallicity of sources relative to Zsolar: ","(a,f6.1,1x,a)", 1., ok, .false.)

    if (statisticalEquilibrium.and.(.not.(molecularPhysics.or.atomicPhysics))) then
       call writeFatal("Must include either molecularPhysics or atomicPhysics for statistical equilibrium calculation")
    endif

    if (radiativeEquilibrium.and.(.not.dustPhysics)) then
       call writeFatal("Can only perform radiative equilibrium using dust physics")
    endif

    if (photoionEquilibrium.and.(.not.photoionPhysics)) then
       call writeFatal("Can only perform a photoionization calculation using photoionization physics")
    endif

    call getLogical("dosmoothgrid", doSmoothGrid, cLine, fLine, nLines, &
         "Smooth AMR grid: ","(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("logspacegrid", logSpaceGrid, cLine, fLine, nLines, &
         "Try and sample a log spaced  AMR grid (x-dir only): ","(a,1l,1x,a)", .false., ok, .false.)

    if(logspacegrid) then
       call getInteger("npoints", npoints, cLine, fLine, nLines, &
            "Number of points on grid : ","(a,i8,a)", 100, ok, .false.)

       call getInteger("nmag", nMag, cLine, fLine, nLines, &
            "number of orders of magnitude spanned : ","(a,i8,a)", 10, ok, .false.)
    endif
    call getReal("smoothfactor", smoothFactor, 1.0, cLine, fLine, nLines, &
         "Inter-cell maximum ratio before smooth: ","(a,f6.1,1x,a)", 3., ok, .false.)

    call getBigInteger("maxmemory", maxMemoryAvailable, cLine, fLine, nLines, &
         "Maximum memory available (Mb): ","(a,i12,a)", 10000_bigInt, ok, .false.)
    maxMemoryAvailable = maxMemoryAvailable * 1000000

    if (statisticalEquilibrium.and.molecularPhysics) call readMolecularLoopParameters(cLine, fLine, nLines)
    if (statisticalEquilibrium.and.atomicPhysics) call readAtomicLoopParameters(cLine, fLine, nLines)
    if (radiativeEquilibrium) call readRadiativeEquilibriumParameters(cLine, fLine, nLines)
    if (photoionEquilibrium) call readPhotoionEquilibriumParameters(cLine, fLine, nLines)
    if (hydrodynamics) call readHydrodynamicsParameters(cLine, fLine, nLines)
    if (timeDependentRT) call readTimeDependentRTParameters(cLine, fLine, nLines)

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

    call getLogical("photometry", calcPhotometry, cLine, fLine, nLines, &
         "Calculate a data cube: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("image", calcImage, cLine, fLine, nLines, &
         "Calculate an image: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("columnimage", calcColumnImage, cLine, fLine, nLines, &
         "Calculate a column density image: ","(a,1l,1x,a)", .false., ok, .false.)


    call getLogical("movie", calcMovie, cLine, fLine, nLines, &
         "Calculate a movie: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("analysis", doAnalysis, cLine, fLine, nLines, &
         "Perform a grid analysis: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("clusteranalysis", doClusterAnalysis, cLine, fLine, nLines, &
         "Perform a grid analysis for cluster models: ","(a,1l,1x,a)", .false., ok, .false.)
    if (doClusterAnalysis) then
       call getLogical("decouplegasdust", decoupleGasDustTemperature, cLine, fLine, nLines, &
            "Decouple gas and dust temperature: ","(a,1l,1x,a)", .false., ok, .false.)

       call getDouble("edgeradius", edgeRadius, pctocm/1.d10, cLine, fLine, nLines, &
            "Cloud edge radius (pc): ","(a,f5.2,a)", 0.d0, ok, .false.)
!
!       call getString("columnfile", columnImageFilename, cLine, fLine, nLines, &
!            "Output column image filename: ","(a,a,1x,a)","none", ok, .false.)
       call getVector("columndir", columnImageDirection, 1.d0, cLine, fLine, nLines, &
            "Direction for column image: ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 1.d0), ok, .false.)
       call getString("columnaxisunits", columnAxisUnits, cLine, fLine, nLines,&
         "Axis units for column image:", "(a,a,1x,a)", "pc", ok, .false.)
       call getString("columndataunits", columnDataUnits, cLine, fLine, nLines,&
            "Data units for column image:", "(a,a,1x,a)", "g/cm2", ok, .false.)

      call readFitsParameters(cLine, fLine, nLines)

      call getLogical("plotavgtemp", plotAvgTemp, cLine, fLine, nLines, &
            "Plot average temperature pixel-by-pixel: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("globalavgtemp", calculateGlobalAvgTemp, cLine, fLine, nLines, &
            "Calculate global average temperature: ","(a,1l,1x,a)", plotAvgTemp, ok, .false.)
      call getLogical("plotavgtdust", plotAvgTdust, cLine, fLine, nLines, &
            "Plot average tdust pixel-by-pixel: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("globalavgtdust", calculateGlobalAvgTdust, cLine, fLine, nLines, &
            "Calculate global average tdust: ","(a,1l,1x,a)", plotAvgTdust, ok, .false.)
      call getLogical("emissionmeasure", calculateEmissionMeasure, cLine, fLine, nLines, &
            "Plot emission measure: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("nlyman", calculateLymanFlux, cLine, fLine, nLines, &
            "Estimate Lyman flux: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("writelums", writeLums, cLine, fLine, nLines, &
            "Write total luminositys to log: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("nundersampled", findNUndersampled, cLine, fLine, nLines, &
            "Find number of undersampled cells: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("calchabing", findHabing, cLine, fLine, nLines, &
            "Write Habing flux to file: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("calchabingcluster", findHabingCluster, cLine, fLine, nLines, &
            "Write Habing flux to file for clustersinks: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("taufuv", calculateTauFUV, cLine, fLine, nLines, &
            "Calculate tau(FUV) between points: ","(a,1l,1x,a)", .false., ok, .false.)
      call getInteger("nionloops", nClusterIonLoops, cLine, fLine, nLines, &
            "Number of photoionization loops: ","(a,i8,a)", 0, ok, .false.)
      call getInteger("primary", primarySource, cLine, fLine, nLines, &
            "Primary source for G0/tau rays: ","(a,i8,a)", 1, ok, .false.)
      call getLogical("globalavgpres", calculateAvgPressure, cLine, fLine, nLines, &
            "Calculate volume average pressure: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("clusterprofile", clusterRadial, cLine, fLine, nLines, &
            "Dump radial profile for cluster: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("radiustotals", radiusTotals, cLine, fLine, nLines, &
            "Calculate total inside radii: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("hiiradius", findHIIRadius, cLine, fLine, nLines, &
            "Calculate total inside radii: ","(a,1l,1x,a)", .false., ok, .false.)
    endif

    call getLogical("spectrum", calcSpectrum, cLine, fLine, nLines, &
         "Calculate a spectrum: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("sourcehistory", sourceHistory, cLine, fLine, nLines, &
         "Write out the source history: ","(a,1l,1x,a)", .false., ok, .false.)
    call getLogical("subsourcehistory", subsourceHistory, cLine, fLine, nLines, &
         "Write out the subsource history: ","(a,1l,1x,a)", .false., ok, .false.)


    call getLogical("writeradialfile", dowriteRadialFile, cLine, fLine, nLines, &
         "Write out a radial dump: ","(a,1l,1x,a)", .false., ok, .false.)

    if (dowriteRadialFile) then
       call getString("radialfilename", radialFilename, cLine, fLine, nLines, &
            "Radial filename: ","(a,a,1x,a)",  " ", ok, .true.)
    endif

    call getLogical("writescatsurface", scatteringSurface, cLine, fLine, nLines, &
         "Write out a scattering surface file: ","(a,1l,1x,a)", .false., ok, .false.)

    if (scatteringSurface) then
       call getReal("lambdatau", lambdatau, 1.0, cLine, fLine, nLines, &
         "Lambda for tau (angstroms): ","(a,1PE10.3,1x,a)", 5500.0, ok, .true.)

       call getString("scatteringfilename", scatteringSurfaceFilename, cLine, fLine, nLines, &
            "Scattering surface filename: ","(a,a,1x,a)",  " ", ok, .true.)
    endif

    if (sourceHistory) then
       call getString("historyfilename", sourceHistoryFilename, cLine, fLine, nLines, &
            "Source history filename: ","(a,a,1x,a)",  " ", ok, .true.)
    endif

    call getLogical("screened", screened, cLine, fLine, nLines, &
         "Screen stellar photons: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("postsublimate", postsublimate, cLine, fLine, nLines, &
         "sublimate dust for sed/image: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("benchmark", calcBenchmark, cLine, fLine, nLines, &
         "Check a benchmark: ","(a,1l,1x,a)", .false., ok, .false.)

#ifdef FITSCUBE
    if (calcDataCube) call readDataCubeParameters(cLine, fLine, nLines)
#endif
    if (calcImage) call readImageParameters(cLine, fLine, nLines)
    if (calcColumnImage) call readColumnImageParameters(cLine, fLine, nLines)
    if (calcMovie) call readMovieParameters(cLine, fLine, nLines)
    if (calcDataCube.or.calcImage.or.calcMovie.or.calcDustCube.or.calcColumnImage) &
         call readFitsParameters(cLine, fLine, nLines)
    if (calcSpectrum) call readSpectrumParameters(cLine, fLine, nLines)
    if (calcPhotometry) call readPhotometryParameters(cLine, fLine, nLines)

    if (writeoutput) write(*,*) " "

    call writeUnusedKeywords(cLine, fLine, nLines)

    ! Close input file
    close(32)

    deallocate (cLine)

    call sanityCheck()

    if ( parameterCheckMode ) then
       call writeBanner("Parameter check mode: inputs read OK","-",IMPORTANT)
#ifdef MPI
       call mpi_abort(MPI_COMM_WORLD, 1, ierr)
#else
       call torus_stop
#endif

    endif

  end subroutine inputs

  subroutine sanityCheck()
    ! put your sanity checks here

    select case(geometry)

       case("ttauri")
          if (TTauriDisc.and.(.not.DustPhysics).and.(alphaDiscTemp==0.)) then
             call writeFatal("ttauridisc specified but dustphysics is false")
             call torus_abort
          endif

       case DEFAULT

     end select

   end subroutine sanityCheck


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
    character(len=lencLine) :: cLine(:)
    character(len=80) :: message
    integer :: nLines
    logical :: fLine(:)
    logical :: ok
    logical :: setSubRadius, oldWindUnits
    character(len=20) :: heightLabel, betaLabel, dustFracLabel
    character(len=20) :: alphaDiscLabel, betaDiscLabel, hDiscLabel, rinDiscLabel, routDiscLabel
    integer :: i, j
    real(double) :: theta, phi
    character(len=20) :: keyword

    select case(geometry)

       case("WB2014")
          call getLogical("WBvel", WBvel, cLine, fLine, nLines, &
               "Just use cell centered velocity ","(a,1l,1x,a)", .false., ok, .false.)
          call getDouble("ringR", ringR, 1.d0, cLine, fLine, nLines, &
               "Radius of ring (WB2014): ","(a,f7.1,1x,a)", 1.d2, ok, .false.)
          call getDouble("ringdR", ringdR, 1.d0, cLine, fLine, nLines, &
               "Width of ring (WB2014): ","(a,f7.1,1x,a)", 0.d0, ok, .false.)
          call getDouble("WB_Sigma0", WB_Sigma0, 1.d0, cLine, fLine, nLines, &
               "Surface density normalization (WB2014): ","(a,f7.1,1x,a)", 1.d2, ok, .false.)
          call getDouble("WB_gamma", WB_gamma, 1.d0, cLine, fLine, nLines, &
               "Surface density power law (WB2014): ","(a,f7.1,1x,a)", 2.d0, ok, .false.)
          call getDouble("WB_Rc", WB_Rc, 1.d0, cLine, fLine, nLines, &
               "Characteristic radius (WB2014): ","(a,f7.1,1x,a)", 100.d0, ok, .false.)
          call getDouble("WB_q", WB_q, 1.d0, cLine, fLine, nLines, &
               "Temperature power law (WB2014): ","(a,f7.1,1x,a)", 0.d5, ok, .false.)
          call getDouble("Tmid1AU", Tmid1AU, 1.d0, cLine, fLine, nLines, &
               "Mid-plane temperature normalization (WB2014): ","(a,f7.1,1x,a)", 100.d0, ok, .false.)
          call getDouble("Tatm1AU", Tatm1AU, 1.d0, cLine, fLine, nLines, &
               "Atmospheric temperature normalization (WB2014): ","(a,f7.1,1x,a)", 300.d0, ok, .false.)

!  real(double) :: WB_Sigma0
!  real(double) :: WB_gamma
!  real(double) :: WB_q
!  real(double) :: Tmid1AU
!  real(double) :: Tatm1AU
       case("protobin")
       call getReal("beta", beta, 1., cLine, fLine, nLines, &
            "Rotation energy to grav enery: ","(a,f7.0,a)", 1., ok, .true.)
       call getDouble("mass", sphereMass, msol, cLine, fLine, nLines, &
            "Sphere mass (in M_sol): ","(a,f7.1,1x,a)", 1.d-30, ok, .true.)
       call getDouble("radius", sphereRadius, pctocm/1.d10, cLine, fLine, nLines, &
            "Sphere radius (in pc): ","(a,e12.3,1x,a)", 1.d-30, ok, .true.)

       case("mgascii")
       call getString("rhofile", rhofile, cLine, fLine, nLines, &
            "Grid input filename: ","(a,a,1x,a)","none", ok, .true.)

       case("lighthouse")
       call getDouble("cavangle", cavAngle, degToRad, cLine, fLine, nLines, &
            "Cavity angle (deg): ","(a,f5.2,a)", 30.d0, ok, .false.)

       case("slab")
       call getdouble("tauslab", tauSlab, 1.d0, cLine, fLine, nLines, &
            "Optical depth of slab at 5500: ","(a,e12.5,a)", 1.d0, ok, .true.)

       case("bondi")
          call getVector("bondicen", bondiCentre, 1.d0, cLine, fLine, nLines, &
               "Centre of bondi accretion (10^10 cm): ", &
               "(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

       case("fractal")
          call getDouble("rho0", rho0, mhydrogen,cLine, fLine, nLines, &
               "Initial number density: ","(a,f6.1,1x,a)", 100.d0, ok, .true.)


       case("etacar")
          oneKappa = .true.
          fastIntegrate = .false.
          lineEmission = .true.
          monteCarloRT = .true.
          call getReal("lamline", lamLine, 1.,cLine, fLine, nLines, &
               "Line emission wavelength: ","(a,f6.1,1x,a)", 850., ok, .true.)

       case("cmfgen")
          oneKappa = .true.
          fastIntegrate = .false.
          lineEmission = .true.
          monteCarloRT = .true.
          call getReal("lamline", lamLine, 1.,cLine, fLine, nLines, &
               "Line emission wavelength: ","(a,f6.1,1x,a)", 850., ok, .true.)

       call getDouble("CMFGEN_Rmin", CMFGEN_Rmin, 1.d0, cLine, fLine, nLines, &
            "radius of central star  [10^10cm] : ", "(a,es9.3,1x,a)", 1.0d0, ok, .true.)
       rcore = real(CMFGEN_Rmin)      ! [10^10cm]


       call getDouble("CMFGEN_Rmax", CMFGEN_Rmax, 1.d0, cLine, fLine, nLines, &
            "max radius of cmfgen data  [rStar] : ", "(a,es9.3,1x,a)", 1.0d0, ok, .true.)
       CMFGEN_rmax = CMFGEN_Rmin * CMFGEN_Rmax


!       call getString("contflux", contFluxFile, cLine, fLine, nLines, &
!            "Continuum flux filename (primary): ","(a,a,1x,a)","none", ok, .true.)


       call getReal("omega", bigOmega, 1., cLine, fLine, nLines, &
            "Ratio of w/w_c: ","(a,f7.2,1x,a)", 0.0, ok, .true.)

       call getLogical("uniformstar", uniformStar, cLine, fLine, nLines, &
            "Assume a uniform star ","(a,1l,1x,a)", .false., ok, .true.)
       if (uniformStar) then
          call getReal("gamma", eddingtonGamma, 1., cLine, fLine, nLines, &
               "Eddington parameter: ","(a,f7.2,1x,a)", 0.0, ok, .true.)

          call getReal("alpha", alphaCAK, 1., cLine, fLine, nLines, &
               "CAK power-law index alpha: ","(a,f7.2,1x,a)", 0.0, ok, .true.)
          if (bigOmega > sqrt(1.d0-eddingtonGamma)) then
             call writeFatal("Omega limit exceeded")
             call torus_stop
          endif
       endif

    case("triangle", "arbitrary")
       call getDouble("ncol", nCol, 1.d0, cLine, fLine, nLines, &
            "Column density of H_2 (cm^-2): ", "(a,es9.3,1x,a)", 1.0d0, ok, .true.)

       call getDouble("tkinetic", tKinetic, 1.d0, cLine, fLine, nLines, &
            "Kinetic temperature: ", "(a,es9.3,1x,a)", 1.0d0, ok, .true.)

       call getDouble("n2max", n2max, 1.d0, cLine, fLine, nLines, &
            "Maximum N_2 density (cm^-3): ", "(a,es9.3,1x,a)", 1.0d0, ok, .true.)

       if(geometry == "triangle") then
          if (n2max < 0.d0) then
             amrGridSize  = real((nCol/abs(n2max))/1.d10)
             amrGridCentreX = amrGridSize/2.0
          else
             amrGridSize = real((2.d0*ncol/n2max)/1.d10)
             amrGridCentreX = amrGridSize/2.0
             amrGridCentreY = 0.d0
             amrGridCentreZ = 0.d0
          endif
       else
             amrGridCentreX = amrGridSize/2.0
             amrGridCentreY = 0.d0
             amrGridCentreZ = 0.d0
       end if





       case("wrshell")
       call getReal("rcore", rCore, real(rSol/1.e10), cLine, fLine, nLines, &
            "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)


       call getReal("rinner", rInner, rcore, cLine, fLine, nLines, &
            "Inner Radius (stellar radii): ","(a,f5.1,a)", 12., ok, .true.)

       call getReal("router", rOuter, rinner, cLine, fLine, nLines, &
            "Outer Radius (inner radius): ","(a,f5.1,a)", 20., ok, .true.)

       call getReal("vterm", vterm, 1.e5, cLine, fLine, nLines, &
            "Terminal velocity (km/s): ","(a,f7.1,a)", 20., ok, .true.)

       call getReal("mdot", mdot, real(msol/(365.25*24*3600.)), cLine, fLine, nLines, &
            "Mass-loss rate (solar masses/year): ","(a,f7.3,a)", 1., ok, .true.)

       call getReal("beta", beta, 1., cLine, fLine, nLines, &
            "Wind beta law index: ","(a,f7.0,a)", 1., ok, .true.)

     case("empty")
       call getDouble("centralmass", centralMass, msol, cLine, fLine, nLines, &
            "Central mass (in M_sol): ","(a,f7.1,1x,a)", 1.d-30, ok, .true.)

     case("uniformden")
       call getDouble("griddensity", gridDensity, 1.d0, cLine, fLine, nLines, &
            "Initial density of grid (in g cm^-3): ","(a,e12.3,1x,a)", 1.d-18, ok, .false.)

     case("arthur06")
       call getDouble("arthurn", arthurN0, 1.d0, cLine, fLine, nLines, &
            "Density at origin (in cm^-3): ","(a,e12.3,1x,a)", 8.d3 , ok, .false.)

       call getDouble("arthurscaleheight", arthurScaleHeight, pctocm/1.d10, cLine, fLine, nLines, &
            "Scale height (in pc): ","(a,e12.3,1x,a)", 0.05d0 , ok, .false.)

     case("shell")
          call getReal("rinner", rInner, real(rSol)/1.e10, cLine, fLine, nLines, &
               "Inner Radius (solar radii): ","(a,f7.3,a)", 12., ok, .true.)
          call getReal("router", rOuter, real(rsol)/1.e10, cLine, fLine, nLines, &
               "Outer Radius (solar radii): ","(a,f5.1,a)", 20., ok, .true.)
          call getDouble("alpha", shellalpha, 1.d0, cLine, fLine, nLines, &
               "Shell density power law index: ","(a,f7.3,a)", 12.d0, ok, .true.)
          call getDouble("shellmass", shellMass, 1.d0, cLine, fLine, nLines, &
               "Shell mass (g): ","(a,f5.1,a)", 20.d0, ok, .true.)

     case("unisphere","gravtest")
       oneKappa = .true.

       call getDouble("mass", sphereMass, msol, cLine, fLine, nLines, &
            "Sphere mass (in M_sol): ","(a,f7.1,1x,a)", 1.d-30, ok, .true.)

       call getDouble("radius", sphereRadius, pctocm/1.d10, cLine, fLine, nLines, &
            "Sphere radius (in pc): ","(a,e12.3,1x,a)", 1.d-30, ok, .true.)

       call getDouble("spheredensityfactor", sphereSurroundingsFac, 1.0d0, cLine, fLine, nLines, &
            "Density ratio at edge of sphere:: ","(a,e12.3,1x,a)", 1.d-2, ok, .false.)

       call getVector("position", spherePosition, 1.d0, cLine, fLine, nLines, &
            "Sphere position (10^10 cm): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

       call getVector("velocity", sphereVelocity, 1.d5/cspeed, cLine, fLine, nLines, &
            "Sphere velocity (km/s): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

     case("3dgaussian")
       call getDouble("mass", sphereMass, msol, cLine, fLine, nLines, &
            "Sphere mass (in M_sol): ","(a,f7.1,1x,a)", 1.d-30, ok, .true.)

       call getDouble("radius", sphereRadius, pctocm/1.d10, cLine, fLine, nLines, &
            "Sphere radius (in pc): ","(a,e12.3,1x,a)", 1.d-30, ok, .true.)

       call getVector("position", spherePosition, 1.d0, cLine, fLine, nLines, &
            "Sphere position (10^10 cm): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

       call getVector("velocity", sphereVelocity, 1.d5/cspeed, cLine, fLine, nLines, &
            "Sphere velocity (km/s): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

     case("krumholz")
       call getDouble("surfacedensity", surfacedensity, 1.d0, cLine, fLine, nLines, &
            "Sphere surface density (in g cm^-2): ","(a,f7.1,1x,a)", 1.d0, ok, .true.)

       call getDouble("mass", sphereMass, msol, cLine, fLine, nLines, &
            "Sphere mass (in M_sol): ","(a,e12.3,1x,a)", 1.d-30, ok, .true.)

       call getDouble("radius", sphereRadius, pctocm/1.d10, cLine, fLine, nLines, &
            "Sphere radius (in pc): ","(a,e12.3,1x,a)", 1.d-30, ok, .true.)

       call getVector("position", spherePosition, 1.d0, cLine, fLine, nLines, &
            "Sphere position (10^10 cm): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

       call getVector("velocity", sphereVelocity, 1.d5/cspeed, cLine, fLine, nLines, &
            "Sphere velocity (km/s): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

       call getReal("beta", beta, 1.0, cLine, fLine, nLines, &
            "Density power law: ","(a,f7.2,a)",-1.5, ok, .true.)

     case("sphere")
       call getDouble("mass", sphereMass, msol, cLine, fLine, nLines, &
            "Sphere mass (in M_sol): ","(a,f7.1,1x,a)", 1.d-30, ok, .true.)

       call getDouble("radius", sphereRadius, pctocm/1.d10, cLine, fLine, nLines, &
            "Sphere radius (in pc): ","(a,e12.3,1x,a)", 1.d-30, ok, .true.)

       call getDouble("spheredensityfactor", sphereSurroundingsFac, 1.0d0, cLine, fLine, nLines, &
            "Density ratio at edge of sphere:: ","(a,e12.3,1x,a)", 1.d-2, ok, .false.)

       call getVector("position", spherePosition, 1.d0, cLine, fLine, nLines, &
            "Sphere position (10^10 cm): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

       call getVector("velocity", sphereVelocity, 1.d5/cspeed, cLine, fLine, nLines, &
            "Sphere velocity (km/s): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

       call getReal("beta", beta, 1.0, cLine, fLine, nLines, &
            "Density power law: ","(a,f7.2,a)",-2.0, ok, .true.)

       call getDouble("omega", omega, 1.d0, cLine, fLine, nLines, &
            "Angular frequency of rotation: ","(a,f7.2,a)",1.d-13, ok, .true.)

       call getReal("isothermtemp", isoThermTemp, 1., cLine, fLine, nLines, &
            "Isothermal temperature (K): ","(a,f7.1,1x,a)", 20.0, ok, .false.)


     case("plumber")
       call getDouble("mass", plumberMass, msol, cLine, fLine, nLines, &
            "Plumber filament mass (in M_sol): ","(a,f7.1,1x,a)", 1.d-30, ok, .true.)
       call getDouble("radius", plumberRadius, pctocm/1.d10, cLine, fLine, nLines, &
            "Plumber filament radius (in pc): ","(a,e12.3,1x,a)", 1.d-30, ok, .true.)
       call getDouble("exponent", plumberExponent, 1.d0, cLine, fLine, nLines, &
            "Plumber density exponent at large radii: ","(a,e12.3,1x,a)", 2.d0, ok, .true.)
       call getDouble("omega", omega, 1.d0, cLine, fLine, nLines, &
            "Angular frequency of rotation: ","(a,f7.2,a)",1.d-13, ok, .true.)


     case("magstream")
       call getString("magstreamfile", magStreamFile, cLine, fLine, nLines, &
            "Magnetic field stream file: ","(a,a,1x,a)","constant", ok, .true.)
       call getReal("ttaurirstar", TTauriRstar, real(rsol), cLine, fLine, nLines, &
            "T Tauri stellar radius (in R_sol): ","(a,f7.1,1x,a)", 2.0, ok, .true.)
       rcore = TTauriRstar/1.0e10       ! [10^10cm]

       ttauriRouter = 20.0*ttaurirStar
       magStreamFileDegrees = .false.
       isothermStream = .true.
       isothermTemp = 10000.

       call getDouble("holeradius", holeRadius, dble(ttaurirstar), cLine, fLine, nLines, &
            "Size of hole in geometrically thin, optically thick disc (in R_star): ","(a,f7.1,1x,a)", &
            1.d10, ok, .false.)


     case("ttauri")
       call getReal("ttaurirstar", TTauriRstar, real(rsol), cLine, fLine, nLines, &
            "T Tauri stellar radius (in R_sol): ","(a,f7.1,1x,a)", 2.0, ok, .true.)
       rcore = TTauriRstar/1.0e10       ! [10^10cm]
       call getReal("ttauridiskheight", TTauriDiskHeight, real(TTauriRstar), cLine, fLine, nLines, &
            "T Tauri disk height (rStar): ","(a,f7.1,1x,a)", 0.05, ok, .false.)
       rcore = TTauriRstar/1.0e10       ! [10^10cm]
       call getReal("ttaurimstar", TTauriMstar, real(msol), cLine, fLine, nLines, &
            "T Tauri stellar mass (in M_sol): ","(a,f7.1,1x,a)", 0.8, ok, .true.)

       call getLogical("ttaurimag", ttauriMagnetosphere, cLine, fLine, nLines, &
            "Include a T Tauri magnetosphere: ","(a,1l,1x,a)", .false., ok, .true.)

        call getReal("ttaurirouter", TTauriRouter, TTaurirStar, cLine, fLine, nLines, &
             "T Tauri outer flow radius (in R_star): ","(a,f7.1,1x,a)", 3.0, ok, .false.)

        call getReal("ttauririnner", TTauriRinner, TTaurirStar, cLine, fLine, nLines, &
          "T Tauri inner flow radius (in R_star): ","(a,f7.1,1x,a)", 2.2, ok, .false.)




       if (ttauriMagnetosphere) then
          call getReal("ttaurirouter", TTauriRouter, TTaurirStar, cLine, fLine, nLines, &
               "T Tauri outer flow radius (in R_star): ","(a,f7.1,1x,a)", 3.0, ok, .true.)
          call getReal("ttauririnner", TTauriRinner, TTaurirStar, cLine, fLine, nLines, &
            "T Tauri inner flow radius (in R_star): ","(a,f7.1,1x,a)", 2.2, ok, .true.)

          call getReal("thotspot", thotspot, 1., cLine, fLine, nLines, &
               "Hot spot temperature (K): ","(a,f8.1,1x,a)", 0., ok, .false.)
          call getReal("dipoleoffset", dipoleOffset, real(degtorad), cLine, fLine, nLines, &
               "Dipole offset (degrees): ","(a,f7.1,1x,a)", 1.e14, ok, .true.)
          call getLogical("usehartmanntemp", useHartmannTemp, cLine, fLine, nLines, &
               "Use temperatures from Hartmann paper:","(a,1l,1x,a)", .false., ok, .false.)
          call getLogical("isotherm", isoTherm, cLine, fLine, nLines, &
               "Use isothermal temperature :","(a,1l,1x,a)", .false., ok, .false.)
          call getReal("isothermtemp", isoThermTemp, 1., cLine, fLine, nLines, &
               "Isothermal temperature (K): ","(a,f7.1,1x,a)", 6500.0, ok, .false.)
          if (.not.(useHartmannTemp .or. isoTherm)) then
             if (writeoutput)  write(*,'(a)') "WARNING: neither useHartmannTemp nor isoTherm specified!"
             call torus_stop
          end if
       endif

       call getDouble("holeradius", holeRadius, dble(ttaurirstar), cLine, fLine, nLines, &
            "Size of hole in geometrically thin, optically thick disc (in R_star): ","(a,f7.1,1x,a)", &
            1.d10, ok, .false.)



       call getDouble("lxoverlbol", lxOverLBol, 1.d0, cLine, fLine, nLines, &
            "X-ray luminosity  (Bolometric luminosities): ","(a,f7.1,1x,a)", 0.d0, ok, .false.)


       call getDouble("maxcellmass", maxCellMass, 1.d0, cLine, fLine, nLines, &
            "Maximum cell mass after splitting (g): ","(a,e12.5,1x,a)", 1.d30, ok, .false.)

       call getLogical("ttauridisc", ttauriDisc, cLine, fLine, nLines, &
            "Dusty disc around magnetosphere: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("flatdisc", flatDisc, cLine, fLine, nLines, &
            "Geometrically flat blackbody disc: ","(a,1l,1x,a)", .false., ok, .false.)

       if (flatDisc) then
          call getReal("rinner", rInner, ttauriRstar/1.e10, cLine, fLine, nLines, &
               "Inner Radius (stellar radii): ","(a,f7.3,a)", 12., ok, .true.)
          call getReal("router", rOuter, real(autocm/1.e10), cLine, fLine, nLines, &
               "Outer Radius (AU): ","(a,f5.1,a)", 20., ok, .true.)
          call getReal("midplanetemp", midplaneDiscTemp, 1., cLine, fLine, nLines, &
               "Blackbody Disc temperature inside sub radius: ","(a,f8.1,a)", 0., ok, .true.)

          call getReal("midplanepower", midplaneDiscPower, 1., cLine, fLine, nLines, &
               "Blackbody disc power law index inside sub radius: ","(a,f8.1,a)", 0., ok, .true.)
       endif
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

          call getReal("epsilon", epsilonDisc, 1., cLine, fLine, nLines, &
               "Isella and Natta epsilon value: ","(a,f6.4,a)", 1., ok, .false.)

          call getReal("alphadisc", alphaDisc, 1., cLine, fLine, nLines, &
               "Disc alpha parameter: ","(a,f5.3,a)", 2.25, ok, .true.)

          call getReal("betadisc", betaDisc, 1., cLine, fLine, nLines, &
               "Disc beta parameter: ","(a,f5.3,a)", 1.25, ok, .true.)


          call getReal("disctemp", alphaDiscTemp, 1., cLine, fLine, nLines, &
               "Disc temperature inside sub radius: ","(a,f8.1,a)", 0., ok, .true.)

          call getReal("discpower", alphaDiscPower, 1., cLine, fLine, nLines, &
               "Disc temperature power law index inside  sub radius: ","(a,f8.1,a)", 0., ok, .true.)



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
       curtainsPhi1s =    curtainsPhi1s * real(pi/180.0)
       curtainsPhi1e =    curtainsPhi1e * real(pi/180.0)
       curtainsPhi2s =    curtainsPhi2s * real(pi/180.0)
       curtainsPhi2e =    curtainsPhi2e * real(pi/180.0)


       ! The following two are used for "constantcurtans" geometry  (RK)
       call getInteger("curtain_number", curtain_number, cLine, fLine, nLines, &
            "Number of curtains : ","(a,i8,a)", 2, ok, .false.)
       call getReal("curtain_width", curtain_width, 1., cLine, fLine, nLines, &
            "Width of each curtain (degree) : ","(a,f7.1,1x,a)", 120.0, ok, .false.)
       ! converting the curtain width from degrees to radians.
       curtain_width =  curtain_width*real(Pi/180.0)


       call getReal("mdotpar1", MdotParameter1, 1., cLine, fLine, nLines, &
            "1st parameter for accretion rate: ", "(a,e9.3,1x,a)", 1.0, ok, .false.)

       if (tTauriMagnetosphere) then
       call getString("mdottype", mDotType, cLine, fLine, nLines, &
            "T Tauri accretion rate model: ","(a,a,1x,a)","constant", ok, .true.)
       call getReal("mdotpar1", MdotParameter1, 1., cLine, fLine, nLines, &
            "1st parameter for accretion rate: ", "(a,e9.3,1x,a)", 1.0, ok, .true.)
       if (maxCellMass > 1.d29) then
          maxCellMass = 1.d17 * (mdotparameter1/1.d-8)
       endif
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
       endif
       call getLogical("ttauriwarp", ttauriwarp, cLine, fLine, nLines, &
            "Include warped disc around magnetosphere: ","(a,1l,1x,a)", .false., ok, .false.)

       call getDouble("hoverr", hoverr, 1.d0, cLine, fLine, nLines, &
            "Warped disc H/R: ","(a,f7.4,1x,a)", 0.3d0, ok, .false.)

       call getLogical("lte", lte, cLine, fLine, nLines, &
            "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("enhance", enhance, cLine, fLine, nLines, &
            "Accretion enhancement: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("ttauriwind", ttauriWind, cLine, fLine, nLines, &
            "T Tauri disc wind present:","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("ttauristellarwind", ttauriStellarWind, cLine, fLine, nLines, &
            "T Tauri stellar wind present:","(a,1l,1x,a)", .false., ok, .false.)

       if (ttauriwind) then

          call getLogical("oldwindunits", oldWindUnits, cLine, fLine, nLines, &
               "Use old ttauri wind units:","(a,1l,1x,a)", .false., ok, .false.)

          if (oldWindUnits) then
             ! --- parameters for ttauri wind
             call getReal("ttaurirouter", TTauriRouter, TTaurirStar, cLine, fLine, nLines, &
                  "T Tauri outer flow radius (in R_star): ","(a,f7.1,1x,a)", 3.0, ok, .true.)

             call getDouble("DW_Rmin", DW_Rmin, ttaurirOuter/1.d10, cLine, fLine, nLines, &
               "Disc wind:: Inner radius of the wind [magnetospheric radii]: ", &
               "(a,1p,e9.3,1x,a)", 70.0d0, ok, .true.)
             call getDouble("DW_Rmax", DW_rMax, ttaurirOuter/1.d10, cLine, fLine, nLines, &
                  "Disc wind:: Outer radius of the disc [magnetospheric radii]: ", &
                  "(a,1p,e9.3,1x,a)", 700.0d0, ok, .true.)
             call getDouble("DW_d", DW_d, DW_rMin, cLine, fLine, nLines, &
                  "Disc wind:: Wind source displacement [inner wind radii]: ", &
                  "(a,1p,e9.3,1x,a)", 70.0d0, ok, .true.)
             call getDouble("DW_Tmax", DW_Tmax,  1.d0, cLine, fLine, nLines, &
                  "Disc wind:: Temperature of disc at inner radius [K]: ", &
                  "(a,1p,e9.3,1x,a)", 2000.0d0, ok, .true.)
             call getDouble("DW_gamma", DW_gamma,  1.d0, cLine, fLine, nLines, &
                  "Disc wind:: Exponent in the disc temperature power law [-]: ", &
                  "(a,f8.3,1x,a)", -0.5d0, ok, .true.)
             call getDouble("DW_Mdot", DW_Mdot,  1.d0, cLine, fLine, nLines, &
                  "Disc wind:: Total mass-loss rate from disc [Msun/yr]: ", &
                  "(a,1p,e9.3,1x,a)", 1.0d-8, ok, .true.)
             call getDouble("DW_alpha", DW_alpha,  1.d0, cLine, fLine, nLines, &
                  "Disc wind:: Exponent in the mass-loss rate per unit area [-]: ", &
                  "(a,1p,e9.3,1x,a)", 0.5d0, ok, .true.)
             call getDouble("DW_beta", DW_beta,  1.d0, cLine, fLine, nLines, &
                  "Disc wind:: Exponent in the modefied beta-velocity law [-]: ", &
                  "(a,1p,e9.3,1x,a)", 0.5d0, ok, .true.)
             call getDouble("DW_Rs", DW_Rs,  DW_rMin, cLine, fLine, nLines, &
                  "Disc wind:: Effective acceleration length [inner wind radii]: ", &
                  "(a,1p,e9.3,1x,a)", 50.0d0, ok, .true.)
             call getDouble("DW_f", DW_f,  1.d0, cLine, fLine, nLines, &
                  "Disc wind:: Scaling on the terminal velocity [-]: ", &
                  "(a,1p,e9.3,1x,a)", 2.0d0, ok, .true.)
             call getDouble("DW_Twind", DW_Twind,  1.d0, cLine, fLine, nLines, &
                  "Disc wind:: Isothermal temperature of disc wind [K]: ", &
                  "(a,1p,e9.3,1x,a)", 5000.0d0, ok, .true.)
          else
             ! --- parameters for ttauri wind
             call getDouble("DW_Rmin", DW_Rmin, ttaurirstar/1.d10, cLine, fLine, nLines, &
               "Disc wind:: Inner radius of the wind [stellar radii]: ", &
               "(a,1p,e9.3,1x,a)", 70.0d0, ok, .true.)
             call getDouble("DW_Rmax", DW_rMax, ttaurirstar/1.d10, cLine, fLine, nLines, &
                  "Disc wind:: Outer radius of the disc [stellar radii]: ", &
                  "(a,1p,e9.3,1x,a)", 700.0d0, ok, .true.)
             call getDouble("DW_d", DW_d, ttaurirstar/1.d10, cLine, fLine, nLines, &
                  "Disc wind:: Wind source displacement [stellar radii]: ", &
                  "(a,1p,e9.3,1x,a)", 70.0d0, ok, .true.)
             call getDouble("DW_Tmax", DW_Tmax,  1.d0, cLine, fLine, nLines, &
                  "Disc wind:: Temperature of disc at inner radius [K]: ", &
                  "(a,1p,e9.3,1x,a)", 2000.0d0, ok, .true.)
             call getDouble("DW_gamma", DW_gamma,  1.d0, cLine, fLine, nLines, &
                  "Disc wind:: Exponent in the disc temperature power law [-]: ", &
                  "(a,f8.3,1x,a)", -0.5d0, ok, .true.)
             call getDouble("DW_Mdot", DW_Mdot,  1.d0, cLine, fLine, nLines, &
                  "Disc wind:: Total mass-loss rate from disc [Msun/yr]: ", &
                  "(a,1p,e9.3,1x,a)", 1.0d-8, ok, .true.)
             call getDouble("DW_alpha", DW_alpha,  1.d0, cLine, fLine, nLines, &
                  "Disc wind:: Exponent in the mass-loss rate per unit area [-]: ", &
                  "(a,1p,e9.3,1x,a)", 0.5d0, ok, .true.)
             call getDouble("DW_beta", DW_beta,  1.d0, cLine, fLine, nLines, &
                  "Disc wind:: Exponent in the modefied beta-velocity law [-]: ", &
                  "(a,1p,e9.3,1x,a)", 0.5d0, ok, .true.)
             call getDouble("DW_Rs", DW_Rs,  DW_rmin, cLine, fLine, nLines, &
                  "Disc wind:: Effective acceleration length [inner wind radii]: ", &
                  "(a,1p,e9.3,1x,a)", 50.0d0, ok, .true.)
             call getDouble("DW_f", DW_f,  1.d0, cLine, fLine, nLines, &
                  "Disc wind:: Scaling on the terminal velocity [-]: ", &
                  "(a,1p,e9.3,1x,a)", 2.0d0, ok, .true.)
             call getDouble("DW_Twind", DW_Twind,  1.d0, cLine, fLine, nLines, &
                  "Disc wind:: Isothermal temperature of disc wind [K]: ", &
                  "(a,1p,e9.3,1x,a)", 5000.0d0, ok, .true.)
          endif

       endif

       if (TTauristellarWind) then
          ! --- parameters for ttauri wind
          call getDouble("SW_Openangle", SW_openAngle, degToRad, cLine, fLine, nLines, &
               "Stellar wind:: Maximum opening angle of steller wind [deg]: ", &
               "(a,1p,e9.3,1x,a)", 30.0d0, ok, .true.)

          call getDouble("SW_eqGap", SW_eqGap, dble(tTauriRstar), cLine, fLine, nLines, &
               "Stellar wind:: Gap between magnetsphere and wind at equator [ttauriRstar]: ", &
               "(a,1p,e9.3,1x,a)", 3.0d0, ok, .false.)

          call getDouble("SW_rAlfvenMult", SW_alfven, 1.d0, cLine, fLine, nLines, &
               "Stellar wind:: multiple to Max R at which momentum is conserved: ", &
               "(a,1p,e9.3,1x,a)", 1.0d0, ok, .false.)

          call getDouble("SW_Rmin", SW_rMin, ttaurirstar/1.d10, cLine, fLine, nLines, &
               "Stellar wind:: Inner radius of the wind [stellar radii]: ", &
               "(a,1p,e9.3,1x,a)", 70.0d0, ok, .true.)

          call getDouble("SW_Rmax", SW_rMax, ttaurirstar/1.d10, cLine, fLine, nLines, &
               "Stellar wind:: Outer radius of the wind [stellar radii]: ", &
               "(a,1p,e9.3,1x,a)", 100.0d0, ok, .true.)

          call getDouble("SW_Vmin", SW_vMin, 1.d5, cLine, fLine, nLines, &
               "Stellar wind:: Wind base velocity [km/s]: ", &
               "(a,1p,e9.3,1x,a)", 10.0d0, ok, .true.)

          call getDouble("SW_Vmax_xEsc", SW_Vmax, 1.d0, cLine, fLine, nLines, &
               "Stellar wind:: Wind max velocity [Vesc]: ", &
               "(a,1p,e9.3,1x,a)", 2.0d0, ok, .true.)

          call getDouble("SW_beta", SW_beta, 1.d0, cLine, fLine, nLines, &
               "Stellar wind:: beta-velocity index []: ", &
               "(a,1p,e9.3,1x,a)", 2.0d0, ok, .true.)

          call getDouble("SW_Mdot", SW_mdot, dble(mSol / yearsToSecs), cLine, fLine, nLines, &
               "Stellar wind:: mass-loss rate [Msol/yr]: ", &
               "(a,1p,e9.3,1x,a)", 1.0d-7, ok, .true.)

          call getDouble("SW_temp", SW_temperature, 1.d0, cLine, fLine, nLines, &
               "Stellar wind:: temperature [K]: ", &
               "(a,1p,e9.3,1x,a)", 10000.0d0, ok, .true.)

          call getDouble("SW_Veq", SW_veq, 1.d5, cLine, fLine, nLines, &
               "Stellar wind:: stellar equatorial rotation velocity [km/s]: ", &
               "(a,1p,e9.3,1x,a)", 0.d0, ok, .false.)


          if (SW_veq == 0.d0) then
             call getDouble("SW_Prot", SW_protation, 86400.d0, cLine, fLine, nLines, &
                  "Stellar wind:: stellar rotation period [d]: ", &
                  "(a,1p,e9.3,1x,a)", 0.d0, ok, .true.)
          endif

          if ((SW_veq /= 0.d0).and.(SW_Protation /= 0.d0)) then
             call writeFatal("Both stellar rotation period and equatorial rotation velocity set")
             stop
          endif

       endif


       if (useHartmannTemp .and. isoTherm) then
          if (writeoutput)  write(*,'(a)') "WARNING: useHartmannTemp and isoTherm both specified!"
          call torus_stop
       end if


       if (useHartmannTemp) &
            call getReal("maxharttemp", maxHartTemp, 1., cLine, fLine, nLines, &
            "Maximum of Hartmann temperature: ","(a,f7.1,1x,a)", 7436., ok, .false.)
       ! sub options for ttauri geometry
       if (ttau_discwind_on) then   ! commnted out here to make ttaur_turn_off_discwind to work
          ! --- parameters for ttauri wind
          call getDouble("DW_d", DW_d, 1.d0, cLine, fLine, nLines, &
               "Disc wind:: Wind soudce displacement [10^10cm]: ", &
               "(a,1p,e9.3,1x,a)", 70.0d0, ok, .true.)
          call getDouble("DW_Rmin", DW_Rmin,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Inner radius of the disc [magnetosphere radii]: ", &
               "(a,1p,e9.3,1x,a)", 70.0d0, ok, .true.)
          call getDouble("DW_Rmax", DW_Rmax,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Outer radius of the disc [magnetosphere radii]: ", &
               "(a,1p,e9.3,1x,a)", 700.0d0, ok, .true.)
          call getDouble("DW_Tmax", DW_Tmax,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Temperature of disc at inner radius [K]: ", &
               "(a,1p,e9.3,1x,a)", 2000.0d0, ok, .true.)
          call getDouble("DW_gamma", DW_gamma,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Exponent in the disc temperature power law [-]: ", &
               "(a,1p,e9.3,1x,a)", -0.5d0, ok, .true.)
          call getDouble("DW_Mdot", DW_Mdot,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Total mass-loss rate from disc [Msun/yr]: ", &
               "(a,1p,e9.3,1x,a)", 1.0d-8, ok, .true.)
          call getDouble("DW_alpha", DW_alpha,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Exponent in the mass-loss rate per unit area [-]: ", &
               "(a,1p,e9.3,1x,a)", 0.5d0, ok, .true.)
          call getDouble("DW_beta", DW_beta,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Exponent in the modefied beta-velocity law [-]: ", &
               "(a,1p,e9.3,1x,a)", 0.5d0, ok, .true.)
          call getDouble("DW_Rs", DW_Rs,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Effective accerelation length [10^10cm]: ", &
               "(a,1p,e9.3,1x,a)", 50.0d0*DW_Rmin, ok, .true.)
          call getDouble("DW_f", DW_f,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Scaling on the terminal velocity [-]: ", &
               "(a,1p,e9.3,1x,a)", 2.0d0, ok, .true.)
          call getDouble("DW_Twind", DW_Twind,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Isotherma temperature of disc wind [K]: ", &
               "(a,1p,e9.3,1x,a)", 5000.0d0, ok, .true.)
       endif
       if (ttau_jet_on) then  ! commented out here to make ttaur_turn_off_jet to work
          ! --- parameters for ttauri wind
          call getDouble("JET_Rmin", JET_Rmin,  1.d0, cLine, fLine, nLines, &
               "Minmium radius of Jet [10^10 cm]: ", &
               "(a,1p,e9.3,1x,a)", TTauriRouter/1.0d10, ok, .false.)
          call getDouble("JET_theta_j", JET_theta_j,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [deg]  jet opening angle: ", &
               "(a,1p,e9.3,1x,a)", 80.0d0, ok, .true.)
          JET_theta_j = JET_theta_j * (Pi/180.0)  ! converting [deg] to [radians]

          call getDouble("JET_Mdot", JET_Mdot,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [Msun/yr] mass loss rate in the jets: ", &
               "(a,1p,e9.3,1x,a)", 1.0d-9, ok, .true.)
          call getDouble("JET_a_param", JET_a_param,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [-] a parameter in density function: ", &
               "(a,1p,e9.3,1x,a)", 0.8d0, ok, .true.)
          call getDouble("JET_b_param", JET_b_param,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [-] b parameter in density function: ", &
               "(a,1p,e9.3,1x,a)", 2.0d0, ok, .true.)
          call getDouble("JET_Vbase", JET_Vbase,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [km/s] Base velocity of jets: ", &
               "(a,1p,e9.3,1x,a)", 20.0d0, ok, .true.)
          call getDouble("JET_Vinf", JET_Vinf,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [km/s] Terminal velocity of jets: ", &
               "(a,1p,e9.3,1x,a)", 200.0d0, ok, .true.)
          call getDouble("JET_beta", JET_beta,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [-] a parameter in velocity function: ", &
               "(a,1p,e9.3,1x,a)", 0.5d0, ok, .true.)
          call getDouble("JET_gamma", JET_gamma,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [-] a parameter in velocity function: ", &
               "(a,1p,e9.3,1x,a)", 0.05d0, ok, .true.)
          call getDouble("JET_T", JET_T,  1.d0, cLine, fLine, nLines, &
               "TTauri jets:: [K]  Isothermal temperature of jets: ", &
               "(a,1p,e9.3,1x,a)", 1.0d4, ok, .true.)
       endif

       case("wind")
       call getReal("rcore", rCore, real(rSol/1.d10), cLine, fLine, nLines, &
            "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)
       call getReal("teff", teff, 1., cLine, fLine, nLines, &
            "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("twind", twind, teff, cLine, fLine, nLines, &
            "Wind temperature (effective temp)): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("v0", v0, 1.e5, cLine, fLine, nLines, &
            "Wind base velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("vterm", vTerm, 1.e5, cLine, fLine, nLines, &
            "Wind terminal velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("beta", beta, 1., cLine, fLine, nLines, &
            "Wind beta law index: ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("mdot", mdot, real(msol)/(365.25*24.*3600.), cLine, fLine, nLines, &
            "mDot (msol/yr): ","(a,e7.1,a)", 1., ok, .true.)
       call getString("contflux1", contFluxFile, cLine, fLine, nLines, &
            "Continuum flux filename (primary): ","(a,a,1x,a)","none", ok, .true.)




       case("fogel")
       call getString("asciiinfile", intextFilename, cLine, fLine, nLines, &
            "Input ascii file for abundance data: ","(a,a,a)","none", ok, .true.)
       call getString("asciioutfile", outtextFilename, cLine, fLine, nLines, &
            "Output ascii file for abundance data: ","(a,a,a)","none", ok, .true.)
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

       case("radpress")
          call getDouble("rcavity", rCavity, 1.d0, cLine, fLine, nLines, &
               "Cavity radius (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .true.)

          call getDouble("ndensity", ndensity, 1.d0, cLine, fLine, nLines, &
               "Number density (per cm^3): ","(a,1pe8.1,1x,a)", 100.d0, ok, .true.)

       case("benchmark", "RHDDisc", "simpledisc")
          call getReal("radius1", rCore, real(rsol/1.e10), cLine, fLine, nLines, &
               "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)

          call getDouble("extmass", extmass, 1.d0, cLine, fLine, nLines, &
               "Ambient mass: ","(a,f5.1,a)", 1.d-20, ok, .false.)

          call getReal("router", rOuter, real(autocm/1.e10), cLine, fLine, nLines, &
               "Outer Radius (AU): ","(a,f8.2,a)", 20., ok, .true.)

          call getReal("height", height, real(autocm/1.e10), cLine, fLine, nLines, &
               "Scale height (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

          call getReal("rho", rho, 1., cLine, fLine, nLines, &
               "Density: ","(a,e12.5,a)", 1., ok, .true.)

          call getReal("rinner", rInner, real(autocm/1.e10), cLine, fLine, nLines, &
               "Inner Radius for dumpresults (AU): ","(a,1pe8.2,a)", 1.e0, ok, .false.)

       case("molebench")
          call getReal("rinner", rInner, 1., cLine, fLine, nLines, &
               "Inner Radius for dumpresults (10^10cm): ","(a,1pe8.2,a)", 1e4, ok, .true.)

          call getReal("router", rOuter, 1., cLine, fLine, nLines, &
               "Outer Radius (10^10cm): ","(a,1pe8.2,a)", 1e6, ok, .true.)

! -------- Geometries which set up a grid from SPH particles -----------------------
       case("sphfile","molcluster", "theGalaxy", "cluster", "wr104")

          call writeBanner("SPH parameters","#",TRIVIAL)

          call getString("sphdatafilename", sphdatafilename, cLine, fLine, nLines, &
               "Input sph data file: ","(a,a,1x,a)","sph.dat.ascii", ok, .true.)

          call getLogical("dragon", dragon, cLine, fLine, nLines, &
               "SPH file is from dragon SPH: ","(a,1l,a)",.false., ok, .false.)

          call getLogical("sphToGridSimple", sphToGridSimple, cLine, fLine, nLines, &
               "Use simple SPH to grid mapping: ","(a,1l,a)",.false., ok, .false.)

          call getReal("hcritPercentile", hcritPercentile, 1., cLine, fLine, nLines, &
               "Percentile for hcrit: ", "(a,f10.4,1x,a)", 0.80, ok, .false.)

          call getReal("hmaxPercentile", hmaxPercentile, 1., cLine, fLine, nLines, &
               "Percentile for hmax: ", "(a,f10.4,1x,a)", 0.99, ok, .false.)

          call getReal("sphNormLimit", sph_norm_limit, 1., cLine, fLine, nLines, &
               "Limit for SPH normalisation: ", "(a,f10.4,a)", 0.5, ok, .false.)

          call getInteger("kerneltype", kerneltype, cLine, fLine, nLines, &
               "Kernel type (0 is exponential/1 is spline): ","(a,i1,a)",0, ok, .false.)

          call getLogical("variableEta", variableEta, cLine, fLine, nLines, &
               "Allow for variable smoothing eta: ","(a,1l,a)",.false., ok, .false.)

          call getLogical("adddisc", adddisc, cLine, fLine, nLines, &
            "Add inner and outer disc to model: ","(a,1l,a)", .false., ok, .false.)

          call getLogical("guessNe", guessNe, cLine, fLine, nLines, &
            "Guess the electron density based on temperature: ","(a,1l,a)", .false., ok, .false.)

          call getLogical("discardsinks", discardSinks, cLine, fLine, nLines, &
            "Discard sink particles: ","(a,1l,a)", .false., ok, .false.)

          call getLogical("sphonepercell", SphOnePerCell, cLine, fLine, nLines, &
            "Split to one particle per cell:", "(a,1l,1x,a)", .false., ok, .false.)

          call getVector("sphveloffset", sphVelOffset, 1.d0, cLine, fLine, nLines, &
               "SPH offset velocity (km/s): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .false.)

          call getVector("sphposoffset", sphPosOffset, 1.d0, cLine, fLine, nLines, &
               "SPH position offset (au): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .false.)


          
          call getString("inputFileFormat", inputFileFormat, cLine, fLine, nLines, &
               "Input file format: ","(a,a,1x,a)","ascii", ok, .false.)

! Conversion from total density to HI density when chemistry is included
! and grid needs to be loaded with HI density
          call getLogical("convertrhotohi", convertRhoToHI, cLine, fLine, nLines, &
               "Convert density to HI:", "(a,1l,1x,a)", .false., ok, .false.)
          call getLogical("sphwithchem", sphwithchem, cLine, fLine, nLines, &
               "SPH has chemistry information:", "(a,1l,1x,a)", .false., ok, .false.)

! If the dump format is ASCII then we can specify which columns to read H2 and molecular abundance from
          if (inputFileFormat=="ascii" .or. inputFileFormat=="ascii-gadget") then
             ! Keep ih2frac for backwards compatibility
             if (checkPresent("sphh2col", cline, nlines)) then
                call getInteger("sphh2col", sphh2col, cLine, fLine, nLines, &
                     "Column containing H2 fraction: ","(a,1x,i2,a)", 11, ok, .false.)
             else
                call getInteger("ih2frac", sphh2col, cLine, fLine, nLines, &
                     "Column containing H2 fraction: ","(a,1x,i2,a)", 11, ok, .false.)
             endif
             ! Keep iCOfrac for backwards compatibility
             if (checkPresent("sphmolcol", cline, nlines)) then
                call getInteger("sphmolcol", sphmolcol, cLine, fLine, nLines, &
                     "Molecular abundances will be read from column: ","(a,1x,i2,a)", 15, ok, .false.)
             else
                call getInteger("iCOfrac", sphmolcol, cLine, fLine, nLines, &
                     "Column containing CO fraction: ","(a,1x,i2,a)", 15, ok, .false.)
             end if
          end if

! These parameters allow SPH particles within a box to be selected. Particles
! outside the box are discarded
    call getLogical("sphboxcut", sphboxcut, cLine, fLine, nLines, &
         "Select SPH particles within a box: ","(a,1l,1x,a)",.false., ok, .false.)
    call getDouble("sphboxxmin", sphboxxmin, 1.0_db, cLine, fLine, nLines, &
         "SPH box min x: ","(a,1pe8.1,1x,a)", -9.99e99_db, ok, .false.)
    call getDouble("sphboxxmax", sphboxxmax, 1.0_db, cLine, fLine, nLines, &
         "SPH box max x: ","(a,1pe8.1,1x,a)",  9.99e99_db, ok, .false.)
    call getDouble("sphboxymin", sphboxymin, 1.0_db, cLine, fLine, nLines, &
         "SPH box min y: ","(a,1pe8.1,1x,a)", -9.99e99_db, ok, .false.)
    call getDouble("sphboxymax", sphboxymax, 1.0_db, cLine, fLine, nLines, &
         "SPH box max y: ","(a,1pe8.1,1x,a)",  9.99e99_db, ok, .false.)
    call getDouble("sphboxzmin", sphboxzmin, 1.0_db, cLine, fLine, nLines, &
         "SPH box min z: ","(a,1pe8.1,1x,a)", -9.99e99_db, ok, .false.)
    call getDouble("sphboxzmax", sphboxzmax, 1.0_db, cLine, fLine, nLines, &
         "SPH box max z: ","(a,1pe8.1,1x,a)",  9.99e99_db, ok, .false.)

! These parameters allow SPH particles within a sphere to be selected. Particles
! outside the sphere are discarded
    call getLogical("sphspherecut", sphspherecut, cLine, fLine, nLines, &
         "Select SPH particles within a sphere: ","(a,1l,1x,a)",.false., ok, .false.)
    call getDouble("sphspherex", sphspherex, 1.0_db, cLine, fLine, nLines, &
         "SPH sphere x: ","(a,1pe8.1,1x,a)", amrgridcentrex, ok, .false.)
    call getDouble("sphspherey", sphspherey, 1.0_db, cLine, fLine, nLines, &
         "SPH sphere y: ","(a,1pe8.1,1x,a)", amrgridcentrey, ok, .false.)
    call getDouble("sphspherez", sphspherez, 1.0_db, cLine, fLine, nLines, &
         "SPH sphere z: ","(a,1pe8.1,1x,a)", amrgridcentrez, ok, .false.)
    call getDouble("sphsphereradius", sphsphereradius, 1.0_db, cLine, fLine, nLines, &
         "SPH sphere radius: ","(a,1pe8.1,1x,a)", real(0.5*amrgridsize,db), ok, .false.)

! Select SPH particles above a minimum density threshold
        call getDouble("sphdensitycut", sphdensitycut, 1.0_db, cLine, fLine, nLines, &
             "SPH density cut: ","(a,1pe8.1,1x,a)", -1.0_db, ok, .false.)

          if ( geometry == "cluster" ) then
             call getReal("removeradius", rCore, 1.0, cLine, fLine, nLines, &
                  "Clearing radius (Torus units): ","(a,f7.3,a)", 2000.0, ok, .false.)
             call getDouble("dphiref", dphiRefine, 1.d0, cLine, fLine, nLines, &
                  "Level of azimuthal refinement (degrees): ","(a,f5.1,a)", 45.d0, ok, .false.)
          end if

          if (geometry == "wr104") then

             call getReal("distance", gridDistance, 1., cLine, fLine, nLines, &
                  "Grid distance (pc): ","(a,f8.1,1x,a)", 100., ok, .true.)
             call getReal("massenv", massEnvelope, real(mMoon), cLine, fLine,  nLines, &
                  "Envelope dust mass (Moon masses): ","(a,f8.2,a)", 10., ok, .true.)
             call getReal("tthresh", tthresh, 1., cLine, fLine,  nLines, &
                  "Dust sublimation temperature (K): ","(a,f8.2,a)", 10., ok, .true.)
          endif
          call getDouble("rcut", rCut, autocm/1.d10, cLine, fLine, nLines, &
                  "Cut off inner radius (au): ","(a,f5.1,a)", -1.d0, ok, .false.)
          call getLogical("discsplit", doDiscSplit, cLine, fLine, nLines, &
               "Split AMR mesh for disc: ","(a,1l,a)", .false., ok, .false.)

    case("spiral")
       call getReal("tthresh", tthresh, 1., cLine, fLine,  nLines, &
            "Dust sublimation temperature (K): ","(a,f8.2,a)", 0., ok, .true.)
       call getReal("vterm", vterm, 1.e5, cLine, fLine, nLines, &
            "Terminal velocity (km/s): ","(a,f7.1,a)", 20., ok, .true.)
       call getReal("period", period, 24.*3600., cLine, fLine, nLines, &
            "Period (days): ","(a,f7.1,a)", 20., ok, .true.)
       call getReal("massenv", massEnvelope, real(mSol), cLine, fLine,  nLines, &
            "Envelope dust mass (solar masses): ","(a,1pe12.3,a)", 10., ok, .true.)


    case("envelope")
       call getReal("rinner", rInner, real(autocm)/1.e10, cLine, fLine, nLines, &
            "Inner radius (AU): ","(a,f7.1,a)", 20., ok, .true.)
       call getReal("router", rOuter, real(autocm)/1.e10, cLine, fLine, nLines, &
            "Outer radius (AU): ","(a,f7.1,a)", 20., ok, .true.)
       call getDouble("n2max", n2max, 1.d0, cLine, fLine, nLines, &
            "Maximum N_2 density (cm^-3): ", "(a,es9.3,1x,a)", 1.0d0, ok, .true.)
       call getReal("beta", beta, 1., cLine, fLine, nLines, &
            "Power law index: ","(a,f7.0,a)", 1., ok, .true.)

!       call getReal("massenv", massEnvelope, real(mSol), cLine, fLine,  nLines, &
!            "Envelope dust mass (solar masses): ","(a,1pe12.3,a)", 10., ok, .true.)


    case("cassandra")

       call getDouble("mstar", mStar, 1.d0, cLine, fLine, nLines, &
            "Mass of central star (solar masses): ","(a,f7.1,a)", 1.d0, ok, .true.)
       mcore = real(mStar)
       call getDouble("metallicity", metallicity, 1.d0, cLine, fLine, nLines, &
            "Mass of central star (solar masses): ","(a,f7.1,a)", 1.d0, ok, .true.)
       call getreal("mdot", mDot, 1., cLine, fLine, nLines, &
            "Mass accretion rate (solar masses per year): ","(a,1p,e8.2,a)", 1., ok, .true.)
       call getReal("router", rOuter, real(autocm/1.d10), cLine, fLine, nLines, &
            "Outer Radius (AU): ","(a,f5.1,a)", 20., ok, .true.)

       call getReal("rinner", rInner, real(autocm)/1.e10, cLine, fLine, nLines, &
            "Inner radius (AU): ","(a,f7.1,a)", 20., ok, .true.)

       call getDouble("minphi", minPhiResolution, degtorad, cLine, fLine, nLines, &
            "Level of azimuthal refinement (degrees): ","(a,f5.1,a)", 30.d0, ok, .false.)
       !Cass change
       call getReal("Q_irr", Q_irr, 1.0, cLine, fLine, nLines, &
            "Irradiation at fixed Q_irr: ","(a,f5.1,a)", 1., ok, .true.)
       call getReal("T_irr", T_irr, 1.0, cLine, fLine, nLines, &
            "Irradiation at fixed Q_irr: ","(a,f5.1,a)", 1., ok, .true.)
       !call getReal("SpiralA", spiralA, 1.0, cLine, fLine, nLines, &
       !     "A in spiral equation: ","(a,f5.1,a)", 1., ok, .true.)
       !call getReal("SpiralB", spiralB, 1.0, cLine, fLine, nLines, &
       !     "B in spiral equation ","(a,f5.1,a)", 1., ok, .true.)
       call getInteger("irrchoice", irrchoice, cLine, fLine, nLines, &
            "irradiation choice","(a,i8,a)", 1, ok, .true.)
       call getInteger("nspiralarms", nspiralarms, cLine, fLine, nLines, &
            "Number of spiral arms:","(a,i8,a)", 1, ok, .true.)
       call getLogical("fixatalphasat", fixatalphasat, cLine, fLine, &
            &nLines,"Fix alpha at saturation for the disc: ",&
            &"(a,1l,1x,a)", .false., ok, .true.)
       call getLogical("fixatQcrit", fixatQcrit, cLine, fLine, &
            &nLines,"Fix Toomre parameter Q at critical value: ",&
            &"(a,1l,1x,a)", .false., ok, .true.)

   case("whitney")

       call getDouble("erinner", erinner, autocm, cLine, fLine, nLines, &
            "Inner radius of envelope (AU): ","(a,f5.1,a)", 180.d0, ok, .true.)

       call getDouble("erouter", erouter, autocm, cLine, fLine, nLines, &
            "Outer radius of envelope (AU): ","(a,f5.1,a)", 180.d0, ok, .true.)

       call getDouble("mdotenv", mdotenv, msol/(365.25d0*24.d0*3600.d0), cLine, fLine, nLines, &
            "Accretion of envelope (solar/year): ","(a,f5.1,a)", 180.d0, ok, .true.)

       call getDouble("drinner", drinner, autocm, cLine, fLine, nLines, &
            "Inner radius of disc (AU): ","(a,f5.1,a)", 180.d0, ok, .true.)

       call getDouble("drouter", drouter, autocm, cLine, fLine, nLines, &
            "Outer radius of disc (AU): ","(a,f5.1,a)", 180.d0, ok, .true.)

       call getReal("mdisc", mdisc, real(msol), cLine, fLine, nLines, &
            "Mass of disc (solar): ","(a,f5.1,a)", 180.0, ok, .true.)

       call getDouble("cavangle", cavAngle, degToRad, cLine, fLine, nLines, &
            "Cavity angle (deg): ","(a,f5.2,a)", 40.d0, ok, .false.)

       call getDouble("cavdens", cavDens, 1.d0, cLine, fLine, nLines, &
            "Cavity density (g/cc): ","(a,e12.2,a)", 1d-30, ok, .false.)

       call getReal("mass1", mcore, real(msol), cLine, fLine, nLines, &
            "Mass of core (solar): ","(a,f5.1,a)", 180.0, ok, .true.)

    call getDouble("rhofloor", rhoFloor, 1.d0, cLine, fLine, nLines, &
         "Minimum density in advection:  ","(a,e12.3,1x,a)", 1.d-30, ok, .false.)

    case("katie")

       call getRealWithUnits("rinner", rInner, "au", "codeunits", cLine, fLine, nLines, &
            "Inner Radius (stellar radii): ", 12., ok, .true.)

       call getRealWithUnits("router", rOuter, "au", "codeunits", cLine, fLine, nLines, &
            "Outer Radius (AU): ", 20., ok, .true.)

       call getRealWithUnits("height", height, "au", "codeunits", cLine, fLine, nLines, &
            "Scale height (AU): ",1.e0,ok,.true.)

       call getReal("alphadisc", alphaDisc, 1., cLine, fLine, nLines, &
            "Disc alpha parameter: ","(a,f5.3,a)", 2.25, ok, .true.)

       call getReal("betadisc", betaDisc, 1., cLine, fLine, nLines, &
            "Disc beta parameter: ","(a,f5.3,a)", 1.25, ok, .true.)

       call getReal("mdisc", mDisc, real(msol), cLine, fLine, nLines, &
            "Disc mass (solar masses): ","(a,f6.4,a)", 1.e-4, ok, .true.)

       call getReal("heightsplitfac", heightSplitFac, 1., cLine, fLine, nLines, &
            "Splitting factor for scale height (local scale heights): ","(a,f5.2,a)", 0.2, ok, .false.)

    case("shakara", "bec")

       oneKappa = .true.
       call getLogical("gasopacity", includeGasOpacity, cLine, fLine, nLines, &
            "Include gas opacity: ","(a,1l,a)", .false., ok, .false.)

       call getReal("mdot", mdot, real(msol * secstoyears), cLine, fLine, nLines, &
            "Mass accretion  rate (msol/yr): ","(a,1p,e12.5,a)", 0.0,  ok, .false.)

       call getLogical("noscat", noScattering, cLine, fLine, nLines, &
            "No scattering opacity in model: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("smoothinneredge", smoothInnerEdge, cLine, fLine, nLines, &
            "Smooth density drop at inner edge: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("curvedinneredge", curvedInnerEdge, cLine, fLine, nLines, &
            "Curved inner edge: ","(a,1l,1x,a)", .false., ok, .false.)

       call getRealWithUnits("radius1", rCore, "rsol", "codeunits", cLine, fLine, nLines, &
            "Core radius: ", 10., ok, .true.)

       call getRealWithUnits("mdisc", mdisc, "msol", "gram", cLine, fLine, nLines, &
            "Disc mass: ", 12., ok, .true.)
       
       call getRealWithUnits("rinner", rInner, "au", "codeunits", cLine, fLine, nLines, &
            "Inner Radius: ", 12., ok, .true.)

       call getRealWithUnits("router", rOuter, "au", "codeunits", cLine, fLine, nLines, &
            "Outer Radius: ", 20., ok, .true.)

       call getRealWithUnits("height", height, "au", "codeunits", cLine, fLine, nLines, &
            "Scale height: ",1.e0,ok,.true.)

       call getReal("alphadisc", alphaDisc, 1., cLine, fLine, nLines, &
            "Disc alpha parameter: ","(a,f5.3,a)", 2.25, ok, .true.)

       call getReal("betadisc", betaDisc, 1., cLine, fLine, nLines, &
            "Disc beta parameter: ","(a,f5.3,a)", 1.25, ok, .true.)


       call getReal("teff1", teff, 1., cLine, fLine, nLines, &
            "Source effective temperature: ","(a,f7.1,a)", 12., ok, .true.)

       call getReal("rsub", rSublimation, rcore, cLine, fLine, nLines, &
               "Sublimation radius (rstar): ","(a,f5.1,a)", 0., ok, .false.)

       call getLogical("setsubradius", setSubRadius, cLine, fLine, nLines, &
            "Set sublimation radius empirically (Whitney et al 2004): ","(a,1l,1x,a)", .false., ok, .false.)

       call getDouble("rhofloor", rhoFloor, 1.d0, cLine, fLine, nLines, &
            "Minimum density:  ","(a,e12.3,1x,a)", 1.d-30, ok, .false.)

       if (.not.setSubRadius) then
          if (rSublimation == 0.) then
             rSublimation = rInner
             write(message, '(a,f7.1,a)') "Dust sublimation radius is ",rSublimation/rcore, " stellar radii"
             call writeInfo(message,TRIVIAL)
          endif
       else
          rsublimation = rCore * (1600./teff)**(-2.1) ! Robitaille 2006 equation 7
          rinner = rsublimation
          write(message, '(a,f7.1,a)') "Dust sublimation radius is ",rSublimation/rcore, " stellar radii"
          call writeInfo(message,TRIVIAL)
       endif



       call getReal("rgapinner", rGapInner, real(autocm/1.d10), cLine, fLine, nLines, &
            "Inner gap  radius (AU): ","(a,f5.1,a)", 1.e30, ok, .false.)

       call getReal("rgapouter", rGapOuter, real(autocm/1.d10), cLine, fLine, nLines, &
            "Outer gap  radius (AU): ","(a,f5.1,a)", 1.e30, ok, .false.)

       call getReal("rhogap", rhoGap, 1., cLine, fLine, nLines, &
            "Density in gap (g/cc): ","(a,f5.1,a)", 1.e-30, ok, .false.)

       call getReal("deltacav", deltaCav, 1., cLine, fLine, nLines, &
            "Scaling factor for inner disc: ","(a,1p,e9.3,a)", 1., ok, .false.)


       call getDouble("phiref", phiRefine, 1.d0, cLine, fLine, nLines, &
            "Range of azimuthal refinement (degrees): ","(a,f5.1,a)", 180.d0, ok, .false.)

       call getDouble("dphiref", dphiRefine, 1.d0, cLine, fLine, nLines, &
            "Level of azimuthal refinement (degrees): ","(a,f5.1,a)", 10.d0, ok, .false.)

       call getDouble("minphi", minPhiResolution, degtorad, cLine, fLine, nLines, &
            "Level of azimuthal refinement (degrees): ","(a,f5.1,a)", 1.d30, ok, .false.)


       call getReal("mass1", mCore, real(msol), cLine, fLine, nLines, &
            "Core mass (solar masses): ","(a,f8.4,a)", 0.5, ok, .true.)

       call getReal("alphadisc", alphaDisc, 1., cLine, fLine, nLines, &
            "Disc alpha parameter: ","(a,f5.3,a)", 2.25, ok, .true.)

       call getReal("betadisc", betaDisc, 1., cLine, fLine, nLines, &
            "Disc beta parameter: ","(a,f5.3,a)", 1.25, ok, .true.)

       call getLogical("hydro", solveVerticalHydro, cLine, fLine, nLines, &
            "Solve vertical hydrostatical equilibrium: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("opaquecore", opaqueCore, cLine, fLine, nLines, &
            "Opaque Core: ","(a,1l,a)", .true., ok, .false.)

       call getLogical("dospiral", dospiral, cLine, fLine, nLines, &
               "Add a spiral density wave: : ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("planetdisc", planetDisc, cLine, fLine, nLines, &
               "Add a disc around the planet: : ","(a,1l,1x,a)", .false., ok, .false.)

       if (solveVerticalHydro) then
          call getInteger("nhydro", nhydro,  cline, fLine, nLines, &
               "Max number of hydro iterations : ","(a,i4,a)", 5, ok, .true.)
       endif


       call getReal("heightsplitfac", heightSplitFac, 1., cLine, fLine, nLines, &
            "Splitting factor for scale height (local scale heights): ","(a,f5.2,a)", 0.2, ok, .false.)


       call getDouble("erinner", erInner, autocm, cLine, fLine, nLines, &
            "Envelope inner radius (AU): ","(a,f10.2,a)", 100.d0, ok, .false.)

       call getDouble("erouter", erOuter, autocm, cLine, fLine, nLines, &
            "Envelope outer radius (AU): ","(a,f10.2,a)", 1.d5, ok, .false.)


       call getDouble("mdotenv", mDotEnv, msol * secstoyears, cLine, fLine, nLines, &
            "Envelope accretion rate (msol/yr): ","(a,f5.2,a)", 1.d-30, ok, .false.)


       call getDouble("cavangle", cavAngle, degToRad, cLine, fLine, nLines, &
            "Cavity angle (deg): ","(a,f5.2,a)", 40.d0, ok, .false.)

       call getDouble("cavdens", cavDens, 1.d0, cLine, fLine, nLines, &
            "Cavity density (g/cc): ","(a,e12.2,a)", 1d-30, ok, .false.)

       call getDouble("ambientdens", ambientDens, 1.d0, cLine, fLine, nLines, &
            "Ambient density (g/cc): ","(a,e12.2,a)", 1d-30, ok, .false.)

       call getDouble("ramb", ramb, autocm, cLine, fLine, nLines, &
            "Ambient density outer radius (AU): ","(a,e12.2,a)", 1d30, ok, .false.)

       call getDouble("rcavity", rCavity, autocm, cLine, fLine, nLines, &
               "Cavity radius (AU): ","(a,1pe8.1,1x,a)", erOuter/autocm, ok, .false.)

       call getLogical("discwind", discWind, cLine, fLine, nLines, &
               "Include disc wind: : ","(a,1l,1x,a)", .false., ok, .false.)



       if (discwind) then
          ! --- parameters for ttauri wind
          call getDouble("DW_Rmin", DW_Rmin, autocm/1.d10, cLine, fLine, nLines, &
               "Disc wind:: Inner radius of the wind [magnetospheric radii]: ", &
               "(a,1p,e9.3,1x,a)", 70.0d0, ok, .true.)
          call getDouble("DW_Rmax", DW_rMax, autocm/1.d10, cLine, fLine, nLines, &
               "Disc wind:: Outer radius of the disc [magnetospheric radii]: ", &
               "(a,1p,e9.3,1x,a)", 700.0d0, ok, .true.)
          call getDouble("DW_d", DW_d, DW_rMin, cLine, fLine, nLines, &
               "Disc wind:: Wind source displacement [inner wind radii]: ", &
               "(a,1p,e9.3,1x,a)", 70.0d0, ok, .true.)
          call getDouble("DW_Tmax", DW_Tmax,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Temperature of disc at inner radius [K]: ", &
               "(a,1p,e9.3,1x,a)", 2000.0d0, ok, .true.)
          call getDouble("DW_gamma", DW_gamma,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Exponent in the disc temperature power law [-]: ", &
               "(a,f8.3,1x,a)", -0.5d0, ok, .true.)
          call getDouble("DW_Mdot", DW_Mdot,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Total mass-loss rate from disc [Msun/yr]: ", &
               "(a,1p,e9.3,1x,a)", 1.0d-8, ok, .true.)
          call getDouble("DW_alpha", DW_alpha,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Exponent in the mass-loss rate per unit area [-]: ", &
               "(a,1p,e9.3,1x,a)", 0.5d0, ok, .true.)
          call getDouble("DW_beta", DW_beta,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Exponent in the modefied beta-velocity law [-]: ", &
               "(a,1p,e9.3,1x,a)", 0.5d0, ok, .true.)
          call getDouble("DW_Rs", DW_Rs,  DW_rMin, cLine, fLine, nLines, &
               "Disc wind:: Effective acceleration length [inner wind radii]: ", &
               "(a,1p,e9.3,1x,a)", 50.0d0, ok, .true.)
          call getDouble("DW_f", DW_f,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Scaling on the terminal velocity [-]: ", &
               "(a,1p,e9.3,1x,a)", 2.0d0, ok, .true.)
          call getDouble("DW_Twind", DW_Twind,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Isothermal temperature of disc wind [K]: ", &
               "(a,1p,e9.3,1x,a)", 5000.0d0, ok, .true.)
       endif



       call getInteger("ndusttype", nDustType, cLine, fLine, nLines,"Number of different dust types: ","(a,i12,a)",1,ok,.false.)

       call getLogical("dustsettling", dustSettling, cLine, fLine, nLines, &
               "Dust settling model: : ","(a,1l,1x,a)", .false., ok, .false.)


       do i = 1, nDustType


          write(heightLabel, '(a,i1.1)') "fracheight",i
          call getDouble(heightLabel, fracdustHeight(i), 1.d0, cLine, fLine, nLines, &
               "Dust scale height as fraction of gas scale height: ","(a,f10.5,1x,a)", 1.d0, ok, .false.)


          write(heightLabel, '(a,i1.1)') "dustheight",i
          call getDouble(heightLabel, dustHeight(i), autocm/1.d10, cLine, fLine, nLines, &
               "Dust scale height at 100 AU (AU): ","(a,f10.5,1x,a)", dble(height)*1.d10/autocm, ok, .false.)

          write(betaLabel, '(a,i1.1)') "dustbeta",i
          call getDouble(betaLabel, dustBeta(i), 1.d0, cLine, fLine, nLines, &
               "Dust beta law (AU): ","(a,f10.5,1x,a)", dble(betaDisc), ok, .false.)

       enddo

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

       rho0  = real(dble(mDisc)*(betaDisc-alphaDisc+2.) / ( twoPi**1.5 * (dble(height)*1.d10)/dble(100.d0*autocm)**betaDisc  &
            * (dble(rInner)*1.d10)**alphaDisc * &
            (((dble(rOuter)*1.d10)**(betaDisc-alphaDisc+2.)-(dble(rInner)*1.d10)**(betaDisc-alphaDisc+2.))) ))
       if (Writeoutput) write(*,*) "rho0: ",rho0


       if (geometry == "bec") then
          call getInteger("nblobs", nBlobs,  cline, fLine, nLines, &
               "Number of blobs: ","(a,i4,a)", 1, ok, .true.)
          do i = 1, nBlobs
             write(keyword, '(a,i1.1)') "blobpos",i
             call getVector(keyword, blobPos(i), autocm/1.d10, cLine, fLine, nLines, &
               "Blob position (au): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)
             write(keyword, '(a,i1.1)') "aradius",i
             call getDouble(keyword, aRadius(i), autocm/1.d10, cLine, fLine, nLines, &
                  "a radius (au):  ","(a,e12.3,1x,a)", 1.d-30, ok, .true.)
             write(keyword, '(a,i1.1)') "cradius",i             
             call getDouble(keyword, cRadius(i), autocm/1.d10, cLine, fLine, nLines, &
                  "c radius (au):  ","(a,e12.3,1x,a)", 1.d-30, ok, .true.)
             write(keyword, '(a,i1.1)') "theta",i             
             call getDouble(keyword, theta, pi/180.d0, cLine, fLine, nLines, &
                  "Blob theta angle (deg):  ","(a,e12.3,1x,a)", 1.d-30, ok, .true.)
             write(keyword, '(a,i1.1)') "phi",i             
             call getDouble(keyword, phi, pi/180.d0, cLine, fLine, nLines, &
                  "Blob phi angle (deg):  ","(a,e12.3,1x,a)", 1.d-30, ok, .true.)
             blobZVec(i) = VECTOR(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta))
          enddo
       endif


    CASE("spiraldisc")

       oneKappa = .true.
       call getLogical("gasopacity", includeGasOpacity, cLine, fLine, nLines, &
            "Include gas opacity: ","(a,1l,a)", .false., ok, .false.)

       call getReal("mdot", mdot,  real(msol * secstoyears), cLine, fLine, nLines, &
            "Mass accretion  rate (msol/yr): ","(a,1p,e12.5,a)", 0.0,  ok, .false.)

       call getLogical("noscat", noScattering, cLine, fLine, nLines, &
            "No scattering opacity in model: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("smoothinneredge", smoothInnerEdge, cLine, fLine, nLines, &
            "Smooth density drop at inner edge: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("curvedinneredge", curvedInnerEdge, cLine, fLine, nLines, &
            "Curved inner edge: ","(a,1l,1x,a)", .false., ok, .false.)

       call getReal("radius1", rCore, real(rsol/1.e10), cLine, fLine, nLines, &
            "Core radius (solar radii): ","(a,f7.3,a)", 10., ok, .true.)

       call getReal("rinner", rInner, real(autocm/1.d10), cLine, fLine, nLines, &
            "Inner Radius (AU): ","(a,f7.3,a)", 12., ok, .true.)

       call getReal("rspiral", rSpiral, real(autocm/1.d10), cLine, fLine, nLines, &
            "Spiral inner radius (AU): ","(a,f7.3,a)", 12., ok, .true.)

       
       call getReal("teff1", teff, 1., cLine, fLine, nLines, &
            "Source effective temperature: ","(a,f7.1,a)", 12., ok, .true.)

          call getReal("rsub", rSublimation, rcore, cLine, fLine, nLines, &
               "Sublimation radius (rstar): ","(a,f5.1,a)", 0., ok, .false.)

       call getLogical("setsubradius", setSubRadius, cLine, fLine, nLines, &
            "Set sublimation radius empirically (Whitney et al 2004): ","(a,1l,1x,a)", .false., ok, .false.)

       call getDouble("rhofloor", rhoFloor, 1.d0, cLine, fLine, nLines, &
            "Minimum density:  ","(a,e12.3,1x,a)", 1.d-30, ok, .false.)

       if (.not.setSubRadius) then
          if (rSublimation == 0.) then
             rSublimation = rInner
             write(message, '(a,f7.1,a)') "Dust sublimation radius is ",rSublimation/rcore, " stellar radii"
             call writeInfo(message,TRIVIAL)
          endif
       else
          rsublimation = rCore * (1600./teff)**(-2.1) ! Robitaille 2006 equation 7
          write(message, '(a,f7.1,a)') "Dust sublimation radius is ",rSublimation/rcore, " stellar radii"
          call writeInfo(message,TRIVIAL)
       endif



       call getReal("router", rOuter, real(autocm/1.d10), cLine, fLine, nLines, &
            "Outer Radius (AU): ","(a,f5.1,a)", 20., ok, .true.)

       call getReal("rgapinner", rGapInner, real(autocm/1.d10), cLine, fLine, nLines, &
            "Inner gap  radius (AU): ","(a,f5.1,a)", 1.e30, ok, .false.)

       call getReal("rgapouter", rGapOuter, real(autocm/1.d10), cLine, fLine, nLines, &
            "Outer gap  radius (AU): ","(a,f5.1,a)", 1.e30, ok, .false.)

       call getReal("rhogap", rhoGap, 1., cLine, fLine, nLines, &
            "Density in gap (g/cc): ","(a,f5.1,a)", 1.e-30, ok, .false.)

       call getReal("deltacav", deltaCav, 1., cLine, fLine, nLines, &
            "Scaling factor for inner disc: ","(a,1p,e9.3,a)", 1., ok, .false.)


       call getDouble("phiref", phiRefine, 1.d0, cLine, fLine, nLines, &
            "Range of azimuthal refinement (degrees): ","(a,f5.1,a)", 180.d0, ok, .false.)

       call getDouble("dphiref", dphiRefine, 1.d0, cLine, fLine, nLines, &
            "Level of azimuthal refinement (degrees): ","(a,f5.1,a)", 10.d0, ok, .false.)

       call getDouble("minphi", minPhiResolution, degtorad, cLine, fLine, nLines, &
            "Level of azimuthal refinement (degrees): ","(a,f5.1,a)", 1.d30, ok, .false.)


       call getReal("height", height, real(autocm/1.d10), cLine, fLine, nLines, &
            "Scale height (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

       call getReal("mass1", mCore, real(msol), cLine, fLine, nLines, &
            "Core mass (solar masses): ","(a,f8.4,a)", 0.5, ok, .true.)

       call getReal("mdisc", mDisc, real(msol), cLine, fLine, nLines, &
            "Disc mass (solar masses): ","(a,f8.4,a)", 1.e-4, ok, .true.)

       call getReal("alphadisc", alphaDisc, 1., cLine, fLine, nLines, &
            "Disc alpha parameter: ","(a,f5.3,a)", 2.25, ok, .true.)

       call getReal("betadisc", betaDisc, 1., cLine, fLine, nLines, &
            "Disc beta parameter: ","(a,f5.3,a)", 1.25, ok, .true.)

       call getLogical("hydro", solveVerticalHydro, cLine, fLine, nLines, &
            "Solve vertical hydrostatical equilibrium: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("opaquecore", opaqueCore, cLine, fLine, nLines, &
            "Opaque Core: ","(a,1l,a)", .true., ok, .false.)

       call getLogical("dospiral", dospiral, cLine, fLine, nLines, &
               "Add a spiral density wave: : ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("planetdisc", planetDisc, cLine, fLine, nLines, &
               "Add a disc around the planet: : ","(a,1l,1x,a)", .false., ok, .false.)

       if (solveVerticalHydro) then
          call getInteger("nhydro", nhydro,  cline, fLine, nLines, &
               "Max number of hydro iterations : ","(a,i4,a)", 5, ok, .true.)
       endif


       call getReal("heightsplitfac", heightSplitFac, 1., cLine, fLine, nLines, &
            "Splitting factor for scale height (local scale heights): ","(a,f5.2,a)", 0.2, ok, .false.)


       call getDouble("erinner", erInner, autocm, cLine, fLine, nLines, &
            "Envelope inner radius (AU): ","(a,f10.2,a)", 100.d0, ok, .false.)

       call getDouble("erouter", erOuter, autocm, cLine, fLine, nLines, &
            "Envelope outer radius (AU): ","(a,f10.2,a)", 1.d5, ok, .false.)


       call getDouble("mdotenv", mDotEnv, msol * secstoyears, cLine, fLine, nLines, &
            "Envelope accretion rate (msol/yr): ","(a,f5.2,a)", 1.d-30, ok, .false.)


       call getDouble("cavangle", cavAngle, degToRad, cLine, fLine, nLines, &
            "Cavity angle (deg): ","(a,f5.2,a)", 40.d0, ok, .false.)

       call getDouble("cavdens", cavDens, 1.d0, cLine, fLine, nLines, &
            "Cavity density (g/cc): ","(a,e12.2,a)", 1d-30, ok, .false.)


       call getLogical("discwind", discWind, cLine, fLine, nLines, &
               "Include disc wind: : ","(a,1l,1x,a)", .false., ok, .false.)




       if (discwind) then
          ! --- parameters for ttauri wind
          call getDouble("DW_Rmin", DW_Rmin, autocm/1.d10, cLine, fLine, nLines, &
               "Disc wind:: Inner radius of the wind [magnetospheric radii]: ", &
               "(a,1p,e9.3,1x,a)", 70.0d0, ok, .true.)
          call getDouble("DW_Rmax", DW_rMax, autocm/1.d10, cLine, fLine, nLines, &
               "Disc wind:: Outer radius of the disc [magnetospheric radii]: ", &
               "(a,1p,e9.3,1x,a)", 700.0d0, ok, .true.)
          call getDouble("DW_d", DW_d, DW_rMin, cLine, fLine, nLines, &
               "Disc wind:: Wind source displacement [inner wind radii]: ", &
               "(a,1p,e9.3,1x,a)", 70.0d0, ok, .true.)
          call getDouble("DW_Tmax", DW_Tmax,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Temperature of disc at inner radius [K]: ", &
               "(a,1p,e9.3,1x,a)", 2000.0d0, ok, .true.)
          call getDouble("DW_gamma", DW_gamma,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Exponent in the disc temperature power law [-]: ", &
               "(a,f8.3,1x,a)", -0.5d0, ok, .true.)
          call getDouble("DW_Mdot", DW_Mdot,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Total mass-loss rate from disc [Msun/yr]: ", &
               "(a,1p,e9.3,1x,a)", 1.0d-8, ok, .true.)
          call getDouble("DW_alpha", DW_alpha,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Exponent in the mass-loss rate per unit area [-]: ", &
               "(a,1p,e9.3,1x,a)", 0.5d0, ok, .true.)
          call getDouble("DW_beta", DW_beta,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Exponent in the modefied beta-velocity law [-]: ", &
               "(a,1p,e9.3,1x,a)", 0.5d0, ok, .true.)
          call getDouble("DW_Rs", DW_Rs,  DW_rMin, cLine, fLine, nLines, &
               "Disc wind:: Effective acceleration length [inner wind radii]: ", &
               "(a,1p,e9.3,1x,a)", 50.0d0, ok, .true.)
          call getDouble("DW_f", DW_f,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Scaling on the terminal velocity [-]: ", &
               "(a,1p,e9.3,1x,a)", 2.0d0, ok, .true.)
          call getDouble("DW_Twind", DW_Twind,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Isothermal temperature of disc wind [K]: ", &
               "(a,1p,e9.3,1x,a)", 5000.0d0, ok, .true.)
       endif



       call getInteger("ndusttype", nDustType, cLine, fLine, nLines,"Number of different dust types: ","(a,i12,a)",1,ok,.false.)

       call getLogical("dustsettling", dustSettling, cLine, fLine, nLines, &
               "Dust settling model: : ","(a,1l,1x,a)", .false., ok, .false.)


       do i = 1, nDustType


          write(heightLabel, '(a,i1.1)') "fracheight",i
          call getDouble(heightLabel, fracdustHeight(i), 1.d0, cLine, fLine, nLines, &
               "Dust scale height as fraction of gas scale height: ","(a,f10.5,1x,a)", 1.d0, ok, .false.)


          write(heightLabel, '(a,i1.1)') "dustheight",i
          call getDouble(heightLabel, dustHeight(i), autocm/1.d10, cLine, fLine, nLines, &
               "Dust scale height at 100 AU (AU): ","(a,f10.5,1x,a)", dble(height)*1.d10/autocm, ok, .false.)

          write(betaLabel, '(a,i1.1)') "dustbeta",i
          call getDouble(betaLabel, dustBeta(i), 1.d0, cLine, fLine, nLines, &
               "Dust beta law (AU): ","(a,f10.5,1x,a)", dble(betaDisc), ok, .false.)

       enddo

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

       rho0  = real(dble(mDisc)*(betaDisc-alphaDisc+2.) / ( twoPi**1.5 * (dble(height)*1.d10)/dble(100.d0*autocm)**betaDisc  &
            * (dble(rInner)*1.d10)**alphaDisc * &
            (((dble(rOuter)*1.d10)**(betaDisc-alphaDisc+2.)-(dble(rInner)*1.d10)**(betaDisc-alphaDisc+2.))) ))
       if (Writeoutput) write(*,*) "rho0: ",rho0

    case("HD169142")

       rhoAmbient = 1.d-30
       oneKappa = .true.
       call getLogical("dustsettling", dustSettling, cLine, fLine, nLines, &
               "Dust settling model: : ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("gasopacity", includeGasOpacity, cLine, fLine, nLines, &
            "Include gas opacity: ","(a,1l,a)", .false., ok, .false.)

       call getReal("mdot", mdot,  real(msol * secstoyears), cLine, fLine, nLines, &
            "Mass accretion  rate (msol/yr): ","(a,1p,e12.5,a)", 0.0,  ok, .false.)

       call getLogical("noscat", noScattering, cLine, fLine, nLines, &
            "No scattering opacity in model: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("smoothinneredge", smoothInnerEdge, cLine, fLine, nLines, &
            "Smooth density drop at inner edge: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("curvedinneredge", curvedInnerEdge, cLine, fLine, nLines, &
            "Curved inner edge: ","(a,1l,1x,a)", .false., ok, .false.)

       call getReal("radius1", rCore, real(rsol/1.e10), cLine, fLine, nLines, &
            "Core radius (solar radii): ","(a,f7.3,a)", 10., ok, .true.)

       call getReal("rinner", rInner, real(autocm/1.d10), cLine, fLine, nLines, &
            "Inner Radius (AU): ","(a,f7.3,a)", 12., ok, .true.)

       call getReal("teff1", teff, 1., cLine, fLine, nLines, &
            "Source effective temperature: ","(a,f7.1,a)", 12., ok, .true.)

          call getReal("rsub", rSublimation, rcore, cLine, fLine, nLines, &
               "Sublimation radius (rstar): ","(a,f5.1,a)", 0., ok, .false.)

       call getLogical("setsubradius", setSubRadius, cLine, fLine, nLines, &
            "Set sublimation radius empirically (Whitney et al 2004): ","(a,1l,1x,a)", .false., ok, .false.)

       if (.not.setSubRadius) then
          if (rSublimation == 0.) then
             rSublimation = rInner
             write(message, '(a,f7.1,a)') "Dust sublimation radius is ",rSublimation/rcore, " stellar radii"
             call writeInfo(message,TRIVIAL)
          endif
       else
          rsublimation = rCore * (1600./teff)**(-2.1) ! Robitaille 2006 equation 7
          write(message, '(a,f7.1,a)') "Dust sublimation radius is ",rSublimation/rcore, " stellar radii"
          call writeInfo(message,TRIVIAL)
       endif

          call getReal("height", height, real(autocm/1e10), cLine, fLine, nLines, &
               "Scale height (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)


       call getReal("router", rOuter, real(autocm/1.d10), cLine, fLine, nLines, &
            "Outer Radius (AU): ","(a,f5.1,a)", 20., ok, .true.)

       call getReal("rgapinner1", rGapInner1, real(autocm/1.d10), cLine, fLine, nLines, &
            "Inner gap inner radius (AU): ","(a,f5.1,a)", 1.e30, ok, .false.)

       call getReal("rgapouter1", rGapOuter1, real(autocm/1.d10), cLine, fLine, nLines, &
            "Inner gap outer radius (AU): ","(a,f5.1,a)", 1.e30, ok, .false.)

       call getReal("rgapinner2", rGapInner2, real(autocm/1.d10), cLine, fLine, nLines, &
            "Outer gap inner radius (AU): ","(a,f5.1,a)", 1.e30, ok, .false.)

       call getReal("rgapouter2", rGapOuter2, real(autocm/1.d10), cLine, fLine, nLines, &
            "Outer gap outer radius (AU): ","(a,f5.1,a)", 1.e30, ok, .false.)

       call getReal("rhogap1", rhoGap1, 1., cLine, fLine, nLines, &
            "Density in gap 1 (g/cc): ","(a,f5.1,a)", 1.e-30, ok, .false.)

       call getReal("rhogap2", rhoGap2, 1., cLine, fLine, nLines, &
            "Density in gap 2 (g/cc): ","(a,f5.1,a)", 1.e-30, ok, .false.)

       call getReal("deltacav", deltaCav, 1., cLine, fLine, nLines, &
            "Scaling factor for inner disc: ","(a,1p,e9.3,a)", 1., ok, .false.)


       call getDouble("phiref", phiRefine, 1.d0, cLine, fLine, nLines, &
            "Range of azimuthal refinement (degrees): ","(a,f5.1,a)", 180.d0, ok, .false.)

       call getDouble("dphiref", dphiRefine, 1.d0, cLine, fLine, nLines, &
            "Level of azimuthal refinement (degrees): ","(a,f5.1,a)", 10.d0, ok, .false.)

       call getDouble("minphi", minPhiResolution, degtorad, cLine, fLine, nLines, &
            "Level of azimuthal refinement (degrees): ","(a,f5.1,a)", 1.d30, ok, .false.)


       call getReal("heightinner", heightinner, real(autocm/1.d10), cLine, fLine, nLines, &
            "Scale height of inner disc (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

       call getReal("heightouter", heightOuter, real(autocm/1.d10), cLine, fLine, nLines, &
            "Scale height of inner disc (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

       call getReal("ringheight", ringHeight, real(autocm/1.d10), cLine, fLine, nLines, &
            "Scale height of inner disc (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)



       call getDouble("envangle", envAngle, degToRad, cLine, fLine, nLines, &
            "Half opening angle of envelope (deg): ","(a,1pe8.2,a)",1.d0,ok,.true.)

       call getDouble("envrho", envRho, 1.d0, cLine, fLine, nLines, &
            "Envelope density (g/cc): ","(a,1pe8.2,a)",1.d0,ok,.true.)

       call getReal("mass1", mCore, real(msol), cLine, fLine, nLines, &
            "Core mass (solar masses): ","(a,f8.4,a)", 0.5, ok, .true.)

       call getReal("mdisc", mDisc, real(msol), cLine, fLine, nLines, &
            "Disc mass (solar masses): ","(a,f8.4,a)", 1.e-4, ok, .true.)


       call getReal("betadisc", betaDisc, 1., cLine, fLine, nLines, &
            "Disc beta parameter: ","(a,f5.3,a)", 1.25, ok, .true.)

       call getReal("alphadisc", alphaDisc, 1., cLine, fLine, nLines, &
            "Disc alpha parameter: ","(a,f5.3,a)", 1.25, ok, .true.)

       call getLogical("hydro", solveVerticalHydro, cLine, fLine, nLines, &
            "Solve vertical hydrostatical equilibrium: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("opaquecore", opaqueCore, cLine, fLine, nLines, &
            "Opaque Core: ","(a,1l,a)", .true., ok, .false.)

       call getLogical("dospiral", dospiral, cLine, fLine, nLines, &
               "Add a spiral density wave: : ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("planetdisc", planetDisc, cLine, fLine, nLines, &
               "Add a disc around the planet: : ","(a,1l,1x,a)", .false., ok, .false.)

       if (solveVerticalHydro) then
          call getInteger("nhydro", nhydro,  cline, fLine, nLines, &
               "Max number of hydro iterations : ","(a,i4,a)", 5, ok, .true.)
       endif


       call getReal("heightsplitfac", heightSplitFac, 1., cLine, fLine, nLines, &
            "Splitting factor for scale height (local scale heights): ","(a,f5.2,a)", 0.2, ok, .false.)


       call getDouble("erinner", erInner, autocm, cLine, fLine, nLines, &
            "Envelope inner radius (AU): ","(a,f10.2,a)", 100.d0, ok, .false.)

       call getDouble("erouter", erOuter, autocm, cLine, fLine, nLines, &
            "Envelope inner radius (AU): ","(a,f10.2,a)", 1.d5, ok, .false.)


       call getDouble("mdotenv", mDotEnv, msol * secstoyears, cLine, fLine, nLines, &
            "Envelope accretion rate (AU): ","(a,f5.2,a)", 1.d-30, ok, .false.)


       call getDouble("cavangle", cavAngle, degToRad, cLine, fLine, nLines, &
            "Cavity angle (deg): ","(a,f5.2,a)", 40.d0, ok, .false.)

       call getDouble("cavdens", cavDens, 1.d0, cLine, fLine, nLines, &
            "Cavity density (g/cc): ","(a,e12.2,a)", 1d-30, ok, .false.)


       call getLogical("discwind", discWind, cLine, fLine, nLines, &
               "Include disc wind: : ","(a,1l,1x,a)", .false., ok, .false.)
       if (discWind) then
          call getDouble("DW_Rmin", DW_Rmin, autocm/1.d10, cLine, fLine, nLines, &
               "Disc wind:: Inner radius of the wind [AU]: ", &
               "(a,1p,e9.3,1x,a)", 70.0d0, ok, .true.)
          call getDouble("DW_Rmax", DW_rMax, autocm/1.d10, cLine, fLine, nLines, &
               "Disc wind:: Outer radius of the disc [AU]: ", &
               "(a,1p,e9.3,1x,a)", 700.0d0, ok, .true.)
          call getDouble("DW_Mdot", DW_Mdot,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Total mass-loss rate from disc [Msun/yr]: ", &
               "(a,1p,e9.3,1x,a)", 1.0d-8, ok, .true.)
       endif


       call getInteger("ndusttype", nDustType, cLine, fLine, nLines,"Number of different dust types: ","(a,i12,a)",1,ok,.false.)

       call getLogical("dustsettling", dustSettling, cLine, fLine, nLines, &
               "Dust settling model: : ","(a,1l,1x,a)", .false., ok, .false.)


       do i = 1, nDustType


          write(heightLabel, '(a,i1.1)') "dustheight",i
          call getDouble(heightLabel, dustHeight(i), autocm/1.d10, cLine, fLine, nLines, &
               "Dust scale height at 100 AU (AU): ","(a,f10.5,1x,a)", dble(height)*1.d10/autocm, ok, .false.)

          write(betaLabel, '(a,i1.1)') "dustbeta",i
          call getDouble(betaLabel, dustBeta(i), 1.d0, cLine, fLine, nLines, &
               "Dust beta law (AU): ","(a,f10.5,1x,a)", dble(betaDisc), ok, .false.)

       enddo

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

       rho0  = 1.d-10


    case("MWC275")

       oneKappa = .true.
       call getLogical("gasopacity", includeGasOpacity, cLine, fLine, nLines, &
            "Include gas opacity: ","(a,1l,a)", .false., ok, .false.)

       call getDouble("hoverr", hoverr, 1.d0, cLine, fLine, nLines, &
            "Inner scale height enhancement  H/R: ","(a,f7.4,1x,a)", 0.3d0, ok, .true.)

       call getReal("mdot", mdot,  real(msol * secstoyears), cLine, fLine, nLines, &
            "Mass accretion  rate (msol/yr): ","(a,1p,e12.5,a)", 0.0,  ok, .false.)

       call getLogical("noscat", noScattering, cLine, fLine, nLines, &
            "No scattering opacity in model: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("smoothinneredge", smoothInnerEdge, cLine, fLine, nLines, &
            "Smooth density drop at inner edge: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("curvedinneredge", curvedInnerEdge, cLine, fLine, nLines, &
            "Curved inner edge: ","(a,1l,1x,a)", .false., ok, .false.)

       call getReal("radius1", rCore, real(rsol/1.e10), cLine, fLine, nLines, &
            "Core radius (solar radii): ","(a,f7.3,a)", 10., ok, .true.)

       call getReal("rinner", rInner, real(AUtocm/1.d10), cLine, fLine, nLines, &
            "Inner Radius (AU): ","(a,f7.3,a)", 12., ok, .true.)

       call getReal("teff1", teff, 1., cLine, fLine, nLines, &
            "Source effective temperature: ","(a,f7.1,a)", 12., ok, .true.)

          call getReal("rsub", rSublimation, real(autocm/1.d10), cLine, fLine, nLines, &
               "Sublimation radius (rstar): ","(a,f5.1,a)", 0., ok, .false.)

       call getLogical("setsubradius", setSubRadius, cLine, fLine, nLines, &
            "Set sublimation radius empirically (Whitney et al 2004): ","(a,1l,1x,a)", .false., ok, .false.)

       if (.not.setSubRadius) then
          if (rSublimation == 0.) then
             rSublimation = rInner
             write(message, '(a,f7.1,a)') "Dust sublimation radius is ",rSublimation/rcore, " stellar radii"
             call writeInfo(message,TRIVIAL)
          endif
       else
          rsublimation = rCore * (1600./teff)**(-2.1) ! Robitaille 2006 equation 7
          write(message, '(a,f7.1,a)') "Dust sublimation radius is ",rSublimation/rcore, " stellar radii"
          call writeInfo(message,TRIVIAL)
       endif



       call getReal("router", rOuter, real(autocm/1.d10), cLine, fLine, nLines, &
            "Outer Radius (AU): ","(a,f5.1,a)", 20., ok, .true.)

       call getReal("rgapinner1", rGapInner1, real(autocm/1.d10), cLine, fLine, nLines, &
            "Inner gap inner radius (AU): ","(a,f5.1,a)", 1.e30, ok, .false.)

       call getReal("rgapouter1", rGapOuter1, real(autocm/1.d10), cLine, fLine, nLines, &
            "Inner gap outer radius (AU): ","(a,f5.1,a)", 1.e30, ok, .false.)

       call getReal("rgapinner2", rGapInner2, real(autocm/1.d10), cLine, fLine, nLines, &
            "Outer gap inner radius (AU): ","(a,f5.1,a)", 1.e30, ok, .false.)

       call getReal("rgapouter2", rGapOuter2, real(autocm/1.d10), cLine, fLine, nLines, &
            "Outer gap outer radius (AU): ","(a,f5.1,a)", 1.e30, ok, .false.)

       call getReal("rhogap", rhoGap, 1., cLine, fLine, nLines, &
            "Density in gap (g/cc): ","(a,f5.1,a)", 1.e-30, ok, .false.)

       call getReal("deltacav", deltaCav, 1., cLine, fLine, nLines, &
            "Scaling factor for inner disc: ","(a,1p,e9.3,a)", 1., ok, .false.)


       call getDouble("phiref", phiRefine, 1.d0, cLine, fLine, nLines, &
            "Range of azimuthal refinement (degrees): ","(a,f5.1,a)", 180.d0, ok, .false.)

       call getDouble("dphiref", dphiRefine, 1.d0, cLine, fLine, nLines, &
            "Level of azimuthal refinement (degrees): ","(a,f5.1,a)", 10.d0, ok, .false.)

       call getDouble("minphi", minPhiResolution, degtorad, cLine, fLine, nLines, &
            "Level of azimuthal refinement (degrees): ","(a,f5.1,a)", 1.d30, ok, .false.)


       call getReal("height", height, real(autocm/1.d10), cLine, fLine, nLines, &
            "Scale height (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

       call getReal("height2", height2, real(autocm/1.d10), cLine, fLine, nLines, &
            "Scale height of outer disc (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)


       call getReal("mass1", mCore, real(msol), cLine, fLine, nLines, &
            "Core mass (solar masses): ","(a,f8.4,a)", 0.5, ok, .true.)

       call getReal("mdisc", mDisc, real(msol), cLine, fLine, nLines, &
            "Disc mass (solar masses): ","(a,f8.4,a)", 1.e-4, ok, .true.)

       call getReal("alphadisc", alphaDisc, 1., cLine, fLine, nLines, &
            "Disc alpha parameter: ","(a,f5.3,a)", 2.25, ok, .true.)

       call getReal("betadisc", betaDisc, 1., cLine, fLine, nLines, &
            "Disc beta parameter: ","(a,f5.3,a)", 1.25, ok, .true.)

       call getLogical("hydro", solveVerticalHydro, cLine, fLine, nLines, &
            "Solve vertical hydrostatical equilibrium: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("opaquecore", opaqueCore, cLine, fLine, nLines, &
            "Opaque Core: ","(a,1l,a)", .true., ok, .false.)

       call getLogical("dospiral", dospiral, cLine, fLine, nLines, &
               "Add a spiral density wave: : ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("planetdisc", planetDisc, cLine, fLine, nLines, &
               "Add a disc around the planet: : ","(a,1l,1x,a)", .false., ok, .false.)

       if (solveVerticalHydro) then
          call getInteger("nhydro", nhydro,  cline, fLine, nLines, &
               "Max number of hydro iterations : ","(a,i4,a)", 5, ok, .true.)
       endif


       call getReal("heightsplitfac", heightSplitFac, 1., cLine, fLine, nLines, &
            "Splitting factor for scale height (local scale heights): ","(a,f5.2,a)", 0.2, ok, .false.)


       call getDouble("erinner", erInner, autocm, cLine, fLine, nLines, &
            "Envelope inner radius (AU): ","(a,f10.2,a)", 100.d0, ok, .false.)

       call getDouble("erouter", erOuter, autocm, cLine, fLine, nLines, &
            "Envelope inner radius (AU): ","(a,f10.2,a)", 1.d5, ok, .false.)


       call getDouble("mdotenv", mDotEnv, msol * secstoyears, cLine, fLine, nLines, &
            "Envelope accretion rate (AU): ","(a,f5.2,a)", 1.d-30, ok, .false.)


       call getDouble("cavangle", cavAngle, degToRad, cLine, fLine, nLines, &
            "Cavity angle (deg): ","(a,f5.2,a)", 40.d0, ok, .false.)

       call getDouble("cavdens", cavDens, 1.d0, cLine, fLine, nLines, &
            "Cavity density (g/cc): ","(a,e12.2,a)", 1d-30, ok, .false.)


       call getLogical("discwind", discWind, cLine, fLine, nLines, &
               "Include disc wind: : ","(a,1l,1x,a)", .false., ok, .false.)
       if (discWind) then
          call getDouble("DW_Rmin", DW_Rmin, autocm/1.d10, cLine, fLine, nLines, &
               "Disc wind:: Inner radius of the wind [AU]: ", &
               "(a,1p,e9.3,1x,a)", 70.0d0, ok, .true.)
          call getDouble("DW_Rmax", DW_rMax, autocm/1.d10, cLine, fLine, nLines, &
               "Disc wind:: Outer radius of the disc [AU]: ", &
               "(a,1p,e9.3,1x,a)", 700.0d0, ok, .true.)
          call getDouble("DW_Mdot", DW_Mdot,  1.d0, cLine, fLine, nLines, &
               "Disc wind:: Total mass-loss rate from disc [Msun/yr]: ", &
               "(a,1p,e9.3,1x,a)", 1.0d-8, ok, .true.)
       endif


       call getInteger("ndusttype", nDustType, cLine, fLine, nLines,"Number of different dust types: ","(a,i12,a)",1,ok,.false.)

       call getLogical("dustsettling", dustSettling, cLine, fLine, nLines, &
               "Dust settling model: : ","(a,1l,1x,a)", .false., ok, .false.)


       do i = 1, nDustType


          write(heightLabel, '(a,i1.1)') "dustheight",i
          call getDouble(heightLabel, dustHeight(i), autocm/1.d10, cLine, fLine, nLines, &
               "Dust scale height at 100 AU (AU): ","(a,f10.5,1x,a)", dble(height)*1.d10/autocm, ok, .false.)

          write(betaLabel, '(a,i1.1)') "dustbeta",i
          call getDouble(betaLabel, dustBeta(i), 1.d0, cLine, fLine, nLines, &
               "Dust beta law (AU): ","(a,f10.5,1x,a)", dble(betaDisc), ok, .false.)

       enddo

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

       rho0  = real(mDisc *(betaDisc-alphaDisc+2.) / ( twoPi**1.5 * (height*1.e10)/real(100.d0*autocm)**betaDisc  &
            * (rInner*1.e10)**alphaDisc * &
            (((rOuter*1.e10)**(betaDisc-alphaDisc+2.)-(rInner*1.e10)**(betaDisc-alphaDisc+2.))) ))

    case("modular")

       oneKappa = .true.
! -- Define the disc-module independent geometry-specific parameters:
       call getReal("radius1", rCore, real(rsol/1.e10), cLine, fLine, nLines, &
            "Core radius (solar radii): ","(a,f7.3,a)", 10., ok, .true.)
!      Stellar radius (as for shakara)

       call getReal("teff1", teff, 1., cLine, fLine, nLines, &
            "Source effective temperature: ","(a,f7.1,a)", 12., ok, .true.)
!      Stellar effective temperature (as for shakara)

       call getReal("rsub", rSublimation, rCore, cLine, fLine, nLines, &
            "Sublimation radius (rstar): ","(a,f5.1,a)", rCore*(1600./teff)**(-2.1), ok, .true.)
!      Sublimation radius (defaults to Whitney+04 prescription if not defined)

       call getReal("mass1", mCore, real(msol), cLine, fLine, nLines, &
            "Core mass (solar masses): ","(a,f8.4,a)", 0.5, ok, .true.)
!      Stellar mass (as for shakara)

       call getReal("mdisc", mDisc, real(msol), cLine, fLine, nLines, &
            "Mass of disc contained in gas (solar masses): ","(a,f8.4,a)", 1.e-4, ok, .true.)
!      Disc (gas plus dust) mass (as for shakara)

! -- Define the AMR parameters
       call getLogical("smoothgridtau", doSmoothGridtau, cLine, fLine, nLines, &
            "Smooth AMR grid using tau: ","(a,1l,1x,a)", .false., ok, .true.)

       call getLogical("dosmoothgrid", doSmoothGrid, cLine, fLine, nLines, &
            "Smooth AMR grid: ","(a,1l,1x,a)", .false., ok, .true.)

       call getReal("heightsplitfac", heightSplitFac, 1., cLine, fLine, nLines, &
            "Splitting factor for scale height (local scale heights): ","(a,f5.2,a)", 0.2, ok, .true.)

! -- Count the number of disc modules
       call getInteger("ndiscmodule", nDiscModule, cLine, fLine, nLines, &
            "Number of different disc modules: ","(a,i12,a)", 1, ok, .true.)

! -- Count the number of dust species and check whether dust settling is specified
       call getInteger("ndusttype", nDustType, cLine, fLine, nLines, &
            "Number of different dust types: ","(a,i12,a)", 1, ok, .true.)

       call getLogical("dustsettling", dustSettling, cLine, fLine, nLines, &
            "Dust settling model: : ","(a,1l,1x,a)", .false., ok, .false.)

! -- Define the geometry parameters specific to each disc module:
       do i = 1, nDiscModule

          write(alphaDiscLabel, '(a,i1.1)') "tiltangle",i
          call getDouble(alphaDiscLabel, tiltAngleMod(i), degtoRad, cLine, fLine, nLines, &
               "Tilt angle for disc module: ","(a,f8.3,a)", 0.d0, ok, .false.)
!         tilt angle for each module of the disc

          write(alphaDiscLabel, '(a,i1.1)') "alphamod",i
          call getDouble(alphaDiscLabel, alphaMod(i), 1.d0, cLine, fLine, nLines, &
               "Alpha parameter for disc module: ","(a,f8.3,a)", 2.25d0, ok, .true.)
!         alpha parameter for each module of the disc

          write(betaDiscLabel, '(a,i1.1)') "betamod",i
          call getDouble(betaDiscLabel, betaMod(i), 1.d0, cLine, fLine, nLines, &
               "Beta parameter for disc module: ","(a,f8.3,a)", 1.25d0, ok, .true.)
!         beta parameter for each module of the disc

          write(hDiscLabel, '(a,i1.1)') "heightmod",i
          call getDouble(hDiscLabel, heightMod(i), autocm/1.d10, cLine, fLine, nLines, &
               "Gas scale height at inner radius of disc module (au): ","(a,1pe8.2,a)", 1.d0, ok, .true.)
!         scale height at inner radius of each module of the disc.

          write(rinDiscLabel, '(a,i1.1)') "rinnermod",i
          call getDouble(rinDiscLabel, rInnerMod(i), autocm/1.d10, cLine, fLine, nLines, &
               "Inner radius of disc module (au): ","(a,f7.3,a)", 12.d0, ok, .true.)
!         inner radius of each module of the disc

          write(routDiscLabel, '(a,i1.1)') "routermod",i
          call getDouble(routDiscLabel, rOuterMod(i), autocm/1.d10, cLine, fLine, nLines, &
               "Outer radius of disc module (au): ","(a,f7.3,a)", 20.d0, ok, .true.)
!         outer radius of each module of the disc

          do j = 1, nDustType
             write(dustFracLabel, '(a,i1.1,i1.1)') "dustfrac",i,j
             call getDouble(dustFracLabel, dustFracMod(i,j), 1.d0, cLine, fLine, nLines, &
                  "Fraction of dust species in disc module :","(a,f10.5,1x,a)", dble(0.01), ok, .true.)
!            fraction of each dust species in each modular section of the disc

             write(betaLabel, '(a,i1.1,i1.1)') "settlebeta",i,j
             call getDouble(betaLabel, dustBetaMod(i,j), 1.d0, cLine, fLine, nLines, &
                  "Dust beta parameter for disc module: ","(a,f10.5,1x,a)", betaMod(i), ok, .false.)
!            beta parameter for each dust prescription (default is that of the gas in that disc module)

             write(heightLabel, '(a,i1.1,i1.1)') "settleheight",i,j
             call getDouble(heightLabel, dustHeightMod(i,j), autocm/1.d10, cLine, fLine, nLines, &
                  "Dust scale height at inner radius of disc module (au): ","(a,f10.5,1x,a)", &
                  heightMod(i), ok, .false.)
!            scale height for each dust prescription (default is that of the gas in that disc module)

          enddo
       if (writeoutput) write(*,*)
       enddo
!      Calculate the density normalisation assuming that the midplane density at r_out,i which is prescribed by disk
!      module 'i' is the same as that at r_in,i+1 which is prescribed by the next adjacent disk module, 'i+1'.
       prod(1) = 1.d0
       if (nDiscModule > 1) then
          do i = 2, nDiscModule
             prod(i) = prod(i-1)*((rInnerMod(i)*1.d10)**(alphaMod(i)-alphaMod(i-1)))
          enddo
       endif
       do i = 1, nDiscModule
          rhoNought(i) = heightMod(i)*1.d10*(rInnerMod(i)*1.d10)**(betaMod(i)*(-1.d0))*prod(i)* &
                         ( ((rOuterMod(i)*1.d10)**(betaMod(i)-alphaMod(i)+2.d0) - &
                         (rInnerMod(i)*1.d10)**(betaMod(i)-alphaMod(i)+2.d0) ) / &
                         (betaMod(i)-alphaMod(i)+2.d0) )
!         density integration constant for each module of the disc
       enddo

       rho0  = dble( (mDisc / twoPi**1.5) * (1/ERF(1.d0/SQRT(2.d0))) * (rInnerMod(1)*1.d10)**((-1.d0)*alphaMod(1)) * &
               (1/SUM(rhoNought)) )
!      density integration constant for full disc

    end select
  end subroutine readGeometrySpecificParameters




  subroutine readGridInitParameters(cLine, fLine, nLines)
    use octal_mod, only: cart2d
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok


    gridUsesAMR = .true.

    call getLogical("setupamr", doSetupAMRgrid, cLine, fLine, nLines, &
         "Set up or readin and AMR grid: ","(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("cylindrical", cylindrical, cLine, fLine, nLines, &
         "Grid uses 3D cylindical  coords: ","(a,1l,1x,a)", .false., ok, .false.)

    call getDouble("dphiref", dphiRefine, 1.d0, cLine, fLine, nLines, &
         "Level of azimuthal refinement (degrees): ","(a,f5.1,a)", 10.d0, ok, .false.)


    call getDouble("minphi", minPhiResolution, degtorad, cLine, fLine, nLines, &
         "Level of azimuthal refinement (degrees): ","(a,f5.1,a)", 1.d30, ok, .false.)



    call getLogical("amr1d", amr1d, cLine, fLine, nLines, &
         "AMR grid is in one-dimensions only: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("cart2d", cart2d, cLine, fLine, nLines, &
         "2d RHD in cartesians: ","(a,1l,1x,a)", .false., ok, .false.)

    spherical = .false.
    if (amr1d) then
       call getLogical("spherical", spherical, cLine, fLine, nLines, &
            "AMR grid is one-dimension and spherical: ","(a,1l,1x,a)", .true., ok, .false.)
    endif

    call getLogical("amr2d", amr2d, cLine, fLine, nLines, &
         "AMR grid is in two-dimensions only: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("amr3d", amr3d, cLine, fLine, nLines, &
         "AMR grid is in three-dimensions: ","(a,1l,1x,a)", .false., ok, .false.)

    if ((.not.(amr1d.or.amr2d.or.amr3d)).and.doSetupAMRgrid) then
       call writeFatal("AMR dimensionality must be specified")
       call torus_stop
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

    call getLogical("dorefine", dorefine, cLine, fLine, nLines, &
         "Adaptively refine AMR grid?: ","(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("dounrefine", dounrefine, cLine, fLine, nLines, &
         "Adaptively unrefine AMR grid?: ","(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("dophotorefine", dophotorefine, cLine, fLine, nLines, &
         "Adaptively refine AMR grid in photoloop?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("refineonmass", refineonmass, cLine, fLine, nLines, &
         "Refine grid using mass in cell?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("refineondensity", refineondensity, cLine, fLine, nLines, &
         "Refine on density gradient?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getDouble("densitygradient", maxDensityGradient, 1.d0 , cLine, fLine, nLines, &
         "Maximum density gradient allowed before AMR grid refines: ","(a,1p,e9.3,1x,a)", 0.1d0, ok, .false.)

    call getLogical("refineonjeans", refineonJeans, cLine, fLine, nLines, &
         "Refine grid using local jeans mass?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getUnitDouble("masstol", masstol, "mass" , cLine, fLine, nLines, &
         "Cell mass tolerance: ","(a,1p,e9.3,1x,a)", 1.d-5*mSol, ok, .false.)

    call getLogical("refineontemperature", refineontemperature, cLine, fLine, nLines, &
         "Refine grid using temperature gradient?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("refineonrhoe", refineonrhoe, cLine, fLine, nLines, &
         "Refine grid using rhoe?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("refineonspeed", refineonspeed, cLine, fLine, nLines, &
         "Refine grid using speed?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("refineonionization", refineonionization, cLine, fLine, nLines, &
         "Refine grid using ionization gradient?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("captureshocks", captureShocks, cLine, fLine, nLines, &
         "Capture shocks?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getDouble("amrtolerance", amrtolerance, 1.d0 , cLine, fLine, nLines, &
         "Maximum gradient allowed before AMR grid refines: ","(a,1p,e9.3,1x,a)", 5.d-2, ok, .false.)

    call getDouble("unreftol", amrUnrefinetolerance, 1.d0 , cLine, fLine, nLines, &
         "Minimum gradient allowed before AMR grid unrefines: ","(a,1p,e9.3,1x,a)", 1.d-2, ok, .false.)

    if(refineontemperature) then
       call getDouble("amrtemperaturetol", amrTemperatureTol, 1.d0 , cLine, fLine, nLines, &
            "Maximum temperature gradient allowed before AMR grid refines: ","(a,1p,e9.3,1x,a)", 1.d-1, ok, .false.)
    end if

    if(refineonspeed) then
       call getDouble("amrspeedtol", amrSpeedTol, 1.d0 , cLine, fLine, nLines, &
            "Maximum speed gradient allowed before AMR grid refines: ","(a,1p,e9.3,1x,a)", 1.d-1, ok, .false.)
    end if

    if(refineonionization) then
       call getDouble("amrionfractol", amrIonFracTol, 1.d0 , cLine, fLine, nLines, &
            "Maximum ionization fraction gradient allowed before AMR grid refines: ","(a,1p,e9.3,1x,a)", 1.d-1, ok, .false.)
    end if

    if(refineonrhoe) then
       call getDouble("amrrhoetol", amrRhoeTol, 1.d0 , cLine, fLine, nLines, &
            "Maximum rhoe gradient allowed before AMR grid refines: ","(a,1p,e9.3,1x,a)", 1.d-1, ok, .false.)
    end if


!    call getReal("amrgridsize", amrGridSize, 1.0, cLine, fLine, nLines, &
!         "Size of adaptive mesh grid: ","(a,1pe8.1,1x,a)", 1000.0, ok, .false.)

    call getRealwithUnits("amrgridsize", amrGridSize, "codeunits", "codeunits", cLine, fLine, nLines, &
         "Size of adaptive mesh grid: ", 1000.0, ok, .false.)
    call getDoubleWithUnits("amrgridcentrex", amrGridCentreX, "codeunits" , "codeunits", cLine, fLine, nLines, &
         "Grid centre X-coordinate: ",0.0d0, ok, .false.)
    call getDoubleWithUnits("amrgridcentrey", amrGridCentreY, "codeunits" , "codeunits", cLine, fLine, nLines, &
         "Grid centre Y-coordinate: ", 0.0d0, ok, .false.)
    call getDoubleWithUnits("amrgridcentrez", amrGridCentreZ, "codeunits", "codeunits",  cLine, fLine, nLines, &
         "Grid centre Z-coordinate: ", 0.0d0, ok, .false.)

    if (amr2d.and.(.not.checkPresent("amrgridcentrex", cline, nlines))) &
         amrGridCentrex = dble(amrGridSize)/2.d0

    if (amr3d.and.cylindrical.and.(.not.checkPresent("amrgridcentrex", cline, nlines))) &
         amrGridCentrex = amrGridSize/2.d0

    if (amr1d.and.(.not.checkPresent("amrgridcentrex", cline, nlines))) &
         amrGridCentrex = amrGridSize/2.d0

    if (spherical) then
       amrGridCentreX = amrgridsize/2.
    endif

    write(*,*) "amrgridsize ",amrgridsize,amrgridcentrex,amrgridcentrey,amrgridcentrez
    call getDouble("limitscalar", limitScalar, 1.d0, cLine, fLine, nLines, &
         "Scalar limit for subcell division: ","(a,1p,e9.3,1x,a)", 1000._db, ok, .false.)

    call getDouble("limittwo", limitScalar2, 1.d0, cLine, fLine, nLines, &
         "Second scalar limit for subcell division: ","(a,1p,e9.3,1x,a)", 0._db, ok, .false.)

    call getLogical("doVelocitySplit", doVelocitySplit, cLine, fLine, nLines, &
         "Use velocity splitting condition for SPH to grid: ","(a,1l,1x,a)",.false., ok, .false.)

    call getString("geometry", geometry, cLine, fLine, nLines, &
         "Geometry: ","(a,a,a)","none",ok, .false.)

! For setting up a grid from a Flash HDF5 dump
    if (checkPresent("flashfilename", cline, nlines)) call readFlashParameters

! For setting up a grid from a VH-1 dump
    if (checkPresent("vh1filename", cline, nlines)) call readVh1Parameters

#ifdef USECFITSIO
! For setting up a grid from a FITS file
     if (checkPresent("fitsgridfile", cline, nlines)) call readFitsGridParameters
#endif

! For setting up a grid from a Ramses file
     if (checkPresent("ramsesfilename", cline, nlines)) call readRamsesParameters

  contains

    subroutine readFlashParameters
      use gridFromFlash, only: setGridFromFlashParameters

      character(len=80) :: flashfilename
      integer :: numBlocks
      real(double) :: slice
      logical :: doReflectY

      call getString("flashfilename", flashfilename, cLine, fLine, nLines, &
           "Flash file name: ","(a,a,a)","default",ok,.true.)

      call getInteger("flashnumblocks", numblocks, cLine, fLine, nLines, &
           "Number of blocks in flash file: ", "(a,i7,a)",1000,ok,.true.)

      call getDouble("flashslice", slice, 1.d0, cLine, fLine, nLines, &
           "Position of slice for 2D cut ","(a,1p,e9.3,1x,a)", 0.8e8_db, ok, .false.)

      call getLogical("flashreflecty", doReflectY, cLine, fLine, nLines, &
           "Reflect y=axis: ","(a,1l,1x,a)", .false., ok, .false.)

      call setGridFromFlashParameters(flashfilename, numblocks, slice, doReflectY)

    end subroutine readFlashParameters

    subroutine readVh1Parameters
      use vh1_mod, only: setGridFromVh1Parameters

      character(len=80) :: vh1filename
      real(double) :: vh1xoffset, vh1yoffset

      call getString("vh1filename", vh1filename, cLine, fLine, nLines, &
           "VH1: file name: ","(a,a,a)","default",ok, .true.)

      call getDouble("vh1xoffset", vh1xoffset, 1.d0, cLine, fLine, nLines, &
           "VH1: offset applied to x-axis: ","(a,1p,e10.3,1x,a)", 0.0e0_db, ok, .false.)

      call getDouble("vh1yoffset", vh1yoffset, 1.d0, cLine, fLine, nLines, &
           "VH1: offset applied to y-axis: ","(a,1p,e10.3,1x,a)", 0.0e0_db, ok, .false.)

      call setGridFromVh1Parameters(vh1filename, vh1xoffset, vh1yoffset)

    end subroutine readVh1Parameters

#ifdef USECFITSIO
    subroutine readFitsGridParameters
      use gridFromFitsFile, only: setGridFromFitsParameters, setGridFromFitsParametersPionAMR
      character(len=80) :: filename
      character(len=20) :: fileLabel
      integer :: i

      call getLogical("pionAMR", pionAMR, cLine, fLine, nLines, &
           "Pion grid is non-uniform?: : ","(a,1l,1x,a)", .false., ok, .false.)

      if(pionAMR) then
         call getInteger("nFitsLevels", nFitsLevels, cLine, fLine, nLines,"Number of fits levels: ","(a,i12,a)",1,ok,.false.)

         do i = 0, nFitsLevels-1
            write(fileLabel, '(a,i1.1)') "fitsgridfile",i

            call getString(fileLabel, pionFitsAMRfile(i+1), cLine, fLine, nLines, &
                 "FITS file for making grid: ","(a,a,a)","default",ok, .true.)

         enddo
         call setGridFromFitsParametersPionAMR(nFitsLevels, pionFitsAMRfile,amr2d,amr3d, pionAMR,maxdepthamr)
      else

         call getString("fitsgridfile", filename, cLine, fLine, nLines, &
              "FITS file for making grid: ","(a,a,a)","default",ok, .true.)
         call setGridFromFitsParameters(filename,amr2d,amr3d, pionAMR)
      endif



    end subroutine readFitsGridParameters
#endif

    subroutine readRamsesParameters
      use ramses_mod, only: setGridFromRamsesParameters
      character(len=80) :: ramsesfilename

      call getString("ramsesfilename", ramsesfilename, cLine, fLine, nLines, &
           "Ramses file name: ","(a,a,a)","none",ok, .false.)

      call setGridFromRamsesParameters(ramsesfilename)
    end subroutine readRamsesParameters


  end subroutine readGridInitParameters

  subroutine readTimeDependentRTParameters(cLine, fLine, nLines)
    use sed_mod, only:  setSedParameters,  SEDlamMin, SEDlamMax, SEDnumLam
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call writeBanner("Time dependendent RT",".",TRIVIAL)

    call getDouble("timestart", timeStart, yearsToSecs, cLine, fLine, nLines, &
         "Start time of calculation (years):  ","(a,e12.3,1x,a)", 0.d0, ok, .false.)

    call getDouble("timeend", timeEnd, yearsToSecs, cLine, fLine, nLines, &
         "End time of calculation (years):  ","(a,e12.3,1x,a)", 1500.d0, ok, .false.)

    call getDouble("varystart", varystart, yearsToSecs, cLine, fLine, nLines, &
         "Start time of source variability (years):  ","(a,e12.3,1x,a)", 1500.d0, ok, .true.)

    call getDouble("varyend", varyend, yearsToSecs, cLine, fLine, nLines, &
         "End time of source variability (years):  ","(a,e12.3,1x,a)", 1.d30, ok, .false.)


    call getDouble("lumfactor", lumFactor, 1.d0, cLine, fLine, nLines, &
         "Increase in luminosity as a factor of standard luminosity:  ","(a,e12.3,1x,a)", 1500.d0, ok, .true.)

    call getDouble("lumdecaytime", lumDecayTime, yearsToSecs, cLine, fLine, nLines, &
         "Luminosity decay timescale (e-folding time in years):  ","(a,e12.3,1x,a)", 1500.d0, ok, .true.)


    call getBigInteger("nphotons", nPhotons, cLine, fLine, nLines,"Number of photons per timestep: ",&
         "(a,i12,a)",1_bigint,ok,.false.)

    call getLogical("quicksublimate", quickSublimate, cLine, fLine, nLines, &
         "Quick and dirty dust sublimation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getInteger("ntime", nTime, cLine, fLine, nLines,"Number of timesteps: ","(a,i12,a)",1,ok,.false.)

    call getReal("inclination", thisinclination, real(degtorad), cLine, fLine, nLines, &
         "Inclination angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)

    call getString("inputfile", gridInputFilename, cLine, fLine, nLines, &
         "Grid input filename: ","(a,a,1x,a)","none", ok, .true.)

    call getReal("distance", gridDistance, real(pctocm), cLine, fLine, nLines, &
         "Grid distance (pc): ","(a,f6.1,1x,a)", 100., ok, .true.)

    call getReal("sedlammin", SEDlamMin, 1.0e4, cLine, fLine, nLines, &
         "Minimum wavelength output to SED (microns)","(a,1PE10.3,1x,a)", 0.1, ok, .false.)

    call getReal("sedlammax", SEDlamMax, 1.0e4, cLine, fLine, nLines, &
         "Maximum wavelength output to SED (microns)","(a,1PE10.3,1x,a)", 2000.0, ok, .false.)


    call getInteger("sednumlam", SEDnumLam, cLine, fLine, nLines, &
         "Number of SED points: ", "(a,i3,1x,a)", 50, ok, .false.)

    call getInteger("nsedradius", nSedRadius, cLine, fLine, nLines, &
         "Number of SED apertures: ", "(a,i3,1x,a)", 1, ok, .false.)

    call getDoubleArray("sedradius", sedRadius(1:nSedRadius), autocm/1.d10, cLine, fLine, nLines, &
         "Radius of aperture for SED (au): ", 1.d30, ok, .false.)

    
    call getInteger("npix", npix, cLine, fLine, nLines, &
         "Number of npixels across image: ", "(a,i3,1x,a)", 50, ok, .true.)

      call getLogical("timedepimage", timedepImage, cLine, fLine, nLines, &
           "Calculate time dependent images: ","(a,1l,1x,a)", .false., ok, .false.)

      if (timedepImage) then
         call getDouble("imagesize", imageSize, autocm/1.d10, cLine, fLine, nLines, &
              "Image size (au):  ","(a,f5.1,1x,a)", 1500.d0, ok, .true.)

         call getDouble("lamstart", lamStart, 1.d4, cLine, fLine, nLines, &
              "Start wavelength for image (microns):  ","(a,f5.1,1x,a)", 1500.d0, ok, .true.)

         call getDouble("lamend", lamEnd, 1.d4, cLine, fLine, nLines, &
              "End wavelength (microns):  ","(a,f5.1,1x,a)", 1500.d0, ok, .true.)

         call getDouble("imagesize", imageSize, autocm/1.d10, cLine, fLine, nLines, &
              "Image size (au):  ","(a,f5.1,1x,a)", 1500.d0, ok, .true.)
      endif

  end subroutine readTimeDependentRTParameters



  subroutine readDustPhysicsParameters(cLine, fLine, nLines)
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok
    real :: grainFracTotal
    integer :: i
    character(len=20) :: grainTypeLabel, grainFracLabel, aMinLabel, grainDensityLabel, &
         aMaxLabel, qDistLabel, pdistLabel, a0label, fillingFactorLabel, tsubLabel, kappaFileLabel, &
         dustFileLabel, heightLabel, betaLabel, tsubPowerLabel

       oneKappa = .true.

       call writeBanner("Dust parameters","#",TRIVIAL)


       call writeBanner("NOTE THAT DUST TO GAS RATIO INFORMATION SHOULD BE PLACED IN GRAINFRAC. ","!",TRIVIAL)
       call writeBanner("DUSTTOGAS IS NOW USED TO DENOTE THE TOTAL DUST MASS IN THE SYSTEM.","!",TRIVIAL)

       call getReal("dusttogas", dusttoGas, 1., cLine, fLine, nLines, &
            "Dust to gas ratio: ","(a,f5.3,a)",0.01,ok,.false.)


       call getLogical("dustsettling", dustSettling, cLine, fLine, nLines, &
               "Dust settling model: : ","(a,1l,1x,a)", .false., ok, .false.)


       call getLogical("pah", usePAH, cLine, fLine, nLines, &
               "Include PAH/VSG physics : ","(a,1l,1x,a)", .false., ok, .false.)

       if (usePAH) then
          call getString("pahtype", pahType, cLine, fLine, nLines, &
               "PAH/VSG dust type: ","(a,a,1x,a)","MW3.1_60", ok, .true.)
          call getDouble("pahscale", PAHscale, 1.d0, cLine, fLine, nLines, &
               "Scaling factor for PAH emissivities and opacities: ","(a,f7.3,a)",1.d0, ok, .false.)
       endif

       call getInteger("ndusttype", nDustType, cLine, fLine, nLines,"Number of different dust types: ","(a,i12,a)",1,ok,.false.)
       if (nDustType .gt. maxDustTypes) then
          if (writeoutput) write (*,*) "Max dust types exceeded: ", maxDustTypes
          call torus_stop
       end if

       call getLogical("readdust", readDustFromFile, cLine, fLine, nLines, &
         "Read dust opacities and phase matrices from file: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("writedust", writeDustToFile, cLine, fLine, nLines, &
         "Write dust opacities and phase matrices to file: ","(a,1l,1x,a)", .false., ok, .false.)


       call getLogical("writemiephase", writeMiePhase, cLine, fLine, nLines, &
            "Write mie scattering phase file: ","(a,1l,1x,a)",.false., ok, .false.)

       call getLogical("readmiephase", readMiePhase, cLine, fLine, nLines, &
            "Read mie scattering phase file: ","(a,1l,1x,a)",.false., ok, .false.)

       call getReal("subrange", subrange, 1., cLine, fLine, nLines, &
            "Sublimation temperature e-folding temperature: ","(a,f5.3,a)",10.,ok,.false.)


          grainFracTotal = 0.
          do i = 1, nDustType
             write(dustFileLabel, '(a,i1.1)') "dustfile",i
             write(kappaFileLabel, '(a,i1.1)') "kappafile",i

             write(grainTypeLabel, '(a,i1.1)') "graintype",i
             write(grainFracLabel, '(a,i1.1)') "grainfrac",i
             write(tsubLabel, '(a,i1.1)') "tsub",i
             write(tsubPowerLabel, '(a,i1.1)') "tsubpower",i
             write(grainDensityLabel, '(a,i1.1)') "graindensity",i
             write(aMinLabel, '(a,i1.1)') "amin",i
             write(aMaxLabel, '(a,i1.1)') "amax",i
             write(qDistLabel, '(a,i1.1)') "qdist",i
             write(pDistLabel, '(a,i1.1)') "pdist",i
             write(a0Label, '(a,i1.1)') "a0",i
             write(fillingFactorLabel, '(a,i1.1)') "porosity",i






          !       if (writeoutput) write(*,'(a,i1.1)') "Dust properties for grain ",i
             !       if (writeoutput) write(*,'(a,i1.1)') "-------------------------------"
             !       if (writeoutput) write(*,*)

             call getLogical(dustFileLabel, dustFile(i), cLine, fLine, nLines, &
                  "Read dust from kappa file: ","(a,a,1x,a)",.false., ok, .false.)

             if (dustFile(i)) then
                call getString(kappaFileLabel, kappaFilename(i), cLine, fLine, nLines, &
                     "Kappa file: ","(a,a,1x,a)","test", ok, .true.)
             endif

             if (dustfile(i)) then
                call getReal(grainFracLabel, grainFrac(i), 1., cLine, fLine, nLines, &
                     "Grain fractional abundance: ","(a,f8.5,1x,a)",0.01 , ok, .false.)
                grainFracTotal = grainFracTotal + grainFrac(i)
             endif


             if (.not.dustFile(i)) then

             call getString(grainTypeLabel, grainType(i), cLine, fLine, nLines, &
                  "Grain type: ","(a,a,1x,a)","sil_dl", ok, .true.)

             call getReal(grainFracLabel, grainFrac(i), 1., cLine, fLine, nLines, &
                  "Grain fractional abundance: ","(a,f8.5,1x,a)",0.01 , ok, .false.)
             grainFracTotal = grainFracTotal + grainFrac(i)

             call getReal(grainDensityLabel, grainDensity(i), 1., cLine, fLine, nLines, &
                  "Grain density (g/cc): ","(a,f8.5,1x,a)",3.6 , ok, .false.)

             call getReal(fillingFactorLabel, porousFillingFactor(i), 1., cLine, fLine, nLines, &
                  "Porosity: ","(a,f5.2,1x,a)", 0. , ok, .false.)


             call getDouble(tsublabel, tsub(i), 1.d0, cLine, fLine, nLines, &
                  "Temperature for dust sublimation (K):  ","(a,e12.3,1x,a)", 1500.d0, ok, .false.)

             call getDouble(tsubpowerlabel, tsubpower(i), 1.d0, cLine, fLine, nLines, &
                  "Density power-law index  for dust sublimation (K):  ","(a,e12.3,1x,a)", 0.d0, ok, .false.)

             if (.not. readDustFromFile) &
                  call getRealwithUnits(aminLabel, aMin(i), "micron", "micron", cLine, fLine, nLines, &
                  "Min grain size: ", 1000.0, ok, .true.)

             if (.not. readDustFromFile) &
                  call getRealwithUnits(amaxLabel, aMax(i), "micron", "micron", cLine, fLine, nLines, &
                  "Min grain size: ", 1000.0, ok, .true.)


             if (.not. readDustFromFile) &
                  call getReal(qDistLabel, qdist(i), 1., cLine, fLine, nLines, &
                  "Grain power law: ","(a,f4.1,1x,a)", 3.5, ok, .true. )

             if (.not. readDustFromFile) &
                  call getReal(a0Label, a0(i), 1., cLine, fLine, nLines, &
                  "Scale length of grain size (microns): ","(a,f8.5,1x,a)", 1.0e20, ok, .false.)


             if (.not. readDustFromFile) &
          call getReal(pdistLabel, pdist(i), 1., cLine, fLine, nLines, &
          "Exponent for exponential cut off: ","(a,f4.1,1x,a)", 1.0, ok, .false. )



          write(heightLabel, '(a,i1.1)') "dustheight",i
          call getDouble(heightLabel, dustHeight(i), autocm/1.d10, cLine, fLine, nLines, &
               "Dust scale height at 100 AU (AU): ","(a,f10.5,1x,a)", dble(height)*1.d10/autocm, ok, .false.)

          write(betaLabel, '(a,i1.1)') "dustbeta",i
          call getDouble(betaLabel, dustBeta(i), 1.d0, cLine, fLine, nLines, &
               "Dust beta law (AU): ","(a,f10.5,1x,a)", dble(betaDisc), ok, .false.)

          endif
             if (writeoutput) write(*,*)
          enddo



          call getLogical("iso_scatter", isotropicScattering, cLine, fLine, nLines, &
               "Isotropic scattering: ","(a,1l,1x,a)", .false., ok, .false.)

         call getLogical("henyey", henyeyGreensteinPhaseFunction, cLine, fLine, nLines, &
              "Use Henyey-Greenstein phase function: ","(a,1l,1x,a)", .false., ok, .false.)

         if (henyeyGreensteinPhaseFunction) &
         call getDouble("gfac", inputgFac, 1.d0, cLine, fLine, nLines, &
            "Henyey-Greenstein g-factor: ","(a,f7.3,a)",0.d0, ok, .false.)


         call getLogical("writepolar", writePolar, cLine, fLine, nLines, &
              "Write polarizability file: ","(a,1l,1x,a)", .false., ok, .false.)
         if (writepolar) then
            call getReal("polarwav", polarWavelength, real(micronsToAngs), cLine, fLine, nLines, &
           "Polarizability wavelength (microns): ","(a,f7.2,1x,a)", 0.0, ok, .true.)
            write(*,*) "p ",polarwavelength
             call getString("polarfile", polarFilename, cLine, fLine, nLines, &
                  "Polarizability filename: ","(a,a,1x,a)","sil_dl", ok, .true.)
         endif
  end subroutine readDustPhysicsParameters

  subroutine readAtomicPhysicsParameters(cLine, fLine, nLines)
    use atom_mod, only: setVoigtParams
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    integer :: i
    character(len=20) :: keyword
    logical :: ok
    real :: C_rad, C_vdw, C_stark

    call writeBanner("Atomic physics data","#",TRIVIAL)

    call getLogical("sei", sei, cLine, fLine, nLines, &
         "Perform Sobolev with Exact integration ","(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("cmf", cmf, cLine, fLine, nLines, &
         "Perform co-moving frame calculation ","(a,1l,1x,a)", .false., ok, .true.)


    ! if (cmf) then
    !    call getInteger("natom", nAtom, cLine, fLine, nLines, &
    !         "Number of model atoms to solve for: ","(a,i1,a)",1, ok, .true.)
    !    do i = 1, nAtom
    !       write(keyword, '(a,i1)') "atom",i
    !       call getString(keyword, atomFileName(i), cLine, fLine, nLines, &
    !            "Use atom filename: ","(a,a,1x,a)","none", ok, .true.)
    !    enddo
    ! endif

    call getInteger("natom", nAtom, cLine, fLine, nLines, &
            "Number of model atoms to solve for: ","(a,i1,a)",1, ok, .true.)
       do i = 1, nAtom
          write(keyword, '(a,i1)') "atom",i
          call getString(keyword, atomFileName(i), cLine, fLine, nLines, &
               "Use atom filename: ","(a,a,1x,a)","none", ok, .true.)
       enddo

    call getLogical("thickcont", opticallyThickContinuum, cLine, fLine, nLines, &
         "Continuum is optically thick: ","(a,1l,a)", .false., ok, .false.)

    call getDouble("xabundance", Xabundance, 1.d0, cLine, fLine, nLines, &
            "Hydrogen abundance (by number): ","(a,f7.3,a)",1.d0, ok, .true.)

    call getDouble("yabundance", Yabundance, 1.d0, cLine, fLine, nLines, &
            "Helium abundance (by number): ","(a,f7.3,a)",1.d0, ok, .true.)

    call getReal("lamline", lamLine, 1.,cLine, fLine, nLines, &
               "Line emission wavelength: ","(a,f8.1,1x,a)", 850., ok, .true.)

    call getReal("vturb", vturb, real(kmstoc), cLine, fLine, nLines, &
               "Turbulent velocity (km/s):","(a,f6.1,1x,a)", 0., ok, .true.)

    call getLogical("starkbroaden",starkBroaden, cline, fline, nlines, &
                "Use Stark Broadening", "(a,1l,1x,a)", .false., ok, .true.)

    !
    ! Voigt profile prameters
    !
    IF (starkBroaden) THEN
      call getReal("C_rad", C_rad, 1., cLine, fLine, nLines, &
           "Damping constant (radiation)     in [A]: ","(a,1PE10.3,1x,a)", 8.16e-3, ok, .false.)
      call getReal("C_vdw", C_vdw, 1., cLine, fLine, nLines, &
           "Damping constant (van der Waals) in [A]: ","(a,1PE10.3,1x,a)", 5.53e-3, ok, .false.)
      call getReal("C_stark", C_stark, 1., cLine, fLine, nLines, &
           "Damping constant (Stark)         in [A]: ","(a,1PE10.3,1x,a)", 1.47e-2, ok, .false.)
      call setVoigtParams(C_rad, C_vdw, C_stark)
    END IF


  end subroutine readAtomicPhysicsParameters

  subroutine readSourceParameters(cLine, fLine, nLines)
    character(len=lencLine) :: cLine(:)
    character(len=80) :: message
    logical :: fLine(:)
    integer :: nLines
    integer :: i
    character(len=20) :: keyword
    logical :: ok
    logical :: teffPresent, radiusPresent, lumPresent





    call writeBanner("Photon source data","#",TRIVIAL)



    call getInteger("nspheresurface", nSphereSurface, cLine, fLine, nLines, &
         "Number of points on sphere surface: ","(a,i2,a)",100,ok,.false.)


    call getInteger("nsource", inputNSource, cLine, fLine, nLines, &
         "Number of sources: ","(a,i2,a)",0,ok,.false.)


    do i = 1, inputnSource
       if (writeoutput) write(*,*) " "
       write(message,'(a,i1)') "Source number: ",i
       call writeInfo(message,TRIVIAL)

       write(keyword, '(a,i1)') "stellar",i
       call getLogical(keyword, stellarSource(i), cLine, fLine, nLines, &
            "Stellar source: ","(a,xxxx,a)",.true., ok, .false.)

       if (doNbodyOnly) then
          write(keyword, '(a,i1)') "sourcepos",i
          call getVector(keyword, sourcePos(i), 1.d0, cLine, fLine, nLines, &
               "Source position: ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

          write(keyword, '(a,i1)') "velocity",i
          call getVector(keyword, sourceVel(i), 1.d0, cLine, fLine, nLines, &
               "Source velocity: ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .false.)
          write(keyword, '(a,i1)') "mass",i
          call getDoubleWithUnits(keyword, sourceMass(i), "msol", "gram", cLine, fLine, nLines, &
               "Source mass: ",1.d0, ok, .true.)


       else

          if (stellarSource(i)) then
             write(keyword, '(a,i1)') "radius",i
             radiusPresent = checkPresent(keyword, cline, nlines)
             write(keyword, '(a,i1)') "teff",i
             teffPresent = checkPresent(keyword, cline, nlines)
             write(keyword, '(a,i1)') "lum",i
             lumPresent = checkPresent(keyword, cline, nlines)

             if (radiusPresent.and.teffPresent.and.lumPresent) then
                call writeFatal("Source overspecified - choose 2 of teff, R, and Lum")
                call torus_stop
             endif

             if (radiusPresent.and.teffPresent) then
                write(keyword, '(a,i1)') "radius",i
                call getDoubleWithUnits(keyword, sourceRadius(i), "rsol","codeunits", cLine, fLine, nLines, &
                     "Source radius: ",1.d0, ok, .true.)

                write(keyword, '(a,i1)') "teff",i
                call getDoubleWithUnits(keyword, sourceTeff(i), "K", "K",cLine, fLine, nLines, &
                     "Source temperature: ",1.d0, ok, .true.)

                sourceLum(i) = fourPi * sourceRadius(i)**2 *1.d20 * stefanBoltz * sourceTeff(i)**4
             else if (radiusPresent.and.lumPresent) then

                write(keyword, '(a,i1)') "radius",i
                call getDoubleWithUnits(keyword, sourceRadius(i), "rsol", "codeunits", cLine, fLine, nLines, &
                     "Source radius: ",1.d0, ok, .true.)

                write(keyword, '(a,i1)') "lum",i
                call getDouble(keyword, sourceLum(i), lSol, cLine, fLine, nLines, &
                     "Source luminosity (solar luminosities) : ","(a,f7.2,a)",1.d0, ok, .true.)
                sourceTeff(i) = (sourceLum(i) / (fourPi*sourceRadius(i)**2*1.d20 * stefanBoltz))**0.25d0
                if (writeoutput) write(*,'(a,f7.1)') "Source temperature (K): ",sourceTeff(i)
             else if (teffPresent.and.lumPresent) then

                write(keyword, '(a,i1)') "lum",i
                call getDouble(keyword, sourceLum(i), lSol, cLine, fLine, nLines, &
                     "Source luminosity (solar luminosities) : ","(a,1pe12.5,a)",1.d0, ok, .true.)

                write(keyword, '(a,i1)') "teff",i
                call getUnitDouble(keyword, sourceTeff(i), "temperature", cLine, fLine, nLines, &
                     "Source temperature (K) : ","(a,f8.0,a)",1.d0, ok, .true.)

                sourceRadius(i) = sqrt(sourceLum(i) / (fourPi * sourceTeff(i)**4 * stefanBoltz))/1.d10
                if (writeoutput) write(*,'(a,f7.1)') "Source radius (rsol): ",sourceRadius(i)*1.d10/rsol
             else
                if (.not.nbodyPhysics) then
                   call writeFatal("Source underspecified - need two of radius, teff and lum")
                   call torus_stop
                endif
             endif

             write(keyword, '(a,i1)') "mass",i
             call getDoubleWithUnits(keyword, sourceMass(i), "msol", "gram", cLine, fLine, nLines, &
                  "Source mass: ",1.d0, ok, .true.)

             write(keyword, '(a,i1)') "mdot",i
             call getDouble(keyword, sourceMdot(i), msol/(365.25d0 * 24.d0 * 3600.d0), cLine, fLine, nLines, &
                  "Source mass-loss rate (solar masses/year) : ","(a,f9.6,a)",0.d0, ok, .false.)

             write(keyword, '(a,i1)') "contflux",i
             call getString(keyword, inputcontFluxFile(i), cLine, fLine, nLines, &
                  "Continuum flux file: ","(a,a,a)","none", ok, .true.)

             write(keyword, '(a,i1)') "sourcepos",i
             call getVector(keyword, sourcePos(i), 1.d0, cLine, fLine, nLines, &
                  "Source position (10^10 cm): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .false.)

             write(keyword, '(a,i1)') "velocity",i
             call getVector(keyword, sourceVel(i), 1.d5, cLine, fLine, nLines, &
                  "Source velocity (km/s): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .false.)

             write(keyword, '(a,i1)') "probsource",i
             call getDouble(keyword, SourceProb(i), 1.d0, cLine, fLine, nLines, &
                  "Probability of photon packet from source: ","(a,f4.2,a)",0.d0, ok, .false.)

             write(keyword, '(a,i1)') "pointsource",i
             call getLogical(keyword, pointsourcearray(i), cLine, fLine, nLines, &
                  "Point source: ","(a,1l,1x,a)", .false., ok, .false.)

          else
             write(keyword, '(a,i1)') "diffusetype",i
             call getString(keyword, diffuseType(i), cLine, fLine, nLines, &
                  "Type of diffuse radiation field: ","(a,a,1x,a)","none", ok, .true.)
          endif
       endif
    enddo

    call getDouble("biasphidir", biasPhiDirection, degtorad, cLine, fLine, nLines, &
         "Azimuthal direction of photon bias: ","(a,f5.0,a)",-1.d0, ok, .false.)

    if (biasPhiDirection > 0.d0) then

       call getDouble("biasphiprob", biasPhiProb, 1.d0, cLine, fLine, nLines, &
            "Probability of photon in bias direction: ","(a,f5.1,a)",-1.d0, ok, .true.)

       call getDouble("biasphiint", biasPhiInterval, degtorad, cLine, fLine, nLines, &
            "Azimuthal interval for biased direction: ","(a,f5.0,a)",-1.d0, ok, .false.)
    endif


    call getLogical("hotspot", hotSpot, cLine, fLine, nLines, &
         "Add a hotspot to the source: ","(a,1l,1x,a)",.false., ok, .false.)


    call getLogical("pulsating", pulsatingStar, cLine, fLine, nLines, &
         "Star is pulsating: ","(a,1l,1x,a)",.false., ok, .false.)
    if (pulsatingStar) then
       call getInteger("nmodes", nModes,  cLine, fLine, nLines, &
            "Number of modes: ","(a,i2,a)",1, ok, .true.)


       if (.not.allocated(lMode)) allocate(lMode(1:nModes))
       if (.not.allocated(mMode)) allocate(mMode(1:nModes))
       if (.not.allocated(periodMode)) allocate(periodMode(1:nModes))
       if (.not.allocated(fracMode)) allocate(fracMode(1:nModes))


       call getIntegerArray("lmode", lMode, cLine, fLine, nLines, &
            "l mode: ",1, ok, .true.)

       call getIntegerArray("mmode", mMode, cLine, fLine, nLines, &
            "m mode: ",1, ok, .true.)

       call getDoubleArray("period", periodMode, 1.d0, cLine, fLine, nLines, &
            "Pulsation period (seconds): ",1.d0, ok, .true.)

       call getDoubleArray("fracmode", fracMode, 1.d0, cLine, fLine, nLines, &
            "Pulsation period (seconds): ",1.d0, ok, .true.)
    endif



  end subroutine readSourceParameters

  subroutine readMolecularPhysicsParameters(cLine, fLine, nLines)
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

       call getLogical("ALMA", ALMA, cLine, fLine, nLines, &
            "Generate a datacube with ALMA (casa)-friendly headers: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("lte", lte, cLine, fLine, nLines, &
            "Read in LTE grid: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("restart", restart, cLine, fLine, nLines, &
            "Restart molecular calculation : ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("molRestartTest",  molRestartTest, cLine, fLine, nLines, &
            "Save first molecular restart dump for testing : ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("densitysubsample", densitysubsample, cLine, fLine, nLines, &
               "Use density interpolation: ","(a,1l,a)", .false., ok, .false.)
       call getLogical("suppresswarnings", suppressWarnings, cLine, fLine, nLines, &
            "Suppress Warnings: ","(a,l,1x,a)",.false., ok, .false.)
! Not currenlty used
!       call getLogical("addnewmoldata", addnewmoldata, cLine, fLine, nLines, &
!            "Add new molecular data to non-molecular grid: ","(a,l,1x,a)",.false., ok, .false.)
       call getString("moleculefile", moleculefile, cLine, fLine, nLines, &
            "Input molecule filename: ","(a,a,1x,a)","none", ok, .true.)
       call getReal("distance", gridDistance, 1., cLine, fLine, nLines, &
            "Grid distance (pc): ","(a,f8.1,1x,a)", 1., ok, .true.)
       call getInteger("initnray", initnray, cLine, fLine, nLines, &
               "Number of fixed rays for stage 1: ","(a,i4,a)", 1024, ok, .false.)
       call getLogical("dongstep", dongstep, cLine, fLine, nLines, &
               "Use Ng Acceleration: ","(a,1l,a)", .false., ok, .false.)
       call getLogical("quasi", quasi, cLine, fLine, nLines, &
               "Use Quasirandom numbers: ","(a,1l,a)", .false., ok, .false.)
       call getLogical("modelwashydro", modelwashydro, cLine, fLine, nLines, &
               "Grid is based on hydro or radiation hydro calculation: ","(a,1l,a)", .false., ok, .false.)
       !thaw -
       call getDouble("lineunrefinethresh", lineUnrefineThresh, 1.0d0 , cLine, fLine, nLines, &
               "Derefine cells which are below this threshold density prior to line calc: ","(a,f4.1,1x,a)", 0.0d0, ok, .false.)
       if(modelWasHydro) then
          call getLogical("zeroghosts", zeroghosts, cLine, fLine, nLines, &
               "Zero the contriubtion from ghost cells to molecular line calc: ","(a,1l,a)", .false., ok, .false.)
       else
          zeroghosts=.false.
       end if

       call getLogical("renewinputrays", renewinputrays, cLine, fLine, nLines, &
            "Restart with 1000 rays: ","(a,1l,a)", .false., ok, .false.)

       call getLogical("forceiniray", forceIniRay, cLine, fLine, nLines, &
               "Force the first ray towards star in stateq calculations: ","(a,1l,a)", .false., ok, .false.)
       if(forceIniRay) then
          call readSourceParameters(cLine, fLine, nLines)
       end if
!-
       call getReal("tolerance", tolerance, 1., cLine, fLine, nLines, &
            "Maximum Fractional Change in level populations:","(a,f4.1,1x,a)", 0.01, ok, .false.)
       call getReal("vturb", vturb, 1., cLine, fLine, nLines, &
            "Subsonic turbulent velocity (km/s):","(a,f4.1,1x,a)", 0.3, ok, .false.)
       call getLogical("noturb", noturb, cLine, fLine, nLines, &
            "No microturbulence","(a,1l,a)",.false., ok, .false.)
       call getInteger("setmaxlevel", setmaxlevel, cLine, fLine, nLines, &
            "Maximum molecular level to be considered:","(a,i2,1x,a)", 0, ok, .false.)
       call getReal("setmaxleveltemp", setmaxleveltemp, 1.0, cLine, fLine, nLines, &
            "Temperature of maximum molecular level to be considered:","(a,f8.3,1x,a)", &
            1.0e30, ok, .false.)
       call getLogical("constantabundance", constantAbundance, cLine, fLine, nLines, &
            "Use a constant abundance: ", "(a,1l,1x,a)", .false., ok, .true.)
       call getReal("molAbundance", molAbundance, 1., cLine, fLine, nLines, &
            "Molecular Abundance:","(a,e12.5,1x,a)", 1e-9, ok, .false.)

       call getLogical("usekromeabundance", useKromeAbundance, cLine, fLine, nLines, &
            "Use previously calculated KROME abundances: ", "(a,1l,1x,a)", .false., ok, .false.)

       if (useKromeAbundance) then
          call getDouble("isotopologuefrac", isotopologueFraction, 1.d0, cLine, fLine, nLines, &
               "Fractional abundance of this isotopologue: ", "(a,f8.3,1x,a)", 0.d0, ok, .false.)
       endif

       call getLogical("removehotmolecular", removeHotMolecular, cLine, fLine, nLines, &
            "Remove molecular material above 100K: ", "(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("isinlte", isinlte, cLine, fLine, nLines, &
            "Assume LTE: ", "(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("setupMolecularLteOnly", setupMolecularLteOnly, cLine, fLine, nLines, &
            "Set up LTE in molecular calculation then exit: ", "(a,1l,1x,a)", .false., ok, .false.)
       call getReal("dusttogas", dusttoGas, 1., cLine, fLine, nLines, &
            "Dust to gas ratio: ","(a,f5.3,a)",0.01,ok,.false.)
       call getLogical("plotlevels", plotlevels, cLine, fLine, nLines, &
            "Plot Molecular Levels ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("realdust", realdust, cLine, fLine, nLines, &
            "Use realistic dust model: ", "(a,1l,1x,a)", .true., ok, .false.)

       call getLogical("outputconvergence", outputconvergence, cLine, fLine, nLines, &
            "Write out convergence data : ", "(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("dotune", dotune, cLine, fLine, nLines, &
            "Write out tuning data : ", "(a,1l,1x,a)", .false., ok, .false.)


       call getLogical("doCOchemistry", doCOchemistry, cLine, fLine, nLines, &
            "Use drop profile to model CO depletion: ", "(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("useDust", useDust, cLine, fLine, nLines, &
            "Calculate continuum emission from dust:", "(a,1l,1x,a)", .false., ok, .false.)

       if (useDust.and.(.not.dustPhysics)) then
          call writeFatal("useDust specified as true but dustphysics is F")
          stop
       endif

       if(doCOchemistry) then

          call getReal("fracCOdepletion", x_D, 1., cLine, fLine, nLines, &
               "Fraction of CO depletion", "(a,f5.3,a)", 0.1, ok, .true.)

          call getReal("fracCOnormal", x_0, 1., cLine, fLine, nLines, &
               "Fraction of CO normal level", "(a,f5.3,a)", 0.1, ok, .true.)

       endif

! Do some sanity checks
       if (restart.and.(.not.readGrid)) then
          call writeFatal("Restart option can only be used when reading a grid")
       endif

       if (restart) then
          inquire(file="restart.dat", exist=ok)
          if (.not.ok) then
             call writeFatal("Failed to find the file restart.dat")
          endif
       endif

       if (setupMolecularLteOnly.and.(.not.isinlte)) &
            call writeWarning("setupMolecularLteOnly only works when isinlte is true")

  end subroutine readMolecularPhysicsParameters

#ifdef CHEMISTRY
  subroutine readChemistryParameters(cline, fline ,nLines)
    use krome_main
    use krome_user
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines,i
    character(len=20) :: keyword
    character(len=16) :: kromeNames(1:krome_nmols)
    character(len=50) :: outputString
    logical :: ok

    call getDouble("timechemistry", timeChemistry, yearstoSecs, cLine, fLine, nLines, &
            "Time for chemistry to evolve: ","(a,1p,e12.3,1x,a)", 1.d6, ok, .true.)


    kromeNames = krome_get_names()
    allocate(kromeInitialAbundances(1:krome_nmols))
    kromeInitialAbundances = 1.d-20
    do i = 1, krome_nmols
       write(keyword,'(a)') "krome_"//trim(kromeNames(i))
       write(outputString,'(a,a,a)') "Initial krome abundance for ",trim(kromenames(i))," :"
       call getDouble(keyword, kromeInitialAbundances(i), 1.d0, cLine, fLine, nLines, &
            trim(outputString),"(a,1p,e12.3,1x,a)", 1.d-20, ok, .false.)
    enddo
  end subroutine readChemistryParameters
#endif

#ifdef PHOTOION
  subroutine readPhotoionPhysicsParameters(cLine, fLine, nLines)
    use ion_mod, only: setAbundances
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines, thisStackLim, numMPIthreads, i
    character(len=4)  :: iChar
    character(len=20) :: keyword
    logical :: ok
    real :: h_abund, he_abund, c_abund, n_abund, o_abund, ne_abund, s_abund

    call getLogical("timedep", timeDependentRT, cLine, fLine, nLines, &
         "Time-dependent RT: ", "(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("quickthermal", quickThermal, cLine, fLine, nLines, &
         "Compute photoionization equilibrium: ","(a,1l,a)", .false., ok, .false.)

    decoupleGasDustTemperature = .false.
    if (.not. quickthermal) then
       call getLogical("decouplegasdust", decoupleGasDustTemperature, cLine, fLine, nLines, &
            "Decouple gas and dust temperature: ","(a,1l,1x,a)", .false., ok, .false.)
    endif

    call getLogical("isothermal", isoThermal, cLine, fLine, nLines, &
         "Isothermal (no photoionization loop): ","(a,1l,1x,a)", .false., ok, .false.)


    call getLogical("usemetals", usemetals, cLine, fLine, nLines, &
         "Use metals in photoionization calculation: ","(a,1l,a)", .true., ok, .false.)

    call getDouble("gasmetallicity", gasMetallicity, 1.d0, cLine, fLine, nLines, &
         "Metallicity of gas relative to Zsolar: ","(a,f6.1,1x,a)", 1.d0, ok, .false.)

    call getLogical("xraymetals", usexraymetals, cLine, fLine, nLines, &
         "Use x-ray metals in photoionization calculation: ","(a,1l,a)", .false., ok, .false.)

    call getLogical("vardustsub", variableDustSublimation, cLine, fLine, nLines, &
         "Variable dust sublimation temperature: ", "(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("nodiffuse", noDiffuseField, cLine, fLine, nLines, &
         "Ignore diffuse ionizing radiation field: ","(a,1l,a)", .false., ok, .false.)

    call getLogical("noionization", noIonization, cLine, fLine, nLines, &
         "Don't ionize (set ionization potentials to high number): ","(a,1l,a)", .false., ok, .false.)

    call getLogical("radpresstest", radPressureTest, cLine, fLine, nLines, &
         "Radiation-pressure test (on-the-spot absorption): ","(a,1l,a)", .false., ok, .false.)

    call getLogical("UV_vector", UV_vector, cLine, fLine, nLines, &
         "Dump UV vectors to file for 3D-PDR: ","(a,1l,a)", .false., ok, .false.)

    call getLogical("dummyuv", dummyuv, cLine, fLine, nLines, &
         "Assign a dummy uv vector field: ","(a,1l,a)", .false., ok, .false.)

    call getLogical("drainefromuv", drainefromuv, cLine, fLine, nLines, &
         "Create draine field from star: ","(a,1l,a)", .false., ok, .false.)

    if(dummyuv) then
       call getDouble("dummyval", dummyval, 1.d0, cLine, fLine, nLines, &
            "Assumed UV field (draines):  ","(a,e12.3,1x,a)", 100.d0, ok, .false.)
    endif

    if(UV_vector) then
       call getDouble("UV_low", UV_low, 1.d0, cLine, fLine, nLines, &
            "Lower bound of UV energy field (eV):  ","(a,e12.3,1x,a)", 5.d0, ok, .false.)

       call getDouble("UV_high", UV_high, 1.d0, cLine, fLine, nLines, &
            "Lower bound of UV energy field (eV):  ","(a,e12.3,1x,a)", 50.d0, ok, .false.)
    endif

    if(nodiffusefield .or. UV_vector) then
       call getLogical("monochromatic", monochromatic, cLine, fLine, nLines, &
            "Use a monochromatic source:", "(a,1l,1x,a)", .false., ok, .false.)
       if(monochromatic) then
          call getDouble("inputEV", inputEV, 1.d0, cLine, fLine, nLines, &
               "Energy of monochromatic photons (eV):  ","(a,e12.3,1x,a)", 5.d0, ok, .false.)
       endif
    end if


    call getLogical("biasToLyman", biasToLyman, cLine, fLine, nLines, &
         "Variance reduction, higher sampling of Lyman photons: ","(a,1l,1x,a)", .false., ok, .false.)

    if (biasToLyman) then
       call getDouble("biasMagnitude", biasMagnitude, 1.d0, cLine, fLine, nLines, &
               "Variance reduction, extent of bias: ", "(a,1p,e9.3,1x,a)", 100.d0, ok, .false.)
    endif

    call getLogical("binPhotons", binPhotons, cLine, fLine, nLines, &
         "Bin and dump photons as a function of wavelength: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("hOnly", hOnly, cLine, fline, nLines, &
         "Hydrogen-only calculation: ", "(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("simplemu", simpleMu, cLine, fline, nLines, &
         "Use simple mean particle mass cals: ", "(a,1l,1x,a)", .false., ok, .false.)
    call getLogical("caseB", caseB, cLine, fline, nLines, &
         "Use case B recombination coefficient: ", "(a,1l,1x,a)", .false., ok, .false.)

    if(caseB .and. .not. nodiffusefield) then
       call writeWarning("Using a case B recombination coefficient with the diffuse field!")
    end if

    call getLogical("startfromneutral", startFromNeutral, cLine, fline, nLines, &
         "Start photoionization loop from neurtral: ", "(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("stellarwinds", stellarWinds, cLine, fline, nLines, &
         "Include stellar wind outflow feedback: ", "(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("supernovae", supernovae, cLine, fline, nLines, &
         "Include supernova feedback: ", "(a,1l,1x,a)", .false., ok, .false.)

    if (stellarWinds .or. supernovae) then
      call getDouble("windradius", stellarWindRadius, 1.d0, cLine, fLine, nLines, &
            "Stellar wind radius (smallestCellSize): ","(a,f5.2,1x,a)", 5.d0, ok, .false.)
    endif

    call getReal("tminglobal", TMinGlobal, 1., cLine, fLine, nLines, &
         "Minimum Temperature (K): ","(a,f4.1,1x,a)", 10., ok, .false.)
    call getReal("tmaxglobal", TMaxGlobal, 1., cLine, fLine, nLines, &
         "Maximum Temperature (K): ","(a,f7.1,1x,a)", 2.e4, ok, .false.)

    if(splitovermpi) then
       call getLogical("optimizeStack", optimizeStack, cLine, fLine, nLines, &
            "Perform a bundle size optimization calculation: ","(a,1l,1x,a)", .false., ok, .false.)

       call getInteger("stacklimit", stacklimit, cLine, fLine, nLines, &
            "Maximum number of photons in a stack","(a,i8,a)", 10, ok, .false.)

       call getInteger("dstack", dstack, cLine, fLine, nLines, &
            "Optimization tweak level","(a,i8,a)", 100, ok, .false.)

       call getLogical("bufferedsend", bufferedSend, cLine, fLine, nLines, &
            "Use buffered send for photon stacks: ","(a,1l,1x,a)", .true., ok, .false.)

       call getLogical("customstacks", customStacks, cLine, fLine, nLines, &
            "Specify stack limits on a thread by thread basis: ","(a,1l,1x,a)", .false., ok, .false.)

       if(customStacks) then
          call getInteger("numMPIthreads", numMPIthreads, cLine, fLine, nLines, &
               "Number of MPI threads, for custom stack specification","(a,i8,a)", 100, ok, .false.)

          allocate(stackLimitArray(numMPIthreads))

          do i = 0, numMPIthreads-1
             ! Set up a left adjusted string containing the image number and trailing spaces
             iChar="    "
             write(iChar,'(i4)') i
             iChar = adjustl(iChar)

             write(keyword,'(a)') "stacklimarray"//iChar
             call getInteger(keyword, thisStackLim, cLine, fLine, nLines, &
                  "Stack limit array element size","(a,i8,a)", stackLimit, ok, .false.)

             stackLimitArray(i+1) = thisStackLim
          end do
       end if

       call getLogical("periodicX", periodicX, cLine, fLine, nLines, &
            "Use periodic photon boundary conditions in x direction:", "(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("periodicY", periodicY, cLine, fLine, nLines, &
            "Use periodic photon boundary conditions in y direction:", "(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("periodicZ", periodicZ, cLine, fLine, nLines, &
            "Use periodic photon boundary conditions in z direction:", "(a,1l,1x,a)", .false., ok, .false.)
       call getInteger("bufferCap", bufferCap, cLine, fLine, nLines, &
!            "Number of photon stacks allowed in the buffer ","(a,i8,a)", 40000, ok, .false.)
            "Number of photon stacks allowed in the buffer ","(a,i8,a)", &
            max(40000,nint(dble(inputnmonte)/dble(stacklimit))), ok, .false.)
    end if

!
! Abundances
!
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

    call setAbundances(h_abund, he_abund, c_abund, n_abund, o_abund, ne_abund, s_abund)

  end subroutine readPhotoionPhysicsParameters
#endif

  subroutine readMolecularLoopParameters(cLine, fLine, nLines)
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call getLogical("hydrovelconv", hydroVelocityConv, cLine, fLine, nLines, &
         "Take velocity data from hydrodynamics: ","(a,1l,1x,a)", .false., ok, .false.)

  end subroutine readMolecularLoopParameters

  subroutine readAtomicLoopParameters(cLine, fLine, nLines)
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call getLogical("lte", lte, cLine, fLine, nLines, &
         "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)

    call getReal("vturb", vturb, real(kmstoc), cLine, fLine, nLines, &
         "Turbulent velocity (km/s):","(a,f6.1,1x,a)", 50., ok, .true.)


  end subroutine readAtomicLoopParameters

  subroutine readRadiativeEquilibriumParameters(cLine, fLine, nLines)
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok
    call getInteger("nlucy", nLucy, cLine, fLine, nLines,"Number of photons per lucy iteration: ","(a,i12,a)",0,ok,.false.)

    call getInteger("iterLucy", iterLucy, cline, fLine, nlines, "Minimum number of Lucy iterations: ", "(a,i3,a)",3,ok,.false.)

    call getInteger("maxiter", maxIterLucy, cline, fLine, nlines, "Maximum number of Lucy iterations: ", "(a,i3,a)",20,ok,.false.)

    call getInteger("maxgaussiter", maxGaussIter, cline, fLine, nlines, "Maximum number of G-S iterations: ", "(a,i3,a)", &
         10000,ok,.false.)

    call getLogical("forceLucyConv", forceLucyConv, cLine, fLine, nLines, &
         "Force convergence of Lucy algorithm: ","(a,1l,1x,a)", .false., ok, .false.)


    call getLogical("restart", restartLucy, cLine, fLine, nLines, &
         "Restart Lucy algorithm from save point: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("solvediffusion", solveDiffusionZone, cLine, fLine, nLines, &
         "Solve diffusion region using Gauss-Seidel: ","(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("quicksublimate", quickSublimate, cLine, fLine, nLines, &
         "Quick and dirty dust sublimation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getReal("lucy_undersampled", lucy_undersampled, 1., cLine, fLine, nLines, &
         "Minimum percentage of undersampled cell in lucy iteration: ", &
         "(a,f6.1,a)",0.0,ok,.false.)

    call getLogical("convergeonundersampled", convergeOnUndersampled, cLine, fLine, nLines, &
         "Prevent convergence when too many cells undersampled: ", "(a,1l,1x,a)", .true., ok, .false.)


    call getLogical("writelucytmpfile", writeLucyTmpFile, cLine, fLine, nLines, &
         "Write temporary lucy grid file: ", "(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("vardustsub", variableDustSublimation, cLine, fLine, nLines, &
         "Variable dust sublimation temperature: ", "(a,1l,1x,a)", .false., ok, .false.)


    call getReal("tthresh", tthresh, 1., cLine, fLine,  nLines, &
         "Dust sublimation temperature (K): ","(a,f8.2,a)", 0., ok, .false.)


    call getInteger("mincrossings", minCrossings, cLine, fLine, nLines, &
         "Minimum crossings required for cell to be sampled: ","(a,i12,a)",200,ok,.false.)

    call getReal("tminglobal", TMinGlobal, 1., cLine, fLine, nLines, &
         "Minimum Temperature (K): ","(a,f4.1,1x,a)", 10., ok, .false.)

    call getReal("taudiff", tauDiff, 1., cLine, fLine, nLines, &
         "Mininum optical depth of cell to be in diffusion approx : ","(a,f7.1,a)",100., ok, .false.)

    call getReal("tauforce", tauForce, 1., cLine, fLine, nLines, &
         "Forced optical depth of cell to be in diffusion approx : ","(a,f7.1,a)",10., ok, .false.)


    call getLogical("resetdiffusion", resetDiffusion, cLine, fLine, nLines, &
         "Reset diffusion zones to false if thin: ","(a,1l,a)",.true., ok, .false.)

    call getLogical("dosmoothgrid", doSmoothGrid, cLine, fLine, nLines, &
         "Smooth AMR grid: ","(a,1l,1x,a)", .false., ok, .false.)

    call getReal("smoothfactor", smoothFactor, 1.0, cLine, fLine, nLines, &
         "Inter-cell maximum ratio before smooth: ","(a,f6.1,1x,a)", 3., ok, .false.)

    call getLogical("smoothgridtau", doSmoothGridtau, cLine, fLine, nLines, &
         "Smooth AMR grid using tau: ","(a,1l,1x,a)", .false., ok, .false.)

      call getLogical("multilucyfiles", multilucyfiles, cLine, fLine, nLines, &
            "Multiple lucy vtk files: ","(a,1l,1x,a)", .false., ok, .false.)

! Also used for calculating tau in VTK output
    call getReal("lambdasmooth", lambdasmooth, 1.0, cLine, fLine, nLines, &
         "Lambda for tau smoothing: ","(a,1PE10.3,1x,a)", 5500.0, ok, .false.)

    if (dosmoothgridtau) then
       call getReal("taumax", tauSmoothMax, 1.0, cLine, fLine, nLines, &
            "Maximum tau for smoothing: ","(a,f10.1,1x,a)", 1., ok, .false.)
       call getReal("taumin", tauSmoothMin, 1.0, cLine, fLine, nLines, &
            "Minimum tau for smoothing: ","(a,f10.1,1x,a)", 0.1, ok, .false.)
    endif

    call getReal("edenstol", eDensTol, 1., cLine, fLine, nLines, &
         "Fractional change in energy density for convergence: ","(a,f7.1,a)",0.001, ok, .false.) ! used for gauss-seidel sweep also


    call getLogical("storescattered",storeScattered, cLine, fLine, nLines, &
         "Store scattered light: ", "(a,1l,1x,a)", .false., ok, .false.)

    call getReal("scatteredlightwavelength", scatteredLightWavelength, 1., cLine, fLine, nLines, &
         "Wavelength of scattered light","(a,f7.1,a)",1.e4, ok, .false.)


  end subroutine readRadiativeEquilibriumParameters

  subroutine readPhotoionEquilibriumParameters(cLine, fLine, nLines)
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call getInteger("nlucy", nLucy, cLine, fLine, nLines,"Number of photons per lucy iteration: ","(a,i12,a)",0,ok,.false.)

    call getbigInteger("nmonte", inputnMonte, cLine, fLine, nLines, &
         "Number of photons in image","(a,i12,a)", 0_bigInt, ok, .false.)

    call getLogical("massiveStars", massiveStars, cLine, fLine, nLines, &
         "Only include stars over 20Msol in photoion calc: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("forceminrho", forceMinRHo, cLine, fLine, nLines, &
         "Enforce a minimum density of 1mHydrogen cm^-3: ","(a,1l,1x,a)", .false., ok, .false.)

    call getInteger("mincrossings", minCrossings, cLine, fLine, nLines, &
         "Minimum crossings required for cell to be sampled: ","(a,i12,a)",200,ok,.false.)

    call getReal("taudiff", tauDiff, 1., cLine, fLine, nLines, &
         "Mininum optical depth of cell to be in diffusion approx : ","(a,f7.1,a)",100., ok, .false.)

    call getReal("edenstol", eDensTol, 1., cLine, fLine, nLines, &
         "Fractional change in energy density for convergence: ","(a,f7.1,a)",0.001, ok, .false.) ! used for gauss-seidel sweep also

    call getLogical("dosmoothgrid", doSmoothGrid, cLine, fLine, nLines, &
         "Smooth AMR grid: ","(a,1l,1x,a)", .false., ok, .false.)

    call getReal("smoothfactor", smoothFactor, 1.0, cLine, fLine, nLines, &
         "Inter-cell maximum ratio before smooth: ","(a,f6.1,1x,a)", 3., ok, .false.)

  end subroutine readPhotoionEquilibriumParameters

  subroutine readHydrodynamicsParameters(cLine, fLine, nLines)


    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    real(double) :: dx, cs
    logical :: ok

    call getReal("cfl", cflNumber, 1., cLine, fLine, nLines, &
         "Courant number:","(a,f4.1,1x,a)", 0.3, ok, .false.)

    call getLogical("gascourant", forcegascourant, cLine, fLine, nLines, &
         "Only use the gas courant condition: ","(a,1l,a)", .false., ok, .false.)

    call getLogical("advecthydro", advectHydro, cLine, fLine, nLines, &
         "Perform hydrodynamic advection: ","(a,1l,a)", .true., ok, .false.)


    call getLogical("modelwashydro", modelwashydro, cLine, fLine, nLines, &
         "Grid is based on hydro or radiation hydro calculation: ","(a,1l,a)", .false., ok, .false.)

    call getDouble("mu", mu, 1.d0, cLine, fLine, nLines, &
         "Mean molecular weight:  ","(a,f5.2,1x,a)", 2.d33, ok, .false.)

    call getDouble("etaviscosity", etaViscosity, 1.d0, cLine, fLine, nLines, &
         "Viscosity eta parameter:  ","(a,e12.3,1x,a)", 3.d0, ok, .false.)

    call getInteger("CD_version", CD_version, cLine, fLine, nLines, &
         "Contact discontinuity test version: ","(a,1x,i4,a)", 1, ok, .false.)


    call getDouble("rhofloor", rhoFloor, 1.d0, cLine, fLine, nLines, &
         "Minimum density in advection:  ","(a,e12.3,1x,a)", 1.d-30, ok, .false.)

    call getLogical("tensorviscosity", useTensorViscosity, cLine, fLine, nLines, &
         "Use tensor form for artificial viscosity: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("cylindricalhydro", cylindricalHydro, cLine, fLine, nLines, &
         "Hydrodynamics in cylindrical coordinates: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("vtkincludeghosts", vtkIncludeGhosts, cLine, fLine, nLines, &
         "Include ghost cells in vtk files: ","(a,1l,1x,a)", .true., ok, .false.)

    call getInteger("vtuToGrid", vtuToGrid, cLine, fLine, nLines, &
         "specify how many vtu files to dump for each grid file: ","(a,1x,i4,a)", 1, ok, .false.)

    call getLogical("imposefixing", imposeFixing, cLine, fLine, nLines, &
         "Impose fixed regions on the hydro models: ","(a,1l,1x,a)", .false., ok, .false.)



    if (cylindricalHydro) then
       !       if(.not. useionparam) then
       amrGridCentreX = amrgridsize/2.
       dx = dble(amrgridSize)/dble(2**5-4)
       amrGridSize = real(dble(amrGridsize) + 4.0d0*dx)
       if (.not.checkPresent("amrgridcentrez", cline, nlines)) then
          amrGridCentrez =  0.49999999*amrGridSize / real(2**maxDepthAMR)! can mess up locate functions if centre is exactly half a cell from a boundary
          amrGridCentrez =  amrGridCentreZ * 1.000001d0
       endif
       vtkIncludeGhosts = .false.
       !       endif
       call getDouble("alpha", alphaViscosity, 1.d0, cLine, fLine, nLines, &
            "Alpha Viscosity: ","(a,f7.2,1x,a)", 0.3d0, ok, .false.)
    endif
!

    if (spherical) then
       amrGridCentreX = amrgridsize/2.
       dx = dble(amrgridSize)/dble(2**4-4)
       amrGridSize = real(dble(amrGridsize) + 4.0d0*dx)
       vtkIncludeGhosts = .false.
    endif

    call getDouble("griddistancescale", gridDistanceScale, 1.d0, cLine, fLine, nLines, &
         "Distance grid scale: ","(a,e12.3,1x,a)", 1.d10, ok, .true.)

    call getUnitDouble("tstart", tStart, "time", cLine, fLine, nLines, &
         "Start time for hydrodynamical calculatione: ","(a,e12.3,1x,a)", 0.d0, ok, .false.)

    call getUnitDouble("tend", tEnd, "time", cLine, fLine, nLines, &
         "End time for calculation: ","(a,e12.3,1x,a)", 1.d10, ok, .true.)

    call getUnitDouble("tdump", tDump, "time", cLine, fLine, nLines, &
         "Time between dump files: ","(a,e12.3,1x,a)", 0.d0, ok, .true.)

    call getLogical("rhieChow", rhieChow, cLine, fLine, nLines, &
         "Use Rhie-Chow interpolation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("radpressure", radiationPressure, cLine, fLine, nLines, &
         "Include radiation pressure terms: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("gasgravity", doGasGravity, cLine, fLine, nLines, &
         "Include gas gravity: ","(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("fluxinterp", fluxinterp, cLine, fLine, nLines, &
         "Interpolate flux at fine coarse boundaries: ","(a,1l,1x,a)", .false., ok, .false.)

    call getString("limitertype", limiterType, cLine, fLine, nLines, &
         "Flux limiter to use: ","(a,a,1x,a)","superbee", ok, .false.)

    call getLogical("includepressure", includePressureTerms, cLine, fLine, nLines, &
         "Include pressure source terms: ","(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("readturb", readTurb, cLine, fLine, nLines, &
         "Read in a turbulent velocity field from file: ","(a,1l,1x,a)", .false., ok, .false.)
    if (readTurb) then
       call getReal("virialalpha", virialAlpha, 1.0, cLine, fLine, nLines, &
            "Virial parameter alpha: ","(a,f7.2,a)", 2.0, ok, .true.)
     endif

    call getLogical("useviscosity", useViscosity, cLine, fLine, nLines, &
         "Use viscosity?: ","(a,1l,1x,a)", .true., ok, .false.)


    call getLogical("fixedrhobound", fixedRhoBound, cLine, fLine, nLines, &
         "Use fixed density boundary conditions?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getDouble("rhoconst", rho_const, 1.d0, cLine, fLine, nLines, &
         "Density for fixed density boundary conditions: ","(a,e12.3,1x,a)", 10.d0*mHydrogen, ok, .false.)


    xplusboundstring = "null"
    xminusboundstring = "null"
    yplusboundstring = "null"
    yminusboundstring = "null"
    zplusboundstring = "null"
    zminusboundstring = "null"

    call getString("xplusboundstring", xplusboundString, cLine, fLine, nLines, &
         "positive x boundary condition:  ","(a,a,a)","null",ok, .false.)
    if(xplusboundstring(1:4) /= "null")  xplusbound = getBoundaryCode(xplusboundString)

    call getString("xminusboundstring", xminusboundString, cLine, fLine, nLines, &
         "negative x boundary condition:  ","(a,a,a)","null",ok, .false.)
    if(xminusboundString(1:4) /= "null") xminusbound = getBoundaryCode(xminusboundString)

    call getString("yplusboundstring", yplusboundString, cLine, fLine, nLines, &
         "positive y boundary condition:  ","(a,a,a)","null",ok, .false.)
    if(yplusboundString(1:4) /= "null") yplusbound = getBoundaryCode(yplusboundString)

    call getString("yminusboundstring", yminusboundString, cLine, fLine, nLines, &
         "negative y boundary condition:  ","(a,a,a)","null",ok, .false.)
    if(yminusboundString(1:4) /= "null") yminusbound = getBoundaryCode(yminusboundString)

    call getString("zplusboundstring", zplusboundString, cLine, fLine, nLines, &
         "positive z boundary condition:  ","(a,a,a)","null",ok, .false.)
    if(zplusboundString(1:4) /= "null") zplusbound = getBoundaryCode(zplusboundString)

    call getString("zminusboundstring", zminusboundString, cLine, fLine, nLines, &
         "negative z boundary condition:  ","(a,a,a)","null",ok, .false.)
    if(zminusboundString(1:4) /= "null") zminusbound = getBoundaryCode(zminusboundString)

    call getLogical("xslope", xslope, cLine, fLine, nLines, &
         "Inflow gradient varies along x-direction: ","(a,1l,a)", .false., ok, .false.)

    call getLogical("yslope", yslope, cLine, fLine, nLines, &
         "Inflow gradient varies along y-direction: ","(a,1l,a)", .false., ok, .false.)

    call getLogical("zslope", zslope, cLine, fLine, nLines, &
         "Inflow gradient varies along z-direction: ","(a,1l,a)", .false., ok, .false.)

    call getDouble("inflowrho", inflowRho, 1.d0, cLine, fLine, nLines, &
         "Density for inflow boundary conditions: ","(a,e12.3,1x,a)", 1.d-23, ok, .false.)

    call getDouble("inflowtemp", inflowtemp, 1.d0, cLine, fLine, nLines, &
         "Temperature for inflow boundary conditions: ","(a,e12.3,1x,a)", 10.d0, ok, .false.)

    inflowPressure = (inflowrho / (2.33d0 * mHydrogen)) * kerg * inflowTemp
    cs = sqrt(inflowPressure/inflowRho)

    call getDouble("inflowspeed", inflowSpeed, 1.d0, cLine, fLine, nLines, &
         "Mach number of speed for inflow boundary conditions: ","(a,e12.3,1x,a)", 3.d0, ok, .false.)

    inflowSpeed = inflowSpeed * cs

    inflowMomentum = inflowSpeed * inflowRho
    inflowEnergy = kerg * inflowTemp/(2.33d0*mHydrogen)
    inflowEnergy = inflowEnergy + 0.5d0*(inflowSpeed)**2
    inflowRhoE = inflowEnergy * inflowRho

  end subroutine readHydrodynamicsParameters


  subroutine readPhotometryParameters(cLine, fLine, nLines)
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call getReal("inclination", thisinclination, real(degtorad), cLine, fLine, nLines, &
         "Inclination angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)
    call getReal("distance", gridDistance, real(pcToCm), cLine, fLine, nLines, &
         "Grid distance (pc): ","(a,f4.1,1x,a)", 100., ok, .false.)
    call getDouble("limbaB", sourcelimbaB, 1.d0, cLine, fLine, nLines, &
         "Limb darkening a coefficient (B): ","(a,e9.3,a)",0.d0, ok, .false.)
    call getDouble("limbbB", sourcelimbbB, 1.d0, cLine, fLine, nLines, &
         "Limb darkening b coefficient (B): ","(a,e9.3,a)",0.d0, ok, .false.)
    call getDouble("limbaV", sourcelimbaV, 1.d0, cLine, fLine, nLines, &
         "Limb darkening a coefficient (V): ","(a,e9.3,a)",0.d0, ok, .false.)
    call getDouble("limbbV", sourcelimbbV, 1.d0, cLine, fLine, nLines, &
         "Limb darkening b coefficient (V): ","(a,e9.3,a)",0.d0, ok, .false.)

  end subroutine readPhotometryParameters

#ifdef FITSCUBE
  subroutine readDataCubeParameters(cLine, fLine, nLines)
    use datacube_mod, only: cubePositionAngle, datacubeFilename, setCubeParams, npixels, RA, DEC
    use angularImage_utils
    use h21cm_mod, only: h21cm

    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok
    TYPE(VECTOR) :: gridCentre
    character(len=100) :: message
    real :: cubeAspectRatio ! Aspect ratio of spatial axes
    real :: WV_background ! background level for moment maps

    call getString("datacubefile", datacubeFilename, cLine, fLine, nLines, &
         "Output datacube  filename: ","(a,a,1x,a)","none", ok, .true.)
    call getReal("cubeaspectratio", cubeAspectRatio, 1.0, cLine, fLine, nLines, &
         "Data cube spatial aspect ratio: ","(a,f4.1,1x,a)", 1.0, ok, .false.)
    call getReal("distance", gridDistance, 1., cLine, fLine, nLines, &
         "Grid distance (pc): ","(a,f9.1,1x,a)", 100., ok, .false.)
    call getInteger("npixels", npixels, cLine, fLine, nLines, &
         "Number of pixels per row: ","(a,1x,i4,a)", 50, ok, .true.)
    call getInteger("ncubes", ncubes, cLine, fLine, nLines, &
         "Number of data cubes: ","(a,i4,a)", 1, ok, .false.)
    call getInteger("nv", nv, cLine, fLine, nLines, &
         "Number of velocity bins: ","(a,i4,a)", 50, ok, .true.)
    call getDouble("maxVel", maxVel, 1.d0, cLine, fLine, nLines, &
         "Maximum Velocity Channel (km/s): ","(a,f8.1,1x,a)", 1.0d0, ok, .true.)
    call getDouble("minVel", minVel, 1.0_db, cLine, fLine, nLines, &
         "Minimum Velocity Channel (km/s): ","(a,f8.1,1x,a)", -1.0d0*maxVel, ok, .false.)
    call getLogical("h21cm", h21cm, cLine, fLine, nLines, &
         "Calculate data cube of 21cm emission: ","(a,1l,a)", .false., ok, .false.)
    call getReal("wvbackground", WV_background, 1.0, cLine, fLine, nLines, &
         "Background to subtract for moment maps: ","(a,f4.1,1x,a)", -1.0, ok, .false.)

    call getDouble("DEC", dec, 1.d0, cLine, fLine, nLines, &
         "Reference declination: ","(a,f8.1,1x,a)", 0.0d0, ok, .false.)

    call getDouble("RA", ra, 1.d0, cLine, fLine, nLines, &
         "Reference right ascension: ","(a,f8.1,1x,a)", 0.0d0, ok, .false.)


    if (atomicPhysics) then
       call getReal("inclination", thisinclination, real(degtorad), cLine, fLine, nLines, &
            "Inclination angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)
       call getReal("positionangle", cubePositionAngle, real(degtorad), cLine, fLine, nLines, &
            "Position angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)
       call getInteger("ninc", nDataCubeInclinations,cLine, fLine, nLines, &
            "Number of inclinations: ","(a,i2,1x,a)", 1, ok, .true.)
       allocate(datacubeInclinations(nDataCubeInclinations))
       call getRealArray("inclinations", datacubeInclinations, real(degtorad), cLine, fLine, nLines, &
            "Inclinations (deg): ",90., ok, .true.)
       call  getInteger("nlamline", nlamLine,cLine, fLine, nLines, &
            "Number of line emission wavelength: ","(a,i2,1x,a)", 1, ok, .true.)
       allocate(lamLineArray(1:nLamLine))
       call getRealArray("lamline", lamLineArray, 1.,cLine, fLine, nLines, &
            "Line emission wavelengths: ", 850., ok, .true.)

       call getReal("vturb", vturb, real(kmstoc), cLine, fLine, nLines, &
            "Turbulent velocity (km/s):","(a,f6.1,1x,a)", 50., ok, .true.)

       call  getInteger("nr1", nr1,cLine, fLine, nLines, &
            "Number of radial rays to capture stellar photosphere: ","(a,i2,1x,a)", 200, ok, .false.)
       call  getInteger("nphi1", nphi1,cLine, fLine, nLines, &
            "Number of azimuthal rays to capture stellar photosphere: ","(a,i2,1x,a)", 100, ok, .false.)
       nr2 = 0
       nphi2 = 0
       if (ttauriMagnetosphere) then
          call  getInteger("nr2", nr2,cLine, fLine, nLines, &
               "Number of radial rays to capture stellar magnetosphere: ","(a,i2,1x,a)", 200, ok, .false.)
          call  getInteger("nphi2", nphi2,cLine, fLine, nLines, &
               "Number of azimuthal rays to capture stellar magnetosphere: ","(a,i2,1x,a)", 100, ok, .false.)
       endif
       nr3 = 0
       nphi3 = 0
       if (ttauriWind) then
          call  getInteger("nr3", nr3,cLine, fLine, nLines, &
               "Number of radial rays to capture disc wind: ","(a,i2,1x,a)", 200, ok, .false.)
          call  getInteger("nphi3", nphi3,cLine, fLine, nLines, &
               "Number of azimuthal rays to capture disc wind: ","(a,i2,1x,a)", 100, ok, .false.)
       endif
       nr4 = 0
       nphi4 = 0

       if (ttauristellarWind) then
          call  getInteger("nr4", nr4,cLine, fLine, nLines, &
               "Number of radial rays to capture stellar wind: ","(a,i2,1x,a)", 200, ok, .false.)
          call  getInteger("nphi4", nphi4,cLine, fLine, nLines, &
               "Number of azimuthal rays to capture stellar wind: ","(a,i2,1x,a)", 100, ok, .false.)
       endif

    endif

    call getLogical("internalView", internalView, cLine, fLine, nLines, &
         "View as our Galaxy:", "(a,1l,1x,a)", .false., ok, .false.)

    ! Read parameters used by galactic plane survey
    if ( internalView ) then

       call getReal("imageside", imageside, 1., cLine, fLine, nLines, &
            "Image size (degrees):","(a,1p,f7.2,1x,a)", amrGridSize, ok, .false.)

       call getLogical("splitCubes", splitCubes, cLine, fLine, nLines, &
            "Split intensity into +/- components: ","(a,1l,a)", .false., ok, .false.)

       call getLogical("refineQ2Only", refineQ2Only, cLine, fLine, nLines, &
            "Limit refinement to 2nd quadrant:", "(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("obsVelFromGrid", obsVelFromGrid, cLine, fLine, nLines, &
            "Set observer velocity from grid:", "(a,1l,1x,a)", .false., ok, .true.)

       call getDouble("intPosX", intPosX,  1.0_db, cLine, fLine, nLines, "Observer x position (x10^10cm):", &
            "(a,e10.4,1x,a)", 0.d0, ok, .false.)
       call getDouble("intPosY", intPosY,  1.0_db, cLine, fLine, nLines, "Observer y position (x10^10cm):", &
            "(a,e10.4,1x,a)", 2.2e12_db, ok, .false.)
       call getDouble("intPosZ", intPosZ, 1.0_db,  cLine, fLine, nLines, "Observer z position (x10^10cm):", &
            "(a,e10.4,1x,a)", 0.d0, ok, .false.)

       call getDouble("intDeltaVx", intDeltaVx, 1.0_db,  cLine, fLine, nLines, "Observer x velocity boost (km/s):", &
            "(a,f8.2,1x,a)", 0.d0, ok, .false.)
       call getDouble("intDeltaVy", intDeltaVy, 1.0_db,  cLine, fLine, nLines, "Observer y velocity boost (km/s):", &
            "(a,f8.2,1x,a)", 0.d0, ok, .false.)
       call getDouble("intDeltaVz", intDeltaVz, 1.0_db,  cLine, fLine, nLines, "Observer z velocity boost (km/s):", &
               "(a,f8.2,1x,a)", 0.d0, ok, .false.)

       call getReal("vturb", vturb, 1., cLine, fLine, nLines, &
            "Subsonic turbulent velocity (km/s):","(a,f4.1,1x,a)", 0.0, ok, .false.)

       ! This is only implemented in angularImage_mod
       call getInteger("dssminsample", dssminsample, cLine, fLine, nLines, &
            "Minimum number of density subsamples: ","(a,i4,a)", 5, ok, .false.)

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

       if ( readgrid ) then
       ! When restarting particles are required for map_dI_to_particles
          call getString("sphdatafilename", sphdatafilename, cLine, fLine, nLines, &
               "Input sph data file: ","(a,a,1x,a)","sph.dat.ascii", ok, .true.)

          call getString("inputFileFormat", inputFileFormat, cLine, fLine, nLines, &
               "Input file format: ","(a,a,1x,a)","binary", ok, .false.)
       end if

    else
       ! Switchable units are only in molecular mod
       call getString("datacubeunits", datacubeunits, cLine, fLine, nLines, &
            "Output datacube units: ","(a,a,1x,a)","intensity", ok, .false.)
       call getString("datacubeaxisunits", datacubeaxisunits, cLine, fLine, nLines, &
            "Data cube axis units: ","(a,a,1x,a)","torus", ok, .false.)

       select case (datacubeaxisunits)
          case ("torus")
             call getReal("imageside", imageside, 1., cLine, fLine, nLines, &
                  "Image size (x10^10):","(a,1p,e9.3,1x,a)", amrGridSize, ok, .false.)
          case ("pc")
             call getReal("imageside", imageside, real(pctocm/1.0e10), cLine, fLine, nLines, &
                  "Image size (pc):","(a,1p,e9.3,1x,a)", amrGridSize, ok, .false.)
          case ("kpc")
             call getReal("imageside", imageside, real(kpctocm/1.0e10), cLine, fLine, nLines, &
                  "Image size (kpc):","(a,1p,e9.3,1x,a)", amrGridSize, ok, .false.)
          case('au')
             call getReal("imageside", imageside, real(autocm/1.0e10), cLine, fLine, nLines, &
                  "Image size (au):","(a,1p,e9.3,1x,a)", amrGridSize, ok, .false.)
          case('cm')
             call getReal("imageside", imageside, 1.0e-10, cLine, fLine, nLines, &
                  "Image size (cm):","(a,1p,e9.3,1x,a)", amrGridSize, ok, .false.)
          case('m')
             call getReal("imageside", imageside, 1.0e-8, cLine, fLine, nLines, &
                  "Image size (m):","(a,1p,e9.3,1x,a)", amrGridSize, ok, .false.)
          end select

    end if

    if ( h21cm ) then

       call getReal("tminglobal", TMinGlobal, 1., cLine, fLine, nLines, &
            "Minimum Temperature (K): ","(a,f4.1,1x,a)", 10., ok, .false.)
       call getInteger("nSubpixels", nSubpixels, cLine, fLine, nLines, &
            "Subpixel splitting (0 denotes adaptive)","(a,i4,a)", 1, ok, .false.)
       call getLogical("densitysubsample", densitysubsample, cLine, fLine, nLines, &
            "Use density interpolation: ","(a,1l,a)", .false., ok, .false.)

       if (internalView) then
          call getDouble("centrevecx", centrevecx, 1.d0, cLine, fLine, nLines, &
               "Image Centre Coordinate (lon): ","(a,f7.2,1x,a)", 0.d0, ok, .true.)
          call getDouble("centrevecy", centrevecy, 1.d0, cLine, fLine, nLines, &
               "Image Centre Coordinate (lat): ","(a,f7.2,1x,a)", 0.d0, ok, .true.)
       else
          if (amr3d) then
             if (cylindrical) then
                call getDouble("centrevecx", centrevecx, 1.d0, cLine, fLine, nLines, &
                  "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .false.)
                call getDouble("centrevecy", centrevecy, 1.d0, cLine, fLine, nLines, &
                  "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .false.)
             else
                call getDouble("centrevecx", centrevecx, 1.d0, cLine, fLine, nLines, &
                     "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", amrGridCentreX, ok, .false.)
                call getDouble("centrevecy", centrevecy, 1.d0, cLine, fLine, nLines, &
                     "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", amrGridCentreY, ok, .false.)
             endif
          else
             call getDouble("centrevecx", centrevecx, 1.d0, cLine, fLine, nLines, &
                  "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .false.)
             call getDouble("centrevecy", centrevecy, 1.d0, cLine, fLine, nLines, &
                  "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .false.)
          endif
       endif

       call getLogical("wanttau", wanttau, cLine, fLine, nLines, &
            "Write Tau information to datacube: ","(a,1l,1x,a)", .false., ok, .false.)
       itrans = 1 ! always 1 for the 21cm line

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

       call getReal("inclination", thisinclination, real(degtorad), cLine, fLine, nLines, &
            "Inclination angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)
       rotateViewAboutX = 90.0 - thisInclination*radtodeg


       call getDouble("rotateviewaboutz", rotateViewAboutZ, 1.d0, cLine, fLine, nLines, &
            "Angle to rotate about Z (deg): ","(a,f4.1,1x,a)", 0.d0, ok, .false.)
       if (internalView) then
          call getDouble("centrevecx", centrevecx, 1.d0, cLine, fLine, nLines, &
               "Image Centre Coordinate (lon): ","(a,f7.2,1x,a)", 0.d0, ok, .true.)
          call getDouble("centrevecy", centrevecy, 1.d0, cLine, fLine, nLines, &
               "Image Centre Coordinate (lat): ","(a,f7.2,1x,a)", 0.d0, ok, .true.)
       else

          if (amr3d) then

             if (cylindrical) then
                call getDouble("centrevecx", centrevecx, 1.d0, cLine, fLine, nLines, &
                  "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .false.)
                call getDouble("centrevecy", centrevecy, 1.d0, cLine, fLine, nLines, &
                  "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .false.)
             else
                call getDouble("centrevecx", centrevecx, 1.d0, cLine, fLine, nLines, &
                     "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", amrGridCentreX, ok, .false.)
                call getDouble("centrevecy", centrevecy, 1.d0, cLine, fLine, nLines, &
                     "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", amrGridCentreY, ok, .false.)
             endif
          else
             call getDouble("centrevecx", centrevecx, 1.d0, cLine, fLine, nLines, &
               "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .false.)
             call getDouble("centrevecy", centrevecy, 1.d0, cLine, fLine, nLines, &
                  "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .false.)
          endif
          call getDouble("centrevecz", centrevecz, 1.d0, cLine, fLine, nLines, &
               "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", amrGridCentreZ, ok, .false.)
       end if
       call getLogical("wanttau", wanttau, cLine, fLine, nLines, &
            "Write Tau information to datacube: ","(a,1l,1x,a)", .false., ok, .false.)
    endif

! Set observer orientation for molecular_mod. This is only required for the far field case
! when using molecular or 21cm physics.
molecular_orientation: if ( .not.internalView .and. (molecularPhysics.or.h21cm)) then

! We should have either rotateViewAbout? parameters galaxy* parameters but not a mixture.
       if ((checkPresent("rotateViewAboutX", cLine, nlines).or. &
            checkPresent("rotateViewAboutY", cLine, nlines).or. &
            checkPresent("rotateViewAboutZ", cLine, nlines)) .and. &
           (checkPresent("galaxyInclination", cLine, nlines).or. &
            checkPresent("galaxyPositionAngle", cLine, nlines)) ) then
          call writeFatal("Found rotateViewAbout parameter and a galaxyInclination/PositionAngle")
          stop
       end if

! Read either galaxy parameters or rotate parameters. Default will set all rotations to zero.
       if ( checkPresent("galaxyInclination", cLine, nlines).or. &
            checkPresent("galaxyPositionAngle", cLine, nlines)) then
          call getDouble("galaxyInclination", galaxyInclination, 1.0_db, cLine, fLine, nLines, &
               "Galaxy Inclination:", "(a,f5.1,1x,a)", 0.d0, ok, .false.)
          call getDouble("galaxyPositionAngle", galaxyPositionAngle, 1.0_db, cLine, fLine, nLines, &
               "Galaxy position angle:", "(a,f5.1,1x,a)", 0.d0, ok, .false.)
          rotateViewAboutX    = 90.0 - galaxyInclination
          imageBasisPrerotate = 90.0 - galaxyPositionAngle
       else
          call getDouble("rotateviewaboutx", rotateViewAboutX, 1.d0, cLine, fLine, nLines, &
               "Angle to rotate about X (deg): ","(a,f4.1,1x,a)", 0.d0, ok, .false.)
          call getDouble("rotateviewabouty", rotateViewAboutY, 1.d0, cLine, fLine, nLines, &
               "Angle to rotate about Y (deg): ","(a,f4.1,1x,a)", 0.d0, ok, .false.)
          call getDouble("rotateviewaboutz", rotateViewAboutZ, 1.d0, cLine, fLine, nLines, &
               "Angle to rotate about Z (deg): ","(a,f4.1,1x,a)", 0.d0, ok, .false.)
          imageBasisPrerotate = 0.0
       endif

    end if molecular_orientation

! Set up values in datacube_mod
    call setCubeParams(npixels, cubeAspectRatio, WV_background)

  end subroutine readDataCubeParameters
#endif
  subroutine readDetectorParameters(cLine, fLine, nLines)
    use constants_mod
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call getString("detectortype", detType, cLine, fLine, nLines,&
         "Type of detector:", "(a,a,1x,a)", "au", ok, .true.)

    call getVector("detectorpos", detectorPosition, 1.d0, cLine, fLine, nLines, &
                  "Detector position: ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

    call getDouble("detectortheta", detectorTheta, degtorad, cLine, fLine, nLines, &
                  "Detector polar angle (deg): ","(a,f6.1,a)",0.d0, ok, .true.)

    call getDouble("detectorphi", detectorPhi, degtorad, cLine, fLine, nLines, &
         "Detector azimuthal angle (deg): ","(a,f6.1,a)",0.d0, ok, .true.)

    select case(detType)
       case("ccd")
          call getDouble("detectorxsize", detectorXsize, 1.d0, cLine, fLine, nLines, &
                  "Detector X size: ","(a,f6.1,a)",0.d0, ok, .true.)

          call getDouble("detectorysize", detectorYsize, 1.d0, cLine, fLine, nLines, &
               "Detector Y size: ","(a,f6.1,a)",0.d0, ok, .true.)

          call getInteger("detectornx", detectornx, cLine, fLine, nLines, &
               "Detector X pixels ","(a,i4,a)",0, ok, .true.)

          call getInteger("detectorny", detectorny, cLine, fLine, nLines, &
                  "Detector Y pixels ","(a,i4,a)",0, ok, .true.)
       case("fibre")
          call getDouble("detradius", detRadius, microntocm, cLine, fLine, nLines, &
                  "Detector fibre radius: ","(a,f6.1,a)",0.d0, ok, .true.)
          call getDouble("detectorna", detectorNA, 1.d0, cLine, fLine, nLines, &
                  "Detector fibre numerical aperture: ","(a,f6.1,a)",0.d0, ok, .true.)


    end select


    call getString("detfilename", detFilename, cLine, fLine, nLines,&
         "Filename for detector output:", "(a,a,1x,a)", "au", ok, .true.)
  end subroutine readDetectorParameters

  subroutine readColumnImageParameters(cLine, fLine, nLines)
    use constants_mod
    use image_utils_mod
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call getString("columnfile", columnImageFilename, cLine, fLine, nLines, &
         "Output column image filename: ","(a,a,1x,a)","none", ok, .true.)

    call getVector("columndir", columnImageDirection, 1.d0, cLine, fLine, nLines, &
         "Direction for column image: ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

    call getString("columnaxisunits", columnAxisUnits, cLine, fLine, nLines,&
         "Axis units for column image:", "(a,a,1x,a)", "pc", ok, .false.)
  end subroutine readColumnImageParameters

  subroutine readImageParameters(cLine, fLine, nLines)
    use constants_mod
    use image_utils_mod

    character(len=lencLine) :: cLine(:)
    character(len=80) :: message
    logical :: fLine(:)
    integer :: nLines
    logical :: ok
    integer :: i
    character(len=20) :: keyword
    character(len=10) :: axisUnits
    character(len=12) :: fluxUnits, globalFluxUnits
    character(len=80) :: outputImageType, imageFilename
    character(len=4)  :: iChar
    integer :: thisnpixels, npixels, nimage
    real :: lambdaImage, thisLambdaImage
    real :: imagesize, thisImageSize, wholeGrid
    real :: inclination, positionAngle, thisPA, thisInc
    real :: offsetX, offsetY, thisOffsetX, thisOffsetY
    real :: aspectRatio, thisAspectRatio
    type(VECTOR) :: viewvec

    call getInteger("nphase", nPhase, cLine, fLine, nLines, &
         "Number of phases: ", "(a,i3,1x,a)", 1, ok, .false.)
    nStartPhase = 1
    nEndPhase = nPhase

    ! Used for optical depth tests in phaseloop
    call getReal("lambdatau", lambdatau, 1.0, cLine, fLine, nLines, &
         "Lambda for tau test: ","(a,1PE10.3,1x,a)", 5500.0, ok, .false.)

    call getBigInteger("nphotons", nphotons, cLine, fLine, nLines, &
         "Number of photons in image: ","(a,i9,a)", 10000_bigInt, ok, .false.)

    call getBigInteger("nphotimage", nphotimage, cLine, fLine, nLines, &
         "Number of photons in image: ","(a,i9,a)", nPhotons, ok, .false.)

    call getInteger("maxscat", maxScat, cLine, fLine, nLines, &
         "Maximum number of scatterings: ","(a,i9,a)", 10000, ok, .false.)

    call getReal("distance", gridDistance, 1., cLine, fLine, nLines, &
         "Grid distance (pc): ","(a,f10.1,1x,a)", 100., ok, .false.)

    call getInteger("nimage", nimage, cLine, fLine, nLines, &
         "Number of images to calculate: ","(a,i4,a)", 1, ok, .false.)
    if (nimage > 9999) call writeWarning("Too many images. Need nimage <= 9999.")

    call setNumImages(nimage)

    call getString("imageaxisunits", axisUnits, cLine, fLine, nLines,&
         "Axis units for image:", "(a,a,1x,a)", "au", ok, .false.)

    ! Size of the whole grid in image axis units
    ! This is used as the default image size
    select case (axisunits)
    case ("au","AU")
       wholeGrid = (amrGridSize*1.0e10) / real(auTocm)
    case ("pc","PC")
       wholeGrid = (amrGridSize*1.0e10) / real(pctocm)
    case ("cm")
       wholeGrid = amrGridSize*1.0e10
    case ("arcsec")
       wholeGrid = (amrGridSize*1.0e10) / real(pctocm) ! pc
       wholeGrid = wholeGrid / real(gridDistance)            ! radians
       wholeGrid = wholeGrid * (180.0/real(pi)) * 3600.0     ! arcsec
    case default
       wholeGrid = amrGridSize
       write(message,*) "Unrecognised units in readImageParameters: ", trim(axisunits)
       call writeWarning(message)
    end select

    write(message,*) "Image size ("//trim(axisUnits)//"): "
    call getReal("imagesize", imageSize, 1.0, cLine, fLine, nLines, &
         trim(message), "(a,1pe10.2,1x,a)", wholeGrid, ok, .false.)

    call getVector("imageorigin", imageOrigin, 1.d0, cLine, fLine, nLines, &
                  "Origin for image (10^10 cm): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .false.)


    call getString("imagefluxunits", globalFluxUnits, cLine, fLine, nLines,&
         "Flux units for image:", "(a,a,1x,a)", "MJy/str", ok, .false.)
    if (trim(globalfluxUnits) == "Jy/beam") then
       call getDouble("beamarea", beamarea, arcsecstoradians**2,cLine, fLine, nLines, &
            "Beam area in square arcseconds: ", "(a,f10.2,1x,a)", 1.d0, ok, .false.)
    endif


    !Thaw for image cut dumps
    call getLogical("dumpCut", dumpCut, cLine, fLine, nLines, &
         "Dump a file of pixel values in a cut across the image: ","(a,1l,1x,a)", .false., ok, .false.)

    call getString("cutType", cutType, cLine, fLine, nLines, &
         "Direction of cut: ","(a,a,1x,a)","horizontal", ok, .false.)


    call getInteger("sliceIndex", sliceIndex, cLine, fLine, nLines, &
         "Cell index axis of image cut","(a,i8,a)", 101, ok, .false.)

    call getLogical("dustcube", calcDustCube, cLine, fLine, nLines, &
         "Calculate a dust data cube instead of an image: ","(a,1l,1x,a)", .false., ok, .false.)

    if (calcDustCube) then
       call getString("lambdafile", lambdaFilename, cLine, fLine, nLines, &
            "Wavelength file for dust cube: ","(a,a,1x,a)","none", ok, .true.)
    endif

    call getInteger("npixels", npixels, cLine, fLine, nLines, &
         "Number of pixels per side in image","(a,i8,a)", 200, ok, .false.)

    call getReal("imageaspect", aspectRatio, 1.0, cLine, fLine, nLines, &
         "Image aspect ratio: ", "(a,f4.1,1x,a)", 1.0, ok, .false.)

! Inclination and position angle for single image, also used as default value for multiple images
    call getReal("inclination", inclination, real(degtorad), cLine, fLine, nLines, &
         "Inclination angle (deg): ","(a,f6.1,1x,a)", 0., ok, .false.)
    call getReal("positionangle", positionAngle, real(degtorad), cLine, fLine, nLines, &
         "Position angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)

! Position of image centre
    call getReal("imagecentrex", offsetx, 1.0, cLine, fLine, nLines, &
         "Image centre x position (10^10 cm): ", "(a,e10.1,1x,a)", 0.0, ok, .false.)

    call getReal("imagecentrey", offsety, 1.0, cLine, fLine, nLines, &
         "Image centre y position (10^10 cm): ", "(a,e10.1,1x,a)", 0.0, ok, .false.)

    call getDouble("fwhmpixels", fwhmpixels, 1.d0, cLine, fLine, nLines, &
         "FHWM of imaging resolution in pixels: ", "(a,f7.2,1x,a)", 0.d0, ok, .false.)


    call getLogical("polimage", polarizationImages, cLine, fLine, nLines, &
         "Write polarization images: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("forcefirstscat", forceFirstScat, cLine, fLine, nLines, &
         "Force first scattering: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("noscatlight", excludeScatteredLight, cLine, fLine, nLines, &
         "Exclude dust-scattered light from image: ","(a,1l,1x,a)", .false., ok, .false.)


    call getRealwithUnits("lambdaimage", lambdaImage, "angstrom", "angstrom", cLine, fLine, nLines, &
         "Wavelength for monochromatic image: ", 1000.0, ok, .false.)

    call getDouble("imagepa", thisimagePA, degtorad, cLine, fLine, nLines, &
         "Position angle for image :","(a,f12.2,1x,a)", 0.d0, ok, .false.)

    if (nimage == 1) then

       call getString("imagefile", imageFilename, cLine, fLine, nLines, &
            "Output image  filename: ","(a,a,1x,a)","none", ok, .true.)

       call getVector("viewvec", viewVec, 1.d0, cLine, fLine, nLines, &
                  "View vector for image: ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .false.)

       
       if (photoionPhysics) then
          call getString("imagetype", outputimageType, cLine, fLine, nLines, &
               "Type of output image: ","(a,a,1x,a)","none", ok, .true.)
          ! Reset default wavelength to 20cm if this is a free-free image
          if (outputimagetype == "freefree" .and. &
               .not. checkPresent("lambdaimage", cline, nlines)) then
             lambdaImage = 20.0e8 ! 20cm in Anstroms
             call writeInfo("Setting default wavelength to 20cm for free-free image")
          endif

       elseif(dustPhysics) then
          outputimageType(1:10) = "stokes    "
             call writeInfo("Type of output image: Stokes image (default for dust physics)")
       elseif (atomicPhysics) then
          outputimageType(1:10) = "stokes    "
          call writeInfo("Type of output image: Stokes image (default for atomic physics)")
       endif

       if (modulus(viewVec) < 0.5d0) then
          call setImageParams(1, LambdaImage, outputimageType,imageFilename, npixels, axisUnits, globalfluxUnits, &
               ImageSize, aspectRatio, inclination,  positionangle, Offsetx, Offsety, gridDistance)
       else
          call setImageParams(1, LambdaImage, outputimageType,imageFilename, npixels, axisUnits, globalfluxUnits, &
               ImageSize, aspectRatio, inclination,  positionangle, Offsetx, Offsety, gridDistance, viewvec=viewvec)
       endif
    else
       do i = 1, nImage

          ! Set up a left adjusted string containing the image number and trailing spaces
          iChar="    "
          write(iChar,'(i4)') i
          iChar = adjustl(iChar)

          call writeInfo(" ")
          write(message,'(a)') "Details for image: "//iChar
          call writeInfo(message)
          call writeInfo(" ")

          write(keyword,'(a)') "imagefile"//iChar
          call getString(keyword, imageFilename, cLine, fLine, nLines, &
               "Output image  filename: ","(a,a,1x,a)","none", ok, .true.)
          write(keyword,'(a)') "lambdaimage"//iChar
          call getRealwithUnits(keyword, thislambdaImage, "angstrom", "angstrom", cLine, fLine, nLines, &
               "Wavelength of monochromatic image: ", 1000.0, ok, .true.)

          if (photoionPhysics) then
             write(keyword,'(a)') "imagetype"//iChar
             call getString(keyword, outputimageType, cLine, fLine, nLines, &
                  "Type of output image: ","(a,a,1x,a)","none", ok, .true.)
          ! Reset default wavelength to 20cm if this is a free-free image and
          ! no wavelength has been specified in the parameters file
             write(keyword,'(a)') "lambdaimage"//iChar
             if ( outputimagetype == "freefree" .and. &
                  .not. checkPresent(keyword, cline, nlines) .and. &
                  .not. checkPresent("lambdaimage", cline, nlines)) then
                thisLambdaImage = 20.0e8 ! 20cm in Anstroms
                call writeInfo("Setting default wavelength to 20cm for free-free image")
             endif

          elseif(dustPhysics) then
             outputimageType(1:10) = "stokes    "
             call writeInfo("Type of output image: Stokes image (default for dust physics)")
          elseif (atomicPhysics) then
             outputimageType(1:10) = "stokes    "
             call writeInfo("Type of output image: Stokes image (default for atomic physics)")
          endif

          write(keyword,'(a)') "npixels"//iChar
          call getInteger(keyword, thisnpixels, cLine, fLine, nLines, &
               "Number of pixels per side in image","(a,i8,a)", npixels, ok, .false.)

          ! Size of this image
          write(keyword,'(a)') "imagesize"//iChar
          write(message,*) "Image size ("//trim(axisUnits)//"): "
          call getReal(keyword, thisImageSize, 1.0, cLine, fLine, nLines, &
            trim(message), "(a,1pe10.2,1x,a)", imageSize, ok, .false.)

          ! Aspect ratio
          write(keyword,'(a)') "imageaspect"//iChar
          call getReal(keyword, thisAspectRatio, 1.0, cLine, fLine, nLines, &
               "Image aspect ratio: ", "(a,f4.1,1x,a)", aspectRatio, ok, .false.)

          write(keyword,'(a)') "imagefluxunits"//iChar
          call getString(keyword,fluxUnits, cLine, fLine, nLines,&
               "Flux units for image:", "(a,a,1x,a)", globalfluxUnits, ok, .false.)


          ! Inclination and position angle
          write(keyword,'(a)') "inclination"//iChar
          call getReal(keyword, thisInc, real(degtorad), cLine, fLine, nLines, &
               "Inclination of image: ","(a,f6.1,1x,a)",inclination*real(radtodeg), ok, .false.)
          write(keyword,'(a)') "positionangle"//iChar
          call getReal(keyword, thisPA, real(degtorad), cLine, fLine, nLines, &
               "Position angle (deg): ","(a,f6.1,1x,a)", positionAngle*real(radtodeg), ok, .false.)

          write(keyword,'(a)') "viewvec"//iChar
          call getVector(keyword, viewVec, 1.d0, cLine, fLine, nLines, &
                  "View vector for image: ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .false.)


          ! Position of image centre
          write(keyword,'(a)') "imagecentrex"//iChar
          call getReal(keyword, thisOffsetx, 1.0, cLine, fLine, nLines, &
               "Image centre x position (10^10 cm): ", "(a,e10.1,1x,a)", offsetx, ok, .false.)
          write(keyword,'(a)') "imagecentrey"//iChar
          call getReal(keyword, thisOffsety, 1.0, cLine, fLine, nLines, &
               "Image centre y position (10^10 cm): ", "(a,e10.1,1x,a)", offsety, ok, .false.)

          if (modulus(viewVec) < 0.5d0) then
             call setImageParams(i, thisLambdaImage, outputimageType,imageFilename, thisnpixels, axisUnits, fluxUnits, &
                  thisImageSize, thisaspectRatio, thisInc, thisPA, thisOffsetx, thisOffsety, gridDistance)
          else
             call setImageParams(i, thisLambdaImage, outputimageType,imageFilename, thisnpixels, axisUnits, fluxUnits, &
                  thisImageSize, thisaspectRatio, thisInc, thisPA, thisOffsetx, thisOffsety, gridDistance,viewvec=viewvec)
          endif
       enddo

       call writeInfo(" ")

    end if

  end subroutine readImageParameters

  subroutine readMovieParameters(cLine, fLine, nLines)
    use image_utils_mod

    character(len=lencLine) :: cLine(:)
    character(len=80) :: message
    logical :: fLine(:)
    integer :: nLines
    logical :: ok
    integer :: i
    character(len=20) :: keyword
    character(len=10) :: axisUnits
    character(len=12) :: fluxUnits
    character(len=80) :: outputImageType, imageFilename

    integer :: thisnpixels, npixels, nimage
    real :: lambdaImage, thisLambdaImage
    real :: imagesize, thisImageSize, wholeGrid
    real :: inclination, positionAngle, thisPA, thisInc
    real :: offsetX, offsetY, thisOffsetX, thisOffsetY
    real :: aspectRatio, thisAspectRatio, thisTime
    type(VECTOR) :: viewVec
    real(double) :: ang

    call getBigInteger("nphotons", nphotons, cLine, fLine, nLines, &
         "Number of photons in image: ","(a,i9,a)", 10000_bigInt, ok, .false.)

    call getBigInteger("nphotimage", nphotimage, cLine, fLine, nLines, &
         "Number of photons in image: ","(a,i9,a)", nPhotons, ok, .false.)

    call getInteger("maxscat", maxScat, cLine, fLine, nLines, &
         "Maximum number of scatterings: ","(a,i9,a)", 10000, ok, .false.)

    call getReal("distance", gridDistance, 1., cLine, fLine, nLines, &
         "Grid distance (pc): ","(a,f6.1,1x,a)", 100., ok, .false.)

    call getVector("imageorigin", imageOrigin, 1.d0, cLine, fLine, nLines, &
                  "Origin for image (10^10 cm): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .false.)

    call getInteger("nimage", nimage, cLine, fLine, nLines, &
         "Number of images to calculate: ","(a,i4,a)", 1, ok, .false.)
    if (nimage > 9999) call writeWarning("Too many images. Need nimage <= 9999.")

    call setNumImages(nimage)

    call getString("imageaxisunits", axisUnits, cLine, fLine, nLines,&
         "Axis units for image:", "(a,a,1x,a)", "au", ok, .false.)

    ! Size of the whole grid in image axis units
    ! This is used as the default image size
    select case (axisunits)
    case ("au","AU")
       wholeGrid = (amrGridSize*1.0e10) / real(auTocm)
    case ("pc","PC")
       wholeGrid = (amrGridSize*1.0e10) / real(pctocm)
    case ("cm")
       wholeGrid = amrGridSize*1.0e10
    case ("arcsec")
       wholeGrid = (amrGridSize*1.0e10) / real(pctocm) ! pc
       wholeGrid = wholeGrid / real(gridDistance)            ! radians
       wholeGrid = wholeGrid * (180.0/real(pi)) * 3600.0     ! arcsec
    case default
       wholeGrid = amrGridSize
       write(message,*) "Unrecognised units in readImageParameters: ", trim(axisunits)
       call writeWarning(message)
    end select

    write(message,*) "Image size ("//trim(axisUnits)//"): "
    call getReal("imagesize", imageSize, 1.0, cLine, fLine, nLines, &
         trim(message), "(a,1pe10.2,1x,a)", wholeGrid, ok, .false.)

    call getString("imagefluxunits", fluxUnits, cLine, fLine, nLines,&
         "Flux units for image:", "(a,a,1x,a)", "MJy/str", ok, .false.)

    call getInteger("npixels", npixels, cLine, fLine, nLines, &
         "Number of pixels per side in image","(a,i8,a)", 200, ok, .false.)

    call getReal("imageaspect", aspectRatio, 1.0, cLine, fLine, nLines, &
         "Image aspect ratio: ", "(a,f4.1,1x,a)", 1.0, ok, .false.)

! Inclination and position angle for single image, also used as default value for multiple images
    call getReal("inclination", inclination, real(degtorad), cLine, fLine, nLines, &
         "Inclination angle (deg): ","(a,f6.1,1x,a)", 0., ok, .false.)

    call getReal("positionangle", positionAngle, real(degtorad), cLine, fLine, nLines, &
         "Position angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)

! Position of image centre
    call getReal("imagecentrex", offsetx, 1.0, cLine, fLine, nLines, &
         "Image centre x position (10^10 cm): ", "(a,e10.1,1x,a)", 0.0, ok, .false.)
    call getReal("imagecentrey", offsety, 1.0, cLine, fLine, nLines, &
         "Image centre y position (10^10 cm): ", "(a,e10.1,1x,a)", 0.0, ok, .false.)

    call getLogical("polimage", polarizationImages, cLine, fLine, nLines, &
         "Write polarization images: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("forcefirstscat", forceFirstScat, cLine, fLine, nLines, &
         "Force first scattering: ","(a,1l,1x,a)", .false., ok, .false.)

    call getRealwithUnits("lambdaimage", lambdaImage, "angstrom", "angstrom", cLine, fLine, nLines, &
         "Wavelength of monochromatic image: ", 1000.0, ok, .true.)

    call writeInfo(" ")
    write(message,'(a)') "Details for movie: "
    call writeInfo(message)
    call writeInfo(" ")

    call getString("moviefile", imageFilename, cLine, fLine, nLines, &
         "Output movie filename: ","(a,a,1x,a)","none", ok, .true.)


    call getRealwithUnits("movielambdae", thislambdaImage, "angstrom", "angstrom", cLine, fLine, nLines, &
         "Wavelength of monochromatic image: ", 1000.0, ok, .true.)

    outputimageType = "stokes"

    if (photoionPhysics) then
       call getString("movietype", outputimageType, cLine, fLine, nLines, &
            "Type of output image: ","(a,a,1x,a)","none", ok, .true.)
          ! Reset default wavelength to 20cm if this is a free-free image and
          ! no wavelength has been specified in the parameters file
    endif
    call getInteger("npixels", thisnpixels, cLine, fLine, nLines, &
         "Number of pixels per side in image","(a,i8,a)", npixels, ok, .false.)

    ! Size of this image
    write(keyword,'(a)') "imagesize"
    write(message,*) "Image size ("//trim(axisUnits)//"): "
    call getReal(keyword, thisImageSize, 1.0, cLine, fLine, nLines, &
         trim(message), "(a,1pe10.2,1x,a)", imageSize, ok, .false.)

    ! Aspect ratio
    write(keyword,'(a)') "imageaspect"
    call getReal(keyword, thisAspectRatio, 1.0, cLine, fLine, nLines, &
         "Image aspect ratio: ", "(a,f4.1,1x,a)", aspectRatio, ok, .false.)

    ! Inclination and position angle
    write(keyword,'(a)') "inclination"
    call getReal(keyword, thisInc, real(degtorad), cLine, fLine, nLines, &
         "Inclination of image: ","(a,f6.1,1x,a)",inclination*real(radtodeg), ok, .false.)
    write(keyword,'(a)') "positionangle"
    call getReal(keyword, thisPA, real(degtorad), cLine, fLine, nLines, &
         "Position angle (deg): ","(a,f6.1,1x,a)", positionAngle*real(radtodeg), ok, .false.)

    ! Position of image centre
    write(keyword,'(a)') "imagecentrex"
    call getReal(keyword, thisOffsetx, 1.0, cLine, fLine, nLines, &
         "Image centre x position (10^10 cm): ", "(a,e10.1,1x,a)", offsetx, ok, .false.)
    write(keyword,'(a)') "imagecentrey"
    call getReal(keyword, thisOffsety, 1.0, cLine, fLine, nLines, &
         "Image centre y position (10^10 cm): ", "(a,e10.1,1x,a)", offsety, ok, .false.)

    do i = 1, nImage

       write(imageFilename,'(a,i4.4,a)') "movie",i,".fits"

       thisTime = 0.d0
       ang = twoPi * dble(i-1)/dble(nImage)
       viewVec = VECTOR(sin(thisInc),0.d0, cos(thisInc))
       viewVec = rotateZ(viewVec, dble(ang))
       if (writeoutput) write(*,*) i, viewvec
       call setImageParams(i, thisLambdaImage, outputimageType,imageFilename, thisnpixels, axisUnits, fluxUnits, &
            thisImageSize, thisaspectRatio, thisInc, thisPA, thisOffsetx, thisOffsety, gridDistance, &
            viewVec=viewVec, imageTime = thisTime)
    enddo

    call writeInfo(" ")

  end subroutine readMovieParameters


  subroutine readBioPhysicsParameters(cLine, fLine, nLines)
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok
    integer :: i
    character(len=20) :: keyword

    call getVector("sourcepos", sourcePosition, 1.d0, cLine, fLine, nLines, &
         "Source position: ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

    call getDouble("sourcetheta", sourceTheta, degtorad, cLine, fLine, nLines, &
         "Source theta angle: ","(a,f5.1,a)",180.d0, ok, .false.)

    call getDouble("sourcephi", sourcePhi, degtorad, cLine, fLine, nLines, &
         "Source phi angle: ","(a,f5.1,a)",0.d0, ok, .false.)

    call getInteger("ncomponent", nComponent, cLine, fLine, nLines, &
         "Number of source components","(a,i2,a)", 1, ok, .false.)

    do i = 1, nComponent

       write(keyword, '(a,i1)') "componenttype",i
       call getString(keyword, componentType(i), cLine, fLine, nLines, &
            "Type of source component: ","(a,a,1x,a)","spectrum.dat", ok, .true.)

       write(keyword, '(a,i1)') "componentpos",i
       call getVector(keyword, componentPosition(i), 1.d0, cLine, fLine, nLines, &
            "Component position: ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

       write(keyword, '(a,i1)') "power",i
       call getDouble(keyword, componentPower(i), 1.d0, cLine, fLine, nLines, &
            "Component power (mW): ","(a,f5.1,1x,a)", 1.0d0, ok, .true.)

       write(keyword, '(a,i1)') "wavelength",i
       call getDouble(keyword, componentWavelength(i), 1.d0, cLine, fLine, nLines, &
            "Component monochromatic wavelength (nm): ","(a,f5.1,1x,a)", 1.0d0, ok, .false.)

       write(keyword, '(a,i1)') "radius",i
       call getDouble(keyword, componentRadius(i), microntocm, cLine, fLine, nLines, &
            "Component radius (microns): ","(a,f5.1,1x,a)", 1.0d0, ok, .true.)

       if (componentType(i) == "fibre") then
          write(keyword, '(a,i1)') "na",i
          call getDouble(keyword, componentNa(i), 1.d0, cLine, fLine, nLines, &
               "Numerical aperture: ","(a,f5.1,1x,a)", 1.0d0, ok, .true.)
          call getString("fibrespectrum", fibreSpectrum, cLine, fLine, nLines, &
                  "Source fibre spectrum: ","(a,a,a)", "none", ok, .true.)
       endif
    enddo

    call getString("wavefrontfile", wavefrontfile, cLine, fLine, nLines, &
            "Wavefront file: ","(a,a,1x,a)","spectrum.dat", ok, .true.)

    call getBigInteger("nphotons", nPhotons, cLine, fLine, nLines, &
         "Number of photons in loop: ", "(a,i15,1x,a)", 100000_bigInt, ok, .false.)


  end subroutine readBioPhysicsParameters


  subroutine readFitsParameters(cLine, fLine, nLines)
    use fits_utils_mod, only: fitsBitpix
    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call getInteger("fitsbitpix", fitsBitpix, cLine, fLine, nLines, &
         "FITS file BITPIX ","(a,i2,a)", -32, ok, .false.)

  end subroutine readFitsParameters

  subroutine readSpectrumParameters(cLine, fLine, nLines)
    use sed_mod, only:  setSedParameters,  SEDlamMin, SEDlamMax, SEDwavLin, SEDnumLam

    character(len=lencLine) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok, isRange
    character(len=80) :: message
    logical :: sed, jansky, SIsed, uniformInCos, lambdaInMicrons
    integer :: nInclination
    real    :: firstInclination, firstPA=0.0
    real    :: lastInclination=80.0, lastPA=0.0
    real, allocatable :: inclinations(:), posangs(:)
    character(len=80) :: outFile


    call getBigInteger("nphotons", nPhotons, cLine, fLine, nLines, &
         "Number of photons in SED: ", "(a,i15,1x,a)", 100000_bigInt, ok, .false.)

    call getBigInteger("nphotspec", nPhotSpec, cLine, fLine, nLines, &
         "Number of photons in SED: ", "(a,i15,1x,a)", nphotons, ok, .false.)

    call getLogical("forcefirstscat", forceFirstScat, cLine, fLine, nLines, &
         "Force first scattering: ","(a,1l,1x,a)", .false., ok, .false.)


    call getInteger("maxscat", maxScat, cLine, fLine, nLines, &
         "Maximum number of scatterings: ","(a,i9,a)", 10000, ok, .false.)

    call getLogical("forcefirstscat", forceFirstScat, cLine, fLine, nLines, &
         "Force first scattering: ","(a,1l,1x,a)", .false., ok, .false.)

    call getInteger("ninc", nInclination, cLine, fLine, nLines, &
         "Number of inclination angles: ", "(a,i3,1x,a)", 1, ok, .false.)

incStyle: if (checkPresent("inclinations", cline, nlines)) then
! Inclinations and postion angles range style
       call parameterDefinesRange("inclinations", firstInclination, lastInclination, isRange, cLine, fLine, nLines)
       if (isRange) then
          write(message,'(a,f6.2,a,f6.2,a)') &
               "SED inclination range is from ", firstInclination, " to ", lastInclination, " degrees"
          call writeInfo(message,TRIVIAL)
          firstInclination = firstInclination * real(degToRad)
          lastInclination = lastInclination * real(degToRad)

          if (checkPresent("posangs", cline, nlines)) then
             call parameterDefinesRange("posangs", firstPA, lastPA, isRange, cLine, fLine, nLines)
             if (isRange) then
                write(message,'(a,f6.2,a,f6.2,a)') &
                     "SED PA range is from ", firstPA, " to ", lastPA, " degrees"
                call writeInfo(message,TRIVIAL)
                firstPA = firstPA * real(degToRad)
                lastPA = lastPA * real(degToRad)
             else
                call writewarning("Position angles and inclinations should be specified in the same format")
             endif
          else
             firstPA=0.0; lastPA=0.0
          endif

       else
! Inclinations and position angles list style
          allocate(inclinations(nInclination))
          call getRealArray("inclinations", inclinations, 1.0, cLine, fLine, nLines, &
               "Inclinations (deg): ",90., ok, .false.)
          inclinations(:) = inclinations(:) * real(degToRad)
          allocate(posangs(nInclination))
          call getRealArray("posangs", posangs, 1.0, cLine, fLine, nLines, &
               "Position angles (deg): ",0., ok, .false.)
          posangs(:) = posangs(:) * real(degToRad)
       end if

    else

       ! Inclinations and postion angles "first and last" style
! lasstinc is required if inclinations keyword is not present and there is >1 inclination to do
       if (nInclination > 1) then
! firstinc is required if inclinations keyword is not present
            call getReal("firstinc", firstInclination, real(degToRad), cLine, fLine, nLines, &
            "First inclination angle (deg): ","(a,f4.1,1x,a)", 10., ok, .true.)
            call getReal("lastinc", lastInclination, real(degToRad), cLine, fLine, nLines, &
            "Last inclination angle (deg): ","(a,f5.1,1x,a)", 80., ok, .true.)
         endif
! Position angles default to zero, we don't need these for 2D geometries
       call getReal("firstPA", firstPA, real(degToRad), cLine, fLine, nLines, &
            "First position angle (deg): ","(a,f4.1,1x,a)", 0.0, ok, .false.)
       call getReal("lastPA", lastPA, real(degToRad), cLine, fLine, nLines, &
            "Last position angle (deg): ","(a,f5.1,1x,a)", 0.0, ok, .false.)

    end if incStyle


    call getLogical("sedincseq", sedIncSeq, cLine, fLine, nLines, &
         "Spectrum filename named in sequence by inclination: ","(a,1l,1x,a)", .false., ok, .false.)

    call getReal("inclination", thisinclination, real(degtorad), cLine, fLine, nLines, &
         "Inclination angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)

    call getReal("positionangle", thisPA, real(degtorad), cLine, fLine, nLines, &
         "Inclination angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)

    call getDouble("ism_av", ism_av, 1.d0, cLine, fLine, nLines, &
         "Interstellar extinction (A_V): ","(a,f5.1,1x,a)", 0.d0, ok, .false.)

    call getDouble("ism_rv", ism_rv, 1.d0, cLine, fLine, nLines, &
         "Interstellar R_V: ","(a,f5.1,1x,a)", 3.1d0, ok, .false.)



    call getDouble("maxVel", maxVel, 1.d0, cLine, fLine, nLines, &
         "Maximum Velocity Channel (km/s): ","(a,f5.1,1x,a)", 1.0d0, ok, .false.)

    call getDouble("minVel", minVel, 1.0_db, cLine, fLine, nLines, &
         "Minimum Velocity Channel (km/s): ","(a,f4.1,1x,a)", -1.0d0*maxVel, ok, .false.)

    call getReal("vmin", vMinSpec, 1.0, cLine, fLine, nLines, &
         "Minimum velocity output to spectrum (km/s)","(a,1PE10.3,1x,a)", -750.0, ok, .false.)

    call getReal("vmax", vMaxSpec, 1.0, cLine, fLine, nLines, &
         "Maximum velocity output to spectrum (km/s)","(a,1PE10.3,1x,a)", 750.0, ok, .false.)

    call getInteger("nv", nv, cLine, fLine, nLines, &
         "Number of velocity bins: ", "(a,i3,1x,a)", 1, ok, .false.)

    call getReal("distance", gridDistance, 1., cLine, fLine, nLines, &
         "Grid distance (pc): ","(a,f6.1,1x,a)", 100., ok, .false.)

    call getInteger("nphase", nPhase, cLine, fLine, nLines, &
         "Number of phases: ", "(a,i3,1x,a)", 1, ok, .false.)

!    call getInteger("nstart", nStartPhase, cLine, fLine, nLines, &
!         "Start at phase: ", "(a,i3,1x,a)", 1, ok, .false.)

!    call getInteger("nend", nEndPhase, cLine, fLine, nLines, &
!         "End at phase: ", "(a,i3,1x,a)", nPhase, ok, .false.)

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

    call getLogical("lambdainmicrons", lambdainmicrons, cLine, fLine, nLines, &
         "Write wavelength array in microns: ","(a,1l,1x,a)", .false., ok, .false.)

    call getReal("sedlammin", SEDlamMin, 1.0e4, cLine, fLine, nLines, &
         "Minimum wavelength output to SED (microns)","(a,1PE10.3,1x,a)", 0.1, ok, .false.)

    call getReal("sedlammax", SEDlamMax, 1.0e4, cLine, fLine, nLines, &
         "Maximum wavelength output to SED (microns)","(a,1PE10.3,1x,a)", 2000.0, ok, .false.)


    call getLogical("sedwavlin", SEDwavLin, cLine, fLine, nLines, &
         "Linear wavelength spacing in SED: ","(a,1l,1x,a)", .false., ok, .false.)

    call getInteger("sednumlam", SEDnumLam, cLine, fLine, nLines, &
         "Number of SED points: ", "(a,i3,1x,a)", 200, ok, .false.)

    call getLogical("lambdafile", lamfile, cLine, fLine, nLines, &
         "Wavelength grid from file: ","(a,1l,1x,a)", .false., ok, .false.)

    if (lamFile) then
       call getString("lamfilename", lamFilename, cLine, fLine, nLines, &
            "Wavelength grid filename: ","(a,a,1x,a)","none", ok, .true.)
    endif

    call getLogical("resolvesilicate", resolveSilicateFeature, cLine, fLine, nLines, &
         "Add wavelength points to resolve silicate feature in SED: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("freefreesed", freefreeSED, cLine, fLine, nLines, &
         "Include free-free emission in SED: ","(a,1l,1x,a)", .false., ok, .false.)


    call getLogical("recombinationsed", recombinationSED, cLine, fLine, nLines, &
         "Include recombination emission in SED: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("forbiddensed", forbiddenSED, cLine, fLine, nLines, &
         "Include forbidden line emission in SED: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("sedcosspaced", UniformInCos, cLine, fLine, nLines, &
         "SED inclinations are uniform in cos(inc): ","(a,1l,1x,a)", .true., ok, .false.)

    if (allocated(inclinations)) then
       call setSedParameters(outFile,jansky,SIsed,sed, lambdaInMicrons, firstPA, lastPA, incList=inclinations, PAlist=posangs)
       deallocate(inclinations)
       deallocate(posangs)
    else
       if (nInclination > 1) then
          call setSedParameters(outFile,jansky,SIsed,sed, lambdaInMicrons, firstPA, lastPA, nInclination=nInclination,&
               firstInc=firstInclination,LastInc=LastInclination, cosSpacing=UniformInCos)
       else
          call setSedParameters(outFile,jansky,SIsed,sed, lambdaInMicrons, firstPA, lastPA, nInclination=nInclination,&
               firstinc=FirstInclination,lastInc=LastInclination,thisInclination=thisInclination, thisPA=thisPA)
       endif
    end if

  end subroutine readSpectrumParameters


!-----------------------------------------------------------------------------------------


subroutine findReal(name, value, unitString, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 real :: value
 character(len=lencLine) :: cLine(:)
 character(len=80) :: restOfLine
 character(len=*) :: unitString
 logical :: fLine(:)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  k = index(cline(i)," ")-1
  j = len_trim(name)
  if ((trim(cLine(i)(1:k)) .eq. name(1:j)).and.(j==k)) then
     if (.not.ok) then
        ok = .true.
        read(cLine(i)(j+1:lencline),*) value
        restOfLine = adjustl(cLine(i)(j+1:))
        k = index(restOfLine," ")
        restOfLine = trim(restOfLine(k+1:))
        read(restOfLine,'(a)') unitString
        fLine(i) = .true.
     else
        call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
        call torus_stop
     endif
  endif
 end do
 end subroutine findReal

 subroutine findDouble(name, value, unitString, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 real(double) :: value
 character(len=lencLine) :: cLine(:)
 character(len=80) :: restOfLine
 character(len=*) :: unitString
 logical :: fLine(:)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  k = index(cline(i)," ")-1
  j = len_trim(name)
  if ((trim(cLine(i)(1:k)) .eq. name(1:j)).and.(j==k)) then
     if (.not.ok) then
        ok = .true.
        read(cLine(i)(j+1:lencline),*) value
        restOfLine = adjustl(cLine(i)(j+1:))
        k = index(restOfLine," ")
        restOfLine = trim(restOfLine(k+1:))
        read(restOfLine,'(a)') unitString
        fLine(i) = .true.
     else
        call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
        call torus_stop
     endif
  endif
 end do
 end subroutine findDouble

 subroutine findUnitDouble(name, value, unit, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 real(double) :: value
 character(len=lencLine) :: cLine(:)
 character(len=*) :: unit
 logical :: fLine(:)
 integer :: nLines
 logical :: ok
 integer :: i, j, k, readstat

 ok = .false.
 do i = 1, nLines
  j = len_trim(name)
  k = index(cline(i)," ")-1
  if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
     if (.not.ok) then
        ok = .true.
        read(cLine(i)(j+1:lencline),*,iostat=readstat) value, unit
        if(readstat /= 0) then
           read(cLine(i)(j+1:lencline),*,iostat=readstat) value
           if(readstat /= 0) then
              write(*,*) "parameters.dat entry with no value given!!" //name
              call torus_stop
           else
              unit = "default"
           end if
        end if
        fLine(i) = .true.
     else
     call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
     call torus_stop
  endif
endif
end do
 end subroutine findUnitDouble

 subroutine findUnitVECTOR(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 type(VECTOR) :: value
 character(len=10) :: unit
 integer :: readStat
 character(len=lencLine) :: cLine(:)
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
          read(cLine(i)(j+1:lencline),*,iostat=readstat) value, unit
          if(readstat /= 0) then
           read(cLine(i)(j+1:lencline),*,iostat=readstat) value
           if(readstat /= 0) then
              write(*,*) "parameters.dat entry with no value given!!" //name
              call torus_stop
           else
              unit = "default"
           end if
        end if
        fLine(i) = .true.
     else
        call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
        call torus_stop
     endif
  endif
end do
end subroutine findUnitVECTOR

 subroutine findVECTOR(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 type(VECTOR) :: value
 character(len=lencLine) :: cLine(:)
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
        read(cLine(i)(j+1:lencline),*) value%x, value%y, value%z
        fline(i) = .true.
  else
     call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
     call torus_stop
  endif
     endif


 end do
end subroutine findVECTOR

subroutine findInteger(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 integer :: value
 character(len=lencLine) :: cLine(:)
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
        call torus_stop
     endif
  endif
end do
end subroutine findInteger

 subroutine findBigInteger(name, value, cLine, fLine, nLines, ok)
   implicit none
   character(len=*) :: name
   logical :: fLine(:)
   integer(kind=bigInt) :: value
   character(len=lencLine) :: cLine(:)
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
        call torus_stop
     endif
  endif
   end do
 end subroutine findBigInteger

subroutine findLogical(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 logical :: value
 character(len=lencLine) :: cLine(:)
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
        call torus_stop
     endif
  endif
end do
 end subroutine findLogical

subroutine findString(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 character(len=*) :: value
 character(len=lencLine) :: cLine(:), tmp
 logical :: fLine(:)
 integer :: nLines
 logical :: ok
 integer :: i, j, k, n

 ok = .false.
 do i = 1, nLines
  j = len_trim(name)
  k = index(cline(i)," ")-1
     if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
  if (.not.ok) then
        ok = .true.
        read(cLine(i)(j+1:),'(a)') tmp
        tmp = trim(adjustl(tmp))
        n = index(tmp, " ")
        if (n == 0) then
           value = tmp
        else
           value = tmp(1:n)
        endif
        fLine(i) = .true.
  else
     call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
     call torus_stop
  endif
endif
 end do
 value = trim(value)
 end subroutine findString

! Return the whole line containing a given keyword
subroutine findLine(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 character(len=lencLine) :: value
 character(len=lencLine) :: cLine(:)
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
        value=(cLine(i))
        fLine(i) = .true.
  else
     call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
     call torus_stop
  endif
endif
 end do
 value = trim(value)
 end subroutine findLine

subroutine findRealArray(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 real :: value(:)
 character(len=lencLine) :: cLine(:)
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
     call torus_stop
  endif
endif
 end do
 end subroutine findRealArray

subroutine findIntegerArray(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 integer :: value(:)
 character(len=lencLine) :: cLine(:)
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
     call torus_stop
  endif
endif
 end do
 end subroutine findIntegerArray

subroutine findDoubleArray(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 real(double) :: value(:)
 character(len=lencLine) :: cLine(:)
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
     call torus_stop
  endif
endif
 end do
end subroutine findDoubleArray

 subroutine getInteger(name, ival, cLine, fLine, nLines, message, format, idef, ok, &
                       musthave)
  character(len=*) :: name
  integer :: ival
  character(len=lencLine) :: cLine(:)
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
       call torus_stop
    endif
    ival = idef
    default = " (default)"
  endif
  if (ok) then
     write(output,format) trim(message),ival,default
     call writeInfo(output, TRIVIAL)
  endif
 end subroutine getInteger

 subroutine getBigInteger(name, ival, cLine, fLine, nLines, message, format, idef, ok, &
                       musthave)
  character(len=*) :: name
  integer(kind=bigInt) :: ival,idef
  character(len=lencLine) :: cLine(:)
  logical :: fLine(:)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default
  logical :: musthave
  logical :: ok
  ok = .true.
  default = " "
  call findBigInteger(name, ival, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput) write(*,'(a,a)') name, " must be defined"
       call torus_stop
    endif
    ival = idef
    default = " (default)"
  endif
  if (ok) then
     write(output,format) trim(message),ival,default
     call writeInfo(output, TRIVIAL)
  endif
end subroutine getBigInteger

 subroutine getRealWithUnits(name, rval, defaultUnits, requiredUnits, cLine, fLine, nLines, message, rdef, ok, &
                    musthave)
  character(len=*) :: name, requiredUnits, defaultUnits
  real :: rval
  real :: thisUnitConversion
  logical :: musthave
  character(len=lencLine) :: cLine(:)
  logical :: fLine(:)
  character(len=100) :: output
  character(len=20) :: unitString, cval
  integer :: nLines
  character(len=*) :: message
  character(len=10) :: default
  real :: rdef
  logical :: ok
  integer :: i, j
  ok = .true.
  thisUnitConversion = 1.
  default = " "
  call findReal(name, rval, unitString, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput) write(*,'(a,a)') name, " must be defined"
       call torus_stop
    endif
    rval = rdef
    default = " (default)"
 else
    if (len(trim(unitString)) .ne. 0) then
       i = unitNumber(unitString)
       j = unitNumber(requiredUnits)
       thisUnitConversion = unitList(i)%unitInCGS/unitList(j)%unitInCGS
    else
       i = unitNumber(defaultUnits)
       j = unitNumber(requiredUnits)
       thisUnitConversion = unitList(i)%unitInCGS/unitList(j)%unitInCGS
    endif
 endif 
 if (ok) then
    call prettyNumber(rval, cval)
    write(output,'(a)') trim(message)//" "//trim(cval)//" "//trim(unitList(i)%longUnitName) // trim(default)
    call writeInfo(output, TRIVIAL)
 endif
 rval = rval * thisunitConversion
end subroutine getRealWithUnits

 subroutine getDoubleWithUnits(name, rval, defaultUnits, requiredUnits, cLine, fLine, nLines, message, rdef, ok, &
                    musthave)
  character(len=*) :: name, requiredUnits, defaultUnits
  real(double) :: rval
  real(double) :: thisUnitConversion
  logical :: musthave
  character(len=lencLine) :: cLine(:)
  logical :: fLine(:)
  character(len=100) :: output
  character(len=20) :: unitString, cval
  integer :: nLines
  character(len=*) :: message
  character(len=10) :: default
  real(double) :: rdef
  logical :: ok
  integer :: i, j
  ok = .true.
  thisUnitConversion = 1.
  default = " "
  call findDouble(name, rval, unitString, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput) write(*,'(a,a)') name, " must be defined"
       call torus_stop
    endif
    rval = rdef
    default = " (default)"
 else
    if (len(trim(unitString)) .ne. 0) then
       i = unitNumber(unitString)
       j = unitNumber(requiredUnits)
       thisUnitConversion = unitList(i)%unitInCGS/unitList(j)%unitInCGS
    else
       i = unitNumber(defaultUnits)
       j = unitNumber(requiredUnits)
       thisUnitConversion = unitList(i)%unitInCGS/unitList(j)%unitInCGS
    endif
 endif
 if (ok) then
    call prettyNumber(rval, cval)
    write(output,'(a)') trim(message)//" "//trim(cval)//" "//trim(unitList(i)%longUnitName) // trim(default)
    call writeInfo(output, TRIVIAL)
 endif
 rval = rval * thisunitConversion
end subroutine getDoubleWithUnits

 logical function checkPresent(keyword, cline, nlines)
   character(len=*) :: keyword
   character(len=lencLine) :: cLine(:)
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


 subroutine getDouble(name, dval, unitConversion, cLine, fLine, nLines, message, &
   format, ddef, ok, musthave)
  character(len=*) :: name
  real(double) :: dval, unitConversion
  logical :: musthave
  character(len=lencLine) :: cLine(:)
  logical :: fLine(:)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default, units
  real(double) :: ddef
  logical :: ok
  ok = .true.
  default = " "
  call findDouble(name, dval, units, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput) write(*,'(a,a)') name, " must be defined"
       call torus_stop
    endif
    dval = ddef
    default = " (default)"
 endif
 if (ok) then
    write(output,format) trim(message)//" ",dval,default
    call writeInfo(output, TRIVIAL)
 endif
 dval = dval  * unitConversion
 end subroutine getDouble


  subroutine getReal(name, val, unitConversion, cLine, fLine, nLines, message, &
   format, def, ok, musthave)
  character(len=*) :: name
  real :: val, unitConversion
  logical :: musthave
  character(len=lencLine) :: cLine(:)
  logical :: fLine(:)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default, units
  real :: def
  logical :: ok
  ok = .true.
  default = " "
  call findReal(name, val, units, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput) write(*,'(a,a)') name, " must be defined"
       call torus_stop
    endif
    val = def
    default = " (default)"
 endif
 if (ok) then
    write(output,format) trim(message)//" ",val,default
    call writeInfo(output, TRIVIAL)
 endif
 val = val  * unitConversion
 end subroutine getReal


 subroutine getUnitDouble(name, dval, unitType, cLine, fLine, nLines, message, &
   format, ddef, ok, musthave)
  character(len=*) :: name
  real(double) :: dval
  logical :: musthave
  character(len=lencLine) :: cLine(:)
  logical :: fLine(:)
  character(len=100) :: output
  character(len=10) :: unitString
  character(len=*) :: unitType
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default
  real(double) :: ddef
  logical :: ok
  ok = .true.
  default = " "
  call findUnitDouble(name, dval, unitString, cLine, fLine, nLines, ok)
  if(unitString == "default") then
!     write(*,*) "Unit not given, proceeding with default for ", unitType
  end if
  if (.not. ok) then
    if (musthave) then
       if (writeoutput) write(*,'(a,a)') name, " must be defined"
       call torus_stop
    endif
    dval = ddef
    default = " (default)"
 endif
 call convertToTorusUnits(unitString, unitType,  dval)
 if (ok) then
    write(output,format) trim(message),dval,default
    call writeInfo(output, TRIVIAL)
 endif
 end subroutine getUnitDouble

 subroutine getVector(name, dval, unitConversion, cLine, fLine, nLines, message, format, ddef, ok, &
                    musthave)
  character(len=*) :: name
  real(double) :: unitConversion
  type(VECTOR) :: dval, ddef
  logical :: musthave
  character(len=lencLine) :: cLine(:)
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
       call torus_stop
    endif
    dval = ddef
    default = " (default)"
 endif
 if (ok) then
    write(output,format) trim(message),dval,default
    call writeInfo(output, TRIVIAL)
 endif
 dval = dval  * unitConversion
end subroutine getVector

 subroutine getVectorWithUnits(name, dval, defaultUnits, requiredUnits, cLine, fLine, nLines, message, format, ddef, ok, &
                    musthave)
  character(len=*) :: name, requiredUnits, defaultUnits
  real(double) :: unitConversion, thisUnitConversion
  type(VECTOR) :: dval, ddef
  logical :: musthave
  character(len=10) :: unitString
  character(len=lencLine) :: cLine(:)
  logical :: fLine(:)
  integer :: i, j
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
       call torus_stop
    endif
    dval = ddef
    default = " (default)"
 else
    if (len(trim(unitString)) .ne. 0) then
       i = unitNumber(unitString)
       j = unitNumber(requiredUnits)
       thisUnitConversion = unitList(i)%unitInCGS/unitList(j)%unitInCGS
    else
       i = unitNumber(defaultUnits)
       j = unitNumber(requiredUnits)
       thisUnitConversion = unitList(i)%unitInCGS/unitList(j)%unitInCGS
    endif
 endif

 if (ok) then
    write(output,format) trim(message),dval,default
    call writeInfo(output, TRIVIAL)
 endif
 dval = dval  * unitConversion
end subroutine getVectorWithUnits


 subroutine getString(name, rval, cLine, fLine, nLines, message, format, rdef, ok, &
                      musthave)
  character(len=*) :: name
  character(len=*) :: rval
  logical :: musthave
  character(len=lencLine) :: cLine(:)
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
       call torus_stop
    endif
    rval = rdef
    default = " (default)"
  endif
  if (ok) then
     write(output,format) trim(message)//" ",trim(rval),default
     call writeInfo(output, TRIVIAL)
  endif
 end subroutine getString


 subroutine getLogical(name, rval, cLine, fLine, nLines, message, cformat, rdef, ok, &
                       musthave)
  character(len=*) :: name
  logical :: rval
  logical :: musthave
  character(len=lencLine) :: cLine(:)
  logical :: fLine(:)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, cformat
  character(len=10) :: default
  logical :: rdef
  character(len=6) :: trueOrFalse, cf
  logical :: ok, thisIsDefault
  character (len=80) :: errorMessage

  cf = cformat
  ok = .true.
  default = " "
  thisIsDefault = .false.
  call findLogical(name, rval, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       write(errorMessage,'(a,a)') name, " must be defined"
       call writeFATAL(errorMessage)
       call torus_stop
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
  character(len=lencLine) :: cLine(:)
  logical :: fLine(:)
  character(len=lencLine+200) :: output
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
       call torus_stop
    endif
    rval(:) = rdef
    default = " (default)"
  endif
  write(output,*) trim(message)//" ",rval,default
  call writeInfo(output, TRIVIAL)
  rVal = rVal * unitConversion
 end subroutine getRealArray

 subroutine getDoubleArray(name, rval, unitConversion, cLine, fLine, nLines, message, rdef, ok, &
                      musthave)
  character(len=*) :: name
  real(double) :: rval(:), unitConversion
  logical :: musthave
  character(len=lencLine) :: cLine(:)
  logical :: fLine(:)
  character(len=lencLine+200) :: output
  integer :: nLines
  character(len=*) :: message
  character(len=10) :: default
  real(double) :: rdef
  logical :: ok
  ok = .true.
  default = " "
  call findDoubleArray(name, rval, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput)  write(*,'(a,a)') name, " must be defined"
       call torus_stop
    endif
    rval(:) = rdef
    default = " (default)"
  endif
  write(output,*) trim(message)//" ",rval,default
  call writeInfo(output, TRIVIAL)
  rVal = rVal * unitConversion
 end subroutine getDoubleArray

 subroutine getIntegerArray(name, rval,  cLine, fLine, nLines, message, rdef, ok, &
                      musthave)
  character(len=*) :: name
  integer :: rval(:)
  logical :: musthave
  character(len=lencLine) :: cLine(:)
  logical :: fLine(:)
  character(len=lencLine+200) :: output
  integer :: nLines
  character(len=*) :: message
  character(len=10) :: default
  integer :: rdef
  logical :: ok
  ok = .true.
  default = " "
  call findIntegerArray(name, rval, cLine, fLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput)  write(*,'(a,a)') name, " must be defined"
       call torus_stop
    endif
    rval(:) = rdef
    default = " (default)"
  endif
  write(output,*) trim(message)//" ",rval,default
  call writeInfo(output, TRIVIAL)
end subroutine getIntegerArray

!Return boundary code for input boundary strings
 integer function getBoundaryCode(boundaryString)
   character(len=20) :: boundaryString

   select case (boundaryString)
      case("reflecting")
         getBoundaryCode = 1
      case("periodic")
         getBoundaryCode = 2
      case("shock")
         getBoundaryCode = 3
      case("freeOutNoIn")
         getBoundaryCode = 4
      case("constDenNoVel")
         getBoundaryCode = 5
      case("inflow")
         getBoundaryCode = 6
      case("inflowGrad")
         getBoundaryCode = 7
      case("gasmix")
         getBoundaryCode = 9
      case("slip-wall")
         getBoundaryCode = 10
      case DEFAULT
         print *, "Unrecognised boundary string:", boundaryString
         print *, "Halting."
         call torus_stop
   end select


 end function getBoundaryCode

! Count the number of lines in a file.
! Should return the same answer as wc -l i.e. includes blank lines in the total
! D. Acreman, June 2010
  integer function file_line_count(filename)

    character(len=*) :: filename
    character(len=1) :: dummy
    integer :: status

    file_line_count = 0
    open(unit=30, status="old", file=filename)
    do
       read(30,'(a1)',iostat=status) dummy
       if ( status /= 0 ) exit
       file_line_count = file_line_count + 1
    end do
    close(30)

  end function file_line_count

  subroutine testFileExists(thisFile)
    character(len=*) :: thisFile
    integer :: error
    open(53, file=thisFile, status="old",form="formatted", iostat=error)

    if (error /= 0) then
       call writeFatal("Error opening file: "//trim(thisFile))
       call torus_stop
    else
       close(53)
    endif
  end subroutine testFileExists

  subroutine parameterDefinesRange(name, firstVal, lastVal, isRange, cLine, fLine, nLines)

    character(len=*), intent(in)  :: name
    real, intent(out)             :: firstVal
    real, intent(out)             :: lastVal
    logical, intent(out)          :: isRange
    character(len=lencLine), intent(in) :: cLine(:)
    logical, intent(in)           :: fLine(:)
    integer, intent(in)           :: nLines

    character(len=lencline) :: thisString
    logical           :: ok
    integer           :: i, j

    call findLine(name, thisString, cLine, fLine, nLines, ok)

! Look for .. to determine if this is a range of values
! The parameters file is lencLine in length
    isRange = .false.
    do i=1, lencLine-4
       if (thisString(i:i+1) == "..") then
          isRange = .true.
          thisString(i:i+1) = "  "
          exit
       end if
    end do

    if (isRange) then
       j = len_trim(name)
       read(thisString(j+1:lencLine),*) firstVal, lastVal
    else
       firstVal = -1.0e-33
       lastVal  =  -1.0e-33
       return
    end if

  end subroutine parameterDefinesRange




  
end module inputs_mod
