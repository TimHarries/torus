module inputs_mod

  use vector_mod
  use unix_mod
  use messages_mod
  use kind_mod
  use constants_mod
  use units_mod
  use zlib_mod, only : uncompressedDumpFiles, buffer, nbuffer
  implicit none

  include "input_variables.f90"

contains
  
  subroutine  inputs()

    implicit none

    integer :: nLines
    logical :: ok
    logical :: compresseddumpfiles
    character(len=80), allocatable :: cLine(:) 
    character(len=80) :: message, thisLine
    logical :: done
    logical, allocatable :: fLine(:)

    ! character vars for unix environ

    character(len=80) :: dataDirectory

    character(len=100) :: paramFile

    integer :: error

    datadirectory = " "
    done = .false.
    ok = .true.
    oneKappa = .false.
    monteCarloRT = .false.

    biasToLyman = .false.

    nDustType = 1
    tMinGlobal = 3.
    filter_set_name = "natural"
    noscattering = .false.
    forceFirstScat = .false.
    usebias = .true.
    
    intProFilename = "none"
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

    call testFileExists(paramFile)

    nLines = file_line_count(paramfile)
    ! nLines+1 element of cline is used in loop which reads parameter file
    ! Needs to be +2 in case the last line doesn't have a newline at the end
    allocate ( cLine(nLines+2) )
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


    call getLogical("binaryxml", useBinaryXMLVTKfiles, cLine, fLine, nLines, &
         "Use binary XML VTK files: ","(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("compressdumps", compressedDumpFiles, cLine, fLine, nLines, &
         "Use compressed dump files: ","(a,1l,1x,a)", .false., ok, .false.)
    uncompressedDumpFiles = .not.compressedDumpFiles
    nbuffer = 0
    buffer = 0

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
            "Number of threads for domain decomposition: ","(a,i3,a)", 0, ok, .false.)
    endif



    call getLogical("debug", debug, cLine, fLine, nLines, &
         "Output debug information: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("setupamr", doSetupAMRgrid, cLine, fLine, nLines, &
         "Set up or readin and AMR grid: ","(a,1l,1x,a)", .true., ok, .false.)


    call getLogical("multimodels", multimodels, cLine, fLine, nLines, &
         "Perform calculation on multiple models: ","(a,1l,1x,a)", .false., ok, .false.)

    nModelStart = 1
    nModelEnd = 1
    iModel = 1
    if (multimodels) then
       call getInteger("nstart", nModelStart, cLine, fLine, nLines, &
            "Starting number for model: ","(a,i7,a)", 1, ok, .true.)
       call getInteger("nend", nModelEnd, cLine, fLine, nLines, &
            "Starting number for model: ","(a,i7,a)", 1, ok, .true.)
    endif


    call getLogical("readgrid", readGrid, cLine, fLine, nLines, &
         "Read grid file: ","(a,1l,1x,a)", .false., ok, .true.)

!    if (multimodels.and.(.not.readgrid)) then
!       write(message,*) "Multiple models specified by readgrid set to false"
!       call writeFatal(message)
!       stop
!    endif

    
    call getLogical("singlemegaphoto", singleMegaPhoto, cLine, fLine, nLines, &
         "Do a long photoionization calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    if(readgrid) call getLogical("justdump", justDump, cLine, fLine, nLines, &
         "Dump a vtk file and abort: ","(a,1l,1x,a)", .false., ok, .false.)

    if (readgrid) call getString("inputfile", gridInputFilename, cLine, fLine, nLines, &
                  "Grid input filename: ","(a,a,1x,a)","none", ok, .true.)
    
!    if (multimodels.and.(index(gridInputFilename, "*") == 0)) then
!       write(message,*) "Multiple models but input filename has no asterixes"
!       call writeFatal(message)
!       stop
!    endif

    call readGridInitParameters(cLine, fLine, nLines)

    call readGeometrySpecificParameters(cLine, fLine, nLines)



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

    call getLogical("nbodyphysics", nBodyPhysics, cLine, fLine, nLines, &
         "Include n-body physics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getDouble("timeunit", timeUnit, 1.d0, cLine, fLine, nLines, &
         "Code unit of time: ","(a,e12.3,1x,a)", 1.d0, ok, .false.)

    call getDouble("massunit", massUnit, 1.d0, cLine, fLine, nLines, &
         "Code unit of mass: ","(a,e12.3,1x,a)", 1.d0, ok, .false.)

    call getDouble("lengthunit", lengthUnit, 1.d0, cLine, fLine, nLines, &
         "Code unit of length: ","(a,e12.3,1x,a)", 1.d0, ok, .false.)

    call getInteger("nmonte", inputnMonte, cLine, fLine, nLines, &
         "Number of photon packets","(a,i8,a)", 0, ok, .false.)

    call getInteger("maxiter", maxPhotoionIter, cLine, fLine, nLines, &
         "Maximum number of iterations","(a,i8,a)", 10, ok, .false.)

    if (nBodyPhysics) then
       call getUnitDouble("tend", tEnd, "time", cLine, fLine, nLines, &
            "End time for calculation: ","(a,e12.3,1x,a)", 1.d10, ok, .false.)

       call getUnitDouble("tdump", tDump, "time", cLine, fLine, nLines, &
            "Time between dump files: ","(a,e12.3,1x,a)", 0.d0, ok, .false.)

       call getDouble("eps", inputEps, 1.d0, cLine, fLine, nLines, &
            "Gravity softening length (cm): ","(a,e12.3,1x,a)", 0.d0, ok, .false.)

       call getLogical("addsinks", addSinkParticles, cLine, fLine, nLines, &
            "Add sink particles: ", "(a,1l,1x,a)", .true., ok, .false.)

       call getLogical("hosokawa", hosokawaTracks, cLine, fLine, nLines, &
            "Use Hosokawa evolutionary tracks: ", "(a,1l,1x,a)", .true., ok, .false.)

       call getLogical("mergesinks", mergeBoundSinks, cLine, fLine, nLines, &
            "Merge bound sinks: ", "(a,1l,1x,a)", .true., ok, .false.)

       call getLogical("movesinks", moveSources, cLine, fLine, nLines, &
            "Allow sources to move: ", "(a,1l,1x,a)", .true., ok, .false.)

    endif

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

    call getLogical("hydrodynamics", hydrodynamics, cLine, fLine, nLines, &
         "Perform a hydrodynamics calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("nbody", donBodyOnly, cLine, fLine, nLines, &
         "Perform an n-body (bigG=1) calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("radiationHydrodynamics", radiationHydrodynamics, cLine, fLine, nLines, &
         "Perform a radiation-hydrodynamics calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("doselfgrav", doselfgrav, cLine, fLine, nLines, &
         "Use self gravity: ","(a,1l,1x,a)", .false., ok, .false.)

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

    call getLogical("readsources", readsources, cLine, fLine, nLines, &
         "Read sources from a file: ","(a,1l,1x,a)", .false., ok, .false.)

    if (.not. readsources) then
       if (checkPresent("nsource", cLine, nLines)) then
          call readSourceParameters(cLine, fLine, nLines)
       endif
    else
       call getString("sourcefile", sourceFilename, cLine, fLine, nLines, &
                  "Source filename: ","(a,a,1x,a)","none", ok, .true.)
    endif


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

    call getReal("smoothfactor", smoothFactor, 1.0, cLine, fLine, nLines, &
         "Inter-cell maximum ratio before smooth: ","(a,f6.1,1x,a)", 3., ok, .false.)

    call getBigInteger("maxmemory", maxMemoryAvailable, cLine, fLine, nLines, &
         "Maximum memory available (Mb): ","(a,i12,a)", 1500, ok, .false.)
    maxMemoryAvailable = maxMemoryAvailable * 1000000

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

    call getLogical("photometry", calcPhotometry, cLine, fLine, nLines, &
         "Calculate a data cube: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("image", calcImage, cLine, fLine, nLines, &
         "Calculate an image: ","(a,1l,1x,a)", .false., ok, .false.)


    call getLogical("spectrum", calcSpectrum, cLine, fLine, nLines, &
         "Calculate a spectrum: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("sourcehistory", sourceHistory, cLine, fLine, nLines, &
         "Write out the source history: ","(a,1l,1x,a)", .false., ok, .false.)


    call getLogical("writeradialfile", dowriteRadialFile, cLine, fLine, nLines, &
         "Write out a radial dump: ","(a,1l,1x,a)", .false., ok, .false.)

    if (dowriteRadialFile) then
       call getString("radialfilename", radialFilename, cLine, fLine, nLines, &
            "Radial filename: ","(a,a,1x,a)",  " ", ok, .true.)
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
    if (calcDataCube.or.calcImage) call readFitsParameters(cLine, fLine, nLines)
    if (calcSpectrum) call readSpectrumParameters(cLine, fLine, nLines)
    if (calcPhotometry) call readPhotometryParameters(cLine, fLine, nLines)

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
    character(len=80) :: cLine(:), message
    integer :: nLines
    logical :: fLine(:)
    logical :: ok
    logical :: setSubRadius
    character(len=20) :: heightLabel, betaLabel
    integer :: i

    select case(geometry)

       case("bondi")
          call getVector("bondicen", bondiCentre, 1.d0, cLine, fLine, nLines, &
               "Centre of bondi accretion (10^10 cm): ", &
               "(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

       case("fractal")
          call getReal("rho0", rho0, real(mhydrogen),cLine, fLine, nLines, &
               "Initial number density: ","(a,f6.1,1x,a)", 100., ok, .true.)


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
             stop
          endif
       endif


    case("starburst")
       call getDouble("mstarburst", mStarburst, 1.d0, cLine, fLine, nLines, &
            "Starburst mass (solar masses): ","(a,f6.1,a)", 1000.d0, ok, .true.)

       call getDouble("clusterradius", clusterRadius, pctocm, cLine, fLine, nLines, &
            "Cluster radius (pc): ","(a,f5.1,a)", 1.d0, ok, .true.)


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

     case("unisphere","gravtest")
       call getDouble("mass", sphereMass, msol, cLine, fLine, nLines, &
            "Sphere mass (in M_sol): ","(a,f7.1,1x,a)", 1.d-30, ok, .true.)

       call getDouble("radius", sphereRadius, pctocm/1.d10, cLine, fLine, nLines, &
            "Sphere radius (in pc): ","(a,e12.3,1x,a)", 1.d-30, ok, .true.)

       call getVector("position", spherePosition, 1.d0, cLine, fLine, nLines, &
            "Sphere position (10^10 cm): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

       call getVector("velocity", sphereVelocity, 1.d5/cspeed, cLine, fLine, nLines, &
            "Sphere velocity (km/s): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

     case("sphere")
       call getDouble("mass", sphereMass, msol, cLine, fLine, nLines, &
            "Sphere mass (in M_sol): ","(a,f7.1,1x,a)", 1.d-30, ok, .true.)

       call getDouble("radius", sphereRadius, pctocm/1.d10, cLine, fLine, nLines, &
            "Sphere radius (in pc): ","(a,e12.3,1x,a)", 1.d-30, ok, .true.)

       call getVector("position", spherePosition, 1.d0, cLine, fLine, nLines, &
            "Sphere position (10^10 cm): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

       call getVector("velocity", sphereVelocity, 1.d5/cspeed, cLine, fLine, nLines, &
            "Sphere velocity (km/s): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

       call getReal("beta", beta, 1.0, cLine, fLine, nLines, &
            "Density power law: ","(a,f7.2,a)",-2.0, ok, .true.)

       call getDouble("omega", omega, 1.d0, cLine, fLine, nLines, &
            "Angular frequency of rotation: ","(a,f7.2,a)",1.d-13, ok, .true.)

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
       if (ttauriMagnetosphere) then
          call getReal("ttaurirouter", TTauriRouter, TTaurirStar, cLine, fLine, nLines, &
               "T Tauri outer flow radius (in R_star): ","(a,f7.1,1x,a)", 3.0, ok, .true.)
          call getReal("ttauririnner", TTauriRinner, TTaurirStar, cLine, fLine, nLines, &
            "T Tauri inner flow radius (in R_star): ","(a,f7.1,1x,a)", 2.2, ok, .true.)

          call getReal("thotspot", thotspot, 1., cLine, fLine, nLines, &
               "Hot spot temperature (K): ","(a,f8.1,1x,a)", 0., ok, .false.)
       endif

       call getDouble("holeradius", holeRadius, dble(ttaurirstar), cLine, fLine, nLines, &
            "Size of hole in geometrically thin, optically thick disc (in R_star): ","(a,f7.1,1x,a)", &
            0.d0, ok, .true.)



       call getDouble("lxoverlbol", lxOverLBol, 1.d0, cLine, fLine, nLines, &
            "X-ray luminosity  (Bolometric luminosities): ","(a,f7.1,1x,a)", 0.d0, ok, .false.)


       call getDouble("maxcellmass", maxCellMass, 1.d0, cLine, fLine, nLines, &
            "Maximum cell mass after splitting (g): ","(a,e12.5,1x,a)", 1.d30, ok, .false.)

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

          call getReal("epsilon", epsilonDisc, 1., cLine, fLine, nLines, &
               "Isella and Natta epsilon value: ","(a,f6.4,a)", 1., ok, .false.)

          call getReal("alphadisc", alphaDisc, 1., cLine, fLine, nLines, &
               "Disc alpha parameter: ","(a,f5.3,a)", 2.25, ok, .true.)

          call getReal("betadisc", betaDisc, 1., cLine, fLine, nLines, &
               "Disc beta parameter: ","(a,f5.3,a)", 1.25, ok, .true.)

          call getReal("disctemp", alphaDiscTemp, 1., cLine, fLine, nLines, &
               "Disc temperature inside sub radius: ","(a,f8.1,a)", 10000., ok, .true.)

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


       if (tTauriMagnetosphere) then
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
       endif
       call getLogical("ttauriwarp", ttauriwarp, cLine, fLine, nLines, &
            "Include warped disc around magnetosphere: ","(a,1l,1x,a)", .false., ok, .false.)

       call getDouble("hoverr", hoverr, 1.d0, cLine, fLine, nLines, &
            "Warped disc H/R: ","(a,f7.4,1x,a)", 0.3d0, ok, .false.)

       call getLogical("lte", lte, cLine, fLine, nLines, &
            "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
       call getReal("dipoleoffset", dipoleOffset, real(degtorad), cLine, fLine, nLines, &
            "Dipole offset (degrees): ","(a,f7.1,1x,a)", 1.e14, ok, .true.)
       call getLogical("enhance", enhance, cLine, fLine, nLines, &
            "Accretion enhancement: ","(a,1l,1x,a)", .false., ok, .false.)
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
               "Disc wind:: Outer radius of the disc wind [disc wind inner radii]: ", &
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
          DW_theta = DW_theta * degtoRad
       endif


       if (useHartmannTemp .and. isoTherm) then 
          if (writeoutput)  write(*,'(a)') "WARNING: useHartmannTemp and isoTherm both specified!"
          stop
       end if

       if (.not.(useHartmannTemp .or. isoTherm)) then 
          if (writeoutput)  write(*,'(a)') "WARNING: neither useHartmannTemp nor isoTherm specified!"
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


       case("molcluster", "theGalaxy", "cluster", "wr104")
          call getString("sphdatafilename", sphdatafilename, cLine, fLine, nLines, &
               "Input sph data file: ","(a,a,1x,a)","sph.dat.ascii", ok, .true.)

          call getReal("hcritPercentile", hcritPercentile, 1., cLine, fLine, nLines, &
               "Percentile for hcrit: ", "(a,f10.4,1x,f10.4)", 0.80, ok, .false.)

          call getReal("hmaxPercentile", hmaxPercentile, 1., cLine, fLine, nLines, &
               "Percentile for hmax: ", "(a,f10.4,1x,f10.4)", 0.99, ok, .false.)

          call getReal("sphNormLimit", sph_norm_limit, 1., cLine, fLine, nLines, &
               "Limit for SPH normalisation: ", "(a,f10.4,a)", 0.5, ok, .false.)

          call getInteger("kerneltype", kerneltype, cLine, fLine, nLines, &
               "Kernel type (0 is exponential/1 is spline): ","(a,i1,a)",0, ok, .false.)

          call getLogical("useHull", useHull, cLine, fLine, nLines, &
            "Use hull particle method: ","(a,1l,a)", .false., ok, .false.)

          call getString("inputFileFormat", inputFileFormat, cLine, fLine, nLines, &
               "Input file format: ","(a,a,1x,a)","binary", ok, .false.)

! Conversion from total density to HI density when chemistry is included
! and grid needs to be loaded with HI density
          call getLogical("convertrhotohi", convertRhoToHI, cLine, fLine, nLines, &
               "Convert density to HI:", "(a,1l,1x,a)", .false., ok, .false.)
          call getInteger("ih2frac", ih2frac, cLine, fLine, nLines, &
               "Column containing H2 fraction ","(a,i2,a)", 11, ok, .false.)
          call getLogical("sphwithchem", sphwithchem, cLine, fLine, nLines, &
               "SPH has chemistry information:", "(a,1l,1x,a)", .false., ok, .false.)

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

    case("shakara")

       oneKappa = .true.
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

       call getReal("teff1", teff, 1., cLine, fLine, nLines, &
            "Source effective temperature: ","(a,f7.1,a)", 12., ok, .true.)


       call getLogical("setsubradius", setSubRadius, cLine, fLine, nLines, &
            "Set sublimation radius empirically (Whitney et al 2004): ","(a,1l,1x,a)", .false., ok, .false.)

       if (.not.setSubRadius) then
          rSublimation = rInner
          write(message, '(a,f7.1,a)') "Dust sublimation radius is ",rSublimation/rcore, " stellar radii"
          call writeInfo(message,TRIVIAL)
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
            "Core mass (solar masses): ","(a,f6.4,a)", 0.5, ok, .true.)

       call getReal("mdisc", mDisc, real(msol), cLine, fLine, nLines, &
            "Disc mass (solar masses): ","(a,f6.4,a)", 1.e-4, ok, .true.)

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

       if (solveVerticalHydro) then
          call getInteger("nhydro", nhydro,  cline, fLine, nLines, &
               "Max number of hydro iterations : ","(a,i4,a)", 5, ok, .true.)
       endif


       call getReal("heightsplitfac", heightSplitFac, 1., cLine, fLine, nLines, &
            "Splitting factor for scale height (local scale heights): ","(a,f5.2,a)", 0.2, ok, .false.)


       call getInteger("ndusttype", nDustType, cLine, fLine, nLines,"Number of different dust types: ","(a,i12,a)",1,ok,.false.)

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

    call getLogical("dorefine", dorefine, cLine, fLine, nLines, &
         "Adaptively refine AMR grid?: ","(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("dounrefine", dounrefine, cLine, fLine, nLines, &
         "Adaptively unrefine AMR grid?: ","(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("dophotorefine", dophotorefine, cLine, fLine, nLines, &
         "Adaptively refine AMR grid in photoloop?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("refineonmass", refineonmass, cLine, fLine, nLines, &
         "Refine grid using mass in cell?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("refineonjeans", refineonJeans, cLine, fLine, nLines, &
         "Refine grid using local jeans mass?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getUnitDouble("masstol", masstol, "mass" , cLine, fLine, nLines, &
         "Cell mass tolerance: ","(a,es9.3,1x,a)", 1.d-5*mSol, ok, .false.)

    call getLogical("refineontemperature", refineontemperature, cLine, fLine, nLines, &
         "Refine grid using temperature gradient?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("refineonrhoe", refineonrhoe, cLine, fLine, nLines, &
         "Refine grid using rhoe?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("refineonspeed", refineonspeed, cLine, fLine, nLines, &
         "Refine grid using speed?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("refineonionization", refineonionization, cLine, fLine, nLines, &
         "Refine grid using ionization gradient?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("captureshocks", captureShocks, cLine, fLine, nLines, &
         "Captue shocks?: ","(a,1l,1x,a)", .false., ok, .false.)

    call getDouble("amrtolerance", amrtolerance, 1.d0 , cLine, fLine, nLines, &
         "Maximum gradient allowed before AMR grid refines: ","(a,es9.3,1x,a)", 5.d-2, ok, .false.) 

    call getDouble("unreftol", amrUnrefinetolerance, 1.d0 , cLine, fLine, nLines, &
         "Minimum gradient allowed before AMR grid unrefines: ","(a,es9.3,1x,a)", 5.d-2, ok, .false.) 

    if(refineontemperature) then
       call getDouble("amrtemperaturetol", amrTemperatureTol, 1.d0 , cLine, fLine, nLines, &
            "Maximum temperature gradient allowed before AMR grid refines: ","(a,es9.3,1x,a)", 1.d-1, ok, .false.)
    end if

    if(refineonspeed) then
       call getDouble("amrspeedtol", amrSpeedTol, 1.d0 , cLine, fLine, nLines, &
            "Maximum speed gradient allowed before AMR grid refines: ","(a,es9.3,1x,a)", 1.d-1, ok, .false.)
    end if
    
    if(refineonionization) then
       call getDouble("amrionfractol", amrIonFracTol, 1.d0 , cLine, fLine, nLines, &
            "Maximum ionization fraction gradient allowed before AMR grid refines: ","(a,es9.3,1x,a)", 1.d-1, ok, .false.)
    end if

    if(refineonrhoe) then
       call getDouble("amrrhoetol", amrRhoeTol, 1.d0 , cLine, fLine, nLines, &
            "Maximum rhoe gradient allowed before AMR grid refines: ","(a,es9.3,1x,a)", 1.d-1, ok, .false.)
    end if
  
    call getReal("amrgridsize", amrGridSize, 1., cLine, fLine, nLines, &
         "Size of adaptive mesh grid: ","(a,1pe8.1,1x,a)", 1000., ok, .true.) 
    call getUnitDouble("amrgridcentrex", amrGridCentreX, "distance" , cLine, fLine, nLines, &
         "Grid centre X-coordinate: ","(a,es9.3,1x,a)", 0.0d0, ok, .false.) 
    call getUnitDouble("amrgridcentrey", amrGridCentreY, "distance" , cLine, fLine, nLines, &
         "Grid centre Y-coordinate: ","(a,es9.3,1x,a)", 0.0d0, ok, .false.) 
    call getUnitDouble("amrgridcentrez", amrGridCentreZ, "distance", cLine, fLine, nLines, &
         "Grid centre Z-coordinate: ","(a,es9.3,1x,a)", 0.0d0, ok, .false.) 

    if (amr2d.and.(.not.checkPresent("amrgridcentrex", cline, nlines))) &
         amrGridCentrex = amrGridSize/2.d0

    if (amr3d.and.cylindrical.and.(.not.checkPresent("amrgridcentrex", cline, nlines))) &
         amrGridCentrex = amrGridSize/2.d0

    if (amr1d.and.(.not.checkPresent("amrgridcentrex", cline, nlines))) &
         amrGridCentrex = amrGridSize/2.d0


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

       call writeBanner("Dust parameters","#",TRIVIAL)


       call writeBanner("NOTE THAT DUST TO GAS RATIO INFORMATION SHOULD BE PLACED IN GRAINFRAC. ","!",TRIVIAL)
       call writeBanner("DUSTTOGAS IS NOW USED TO DENOTE THE TOTAL DUST MASS IN THE SYSTEM.","!",TRIVIAL)

       call getReal("dusttogas", dusttoGas, 1., cLine, fLine, nLines, &
            "Dust to gas ratio: ","(a,f5.3,a)",0.01,ok,.false.)

       call getInteger("ndusttype", nDustType, cLine, fLine, nLines,"Number of different dust types: ","(a,i12,a)",1,ok,.false.)
       if (nDustType .gt. maxDustTypes) then
          if (writeoutput) write (*,*) "Max dust types exceeded: ", maxDustTypes
          stop
       end if
       
       call getLogical("readdust", readDustFromFile, cLine, fLine, nLines, &
         "Read dust opacities and phase matrices from file: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("writedust", writeDustToFile, cLine, fLine, nLines, &
         "Write dust opacities and phase matrices to file: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("dustfile", dustfile, cLine, fline, nLines, &
            "Get dust properties from file: ","(a,1l,1x,a)", .false., ok, .false.)
       if (dustfile) then
          call getString("kappafile", dustFilename(1), cline, fLine, nLines, &
               "Dust properties filename: ","(a,a,1x,a)","none", ok, .true.)
          nDustType = 1
          grainfrac = 1.
          isoTropicScattering = .true.
          

       else

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
             
             if (.not. readDustFromFile) &
                  call getReal(aminLabel, aMin(i), 1., cLine, fLine, nLines, &
                  "Min grain size (microns): ","(a,f8.5,1x,a)", 0.005, ok,  .true.)
             
             if (.not. readDustFromFile) &
                  call getReal(amaxLabel, aMax(i), 1., cLine, fLine, nLines, &
                  "Max grain size (microns): ","(a,f10.5,1x,a)", 0.25, ok, .true.)
             
             if (.not. readDustFromFile) &
                  call getReal(qDistLabel, qdist(i), 1., cLine, fLine, nLines, &
                  "Grain power law: ","(a,f4.1,1x,a)", 3.5, ok, .true. )
             
             if (.not. readDustFromFile) &
                  call getReal(a0Label, a0(i), 1., cLine, fLine, nLines, &
                  "Scale length of grain size (microns): ","(a,f8.5,1x,a)", 1.0e20, ok, .false.)
             
             
             if (.not. readDustFromFile) &
          call getReal(pdistLabel, pdist(i), 1., cLine, fLine, nLines, &
          "Exponent for exponential cut off: ","(a,f4.1,1x,a)", 1.0, ok, .false. )




             if (writeoutput) write(*,*)
          enddo
          if (.not. readDustFromFile) &
               call getLogical("iso_scatter", isotropicScattering, cLine, fLine, nLines, &
               "Isotropic scattering: ","(a,1l,1x,a)", .false., ok, .false.)
       endif
  end subroutine readDustPhysicsParameters

  subroutine readAtomicPhysicsParameters(cLine, fLine, nLines)
    use atom_mod, only: setVoigtParams
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    integer :: i 
    character(len=20) :: keyword
    logical :: ok
    real :: C_rad, C_vdw, C_stark

    call writeBanner("Atomic physics data","#",TRIVIAL)



    call getLogical("cmf", cmf, cLine, fLine, nLines, &
         "Perform co-moving frame calculation ","(a,1l,1x,a)", .false., ok, .true.)

    if (cmf) then
       call getInteger("natom", nAtom, cLine, fLine, nLines, &
            "Number of model atoms to solve for: ","(a,i1,a)",1,ok,.true.)
       do i = 1, nAtom
          write(keyword, '(a,i1)') "atom",i
          call getString(keyword, atomFileName(i), cLine, fLine, nLines, &
               "Use atom filename: ","(a,a,1x,a)","none", ok, .true.)
       enddo
    endif


    call getLogical("thickcont", opticallyThickContinuum, cLine, fLine, nLines, &
         "Continuum is optically thick: ","(a,1l,a)", .false., ok, .false.)

    call getDouble("xabundance", Xabundance, 1.d0, cLine, fLine, nLines, &
            "Hydrogen abundance (by number): ","(a,f7.3,a)",1.d0, ok, .true.)

    call getDouble("yabundance", Yabundance, 1.d0, cLine, fLine, nLines, &
            "Helium abundance (by number): ","(a,f7.3,a)",1.d0, ok, .true.)

    call getReal("lamline", lamLine, 1.,cLine, fLine, nLines, &
               "Line emission wavelength: ","(a,f8.1,1x,a)", 850., ok, .true.)          

    call getReal("vturb", vturb, real(kmstoc), cLine, fLine, nLines, &
               "Turbulent velocity (km/s):","(a,f6.1,1x,a)", 50., ok, .true.)
    !
    ! Voigt profile prameters
    !
    call getReal("C_rad", C_rad, 1., cLine, fLine, nLines, &
         "Damping constant (radiation)     in [A]: ","(a,1PE10.3,1x,a)", 0.0, ok, .false.)
    call getReal("C_vdw", C_vdw, 1., cLine, fLine, nLines, &
         "Damping constant (van der Waals) in [A]: ","(a,1PE10.3,1x,a)", 0.0, ok, .false.)
    call getReal("C_stark", C_stark, 1., cLine, fLine, nLines, &
         "Damping constant (Stark)         in [A]: ","(a,1PE10.3,1x,a)", 0.0, ok, .false.)
    call setVoigtParams(C_rad, C_vdw, C_stark)

  end subroutine readAtomicPhysicsParameters

  subroutine readSourceParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:), message
    logical :: fLine(:)
    integer :: nLines
    integer :: i 
    character(len=20) :: keyword
    logical :: ok





    call writeBanner("Photon source data","#",TRIVIAL)



    call getInteger("nsource", inputNSource, cLine, fLine, nLines, &
         "Number of sources: ","(a,i2,a)",1,ok,.true.)
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
          call getDouble(keyword, sourceMass(i), 1.d0, cLine, fLine, nLines, &
               "Source mass (solar masses) : ","(a,f7.2,a)",1.d0, ok, .true.)


       else

          if (stellarSource(i)) then
             write(keyword, '(a,i1)') "radius",i
             call getDouble(keyword, sourceRadius(i), rsol/1.d10, cLine, fLine, nLines, &
                  "Source radius (solar radii) : ","(a,f7.2,a)",1.d0, ok, .true.)

             write(keyword, '(a,i1)') "teff",i
             call getUnitDouble(keyword, sourceTeff(i), "temperature", cLine, fLine, nLines, &
                  "Source temperature (K) : ","(a,f8.0,a)",1.d0, ok, .true.)

             write(keyword, '(a,i1)') "mass",i
             call getDouble(keyword, sourceMass(i), mSol, cLine, fLine, nLines, &
                  "Source mass (solar masses) : ","(a,f7.2,a)",1.d0, ok, .true.)

             write(keyword, '(a,i1)') "mdot",i
             call getDouble(keyword, sourceMdot(i), msol/(365.25d0 * 24.d0 * 3600.d0), cLine, fLine, nLines, &
                  "Source mass-loss rate (solar masses/year) : ","(a,f9.6,a)",0.d0, ok, .false.)

             write(keyword, '(a,i1)') "contflux",i
             call getString(keyword, inputcontFluxFile(i), cLine, fLine, nLines, &
                  "Continuum flux file: ","(a,a,a)","none", ok, .true.)

             write(keyword, '(a,i1)') "sourcepos",i
             call getVector(keyword, sourcePos(i), 1.d0, cLine, fLine, nLines, &
                  "Source position (10^10 cm): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)

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
               "Number of fixed rays for stage 1: ","(a,i4,a)", 1024, ok, .false.)
       call getLogical("dongstep", dongstep, cLine, fLine, nLines, &
               "Use Ng Acceleration: ","(a,1l,a)", .false., ok, .false.)
       call getLogical("quasi", quasi, cLine, fLine, nLines, &
               "Use Quasirandom numbers: ","(a,1l,a)", .false., ok, .false.)
       call getLogical("modelwashydro", modelwashydro, cLine, fLine, nLines, &
               "Grid is based on hydro or radiation hydro calculation: ","(a,1l,a)", .false., ok, .false.)
!thaw -        
       if(modelWasHydro) then
          call getLogical("zeroghosts", zeroghosts, cLine, fLine, nLines, &
               "Zero the contriubtion from ghost cells to molecular line calc: ","(a,1l,a)", .false., ok, .false.)
       else
          zeroghosts=.false.
       end if
  
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
       call getLogical("removehotmolecular", removeHotMolecular, cLine, fLine, nLines, &
            "Remove molecular material above 100K: ", "(a,1l,1x,a)", .false., ok, .true.)
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


       call getLogical("doCOchemistry", doCOchemistry, cLine, fLine, nLines, &
            "Use drop profile to model CO depletion: ", "(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("useDust", useDust, cLine, fLine, nLines, &
            "Calculate continuum emission from dust:", "(a,1l,1x,a)", .false., ok, .false.)

       if(doCOchemistry) then

          call getReal("fracCOdepletion", x_D, 1., cLine, fLine, nLines, &
               "Fraction of CO depletion", "(a,f5.3,a)", 0.1, ok, .true.)
          
       endif

       


  end subroutine readMolecularPhysicsParameters

#ifdef PHOTOION
  subroutine readPhotoionPhysicsParameters(cLine, fLine, nLines)
    use ion_mod, only: setAbundances
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines, thisStackLim, numMPIthreads, i
    character(len=4)  :: iChar
    character(len=20) :: keyword
    logical :: ok
    real :: h_abund, he_abund, c_abund, n_abund, o_abund, ne_abund, s_abund  

    call getLogical("quickthermal", quickThermal, cLine, fLine, nLines, &
         "Compute photoionization equilibrium: ","(a,1l,a)", .false., ok, .false.)

    call getLogical("usemetals", usemetals, cLine, fLine, nLines, &
         "Use metals in photoionization calculation: ","(a,1l,a)", .true., ok, .false.)

    call getLogical("vardustsub", variableDustSublimation, cLine, fLine, nLines, &
         "Variable dust sublimation temperature: ", "(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("nodiffuse", noDiffuseField, cLine, fLine, nLines, &
         "Ignore diffuse radiation field: ","(a,1l,a)", .false., ok, .false.)

    call getLogical("radpresstest", radPressureTest, cLine, fLine, nLines, &
         "Radiation-pressure test (on-the-spot absorption): ","(a,1l,a)", .false., ok, .false.)


    if(nodiffusefield) then
       call getLogical("monochromatic", monochromatic, cLine, fLine, nLines, &
            "Use a monochromatic source:", "(a,1l,1x,a)", .false., ok, .false.)
    end if


    call getLogical("biasToLyman", biasToLyman, cLine, fLine, nLines, &
         "Variance reduction, higher sampling of Lyman photons: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("binPhotons", binPhotons, cLine, fLine, nLines, &
         "Bin and dump photons as a function of wavelength: ","(a,1l,1x,a)", .false., ok, .false.)

    call getDouble("biasMagnitude", biasMagnitude, 1.d0, cLine, fLine, nLines, &
            "Variance reduction, extent of bias: ", "(a,es9.3,1x,a)", 100.d0, ok, .false.)

    call getLogical("hOnly", hOnly, cLine, fline, nLines, &
         "Hydrogen-only calculation: ", "(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("startfromneutral", startFromNeutral, cLine, fline, nLines, &
         "Start photoionization loop from neurtral: ", "(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("stellarwinds", stellarWinds, cLine, fline, nLines, &
         "Include stellar wind outflow feedback: ", "(a,1l,1x,a)", .false., ok, .false.)

    call getReal("tminglobal", TMinGlobal, 1., cLine, fLine, nLines, &
         "Minimum Temperature (K): ","(a,f4.1,1x,a)", 10., ok, .false.)

    if(splitovermpi) then
       call getLogical("optimizeStack", optimizeStack, cLine, fLine, nLines, &
            "Perform a bundle size optimization calculation: ","(a,1l,1x,a)", .false., ok, .false.)

       call getInteger("stacklimit", stacklimit, cLine, fLine, nLines, &
            "Maximum number of photons in a stack","(a,i8,a)", 200, ok, .false.)
       
       call getInteger("dstack", dstack, cLine, fLine, nLines, &
            "Optimization tweak level","(a,i8,a)", 100, ok, .false.)

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
            "Number of photon stacks allowed in the buffer ","(a,i8,a)", 40000, ok, .false.)
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
         "Turbulent velocity (km/s):","(a,f6.1,1x,a)", 50., ok, .true.)


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
         "(a,f6.1,a)",0.0,ok,.false.)

    call getLogical("convergeonundersampled", convergeOnUndersampled, cLine, fLine, nLines, &
         "Prevent convergence when too many cells undersampled: ", "(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("vardustsub", variableDustSublimation, cLine, fLine, nLines, &
         "Variable dust sublimation temperature: ", "(a,1l,1x,a)", .false., ok, .false.)

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
            "Maximum tau for smoothing: ","(a,f10.1,1x,a)", 0.5, ok, .false.)
       call getReal("taumin", tauSmoothMin, 0.1, cLine, fLine, nLines, &
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

    call getInteger("nlucy", nLucy, cLine, fLine, nLines,"Number of photons per lucy iteration: ","(a,i12,a)",0,ok,.false.)

    call getInteger("nmonte", inputnMonte, cLine, fLine, nLines, &
         "Number of photons in image","(a,i8,a)", 0, ok, .false.)


    call getInteger("mincrossings", minCrossings, cLine, fLine, nLines, &
         "Minimum crossings required for cell to be sampled: ","(a,i12,a)",100,ok,.false.)

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
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call getReal("cfl", cflNumber, 1., cLine, fLine, nLines, &
         "Courant number:","(a,f4.1,1x,a)", 0.3, ok, .false.)

    call getLogical("modelwashydro", modelwashydro, cLine, fLine, nLines, &
         "Grid is based on hydro or radiation hydro calculation: ","(a,1l,a)", .false., ok, .false.)

    call getDouble("etaviscosity", etaViscosity, 1.d0, cLine, fLine, nLines, &
         "Viscosity eta parameter:  ","(a,e12.3,1x,a)", 3.d0, ok, .false.)

    call getLogical("tensorviscosity", useTensorViscosity, cLine, fLine, nLines, &
         "Use tensor form for artificial viscosity: ","(a,1l,1x,a)", .false., ok, .false.)

    call getDouble("griddistancescale", gridDistanceScale, 1.d0, cLine, fLine, nLines, &
         "Distance grid scale: ","(a,e12.3,1x,a)", 1.d10, ok, .true.)

    call getUnitDouble("tstart", tStart, "time", cLine, fLine, nLines, &
         "Start time for hydrodynamical calculatione: ","(a,e12.3,1x,a)", 0.d0, ok, .false.)

    call getUnitDouble("tend", tEnd, "time", cLine, fLine, nLines, &
         "End time for calculation: ","(a,e12.3,1x,a)", 1.d10, ok, .true.)

    call getUnitDouble("tdump", tDump, "time", cLine, fLine, nLines, &
         "Time between dump files: ","(a,e12.3,1x,a)", 0.d0, ok, .true.)

    call getLogical("rhieChow", rhieChow, cLine, fLine, nLines, &
         "Use Rhie-Chow interpolation: ","(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("radpressure", radiationPressure, cLine, fLine, nLines, &
         "Include radiation pressure terms: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("gasgravity", doGasGravity, cLine, fLine, nLines, &
         "Include gas gravity: ","(a,1l,1x,a)", .true., ok, .false.)

    call getLogical("fluxinterp", fluxinterp, cLine, fLine, nLines, &
         "Interpolate flux at fine coarse boundaries: ","(a,1l,1x,a)", .false., ok, .false.)

    call getString("limitertype", limiterType, cLine, fLine, nLines, &
         "Flux limiter to use: ","(a,a,1x,a)","superbee", ok, .false.)

    call getLogical("readturb", readTurb, cLine, fLine, nLines, &
         "read in a turbulent velocity field from file: ","(a,1l,1x,a)", .false., ok, .false.)

    if(readTurb) then
       call getString("turbvelfilex", turbVelFilex, cLine, fLine, nLines, &
       "Turbulent velocity field file (x): ","(a,a,1x,a)","cube_v1.dat", ok, .false.)
       
       call getString("turbvelfiley", turbVelFiley, cLine, fLine, nLines, &
       "Turbulent velocity field file (y): ","(a,a,1x,a)","cube_v2.dat", ok, .false.)
       
       call getString("turbvelfilez", turbVelFilez, cLine, fLine, nLines, &
       "Turbulent velocity field file (z): ","(a,a,1x,a)","cube_v3.dat", ok, .false.)

       call getInteger("nturblines", nturblines, cLine, fLine, nLines, &
       "Number of lines in turbulent velocity field file: ","(a,i2,a)", 128, ok, .false.)
       
       
    end if


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

  end subroutine readHydrodynamicsParameters


  subroutine readPhotometryParameters(cLine, fLine, nLines)
    character(len=80) :: cLine(:)
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
    use datacube_mod, only: cubePositionAngle, npixels, datacubeFilename
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok
    TYPE(VECTOR) :: gridCentre
    character(len=100) :: message

    call getReal("inclination", thisinclination, real(degtorad), cLine, fLine, nLines, &
         "Inclination angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)
    call getReal("positionangle", cubePositionAngle, real(degtorad), cLine, fLine, nLines, &
         "Position angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)
    call getString("datacubefile", datacubeFilename, cLine, fLine, nLines, &
         "Output datacube  filename: ","(a,a,1x,a)","none", ok, .true.)
    call getReal("imageside", imageside, 1., cLine, fLine, nLines, &
         "Image size (x10^10cm):","(a,es9.3,1x,a)", 5e7, ok, .true.)
    call getReal("distance", gridDistance, real(pcToCm), cLine, fLine, nLines, &
         "Grid distance (pc): ","(a,f4.1,1x,a)", 100., ok, .false.)
    call getInteger("npixels", npixels, cLine, fLine, nLines, &
         "Number of pixels per row: ","(a,1x,i4,a)", 50, ok, .true.)
    call getInteger("ncubes", ncubes, cLine, fLine, nLines, &
         "Number of data cubes ","(a,i4,a)", 1, ok, .false.)
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

    if (atomicPhysics) then
          call getReal("lamline", lamLine, 1.,cLine, fLine, nLines, &
               "Line emission wavelength: ","(a,f8.1,1x,a)", 850., ok, .true.)          
          call getReal("vturb", vturb, real(kmstoc), cLine, fLine, nLines, &
               "Turbulent velocity (km/s):","(a,f6.1,1x,a)", 50., ok, .true.)
    endif

    call getLogical("internalView", internalView, cLine, fLine, nLines, &
         "View as our Galaxy:", "(a,1l,1x,a)", .false., ok, .false.)

    ! Read parameters used by galactic plane survey 
    if ( internalView ) then 

       call getLogical("obsVelFromGrid", obsVelFromGrid, cLine, fLine, nLines, &
            "Set observer velocity from grid:", "(a,1l,1x,a)", .false., ok, .true.)

       call getDouble("intPosX", intPosX,  1.0_db, cLine, fLine, nLines, "Observer x position (x10^10cm)", &
            "(a,e10.4,1x,a)", 0.d0, ok, .false.)
       call getDouble("intPosY", intPosY,  1.0_db, cLine, fLine, nLines, "Observer y position (x10^10cm)", &
            "(a,e10.4,1x,a)", 2.2e12_db, ok, .false.)
       call getDouble("intPosZ", intPosZ, 1.0_db,  cLine, fLine, nLines, "Observer z position (x10^10cm)", &
            "(a,e10.4,1x,a)", 0.d0, ok, .false.)
          
       call getDouble("intDeltaVx", intDeltaVx, 1.0_db,  cLine, fLine, nLines, "Observer x velocity boost (km/s)", &
            "(a,f8.2,1x,a)", 0.d0, ok, .false.)
       call getDouble("intDeltaVy", intDeltaVy, 1.0_db,  cLine, fLine, nLines, "Observer y velocity boost (km/s)", &
            "(a,f8.2,1x,a)", 0.d0, ok, .false.)
       call getDouble("intDeltaVz", intDeltaVz, 1.0_db,  cLine, fLine, nLines, "Observer z velocity boost (km/s)", &
               "(a,f8.2,1x,a)", 0.d0, ok, .false.)

       call getLogical("thermalLineWidth", thermalLineWidth, cLine, fLine, nLines, &
            "Thermal line width:", "(a,1l,1x,a)", .true., ok, .false.)
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
       else
          ! Option for agressive splitting
          call getLogical("sphonepercell", SphOnePerCell, cLine, fLine, nLines, &
            "Split to one particle per cell:", "(a,1l,1x,a)", .false., ok, .false.)
       end if

    end if

    if ( h21cm ) then 

       call getReal("tminglobal", TMinGlobal, 1., cLine, fLine, nLines, &
            "Minimum Temperature (K): ","(a,f4.1,1x,a)", 10., ok, .false.)
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
       itrans = 1 ! always 1 for the 21cm line

       if ( .not. internalView ) then 
          ! Far field h21cm case
          call getDouble("galaxyInclination", galaxyInclination, 1.0_db, cLine, fLine, nLines, &
               "Galaxy Inclination:", "(a,f4.1,1x,a)", 50.d0, ok, .false.)
          call getDouble("galaxyPositionAngle", galaxyPositionAngle, 1.0_db, cLine, fLine, nLines, &
               "Galaxy position angle:", "(a,f4.1,1x,a)", 20.d0, ok, .false.)
! In molecular_mod the observer is looking along the y-axis
          rotateViewAboutX = 90.0 - galaxyInclination
          rotateViewAboutY = galaxyPositionAngle - 90.0
          rotateViewAboutZ = 0.0
       end if
       
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
       if (internalView) then
          call getDouble("centrevecx", centrevecx, 1.d0, cLine, fLine, nLines, &
               "Image Centre Coordinate (lon): ","(a,f7.2,1x,a)", 0.d0, ok, .true.)
          call getDouble("centrevecy", centrevecy, 1.d0, cLine, fLine, nLines, &
               "Image Centre Coordinate (lat): ","(a,f7.2,1x,a)", 0.d0, ok, .true.)
       else
          call getDouble("centrevecx", centrevecx, 1.d0, cLine, fLine, nLines, &
               "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .true.)
          call getDouble("centrevecy", centrevecy, 1.d0, cLine, fLine, nLines, &
               "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .true.)
          call getDouble("centrevecz", centrevecz, 1.d0, cLine, fLine, nLines, &
               "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .true.)
       end if
       call getLogical("wanttau", wanttau, cLine, fLine, nLines, &
            "Write Tau information to datacube: ","(a,1l,1x,a)", .false., ok, .false.)
    endif
       
  end subroutine readDataCubeParameters
#endif

  subroutine readImageParameters(cLine, fLine, nLines)
    use image_utils_mod

    character(len=80) :: cLine(:), message
    logical :: fLine(:)
    integer :: nLines
    logical :: ok
    integer :: i
    character(len=20) :: keyword
    character(len=10) :: axisUnits
    character(len=10) :: fluxUnits
    character(len=80) :: outputImageType, imageFilename
    character(len=4)  :: iChar
    integer :: thisnpixels, npixels, nimage
    real :: lambdaImage, thisLambdaImage
    real :: imagesize, thisImageSize, wholeGrid
    real :: inclination, positionAngle, thisPA, thisInc
    real :: offsetX, offsetY, thisOffsetX, thisOffsetY
    real :: aspectRatio, thisAspectRatio

    call getBigInteger("nphotons", nphotons, cLine, fLine, nLines, &
         "Number of photons in image: ","(a,i8,a)", 10000, ok, .true.)

    call getReal("distance", gridDistance, 1., cLine, fLine, nLines, &
         "Grid distance (pc): ","(a,f4.1,1x,a)", 100., ok, .false.)

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
       wholeGrid = (amrGridSize*1.0e10) / real(autocm)
    case ("pc","PC")
       wholeGrid = (amrGridSize*1.0e10) / real(pctocm)
    case ("cm")
       wholeGrid = amrGridSize*1.0e10
    case ("arcsec")
       wholeGrid = (amrGridSize*1.0e10) / real(pctocm) ! pc
       wholeGrid = wholeGrid / gridDistance            ! radians
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
         "Inclination angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)
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

    call getReal("lambdaimage", lambdaImage, 1., cLine, fLine, nLines, &
         "Wavelength for monochromatic image (A):","(a,f8.1,1x,a)", 6562.8, ok, .false.)

    if (nimage == 1) then

       call getString("imagefile", imageFilename, cLine, fLine, nLines, &
            "Output image  filename: ","(a,a,1x,a)","none", ok, .true.)

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
              
       call setImageParams(1, lambdaImage, outputimageType, imageFilename, npixels, axisUnits, fluxUnits, &
            imageSize, aspectRatio, inclination, positionAngle, offsetx, offsety, gridDistance)
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
          call getReal(keyword, thisLambdaImage, 1., cLine, fLine, nLines, &
               "Wavelength for monochromatic image (A):","(a,f8.1,1x,a)", lambdaImage, ok, .false.)

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

          ! Inclination and position angle
          write(keyword,'(a)') "inclination"//iChar
          call getReal(keyword, thisInc, real(degtorad), cLine, fLine, nLines, &
               "Inclination of image: ","(a,f4.1,1x,a)",inclination*real(radtodeg), ok, .false.)
          write(keyword,'(a)') "positionangle"//iChar
          call getReal(keyword, thisPA, real(degtorad), cLine, fLine, nLines, &
               "Position angle (deg): ","(a,f4.1,1x,a)", positionAngle*real(radtodeg), ok, .false.)

          ! Position of image centre
          write(keyword,'(a)') "imagecentrex"//iChar
          call getReal(keyword, thisOffsetx, 1.0, cLine, fLine, nLines, &
               "Image centre x position (10^10 cm): ", "(a,e10.1,1x,a)", offsetx, ok, .false.)
          write(keyword,'(a)') "imagecentrey"//iChar
          call getReal(keyword, thisOffsety, 1.0, cLine, fLine, nLines, &
               "Image centre y position (10^10 cm): ", "(a,e10.1,1x,a)", offsety, ok, .false.)

          call setImageParams(i, thisLambdaImage, outputimageType,imageFilename, thisnpixels, axisUnits, fluxUnits, &
               thisImageSize, thisaspectRatio, thisInc, thisPA, thisOffsetx, thisOffsety, gridDistance)
       enddo

    end if
   
  end subroutine readImageParameters

  subroutine readFitsParameters(cLine, fLine, nLines)
    use fits_utils_mod, only: fitsBitpix
    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok

    call getInteger("fitsbitpix", fitsBitpix, cLine, fLine, nLines, &
         "FITS file BITPIX ","(a,i2,a)", -32, ok, .false.)

  end subroutine readFitsParameters

  subroutine readSpectrumParameters(cLine, fLine, nLines)
    use sed_mod, only:  setSedParameters,  SEDlamMin, SEDlamMax, SEDwavLin, SEDnumLam

    character(len=80) :: cLine(:)
    logical :: fLine(:)
    integer :: nLines
    logical :: ok, isRange
    character(len=80) :: message
    logical :: sed, jansky, SIsed
    integer :: nInclination
    real    :: firstInclination
    real    :: lastInclination=80.0
    real, allocatable :: inclinations(:)
    character(len=80) :: outFile 


    call getBigInteger("nphotons", nPhotons, cLine, fLine, nLines, &
         "Number of photons in SED: ", "(a,i15,1x,a)", 100000, ok, .false.)

    call getInteger("ninc", nInclination, cLine, fLine, nLines, &
         "Number of inclination angles: ", "(a,i3,1x,a)", 1, ok, .false.)

    if (checkPresent("inclinations", cline, nlines)) then
       call parameterDefinesRange("inclinations", firstInclination, lastInclination, isRange, cLine, fLine, nLines)
       if (isRange) then 
          write(message,'(a,f6.2,a,f6.2,a)') &
               "SED inclination range is from ", firstInclination, " to ", lastInclination, " degrees"
          call writeInfo(message,TRIVIAL)
          firstInclination = firstInclination * real(degToRad)
          lastInclination = lastInclination * real(degToRad)
       else
          allocate(inclinations(nInclination))
          call getRealArray("inclinations", inclinations, 1.0, cLine, fLine, nLines, &
               "Inclinations (deg): ",90., ok, .false.)
          inclinations(:) = inclinations(:) * real(degToRad)
       end if
    else
       call getReal("firstinc", firstInclination, 1.0, cLine, fLine, nLines, &
            "First inclination angle (deg): ","(a,f4.1,1x,a)", 10., ok, .true.)
       firstInclination = firstInclination * real(degToRad)
       if (nInclination > 1) &
            call getReal("lastinc", lastInclination, 1.0, cLine, fLine, nLines, &
            "Last inclination angle (deg): ","(a,f4.1,1x,a)", 80., ok, .true.)
       lastInclination = lastInclination * real(degToRad)
    end if

    call getReal("inclination", thisinclination, real(degtorad), cLine, fLine, nLines, &
         "Inclination angle (deg): ","(a,f4.1,1x,a)", 0., ok, .false.)


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
         "Grid distance (pc): ","(a,f4.1,1x,a)", 100., ok, .false.)

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

    if (allocated(inclinations)) then 
       call setSedParameters(outFile,jansky,SIsed,sed,incList=inclinations)
       deallocate(inclinations)
    else
       call setSedParameters(outFile,jansky,SIsed,sed,nInclination=nInclination,&
            firstInc=firstInclination,LastInc=LastInclination, cosSpacing=.true.)
    end if

  end subroutine readSpectrumParameters


!-----------------------------------------------------------------------------------------


subroutine findReal(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 real :: value
 character(len=80) :: cLine(:)
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
 character(len=80) :: cLine(:)
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

 subroutine findUnitDouble(name, value, unit, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 real(double) :: value
 character(len=80) :: cLine(:)
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
        read(cLine(i)(j+1:80),*,iostat=readstat) value, unit
        if(readstat /= 0) then
           read(cLine(i)(j+1:80),*,iostat=readstat) value
           if(readstat /= 0) then
              write(*,*) "parameters.dat entry with no value given!!" //name
              stop
           else
              unit = "default"
           end if
        end if
        fLine(i) = .true.
     else
     call writeFatal("Keyword "//name(1:j)//" appears more than once in the input deck")
     stop
  endif
endif
end do
 end subroutine findUnitDouble

 subroutine findVECTOR(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 type(VECTOR) :: value
 character(len=80) :: cLine(:)
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
 character(len=80) :: cLine(:)
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
   character(len=80) :: cLine(:)
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
 character(len=80) :: cLine(:)
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
 character(len=80) :: cLine(:), tmp
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
     stop
  endif
endif
 end do
 value = trim(value)
 end subroutine findString

! Return the whole line containing a given keyword
subroutine findLine(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 character(len=80) :: value
 character(len=80) :: cLine(:)
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
     stop
  endif
endif
 end do
 value = trim(value)
 end subroutine findLine

subroutine findRealArray(name, value, cLine, fLine, nLines, ok)
 implicit none
 character(len=*) :: name
 real :: value(:)
 character(len=80) :: cLine(:)
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
  character(len=80) :: cLine(:)
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
  if (ok) then
     write(output,format) trim(message),ival,default
     call writeInfo(output, TRIVIAL)
  endif
 end subroutine getInteger

 subroutine getBigInteger(name, ival, cLine, fLine, nLines, message, format, idef, ok, &
                       musthave)
  character(len=*) :: name
  integer(kind=bigInt) :: ival
  character(len=80) :: cLine(:)
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
  if (ok) then
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
  character(len=80) :: cLine(:)
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
 if (ok) then
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


 subroutine getDouble(name, dval, unitConversion, cLine, fLine, nLines, message, &
   format, ddef, ok, musthave)
  character(len=*) :: name
  real(double) :: dval, unitConversion
  logical :: musthave
  character(len=80) :: cLine(:)
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
 if (ok) then
    write(output,format) trim(message),dval,default
    call writeInfo(output, TRIVIAL)
 endif
 dval = dval  * unitConversion
 end subroutine getDouble


 subroutine getUnitDouble(name, dval, unitType, cLine, fLine, nLines, message, &
   format, ddef, ok, musthave)
  character(len=*) :: name
  real(double) :: dval
  logical :: musthave
  character(len=80) :: cLine(:)
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
       stop
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
  character(len=80) :: cLine(:)
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
 if (ok) then
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
  character(len=80) :: cLine(:)
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
  character(len=80) :: cLine(:)
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
       STOP
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
  character(len=80) :: cLine(:)
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
      case DEFAULT
         print *, "Unrecognised boundary string:", boundaryString
         print *, "Halting."
         stop
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
       stop
    else
       close(53)
    endif
  end subroutine testFileExists

  subroutine parameterDefinesRange(name, firstVal, lastVal, isRange, cLine, fLine, nLines)

    character(len=*), intent(in)  :: name 
    real, intent(out)             :: firstVal
    real, intent(out)             :: lastVal
    logical, intent(out)          :: isRange
    character(len=80), intent(in) :: cLine(:)
    logical, intent(in)           :: fLine(:)
    integer, intent(in)           :: nLines

    character(len=80) :: thisString
    logical           :: ok 
    integer           :: i, j

    call findLine(name, thisString, cLine, fLine, nLines, ok)

! Look for .. to determine if this is a range of values 
! The parameters file is 80 column maximum
    isRange = .false. 
    do i=1, 76
       if (thisString(i:i+1) == "..") then 
          isRange = .true.
          thisString(i:i+1) = "  "
          exit
       end if
    end do

    if (isRange) then 
       j = len_trim(name)
       read(thisString(j+1:80),*) firstVal, lastVal
    else
       firstVal = -1.0e-33
       lastVal  =  -1.0e-33
       return
    end if

  end subroutine parameterDefinesRange

end module inputs_mod

