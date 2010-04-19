module inputs_mod

  use vector_mod
  use unix_mod
  use messages_mod
  use kind_mod
  use input_variables
  use constants_mod

  implicit none

contains
  
  subroutine  inputs()

    implicit none

    integer :: nLines
    integer :: errNo
    logical :: ok
    character(len=80) :: cLine(200) 
    logical :: done

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

    nBlobs = 0
    nLines = 0
    inputnMonte = 0

    call unixGetEnv("TORUS_JOB_DIR",absolutePath)
!   call get_environment_variable("TORUS_JOB_DIR",absolutePath)

    if (writeoutput) write(*,*) absolutePath
    paramFile = trim(absolutePath)//"parameters.dat"

    open(unit=32, file=paramfile, status='old', iostat=error)
    if (error /=0) then
       print *, 'Panic: parameter file open error, file:',trim(paramFile) ; stop
    end if

    do
       nLines = nLines + 1
       read(32,'(a80)',end=10) cLine(nLines)
       if (trim(cLine(nLines)(1:1)) == "%") nLines = nLines - 1   ! % is a comment
    end do
10  continue
    nLines = nLines - 1
    errno = 0

    call getInteger("verbosity", verbosityLevel, cLine, nLines, &
         "Verbosity level: ", "(a,i8,1x,a)", 3, ok, .false.)

! the grid setup. Either we read the grid in or set it up from scratch

    call writeBanner("Grid setup parameters","*",TRIVIAL)

    call getLogical("splitovermpi", splitOverMPI, cLine, nLines, &
         "Grid is domain decomposed over MPI: ","(a,1l,1x,a)", .false., ok, .true.)


    call getLogical("readgrid", readGrid, cLine, nLines, &
         "Read grid file: ","(a,1l,1x,a)", .false., ok, .true.)

    if (.not.readgrid) then
       call readGridInitParameters(cLine, nLines)
       call readGeometrySpecificParameters(cLine, nLines)
    else
       call getString("inputfile", gridInputFilename, cLine, nLines, &
                  "Grid input filename: ","(a,a,1x,a)","none", ok, .true.)
    endif


! the physical ingredients of the model

    call writeBanner("Physical ingredients of model","*",TRIVIAL)

    call getLogical("dustphysics", dustPhysics, cLine, nLines, &
         "Include dust physics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("atomicphysics", atomicPhysics, cLine, nLines, &
         "Include atomic physics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("molecularphysics", molecularPhysics, cLine, nLines, &
         "Include atomic physics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("photoionphysics", photoionPhysics, cLine, nLines, &
         "Include photoionization physics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    if (.not.(dustPhysics.or.atomicPhysics.or.molecularPhysics.or.photoionPhysics)) then
       call writeFatal("Must include one of: dustPhysics, atomicPhysics, molecularPhysics, photoionPhysics")
    endif


    if (dustPhysics) call readDustPhysicsParameters(cLine, nLines)
    if (atomicPhysics) call readAtomicPhysicsParameters(cLine, nLines)
    if (molecularPhysics) call readMolecularPhysicsParameters(cLine, nLines)
    if (photoionPhysics) call readPhotoionPhysicsParameters(cLine, nLines)

    if (checkPresent("nsource", cLine, nLines)) then
       call readSourceParameters(cLine, nLines)
    endif


! the type of calculation

    call writeBanner("Type of model calculation","*",TRIVIAL)

    call getLogical("stateq", statisticalEquilibrium, cLine, nLines, &
         "Perform a statistical equililbrium calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("radeq", radiativeEquilibrium, cLine, nLines, &
         "Perform a radiative equilibrium calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("photoioneq", photoionEquilibrium, cLine, nLines, &
         "Perform a photoionization calculation: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("hydrodynamics", hydrodynamics, cLine, nLines, &
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

    if (statisticalEquilibrium.and.molecularPhysics) call readMolecularLoopParameters(cLine, nLines)
    if (statisticalEquilibrium.and.atomicPhysics) call readAtomicLoopParameters(cLine, nLines)
    if (radiativeEquilibrium) call readRadiativeEquilibriumParameters(cLine, nLines)
    if (photoionEquilibrium) call readPhotoionEquilibriumParameters(cLine, nLines)
    if (hydrodynamics) call readHydrodynamicsParameters(cLine, nLines)

! now do we dump the output grid

    call writeBanner("Output file details","*",TRIVIAL)

    call getLogical("writegrid", writeGrid, cLine, nLines, &
         "Write grid file: ","(a,1l,1x,a)", .false., ok, .true.)

    if (writeGrid) then
       call getString("outputfile", gridOutputFilename, cLine, nLines, &
                  "Grid output filename: ","(a,a,1x,a)","none", ok, .true.)
    endif

! now we select the outputs

    call getLogical("datacube", calcDataCube, cLine, nLines, &
         "Calculate a data cube: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("image", calcImage, cLine, nLines, &
         "Calculate an image: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("spectrum", calcSpectrum, cLine, nLines, &
         "Calculate a spectrum: ","(a,1l,1x,a)", .false., ok, .false.)


    if (calcDataCube) call readDataCubeParameters(cLine, nLines)
    if (calcImage) call readImageParameters(cLine, nLines)
    if (calcSpectrum) call readSpectrumParameters(cLine, nLines)

    if (writeoutput) write(*,*) " "

    ! Close input file
    close(32)

  end subroutine inputs


  subroutine readGeometrySpecificParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
    logical :: ok

    select case(geometry)
       case("clumpyagb")
          call getReal("rinner", rinner, real(rsol/1.e10), cLine, nLines, &
               "Inner radius (solar radii): ","(a,1pe8.1,1x,a)", 1000., ok, .true.) 

          call getReal("router", router, rinner, cLine, nLines, &
               "Outer radius (solar radii): ","(a,1pe8.1,1x,a)", 1000., ok, .true.) 

          call getReal("vterm", vterm, 1.e5, cLine, nLines, &
               "Terminal velocity (km/s): ","(a,1pe8.1,1x,a)", 1000., ok, .true.) 

          call getReal("mdot", mdot, real(mSol) /( 365.25 * 24. * 3600.),  cLine, nLines, &
               "Mass-loss rate (solar masses per year): ","(a,1pe8.1,1x,a)", 1000., ok, .true.) 


    end select
  end subroutine readGeometrySpecificParameters
         

  subroutine readGridInitParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
    logical :: ok


    gridUsesAMR = .true.

    call getLogical("cylindrical", cylindrical, cLine, nLines, &
         "Grid uses 3D cylindical  coords: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("amr1d", amr1d, cLine, nLines, &
         "AMR grid is in one-dimensions only: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("amr2d", amr2d, cLine, nLines, &
         "AMR grid is in two-dimensions only: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("amr3d", amr3d, cLine, nLines, &
         "AMR grid is in three-dimensions: ","(a,1l,1x,a)", .false., ok, .false.)

    call getInteger("mindepthamr", minDepthAMR, cLine, nLines, "Minimum cell depth of AMR grid: ", &
         & "(a,i3,a)",5,ok,.true.)

    call getInteger("maxdepthamr", maxDepthAMR, cLine, nLines, "Maximum cell depth of AMR grid: ", &
         & "(a,i3,a)",31,ok,.true.)

    call getReal("amrgridsize", amrGridSize, 1., cLine, nLines, &
         "Size of adaptive mesh grid: ","(a,1pe8.1,1x,a)", 1000., ok, .true.) 
    call getDouble("amrgridcentrex", amrGridCentreX, 1.d0 , cLine, nLines, &
         "Grid centre X-coordinate: ","(a,es9.3,1x,a)", 0.0d0, ok, .false.) 
    call getDouble("amrgridcentrey", amrGridCentreY, 1.d0 , cLine, nLines, &
         "Grid centre Y-coordinate: ","(a,es9.3,1x,a)", 0.0d0, ok, .false.) 
    call getDouble("amrgridcentrez", amrGridCentreZ, 1.d0, cLine, nLines, &
         "Grid centre Z-coordinate: ","(a,es9.3,1x,a)", 0.0d0, ok, .false.) 

    call getDouble("limitscalar", limitScalar, 1.d0, cLine, nLines, &
         "Scalar limit for subcell division: ","(a,es9.3,1x,a)", 1000._db, ok, .false.) 

    call getDouble("limittwo", limitScalar2, 1.d0, cLine, nLines, &
         "Second scalar limit for subcell division: ","(a,es9.3,1x,a)", 0._db, ok, .false.) 

    call getString("geometry", geometry, cLine, nLines, &
         "Geometry: ","(a,a,a)","sphere",ok, .true.)

  end subroutine readGridInitParameters

  subroutine readDustPhysicsParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
    logical :: ok
    real :: grainFracTotal
    integer :: i
    character(len=20) :: grainTypeLabel, grainFracLabel, aMinLabel, &
         aMaxLabel, qDistLabel, pdistLabel, a0label

       oneKappa = .true.

       call getReal("dusttogas", dusttoGas, 1., cLine, nLines, &
            "Dust to gas ratio: ","(a,f5.3,a)",0.01,ok,.false.)

       call getInteger("ndusttype", nDustType, cLine, nLines,"Number of different dust types: ","(a,i12,a)",1,ok,.false.)
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
          call getString(grainTypeLabel, grainType(i), cLine, nLines, &
               "Grain type: ","(a,a,1x,a)","sil_dl", ok, .true.)

          call getReal(grainFracLabel, grainFrac(i), 1., cLine, nLines, &
               "Grain fractional abundance: ","(a,f8.5,1x,a)",1. , ok, .false.)
          grainFracTotal = grainFracTotal + grainFrac(i)

          call getReal(aminLabel, aMin(i), 1., cLine, nLines, &
               "Min grain size (microns): ","(a,f8.5,1x,a)", 0.005, ok,  .true.)

          call getReal(amaxLabel, aMax(i), 1., cLine, nLines, &
               "Max grain size (microns): ","(a,f10.5,1x,a)", 0.25, ok, .true.)

          call getReal(qDistLabel, qdist(i), 1., cLine, nLines, &
               "Grain power law: ","(a,f4.1,1x,a)", 3.5, ok, .true. )

          call getReal(a0Label, a0(i), 1., cLine, nLines, &
               "Scale length of grain size (microns): ","(a,f8.5,1x,a)", 1.0e20, ok, .false.)


          call getReal(pdistLabel, pdist(i), 1., cLine, nLines, &
               "Exponent for exponential cut off: ","(a,f4.1,1x,a)", 1.0, ok, .false. )
          if (writeoutput) write(*,*)
       enddo

  end subroutine readDustPhysicsParameters

  subroutine readAtomicPhysicsParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
    integer :: i 
    character(len=20) :: keyword
    logical :: ok

    call writeBanner("Atomic physics data","#",TRIVIAL)

    call getInteger("natom", nAtom, cLine, nLines, &
         "Number of model atoms to solve for: ","(a,i1,a)",1,ok,.true.)
    do i = 1, nAtom
       write(keyword, '(a,i1)') "atom",i
       call getString(keyword, atomFileName(i), cLine, nLines, &
            "Use atom filename: ","(a,a,1x,a)","none", ok, .true.)
    enddo
    call getInteger("itrans", itransLine, cLine, nLines, &
         "Index of line transition: ","(a,i4,a)",4,ok,.true.)

    call getInteger("iatom", itransAtom, cLine, nLines, &
         "Atom of line transition: ","(a,i4,a)",4,ok,.true.)



  end subroutine readAtomicPhysicsParameters

  subroutine readSourceParameters(cLine, nLines)
    character(len=80) :: cLine(:), message
    integer :: nLines
    integer :: i 
    character(len=20) :: keyword
    logical :: ok


    call writeBanner("Stellar source  data","#",TRIVIAL)

    call getInteger("nsource", inputNSource, cLine, nLines, &
         "Number of sources: ","(a,i2,a)",1,ok,.true.)
    do i = 1, inputnSource
       if (writeoutput) write(*,*) " "
       write(message,'(a,i1)') "Source number: ",i
       call writeInfo(message,TRIVIAL)
       write(keyword, '(a,i1)') "radius",i
       call getDouble(keyword, sourceRadius(i), rsol/1.d10, cLine, nLines, &
            "Source radius (solar radii) : ","(a,f7.2,a)",1.d0, ok, .true.)

       write(keyword, '(a,i1)') "teff",i
       call getDouble(keyword, sourceTeff(i), 1.d0, cLine, nLines, &
            "Source temperature (K) : ","(a,f8.0,a)",1.d0, ok, .true.)

       write(keyword, '(a,i1)') "mass",i
       call getDouble(keyword, sourceMass(i), mSol, cLine, nLines, &
            "Source mass (solar masses) : ","(a,f7.2,a)",1.d0, ok, .true.)

       write(keyword, '(a,i1)') "contflux",i
       call getString(keyword, inputcontFluxFile(i), cLine, nLines, &
            "Continuum flux file: ","(a,a,a)","none", ok, .true.)

       write(keyword, '(a,i1)') "sourcepos",i
       call getVector(keyword, sourcePos(i), 1.d0, cLine, nLines, &
            "Source position (10^10 cm): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)
    enddo

  end subroutine readSourceParameters

  subroutine readMolecularPhysicsParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
    logical :: ok

       call getLogical("lte", lte, cLine, nLines, &
            "Read in LTE grid: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("densitysubsample", densitysubsample, cLine, nLines, &
               "Use density interpolation: ","(a,1l,a)", .false., ok, .false.)
       call getLogical("suppresswarnings", suppressWarnings, cLine, nLines, &
            "Suppress Warnings: ","(a,l,1x,a)",.false., ok, .false.) 
       call getLogical("addnewmoldata", addnewmoldata, cLine, nLines, &
            "Add new molecular data to non-molecular grid: ","(a,l,1x,a)",.false., ok, .false.)
       call getString("moleculefile", moleculefile, cLine, nLines, &
            "Input molecule filename: ","(a,a,1x,a)","none", ok, .true.)
       call getReal("distance", gridDistance, real(pctocm), cLine, nLines, &
            "Grid distance (pc): ","(a,f6.1,1x,a)", 1., ok, .true.)
       call getInteger("initnray", initnray, cLine, nLines, &
               "Number of fixed rays for stage 1: ","(a,i4,a)", 100, ok, .false.)
       call getLogical("dongstep", dongstep, cLine, nLines, &
               "Use Ng Acceleration: ","(a,1l,a)", .false., ok, .false.)
       call getLogical("quasi", quasi, cLine, nLines, &
               "Use Quasirandom numbers: ","(a,1l,a)", .false., ok, .false.)
       call getReal("tolerance", tolerance, 1., cLine, nLines, &
            "Maximum Fractional Change in level populations:","(a,f4.1,1x,a)", 0.01, ok, .true.)
       call getReal("vturb", vturb, real(kmstoc), cLine, nLines, &
            "Subsonic turbulent velocity (km/s):","(a,f4.1,1x,a)", 0.3, ok, .true.)
       call getLogical("noturb", noturb, cLine, nLines, &
            "No microturbulence","(a,1l,a)",.false., ok, .false.)
       call getInteger("setmaxlevel", setmaxlevel, cLine, nLines, &
            "Maximum molecular level to be considered:","(a,i2,1x,a)", 0, ok, .false.)
       call getReal("molAbundance", molAbundance, 1., cLine, nLines, &
            "Molecular Abundance:","(a,e12.5,1x,a)", 1e-9, ok, .true.)
       call getLogical("isinlte", isinlte, cLine, nLines, &
            "Assume LTE: ", "(a,1l,1x,a)", .false., ok, .false.)
       call getReal("dusttogas", dusttoGas, 1., cLine, nLines, &
            "Dust to gas ratio: ","(a,f5.3,a)",0.01,ok,.false.)
       call getLogical("plotlevels", plotlevels, cLine, nLines, &
            "Plot Molecular Levels ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("realdust", realdust, cLine, nLines, &
            "Use realistic dust model: ", "(a,1l,1x,a)", .true., ok, .false.)

       call getLogical("doCOchemistry", doCOchemistry, cLine, nLines, &
            "Use drop profile to model CO depletion: ", "(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("useDust", useDust, cLine, nLines, &
            "Calculate continuum emission from dust:", "(a,1l,1x,a)", .false., ok, .true.)

       if(doCOchemistry) then

          call getReal("fracCOdepletion", x_D, 1., cLine, nLines, &
               "Fraction of CO depletion", "(a,1l,1x,a)", 0.1, ok, .true.)
          
       endif


       if(geometry .eq. 'molcluster') then
          call getString("sphdatafilename", sphdatafilename, cLine, nLines, &
               "Input sph data file: ","(a,a,1x,a)","sph.dat.ascii", ok, .true.)
       endif


  end subroutine readMolecularPhysicsParameters


  subroutine readPhotoionPhysicsParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
    logical :: ok

    call getLogical("quickthermal", quickThermal, cLine, nLines, &
         "Compute photoionization equilibrium: ","(a,1l,a)", .false., ok, .false.)

    call getLogical("usemetals", usemetals, cLine, nLines, &
         "Use metals in photoionization calculation: ","(a,1l,a)", .true., ok, .false.)


    call getReal("h_abund", h_abund, 1., cLine, nLines, &
         "Hydrogen abdunance: ","(a,1PF8.3,a)", &
         1., ok, .false.)

    call getReal("he_abund", he_abund, 1., cLine, nLines, &
         "Helium abdunance: ","(a,1PF8.3,a)", &
         0.1, ok, .false.)

    call getReal("c_abund", c_abund, 1., cLine, nLines, &
         "Carbon abdunance: ","(a,1PF8.3,a)", &
         22.e-5, ok, .false.)

    call getReal("n_abund", n_abund, 1., cLine, nLines, &
         "Nitrogen abdunance: ","(a,1PF8.3,a)", &
         4.e-5, ok, .false.)

    call getReal("o_abund", o_abund, 1., cLine, nLines, &
         "Oxygen abdunance: ","(a,1PF8.3,a)", &
         33.e-5, ok, .false.)

    call getReal("ne_abund", ne_abund, 1., cLine, nLines, &
         "Neon abdunance: ","(a,1PF8.3,a)", &
         5.e-5, ok, .false.)

    call getReal("s_abund", s_abund, 1., cLine, nLines, &
         "Sulphur abdunance: ","(a,1PF8.3,a)", &
         0.9e-5, ok, .false.)

  end subroutine readPhotoionPhysicsParameters


  subroutine readMolecularLoopParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
  end subroutine readMolecularLoopParameters

  subroutine readAtomicLoopParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
    logical :: ok

    call getLogical("lte", lte, cLine, nLines, &
         "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
    call getReal("vturb", vturb, real(kmsToC), cLine, nLines, &
         "Turbulent velocity (km/s):","(a,f4.1,1x,a)", 50., ok, .true.)


  end subroutine readAtomicLoopParameters

  subroutine readRadiativeEquilibriumParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
  end subroutine readRadiativeEquilibriumParameters

  subroutine readPhotoionEquilibriumParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
    logical :: ok

    call getInteger("nmonte", inputnMonte, cLine, nLines, &
         "Number of photons in image","(a,i8,a)", 0, ok, .false.)

    call getReal("taudiff", tauDiff, 1., cLine, nLines, &
         "Mininum optical depth of cell to be in diffusion approx : ","(a,f7.1,a)",100., ok, .false.)

  end subroutine readPhotoionEquilibriumParameters

  subroutine readHydrodynamicsParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
  end subroutine readHydrodynamicsParameters

  subroutine readDataCubeParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
    logical :: ok

    call getString("datacubefile", datacubeFilename, cLine, nLines, &
         "Output datacube  filename: ","(a,a,1x,a)","none", ok, .true.)
    call getReal("imageside", imageside, 1., cLine, nLines, &
         "Image size (x10^10cm):","(a,es7.2e1,1x,a)", 5e7, ok, .true.)
    call getInteger("npixels", npixels, cLine, nLines, &
         "Number of pixels per row: ","(a,i4,a)", 50, ok, .true.)
    call getInteger("nv", nv, cLine, nLines, &
         "Number of velocity bins ","(a,i4,a)", 50, ok, .true.)
    call getInteger("nSubpixels", nSubpixels, cLine, nLines, &
         "Subpixel splitting (0 denotes adaptive)","(a,i4,a)", 0, ok, .true.)
    call getInteger("itrans", itrans, cLine, nLines, &
         "Molecular Line Transition","(a,i4,a)", 1, ok, .true.)
    call getReal("beamsize", beamsize, 1., cLine, nLines, &
         "Beam size (arcsec): ","(a,f4.1,1x,a)", 1000., ok, .true.)
    call getDouble("rotateviewaboutx", rotateViewAboutX, 1.d0, cLine, nLines, &
         "Angle to rotate about X (deg): ","(a,f4.1,1x,a)", 0.d0, ok, .false.)
    call getDouble("rotateviewaboutz", rotateViewAboutZ, 1.d0, cLine, nLines, &
         "Angle to rotate about Z (deg): ","(a,f4.1,1x,a)", 0.d0, ok, .false.)
    call getDouble("maxVel", maxVel, 1.d0, cLine, nLines, &
         "Maximum Velocity Channel (km/s): ","(a,f4.1,1x,a)", 1.0d0, ok, .true.)
    call getDouble("centrevecx", centrevecx, 1.d0, cLine, nLines, &
         "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .true.)
    call getDouble("centrevecy", centrevecy, 1.d0, cLine, nLines, &
         "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .true.)
    call getDouble("centrevecz", centrevecz, 1.d0, cLine, nLines, &
         "Image Centre Coordinate (10^10cm): ","(a,1pe8.1,1x,a)", 0.d0, ok, .true.)
    call getLogical("wanttau", wanttau, cLine, nLines, &
         "Write Tau information to datacube: ","(a,1l,1x,a)", .false., ok, .false.)


  end subroutine readDataCubeParameters

  subroutine readImageParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
    logical :: ok
    integer :: i
    character(len=20) :: keyword


    call getBigInteger("nphotons", nphotons, cLine, nLines, &
         "Number of photons in image: ","(a,i8,a)", 10000, ok, .true.)

    call getInteger("nimage", nimage, cLine, nLines, &
         "Number of images to calculate: ","(a,i8,a)", 1, ok, .false.)
    if (nimage == 1) then

       call getString("imagefile", imageFilename(1), cLine, nLines, &
            "Output image  filename: ","(a,a,1x,a)","none", ok, .true.)
       call getReal("lambdaimage", lambdaImage(1),1., cLine, nLines, &
         "Wavelength for monochromatic image (A):","(a,f8.1,1x,a)", 6562.8, ok, .false.)
       call getString("imagetype", outputimageType(1), cLine, nLines, &
            "Type of output image: ","(a,a,1x,a)","none", ok, .true.)
       call getInteger("npixels", npixelsArray(1), cLine, nLines, &
            "Number of pixels per side in image","(a,i8,a)", 200, ok, .false.)
    else
       do i = 1, nImage

          write(keyword,'(a,i1.1)') "imagefile",i
          call getString(keyword, imageFilename(i), cLine, nLines, &
               "Output image  filename: ","(a,a,1x,a)","none", ok, .true.)
          write(keyword,'(a,i1.1)') "lambdaimage",i
          call getReal(keyword, lambdaImage(i),1., cLine, nLines, &
               "Wavelength for monochromatic image (A):","(a,f8.1,1x,a)", 6562.8, ok, .false.)
          write(keyword,'(a,i1.1)') "imagetype",i
          call getString(keyword, outputimageType(i), cLine, nLines, &
               "Type of output image: ","(a,a,1x,a)","none", ok, .true.)
          write(keyword,'(a,i1.1)') "npixels",i
          call getInteger(keyword, npixelsArray(i), cLine, nLines, &
               "Number of pixels per side in image","(a,i8,a)", 200, ok, .false.)
    enddo
 end if
       
  end subroutine readImageParameters

  subroutine readSpectrumParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
  end subroutine readSpectrumParameters



!-----------------------------------------------------------------------------------------


subroutine findReal(name, value, cLine, nLines, ok)
 implicit none
 character(len=*) :: name
 real :: value
 character(len=80) :: cLine(*)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  k = index(cline(i)," ")-1
  j = len_trim(name)
  if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
       ok = .true.
       read(cLine(i)(j+1:80),*) value
  endif
 end do
 end subroutine findReal

 subroutine findDouble(name, value, cLine, nLines, ok)
 implicit none
 character(len=*) :: name
 real(double) :: value
 character(len=80) :: cLine(*)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  j = len_trim(name)
  k = index(cline(i)," ")-1
  if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
       ok = .true.
       read(cLine(i)(j+1:80),*) value
  endif
 end do
 end subroutine findDouble

 subroutine findVECTOR(name, value, cLine, nLines, ok)
 implicit none
 character(len=*) :: name
 type(VECTOR) :: value
 character(len=80) :: cLine(*)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  j = len_trim(name)
  k = index(cline(i)," ")-1
  if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
       ok = .true.
       read(cLine(i)(j+1:80),*) value%x, value%y, value%z
  endif
 end do
end subroutine findVECTOR

subroutine findInteger(name, value, cLine, nLines, ok)
 implicit none
 character(len=*) :: name
 integer :: value
 character(len=80) :: cLine(*)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  j = len_trim(name)
  k = index(cline(i)," ")-1
  if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
       ok = .true.
       read(cLine(i)(j+1:),*) value
  endif
 end do
 end subroutine findInteger

subroutine findBigInteger(name, value, cLine, nLines, ok)
 implicit none
 character(len=*) :: name
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
       ok = .true.
       read(cLine(i)(j+1:),*) value
  endif
 end do
end subroutine findBigInteger

subroutine findLogical(name, value, cLine, nLines, ok)
 implicit none
 character(len=*) :: name
 logical :: value
 character(len=80) :: cLine(*)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  k = index(cline(i)," ")-1
  j = len_trim(name)
  if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
       ok = .true.
       read(cLine(i)(j+1:),*) value
  endif
 end do
 end subroutine findLogical

subroutine findString(name, value, cLine, nLines, ok)
 implicit none
 character(len=*) :: name
 character(len=*) :: value
 character(len=80) :: cLine(*)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  j = len_trim(name)
  k = index(cline(i)," ")-1
  if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
       ok = .true.
       read(cLine(i)(j+1:),*) value
  endif
 end do
 value = trim(value)
 end subroutine findString

subroutine findRealArray(name, value, cLine, nLines, ok)
 implicit none
 character(len=*) :: name
 real :: value(:)
 character(len=80) :: cLine(*)
 integer :: nLines
 logical :: ok
 integer :: i, j, k

 ok = .false.
 do i = 1, nLines
  j = len_trim(name)
  k = index(cline(i)," ")-1
  if (trim(cLine(i)(1:k)) .eq. name(1:j)) then
       ok = .true.
       read(cLine(i)(j+1:),*) value
  endif
 end do
 end subroutine findRealArray

 subroutine getInteger(name, ival, cLine, nLines, message, format, idef, ok, &
                       musthave)
  character(len=*) :: name
  integer :: ival
  character(len=80) :: cLine(*)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default
  logical :: musthave
  integer :: idef
  logical :: ok
  ok = .true.
  default = " "
  call findInteger(name, ival, cLine, nLines, ok)
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

 subroutine getBigInteger(name, ival, cLine, nLines, message, format, idef, ok, &
                       musthave)
  character(len=*) :: name
  integer(kind=bigInt) :: ival
  character(len=80) :: cLine(*)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default
  logical :: musthave
  integer :: idef
  logical :: ok
  ok = .true.
  default = " "
  call findBigInteger(name, ival, cLine, nLines, ok)
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

 subroutine getReal(name, rval, unitConversion, cLine, nLines, message, format, rdef, ok, &
                    musthave)
  character(len=*) :: name
  real :: rval
  real :: unitConversion
  logical :: musthave
  character(len=80) :: cLine(*)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default
  real :: rdef
  logical :: ok
  ok = .true.
  default = " "
  call findReal(name, rval, cLine, nLines, ok)
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


 subroutine getDouble(name, dval, unitConversion, cLine, nLines, message, format, ddef, ok, &
                    musthave)
  character(len=*) :: name
  real(double) :: dval, unitConversion
  logical :: musthave
  character(len=80) :: cLine(*)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default
  real(double) :: ddef
  logical :: ok
  ok = .true.
  default = " "
  call findDouble(name, dval, cLine, nLines, ok)
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

 subroutine getVector(name, dval, unitConversion, cLine, nLines, message, format, ddef, ok, &
                    musthave)
  character(len=*) :: name
  real(double) :: unitConversion
  type(VECTOR) :: dval, ddef
  logical :: musthave
  character(len=80) :: cLine(*)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default
  logical :: ok
  ok = .true.
  default = " "
  call findVECTOR(name, dval, cLine, nLines, ok)
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


 subroutine getString(name, rval, cLine, nLines, message, format, rdef, ok, &
                      musthave)
  character(len=*) :: name
  character(len=*) :: rval
  logical :: musthave
  character(len=80) :: cLine(*)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message, format
  character(len=10) :: default
  character(len=*) :: rdef
  logical :: ok
  ok = .true.
  default = " "
  call findString(name, rval, cLine, nLines, ok)
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


 subroutine getLogical(name, rval, cLine, nLines, message, cformat, rdef, ok, &
                       musthave)
  character(len=*) :: name
  logical :: rval
  logical :: musthave
  character(len=80) :: cLine(*)
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
  call findLogical(name, rval, cLine, nLines, ok)
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

 subroutine getRealArray(name, rval, unitConversion, cLine, nLines, message, rdef, ok, &
                      musthave)
  character(len=*) :: name
  real :: rval(:), unitConversion
  logical :: musthave
  character(len=80) :: cLine(*)
  character(len=100) :: output
  integer :: nLines
  character(len=*) :: message
  character(len=10) :: default
  real :: rdef
  logical :: ok
  ok = .true.
  default = " "
  call findRealArray(name, rval, cLine, nLines, ok)
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

