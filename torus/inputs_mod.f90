module inputs_mod

use messages_mod
use utils_mod
use input_variables
use jets_mod

implicit none

public

contains

  subroutine inputs()

    use constants_mod
    use vector_mod
    use unix_mod
    use input_variables         ! variables that would be passed as arguments

    implicit none

    integer :: nLines
    integer :: errNo, i
    logical :: ok
    character(len=80) :: cLine(200) 
    character(len=10) :: default

    character(len=20) :: grainTypeLabel, aminLabel, aMaxLabel, a0label
    character(len=20) :: qdistlabel, pdistlabel, grainFracLabel
    real :: grainFracTotal
    logical :: done

    ! character vars for unix environ

    character(len=80) :: dataDirectory

    character(len=20) :: keyword
!    character(len=80) :: message
    character(len=80) :: paramFile

    integer :: error

    oneKappa = .false.

    contrast = 1.
    grainSize = 1.
    nDustType = 1

    nBlobs = 0
    nLines = 0

    inputOK = .true.

    call unixGetEnv("TORUS_JOB_DIR",absolutePath)
    write(*,*) absolutePath
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
    call unixGetenv("TORUS_DATA",dataDirectory,i,errno)
    if (errNo /= 0) then
       write(*,'(A)') "! warning TORUS_DATA environment variable should be defined"
       stop
    endif


    call getInteger("verbosity", verbosityLevel, cLine, nLines, &
         "Verbosity level: ", "(a,i8,1x,a)", 1, ok, .false.)


    call getInteger("nphotons", nPhotons, cLine, nLines, &
         "Number of photons: ", "(a,i10,1x,a)", 100000, ok, .true.)

    call getInteger("idump", idump, cLine, nLines, &
         "Hydrodynamic dump number: ", "(a,i4,1x,a)", 1, ok, .false.)

    call getLogical("blockhandout", blockHandout, cLine, nLines, &
         "Use blockhandout for parallel computations ", "(a,1l,1x,a)", .true., ok, .false.)


    call getReal("distance", gridDistance, cLine, nLines, &
         "Grid distance (pc): ","(a,f4.1,1x,a)", 100., ok, .false.)

    call getReal("tminglobal", TMinGlobal, cLine, nLines, &
         "Minimum Temperature (K): ","(a,f4.1,1x,a)", 2.8, ok, .false.)

    call getReal("lamstart", lamstart, cLine, nLines, &
         "X-array start (angs or kms): ", "(a,f7.1,1x,a)", 2000., ok, .true.)

    call getReal("lamend", lamend, cLine, nLines, &
         "X-array end (angs or kms): ", "(a,f12.1,1x,a)", 10000., ok, .true.)


    call getLogical("wavlin", lamLinear, cLine, nLines, &
         "Linear wavelength array: ","(a,1l,1x,a)", .true., ok, .true.)

    call getLogical("lambdafile", lamfile, cLine, nLines, &
         "Wavelength grid from file: ","(a,1l,1x,a)", .false., ok, .false.)

    if (lamFile) then
       call getString("lamfilename", lamFilename, cLine, nLines, &
            "Wavelength grid filename: ","(a,a,1x,a)","none", ok, .true.)
    endif

    call getInteger("nlambda", nlambda, cLine, nLines, &
         "Number of wavelength/velocity bins: ", "(a,i4,1x,a)", 20, ok, .true.)


    call getReal("lambdatau", lambdatau, cLine, nLines, &
         "Lambda for tau test: ","(a,1PE10.3,1x,a)", 5500.0, ok, .false.)

    call getLogical("sed", sed, cLine, nLines, &
         "Write spectrum as lambda vs lambda Flambda: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("dustysed", dustysed, cLine, nLines, &
         "Write spectrum as normalized lambda Flambda: ","(a,1l,1x,a)", .false., ok, .false.)


    call getLogical("sised", sised, cLine, nLines, &
         "Write spectrum as lambda (microns) vs lambda F_lambda (microns * W/m^2):","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("jansky", jansky, cLine, nLines, &
         "Write spectrum in janskies: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("inarcsec", imageinArcsec, cLine, nLines, &
         "Write image distances in arcseconds: ","(a,1l,1x,a)", .false., ok, .false.)

    if (imageInArcSec) then
       call getReal("imagearcsec", imageSizeinArcsec, cLine, nLines, &
            "Image size in arcseconds: ","(a,f10.3,1x,a)", 0.130, ok, .true.)
    endif


    call getLogical("sed", sed, cLine, nLines, &
         "Write spectrum as normalized lambda Flambda: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("jansky", jansky, cLine, nLines, &
         "Write spectrum in janskies: ","(a,1l,1x,a)", .false., ok, .false.)


    call getInteger("npix", npix, cLine, nLines, &
         "Number of pixels for polimages: ", "(a,i3,1x,a)", 50, ok, .false.)


    call getString("distortion", distortionType, cLine, nLines, &
         "Distortion type: ","(a,a,1x,a)","none", ok, .false.)

    call getInteger("ninc", nInclination, cLine, nLines, &
         "Number of inclination angles: ", "(a,i3,1x,a)", 1, ok, .false.)

    allocate(inclinations(nInclination))
    call findRealArray("inclinations", inclinations, cLine, nLines, ok)
    if (ok) then
       call getRealArray("inclinations", inclinations, cLine, nLines, &
            "Inclinations (deg): ",'(a,4i4,1x,a)',90., ok, .false.)
       inclinations(:) = inclinations(:) * degToRad
    else
       deallocate(inclinations)
       call getReal("firstinc", firstInclination, cLine, nLines, &
            "First inclination angle (deg): ","(a,f4.1,1x,a)", 10., ok, .true.)
       firstInclination = firstInclination * degToRad
       if (nInclination > 1) &
            call getReal("lastinc", lastInclination, cLine, nLines, &
            "Last inclination angle (deg): ","(a,f4.1,1x,a)", 80., ok, .true.)
       lastInclination = lastInclination * degToRad
    end if

!    call getReal("scale", scale, cLine, nLines, &
!         "Scale (rsolar): ","(a,f6.1,1x,a)", 1000., ok, .false.) 


    call getString("gridtype", gridcoords, cLine, nLines, &
         "Grid type: ","(a,a,a)","cartesian",ok, .true.)

    call getLogical("gridusesamr", gridUsesAMR, cLine, nLines, &
         "Grid uses adaptive mesh refinement: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("cylindrical", cylindrical, cLine, nLines, &
         "Grid uses 3D cylindical  coords: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("amr1d", amr1d, cLine, nLines, &
         "AMR grid is in one-dimensions only: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("amr2d", amr2d, cLine, nLines, &
         "AMR grid is in two-dimensions only: ","(a,1l,1x,a)", .false., ok, .false.)

    call getLogical("amr3d", amr3d, cLine, nLines, &
         "AMR grid is in three-dimensions: ","(a,1l,1x,a)", .false., ok, .false.)

    call getInteger("mindepthamr", minDepthAMR, cLine, nLines, "Minimum cell depth of AMR grid: ", &
         & "(a,i3,a)",5,ok,.false.)

    call getInteger("maxdepthamr", maxDepthAMR, cLine, nLines, "Maximum cell depth of AMR grid: ", &
         & "(a,i3,a)",30,ok,.false.)


    if (gridUsesAMR) then
       call getReal("amrgridsize", amrGridSize, cLine, nLines, &
            "Size of adaptive mesh grid: ","(a,1pe8.1,1x,a)", 1000., ok, .true.) 
       call getDouble("amrgridcentrex", amrGridCentreX, cLine, nLines, &
            "Grid centre X-coordinate: ","(a,es9.3,1x,a)", 0.0d0, ok, .false.) 
       !        if (amrGridCentreX == 0.0) then 
       !           print *, 'WARNING: amrGridCentreX == 0. This may cause numerical problems!'
       !        end if
       call getDouble("amrgridcentrey", amrGridCentreY, cLine, nLines, &
            "Grid centre Y-coordinate: ","(a,es9.3,1x,a)", 0.0d0, ok, .false.) 
       !        if (amrGridCentreY == 0.0) then 
       !           print *, 'WARNING: amrGridCentreY == 0. This may cause numerical problems!'
       !        end if
       call getDouble("amrgridcentrez", amrGridCentreZ, cLine, nLines, &
            "Grid centre Z-coordinate: ","(a,es9.3,1x,a)", 0.0d0, ok, .false.) 
       !        if (amrGridCentreZ == 0.0) then 
       !           print *, 'WARNING: amrGridCentreZ == 0. This may cause numerical problems!'
       !        end if
       call getDouble("limitscalar", limitScalar, cLine, nLines, &
            "Scalar limit for subcell division: ","(a,es9.3,1x,a)", 1000._db, ok, .false.) 
       call getDouble("limittwo", limitScalar2, cLine, nLines, &
            "Second scalar limit for subcell division: ","(a,es9.3,1x,a)", 0._db, ok, .false.) 
       call getLogical("dosmoothgrid", doSmoothGrid, cLine, nLines, &
            "Smooth AMR grid: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("smoothgridtau", doSmoothGridtau, cLine, nLines, &
            "Smooth AMR grid using tau: ","(a,1l,1x,a)", .false., ok, .false.)
       if (dosmoothgridtau) then
          call getReal("lambdasmooth", lambdasmooth, cLine, nLines, &
               "Lambda for tau smoothing: ","(a,1PE10.3,1x,a)", 5500.0, ok, .true.)
          call getReal("taumax", tauSmoothMax, cLine, nLines, &
               "Maximum tau for smoothing: ","(a,f10.1,1x,a)", 1.0, ok, .true.)
          call getReal("taumin", tauSmoothMin, cLine, nLines, &
               "Minimum tau for smoothing: ","(a,f10.1,1x,a)", 1.0, ok, .true.)
       endif

       if (doSmoothGrid) then
          call getReal("smoothfactor", smoothFactor, cLine, nLines, &
               "Inter-cell maximum ratio before smooth: ","(a,f6.1,1x,a)", 5., ok, .true.)
       else
          smoothFactor = 0.0      
       end if
       call getReal("samplefreq", sampleFreq, cLine, nLines, &
            "Max samples per AMR subcell: ","(a,f6.1,1x,a)", 2., ok, .true.) 
       call getString("popfile", popFilename, cLine, nLines, &
            "Grid populations filename: ","(a,a,1x,a)","none", ok, .false.)
       call getLogical("writepops", writePops, cLine, nLines, &
            "Write populations file: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("readpops", readPops, cLine, nLines, &
            "Read populations file: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("writephasepops", writePhasePops, cLine, nLines, &
            "Write populations file at each phase: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("readphasepops", readPhasePops, cLine, nLines, &
            "Read populations file (specific phase): ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("stateq2d", statEq2d, cLine, nLines, &
            "Statistical equilibrium can be 2-D: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("stateq1stoctant", statEq1stOctant, cLine, nLines, &
            "amrStateq is performed only in 1st octant: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("nophaseupdate", noPhaseUpdate, cLine, nLines, &
            "Disable updating AMR grid at each phase: ","(a,1l,1x,a)", .true., ok, .false.)
       call getLogical("2donly", amr2dOnly, cLine, nLines, &
            "Only use 2D plane in AMR grid: ","(a,1l,1x,a)", .false., ok, .false.)
       if (amr2dOnly .and. .not. statEq2d) then
          if (writeoutput) write(*,'(a)') "WARNING: turning on statEq2d for amr2dOnly"
          statEq2d = .true.
       end if
       if (statEq2d .and. ( (amrGridCentreX < 0.0) .or. (amrGridCentreY < 0.0) ))  & 
            print *, 'WARNING: grid centre should probably set (x>0) and (y>0) if ',&
            'doing statEq2d'     
       if (statEq2d .and. statEq1stOctant)  then
          print *, 'Error:: statEq2d and statEq1stOctant cannot be both T. '
          print *, '        Change one or both in your parameter file!'
          stop
       end if

       if (amr2d) then
          if (writeoutput) write(*,*) "WARNING: amr grid is in two-d - switching off stateq2d"
          stateq2d = .false.
          amr2donly = .false.
       endif

       if (readPops .and. writePops) then
          if (writeoutput) write(*,'(a)') "WARNING: both readPops and writePops set."
       endif
       call getLogical("readfileformatted", readFileFormatted, cLine, nLines, &
            "Populations input file is formatted: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("writefileformatted", writeFileFormatted, cLine, nLines, &
            "Populations output file will be formatted: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("forcelinechange", forceLineChange, cLine, nLines, &
            "Recalculate opacities for different line transition: ","(a,1l,1x,a)", .false., ok, .false.)
! Set the default PGPLOT device to be used by plot_AMR_values
       call getString( "pgplotdevice", pgplotDevice, cline, nlines, &
            "Default PGPLOT device for plot_AMR_values: ", "(a,a,1x,a)", "ps/vcps", ok, .false. )

    else    

       if (gridcoords.eq."cartesian") then
          call getInteger("nx", nx, cLine, nLines, &
               "Grid size in x-direction: ", "(a,i3,1x,a)", 100, ok, .true.)
          call getInteger("ny", ny, cLine, nLines, &
               "Grid size in y-direction: ", "(a,i3,1x,a)", 100, ok, .true.)
          call getInteger("nz", nz, cLine, nLines, &
               "Grid size in z-direction: ", "(a,i3,1x,a)", 100, ok, .true.)
       endif

       if (gridcoords.eq."polar") then
          call getInteger("nr", nr, cLine, nLines, &
               "Grid size in r-direction: ", "(a,i3,1x,a)", 100, ok, .true.)
          call getInteger("nmu", nmu, cLine, nLines, &
               "Grid size in mu-direction: ", "(a,i3,1x,a)", 100, ok, .true.)
          call getInteger("nphi", nphi, cLine, nLines, &
               "Grid size in phi-direction: ", "(a,i3,1x,a)", 100, ok, .true.)
       endif

    end if

    call getString("geometry", geometry, cLine, nLines, &
         "Geometry: ","(a,a,a)","sphere",ok, .true.)

    !  call getLogical("raman", doRaman, cLine, nLines, &
    !            "Raman-scattering model: ","(a,1l,1x,a)", .false., ok, .false.)


    !  if (geometry(1:5) .eq. "torus") then
    !     call getReal("rtorus", rTorus, cLine, nLines, &
    !          "Radius of torus (AU): ","(a,f5.1,a)", 0.7, ok, .true.)
    !     call getReal("router", rOuter, cLine, nLines, &
    !          "Radius of torus x-section (AU): ","(a,f5.1,a)", 0.1, ok, .true.)
    !     call getReal("teff", teff, cLine, nLines, &
    !          "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)
    !     rTorus = rTorus * auTocm
    !     rOuter = rOuter * auTocm
    !  endif

    if (geometry == "rolf") then
       call getReal("mdot", mdot, cLine, nLines, &
            "Hot-star mass-loss rate (msol/yr): ","(a,1p,e12.5,a)", 0.7, ok, .true.)
       call getReal("vterm", vterm, cLine, nLines, &
            "Hot-star terminal velocity (km/s): ","(a,f5.1,a)", 0., ok, .true.)
       call getReal("o6width", o6width, cLine, nLines, &
            "FWHM of O6 line (km/s): ","(a,f5.1,a)", 0., ok, .true.)
       o6width = o6width * 1.e5
       vterm = vterm * 1.e5
       mdot = mdot * msol / (365.25 * 24. * 60. * 60.)
       doRaman = .true.
    endif


    if (trim(geometry) .eq. "dustblob") then
       call getReal("blobdistance", dustBlobDistance, cLine, nLines, &
            "Dust blob distance: ","(a,f5.1,a)", 0.5, ok, .true.)
       call getReal("phi", phiDustBlob, cLine, nLines, &
            "Dust blob azimuthal angle: ","(a,f5.1,a)", 0., ok, .true.)
       phiDustBlob = phiDustBlob * degToRad
       call getReal("xsize", xDustBlobSize, cLine, nLines, &
            "Dust blob x-size: ","(a,f5.1,a)", 0.1, ok, .true.)
       call getReal("ysize", yDustBlobSize, cLine, nLines, &
            "Dust blob y-size: ","(a,f5.1,a)", 0.1, ok, .true.)
       call getReal("zsize", zDustBlobSize, cLine, nLines, &
            "Dust blob z-size: ","(a,f5.1,a)", 0.1, ok, .true.)
    endif

    if (trim(geometry) .eq. "collide") then
       call getReal("rho", dustBlobDistance, cLine, nLines, &
            "Electron density: ","(a,1pe7.1,a)", 1.e10, ok, .true.)
       call getReal("binarysep", binarySep, cLine, nLines, &
            "Binary separation (cm): ","(a,1pe7.1,a)", 1.e13, ok, .true.)
       call getReal("momratio", momRatio, cLine, nLines, &
            "Wind momentum ratio (p/s): ","(a,f3.1,a)", 0.1, ok, .true.)
       call getReal("mdot", logMassLossRate, cLine, nLines, &
            "Log mDot: ","(a,f3.1,a)", -5., ok, .true.)
    endif

    if (distortionType.eq."binary") then
       call getReal("momratio", momRatio, cLine, nLines, &
            "Wind momentum ratio (p/s): ","(a,f3.1,a)", 0.1, ok, .true.)
    endif



    if (geometry .eq. "wr137") then
       call getReal("rcore", rCore, cLine, nLines, &
            "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)
       call getReal("temp", temp1, cLine, nLines, &
            "Wind temp (K): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("vterm", vTerm, cLine, nLines, &
            "Wind terminal velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("beta", beta, cLine, nLines, &
            "Wind beta law index: ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("mdot", mdot, cLine, nLines, &
            "mDot (msol/yr): ","(a,e7.1,a)", 1., ok, .true.)

       rCore = rCore * rSol
       vTerm = vTerm * 1.e5
       mdot = mdot*msol/(365.25*24.*3600.)

    endif

    if (geometry .eq. "wr104") then
       call getReal("massenv", massEnvelope, cLine, nLines, &
            "Envelope dust mass (Moon masses): ","(a,f5.2,a)", 10., ok, .true.)
       massEnvelope = massEnvelope * mMoon
    endif


    call getLogical("resonance", resonanceLine, cLine, nLines, &
         "Resonance line: ","(a,1l,1x,a)", .false., ok, .false.)

    if (geometry .eq. "resonance") then
       !     call getReal("rcore", rCore, cLine, nLines, &
       !          "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)
       !     call getReal("temp", temp1, cLine, nLines, &
       !          "Wind temp (K): ","(a,f7.0,a)", 1., ok, .true.)
       !     call getReal("vterm", vTerm, cLine, nLines, &
       !          "Wind terminal velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       !     call getReal("beta", beta, cLine, nLines, &
       !        "Wind beta law index: ","(a,f7.0,a)", 1., ok, .true.)
       !     call getReal("mdot", mdot, cLine, nLines, &
       !          "mDot (msol/yr): ","(a,e7.1,a)", 1., ok, .true.)
       !
       !     rCore = rCore * rSol
       !     vTerm = vTerm * 1.e5
       !     mdot = mdot*msol/(365.25*24.*3600.)

    endif

    if (geometry .eq. "spiralwind") then
       call getReal("rcore", rCore, cLine, nLines, &
            "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)
       call getReal("teff", teff, cLine, nLines, &
            "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("v0", v0, cLine, nLines, &
            "Wind base velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("vterm", vTerm, cLine, nLines, &
            "Wind terminal velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("beta", beta, cLine, nLines, &
            "Wind beta law index: ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("mdot", mdot, cLine, nLines, &
            "mDot (msol/yr): ","(a,e7.1,a)", 1., ok, .true.)
       rCore = rCore * rSol
       v0 = v0 * 1.e5
       vTerm = vTerm * 1.e5
       mdot = mdot*msol/(365.25*24.*3600.)
       vRot = vRot * 1.e5
       call getString("contflux", contFluxFile, cLine, nLines, &
            "Continuum flux filename: ","(a,a,1x,a)","none", ok, .true.)


    endif


    if (geometry .eq. "spiralwind") then
       call getReal("rcore", rCore, cLine, nLines, &
            "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)
       call getReal("teff", teff, cLine, nLines, &
            "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("v0", v0, cLine, nLines, &
            "Wind base velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("vterm", vTerm, cLine, nLines, &
            "Wind terminal velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("beta", beta, cLine, nLines, &
            "Wind beta law index: ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("mdot", mdot, cLine, nLines, &
            "mDot (msol/yr): ","(a,e7.1,a)", 1., ok, .true.)
       rCore = rCore * rSol
       v0 = v0 * 1.e5
       vTerm = vTerm * 1.e5
       mdot = mdot*msol/(365.25*24.*3600.)
       vRot = vRot * 1.e5
       call getString("contflux", contFluxFile, cLine, nLines, &
            "Continuum flux filename: ","(a,a,1x,a)","none", ok, .true.)


    endif


    if (geometry .eq. "spiralwind") then
       call getReal("rcore", rCore, cLine, nLines, &
            "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)
       call getReal("teff", teff, cLine, nLines, &
            "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("v0", v0, cLine, nLines, &
            "Wind base velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("vterm", vTerm, cLine, nLines, &
            "Wind terminal velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("beta", beta, cLine, nLines, &
            "Wind beta law index: ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("mdot", mdot, cLine, nLines, &
            "mDot (msol/yr): ","(a,e7.1,a)", 1., ok, .true.)
       rCore = rCore * rSol
       v0 = v0 * 1.e5
       vTerm = vTerm * 1.e5
       mdot = mdot*msol/(365.25*24.*3600.)
       vRot = vRot * 1.e5
       call getString("contflux", contFluxFile, cLine, nLines, &
            "Continuum flux filename: ","(a,a,1x,a)","none", ok, .true.)


    endif

    if ((geometry .eq. "puls").or.(geometry .eq. "wind")) then
       call getReal("rcore", rCore, cLine, nLines, &
            "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)
       call getReal("teff", teff, cLine, nLines, &
            "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("v0", v0, cLine, nLines, &
            "Wind base velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("vterm", vTerm, cLine, nLines, &
            "Wind terminal velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("vcont", vContrast, cLine, nLines, &
            "Polar/equator wind speed contrast: ","(a,f8.1,a)", 1., ok, .true.)
       call getReal("beta", beta, cLine, nLines, &
            "Wind beta law index: ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("mdot", mdot, cLine, nLines, &
            "mDot (msol/yr): ","(a,e7.1,a)", 1., ok, .true.)
       call getString("coreprofile", intProFilename, cLine, nLines, &
            "Core profile: ","(a,a,1x,a)","none", ok, .false.)
       call getReal("vrot", vRot, cLine, nLines, &
            "Rotational velocity (km/s): ","(a,f5.1,a)", 0., ok, .true.)
       call getLogical("lte", lte, cLine, nLines, &
            "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("lycontthick", LyContThick, cLine, nLines, &
            "Optically thick Lyman continuum:","(a,1l,1x,a)", .false., ok, .false.)
       call getString("contflux", contFluxFile, cLine, nLines, &
            "Continuum flux filename (primary): ","(a,a,1x,a)","none", ok, .true.)
       call getString("popfile", popFilename, cLine, nLines, &
            "Grid populations filename: ","(a,a,1x,a)","none", ok, .true.)
       call getLogical("writepops", writePops, cLine, nLines, &
            "Write populations file: ","(a,1l,1x,a)", .true., ok, .true.)
       call getLogical("readpops", readPops, cLine, nLines, &
            "Read populations file: ","(a,1l,1x,a)", .true., ok, .true.)
       call getInteger("nlower", nLower, cLine, nLines,"Lower level: ","(a,i2,a)",2,ok,.true.)
       call getInteger("nupper", nUpper, cLine, nLines,"Upper level: ","(a,i2,a)",3,ok,.true.)


       call getReal("xfac", xFac, cLine, nLines, &
            "Latitudinal x-fac: ","(a,f6.3,a)", 0., ok, .true.)


       rCore = rCore * rSol
       v0 = v0 * 1.e5
       vTerm = vTerm * 1.e5
       mdot = mdot*msol/(365.25*24.*3600.)
       vRot = vRot * 1.e5

    endif



    if (geometry.eq."binary") then
       call getLogical("lte", lte, cLine, nLines, &
            "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("lycontthick", LyContThick, cLine, nLines, &
            "Optically thick Lyman continuum:","(a,1l,1x,a)", .false., ok, .false.)
       call getReal("mass1", mass1, cLine, nLines, &
            "Primary mass (msol): ","(a,f5.1,a)", 1., ok, .true.)
       call getReal("radius1", radius1, cLine, nLines, &
            "Primary radius (rsol): ","(a,f5.1,a)", 1., ok, .true.)
       call getReal("mdot1", mdot1, cLine, nLines, &
            "Primary mDot (msol/yr): ","(a,f6.3,a)", 1., ok, .true.)
       call getReal("temp1", temp1, cLine, nLines, &
            "Primary wind temp (K): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("v01", vNought1, cLine, nLines, &
            "Primary wind base velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("vterm1", vTerm1, cLine, nLines, &
            "Primary wind terminal velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("beta1", beta1, cLine, nLines, &
            "Primary wind beta law index: ","(a,f7.0,a)", 1., ok, .true.)
       call getString("coreprofile1", intProFilename, cLine, nLines, &
            "Primary core profile: ","(a,a,1x,a)","none", ok, .false.)


       call getReal("mass2", mass2, cLine, nLines, &
            "Secondary mass (msol): ","(a,f5.1,a)", 1., ok, .true.)
       call getReal("radius2", radius2, cLine, nLines, &
            "Secondary radius (rsol): ","(a,f5.1,a)", 1., ok, .true.)
       call getReal("mdot2", mdot2, cLine, nLines, &
            "Secondary mDot (msol/yr): ","(a,f6.3,a)", 1., ok, .true.)
       call getReal("temp2", temp2, cLine, nLines, &
            "Secondary wind temp (K): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("v02", vNought2, cLine, nLines, &
            "Secondary wind base velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("vterm2", vTerm2, cLine, nLines, &
            "Secondary wind terminal velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("beta2", beta2, cLine, nLines, &
            "Secondary wind beta law index: ","(a,f7.0,a)", 1., ok, .true.)
       call getString("coreprofile2", intProFilename2, cLine, nLines, &
            "Secondary core profile: ","(a,a,1x,a)","none", ok, .false.)



       call getReal("period", period, cLine, nLines, &
            "Period (days): ","(a,f7.1,a)", 4., ok, .true.)

       call getReal("deflect", deflectionAngle, cLine, nLines, &
            "Wind-wind deflection angle (degs): ","(a,f7.1,a)", 0., ok, .false.)

       call getReal("shockwidth", shockWidth, cLine, nLines, &
            "Shock width (binary separations): ","(a,f7.1,a)", 0.1, ok, .true.)

       call getReal("shockfac", shockFac, cLine, nLines, &
            "Shock density: ","(a,f7.1,a)", 10., ok, .true.)

       call getString("contflux1", contFluxFile, cLine, nLines, &
            "Continuum flux filename (primary): ","(a,a,1x,a)","none", ok, .true.)
       call getString("contflux2", contFluxFile2, cLine, nLines, &
            "Continuum flux filename (secondary): ","(a,a,1x,a)","none", ok, .true.)
       call getInteger("nlower", nLower, cLine, nLines,"Lower level: ","(a,i2,a)",2,ok,.true.)
       call getInteger("nupper", nUpper, cLine, nLines,"Upper level: ","(a,i2,a)",3,ok,.true.)
       call getString("popfile", popFilename, cLine, nLines, &
            "Grid populations filename: ","(a,a,1x,a)","none", ok, .true.)
       call getLogical("writepops", writePops, cLine, nLines, &
            "Write populations file: ","(a,1l,1x,a)", .true., ok, .true.)
       call getLogical("readpops", readPops, cLine, nLines, &
            "Read populations file: ","(a,1l,1x,a)", .true., ok, .true.)



       period = period * 24. * 3600.
       mass1 = mass1 * mSol
       mass2 = mass2 * mSol
       radius1 = radius1 * rSol
       radius2 = radius2 * rSol
       vNought1 = vNought1 * 1.e5 
       vNought2 = vNought2 * 1.e5 
       vterm1 = vterm1 * 1.e5
       vterm2 = vterm2 * 1.e5
       mdot1 = 10.**mdot1*msol/(365.25*24.*3600.)
       mdot2 = 10.**mdot2*msol/(365.25*24.*3600.)
       deflectionAngle = deflectionAngle * degTorad


       call getReal("rho", rho, cLine, nLines, &
            "Density (xxx): ","(a,1p,e10.2,1p,1x,a)", 0., ok, .false.)
    endif

    if (geometry.eq."gammavel") then

       call getReal("mass1", mass1, cLine, nLines, &
            "Primary mass (msol): ","(a,f5.1,a)", 1., ok, .true.)
       call getReal("rstar1", rstar1, cLine, nLines, &
            "Primary radius (rsol): ","(a,f5.1,a)", 1., ok, .true.)
       call getReal("mdot1", mdot1, cLine, nLines, &
            "Primary log mDot (msol/yr): ","(a,f6.3,a)", 1., ok, .true.)
       call getReal("teff1", teff1, cLine, nLines, &
            "Primary wind temp (K): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("v01", vNought1, cLine, nLines, &
            "Primary wind base velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("vterm1", vTerm1, cLine, nLines, &
            "Primary wind terminal velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("beta1", beta1, cLine, nLines, &
            "Primary wind beta law index: ","(a,f7.0,a)", 1., ok, .true.)

       call getReal("mass2", mass2, cLine, nLines, &
            "Secondary mass (msol): ","(a,f5.1,a)", 1., ok, .true.)
       call getReal("rstar2", rstar2, cLine, nLines, &
            "Secondary radius (rsol): ","(a,f5.1,a)", 1., ok, .true.)
       call getReal("mdot2", mdot2, cLine, nLines, &
            "Secondary log mDot (msol/yr): ","(a,f6.3,a)", 1., ok, .true.)
       call getReal("teff2", teff2, cLine, nLines, &
            "Secondary wind temp (K): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("v02", vNought2, cLine, nLines, &
            "Secondary wind base velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("vterm2", vTerm2, cLine, nLines, &
            "Secondary wind terminal velocity (km/s): ","(a,f7.0,a)", 1., ok, .true.)
       call getReal("beta2", beta2, cLine, nLines, &
            "Secondary wind beta law index: ","(a,f7.0,a)", 1., ok, .true.)


       call getReal("binarysep", binarySep, cLine, nLines, &
            "Binary separation (AU): ","(a,f7.1,a)", 4., ok, .true.)

       call getReal("deflect", deflectionAngle, cLine, nLines, &
            "Wind-wind deflection angle (degs): ","(a,f7.1,a)", 0., ok, .false.)

       call getReal("shockwidth", shockWidth, cLine, nLines, &
            "Shock width (binary separations): ","(a,f7.1,a)", 0.1, ok, .false.)

       call getReal("shockfac", shockFac, cLine, nLines, &
            "Shock density: ","(a,f7.1,a)", 10., ok, .false.)

       call getString("contflux1", contFluxFile1, cLine, nLines, &
            "Continuum flux filename (primary): ","(a,a,1x,a)","none", ok, .true.)
       call getString("contflux2", contFluxFile2, cLine, nLines, &
            "Continuum flux filename (secondary): ","(a,a,1x,a)","none", ok, .true.)



       binarySep = binarySep * autocm / 1.e10
       mass1 = mass1 * mSol
       mass2 = mass2 * mSol
       rstar1 = rstar1 * rSol / 1.d10
       rstar2 = rstar2 * rSol / 1.d10
       vNought1 = vNought1 * 1.e5 
       vNought2 = vNought2 * 1.e5 
       vterm1 = vterm1 * 1.e5
       vterm2 = vterm2 * 1.e5
       mdot1 = mdot1*msol/(365.25*24.*3600.)
       mdot2 = mdot2*msol/(365.25*24.*3600.)
       deflectionAngle = deflectionAngle * degTorad

    endif

    if (geometry .eq. "shell") then
       call getReal("radius", radius, cLine, nLines, &
            "Stellar radius: ","(a,f5.1,a)", 100., ok, .true.)
       call getReal("shellfrac", shellfrac, cLine, nLines, &
            "Shell fraction: ","(a,f5.1,a)", 0.01, ok, .true.)
       call getReal("kfac", kfac, cLine, nLines, &
            "kfac: ", "(a,f5.1,a)", 0.5, ok, .true.)
       call getReal("rho", rho, cLine, nLines, &
            "Density (xxx): ","(a,1p,e10.2,1p,1x,a)", 1.e-6, ok, .true.)
    endif

    if (geometry .eq. "stateq") then

       call getReal("kfac", kfac, cLine, nLines, &
            "kfac: ", "(a,f5.1,a)", 0., ok, .true.)
       call getReal("densfac", scaledensity, cLine, nLines, &
            "Density scaling factor: ", "(a,f5.1,a)", 1., ok, .false.)
       call getString("opacityfile", opacityDataFile, cLine, nLines, &
            "Opacity file: ","(a,a,1x,a)","opacity_data", ok, .true.)
       call getString("coreprofile", intProFilename, cLine, nLines, &
            "Core profile: ","(a,a,1x,a)","none", ok, .false.)
       call getString("contflux", contFluxFile, cLine, nLines, &
            "Continuum flux filename: ","(a,a,1x,a)","none", ok, .true.)

    endif

    call getInteger("nblobs", nBlobs, cLine, nLines, &
         "Number of blobs: ", "(a,i3,1x,a)", 0, ok, .false.)
    if (nBlobs /= 0) then
       call getLogical("freshblobs", freshBlobs, cLine, nLines, &
            "Fresh Blobs: ","(a,1l,1x,a)", .true., ok, .true.)
       call getReal("blobcontrast", blobContrast, cLine, nLines, &
            "Blob contrast: ","(a,f7.1)", 1., ok, .true.)
    endif


    call getInteger("nphase", nPhase, cLine, nLines, &
         "Number of phases: ", "(a,i3,1x,a)", 1, ok, .false.)

    call getInteger("nstart", nStartPhase, cLine, nLines, &
         "Start at phase: ", "(a,i3,1x,a)", 1, ok, .false.)

    call getInteger("nend", nEndPhase, cLine, nLines, &
         "End at phase: ", "(a,i3,1x,a)", nPhase, ok, .false.)

    call getReal("phasetime", phaseTime, cLine, nLines, &
         "Time of each phase of simulation (s): ","(a,f7.1,a)", 30.*60., ok, .false.)

    call getReal("phaseoffset", phaseOffset, cLine, nLines, &
         "Phase offset: ","(a,f3.1,a)", 0., ok, .false.)

    call getLogical("forcerotate", forceRotate, cLine, nLines, &
         "Force view rotation ON","(a,1l,1x,a)", .false., ok, .false.)
    if (.not. forceRotate) then 
       call getLogical("forcenorotate", forceNoRotate, cLine, nLines, &
            "Force view rotation OFF","(a,1l,1x,a)", .false., ok, .false.)
    end if

    call getReal("probcont", probContPhoton, cLine, nLines, &
         "ProbContPhoton: ", "(a,f4.2,a)", 0.2, ok, .true.)


    if (geometry(1:6) .eq. "sphere") then
       call getReal("radius", radius, cLine, nLines, &
            "radius: ","(a,f5.1,a)", 0.5, ok, .true.)
       call getReal("kfac", kfac, cLine, nLines, &
            "kfac: ", "(a,f5.1,a)", 0.5, ok, .true.)
    endif

    if (geometry(1:7) .eq. "ellipse") then
       call getReal("rmaj", rmaj, cLine, nLines, &
            "rMajor (solar): ", "(a,f5.1,a)", 0.5, ok, .true.)
       call getReal("rmin", rmin, cLine, nLines, &
            "rMinor (solar): ", "(a,f5.1,a)", 0.5, ok, .true.)
       call getReal("rinner", rinner, cLine, nLines, &
            "rInner (solar): ", "(a,f5.1,a)", 0.5, ok, .true.)
       call getReal("rho", rho, cLine, nLines, &
            "Density (xxx): ","(a,1p,e10.2,1p,1x,a)", 1.e-6, ok, .true.)
       call getReal("teff", teff, cLine, nLines, &
            "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)
       rmaj = rmaj * rSol / 1.e10
       rmin = rmin * rSol / 1.e10
       rinner = rinner * rSol / 1.e10
    endif

    if (geometry(1:4) .eq. "disk") then

       call getReal("rcore", rCore, cLine, nLines, &
            "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)

       call getReal("rinner", rInner, cLine, nLines, &
            "Inner Radius (solar radii): ","(a,f5.1,a)", 12., ok, .true.)

       call getReal("router", rOuter, cLine, nLines, &
            "Outer Radius (solar radii): ","(a,f5.1,a)", 20., ok, .true.)

       call getReal("mcore", mCore, cLine, nLines, &
            "Core mass (solar masses): ","(a,f5.1,a)", 10., ok, .true.)

       call getReal("height", height, cLine, nLines, &
            "Scale height (solar radii): ","(a,1pe8.2,a)",1.e0,ok,.true.)

       call getReal("teff", teff, cLine, nLines, &
            "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)

       call getReal("disktemp", diskTemp, cLine, nLines, &
            "Disk temperature (K): ","(a,f9.1,a)",10000.,ok,.true.)

       call getReal("rho", rho, cLine, nLines, &
            "Electron Density (cm^-3): ","(a,1p,e10.2,1p,1x,a)", 1.e13, ok, .true.)

       call getString("coreprofile", intProFilename, cLine, nLines, &
            "Core profile: ","(a,a,1x,a)","none", ok, .false.)

       call getString("contflux", contFluxFile, cLine, nLines, &
            "Continuum flux filename: ","(a,a,1x,a)","none", ok, .true.)

       call getInteger("nspot", nSpot, cLine, nLines, &
            "Number of spots: ", "(a,i3,1x,a)", 0, ok, .false.)

       if (nSpot > 0) then

          call getReal("theta", thetaSpot, cLine, nLines, &
               "Theta of spot (degs): ","(a,f6.2,a)", 0., ok, .true.)

          call getReal("phi", phiSpot, cLine, nLines, &
               "Phi of spot (degs): ","(a,f6.2,a)", 0., ok, .true.)

          call getReal("tspot", tSpot, cLine, nLines, &
               "Temp of spot (degs): ","(a,f7.0,a)", 0., ok, .true.)

          call getReal("fspot", fSpot, cLine, nLines, &
               "Fractional coverage of spots: ","(a,f6.2,a)", 0., ok, .true.)

          call getLogical("photline", photLine, cLine, nLines, &
               "Line produced over whole star: ","(a,1l,1x,a)", .false., ok, .true.)


          thetaSpot = thetaSpot * degToRad
          phiSpot = phiSpot * degToRad


       endif


       rCore = rCore * rSol
       rInner = rInner * rSol
       rOuter = rOuter * rSol
       height = height * rSol
       mCore = mCore * mSol

    endif

    if (geometry .eq. "testamr") then

       call getReal("rcore", rCore, cLine, nLines, &
            "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)

       call getReal("rinner", rInner, cLine, nLines, &
            "Inner Radius (solar radii): ","(a,f5.1,a)", 12., ok, .true.)

       call getReal("router", rOuter, cLine, nLines, &
            "Outer Radius (inner radius): ","(a,f5.1,a)", 20., ok, .true.)

       call getReal("teff", teff, cLine, nLines, &
            "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)

       call getReal("rho", rho, cLine, nLines, &
            "Density (xxx): ","(a,1p,e10.2,1p,1x,a)", 0., ok, .true.)

       rCore = rCore * rSol
       rInner = rInner * rSol
       rOuter = rOuter * rInner
    endif

    if (geometry .eq. "proto") then

       call getReal("rcore", rCore, cLine, nLines, &
            "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)

       call getReal("rinner", rInner, cLine, nLines, &
            "Inner Radius (solar radii): ","(a,f5.1,a)", 12., ok, .true.)

       call getReal("router", rOuter, cLine, nLines, &
            "Outer Radius (inner radius): ","(a,f5.1,a)", 20., ok, .true.)

       call getReal("teff", teff, cLine, nLines, &
            "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)

       call getReal("taurad", tauRad, cLine, nLines, &
            "Radial optical depth at specified lambda: ","(a,1p,f10.2,1p,1x,a)", 1., ok, .true.)

       call getReal("rpower", rpower, cLine, nLines, &
            "Radial power-law index (r^-rpower): ","(a,1p,e10.2,1p,1x,a)", 1.5, ok, .true.)

       rho = 1.e-10
       rCore = rCore * rSol
       rInner = rInner * rSol
       rOuter = rOuter * rInner
    endif

    if (geometry .eq. "wrshell") then

       call getReal("rcore", rCore, cLine, nLines, &
            "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)

       call getreal("teff", teff, cLine, nLines, &
            "Source temperature (K): ","(a,f7.0,1x,a)", 3000., ok, .true.)

       call getReal("rinner", rInner, cLine, nLines, &
            "Inner Radius (stellar radii): ","(a,f5.1,a)", 12., ok, .true.)

       call getReal("router", rOuter, cLine, nLines, &
            "Outer Radius (inner radius): ","(a,f5.1,a)", 20., ok, .true.)

       call getReal("vterm", vterm, cLine, nLines, &
            "Terminal velocity (km/s): ","(a,f7.1,a)", 20., ok, .true.)

       call getReal("mdot", mdot, cLine, nLines, &
            "Mass-loss rate (solar masses/year): ","(a,f7.3,a)", 1., ok, .true.)

       call getString("contflux", contFluxFile, cLine, nLines, &
            "Continuum flux filename: ","(a,a,1x,a)","none", ok, .true.)

       call getReal("beta", beta, cLine, nLines, &
            "Wind beta law index: ","(a,f7.0,a)", 1., ok, .true.)

       rCore = rCore * rSol / 1.e10
       rInner = rInner * rCore
       rOuter = rOuter * rInner
       vTerm = vTerm * 1.e5
       mdot = mdot*mSol/(365.25*24.*3600.)
    endif



    ! if (geometry(1:4) .eq. "star") then 
    !   call getReal("radius", radius, cLine, nLines, &
    !       "radius: ","(a,f5.1,a)", 100., ok, .true.)
    !
    !   call getReal("mdot", mdot, cLine, nLines, &
    !       "Mdot (Msol/yr): ","(a,1pe8.2,a)",1.e-6,ok,.true.)
    !
    !
    !   call getReal("velocity", vterm, cLine, nLines, &
    !       "Wind velocity (km/s): ","(a,1pe8.2,a)",1.e1,ok,.true.)
    !
    !   call getReal("kfac", kfac, cLine, nLines, &
    !       "kfac: ", "(a,f5.1,a)", 0.5, ok, .true.)
    !
    !
    ! endif


    if (geometry == "raman") then
       call getReal("mdot", mdot, cLine, nLines, &
            "Mdot (Msol/yr): ","(a,1pe8.2,a)",1.e-6,ok,.true.)
       call getReal("vinf", vterm, cLine, nLines, &
            "Wind velocity (km/s): ","(a,1pe8.2,a)",1.e1,ok,.true.)
       call getReal("beta", beta, cLine, nLines, &
            "Wind velocity law index (beta): ","(a,f3.1,a)",0.,ok,.true.)
       call getReal("vo6", vo6, cLine, nLines, &
            "Vel width of O6 source (km/s): ","(a,1pe8.2,a)",20.,ok,.true.)
       call getReal("ramvel", ramVel, cLine, nLines, &
            "Absolute speed of O6 source (km/s): ","(a,1pe8.2,a)",0.,ok,.true.)
       vo6 = vo6 * 1.e5
       ramVel = ramVel * 1.e5
       vterm = vterm  * 1.e5
       doRaman = .true.
       mdot = mdot*mSol/(365.25*24.*3600.)
    endif

    if (distortionType .eq. "raman") then
       call getString("ramandist", ramanDist, cLine, nLines, &
            "Raman-scattering distortion type: ","(a,a,1x,a)","none", ok, .true.)
    endif

    if (distortionType .eq. "spiral") then 
       call getReal("vrot", vRot, cLine, nLines, &
            "Rotational velocity (km/s): ","(a,f5.1,a)", 100., ok, .true.)
       vRot = vRot * 1.e5
       call getInteger("nspiral", nSpiral, cLine, nLines, &
            "Number of spiral arms: ","(a,i3)",1,ok,.true.)
       call getReal("dipoleoffset", dipoleOffset, cLine, nLines, &
            "Dipole offset (degrees): ","(a,f7.1,1x,a)", 0.0e0, ok, .true.)
       dipoleOffset = dipoleOffset * degToRad
    endif

    if (distortionType .eq. "rotation") then 
       call getReal("vrot", vRot, cLine, nLines, &
            "Rotational velocity (km/s): ","(a,f5.1,a)", 100., ok, .true.)
       vRot = vRot * 1.e5
    endif

    ! 
    call getLogical("formalsol", formalsol, cLine, nLines, &
         "Solve for formal solution (no photon loops)?: ", "(a,1l,1x,a)", .false., ok, .false.)

    if (formalsol) then
       call getInteger("form_nphi", form_nphi, cLine, nLines, &
            "# of angular poins for formal integration : ","(a,i4,a)", 30, ok, .true.)
       call getInteger("form_nr_core", form_nr_core, cLine, nLines, &
            "# of radial poins for formal integration (core) : ","(a,i4,a)", 30, ok, .true.)
       call getInteger("form_nr_acc", form_nr_acc, cLine, nLines, &
            "# of radial poins for formal integration (accretion) : ","(a,i4,a)", 30, ok, .true.)
       call getInteger("form_nr_wind", form_nr_wind, cLine, nLines, &
            "# of radial poins for formal integration (wind) : ","(a,i4,a)", 30, ok, .true.)
       call getLogical("do_pos_disp", do_pos_disp, cLine, nLines, &
            "Position displacement calculation on?: ", "(a,1l,1x,a)", .false., ok, .false.)
    end if


    ! sub option for ttauri geometry
    call getLogical("ttau_acc_on", ttau_acc_on, cLine, nLines, &
         "Include TTauri magnetosphere?: ", "(a,1l,1x,a)", .true., ok, .false.)
    call getLogical("ttau_disc_on", ttau_disc_on, cLine, nLines, &
         "Include TTauri Disc?: ", "(a,1l,1x,a)", .false., ok, .false.)
    call getLogical("ttau_discwind_on", ttau_discwind_on, cLine, nLines, &
         "Include TTauri disc wind?: ", "(a,1l,1x,a)", .false., ok, .false.)
    call getLogical("ttau_jet_on", ttau_jet_on, cLine, nLines, &
         "Include TTauri jets?: ", "(a,1l,1x,a)", .false., ok, .false.)
    call getLogical("ttau_fuzzy_edge", ttau_fuzzy_edge, cLine, nLines, &
         "Uses fuzzy edge for accretion flow: ", "(a,1l,1x,a)", .false., ok, .false.)

    ! Use this paremeter to turn off the alpha disc when the 
    ! grid read in has alpha disc
    call getLogical("ttau_turn_off_disc", ttau_turn_off_disc, cLine, nLines, &
         "Alpha disc may be read in but turned off: ","(a,1l,1x,a)", .false., ok, .false.)
    call getLogical("ttau_turn_off_jet", ttau_turn_off_jet, cLine, nLines, &
         "Jet may be read in but turned off: ","(a,1l,1x,a)", .false., ok, .false.)
    call getLogical("ttau_turn_off_acc", ttau_turn_off_acc, cLine, nLines, &
         "Magnetosphere may be read in but turned off: ","(a,1l,1x,a)", .false., ok, .false.)


    call findLogical("mie", mie, cLine, nLines, ok)
    call findLogical("iso_scatter", isotropicScattering, cLine, nLines, ok)

    if (.not. ok) then
       mie = .true.
       default = " (default)"
       if (writeoutput) write(*,*) "Error:: mie and iso_scatter must be defined in your parameter file."
       stop
    endif


    if (mie) then
       call writeInfo("Scattering phase matrix: Mie",TRIVIAL)
    elseif (geometry == "ttauri".and. ttau_disc_on) then
       call writeInfo("Scattering phase matrix: Mie (disc) and Rayleigh (elsewhere)")
    else
       call writeInfo("Scattering phase matrix: Rayleigh")
    endif

    call getLogical("iso_scatter", isotropicScattering, cLine, nLines, &
         "Isotropic scattering function: ","(a,1l,1x,a)",.false.,ok,.false.)

    if  (isotropicScattering) then
       call writeWarning("ISOTROPIC SCATTERING PHASE MATRIX ENFORCED")
    endif

    call getLogical("readmiephase", readMiePhase, cLine, nLines, &
         "Read miePhase from file: ","(a,1l,1x,a)", .false., ok, .false.)
    call getLogical("writemiephase", writeMiePhase, cLine, nLines, &
         "Write miePhase to file: ","(a,1l,1x,a)", .false., ok, .false.)

    if (mie) then
       call getLogical("forcedwavelength", forcedWavelength, cLine, nLines, &
            "Forced photon wavelength: ","(a,1l,1x,a)", .false., ok, .false.)
       if (forcedWavelength) then
          call getReal("usephotonwavelength", usePhotonWavelength, cLine, nLines, &
               "Force photons to be produced at wavelength: ","(a,f4.2,a)", 5500., ok, .true.)
       endif
    endif
    if (mie) then
       call getLogical("lucyrad", lucyRadiativeEq, cLine, nLines, &
            "Lucy radiative equ.: ","(a,1l,1x,a)", .false., ok, .true.)
       call getString("lucyfilein", lucyFilenameIn, cLine, nLines, &
            "Input Lucy grid filename: ","(a,a,1x,a)","none", ok, .false.)
       call getString("lucyfileout", lucyFilenameOut, cLine, nLines, &
            "Output Lucy grid filename: ","(a,a,1x,a)","none", ok, .false.)
       call getLogical("writelucy", writeLucy, cLine, nLines, &
            "Write lucy grid file: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("readlucy", readLucy, cLine, nLines, &
            "Read lucy grid file: ","(a,1l,1x,a)", .false., ok, .false.)

       call getLogical("redolucy", redolucy, cLine, nLines, &
            "Redo lucy radiative equilibrium routine from a input file: ", &
            "(a,1l,1x,a)", .false., ok, .false.)
       if (writeLucy)  then
          call getLogical("writefileformatted", writeFileFormatted, cLine, nLines, &
               "Grid + lucy temperature data output file will be formatted: ","(a,1l,1x,a)", &
               .false., ok, .false.)
       end if
       if (readLucy) then
          call getLogical("readfileformatted", readFileFormatted, cLine, nLines, &
               "Grid + lucy temperature input file is read as formatted: ","(a,1l,1x,a)",  &
               .false., ok, .false.)
       end if


       if (lucyRadiativeEq) then
          call getInteger("nlucy", nLucy, cLine, nLines,"Number of photons per lucy iteration: ","(a,i12,a)",20000,ok,.false.)
          call getInteger("iterLucy", iterLucy, cline, nlines, "Minimum number of Lucy iterations: ", "(a,i3,a)",5,ok,.false.)

          call getLogical("forceLucyConv", forceLucyConv, cLine, nLines, &
               "Force convergence of Lucy algorithm: ","(a,1l,1x,a)", .false., ok, .false.)

          call getReal("lucy_undersampled", lucy_undersampled, cLine, nLines, &
               "Minimum percentage of undersampled cell in lucy iteration: ", &
               "(a,f4.2,a)",0.0,ok,.false.)

          call getReal("diffdepth", diffDepth, cLine, nLines, &
               "Depth of diffusion zone (in Rosseland optical depths): ", &
               "(a,f5.1,a)",10.0,ok,.false.)

          call getInteger("mincrossings", minCrossings, cLine, nLines, &
               "Minimum crossings required for cell to be sampled: ","(a,i12,a)",5,ok,.false.)
       endif
       call getReal("probdust", probDust, cLine, nLines, &
            "Probability of photon from dusty envelope: ","(a,f4.2,a)", 0.8, ok, .true.)

       call getReal("tthresh", tthresh, cLine, nLines, &
            "Temperature threshold for dust (K): ","(a,f10.2,a)", 1000., ok, .true.)

       oneKappa = .true.

       call getReal("dusttogas", dusttoGas, cLine, nLines, &
            "Dust to gas ratio: ","(a,f5.3,a)",0.01,ok,.false.)

       call getInteger("ndusttype", nDustType, cLine, nLines,"Number of different dust types: ","(a,i12,a)",1,ok,.false.)
       if (nDustType .gt. maxDustTypes) then
          if (writeoutput) write (*,*) "Max dust types exceeded: ", maxDustTypes
          stop
       end if

       call getLogical("hydro", solveVerticalHydro, cLine, nLines, &
            "Solve vertical hydrostatical equilibrium: ","(a,1l,1x,a)", .false., ok, .false.)

       if (solveVerticalHydro) then
          call getInteger("nhydro", nhydro, cLine, nLines, &
               "Max number of hydro iterations : ","(a,i4,a)", 5, ok, .true.)
       endif

       call getLogical("dustfile", dustfile, cLine, nLines, &
            "Get dust properties from file: ","(a,1l,1x,a)", .false., ok, .false.)
       if (dustfile) then
          call getString("kappafile", dustFilename(1), cLine, nLines, &
               "Dust properties filename: ","(a,a,1x,a)","none", ok, .true.)
       endif

       !   if (nDustType > 1) then
       !
       !      do i = 1, nDustType
       !         write(keyword,'(a,i1)') "kappafile",i
       !         write(message,'(a,i1,a)') "Dust Properties filename ",i,": "
       !         call getString(trim(keyword), dustFilename(i), cLine, nLines, &
       !              message,"(a,a,1x,a)","none", ok, .true.)
       !      enddo
       !      
       !   endif
    endif

    !
    ! Now using new grain parameters e.g. grainTypeLabel which will be read in later....
    ! The following should be removed later as they no longer needed.

    ! ! if (mie .or. (geometry == "ttauri" .and. ttau_disc_on)) then
    !  if (mie .or. (geometry == "ttauri")) then
    !      call getString("graintype", grainType, cLine, nLines, &
    !           "Grain type: ","(a,a,1x,a)","sil_dl", ok, .true.)

    !      ! read the relative abundances (which will be normalized later.)
    !      call getReal("x_sil_ow", X_grain(1), cLine, nLines, &
    !           "Abundance(Si-Ow): ","(a,1pe8.2,a)",0.,ok,.false.)
    !      call getReal("x_sil_oc", X_grain(2), cLine, nLines, &
    !           "Abundance(Si-Oc): ","(a,1pe8.2,a)",0.,ok,.false.)
    !      call getReal("x_sil_dl", X_grain(3), cLine, nLines, &
    !           "Abundance(Si-DL): ","(a,1pe8.2,a)",0.,ok,.false.)
    !      call getReal("x_amc_hn", X_grain(4), cLine, nLines, &
    !           "Abundance(amC-Hn): ","(a,1pe8.2,a)",0.,ok,.false.)
    !      call getReal("x_sil_pg", X_grain(5), cLine, nLines, &
    !           "Abundance(Si-Pg): ","(a,1pe8.2,a)",0.,ok,.false.)
    !      call getReal("x_gr1_dl", X_grain(6), cLine, nLines, &
    !           "Abundance(grf1-DL): ","(a,1pe8.2,a)",0.,ok,.false.)
    !      call getReal("x_gr2_dl", X_grain(7), cLine, nLines, &
    !           "Abundance(grf2-DL): ","(a,1pe8.2,a)",0.,ok,.false.)        


    !  endif



    if (.not.(mie)) then
       default = " "
       call findReal("kappasca", inputKappaSca, cLine, nLines, ok)
       if (ok) then
          if (writeoutput) write(*,'(a,1pe12.4)') "Scattering cross-section: ",inputKappaSca
       endif
       default = " "
       call findReal("kappaabs", inputKappaAbs, cLine, nLines, ok)
       if (ok) then
          if (writeoutput) write(*,'(a,1pe12.4)') "Absorption cross-section: ",inputKappaAbs
       endif
    endif


    call getLogical("secondsource", secondSource, cLine, nLines, &
         "Second source: ","(a,1l,a)", .false., ok, .false.)

    if (secondSource) then
       secondSourcePosition = VECTOR(0.,0.,0.)
       call getReal("secondpos", secondSourcePosition%x, cLine, nLines, &
            "Second source distance: ","(a,e7.1,1x,a)", 1.e14, ok, .true.)
    endif

    call getLogical("coreline", coreEmissionLine, cLine, nLines, &
         "Core Emission Line: ","(a,1l,a)", .false., ok, .false.)


    call getLogical("photoionization", photoionization, cLine, nLines, &
         "Compute photoionization equilibrium: ","(a,1l,a)", .false., ok, .false.)

 call getLogical("molecular", molecular, cLine, nLines, &
   "Compute molecular line transport: ","(a,1l,a)", .false., ok, .false.)

 if (molecular) then
    onekappa = .true.
   call getString("molfilein", molFilenameIn, cLine, nLines, &
        "Input Lucy grid filename: ","(a,a,1x,a)","none", ok, .false.)
   call getString("molfileout", molFilenameOut, cLine, nLines, &
        "Output Lucy grid filename: ","(a,a,1x,a)","none", ok, .false.)
   call getLogical("writemol", writeMol, cLine, nLines, &
        "Write lucy grid file: ","(a,1l,1x,a)", .false., ok, .false.)
   call getLogical("readmol", readMol, cLine, nLines, &
          "Read molecular grid file: ","(a,1l,1x,a)", .false., ok, .false.)
   call getLogical("writeLucy", writeLucy, cLine, nLines, &
        "Write lucy grid file: ","(a,1l,1x,a)", .false., ok, .false.)
   call getLogical("readLucy", readLucy, cLine, nLines, &
          "Read molecular grid file: ","(a,1l,1x,a)", .false., ok, .false.)
   call getLogical("openlucy", openLucy, cLine, nLines, &
        "Open Existing lucy file: ","(a,1l,1x,a)", .false., ok, .false.)
    call getReal("distance", gridDistance, cLine, nLines, &
         "Grid distance (pc): ","(a,f6.1,1x,a)", 1., ok, .true.)
    gridDistance = gridDistance * pcTocm   ! cm
    call getReal("tolerance", tolerance, cLine, nLines, &
         "Maximum Fractional Change in level populations:","(a,f4.1,1x,a)", 0.01, ok, .true.)
    call getReal("molAbundance", molAbundance, cLine, nLines, &
         "Molecular Abundance:","(a,es6.2,1x,a)", 1e-9, ok, .true.)
    call getLogical("useDust", useDust, cLine, nLines, &
         "Calculate continuum emission from dust:", "(a,1l,1x,a)", .false., ok, .false.)
    call getLogical("isinLTE", isinlte, cLine, nLines, &
         "Assume LTE: ", "(a,1l,1x,a)", .false., ok, .false.)
    call getReal("dusttogas", dusttoGas, cLine, nLines, &
         "Dust to gas ratio: ","(a,f5.3,a)",0.01,ok,.false.)
    

! Image parameters
    if(readmol) then

       call getReal("imageside", imageside, cLine, nLines, &
            "Image size (x10^10cm):","(a,es7.2e1,1x,a)", 5e7, ok, .true.)
       call getInteger("npixels", npixels, cLine, nLines, &
            "Number of pixels per row: ","(a,i4,a)", 50, ok, .true.)
       call getInteger("nv", nv, cLine, nLines, &
            "Number of velocity bins ","(a,i4,a)", 50, ok, .true.)
       call getInteger("nSubpixels", nSubpixels, cLine, nLines, &
            "Subpixel splitting (0 denotes adaptive)","(a,i4,a)", 0, ok, .true.)
       call getInteger("itrans", itrans, cLine, nLines, &
            "Molecular Line Transition","(a,i4,a)", 1, ok, .true.)
       call getReal("beamsize", beamsize, cLine, nLines, &
            "Beam size (arcsec): ","(a,f4.1,1x,a)", 1000., ok, .true.)
       call getReal("inclination", inc, cLine, nLines, &
            "Inclination angle (deg): ","(a,f4.1,1x,a)", -1., ok, .true.)
       call getDouble("maxVel", maxVel, cLine, nLines, &
            "Maximum Velocity Channel (km/s): ","(a,f4.1,1x,a)", 1.0d0, ok, .true.)

       if(writelucy .or. readlucy) then
          call getString("lucyfilein", lucyFilenameIn, cLine, nLines, &
               "Input Lucy grid filename: ","(a,a,1x,a)","none", ok, .false.)
          call getString("lucyfileout", lucyFilenameOut, cLine, nLines, &
               "Output Lucy grid filename: ","(a,a,1x,a)","none", ok, .false.)
       endif

    endif

    if(geometry .eq. 'molecularFilament') then

       call getReal("r0", r0, cLine, nLines, &
            "Core Radius of Filament:","(a,f4.1,1x,a)", 1., ok, .true.)
       call getReal("rhoC", rhoC, cLine, nLines, &
            "Central density along filament axis:","(a,f4.1,1x,a)", 1., ok, .true.)
    endif

    if((geometry .eq. 'molebench') .or. (geometry .eq. 'h2obench1') .or. &
       (geometry .eq. 'h2obench2') .or. (geometry .eq. 'AGBStar')) then

       call getReal("rinner", rInner, cLine, nLines, &
            "Inner Radius for dumpresults (10^10cm): ","(a,f5.1,a)", 1e4, ok, .true.)
       
       call getReal("router", rOuter, cLine, nLines, &
            "Outer Radius (10^10cm): ","(a,f5.1,a)", 1e6, ok, .true.)
    endif

endif

    call getLogical("hydrodynamics", hydrodynamics, cLine, nLines, &
         "Do hydrodynamics: ","(a,1l,a)", .false., ok, .false.)

    call getReal("cfl", cflNumber, cLine, nLines, &
         "Courant number:","(a,f4.1,1x,a)", 0.3, ok, .false.)


    if (hydrodynamics) then
       call getReal("x1", x1, cLine, nLines, &
         "Start of x-array:","(a,f4.1,1x,a)", 0.01, ok, .true.)
       call getReal("x2", x2, cLine, nLines, &
         "Start of x-array:","(a,f4.1,1x,a)", 0.01, ok, .true.)
    endif


    call getLogical("cmf", cmf, cLine, nLines, &
         "Compute CMF statistical equilibrium: ","(a,1l,a)", .false., ok, .false.)

    if (cmf) then
       call getString("lucyfilein", lucyFilenameIn, cLine, nLines, &
            "Input Lucy grid filename: ","(a,a,1x,a)","none", ok, .false.)
       call getString("lucyfileout", lucyFilenameOut, cLine, nLines, &
            "Output Lucy grid filename: ","(a,a,1x,a)","none", ok, .false.)
       call getLogical("writelucy", writeLucy, cLine, nLines, &
            "Write lucy grid file: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("readlucy", readLucy, cLine, nLines, &
            "Read lucy grid file: ","(a,1l,1x,a)", .false., ok, .false.)
       call getInteger("natom", nAtom, cLine, nLines, &
            "Number of model atoms to solve for: ","(a,i12,a)",1,ok,.true.)
       do i = 1, nAtom
          write(keyword, '(a,i1)') "atom",i
          call getString(keyword, atomFileName(i), cLine, nLines, &
               "Use atom filename: ","(a,a,1x,a)","none", ok, .true.)
       enddo
       call getInteger("itrans", itransLine, cLine, nLines, &
            "Index of line transition: ","(a,i12,a)",4,ok,.true.)

       call getInteger("iatom", itransAtom, cLine, nLines, &
            "Index of line transition: ","(a,i12,a)",4,ok,.true.)

    endif

    call getLogical("debug", debug, cLine, nLines, &
         "Write debug output: ","(a,1l,a)", .false., ok, .false.)


    call getLogical("thickcont", opticallyThickContinuum, cLine, nLines, &
         "Continuum is optically thick: ","(a,1l,a)", .false., ok, .false.)

    if (photoionization) then
       call getInteger("nlucy", nLucy, cLine, nLines,"Number of photons per lucy iteration: ","(a,i12,a)",20000,ok,.false.)
       call getReal("lucy_undersampled", lucy_undersampled, cLine, nLines, &
            "Minimum percentage of undersampled cell in lucy iteration: ", &
            "(a,f4.2,a)",0.0,ok,.false.)

       call getReal("diffdepth", diffDepth, cLine, nLines, &
            "Depth of diffusion zone (in Rosseland optical depths): ", &
            "(a,f5.1,a)",10.0,ok,.false.)

       call getInteger("mincrossings", minCrossings, cLine, nLines, &
            "Minimum crossings required for cell to be sampled: ","(a,i12,a)",5,ok,.true.)
    endif

    if (geometry == "starburst") then
       photoionization = .true.
       onekappa = .true.
    endif


    if (geometry == "lexington") then
       photoionization = .true.
       onekappa = .true.
       call getReal("rinner", rInner, cLine, nLines, &
            "Inner Radius (10^17cm): ","(a,f5.1,a)", 30., ok, .true.)
       rinner = rinner * 1.e7
    endif


    if (geometry == "ttauri" .or. geometry == "magstream") then

       if (geometry == "magstream") then

          call getString("magstreamfile", MagStreamFile, cLine, nLines, &
               "Magnetic stream data filename: ","(a,a,1x,a)","none", ok, .true.)
          call getLogical("magfiledegrees", magStreamFileDegrees, cLine, nLines, &
               "Magnetic stream data angle in degrees: ","(a,1l,a)", .true., ok, .true.)        
          call getLogical("isothermstream", isothermStream, cLine, nLines, &
               "Use isothermal temperature in accretion stream:","(a,1l,1x,a)", .false., ok, .true.)
          call getReal("isothermtemp", isothermTemp, cLine, nLines,"Isothermal temperature (K): ","(a,e9.3,1x,a)",7500.,ok,.true.)

          call getLogical("limitspottemp", limitSpotTemp, cLine, nLines, &
               "Limit on accretion spot temperature: ","(a,1l,1x,a)", .false., ok, .false.)
          if (limitSpotTemp) then
             call getReal("maxspottemp", maxSpotTemp, cLine, nLines, &
                  "Maximum hotspot temperature (K): ","(a,f7.1,1x,a)", 7000., ok, .true.)
          end if
          call getLogical("scaleflowrho", scaleFlowRho, cLine, nLines, &
               "Rescale accretion flow densities: ","(a,1l,1x,a)", .false., ok, .false.)
          if (scaleFlowRho) then
             call getReal("flowrhoscale", flowRhoScale, cLine, nLines, &
                  "Scaling factor for accretion flow densities: ","(a,e14.3,1x,a)", 1., ok, .true.)
          end if

       endif

       call getReal("ttaurirstar", TTauriRstar, cLine, nLines, &
            "T Tauri stellar radius (in R_sol): ","(a,f7.1,1x,a)", 2.0, ok, .true.)
       TTauriRstar = TTauriRstar * rSol ! [cm]
       rcore = TTauriRstar/1.0e10       ! [10^10cm]
       call getReal("ttaurimstar", TTauriMstar, cLine, nLines, &
            "T Tauri stellar mass (in M_sol): ","(a,f7.1,1x,a)", 0.8, ok, .true.)
       TTauriMstar = TTauriMstar * mSol
       call getReal("ttaurirouter", TTauriRouter, cLine, nLines, &
            "T Tauri outer flow radius (in R_star): ","(a,f7.1,1x,a)", 3.0, ok, .true.)
       TTauriRouter = TTauriRouter * TTauriRstar
       call getReal("ttauririnner", TTauriRinner, cLine, nLines, &
            "T Tauri inner flow radius (in R_star): ","(a,f7.1,1x,a)", 2.2, ok, .true.)
       TTauriRinner = TTauriRinner * TTauriRstar
       call getReal("ttauridiskheight", TTauriDiskHeight, cLine, nLines, &
            "T Tauri disk height (in R_star): ","(a,f7.1,1x,a)", 6.e-2, ok, .false.)
       TTauriDiskHeight = TTauriDiskHeight * TTauriRstar
       call getReal("ttauridiskrin", TTauriDiskRin, cLine, nLines, &
            "T Tauri disk inner radius  (in R_star): ","(a,f7.1,1x,a)", TTauriRouter/TTauriRstar, ok, .false.)
       call getReal("thindiskrin", ThinDiskRin, cLine, nLines, &
            "Thin disk inner radius  (in R_star): ","(a,f7.1,1x,a)", TTauriRouter/TTauriRstar, ok, .false.)
       call getReal("curtainsphi1s", curtainsPhi1s, cLine, nLines, &
            "Curtains 1: Phi start: (degrees): ","(a,f7.1,1x,a)", 30.0, ok, .false.)
       call getReal("curtainsphi1e", curtainsPhi1e, cLine, nLines, &
            "Curtains 1: Phi end: (degrees): ","(a,f7.1,1x,a)", 150.0, ok, .false.)
       call getReal("curtainsphi2s", curtainsPhi2s, cLine, nLines, &
            "Curtains 2: Phi start: (degrees): ","(a,f7.1,1x,a)", 210.0, ok, .false.)
       call getReal("curtainsphi2e", curtainsPhi2e, cLine, nLines, &
            "Curtains 2: Phi end: (degrees): ","(a,f7.1,1x,a)", 330.0, ok, .false.)
       !  converting the angles in radians  (RK) 
       curtainsPhi1s =    curtainsPhi1s * (pi/180.0) 
       curtainsPhi1e =    curtainsPhi1e * (pi/180.0) 
       curtainsPhi2s =    curtainsPhi2s * (pi/180.0) 
       curtainsPhi2e =    curtainsPhi2e * (pi/180.0) 

       ! The following two are used for "constantcurtans" geometry  (RK)
       call getInteger("curtain_number", curtain_number, cLine, nLines, &
            "Number of curtains : ","(a,i8,a)", 2, ok, .false.)
       call getReal("curtain_width", curtain_width, cLine, nLines, &
            "Width of each curtain (degree) : ","(a,f7.1,1x,a)", 120.0, ok, .false.)
       ! converting the curtain width from degrees to radians.
       curtain_width =  curtain_width*Pi/180.0


       call getString("mdottype", mDotType, cLine, nLines, &
            "T Tauri accretion rate model: ","(a,a,1x,a)","constant", ok, .true.)
       call getReal("mdotpar1", MdotParameter1, cLine, nLines, &
            "1st parameter for accretion rate: ", "(a,e9.3,1x,a)", 1.0, ok, .true.)
       call getReal("mdotpar2", MdotParameter2, cLine, nLines, &
            "2nd parameter for accretion rate: ", "(a,e9.3,1x,a)", 1.0, ok, .false.)
       call getReal("mdotpar3", MdotParameter3, cLine, nLines, &
            "3rd parameter for accretion rate: ", "(a,e9.3,1x,a)", 1.0, ok, .false.)
       call getReal("mdotpar4", MdotParameter4, cLine, nLines, &
            "4th parameter for accretion rate: ", "(a,e9.3,1x,a)", 1.0, ok, .false.)
       call getReal("mdotpar5", MdotParameter5, cLine, nLines, &
            "5th parameter for accretion rate: ", "(a,e9.3,1x,a)", 1.0, ok, .false.)
       call getReal("mdotpar6", MdotParameter6, cLine, nLines, &
            "6th parameter for accretion rate: ", "(a,e9.3,1x,a)", 1.0, ok, .false.)

       call getLogical("lte", lte, cLine, nLines, &
            "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("lycontthick", LyContThick, cLine, nLines, &
            "Optically thick Lyman continuum:","(a,1l,1x,a)", .false., ok, .false.)
       call getString("contflux", contFluxFile, cLine, nLines, &
            "Continuum flux filename (primary): ","(a,a,1x,a)","none", ok, .true.)
       call getString("popfile", popFilename, cLine, nLines, &
            "Grid populations filename: ","(a,a,1x,a)","none", ok, .true.)
       call getLogical("writepops", writePops, cLine, nLines, &
            "Write populations file: ","(a,1l,1x,a)", .true., ok, .true.)
       call getLogical("readpops", readPops, cLine, nLines, &
            "Read populations file: ","(a,1l,1x,a)", .true., ok, .true.)
       call getLogical("writephasepops", writePhasePops, cLine, nLines, &
            "Write populations file at each phase: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("readphasepops", readPhasePops, cLine, nLines, &
            "Read populations file (specific phase): ","(a,1l,1x,a)", .false., ok, .false.)
       !call getLogical("curtains", curtains, cLine, nLines, &
            !         "Curtains of accretion: ","(a,1l,1x,a)", .false., ok, .false.)
       call getReal("dipoleoffset", dipoleOffset, cLine, nLines, &
            "Dipole offset (degrees): ","(a,f7.1,1x,a)", 1.e14, ok, .true.)
       dipoleOffset = dipoleOffset * degToRad
       call getLogical("enhance", enhance, cLine, nLines, &
            "Accretion enhancement: ","(a,1l,1x,a)", .false., ok, .false.)
       call getInteger("nlower", nLower, cLine, nLines,"Lower level: ","(a,i2,a)",2,ok,.true.)
       call getInteger("nupper", nUpper, cLine, nLines,"Upper level: ","(a,i2,a)",3,ok,.true.)
       call getLogical("usehartmanntemp", useHartmannTemp, cLine, nLines, &
            "Use temperatures from Hartmann paper:","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("isotherm", isoTherm, cLine, nLines, &
            "Use isothermal temperature :","(a,1l,1x,a)", .false., ok, .false.)
       call getReal("isothermtemp", isoThermTemp, cLine, nLines, &
            "Isothermal temperature (K): ","(a,f7.1,1x,a)", 6500.0, ok, .false.)
       if (useHartmannTemp .and. isoTherm) then 
          if (writeoutput)  write(*,'(a)') "WARNING: useHartmannTemp and isoTherm6500 both specified!"
          stop
       end if
       if (useHartmannTemp) &
            call getReal("maxharttemp", maxHartTemp, cLine, nLines, &
            "Maximum of Hartmann temperature: ","(a,f7.1,1x,a)", 7436., ok, .false.)
       ! sub options for ttauri geometry
          if (ttau_discwind_on) then   ! commnted out here to make ttaur_turn_off_discwind to work
       ! --- parameters for ttauri wind
       call getDouble("DW_d", DW_d, cLine, nLines, &
            "Disc wind:: Wind soudce displacement [10^10cm]: ", &
            "(a,es9.3,1x,a)", 70.0d0, ok, .true.) 
       call getDouble("DW_Rmin", DW_Rmin, cLine, nLines, &
            "Disc wind:: Inner radius of the disc [10^10cm]: ", &
            "(a,es9.3,1x,a)", 70.0d0, ok, .true.) 
       call getDouble("DW_Rmax", DW_Rmax, cLine, nLines, &
            "Disc wind:: Outer radius of the disc [10^10cm]: ", &
            "(a,es9.3,1x,a)", 700.0d0, ok, .true.) 
       call getDouble("DW_Tmax", DW_Tmax, cLine, nLines, &
            "Disc wind:: Temperature of disc at inner radius [K]: ", &
            "(a,es9.3,1x,a)", 2000.0d0, ok, .true.) 
       call getDouble("DW_gamma", DW_gamma, cLine, nLines, &
            "Disc wind:: Exponent in the disc temperature power law [-]: ", &
            "(a,es9.3,1x,a)", -0.5d0, ok, .true.) 
       call getDouble("DW_Mdot", DW_Mdot, cLine, nLines, &
            "Disc wind:: Total mass-loss rate from disc [Msun/yr]: ", &
            "(a,es9.3,1x,a)", 1.0d-8, ok, .true.) 
       call getDouble("DW_alpha", DW_alpha, cLine, nLines, &
            "Disc wind:: Exponent in the mass-loss rate per unit area [-]: ", &
            "(a,es9.3,1x,a)", 0.5d0, ok, .true.) 
       call getDouble("DW_beta", DW_beta, cLine, nLines, &
            "Disc wind:: Exponent in the modefied beta-velocity law [-]: ", &
            "(a,es9.3,1x,a)", 0.5d0, ok, .true.) 
       call getDouble("DW_Rs", DW_Rs, cLine, nLines, &
            "Disc wind:: Effective accerelation length [10^10cm]: ", &
            "(a,es9.3,1x,a)", 50.0d0*DW_Rmin, ok, .true.) 
       call getDouble("DW_f", DW_f, cLine, nLines, &
            "Disc wind:: Scaling on the terminal velocity [-]: ", &
            "(a,es9.3,1x,a)", 2.0d0, ok, .true.) 
       call getDouble("DW_Twind", DW_Twind, cLine, nLines, &
            "Disc wind:: Isotherma temperature of disc wind [K]: ", &
            "(a,es9.3,1x,a)", 5000.0d0, ok, .true.) 
       endif
          if (ttau_jet_on) then  ! commented out here to make ttaur_turn_off_jet to work
       ! --- parameters for ttauri wind
       call getDouble("JET_Rmin", JET_Rmin, cLine, nLines, &
            "Minmium radius of Jet [10^10 cm]: ", &
            "(a,es9.3,1x,a)", TTauriRouter/1.0d10, ok, .false.) 
       call getDouble("JET_theta_j", JET_theta_j, cLine, nLines, &
            "TTauri jets:: [deg]  jet opening angle: ", &
            "(a,es9.3,1x,a)", 80.0d0, ok, .true.) 
       JET_theta_j = JET_theta_j * (Pi/180.0)  ! converting [deg] to [radians]

       call getDouble("JET_Mdot", JET_Mdot, cLine, nLines, &
            "TTauri jets:: [Msun/yr] mass loss rate in the jets: ", &
            "(a,es9.3,1x,a)", 1.0d-9, ok, .true.) 
       call getDouble("JET_a_param", JET_a_param, cLine, nLines, &
            "TTauri jets:: [-] a parameter in density function: ", &
            "(a,es9.3,1x,a)", 0.8d0, ok, .true.) 
       call getDouble("JET_b_param", JET_b_param, cLine, nLines, &
            "TTauri jets:: [-] b parameter in density function: ", &
            "(a,es9.3,1x,a)", 2.0d0, ok, .true.) 
       call getDouble("JET_Vbase", JET_Vbase, cLine, nLines, &
            "TTauri jets:: [km/s] Base velocity of jets: ", &
            "(a,es9.3,1x,a)", 20.0d0, ok, .true.) 
       call getDouble("JET_Vinf", JET_Vinf, cLine, nLines, &
            "TTauri jets:: [km/s] Terminal velocity of jets: ", &
            "(a,es9.3,1x,a)", 200.0d0, ok, .true.) 
       call getDouble("JET_beta", JET_beta, cLine, nLines, &
            "TTauri jets:: [-] a parameter in velocity function: ", &
            "(a,es9.3,1x,a)", 0.5d0, ok, .true.) 
       call getDouble("JET_gamma", JET_gamma, cLine, nLines, &
            "TTauri jets:: [-] a parameter in velocity function: ", &
            "(a,es9.3,1x,a)", 0.05d0, ok, .true.) 
       call getDouble("JET_T", JET_T, cLine, nLines, &
            "TTauri jets:: [K]  Isothermal temperature of jets: ", &
            "(a,es9.3,1x,a)", 1.0d4, ok, .true.) 
       endif
    endif

 if (geometry == "ttwind" .or. geometry == "donati") then
   call getLogical("lte", lte, cLine, nLines, &
            "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
   call getLogical("lycontthick", LyContThick, cLine, nLines, &
            "Optically thick Lyman continuum:","(a,1l,1x,a)", .false., ok, .false.)
   call getString("contflux", contFluxFile, cLine, nLines, &
        "Continuum flux filename (primary): ","(a,a,1x,a)","none", ok, .true.)
   call getString("popfile", popFilename, cLine, nLines, &
        "Grid populations filename: ","(a,a,1x,a)","none", ok, .true.)
   call getLogical("writepops", writePops, cLine, nLines, &
            "Write populations file: ","(a,1l,1x,a)", .true., ok, .true.)
   call getLogical("readpops", readPops, cLine, nLines, &
            "Read populations file: ","(a,1l,1x,a)", .true., ok, .true.)
    call getReal("dipoleoffset", dipoleOffset, cLine, nLines, &
         "Dipole offset (degrees): ","(a,f7.1,1x,a)", 1.e14, ok, .true.)
    dipoleOffset = dipoleOffset * degToRad
endif



 if (geometry == "luc_cir3d") then
      call getDouble("CIR_Rstar", CIR_Rstar, cLine, nLines, &
           "radius of central star  [R_sun] : ", "(a,es9.3,1x,a)", 1.0d0, ok, .true.) 
      rcore = (CIR_Rstar*rSol)/1.0e10       ! [10^10cm]
      call getDouble("CIR_Mass", CIR_Mass, cLine, nLines, &
           "Mass of the star [M_sun] : ", "(a,es9.3,1x,a)", 1.0d0, ok, .true.)       
      call getDouble("CIR_Twind", CIR_Twind, cLine, nLines, &
           "Isothemal temperature of the stellar wind [K] : ", &
           "(a,es9.3,1x,a)", 1.0d4, ok, .true.)       
      call getDouble("CIR_Mdot_scale", CIR_Mdot_scale, cLine, nLines, &
           "Scaling factor for CIR density and Mdot [-] : ", &
           "(a,es9.3,1x,a)", 1.0d0, ok, .true.)       
      call getLogical("lte", lte, cLine, nLines, &
           "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("lycontthick", LyContThick, cLine, nLines, &
           "Optically thick Lyman continuum:","(a,1l,1x,a)", .false., ok, .false.)
      call getString("contflux", contFluxFile, cLine, nLines, &
           "Continuum flux filename (primary): ","(a,a,1x,a)","none", ok, .true.)
      call getString("popfile", popFilename, cLine, nLines, &
           "Grid populations filename: ","(a,a,1x,a)","none", ok, .true.)
      call getLogical("writepops", writePops, cLine, nLines, &
           "Write populations file: ","(a,1l,1x,a)", .true., ok, .true.)
      call getLogical("readpops", readPops, cLine, nLines, &
           "Read populations file: ","(a,1l,1x,a)", .true., ok, .true.)
      call getLogical("writephasepops", writePhasePops, cLine, nLines, &
           "Write populations file at each phase: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("readphasepops", readPhasePops, cLine, nLines, &
           "Read populations file (specific phase): ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("enhance", enhance, cLine, nLines, &
           "Accretion enhancement: ","(a,1l,1x,a)", .false., ok, .false.)
      call getInteger("nlower", nLower, cLine, nLines,"Lower level: ",&
           "(a,i2,a)",2,ok,.true.)
      call getInteger("nupper", nUpper, cLine, nLines,"Upper level: ", &
           "(a,i2,a)",3,ok,.true.)
      call getString("coreprofile", intProFilename, cLine, nLines, &
           "Core profile: ","(a,a,1x,a)","none", ok, .false.)

 end if


 if (geometry == "cmfgen") then
      call getDouble("CMFGEN_Rmin", CMFGEN_Rmin, cLine, nLines, &
           "radius of central star  [10^10cm] : ", "(a,es9.3,1x,a)", 1.0d0, ok, .true.) 
      rcore = CMFGEN_Rmin      ! [10^10cm]
      call getLogical("lycontthick", LyContThick, cLine, nLines, &
           "Optically thick Lyman continuum:","(a,1l,1x,a)", .false., ok, .false.)
      call getString("contflux", contFluxFile, cLine, nLines, &
           "Continuum flux filename (primary): ","(a,a,1x,a)","none", ok, .true.)
      call getString("popfile", popFilename, cLine, nLines, &
           "Grid populations filename: ","(a,a,1x,a)","none", ok, .true.)
      call getLogical("writepops", writePops, cLine, nLines, &
           "Write populations file: ","(a,1l,1x,a)", .true., ok, .true.)
      call getLogical("readpops", readPops, cLine, nLines, &
           "Read populations file: ","(a,1l,1x,a)", .true., ok, .true.)
      call getLogical("writephasepops", writePhasePops, cLine, nLines, &
           "Write populations file at each phase: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("readphasepops", readPhasePops, cLine, nLines, &
           "Read populations file (specific phase): ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("enhance", enhance, cLine, nLines, &
           "Accretion enhancement: ","(a,1l,1x,a)", .false., ok, .false.)
 end if


 if (geometry == "romanova") then
      call getDouble("ROM_Rs", ROM_Rs, cLine, nLines, &
           "radius of central star  [R_sun] : ", "(a,es9.3,1x,a)", 1.0d0, ok, .true.) 
      rcore = (ROM_Rs*rSol)/1.0e10       ! [10^10cm]
      call getDouble("ROM_Mass", ROM_Mass, cLine, nLines, &
           "Mass of the star [M_sun] : ", "(a,es9.3,1x,a)", 1.0d0, ok, .true.)       
      call getLogical("ROM_isoT", ROM_isoT, cLine, nLines, &
           "Switch on and off isothermal option: ","(a,1l,1x,a)", .false., ok, .false.)
      call getDouble("ROM_T_flow", ROM_T_flow, cLine, nLines, &
           "Isothemal temperature of the flow [K] : ", &
           "(a,es9.3,1x,a)", 1.0d4, ok, .true.)       
      call getDouble("ROM_r_ref", ROM_r_ref, cLine, nLines, &
           "Reference length value [cm] : ", &
           "(a,es9.3,1x,a)", 1.0d0, ok, .true.)       
      call getDouble("ROM_rho_ref", ROM_rho_ref, cLine, nLines, &
           "Reference density value [g/cm^3] : ", &
           "(a,es9.3,1x,a)", 1.0d0, ok, .true.)       
      call getDouble("ROM_T_ref", ROM_T_ref, cLine, nLines, &
           "Reference temperature value [K] : ", &
           "(a,es9.3,1x,a)", 1.0d0, ok, .true.)       
      call getDouble("ROM_v_ref", ROM_v_ref, cLine, nLines, &
           "Reference speed value [cm/s] : ", &
           "(a,es9.3,1x,a)", 1.0d0, ok, .true.)       
      call getDouble("ROM_tilt", ROM_tilt, cLine, nLines, &
           "Tilt angle of magnetic axis [deg.] : ", &
           "(a,es9.3,1x,a)", 1.0d0, ok, .true.)       
      ROM_tilt = ROM_tilt*piDouble/180.0d0  ! [rad.]
      call getDouble("ROM_Period", ROM_Period, cLine, nLines, &
           "Rotational Period of star [days] : ", &
           "(a,es9.3,1x,a)", 1.0d0, ok, .true.)       
      ! converting it to seconds
      ROM_Period = ROM_Period*24.0d0*60.0d0*60.0d0   ! (sec)
      call getString("ROM_datafile", ROM_datafile, cLine, nLines, &
           "Filename of the data from Romanova: ","(a,a,1x,a)","none", ok, .true.)
      ! The followings are actually common for a line calculations.
      call getLogical("lte", lte, cLine, nLines, &
           "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("lycontthick", LyContThick, cLine, nLines, &
           "Optically thick Lyman continuum:","(a,1l,1x,a)", .false., ok, .false.)
      call getString("contflux", contFluxFile, cLine, nLines, &
           "Continuum flux filename (primary): ","(a,a,1x,a)","none", ok, .true.)
      call getString("popfile", popFilename, cLine, nLines, &
           "Grid populations filename: ","(a,a,1x,a)","none", ok, .true.)
      call getLogical("writepops", writePops, cLine, nLines, &
           "Write populations file: ","(a,1l,1x,a)", .true., ok, .true.)
      call getLogical("readpops", readPops, cLine, nLines, &
           "Read populations file: ","(a,1l,1x,a)", .true., ok, .true.)
      call getLogical("writephasepops", writePhasePops, cLine, nLines, &
           "Write populations file at each phase: ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("readphasepops", readPhasePops, cLine, nLines, &
           "Read populations file (specific phase): ","(a,1l,1x,a)", .false., ok, .false.)
      call getLogical("enhance", enhance, cLine, nLines, &
           "Accretion enhancement: ","(a,1l,1x,a)", .false., ok, .false.)
      call getInteger("nlower", nLower, cLine, nLines,"Lower level: ",&
           "(a,i2,a)",2,ok,.true.)
      call getInteger("nupper", nUpper, cLine, nLines,"Upper level: ", &
           "(a,i2,a)",3,ok,.true.)
      call getString("coreprofile", intProFilename, cLine, nLines, &
           "Core profile: ","(a,a,1x,a)","none", ok, .false.)
      ! for thin disc implementation
      call getReal("thindiskrin", ThinDiskRin, cLine, nLines, &
           "Thin disk inner radius  (10^10cm): ","(a,f7.1,1x,a)", 60.0, ok, .false.)
 end if


 
 if (coreEmissionLine) then
    
    call getReal("velfwhm", velWidthCoreEmissionLine, cLine, nLines, &
   "FWHM Velocity width of core emission line (km/s): ","(a,f7.0,1x,a)", 1000., ok, .true.)
    call getReal("relint", relintCoreEmissionLine, cLine, nLines, &
   "Relative intensity of core emission line: ","(a,f7.0,1x,a)", 10., ok, .true.)
    velWidthCoreEmissionLine = velWidthCoreEmissionLine * 1.e5
 endif


 if (geometry(1:8) == "windtest") then
   call getLogical("lte", lte, cLine, nLines, &
            "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
   call getLogical("lycontthick", LyContThick, cLine, nLines, &
            "Optically thick Lyman continuum:","(a,1l,1x,a)", .false., ok, .false.)
   call getString("contflux", contFluxFile, cLine, nLines, &
        "Continuum flux filename (primary): ","(a,a,1x,a)","none", ok, .true.)
   call getString("popfile", popFilename, cLine, nLines, &
        "Grid populations filename: ","(a,a,1x,a)","none", ok, .true.)
   call getLogical("writepops", writePops, cLine, nLines, &
            "Write populations file: ","(a,1l,1x,a)", .true., ok, .true.)
   call getLogical("readpops", readPops, cLine, nLines, &
            "Read populations file: ","(a,1l,1x,a)", .true., ok, .true.)
   call getInteger("nlower", nLower, cLine, nLines,"Lower level: ","(a,i2,a)",2,ok,.true.)
   call getInteger("nupper", nUpper, cLine, nLines,"Upper level: ","(a,i2,a)",3,ok,.true.)

 endif

 call getLogical("screened", screened, cLine, nLines, &
   "Source screened: ","(a,1l,a)", .false., ok, .false.)

 call getLogical("pencil", pencilBeam, cLine, nLines, &
   "Pencil beam illumination: ","(a,1l,a)", .false., ok, .false.)
 

 call getreal("teff", teff, cLine, nLines, &
   "Source temperature (K): ","(a,f7.0,1x,a)", 3000., ok, .false.)

 call getLogical("plezmodel", plezmodelOn, cLine, nLines, &
   "Use a Plez model atmosphere: ","(a,1l,a)", .false., ok, .false.)

 call getLogical("opaquecore", opaqueCore, cLine, nLines, &
   "Opaque Core: ","(a,1l,a)", .true., ok, .false.)

 call getLogical("thin_disc_on", thin_disc_on, cLine, nLines, &
   "Include geometrically thin Disc?: ", "(a,1l,1x,a)", .true., ok, .false.)

 call getLogical("tio", fillTio, cLine, nLines, &
   "TiO opacity: ","(a,1l,a)", .false., ok, .false.)

 call getLogical("gasopacity", includeGasOpacity, cLine, nLines, &
   "Include gas opacity: ","(a,1l,a)", .false., ok, .false.)

 if (fillTio .or. mie .or. useDust &
      .or. (geometry == "ttauri".and.mie) ) then
!      .or. (geometry == "ttauri" .and. ttau_disc_on) ) then

! amin and amax are left as microns here

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

       call getReal(grainFracLabel, grainFrac(i), cLine, nLines, &
            "Grain fractional abundance: ","(a,f8.5,1x,a)",1. , ok, .false.)
       grainFracTotal = grainFracTotal + grainFrac(i)

       call getReal(aminLabel, aMin(i), cLine, nLines, &
            "Min grain size (microns): ","(a,f8.5,1x,a)", 0.005, ok,  .true.)

       call getReal(amaxLabel, aMax(i), cLine, nLines, &
            "Max grain size (microns): ","(a,f10.5,1x,a)", 0.25, ok, .true.)

       call getReal(qDistLabel, qdist(i), cLine, nLines, &
            "Grain power law: ","(a,f4.1,1x,a)", 3.5, ok, .true. )

       call getReal(a0Label, a0(i), cLine, nLines, &
            "Scale length of grain size (microns): ","(a,f8.5,1x,a)", 1.0e20, ok, .false.)


       call getReal(pdistLabel, pdist(i), cLine, nLines, &
         "Exponent for exponential cut off: ","(a,f4.1,1x,a)", 1.0, ok, .false. )
       if (writeoutput) write(*,*)
    enddo

    if ((grainFracTotal .ne. 1.).and.(geometry .eq. "ppdisk")) then
       if (writeoutput) write (*,*) "Error: dust fractional abundances do not add up to 1."
       stop
    endif 

 endif


 call getLogical("disc_on", disc_on, cLine, nLines, &
   "Include accreation discs in SPH model: ","(a,1l,a)", .true., ok, .false.)
 
 call getInteger("idx_restrict_star", idx_restrict_star, cLine, nLines,  &
      "restrcting a calculation to this star : ","(a,i12,a)", 0, ok, .false.)

 call getLogical("rayleigh", fillRayleighOpacity, cLine, nLines, &
   "Rayleigh scattering opacity: ","(a,1l,a)", .false., ok, .false.)

 call getLogical("line", lineEmission, cLine, nLines, &
   "Line emission: ","(a,1l,a)", .false., ok, .false.)

 call getLogical("forceFirstScat", forceFirstScat, cLine, nLines, &
   "Forces the first scattering? : ","(a,1l,a)", .false., ok, .false.)

 if (lineEmission.or.doRaman.or.geometry.eq."planet") then

    call getReal("lamline", lamLine, cLine, nLines, &
         "Line emission wavelength: ","(a,f6.1,1x,a)", 6563., ok, .true.)

    if (lamStart < 0.) then
       lamStart = (1. + lamStart*1.e5/cSpeed)*lamLine
       lamEnd = (1. + lamEnd*1.e5/cSpeed)*lamLine
    endif

 endif

 call getLogical("thinline", thinLine, cLine, nLines, &
  "Optically thin line: ","(a,1l,a)",.false., ok, .false.)

 call getLogical("emissoff", lineOff, cLine, nLines, &
  "Line emission switched off: ","(a,1l,a)",.false., ok, .false.)

 call getLogical("staroff", starOff, cLine, nLines, &
  "Star emission switched off: ","(a,1l,a)",.false., ok, .false.)


 call getLogical("voigtprof", VoigtProf, cLine, nLines, &
         "Use Voigt profile: ","(a,1l,1x,a)", .false., ok, .false.)

 !
 ! Voigt profile prameters
 !
 call getReal("C_rad", C_rad, cLine, nLines, &
      "Damping constant (radiation)     in [A]: ","(a,1PE10.3,1x,a)", 0.0, ok, .false.)
 call getReal("C_vdw", C_vdw, cLine, nLines, &
      "Damping constant (van der Waals) in [A]: ","(a,1PE10.3,1x,a)", 0.0, ok, .false.)
 call getReal("C_stark", C_stark, cLine, nLines, &
      "Damping constant (Stark)         in [A]: ","(a,1PE10.3,1x,a)", 0.0, ok, .false.)



 call getLogical("thomson", fillThomson, cLine, nLines, &
   "Electron scattering: ","(a,1l,a)", .false., ok, .false.)

 call getLogical("usebias", usebias, cLine, nLines, &
   "Variance reduction: ","(a,1l,a)", .true., ok, .true.)

 call getLogical("interp", useInterp, cLine, nLines, &
   "Use opacity interpolation: ","(a,1l,a)", .true., ok, .false.)

 call getInteger("maxscat", maxScat, cLine, nLines, &
      "Max no of scatters: ","(a,i8,1x,a)", 1, ok, .false.)

 if (maxscat == 1) then
    if (writeoutput) write(*,*) "maxscat: USEAGE CHANGED, now limits maximum no of scatterings"
    if (writeoutput) write(*,*) "maxscat: You don't want to set this to 1."
    stop
 endif




 call getString("filename", outFile, cLine, nLines, &
  "Output spectrum filename: ","(a,a,1x,a)","spectrum.dat", ok, .true.)
 call replaceDots(outFile, done)
 if (done) then
    if (writeoutput) write(*,'(a,a)') "!!! Filename now: ",trim(outFile)
 endif
 

 call getLogical("ndf", useNdf, cLine, nLines, &
  "Use NDF format files: ","(a,1l,a)",.true., ok, .false.)

! Use plot_maps parameter instead (it already exists).
! call getLogical("doplots", doplots, cLine, nLines, &
!  "Plot outputs: ","(a,1l,a)",.true., ok, .false.)

 call getString("device", device, cLine, nLines, &
  "Plot device: ","(a,a,1x,a)","/xs", ok, .false.)

 call getLogical("plot_maps", plot_maps, cLine, nLines, &
  "Plot values on specified plane?: ","(a,1l,a)",.true., ok, .false.)

 call getString("plane_for_plot", plane_for_plot, cLine, nLines, &
  "Plane(s) for plotting: ","(a,a,1x,a)","x-y", ok, .false.)

 call getLogical("show_value_3rd_dim", show_value_3rd_dim, cLine, nLines, &
  "Display the third dimension on plot_AMR_value: ","(a,1l,a)",.false., ok, .false.)

 call getReal("zoomfactor", zoomFactor, cLine, nLines, &
      "Image zoom factor: ", "(a,1pe7.4,1x,a)", 0.01, ok, .false.)

 call getString("misc", misc, cLine, nLines, &
  "Miscallenous rubbish: ","(a,a,1x,a)","junk", ok, .false.)


 call getLogical("movie", movie, cLine, nLines, &
  "Make a movie: ","(a,1l,a)",.false., ok, .false.)

 call getLogical("stokesimage", stokesImage, cLine, nLines, &
      "Output stokes image: ","(a,1l,a)",.false., ok, .false.)

 call getLogical("narrowband", narrowBandImage, cLine, nLines, &
      "Narrow band image: ","(a,1l,a)",.false., ok, .false.)

 call getReal("imagesize", setImageSize, cLine, nLines, &
      "Image size (AU): ", "(a,1pe10.2,1x,a)", 0., ok, .false.)
 setimagesize = setimagesize * auTocm

 call getReal("imageScale", imageScale, cLine, nLines, &
      "Fraction of original image size (-): ", "(a,1pe10.2,1x,a)", 1.0, ok, .false.)



 if (stokesImage) then
    if (.not.narrowBandImage) then
       call getReal("vmin", vmin, cLine, nLines, &
            "Minimum velocity for image (km/s): ", "(a,1pe10.2,1x,a)", -20000., ok, .true.)
       call getReal("vmax", vmax, cLine, nLines, &
            "Maximum velocity for image (km/s): ", "(a,1pe10.2,1x,a)", 20000., ok, .true.)
    endif
  call getString("filter_set_name", filter_set_name, cLine, nLines, &
       "Name of filter set: ","(a,a,1x,a)","step_functions", ok, .false.)
 endif

 call getLogical("pvimage", doPVimage, cLine, nLines, &
      "Output pv image: ","(a,1l,a)",.false., ok, .false.)
 if (doPVimage) then
    call getInteger("nslit", nSlit, cLine, nLines, &
         "Number of slits: ","(a,i3,a)", 1, ok, .true.)
    call getInteger("pvnp", np, cLine, nLines, &
         "Size of PV image in p: ","(a,i3,a)", 1, ok, .true.)
    call getInteger("pvnv", nv, cLine, nLines, &
         "Size of PV image in v: ","(a,i3,a)", 1, ok, .true.)
    call getReal("slitposx1", slitPosition1%x, cLine, nLines, &
     "Slit start position x (arcsec): ", "(a,f7.1,1x,a)", 0., ok, .true.)
    call getReal("slitposy1", slitPosition1%y, cLine, nLines, &
     "Slit start position y (arcsec): ", "(a,f7.1,1x,a)", 0., ok, .true.)
    if (nSlit > 1) then
       call getReal("slitposx2", slitPosition2%x, cLine, nLines, &
            "Slit end position x (arcsec): ", "(a,f7.1,1x,a)", 0., ok, .true.)
       call getReal("slitposy2", slitPosition2%y, cLine, nLines, &
            "Slit end position y (arcsec): ", "(a,f7.1,1x,a)", 0., ok, .true.)
    endif
    call getReal("slitpa", slitPA, cLine, nLines, &
     "Slit position angle (degrees): ", "(a,f7.1,1x,a)", 0., ok, .true.)
    call getReal("slitlen", slitLength, cLine, nLines, &
     "Slit length (arcsec): ", "(a,f7.1,1x,a)", 0., ok, .true.)
    call getReal("slitwid", slitWidth, cLine, nLines, &
     "Slit width (arcsec): ", "(a,f7.1,1x,a)", 0., ok, .true.)
    call getReal("vmin", vmin, cLine, nLines, &
     "Minimum velocity for image (km/s): ", "(a,f7.1,1x,a)", -20000., ok, .true.)
    call getReal("vmax", vmax, cLine, nLines, &
     "Maximum velocity for image (km/s): ", "(a,f7.1,1x,a)", 20000., ok, .true.)

    call getReal("vsys", vsys, cLine, nLines, &
     "Systemic velocity (km/s): ", "(a,f7.1,1x,a)", 20000., ok, .true.)
    slitPA = slitPA * degToRad

    call getReal("vfwhm", vfwhm, cLine, nLines, &
         "Velocity FWHM of pv smoothing (km/s)","(a,f7.1,1x,a)",1.,ok,.true.)

    call getReal("pfwhm", pfwhm, cLine, nLines, &
         "Mean seeing (FWHM) for pv image (arcsec)","(a,f7.1,1x,a)",1.,ok,.true.)

    call getReal("distance", gridDistance, cLine, nLines, &
         "Grid distance (pc): ","(a,f4.1,1x,a)", 1., ok, .true.)
    gridDistance = gridDistance * pcTocm / 1.e10  ! now in 10^10 cm
 endif


 call getLogical("plotvel", plotVelocity, cLine, nLines, &
  "Plot velocity vectors: ","(a,1l,a)",.false., ok, .false.)

 call getLogical("sphericity", sphericitytest, cLine, nLines, &
  "Perform sphericity test: ","(a,1l,a)",.false., ok, .false.)

 call getReal("tauextra", tauExtra, cLine, nLines, &
  "Foreground optical depth : ","(a,f7.1,a)",0., ok, .false.)

 call getReal("tauext2", tauExtra2, cLine, nLines, &
  "Foreground optical depth : ","(a,f7.1,a)",0., ok, .false.)

 call getReal("taudiff", tauDiff, cLine, nLines, &
  "Mininum optical depth of cell to be in diffusion approx : ","(a,f7.1,a)",100., ok, .false.)

 call getReal("tauforce", tauForce, cLine, nLines, &
  "Forced optical depth of cell to be in diffusion approx : ","(a,f7.1,a)",100., ok, .false.)

 call getLogical("resetdiffusion", resetDiffusion, cLine, nLines, &
  "Reset diffusion zones to false if thin: ","(a,1l,a)",.true., ok, .false.)

 call getReal("edenstol", eDensTol, cLine, nLines, &
  "Fractional change in energy density for convergence: ","(a,f7.1,a)",0.01, ok, .false.) ! used for gauss-seidel sweep also


 if (geometry(1:9) .eq. "benchmark") then

   call getReal("rcore", rCore, cLine, nLines, &
       "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)

   call getReal("rinner", rInner, cLine, nLines, &
       "Inner Radius (AU): ","(a,f5.1,a)", 12., ok, .true.)

   call getReal("router", rOuter, cLine, nLines, &
       "Outer Radius (AU): ","(a,f8.2,a)", 20., ok, .true.)

   call getReal("height", height, cLine, nLines, &
       "Scale height (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

   call getReal("teff", teff, cLine, nLines, &
          "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)

   call getReal("rho", rho, cLine, nLines, &
          "Density: ","(a,e12.5,a)", 1., ok, .true.)

   rInner = rInner * auToCm / 1.e10
   rOuter = rOuter * auToCm / 1.e10
   height = height * autoCm / 1.e10
   rCore = rCore * rSol / 1.e10

endif

if (geometry == "whitney") then


   call getReal("rstellar", rStellar, cLine, nLines, &
       "Stellar radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)
   rStellar = rStellar * rSol
   
   call getReal("teff", teff, cLine, nLines, &
          "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)

   call getReal("mcore", mcore, cLine, nLines, &
       "Stellar mass (solar masses): ","(a,f5.1,a)", 10., ok, .true.)
   mcore = mcore * msol

   call getReal("mdotenv", mdotenv, cLine, nLines, &
       "Envelope infall rate (solar masses/year): ","(a,e10.3,a)", 20., ok, .true.)

   mDotEnv = mDotEnv * mSol / (365.25d0*24.d0*3600.d0)


   call getReal("erinner", erInner, cLine, nLines, &
       "Envelope inner radius (stellar radii): ","(a,f5.1,a)", 12., ok, .true.)

   erInner = erInner * rStellar

   call getReal("erouter", erOuter, cLine, nLines, &
       "Envelope outer radius (AU): ","(a,f5.1,a)", 20., ok, .true.)

   erOuter = erOuter * autocm

   call getReal("mdisc", mDisc, cLine, nLines, &
       "Disc mass (solar masses): ","(a,f5.3,a)", 1.e-4, ok, .true.)

   mDisc  = mDisc * mSol

   call getReal("drinner", drInner, cLine, nLines, &
       "Disc inner Radius (stellar radii): ","(a,f5.1,a)", 12., ok, .true.)

   drInner = drInner * rStellar

   call getReal("drouter", drOuter, cLine, nLines, &
       "Disc outer Radius (AU): ","(a,f5.1,a)", 20., ok, .true.)

   drOuter = drOuter * autocm


   call getReal("cavdens", cavdens, cLine, nLines, &
       "Cavity density (nh2/cm^3): ","(a,f5.1,a)", 10., ok, .true.)

   cavDens = 2. * cavDens * mHydrogen

   call getReal("cavangle", cavangle, cLine, nLines, &
       "Cavity opening angle (deg): ","(a,f5.1,a)", 10., ok, .true.)

   cavAngle = cavAngle * degtorad

endif

 if (geometry .eq. "shakara" .or. (geometry .eq. "iras04158")) then

    call getLogical("noscat", noScattering, cLine, nLines, &
         "No scattering opacity in model: ","(a,1l,1x,a)", .false., ok, .false.)

    call getReal("rcore", rCore, cLine, nLines, &
         "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)

    call getReal("rinner", rInner, cLine, nLines, &
         "Inner Radius (stellar radii): ","(a,f5.1,a)", 12., ok, .true.)

    call getReal("rsub", rSublimation, cLine, nLines, &
         "Dust sublimation radius (stellar radii): ","(a,f5.1,a)", 12., ok, .false.)

    call getReal("router", rOuter, cLine, nLines, &
         "Outer Radius (AU): ","(a,f5.1,a)", 20., ok, .true.)

    call getReal("height", height, cLine, nLines, &
         "Scale height (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

    call getReal("teff", teff, cLine, nLines, &
         "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)

    call getReal("rho", rho0, cLine, nLines, &
         "Density (g/cm^3): ","(a,f7.0,a)", 1., ok, .true.)

   call getReal("mcore", mCore, cLine, nLines, &
       "Core mass (solar masses): ","(a,f5.1,a)", 0.5, ok, .true.)

   call getReal("mdisc", mDisc, cLine, nLines, &
       "Disc mass (solar masses): ","(a,f5.3,a)", 1.e-4, ok, .true.)

   call getReal("alphadisc", alphaDisc, cLine, nLines, &
       "Disc alpha parameter: ","(a,f5.3,a)", 2.25, ok, .true.)

   call getReal("betadisc", betaDisc, cLine, nLines, &
       "Disc beta parameter: ","(a,f5.3,a)", 1.25, ok, .true.)

   call getString("contflux", contFluxFile, cLine, nLines, &
        "Continuum flux filename: ","(a,a,1x,a)","none", ok, .true.)

   call getLogical("vardustsub", variableDustSublimation, cLine, nLines, &
        "Variable dust sublimation temperature: ", "(a,1l,1x,a)", .false., ok, .true.)



   rCore = rCore * rSol / 1.e10
   rInner = rInner * rCore
   rSublimation = rSublimation * rCore
   rOuter = rOuter * autoCm / 1.e10
   height = height * autoCm / 1.e10
   mCore = mCore * mSol
   mDisc = mDisc * mSol

endif

if (geometry .eq. "planetgap") then

   call getLogical("planetgap", planetGap, cLine, nLines, &
        "Planet gap present: ", "(a,1l,1x,a)", .false., ok, .true.)

   call getReal("rcore", rCore, cLine, nLines, &
       "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)

   call getReal("rinner", rInner, cLine, nLines, &
       "Inner Radius (stellar radii): ","(a,f5.1,a)", 12., ok, .true.)

   call getReal("router", rOuter, cLine, nLines, &
       "Outer Radius (AU): ","(a,f5.1,a)", 20., ok, .true.)

   call getReal("height", height, cLine, nLines, &
       "Scale height (rcore): ","(a,1pe8.2,a)",1.e0,ok,.true.)

   call getReal("teff", teff, cLine, nLines, &
          "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)

   call getReal("mdisc", mDisc, cLine, nLines, &
       "Disc mass (solar masses): ","(a,f5.3,a)", 1.e-4, ok, .true.)

   call getReal("alphadisc", alphaDisc, cLine, nLines, &
       "Disc alpha parameter: ","(a,f5.3,a)", 2.25, ok, .true.)

   call getReal("betadisc", betaDisc, cLine, nLines, &
       "Disc beta parameter: ","(a,f5.3,a)", 1.25, ok, .true.)

   call getReal("mcore", mCore, cLine, nLines, &
       "Core mass (solar masses): ","(a,f5.1,a)", 0.5, ok, .true.)

   call getReal("rgap", rGap, cLine, nLines, &
       "Radius of gap (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

!   call getReal("mplanet", mPlanet, cLine, nLines, &
!       "Planet mass (Msol): ","(a,1pe8.5,a)",1.e-3,ok,.true.)
   call getReal("gapwidth", gapWidth, cLine, nLines, &
       "Gap width (AU): ","(a,1pe10.4,a)",0.1,ok,.true.)

   call getReal("gapalpha", gapViscAlpha, cLine, nLines, &
       "Alpha-parameter for gap viscosity: ","(a,1pe10.5,a)",1.e-4,ok,.true.)


   rCore = rCore * rSol / 1.e10
   rInner = rInner * rCore
   rOuter = rOuter * autoCm / 1.e10
   mDisc = mDisc * mSol
   mCore = mCore * mSol

endif

 if (geometry .eq. "warpeddisc") then

   call getLogical("noscat", noScattering, cLine, nLines, &
        "No scattering opacity in model: ","(a,1l,1x,a)", .false., ok, .false.)

   call getReal("rcore", rCore, cLine, nLines, &
       "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)

   call getReal("rinner", rInner, cLine, nLines, &
       "Inner Radius (stellar radii): ","(a,f5.1,a)", 12., ok, .true.)


   call getReal("router", rOuter, cLine, nLines, &
       "Outer Radius (AU): ","(a,f5.1,a)", 20., ok, .true.)

   call getReal("height", height, cLine, nLines, &
       "Scale height (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

   call getReal("warpradius", warpradius, cLine, nLines, &
       "Warp radius (inner radii): ","(a,f5.1,a)", 12., ok, .true.)

   call getReal("warpsigma", warpsigma, cLine, nLines, &
       "Warp sigma in radius (inner radii): ","(a,f5.1,a)", 12., ok, .true.)

   call getReal("warpheight", warpfracheight, cLine, nLines, &
       "Fractional height of disc (warp radius): ","(a,f5.1,a)", 12., ok, .true.)

   call getReal("warpangle", warpangle, cLine, nLines, &
       "PA of peak of warp (degrees): ","(a,f5.1,a)", 0., ok, .true.)

   call getReal("teff", teff, cLine, nLines, &
          "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)

   call getReal("mcore", mCore, cLine, nLines, &
       "Core mass (solar masses): ","(a,f5.1,a)", 0.5, ok, .true.)

   call getReal("mdisc", mDisc, cLine, nLines, &
       "Disc mass (solar masses): ","(a,f5.3,a)", 1.e-4, ok, .true.)

   call getReal("alphadisc", alphaDisc, cLine, nLines, &
       "Disc alpha parameter: ","(a,f5.3,a)", 2.25, ok, .true.)

   call getReal("betadisc", betaDisc, cLine, nLines, &
       "Disc beta parameter: ","(a,f5.3,a)", 1.25, ok, .true.)

   call getString("contflux", contFluxFile, cLine, nLines, &
        "Continuum flux filename: ","(a,a,1x,a)","none", ok, .true.)

   call getLogical("vardustsub", variableDustSublimation, cLine, nLines, &
        "Variable dust sublimation temperature: ", "(a,1l,1x,a)", .false., ok, .true.)

   call getLogical("warphydro", hydroWarp, cLine, nLines, &
        "Impose a warp on a completed hydro disc: ","(a,1l,1x,a)", .false., ok, .false.)

   rCore = rCore * rSol / 1.e10
   rInner = rInner * rCore
   rSublimation = rSublimation * rCore
   rOuter = rOuter * autoCm / 1.e10
   height = height * autoCm / 1.e10
   mCore = mCore * mSol
   mDisc = mDisc * mSol

   warpradius = warpradius * rinner
   warpsigma = warpsigma * rinner
   warpangle = warpangle * pi/180.

endif

if (geometry .eq. "ppdisk") then
   call getReal("rsmooth", rSmooth, cLine, nLines, &
       "Inner Smoothing Radius (AU): ","(a,f5.1,a)", 0.4, ok, .true.)

   call getReal("height", height, cLine, nLines, &
       "Scale height (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

   call getReal("rheight", rHeight, cLine, nLines, &
       "Scale height measured at (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

   call getReal("rgap", rGap, cLine, nLines, &
       "Gap is at (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

   call getReal("mplanet", mPlanet, cLine, nLines, &
       "Planet mass (Msol): ","(a,1pe8.5,a)",1.e-3,ok,.true.)

   call getReal("gapalpha", gapViscAlpha, cLine, nLines, &
       "alpha-parameter for gap viscosity: ","(a,1pe8.5,a)",1.e-4,ok,.true.)

   call getReal("rho", rho0, cLine, nLines, &
        "Density scaling factor: ","(a,f8.5,a)", 1., ok, .true.)

   call getReal("flaringpower", flaringPower, cLine, nLines, &
        "Disc flaring power: ","(a,f8.5,a)", 1.25, ok, .true.)

   call getReal("sigmapower", sigmaPower, cLine, nLines, &
        "Disc surface density power: ","(a,f8.5,a)", 0.5, ok, .true.)

   call getReal("rinner", rInner, cLine, nLines, &
       "Inner Radius (10^10 cm): ","(a,f8.2,a)", 1.5E2, ok, .true.)

   call getReal("router", rOuter, cLine, nLines, &
       "Outer Radius (10^10 cm): ","(a,f8.2,a)", 6.0E3, ok, .true.)

   if (solveVerticalHydro) then
      call getReal("mdisc", mDisc, cLine, nLines, &
          "Disc mass (solar masses): ","(a,1pe8.2,a)", 1.e-4, ok, .true.)

      call getReal("mcore", mCore, cLine, nLines, &
          "Core mass (solar masses): ","(a,f5.1,a)", 0.5, ok, .true.)

      mCore = mCore * mSol
      mDisc = mDisc * mSol
   end if

   call getReal("teff", Teff, cLine, nLines, &
        "Effective temp (K): ","(a,f7.0,a)", 5780., ok, .true.)

   call getReal("mbol", Mbol, cLine, nLines, &
        "Bolometric magnitude (mags): ","(a,f7.3,a)", 4.75, ok, .true.)

   call getReal("rcore", rcore, cLine, nLines, &
        "Core radius (solar radii): ","(a,f5.1,a)", 0., ok, .false.)

   rCore = rCore * rSol

! This is not ideal having to calculate the stellar radius here as well as in
! torusMain - if the assumed bolometric magnitude of the sun changes, for
! example, we must change code in both places
   if (rCore .eq. 0.) then
      rCore = sqrt((10.**((4.64-Mbol)/2.5)*lSol)/(4.*pi*stefanBoltz*Teff**4))
   end if

   if (rSmooth .eq. 0.) then
!       source(1)%luminosity = 10.**((4.64-Mbol)/2.5)*lSol
!       source(1)%radius = sqrt(source(1)%luminosity/(4.*pi*stefanBoltz*Teff**4))/1.d10
!     rSmooth = 10 stellar radii if nothing else is specified
!      rSmooth = 10.*sqrt((10.**((4.64-Mbol)/2.5)*lSol)/(4.*pi*stefanBoltz*Teff**4))/auToCm
      rSmooth = 10.*rCore/auToCm
   end if

   ! rCore needs to be in units of 10^10 cm and this is not done elsewhere in this case
   rCore = rCore / 1.e10

   if (rInner .eq. 0.) then
      rInner = 0.5 * rSmooth * auToCm / 1.e10
   end if
endif



 if (geometry .eq. "clumpydisc") then

    call getInteger("nclumps", nClumps, cLine, nLines,"Number of clumps: ", &
         & "(a,i2,a)",1000,ok,.true.)

   call getReal("rcore", rCore, cLine, nLines, &
       "Core radius (solar radii): ","(a,f5.1,a)", 10., ok, .true.)

   call getReal("rinner", rInner, cLine, nLines, &
       "Inner Radius (stellar radii): ","(a,f5.1,a)", 12., ok, .true.)

   call getReal("router", rOuter, cLine, nLines, &
       "Outer Radius (AU): ","(a,f5.1,a)", 20., ok, .true.)

   call getReal("height", height, cLine, nLines, &
       "Scale height (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

   call getReal("teff", teff, cLine, nLines, &
          "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)

   call getReal("rho", rho0, cLine, nLines, &
          "Density (g/cm^3): ","(a,f7.0,a)", 1., ok, .true.)

   call getReal("mcore", mCore, cLine, nLines, &
       "Core mass (solar masses): ","(a,f5.1,a)", 0.5, ok, .true.)

   call getReal("mdisc", mDisc, cLine, nLines, &
       "Disc mass (solar masses): ","(a,f5.1,a)", 1.e-4, ok, .true.)

   rCore = rCore * rSol / 1.e10
   rInner = rInner * rCore
   rOuter = rOuter * autoCm / 1.e10
   height = height * autoCm / 1.e10
   mCore = mCore * mSol
   mDisc = mDisc * mSol

endif

 if (geometry .eq. "aksco") then

   call getReal("rinner", rInner, cLine, nLines, &
       "Inner Radius (AU): ","(a,f5.1,a)", 12., ok, .true.)

   call getReal("router", rOuter, cLine, nLines, &
       "Outer Radius (AU): ","(a,f5.1,a)", 20., ok, .true.)

   call getReal("height", height, cLine, nLines, &
       "Scale height (AU): ","(a,1pe8.2,a)",1.e0,ok,.true.)

   call getReal("rho", rho0, cLine, nLines, &
          "Density (g/cm^3): ","(a,f7.0,a)", 1., ok, .true.)

   rInner = rInner * autoCm / 1.e10
   rOuter = rOuter * autoCm / 1.e10
   height = height * autoCm / 1.e10
   solveVerticalHydro = .false.
   mCore = 0.5 * mSol 
   sigma0 = 1.e-3 ! at one AU

endif

 ! For bipolar jets geometry

 if (geometry(1:4) ==  "jets") then
    call getLogical("lte", lte, cLine, nLines, &
         "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
    call getLogical("lycontthick", LyContThick, cLine, nLines, &
            "Optically thick Lyman continuum","(a,1l,1x,a)", .false., ok, .false.)
    call getString("contflux", contFluxFile, cLine, nLines, &
         "Continuum flux filename (primary): ","(a,a,1x,a)","none", ok, .true.)
    call getString("popfile", popFilename, cLine, nLines, &
         "Grid populations filename: ","(a,a,1x,a)","none", ok, .true.)
    call getLogical("writepops", writePops, cLine, nLines, &
         "Write populations file: ","(a,1l,1x,a)", .true., ok, .true.)
    call getLogical("readpops", readPops, cLine, nLines, &
         "Read populations file: ","(a,1l,1x,a)", .true., ok, .true.)
    call getLogical("curtains", curtains, cLine, nLines, &
            "Curtains of accretion: ","(a,1l,1x,a)", .false., ok, .false.)
    call getReal("dipoleoffset", dipoleOffset, cLine, nLines, &
         "Dipole offset (degrees): ","(a,f7.1,1x,a)", 1.e14, ok, .true.)
    dipoleOffset = dipoleOffset * degToRad

    call getString("ion_name", ion_name, cLine, nLines, &
         "Name of the ion for a line: ","(a,a,1x,a)","H I", ok, .true.)
    call getReal("ion_frac", ion_frac, cLine, nLines, &
          "Inonization fraction of the ion above: ","(a,f7.0,a)", 1.0, ok, .true.)    
    
    call getInteger("nlower", nLower, cLine, nLines,"Lower level: ", &
         & "(a,i2,a)",2,ok,.true.)
    
    call getInteger("nupper", nUpper, cLine, nLines,"Upper level: ", &
         & "(a,i2,a)",3,ok,.true.)

    call getReal("Rmin_bp", Rmin_bp, cLine, nLines, &
         "Radius of central star  [10^10cm] : ","(a,f5.1,a)", 1.0, ok, .true.)
    call getReal("Rmax_bp", Rmax_bp, cLine, nLines, &
         "Cutoff radius in [10^10 cm] : ","(a,1PE7.1,a)", 10.0, ok, .true.)
    call getReal("Vo_bp", Vo_bp, cLine, nLines, &
         "A small offset in V in  [km/s] : ","(a,f5.1,a)", 5.0, ok, .true.)

    call getReal("Vinf_jets", Vinf_jets, cLine, nLines, &
         "Terminal velocity in  [kms] : ","(a,f5.1,a)", 250.0, ok, .true.)
    call getReal("beta_jets", beta_jets, cLine, nLines, &
         "Beta in the beta-belocity law [-] : ","(a,f5.1,a)", 0.5, ok, .true.)
    call getReal("Mdot_jets", Mdot_jets, cLine, nLines, &
         "Mass loss rate in jets  [M_solar/yr] : ","(a,1PE7.1,a)", 1.0e-9, ok, .true.)
    call getReal("theta_o_jets", theta_o_jets, cLine, nLines, &
         "Half opening angle [degrees] : ","(a,f5.1,a)", 30.0, ok, .true.)
    call getReal("Tcore_jets", Tcore_jets, cLine, nLines, &
         "Temperature at the core at Rmin in [10^4 K] : ","(a,f5.1,a)", 5.0, ok, .true.)
    call getReal("e6_jets", e6_jets, cLine, nLines, &
         "Exponet in tempereture eq. [-] : ","(a,f5.1,a)", 1.0, ok, .true.)


    call getReal("Vinf_disk", Vinf_disk, cLine, nLines, &
         "Terminal velocity in  [kms] : ","(a,f5.1,a)", 250.0, ok, .true.)
    call getReal("beta_disk", beta_disk, cLine, nLines, &
         "Beta in the beta-belocity law [-] : ","(a,f5.1,a)", 0.5, ok, .true.)
    call getReal("Mdot_disk", Mdot_disk, cLine, nLines, &
         "Mass loss rate in jets  [M_solar/yr] : ","(a,1PE7.1,a)", 1.0e-9, ok, .true.)
    call getReal("theta_o_disk", theta_o_disk, cLine, nLines, &
         "Half opening angle [degrees] : ","(a,f5.1,a)", 30.0, ok, .true.)
    call getReal("Tcore_disk", Tcore_disk, cLine, nLines, &
         "Temperature at the core at Rmin in [10^4 K] : ","(a,f5.1,a)", 5.0, ok, .true.)
    call getReal("e6_disk", e6_disk, cLine, nLines, &
         "Exponet in tempereture eq. [-] : ","(a,f5.1,a)", 1.0, ok, .true.)

    call getReal("Rdisk_min", Rdisk_min, cLine, nLines, &
         "The minimum radius of the disk. [10^10 cm] : ","(a,f5.1,a)", &
         1.0, ok, .true.)
    call getReal("Rdisk_max", Rdisk_max, cLine, nLines, &
         "The maximum radius of the disk. [10^10 cm] : ","(a,f5.1,a)", &
         100.0, ok, .true.)
    call getReal("h_disk", h_disk, cLine, nLines, &
         "Thickness of the disk [10^10 cm] : ","(a,f5.1,a)", &
         2.0, ok, .true.)
    call getReal("rho_scale", rho_scale, cLine, nLines, &
         "The density in the units of rho max for disk wind : ","(a,1PE7.1,a)", &
         1.0e-5, ok, .true.)



    
    ! save the variables in an object in jets_mod module.
    call set_jets_parameters(Rmin_bp, Rmax_bp, Vo_bp, &
         & Vinf_jets, beta_jets, Mdot_jets, theta_o_jets, Tcore_jets, e6_jets, &
         & Vinf_disk, beta_disk, Mdot_disk, theta_o_disk, Tcore_disk, e6_disk, &
         & Rdisk_min, Rdisk_max, h_disk, rho_scale)
    
 end if
 

    call getReal("h_abund", h_abund, cLine, nLines, &
         "Hydrogen abdunance: ","(a,1PF8.3,a)", &
         1., ok, .false.)

    call getReal("he_abund", he_abund, cLine, nLines, &
         "Helium abdunance: ","(a,1PF8.3,a)", &
         0.1, ok, .false.)

    call getReal("c_abund", c_abund, cLine, nLines, &
         "Carbon abdunance: ","(a,1PF8.3,a)", &
         22.e-5, ok, .false.)

    call getReal("n_abund", n_abund, cLine, nLines, &
         "Nitrogen abdunance: ","(a,1PF8.3,a)", &
         4.e-5, ok, .false.)

    call getReal("o_abund", o_abund, cLine, nLines, &
         "Oxygen abdunance: ","(a,1PF8.3,a)", &
         33.e-5, ok, .false.)

    call getReal("ne_abund", ne_abund, cLine, nLines, &
         "Neon abdunance: ","(a,1PF8.3,a)", &
         5.e-5, ok, .false.)

    call getReal("s_abund", s_abund, cLine, nLines, &
         "Sulphur abdunance: ","(a,1PF8.3,a)", &
         0.9e-5, ok, .false.)

 
 if (writeoutput) write(*,*) " "
 
! Close input file
close(32)

end subroutine inputs

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

 subroutine getReal(name, rval, cLine, nLines, message, format, rdef, ok, &
                    musthave)
  character(len=*) :: name
  real :: rval
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
    write(output,format) trim(message),rval,default
    call writeInfo(output, TRIVIAL)
 endif
 end subroutine getReal
 

 subroutine getDouble(name, dval, cLine, nLines, message, format, ddef, ok, &
                    musthave)
  character(len=*) :: name
  real(double) :: dval
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
 end subroutine getDouble


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
  write(output,format) trim(message),trim(rval),default
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
  character(len=6) :: trueOrFalse, tmp
  logical :: ok, thisIsDefault
  ok = .true.
  default = " "
  thisIsDefault = .false.
  tmp = cformat(1:6)
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

 subroutine getRealArray(name, rval, cLine, nLines, message, format, rdef, ok, &
                      musthave)
  character(len=*) :: name
  real :: rval(:)
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
  call findRealArray(name, rval, cLine, nLines, ok)
  if (.not. ok) then
    if (musthave) then
       if (writeoutput)  write(*,'(a,a)') name, " must be defined"
       stop
    endif
    rval(:) = rdef
    default = " (default)"
  endif
  write(output,*) trim(message),rval,default
  call writeInfo(output, TRIVIAL)
 end subroutine getRealArray

end module inputs_mod
