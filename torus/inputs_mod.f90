module inputs_mod

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
  character(len=80) :: cLine(100) 
  character(len=10) :: default


  logical :: done

  ! character vars for unix environ

  character(len=80) :: dataDirectory



  contrast = 1.
  grainSize = 1.

  nBlobs = 0
  nLines = 0

  inputOK = .true.

  do
     nLines = nLines + 1
     read(*,'(a80)',end=10) cLine(nLines)
     if (trim(cLine(nLines)(1:1)) == "%") nLines = nLines - 1   ! % is a comment
  end do
10 continue
  nLines = nLines - 1

  errno = 0
  call unixGetenv("TORUS_DATA",dataDirectory,i,errno)
  if (errNo /= 0) then
     write(*,'(A)') "! warning TORUS_DATA environment variable should be defined"
     stop
  endif

  write(*,*) " "
  write(*,'(a)') "TORUS: spectropolarimetric scattering model"
  write(*,'(a)') "-------------------------------------------"
  write(*,*) " "


  write(*,'(a)') "Input parameters"
  write(*,'(a)') "----------------"
  write(*,*) " "


  call getInteger("nphotons", nPhotons, cLine, nLines, &
       "Number of photons: ", "(a,i8,1x,a)", 100000, ok, .true.)

 call getReal("lamstart", lamstart, cLine, nLines, &
     "X-array start (angs or kms): ", "(a,f7.1,1x,a)", 2000., ok, .true.)

 call getReal("lamend", lamend, cLine, nLines, &
     "X-array end (angs or kms): ", "(a,f7.1,1x,a)", 10000., ok, .true.)


 call getLogical("wavlin", lamLinear, cLine, nLines, &
         "Linear wavelength array: ","(a,1l,1x,a)", .true., ok, .true.)

 

 call getInteger("nlambda", nlambda, cLine, nLines, &
     "Number of wavelength/velocity bins: ", "(a,i4,1x,a)", 20, ok, .true.)


 if (mie) then
    call getLogical("lucy", lucyRadiativeEq, cLine, nLines, &
         "Lucy radiative equ.: ","(a,1l,1x,a)", .false., ok, .false.)
 endif


  call getString("distortion", distortionType, cLine, nLines, &
       "Distortion type: ","(a,a,1x,a)","none", ok, .false.)

  call getReal("inclination", inclination, cLine, nLines, &
       "Inclination: ","(a,f4.1,1x,a)", 90., ok, .true.)


  call getReal("scale", scale, cLine, nLines, &
       "Scale (rsolar): ","(a,f6.1,1x,a)", 1000., ok, .false.) 


  call getString("gridtype", gridcoords, cLine, nLines, &
       "Grid type: ","(a,a,a)","cartesian",ok, .true.)

  call getLogical("gridusesamr", gridUsesAMR, cLine, nLines, &
       "Grid uses adaptive mesh refinement: ","(a,1l,1x,a)", .false., ok, .false.)
  
  if (gridUsesAMR) then
     call getReal("amrgridsize", amrGridSize, cLine, nLines, &
          "Size of adaptive mesh grid: ","(a,f6.1,1x,a)", 1000., ok, .true.) 
     call getReal("amrgridcentrex", amrGridCentreX, cLine, nLines, &
          "Grid centre X-coordinate: ","(a,f6.1,1x,a)", 0., ok, .false.) 
        if (amrGridCentreX == 0.0) then 
           print *, 'WARNING: amrGridCentreX == 0. This may cause numerical problems!'
        end if
     call getReal("amrgridcentrey", amrGridCentreY, cLine, nLines, &
          "Grid centre Y-coordinate: ","(a,f6.1,1x,a)", 0., ok, .false.) 
        if (amrGridCentreY == 0.0) then 
           print *, 'WARNING: amrGridCentreY == 0. This may cause numerical problems!'
        end if
     call getReal("amrgridcentrez", amrGridCentreZ, cLine, nLines, &
          "Grid centre Z-coordinate: ","(a,f6.1,1x,a)", 0., ok, .false.) 
        if (amrGridCentreZ == 0.0) then 
           print *, 'WARNING: amrGridCentreZ == 0. This may cause numerical problems!'
        end if
     call getReal("limitscalar", limitScalar, cLine, nLines, &
          "Scalar limit for subcell division: ","(a,es9.3,1x,a)", 1000., ok, .true.) 
     call getLogical("dosmoothgrid", doSmoothGrid, cLine, nLines, &
          "Smooth AMR grid: ","(a,1l,1x,a)", .false., ok, .false.)
     if (doSmoothGrid) then
       call getReal("smoothfactor", smoothFactor, cLine, nLines, &
            "Inter-cell maximum ratio before smooth: ","(a,f6.1,1x,a)", 5., ok, .false.)
     else
       smoothFactor = 0.0      
     end if
     call getReal("samplefreq", sampleFreq, cLine, nLines, &
          "Max samples per AMR subcell: ","(a,f6.1,1x,a)", 2., ok, .true.) 
     call getString("popfile", popFilename, cLine, nLines, &
          "Grid populations filename: ","(a,a,1x,a)","none", ok, .true.)
     call getLogical("writepops", writePops, cLine, nLines, &
          "Write populations file: ","(a,1l,1x,a)", .true., ok, .true.)
     call getLogical("readpops", readPops, cLine, nLines, &
          "Read populations file: ","(a,1l,1x,a)", .true., ok, .true.)
     if (readPops .and. writePops) &
        write(*,'(a)') "WARNING: both readPops and writePops set."
     if (readPops) &
        call getLogical("readfileformatted", readFileFormatted, cLine, nLines, &
             "Populations input file is formatted: ","(a,1l,1x,a)", .false., ok, .false.)
     if (writePops) &
        call getLogical("writefileformatted", writeFileFormatted, cLine, nLines, &
             "Populations output file will be formatted: ","(a,1l,1x,a)", .false., ok, .false.)
      
     
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


  if (geometry(1:5) .eq. "torus") then
     call getReal("rtorus", rTorus, cLine, nLines, &
          "Radius of torus (AU): ","(a,f5.1,a)", 0.7, ok, .true.)
     call getReal("router", rOuter, cLine, nLines, &
          "Radius of torus x-section (AU): ","(a,f5.1,a)", 0.1, ok, .true.)
     call getReal("teff", teff, cLine, nLines, &
          "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)
     rTorus = rTorus * auTocm
     rOuter = rOuter * auTocm
  endif

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

  if (geometry .eq. "wr104") then
     call getReal("teff", teff, cLine, nLines, &
          "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)
     call getReal("rho", rho, cLine, nLines, &
          "Density: ","(a,f7.2,a)", 1., ok, .true.)
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

  if (geometry .eq. "resonance") then
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

  if ((geometry .eq. "puls").or.(geometry .eq. "wind")) then
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
     call getString("coreprofile", intProFilename, cLine, nLines, &
          "Core profile: ","(a,a,1x,a)","none", ok, .false.)
     call getReal("vrot", vRot, cLine, nLines, &
          "Rotational velocity (km/s): ","(a,f5.1,a)", 0., ok, .true.)
     call getLogical("lte", lte, cLine, nLines, &
          "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
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
   call getReal("mass1", mass1, cLine, nLines, &
        "Primary mass (msol): ","(a,f5.1,a)", 1., ok, .true.)
   call getReal("radius1", radius1, cLine, nLines, &
        "Primary radius (rsol): ","(a,f5.1,a)", 1., ok, .true.)
   call getReal("mdot1", mdot1, cLine, nLines, &
        "Primary log mDot (msol/yr): ","(a,f6.3,a)", 1., ok, .true.)
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
        "Secondary log mDot (msol/yr): ","(a,f6.3,a)", 1., ok, .true.)
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
   mdot1 = (10.**mdot1)*msol/(365.25*24.*3600.)
   mdot2 = (10.**mdot2)*msol/(365.25*24.*3600.)
   deflectionAngle = deflectionAngle * degTorad

endif

 call getReal("rho", rho, cLine, nLines, &
      "Density (xxx): ","(a,1p,e10.2,1p,1x,a)", 0., ok, .false.)

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

 call getReal("phaseoffset", phaseOffset, cLine, nLines, &
      "Phase offset: ","(a,f3.1,a)", 0., ok, .false.)


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
     "rMajor (AU): ", "(a,f5.1,a)", 0.5, ok, .true.)
   call getReal("rmin", rmin, cLine, nLines, &
     "rMinor (AU): ", "(a,f5.1,a)", 0.5, ok, .true.)
    call getReal("rho", rho, cLine, nLines, &
    "Density (xxx): ","(a,1p,e10.2,1p,1x,a)", 1.e-6, ok, .true.)
    call getReal("teff", teff, cLine, nLines, &
         "Effective temp (K): ","(a,f7.0,a)", 1., ok, .true.)
    rmaj = rmaj * AUtocm / 1.e10
    rmin = rmin * AUtocm / 1.e10
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






 if (geometry(1:4) .eq. "star") then 
   call getReal("radius", radius, cLine, nLines, &
       "radius: ","(a,f5.1,a)", 100., ok, .true.)

   call getReal("mdot", mdot, cLine, nLines, &
       "Mdot (Msol/yr): ","(a,1pe8.2,a)",1.e-6,ok,.true.)


   call getReal("velocity", vterm, cLine, nLines, &
       "Wind velocity (km/s): ","(a,1pe8.2,a)",1.e1,ok,.true.)

   call getReal("kfac", kfac, cLine, nLines, &
       "kfac: ", "(a,f5.1,a)", 0.5, ok, .true.)


 endif


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
       "Dipole offset (degrees): ","(a,f7.1,1x,a)", 1.e14, ok, .true.)
   dipoleOffset = dipoleOffset * degToRad
 endif

 if (distortionType .eq. "rotation") then 
   call getReal("vrot", vRot, cLine, nLines, &
       "Rotational velocity (km/s): ","(a,f5.1,a)", 100., ok, .true.)
   vRot = vRot * 1.e5
 endif


 call findLogical("mie", mie, cLine, nLines, ok)
 if (.not. ok) then
    mie = .true.
    default = " (default)"
 endif
 if (mie) then
   write(*,'(a)') "Scattering phase matrix: Mie"
  else
   write(*,'(a)') "Scattering phase matrix: Rayleigh"
 endif


 if (mie) then
    call getLogical("lucy", lucyRadiativeEq, cLine, nLines, &
         "Lucy radiative equ.: ","(a,1l,1x,a)", .false., ok, .true.)
 endif



 if (mie) then
     call getString("graintype", grainType, cLine, nLines, &
  "Grain type: ","(a,a,1x,a)","silicate", ok, .true.)
 endif

 if (.not.mie) then
 default = " "
 call findReal("kappasca", inputKappaSca, cLine, nLines, ok)
 if (ok) then
   write(*,'(a,1pe12.4)') "Scattering cross-section: ",inputKappaSca
 endif
 default = " "
 call findReal("kappaabs", inputKappaAbs, cLine, nLines, ok)
 if (ok) then
   write(*,'(a,1pe12.4)') "Absorption cross-section: ",inputKappaAbs
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

 if (geometry == "ttauri") then
   call getLogical("lte", lte, cLine, nLines, &
            "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
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
   call getLogical("enhance", enhance, cLine, nLines, &
            "Accretion enhancement: ","(a,1l,1x,a)", .false., ok, .false.)
   call getInteger("nlower", nLower, cLine, nLines,"Lower level: ","(a,i2,a)",2,ok,.true.)
   call getInteger("nupper", nUpper, cLine, nLines,"Upper level: ","(a,i2,a)",3,ok,.true.)
endif

 if (geometry == "ttwind" .or. geometry == "donati") then
   call getLogical("lte", lte, cLine, nLines, &
            "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
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

 
 if (coreEmissionLine) then
    
    call getReal("velfwhm", velWidthCoreEmissionLine, cLine, nLines, &
   "FWHM Velocity width of core emission line (km/s): ","(a,f7.0,1x,a)", 1000., ok, .true.)
    call getReal("relint", relintCoreEmissionLine, cLine, nLines, &
   "Relative intensity of core emission line: ","(a,f7.0,1x,a)", 10., ok, .true.)
    velWidthCoreEmissionLine = velWidthCoreEmissionLine * 1.e5
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

 call getLogical("tio", fillTio, cLine, nLines, &
   "TiO opacity: ","(a,1l,a)", .false., ok, .false.)

 if (fillTio.or.mie) then
    call getReal("amin", aMin, cLine, nLines, &
         "Min grain size (microns): ","(a,f8.5,1x,a)", 0.0005, ok,  .false.)
    aMin = aMin * 1.e-4

    call getReal("amax", aMax, cLine, nLines, &
         "Max grain size (microns): ","(a,f8.5,1x,a)", 0.25, ok, .false.)
    aMax = aMax * 1.e-4

    call getReal("qdist", qdist, cLine, nLines, &
         "Grain power law: ","(a,f4.1,1x,a)", -3.5, ok, .false. )
 endif


 call getLogical("rayleigh", fillRayleighOpacity, cLine, nLines, &
   "Rayleigh scattering opacity: ","(a,1l,a)", .false., ok, .false.)

 call getLogical("line", lineEmission, cLine, nLines, &
   "Line emission: ","(a,1l,a)", .false., ok, .false.)

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


 call getLogical("thomson", fillThomson, cLine, nLines, &
   "Electron scattering: ","(a,1l,a)", .false., ok, .false.)

 call getLogical("usebias", usebias, cLine, nLines, &
   "Variance reduction: ","(a,1l,a)", .true., ok, .false.)

 call getLogical("interp", useInterp, cLine, nLines, &
   "Use opacity interpolation: ","(a,1l,a)", .true., ok, .false.)

 call getInteger("maxscat", maxScat, cLine, nLines, &
      "Max no of scatters: ","(a,i4,1x,a)", 1, ok, .false.)





 call getString("filename", outFile, cLine, nLines, &
  "Output spectrum filename: ","(a,a,1x,a)","spectrum.dat", ok, .true.)
 call replaceDots(outFile, done)
 if (done) then
    write(*,'(a,a)') "!!! Filename now: ",trim(outFile)
 endif
 

 call getLogical("ndf", useNdf, cLine, nLines, &
  "Use NDF format files: ","(a,1l,a)",.true., ok, .false.)

 call getString("device", device, cLine, nLines, &
  "Plot device: ","(a,a,1x,a)","/xs", ok, .false.)

 call getString("misc", misc, cLine, nLines, &
  "Miscallenous rubbish: ","(a,a,1x,a)","junk", ok, .false.)


 call getLogical("movie", movie, cLine, nLines, &
  "Make a movie: ","(a,1l,a)",.false., ok, .false.)

 call getLogical("stokesimage", stokesImage, cLine, nLines, &
      "Output stokes image: ","(a,1l,a)",.false., ok, .false.)

 call getLogical("narrowband", narrowBandImage, cLine, nLines, &
      "Output stokes image: ","(a,1l,a)",.false., ok, .false.)


 if (stokesImage) then
    if (.not.narrowBandImage) then
       call getReal("vmin", vmin, cLine, nLines, &
            "Minimum velocity for image (km/s): ", "(a,1pe10.2,1x,a)", -20000., ok, .true.)
       call getReal("vmax", vmax, cLine, nLines, &
            "Maximum velocity for image (km/s): ", "(a,1pe10.2,1x,a)", 20000., ok, .true.)
    else
       call getReal("lmin", vmin, cLine, nLines, &
            "Minimum wavelength for image (angs): ", "(a,1pe10.2,1x,a)", -20000., ok, .true.)
       call getReal("lmax", vmax, cLine, nLines, &
            "Maximum wavelength for image (angs): ", "(a,1pe10.2,1x,a)", 20000., ok, .true.)
    endif
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
    gridDistance = gridDistance * pcTocm / 1.e10
 endif


 call getLogical("plotvel", plotVelocity, cLine, nLines, &
  "Plot velocity vectors: ","(a,1l,a)",.false., ok, .false.)

 call getLogical("sphericity", sphericitytest, cLine, nLines, &
  "Perform sphericity test: ","(a,1l,a)",.false., ok, .false.)



 ! For bipolar jets geometry

 if (geometry(1:4) ==  "jets") then
    call getLogical("lte", lte, cLine, nLines, &
	 "Statistical equ. in LTE: ","(a,1l,1x,a)", .false., ok, .false.)
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
    call getInteger("nlower", nLower, cLine, nLines,"Lower level: ", &
	 & "(a,i2,a)",2,ok,.true.)
    
    call getInteger("nupper", nUpper, cLine, nLines,"Upper level: ", &
	 & "(a,i2,a)",3,ok,.true.)

    
    call getReal("Rmin_bp", Rmin_bp, cLine, nLines, &
	 "Radius of central star  [10^10cm] : ","(a,f5.1,a)", 1.0, ok, .true.)
    call getReal("Rmax_bp", Rmax_bp, cLine, nLines, &
	 "Cutoff radius in [10^10 cm] : ","(a,f5.1,a)", 10.0, ok, .true.)
    call getReal("Vinf_bp", Vinf_bp, cLine, nLines, &
	 "Terminal velocity in  [kms] : ","(a,f5.1,a)", 250.0, ok, .true.)
    call getReal("beta_bp", beta_bp, cLine, nLines, &
	 "Beta in the beta-belocity law [-] : ","(a,f5.1,a)", 0.5, ok, .true.)
    call getReal("Vo_bp", Vo_bp, cLine, nLines, &
	 "A small offset in V in  [km/s] : ","(a,f5.1,a)", 5.0, ok, .true.)
    call getReal("Mdot_bp", Mdot_bp, cLine, nLines, &
	 "Mass loss rate in jets  [M_solar/yr] : ","(a,1PE7.1,a)", 1.0e-9, ok, .true.)
    call getReal("theta_o_bp", theta_o_bp, cLine, nLines, &
	 "Half opening angle [degrees] : ","(a,f5.1,a)", 30.0, ok, .true.)
    call getReal("Tcore_bp", Tcore_bp, cLine, nLines, &
	 "Temperature at the core at Rmin in [10^4 K] : ","(a,f5.1,a)", 5.0, ok, .true.)
    call getReal("e6_bp", e6_bp, cLine, nLines, &
	 "Exponet in tempereture eq. [-] : ","(a,f5.1,a)", 1.0, ok, .true.)

    ! save the variables in an object in jets_mod module.
    call set_jets_parameters(Rmin_bp, Rmax_bp, Vinf_bp, beta_bp,  Vo_bp, &
	 & Mdot_bp, theta_o_bp, Tcore_bp, e6_bp)
    
 end if
 

!   mdot = mdot*msol/(365.25*24.*3600.)
!   vterm = vterm * 1.e5
!   v0 = v0 * 1.e5


 write(*,*) " "

666 continue
end subroutine inputs

subroutine findReal(name, value, cLine, nLines, ok)
 implicit none
 character(len=*) :: name
 real :: value
 character(len=80) :: cLine(*)
 integer :: nLines
 logical :: ok
 integer :: i, j

 ok = .false.
 do i = 1, nLines
  j = len(name)
  if (trim(cLine(i)(1:j)) .eq. name) then
       ok = .true.
       read(cLine(i)(j+1:80),*) value
  endif
 end do
 end subroutine findReal

 
subroutine findInteger(name, value, cLine, nLines, ok)
 implicit none
 character(len=*) :: name
 integer :: value
 character(len=80) :: cLine(*)
 integer :: nLines
 logical :: ok
 integer :: i, j

 ok = .false.
 do i = 1, nLines
  j = len(name)
  if (trim(cLine(i)(1:j)) .eq. name) then
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
 integer :: i, j

 ok = .false.
 do i = 1, nLines
  j = len(name)
  if (trim(cLine(i)(1:j)) .eq. name) then
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
 integer :: i, j

 ok = .false.
 do i = 1, nLines
  j = len(name)
  if (trim(cLine(i)(1:j)) .eq. name) then
       ok = .true.
       read(cLine(i)(j+1:),*) value
  endif
 end do
 value = trim(value)
 end subroutine findString

 subroutine getInteger(name, ival, cLine, nLines, message, format, idef, ok, &
                       musthave)
  character(len=*) :: name
  integer :: ival
  character(len=80) :: cLine(*)
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
       write(*,'(a,a)') name, " must be defined"
       stop
    endif
    ival = idef
    default = " (default)"
  endif
  if (musthave.or.(ival /= idef)) then
     write(*,format) trim(message),ival,default
  endif
 end subroutine getInteger

 subroutine getReal(name, rval, cLine, nLines, message, format, rdef, ok, &
                    musthave)
  character(len=*) :: name
  real :: rval
  logical :: musthave
  character(len=80) :: cLine(*)
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
       write(*,'(a,a)') name, " must be defined"
       stop
    endif
    rval = rdef
    default = " (default)"
 endif
 if (musthave) then
    write(*,format) trim(message),rval,default
 endif
 end subroutine getReal


 subroutine getString(name, rval, cLine, nLines, message, format, rdef, ok, &
                      musthave)
  character(len=*) :: name
  character(len=*) :: rval
  logical :: musthave
  character(len=80) :: cLine(*)
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
       write(*,'(a,a)') name, " must be defined"
       stop
    endif
    rval = rdef
    default = " (default)"
  endif
  write(*,format) trim(message),trim(rval),default
 end subroutine getString


 subroutine getLogical(name, rval, cLine, nLines, message, cformat, rdef, ok, &
                       musthave)
  character(len=*) :: name
  logical :: rval
  logical :: musthave
  character(len=80) :: cLine(*)
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
       write(*,'(a,a)') name, " must be defined"
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
     write(*,'(a,a,a)') trim(message),trueOrFalse,default
  endif
 end subroutine getLogical

 
end module inputs_mod
