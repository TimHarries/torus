module inputs_mod

use utils_mod

implicit none

public

contains

subroutine inputs(nPhotons, nx, ny, nz, nr, nmu, nphi, &
     inclination, geometry, gridtype, mie, &
     rTorus, rOuter, scale, rho, lamStart, lamEnd, nLambda, &
     outFile, grainSize, radius, kfac, rMin, rMaj, &
     aMin, aMax, qDist, height, kappaSca, kappaAbs, &
     screened, useBias, maxScat, Teff, mdot, vterm, plezModelOn, &
     opaqueCore, fillTio, lineEmission, lamLine, &
     fillThomson, fillRayleigh, device, movie, shellfrac, &
     inputOK, probContPhoton, &
     nPhase, nBlobs, contrast, beta, opacityDataFile, &
     intProFilename, intProFilename2, distortionType, vRot, plotVelocity, &
     nStartPhase, nEndphase, scaleDensity, grainType,nSpiral, thinLine, &
     freshBlobs, stokesImage, vmin, vmax, gridDistance, pencilBeam, &
     dustBlobDistance, phiDustBlob, xDustBlobSize, yDustBlobSize, &
     zDustBlobSize, secondSource, secondSourcePosition, lumRatio, &
     binarySep, momRatio, vo6, rCore, rInner, coreEmissionLine, &
     velWidthCoreEmission, relIntCoreEmission, contFluxFile, lineOff, &
     logMassLossRate, ramVel, ramanDist, phaseOffset, &
     contFluxFile2, mass1, mass2, period, radius1, radius2, mdot1, mdot2, &
     popFilename, readPops, writePops, vNought1, vNought2, beta1, beta2, &
     vterm1, vterm2, temp1, temp2, shockWidth, shockFac, doRaman, deflectionAngle, lte, blobContrast, &
     pvimage, slitPosition1,slitPosition2, slitPA, slitWidth, slitLength, &
     nSlit, np, nv, vfwhm, pfwhm, vSys, useNdf, mCore, diskTemp,curtains, &
     dipoleOffset, enhance, v0, o6width, misc, nLower, nUpper, &
     nSpot, fSpot, tSpot, thetaSpot, phiSpot, photLine)

  use constants_mod
  use vector_mod
  use unix_mod
  implicit none
  logical :: curtains, enhance, done
  real :: dipoleOffset
  logical :: lte, useNdf
  character(len=*) :: misc
  logical pvimage
  real :: mCore, diskTemp, v0, o6width
  real :: slitPA, slitWidth, slitLength
  integer :: nSlit, np, nv
  real :: vfwhm, pfwhm, vSys
  type(VECTOR) :: slitPosition1, slitPosition2
  integer :: nphase, nBlobs, nStartPhase, nEndPhase
  logical :: doRaman
  integer :: nLower, nUpper
  real :: contrast, beta
  real :: radius1, radius2, mdot1, mdot2
  real :: vNought1, vNought2, beta1, beta2
  real :: temp1, temp2
  real :: shockWidth, shockFac
  real :: vterm1, vterm2
  real :: phaseOffset
  logical :: readPops, writePops
  character(len=*) :: popFilename
  real :: logMassLossRate
  real :: period
  real :: mass1, mass2
  logical :: lineOff
  integer :: errNo, i
  integer :: nSpiral
  logical :: freshBlobs
  real :: blobContrast
  logical :: thinLine
  logical :: stokesImage
  logical :: pencilBeam
  real :: vo6
  real :: vMin, vMax
  real :: binarysep, momRatio
  character(len=*) :: grainType
  character(len=*) :: device, opacityDataFile, intProFilename,  intProFilename2
  character(len=*) :: distortionType
  character(len=*) :: contFluxFile, contFluxFile2
  character(len=*) :: ramanDist
  logical :: plotVelocity
  real :: rTorus, rOuter, shellfrac
  real :: scaleDensity
  logical :: fillThomson
  logical :: movie
  real :: probContPhoton
  real :: inclination
  real :: vRot
  logical :: opaqueCore, lineEmission
  real :: lamLine
  real :: grainSize
  real :: Teff
  logical :: fillRayleigh
  logical :: fillTio
  integer :: nx,ny,nz
  integer :: nr, nmu, nphi
  logical :: useBias
  logical :: plezModelOn
  real :: kappaSca, kappaAbs
  real :: scale
  real :: height
  real :: rho
  logical :: screened
  real :: aMin, aMax, qDist
  character(len=10) geometry, gridtype
  integer :: nPhotons
  logical :: mie
  logical :: ok
  real :: mdot, vterm
  real :: lamStart, lamEnd
  real :: radius, kfac
  real :: rMin, rMaj
  integer :: nLines, nLambda
  character(len=80) :: cLine(100) 
  character(len=9) :: default
  character(len=*) :: outFile
  real :: gridDistance, ramVel
  integer :: maxScat
  logical :: inputOK
  real :: dustBlobDistance, phiDustBlob
  real :: xDustBlobSize, yDustBlobSize, zDustBlobSize
  real :: rCore, rInner

  logical :: coreEmissionLine
  real :: velWidthCoreEmission
  real :: relIntCoreEmission
  real :: deflectionAngle

  ! variables for a second source of radiation

  logical :: secondSource                    ! second photon source?
  type(VECTOR) :: secondSourcePosition       ! the position of it
  real :: lumRatio                           ! lumonsity ratio

  ! character vars for unix environ

  character(len=80) :: dataDirectory

  ! Spot stuff
  
  integer :: nSpot                       ! number of spots
  real :: fSpot                          ! factional area coverage of spots
  real :: tSpot                          ! spot temperatures
  real :: thetaSpot, phiSpot             ! spot coords
  logical :: photLine                    ! photospheric line production
  


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

 

 call getInteger("nlambda", nlambda, cLine, nLines, &
     "Number of wavelength/velocity bins: ", "(a,i4,1x,a)", 20, ok, .true.)




  call getString("distortion", distortionType, cLine, nLines, &
       "Distortion type: ","(a,a,1x,a)","none", ok, .false.)

  call getReal("inclination", inclination, cLine, nLines, &
       "Inclination: ","(a,f4.1,1x,a)", 90., ok, .true.)


  call getReal("scale", scale, cLine, nLines, &
       "Scale (rsolar): ","(a,f6.1,1x,a)", 1000., ok, .false.) 


  call getString("gridtype", gridtype, cLine, nLines, &
       "Grid type: ","(a,a,a)","cartesian",ok, .true.)

  if (gridtype.eq."cartesian") then
     call getInteger("nx", nx, cLine, nLines, &
          "Grid size in x-direction: ", "(a,i3,1x,a)", 100, ok, .true.)
     call getInteger("ny", ny, cLine, nLines, &
          "Grid size in y-direction: ", "(a,i3,1x,a)", 100, ok, .true.)
     call getInteger("nz", nz, cLine, nLines, &
          "Grid size in z-direction: ", "(a,i3,1x,a)", 100, ok, .true.)
  endif

  if (gridtype.eq."polar") then
     call getInteger("nr", nr, cLine, nLines, &
          "Grid size in r-direction: ", "(a,i3,1x,a)", 100, ok, .true.)
     call getInteger("nmu", nmu, cLine, nLines, &
          "Grid size in mu-direction: ", "(a,i3,1x,a)", 100, ok, .true.)
     call getInteger("nphi", nphi, cLine, nLines, &
          "Grid size in phi-direction: ", "(a,i3,1x,a)", 100, ok, .true.)
  endif


  call getString("geometry", geometry, cLine, nLines, &
       "Geometry: ","(a,a,a)","sphere",ok, .true.)

!  call getLogical("raman", doRaman, cLine, nLines, &
!            "Raman-scattering model: ","(a,1l,1x,a)", .false., ok, .false.)


  if (geometry(1:5) .eq. "torus") then
     call getReal("rtorus", rTorus, cLine, nLines, &
          "Radius of torus (AU): ","(a,f5.1,a)", 0.7, ok, .true.)
     call getReal("router", rOuter, cLine, nLines, &
          "Radius of torus x-section (AU): ","(a,f5.1,a)", 0.1, ok, .true.)
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
     "rMajor: ", "(a,f5.1,a)", 0.5, ok, .true.)
   call getReal("rmin", rmin, cLine, nLines, &
     "rMinor: ", "(a,f5.1,a)", 0.5, ok, .true.)
    call getReal("rho", rho, cLine, nLines, &
    "Density (xxx): ","(a,1p,e10.2,1p,1x,a)", 1.e-6, ok, .true.)
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
     call getString("graintype", grainType, cLine, nLines, &
  "Grain type: ","(a,a,1x,a)","silicate", ok, .true.)
 endif

 if (.not.mie) then
 default = " "
 call findReal("kappasca", kappaSca, cLine, nLines, ok)
 if (ok) then
   write(*,'(a,1pe12.4)') "Scattering cross-section: ",kappaSca
 endif
 default = " "
 call findReal("kappaabs", kappaAbs, cLine, nLines, ok)
 if (ok) then
   write(*,'(a,1pe12.4)') "Absorption cross-section: ",kappaAbs
 endif
endif


 call getLogical("secondsource", secondSource, cLine, nLines, &
   "Second source: ","(a,1l,a)", .false., ok, .false.)

 if (secondSource) then
    secondSourcePosition = VECTOR(0.,0.,0.)
    call getReal("secondpos", secondSourcePosition%x, cLine, nLines, &
         "Second source distance: ","(a,e7.1,1x,a)", 1.e14, ok, .true.)
    call getReal("lumratio", lumratio, cLine, nLines, &
         "Luminosity Ratio (p/s): ","(a,f7.2,1x,a)", 1., ok, .true.)
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
    
    call getReal("velfwhm", velWidthCoreEmission, cLine, nLines, &
   "FWHM Velocity width of core emission line (km/s): ","(a,f7.0,1x,a)", 1000., ok, .true.)
    call getReal("relint", relintCoreEmission, cLine, nLines, &
   "Relative intensity of core emission line: ","(a,f7.0,1x,a)", 10., ok, .true.)
    velWidthCoreEmission = velWidthCoreEmission * 1.e5
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


 call getLogical("rayleigh", fillRayleigh, cLine, nLines, &
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

 if (stokesImage) then
    call getReal("vmin", vmin, cLine, nLines, &
     "Minimum velocity for image (km/s): ", "(a,f7.1,1x,a)", -20000., ok, .true.)
    call getReal("vmax", vmax, cLine, nLines, &
     "Maximum velocity for image (km/s): ", "(a,f7.1,1x,a)", 20000., ok, .true.)
!    vmax = (vmax*1.e5) / cSpeed
!    vmin = (vmin*1.e5) / cSpeed
    vmin = vmin * 1.e5
    vmax = vmax * 1.e5
    
 endif

 call getLogical("pvimage", pvimage, cLine, nLines, &
      "Output pv image: ","(a,1l,a)",.false., ok, .false.)
 if (pvimage) then
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
  character(len=9) :: default
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
    default = "(default)"
  endif
  if (musthave) then
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
  character(len=9) :: default
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
    default = "(default)"
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
  character(len=9) :: default
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
  character(len=9) :: default
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
