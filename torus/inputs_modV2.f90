module inputs_mod

  use vector_mod
  use unix_mod
  use messages_mod
  use kind_mod
  use input_variables
  use constants_mod

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

    nBlobs = 0
    nLines = 0

    call unixGetEnv("TORUS_JOB_DIR",absolutePath)
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

    call getLogical("readgrid", readGrid, cLine, nLines, &
         "Read grid file: ","(a,1l,1x,a)", .false., ok, .true.)

    if (.not.readgrid) then
       call readGridInitParameters(cLine, nLines)
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
         "Include atomic physics in calculation: ","(a,1l,1x,a)", .false., ok, .false.)

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
         & "(a,i3,a)",5,ok,.false.)

    call getInteger("maxdepthamr", maxDepthAMR, cLine, nLines, "Maximum cell depth of AMR grid: ", &
         & "(a,i3,a)",31,ok,.false.)

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
       call getDouble(keyword, sourceRadius(i), rsol, cLine, nLines, &
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
       call getVector(keyword, sourcePos(i), 1.d10, cLine, nLines, &
            "Source position (10^10 cm): ","(a,3(1pe12.3),a)",VECTOR(0.d0, 0.d0, 0.d0), ok, .true.)
    enddo

  end subroutine readSourceParameters

  subroutine readMolecularPhysicsParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
    logical :: ok

       call getLogical("lte", lte, cLine, nLines, &
            "Read in LTE grid: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("suppresswarnings", suppressWarnings, cLine, nLines, &
            "Suppress Warnings: ","(a,l,1x,a)",.false., ok, .false.) 
       call getLogical("addnewmoldata", addnewmoldata, cLine, nLines, &
            "Add new molecular data to non-molecular grid: ","(a,l,1x,a)",.false., ok, .false.)
       call getString("moleculefile", moleculefile, cLine, nLines, &
            "Input molecule filename: ","(a,a,1x,a)","none", ok, .false.)
       call getString("molfilein", molFilenameIn, cLine, nLines, &
            "Input Lucy grid filename: ","(a,a,1x,a)","none", ok, .false.)
       call getString("molfileout", molFilenameOut, cLine, nLines, &
            "Output Lucy grid filename: ","(a,a,1x,a)","none", ok, .false.)
       call getLogical("writemol", writeMol, cLine, nLines, &
            "Write lucy grid file: ","(a,1l,1x,a)", .false., ok, .true.)
       call getLogical("readmol", readMol, cLine, nLines, &
            "Read molecular grid file: ","(a,1l,1x,a)", .false., ok, .true.)
       call getLogical("writeLucy", writeLucy, cLine, nLines, &
            "Write lucy grid file: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("readLucy", readLucy, cLine, nLines, &
            "Read molecular grid file: ","(a,1l,1x,a)", .false., ok, .false.)
       call getLogical("openlucy", openLucy, cLine, nLines, &
            "Open Existing lucy file: ","(a,1l,1x,a)", .false., ok, .false.)
       call getReal("distance", gridDistance, real(pctocm), cLine, nLines, &
            "Grid distance (pc): ","(a,f6.1,1x,a)", 1., ok, .true.)
       call getLogical("dongstep", dongstep, cLine, nLines, &
               "Use Ng Acceleration: ","(a,1l,a)", .false., ok, .false.)
       call getLogical("quasi", quasi, cLine, nLines, &
               "Use Quasirandom numbers: ","(a,1l,a)", .false., ok, .false.)
       call getReal("tolerance", tolerance, 1., cLine, nLines, &
            "Maximum Fractional Change in level populations:","(a,f4.1,1x,a)", 0.01, ok, .true.)
       call getReal("vturb", vturb, 1., cLine, nLines, &
            "Subsonic turbulent velocity (km/s):","(a,f4.1,1x,a)", 0.3, ok, .true.)
       call getLogical("noturb", noturb, cLine, nLines, &
            "No microturbulence","(a,1l,a)",.false., ok, .false.)
       call getInteger("setmaxlevel", setmaxlevel, cLine, nLines, &
            "Maximum molecular level to be considered:","(a,i2,1x,a)", 0, ok, .false.)
       call getReal("molAbundance", molAbundance, 1., cLine, nLines, &
            "Molecular Abundance:","(a,e12.5,1x,a)", 1e-9, ok, .true.)
       call getLogical("useDust", useDust, cLine, nLines, &
            "Calculate continuum emission from dust:", "(a,1l,1x,a)", .false., ok, .true.)
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
  end subroutine readPhotoionEquilibriumParameters

  subroutine readHydrodynamicsParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
  end subroutine readHydrodynamicsParameters

  subroutine readDataCubeParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
  end subroutine readDataCubeParameters

  subroutine readImageParameters(cLine, nLines)
    character(len=80) :: cLine(:)
    integer :: nLines
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

