module gridio_mod
  use gridtype_mod
  use kind_mod
  use constants_mod                   ! physical constants
  use vector_mod                      ! vector math
  use atom_mod                        ! LTE atomic physics
  use utils_mod
  use octal_mod                       ! octal type for amr
  use amr_mod
  use density_mod                     ! to use generic density function
  use cluster_class
  use cmfgen_class
  use messages_mod
  use ion_mod
  use mpi_global_mod
  use mpi_amr_mod

  implicit none

  public


  interface writeAttributePointerFlexi
     module procedure writeAttributePointerInteger1dFlexi
     module procedure writeAttributePointerReal1dFlexi
     module procedure writeAttributePointerDouble1dFlexi
     module procedure writeAttributePointerDouble2dFlexi
     module procedure writeAttributePointerDouble3dFlexi
     module procedure writeAttributePointerReal2dFlexi
  end interface

  interface writeAttributeStaticFlexi
     module procedure writeAttributeStaticIntegerSingleFlexi
     module procedure writeAttributeStaticRealSingleFlexi
     module procedure writeAttributeStaticCharacterSingleFlexi
     module procedure writeAttributeStaticVectorSingleFlexi
     module procedure writeAttributeStaticInteger1DFlexi
     module procedure writeAttributeStaticLogical1DFlexi
     module procedure writeAttributeStaticVector1DFlexi
     module procedure writeAttributeStaticDouble1dFlexi
     module procedure writeAttributeStaticReal1dFlexi
     module procedure writeAttributeStaticDoubleSingleFlexi
     module procedure writeAttributeStaticLogicalSingleFlexi
  end interface

  interface readSingleFlexi
     module procedure readSingleIntegerFlexi
     module procedure readSingleRealFlexi
     module procedure readSingleDoubleFlexi
     module procedure readSingleCharacterFlexi
     module procedure readSingleLogicalFlexi
     module procedure readSingleVectorFlexi
  end interface

  interface readArrayFlexi
     module procedure readArrayVectorFlexi
     module procedure readArrayRealFlexi
     module procedure readArrayDoubleFlexi
     module procedure readArrayIntegerFlexi
     module procedure readArrayLogicalFlexi
  end interface

  interface readPointerFlexi
     module procedure readDoublePointer1DFlexi
     module procedure readRealPointer1DFlexi
     module procedure readIntegerPointer1DFlexi
     module procedure readVectorPointer1DFlexi
     module procedure readDoublePointer2DFlexi
     module procedure readDoublePointer3DFlexi
     module procedure readRealPointer2DFlexi
  end interface

contains

  subroutine writeAMRgrid(filename,fileFormatted,grid)
    ! writes out the 'grid' for an adaptive mesh geometry  

    implicit none
  
    character(len=*)           :: filename
    logical, intent(in)        :: fileFormatted
    type(GRIDTYPE), intent(in) :: grid
    
    integer, dimension(8) :: timeValues ! system date and time
    integer               :: error      ! error code

    if (fileFormatted) then 
       open(unit=20,iostat=error, file=filename, form="formatted", status="replace")
    else 
       open(unit=20,iostat=error, file=filename, form="unformatted", status="replace")
    end if        
    call writeInfo("Writing AMR grid file to: "//trim(filename),TRIVIAL)
    
    call date_and_time(values=timeValues)
    if (fileFormatted) then
       write(unit=20,fmt=*,iostat=error) timeValues(:)
    else
       write(unit=20,iostat=error) timeValues(:)
    endif

    call writeFileTag(20, "GRIDBEGINS", fileFormatted)
    call writeAttributeStaticFlexi(20,"nLambda", grid%nLambda, fileFormatted)
    call writeAttributeStaticFlexi(20,"flatSpec", grid%flatSpec, fileFormatted)
    call writeAttributeStaticFlexi(20,"adaptive", grid%adaptive, fileFormatted)
    call writeAttributeStaticFlexi(20,"cartesian", grid%cartesian, fileFormatted)
    call writeAttributeStaticFlexi(20,"isotropic", grid%isotropic, fileFormatted)
    call writeAttributeStaticFlexi(20,"hitCore", grid%hitCore, fileFormatted)
    call writeAttributeStaticFlexi(20,"diskRadius", grid%diskRadius, fileFormatted)
    call writeAttributeStaticFlexi(20,"diskNormal", grid%diskNormal, fileFormatted)
    call writeAttributeStaticFlexi(20,"DipoleOffset", grid%DipoleOffset, fileFormatted)
    call writeAttributeStaticFlexi(20,"geometry", grid%geometry, fileFormatted)
    call writeAttributeStaticFlexi(20,"rCore", grid%rCore, fileFormatted)
    call writeAttributeStaticFlexi(20,"lCore", grid%lCore, fileFormatted)
    call writeAttributeStaticFlexi(20,"chanceWind", grid%chanceWindOverTotalContinuum, fileFormatted)
    call writeAttributeStaticFlexi(20,"lineEmission", grid%lineEmission, fileFormatted)
    call writeAttributeStaticFlexi(20,"contEmission", grid%contEmission, fileFormatted)
    call writeAttributeStaticFlexi(20,"doRaman", grid%doRaman, fileFormatted)
    call writeAttributeStaticFlexi(20,"resonanceLine", grid%resonanceLine, fileFormatted)
    call writeAttributeStaticFlexi(20,"rStar1", grid%rStar1, fileFormatted)
    call writeAttributeStaticFlexi(20,"rStar2", grid%rStar2, fileFormatted)
    call writeAttributeStaticFlexi(20,"lumRatio", grid%lumRatio, fileFormatted)
    call writeAttributeStaticFlexi(20,"tempSource", grid%tempSource, fileFormatted)
    call writeAttributeStaticFlexi(20,"starPos1", grid%starPos1, fileFormatted)
    call writeAttributeStaticFlexi(20,"starPos2", grid%starPos2, fileFormatted)
    call writeAttributeStaticFlexi(20,"lambda2", grid%lambda2, fileFormatted)
    call writeAttributeStaticFlexi(20,"maxLevels", grid%maxLevels, fileFormatted)
    call writeAttributeStaticFlexi(20,"maxDepth", grid%maxDepth, fileFormatted)
    call writeAttributeStaticFlexi(20,"halfSmallestSubcell", grid%halfSmallestSubcell, fileFormatted)
    call writeAttributeStaticFlexi(20,"nOctals", grid%nOctals, fileFormatted)
    call writeAttributeStaticFlexi(20,"smoothingFactor", grid%smoothingFactor, fileFormatted)
    call writeAttributeStaticFlexi(20,"oneKappa", grid%oneKappa, fileFormatted)
    call writeAttributeStaticFlexi(20,"rInner", grid%rInner, fileFormatted)
    call writeAttributeStaticFlexi(20,"rOuter", grid%rOuter, fileFormatted)
    call writeAttributeStaticFlexi(20,"amr2dOnly", grid%amr2dOnly, fileFormatted)
    call writeAttributeStaticFlexi(20,"photoionization", grid%photoionization, fileFormatted)
    call writeAttributeStaticFlexi(20,"iDump", grid%iDump, fileFormatted)
    call writeAttributeStaticFlexi(20,"currentTime", grid%currentTime, fileFormatted)
    call writeAttributeStaticFlexi(20,"lamarray", grid%lamarray,fileFormatted)
    call writeAttributePointerFlexi(20,"oneKappaAbs", grid%oneKappaAbs,fileFormatted)
    call writeAttributePointerFlexi(20,"oneKappaSca", grid%oneKappaSca,fileFormatted)
    call writeFileTag(20, "GRIDENDS", fileFormatted)
    call writeOctreePrivateFlexi(grid%octreeRoot,fileFormatted, grid)
    close(unit=20)
    
  contains
  
    recursive subroutine writeOctreePrivateFlexi(thisOctal,fileFormatted, grid)
       ! writes out an octal from the grid octree

       type(octal), intent(in), target :: thisOctal
       logical, intent(in)             :: fileFormatted
       type(gridtype) :: grid
       type(octal), pointer :: thisChild
       integer              :: iChild

       call writeFileTag(20, "OCTALBEGINS", fileFormatted)
       call writeAttributeStaticFlexi(20, "nDepth", thisOctal%nDepth, fileFormatted)
       call writeAttributeStaticFlexi(20, "nChildren", thisOctal%nChildren, fileFormatted)
       call writeAttributeStaticFlexi(20, "indexChild", thisOctal%indexChild, fileFormatted)
       call writeAttributeStaticFlexi(20, "hasChild", thisOctal%hasChild, fileFormatted)
       call writeAttributeStaticFlexi(20, "centre", thisOctal%centre, fileFormatted)
       call writeAttributeStaticFlexi(20, "rho", thisOctal%rho, fileFormatted)
       call writeAttributeStaticFlexi(20, "temperature", thisOctal%temperature, fileFormatted)
       call writeAttributeStaticFlexi(20, "label", thisOctal%label, fileFormatted)
       call writeAttributeStaticFlexi(20, "subcellSize", thisOctal%subcellSize, fileFormatted)
       call writeAttributeStaticFlexi(20, "threeD", thisOctal%threeD, fileFormatted)
       call writeAttributeStaticFlexi(20, "twoD", thisOctal%twoD, fileFormatted)
       call writeAttributeStaticFlexi(20, "maxChildren", thisOctal%maxChildren, fileFormatted)
       call writeAttributeStaticFlexi(20, "cylindrical", thisOctal%cylindrical, fileFormatted)
       call writeAttributeStaticFlexi(20, "splitAzimuthally", thisOctal%splitAzimuthally, fileFormatted)
       call writeAttributeStaticFlexi(20, "phi", thisOctal%phi, fileFormatted)
       call writeAttributeStaticFlexi(20, "dphi", thisOctal%dphi, fileFormatted)
       call writeAttributeStaticFlexi(20, "r", thisOctal%r, fileFormatted)
       call writeAttributeStaticFlexi(20, "parentSubcell", thisOctal%parentSubcell, fileFormatted)
       call writeAttributeStaticFlexi(20, "inStar", thisOctal%inStar, fileFormatted)
       call writeAttributeStaticFlexi(20, "inFlow", thisOctal%inFlow, fileFormatted)
       call writeAttributeStaticFlexi(20, "velocity", thisOctal%velocity, fileFormatted)
       call writeAttributeStaticFlexi(20, "cornervelocity", thisOctal%cornervelocity, fileFormatted)
       call writeAttributePointerFlexi(20, "chiLine", thisOctal%chiLine, fileFormatted)
       call writeAttributePointerFlexi(20, "etaLine", thisOctal%etaLine, fileFormatted)
       call writeAttributePointerFlexi(20, "etaCont", thisOctal%etaCont, fileFormatted)
       call writeAttributePointerFlexi(20, "biasLine3D", thisOctal%biasLine3D, fileFormatted)
       call writeAttributePointerFlexi(20, "biasCont3D", thisOctal%biasCont3D, fileFormatted)
       call writeAttributePointerFlexi(20, "probDistLine", thisOctal%probDistLine, fileFormatted)
       call writeAttributePointerFlexi(20, "probDistCont", thisOctal%probDistCont, fileFormatted)
       call writeAttributePointerFlexi(20, "ne", thisOctal%ne, fileFormatted)
       call writeAttributePointerFlexi(20, "nH", thisOctal%nH, fileFormatted)
       call writeAttributePointerFlexi(20, "nTot", thisOctal%nTot, fileFormatted)
       call writeAttributePointerFlexi(20, "dustType", thisOctal%dustType, fileFormatted)


       call writeAttributePointerFlexi(20, "kappaAbs", thisOctal%kappaAbs, fileFormatted)
       call writeAttributePointerFlexi(20, "kappaSca", thisOctal%kappaSca, fileFormatted)

       call writeAttributePointerFlexi(20, "ionFrac", thisOctal%ionFrac, fileFormatted)
       call writeAttributePointerFlexi(20, "photoIonCoeff", thisOctal%photoIonCoeff, fileFormatted)

       call writeAttributePointerFlexi(20, "molecularLevel", thisOctal%molecularLevel, fileFormatted)
       call writeAttributePointerFlexi(20, "jnu", thisOctal%jnu, fileFormatted)
       call writeAttributePointerFlexi(20, "bnu", thisOctal%bnu, fileFormatted)
       call writeAttributePointerFlexi(20, "molAbundance", thisOctal%molAbundance, fileFormatted)
       call writeAttributePointerFlexi(20, "nh2", thisOctal%nh2, fileFormatted)
       call writeAttributePointerFlexi(20, "microTurb", thisOctal%microTurb, fileFormatted)

       call writeAttributePointerFlexi(20, "N", thisOctal%n, fileFormatted)
       call writeAttributePointerFlexi(20, "departCoeff", thisOctal%departCoeff, fileFormatted)
       call writeAttributePointerFlexi(20, "dustTypeFraction", thisOctal%dustTypeFraction, fileFormatted)


       call writeAttributePointerFlexi(20, "atomAbundance", thisOctal%atomAbundance, fileFormatted)
       call writeAttributePointerFlexi(20, "atomLevel", thisOctal%atomLevel, fileFormatted)
       call writeAttributePointerFlexi(20, "jnuCont", thisOctal%jnuCont, fileFormatted)
       call writeAttributePointerFlexi(20, "jnuLine", thisOctal%jnuLine, fileFormatted)

       call writeAttributeStaticFlexi(20, "mpiThread", thisOctal%mpiThread, fileFormatted)
       call writeFileTag(20, "OCTALENDS", fileFormatted)
       
       if (thisOctal%nChildren > 0) then 
          do iChild = 1, thisOctal%nChildren, 1
             thisChild => thisOctal%child(iChild)
             call writeOctreePrivateFlexi(thisChild,fileFormatted,grid)
          end do
       end if

     end subroutine writeOctreePrivateFlexi
    
    
   end subroutine writeAMRgrid

  subroutine readAMRgrid(filename,fileFormatted,grid)
    ! reads in a previously saved 'grid' for an adaptive mesh geometry  
    implicit none

    character(len=*)            :: filename
    character(len=80)            :: message
    logical, intent(in)         :: fileFormatted
    type(GRIDTYPE), intent(inout) :: grid
    
    integer, dimension(8) :: timeValues    ! system date and time
    integer               :: dummy         
    integer               :: error         ! status code
    integer :: nOctal
    character(len=80) :: absolutePath, inFile
    character(len=20) :: tag

    error = 0

  call unixGetEnv("TORUS_JOB_DIR",absolutePath)
  inFile = trim(absolutePath)//trim(filename)

    if (fileFormatted) then
       open(unit=20, iostat=error, file=inFile, form="formatted", status="old")
       if (error /=0) then
         print *, 'Panic: file open error in readAMRgrid, file:',trim(inFile) ; stop
       end if
       ! read the file's time stamp
       read(unit=20,fmt=*,iostat=error) timeValues 
       if (error /=0) then
         print *, 'Panic: read error in readAMRgrid (formatted timeValues)' ; stop
       end if
    else
       open(unit=20, iostat=error, file=inFile, form="unformatted", status="old")
       if (error /=0) then
         print *, 'Panic: file open error in readAMRgrid, file:',trim(inFile) ; stop
       end if
       ! read the file's time stamp
       read(unit=20,iostat=error) timeValues
       if (error /=0) then
         print *, 'Panic: read error in readAMRgrid (unformatted timeValues)' ; stop
       end if
    end if

    write(message,'(a,a)') "Reading flexible AMR file from: ",trim(filename)
    call writeInfo(message,TRIVIAL)

    write(message,'(a,i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') ' - data file written at: ', &
                          timeValues(1),'/',timeValues(2),'/',&
                          timeValues(3),'  ',timeValues(5),':',timeValues(6)
    call writeInfo(message,TRIVIAL)


    do while(.true.)
       if (fileFormatted) then 
         read(20, '(a)') tag
       else
          read(20) tag
       endif 
       tag = ADJUSTL(tag)
       if (tag == "GRIDBEGINS") cycle
       if (tag == "GRIDENDS") exit

       select case (trim(tag))
          
          case("nLambda")
             call readSingleFlexi(20, grid%nLambda, fileFormatted)
          case("flatSpec")
             call readSingleFlexi(20,grid%flatSpec, fileFormatted)
          case("adaptive")
             call readSingleFlexi(20,grid%adaptive, fileFormatted)
          case("cartesian")
             call readSingleFlexi(20,grid%cartesian, fileFormatted)
          case("isotropic")
             call readSingleFlexi(20,grid%isotropic, fileFormatted)
          case("hitCore")
             call readSingleFlexi(20,grid%hitCore, fileFormatted)
          case("diskRadius")
             call readSingleFlexi(20,grid%diskRadius, fileFormatted)
          case("diskNormal")
             call readSingleFlexi(20,grid%diskNormal, fileFormatted)
          case("DipoleOffset")
             call readSingleFlexi(20,grid%DipoleOffset, fileFormatted)
          case("geometry")
             call readSingleFlexi(20,grid%geometry, fileFormatted)
          case("rCore")
             call readSingleFlexi(20,grid%rCore, fileFormatted)
          case("lCore")
             call readSingleFlexi(20,grid%lCore, fileFormatted)
          case("chanceWind")
             call readSingleFlexi(20,grid%chanceWindOverTotalContinuum, fileFormatted)
          case("lineEmission")
             call readSingleFlexi(20,grid%lineEmission, fileFormatted)
          case("contEmission")
             call readSingleFlexi(20,grid%contEmission, fileFormatted)
          case("doRaman")
             call readSingleFlexi(20,grid%doRaman, fileFormatted)
          case("resonanceLine")
             call readSingleFlexi(20,grid%resonanceLine, fileFormatted)
          case("rStar1")
             call readSingleFlexi(20,grid%rStar1, fileFormatted)
          case("rStar2")
             call readSingleFlexi(20,grid%rStar2, fileFormatted)
          case("lumRatio")
             call readSingleFlexi(20,grid%lumRatio, fileFormatted)
          case("tempSource")
             call readSingleFlexi(20,grid%tempSource, fileFormatted)
          case("starPos1")
             call readSingleFlexi(20,grid%starPos1, fileFormatted)
          case("starPos2")
             call readSingleFlexi(20,grid%starPos2, fileFormatted)
          case("lambda2")
             call readSingleFlexi(20,grid%lambda2, fileFormatted)
          case("maxLevels")
             call readSingleFlexi(20,grid%maxLevels, fileFormatted)
          case("maxDepth")
             call readSingleFlexi(20,grid%maxDepth, fileFormatted)
          case("halfSmallestSubcell")
             call readSingleFlexi(20,grid%halfSmallestSubcell, fileFormatted)
          case("nOctals")
             call readSingleFlexi(20,grid%nOctals, fileFormatted)
          case("smoothingFactor")
             call readSingleFlexi(20,grid%smoothingFactor, fileFormatted)
          case("oneKappa")
             call readSingleFlexi(20,grid%oneKappa, fileFormatted)
          case("rInner")
             call readSingleFlexi(20,grid%rInner, fileFormatted)
          case("rOuter")
             call readSingleFlexi(20,grid%rOuter, fileFormatted)
          case("amr2dOnly")
             call readSingleFlexi(20,grid%amr2dOnly, fileFormatted)
          case("photoionization")
             call readSingleFlexi(20,grid%photoionization, fileFormatted)
          case("iDump")
             call readSingleFlexi(20,grid%iDump, fileFormatted)
          case("currentTime")
             call readSingleFlexi(20,grid%currentTime, fileFormatted)
          case("lamarray")
             call readPointerFlexi(20,grid%lamarray,fileFormatted)
          case("oneKappaAbs")
             call readPointerFlexi(20,grid%oneKappaAbs,fileFormatted)
          case("oneKappaSca")
             call readPointerFlexi(20,grid%oneKappaSca,fileFormatted)
          case DEFAULT
             write(message,'(a,a)') "Unrecognised grid attribute: "//trim(tag)
        end select
     end do




       allocate(grid%octreeRoot)
       nOctal = 0
       call readOctreePrivateFlexi(grid%octreeRoot,null(),fileFormatted, nOctal, grid)

    ! check that we are at the end of the file
    if (fileFormatted) then
       read(unit=20,fmt=*, iostat=error) dummy
    else
       read(unit=20, iostat=error) dummy
    end if
    if (error == 0) then
       print *, 'Panic: read error (expected end of file) in readAMRgridflexi' ; stop
    end if
    
    close(unit=20)

  contains
   
    recursive subroutine readOctreePrivateFlexi(thisOctal,parent,fileFormatted, noctal, grid)
      ! read in an octal to the grid octree

      implicit none
      type(octal), pointer :: thisOctal
      type(octal), pointer :: parent
      type(gridtype) :: grid

      logical, intent(in)  :: fileFormatted
      integer :: nOctal
      character(len=20) :: tag

      type(octal), pointer :: thisChild
      integer              :: iChild

      nOctal = nOctal+1
      thisOctal%parent => parent

      do while (.true.)



         if (fileFormatted) then
            read(20, *) tag
         else
            read(20) tag
         endif
         tag = ADJUSTL(tag)
         if (tag == "OCTALBEGINS") cycle
         if (tag == "OCTALENDS") exit

         select case (tag)
         case("nDepth")
            call readSingleFlexi(20, thisOctal%nDepth, fileFormatted)
         case("nChildren")
            call readSingleFlexi(20, thisOctal%nChildren, fileFormatted)
         case("indexChild")
            call readArrayFlexi(20, thisOctal%indexChild, fileFormatted)
         case("hasChild")
            call readArrayFlexi(20, thisOctal%hasChild, fileFormatted)
         case("centre")
            call readSingleFlexi(20, thisOctal%centre, fileFormatted)
         case("rho")
            call readArrayFlexi(20, thisOctal%rho, fileFormatted)
         case("temperature")
            call readArrayFlexi(20, thisOctal%temperature, fileFormatted)
         case("label")
            call readArrayFlexi(20, thisOctal%label, fileFormatted)
         case("subcellSize")
            call readSingleFlexi(20, thisOctal%subcellSize, fileFormatted)
         case("threeD")
            call readSingleFlexi(20, thisOctal%threeD, fileFormatted)
         case("twoD")
            call readSingleFlexi(20, thisOctal%twoD, fileFormatted)
         case("maxChildren")
            call readSingleFlexi(20, thisOctal%maxChildren, fileFormatted)
         case("cylindrical")
            call readSingleFlexi(20, thisOctal%cylindrical, fileFormatted)
         case("splitAzimuthally")
            call readSingleFlexi(20, thisOctal%splitAzimuthally, fileFormatted)
         case("phi")
            call readSingleFlexi(20, thisOctal%phi, fileFormatted)
         case("dphi")
            call readSingleFlexi(20, thisOctal%dphi, fileFormatted)
         case("r")
            call readSingleFlexi(20, thisOctal%r, fileFormatted)
         case("parentSubcell")
            call readSingleFlexi(20, thisOctal%parentSubcell, fileFormatted)
         case("inStar")
            call readArrayFlexi(20, thisOctal%inStar, fileFormatted)
         case("inFlow")
            call readArrayFlexi(20, thisOctal%inFlow, fileFormatted)
         case("velocity")
            call readArrayFlexi(20, thisOctal%velocity, fileFormatted)
         case("cornervelocity")
            call readArrayFlexi(20, thisOctal%cornervelocity, fileFormatted)
         case("chiLine")
            call readPointerFlexi(20, thisOctal%chiLine, fileFormatted)
         case("etaLine")
            call readPointerFlexi(20, thisOctal%etaLine, fileFormatted)
         case("etaCont")
            call readPointerFlexi(20, thisOctal%etaCont, fileFormatted)
         case("biasLine3D")
            call readPointerFlexi(20, thisOctal%biasLine3D, fileFormatted)
         case("biasCont3D")
            call readPointerFlexi(20, thisOctal%biasCont3D, fileFormatted)
         case("probDistLine")
            call readPointerFlexi(20, thisOctal%probDistLine, fileFormatted)
         case("probDistCont")
            call readPointerFlexi(20, thisOctal%probDistCont, fileFormatted)
         case("ne")
            call readPointerFlexi(20, thisOctal%ne, fileFormatted)
         case("nH")
            call readPointerFlexi(20, thisOctal%nH, fileFormatted)
         case("nTot")
            call readPointerFlexi(20, thisOctal%nTot, fileFormatted)
         case("dustType")
            call readPointerFlexi(20, thisOctal%dustType, fileFormatted)
         case("dustTypeFraction")
            call readPointerFlexi(20, thisOctal%dustTypeFraction, fileFormatted)
          case("mpiThread")
             call readArrayFlexi(20, thisOctal%mpiThread,fileFormatted)
          case("kappaAbs")
             call readPointerFlexi(20, thisOctal%kappaAbs, fileFormatted)
          case("kappaSca")
             call readPointerFlexi(20, thisOctal%kappaSca, fileFormatted)
          case("ionFrac")
             call readPointerFlexi(20, thisOctal%ionFrac, fileFormatted)
          case("photoIonCoeff")
             call readPointerFlexi(20, thisOctal%photoIonCoeff, fileFormatted)
          case("molecularLevel")
             call readPointerFlexi(20, thisOctal%molecularLevel, fileFormatted)
          case("jnu")
             call readPointerFlexi(20, thisOctal%jnu, fileFormatted)
          case("bnu")
             call readPointerFlexi(20, thisOctal%bnu, fileFormatted)
          case("molAbundance")
             call readPointerFlexi(20, thisOctal%molAbundance, fileFormatted)
          case("nh2")
             call readPointerFlexi(20, thisOctal%nh2, fileFormatted)
          case("microTurb")
             call readPointerFlexi(20, thisOctal%microTurb, fileFormatted)
          case("N")
             call readPointerFlexi(20, thisOctal%n, fileFormatted)
          case("departCoeff")
             call readPointerFlexi(20, thisOctal%departCoeff, fileFormatted)
          case("atomAbundance")
             call readPointerFlexi(20, thisOctal%atomAbundance, fileFormatted)
          case("atomLevel")
             call readPointerFlexi(20, thisOctal%atomLevel, fileFormatted)
          case("jnuCont")
             call readPointerFlexi(20, thisOctal%jnuCont, fileFormatted)
          case("jnuLine")
             call readPointerFlexi(20, thisOctal%jnuLine, fileFormatted)

         case DEFAULT
            write(message,*) "Unrecognised tag on read: "//trim(tag)
            call writeWarning(message)
            call readDummyData(20, fileFormatted)
         end select

      end do


      if (thisOctal%nChildren > 0) then 
         allocate(thisOctal%child(1:thisOctal%nChildren)) 
         do iChild = 1, thisOctal%nChildren, 1
            thisChild => thisOctal%child(iChild)
            call readOctreePrivateFlexi(thisChild,thisOctal,fileFormatted, nOctal, grid)
         end do
      end if

    end subroutine readOctreePrivateFlexi
 
  end subroutine readAMRgrid




  subroutine writeAttributeStaticIntegerSingleFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    integer :: value
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "isingle"

    attributeName = name
    if (fileFormatted) then
       write(lUnit,*) attributeName
       write(lUnit,*) dataType
       write(lUnit,*) value
    else
       write(lUnit) attributeName
       write(lUnit) dataType
       write(lUnit) value
    endif
  end subroutine writeAttributeStaticIntegerSingleFlexi

  subroutine writeAttributeStaticRealSingleFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    real :: value
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "rsingle"

    attributeName = name
    if (fileFormatted) then
       write(lUnit,*) attributeName
       write(lUnit,*) dataType
       write(lUnit,*) value
    else
       write(lUnit) attributeName
       write(lUnit) dataType
       write(lUnit) value
    endif
  end subroutine writeAttributeStaticRealSingleFlexi

  subroutine writeAttributeStaticCharacterSingleFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    character(len=*) :: value
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "csingle"

    attributeName = name
    if (fileFormatted) then
       write(lUnit,*) attributeName
       write(lUnit,*)  dataType
       write(lUnit,*) len(trim(value))
       write(lUnit,*) trim(value)
    else
       write(lUnit) attributeName
       write(lUnit)  dataType
       write(lUnit) len(trim(value))
       write(lUnit) trim(value)
    endif
  end subroutine writeAttributeStaticCharacterSingleFlexi

  subroutine writeAttributeStaticInteger1dFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    integer :: value(:)
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "i1darray"

    attributeName = name
    if (fileFormatted) then
       write(lUnit,*) attributeName
       write(lUnit,*) dataType
       write(lUnit,*) size(value)
       write(lUnit,*) value(1:size(value))
    else
       write(lUnit) attributeName
       write(lUnit) dataType
       write(lUnit) size(value)
       write(lUnit) value(1:size(value))
    endif
  end subroutine writeAttributeStaticInteger1dFlexi

  subroutine writeAttributeStaticLogical1dFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    logical :: value(:)
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "l1darray"

    attributeName = name
    if (fileFormatted) then
       write(lUnit,*) attributeName
       write(lUnit,*) dataType
       write(lUnit,*) size(value)
       write(lUnit,*) value(1:size(value))
    else
       write(lUnit) attributeName
       write(lUnit) dataType
       write(lUnit) size(value)
       write(lUnit) value(1:size(value))
    endif
  end subroutine writeAttributeStaticLogical1dFlexi

  subroutine writeAttributeStaticDouble1dFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    real(double) :: value(:)
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "d1darray"

    attributeName = name
    if (fileFormatted) then
       write(lUnit,*) attributeName
       write(lUnit,*) dataType
       write(lUnit,*) size(value)
       write(lUnit,*) value(1:size(value))
    else
       write(lUnit) attributeName
       write(lUnit) dataType
       write(lUnit) size(value)
       write(lUnit) value(1:size(value))
    endif
  end subroutine writeAttributeStaticDouble1dFlexi

  subroutine writeAttributeStaticDoubleSingleFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    real(double) :: value
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "dsingle"

    attributeName = name
    if (fileFormatted) then
       write(lUnit,*) attributeName
       write(lUnit,*) dataType
       write(lUnit,*) value
    else
       write(lUnit) attributeName
       write(lUnit) dataType
       write(lUnit) value
    endif
  end subroutine writeAttributeStaticDoubleSingleFlexi

  subroutine writeAttributeStaticLogicalSingleFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    logical :: value
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "lsingle"

    attributeName = name
    if (fileFormatted) then
       write(lUnit,*) attributeName
       write(lUnit,*) dataType
       write(lUnit,*) value
    else
       write(lUnit) attributeName
       write(lUnit) dataType
       write(lUnit) value
    endif
  end subroutine writeAttributeStaticLogicalSingleFlexi

  subroutine writeAttributeStaticReal1dFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    real :: value(:)
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "r1darray"

    attributeName = name
    if (fileFormatted) then
       write(lUnit,*) attributeName
       write(lUnit,*) dataType
       write(lUnit,*) size(value)
       write(lUnit,*) value(1:size(value))
    else
       write(lUnit) attributeName
       write(lUnit) dataType
       write(lUnit) size(value)
       write(lUnit) value(1:size(value))
    endif
  end subroutine writeAttributeStaticReal1dFlexi

  subroutine writeAttributeStaticVector1dFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    TYPE(vector) :: value(:)
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "v1darray"

    attributeName = name
    if (fileFormatted) then
       write(lUnit,*) attributeName
       write(lUnit,*) dataType
       write(lUnit,*) size(value)
       write(lUnit,*) value(1:size(value))
    else
       write(lUnit) attributeName
       write(lUnit) dataType
       write(lUnit) size(value)
       write(lUnit) value(1:size(value))
    endif
  end subroutine writeAttributeStaticVector1dFlexi

  subroutine writeAttributeStaticVectorSingleFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    TYPE(vector) :: value
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "vsingle"

    attributeName = name
    if (fileFormatted) then
       write(lUnit,*) attributeName
       write(lUnit,*) dataType
       write(lUnit,*) value
    else
       write(lUnit) attributeName
       write(lUnit) dataType
       write(lUnit) value
    endif
  end subroutine writeAttributeStaticVectorSingleFlexi

  subroutine writeAttributePointerDouble1DFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    real(double), pointer :: value(:)
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "d1darray"

    attributeName = name
    if (associated(value)) then
       if (fileFormatted) then
          write(lUnit,*) attributeName
          write(lUnit,*) dataType
          write(lUnit,*) SIZE(value)
          write(lUnit,*) value(1:SIZE(value))
       else
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) SIZE(value)
          write(lUnit) value(1:SIZE(value))
       endif
    endif
  end subroutine writeAttributePointerDouble1DFlexi

  subroutine writeAttributePointerInteger1DFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    integer, pointer :: value(:)
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "i1darray"

    attributeName = name
    if (associated(value)) then
       if (fileFormatted) then
          write(lUnit,*) attributeName
          write(lUnit,*) dataType
          write(lUnit,*) SIZE(value)
          write(lUnit,*) value(1:SIZE(value))
       else
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) SIZE(value)
          write(lUnit) value(1:SIZE(value))
       endif
    endif
  end subroutine writeAttributePointerInteger1DFlexi

  subroutine writeAttributePointerReal1DFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    real, pointer :: value(:)
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "r1darray"

    attributeName = name
    if (associated(value)) then
       if (fileFormatted) then
          write(lUnit,*) attributeName
          write(lUnit,*) dataType
          write(lUnit,*) SIZE(value)
          write(lUnit,*) value(1:SIZE(value))
       else
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) SIZE(value)
          write(lUnit) value(1:SIZE(value))
       endif
    endif
  end subroutine writeAttributePointerReal1DFlexi

  subroutine writeAttributePointerDouble2DFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    real(double), pointer :: value(:,:)
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "d2darray"

    attributeName = name
    if (associated(value)) then
       if (fileFormatted) then
          write(lUnit,*) attributeName
          write(lUnit,*) dataType
          write(lUnit,*) SIZE(value,1),SIZE(value,2)
          write(lUnit,*) value(1:SIZE(value,1),1:SIZE(value,2))
       else
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) SIZE(value,1),SIZE(value,2)
          write(lUnit) value(1:SIZE(value,1),1:SIZE(value,2))
       endif
    endif
  end subroutine writeAttributePointerDouble2DFlexi

  subroutine writeAttributePointerDouble3DFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    real(double), pointer :: value(:,:,:)
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "d3darray"

    attributeName = name
    if (associated(value)) then
       if (fileFormatted) then
          write(lUnit,*) attributeName
          write(lUnit,*) dataType
          write(lUnit,*) SIZE(value,1), SIZE(value,2), SIZE(value,3)
          write(lUnit,*) value(1:SIZE(value,1),1:SIZE(value,2),1:SIZE(value,3))
       else
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) SIZE(value,1),SIZE(value,2), SIZE(value,3)
          write(lUnit) value(1:SIZE(value,1),1:SIZE(value,2),1:SIZE(value,3))
       endif
    endif
  end subroutine writeAttributePointerDouble3DFlexi

  subroutine writeAttributePointerReal2DFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    real, pointer :: value(:,:)
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "r2darray"

    attributeName = name
    if (associated(value)) then
       if (fileFormatted) then
          write(lUnit,*) attributeName
          write(lUnit,*) dataType
          write(lUnit,*) SIZE(value,1),SIZE(value,2)
          write(lUnit,*) value(1:SIZE(value,1),1:SIZE(value,2))
       else
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) SIZE(value,1),SIZE(value,2)
          write(lUnit) value(1:SIZE(value,1),1:SIZE(value,2))
       endif
    endif
  end subroutine writeAttributePointerReal2DFlexi



    subroutine writeFileTag(lUnit, tag, fileFormatted)
      integer :: lUnit
      character(len=*) tag
      logical :: fileFormatted
      character(len=20) :: wtag
      wtag = tag

      if (fileFormatted) then
         write(lUnit,*) wtag
      else
         write(lUnit) wtag
      endif
    end subroutine writeFileTag

    subroutine readSingleIntegerFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      integer :: value
      logical :: fileFormatted

      call testDataType("isingle", fileFormatted)
      if (fileFormatted) then
         read(lUnit,*) value
      else
         read(lUnit) value
      endif
    end subroutine readSingleIntegerFlexi

    subroutine readSingleLogicalFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      logical :: value
      logical :: fileFormatted

      call testDataType("lsingle", fileFormatted)
      if (fileFormatted) then
         read(lUnit,*) value
      else
         read(lUnit) value
      endif
    end subroutine readSingleLogicalFlexi

    subroutine readSingleRealFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      real :: value
      logical :: fileFormatted

      call testDataType("rsingle", fileFormatted)
      if (fileFormatted) then
         read(lUnit,*) value
      else
         read(lUnit) value
      endif
    end subroutine readSingleRealFlexi

    subroutine readSingleVectorFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      type(VECTOR) :: value
      logical :: fileFormatted

      call testDataType("vsingle", fileFormatted)
      if (fileFormatted) then
         read(lUnit,*) value
      else
         read(lUnit) value
      endif
    end subroutine readSingleVectorFlexi

    subroutine readSingleDoubleFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      real(double) :: value
      logical :: fileFormatted

      call testDataType("dsingle", fileFormatted)
      if (fileFormatted) then
         read(lUnit,*) value
      else
         read(lUnit) value
      endif
    end subroutine readSingleDoubleFlexi

    subroutine readSingleCharacterFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      integer :: i
      character(len=*) :: value
      logical :: fileFormatted

      call testDataType("csingle", fileFormatted)
      if (fileFormatted) then
         read(lUnit,*) i
         read(lUnit,*) value(1:i)
      else
         read(lUnit) i
         read(lUnit) value(1:i)
      endif
    end subroutine readSingleCharacterFlexi
    
    subroutine readDoublePointer1DFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      real(double), pointer :: value(:)
      logical :: fileFormatted
      integer :: n

      call testDataType("d1darray", fileFormatted)
      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      if (fileFormatted) then
         read(lUnit,*) n
         allocate(value(1:n))
         read(lUnit,*) value(1:n)
      else
         read(lUnit) n
         allocate(value(1:n))
         read(lUnit) value(1:n)
      endif
    end subroutine readDoublePointer1DFlexi

    subroutine readDoublePointer2DFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      real(double), pointer :: value(:,:)
      logical :: fileFormatted
      integer :: n, m

      call testDataType("d2darray", fileFormatted)
      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      if (fileFormatted) then
         read(lUnit,*) n,m
         allocate(value(1:n,1:m))
         read(lUnit,*) value(1:n,1:m)
      else
         read(lUnit) n, m
         allocate(value(1:n,1:m))
         read(lUnit) value(1:n,1:m)
      endif
    end subroutine readDoublePointer2DFlexi

    subroutine readDoublePointer3DFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      real(double), pointer :: value(:,:,:)
      logical :: fileFormatted
      integer :: n, m, j

      call testDataType("d3darray", fileFormatted)
      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      if (fileFormatted) then
         read(lUnit,*) n,m,j
         allocate(value(1:n,1:m,1:j))
         read(lUnit,*) value(1:n,1:m,1:j)
      else
         read(lUnit) n, m, j
         allocate(value(1:n,1:m,1:j))
         read(lUnit) value(1:n,1:m,1:j)
      endif
    end subroutine readDoublePointer3DFlexi

    subroutine readRealPointer2DFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      real, pointer :: value(:,:)
      logical :: fileFormatted
      integer :: n, m

      call testDataType("r2darray", fileFormatted)
      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      if (fileFormatted) then
         read(lUnit,*) n,m
         allocate(value(1:n,1:m))
         read(lUnit,*) value(1:n,1:m)
      else
         read(lUnit) n, m
         allocate(value(1:n,1:m))
         read(lUnit) value(1:n,1:m)
      endif
    end subroutine readRealPointer2DFlexi

    subroutine readRealPointer1DFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      real, pointer :: value(:)
      logical :: fileFormatted
      integer :: n

      call testDataType("r1darray", fileFormatted)
      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      if (fileFormatted) then
         read(lUnit,*) n
         allocate(value(1:n))
         read(lUnit,*) value(1:n)
      else
         read(lUnit) n
         allocate(value(1:n))
         read(lUnit) value(1:n)
      endif
    end subroutine readRealPointer1DFlexi

    subroutine readIntegerPointer1DFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      integer, pointer :: value(:)
      logical :: fileFormatted
      integer :: n

      call testDataType("i1darray", fileFormatted)
      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      if (fileFormatted) then
         read(lUnit,*) n
         allocate(value(1:n))
         read(lUnit,*) value(1:n)
      else
         read(lUnit) n
         allocate(value(1:n))
         read(lUnit) value(1:n)
      endif
    end subroutine readIntegerPointer1DFlexi

    subroutine readVectorPointer1DFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      type(VECTOR), pointer :: value(:)
      logical :: fileFormatted
      integer :: n

      call testDataType("v1darray", fileFormatted)
      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      if (fileFormatted) then
         read(lUnit,*) n
         allocate(value(1:n))
         read(lUnit,*) value(1:n)
      else
         read(lUnit) n
         allocate(value(1:n))
         read(lUnit) value(1:n)
      endif
    end subroutine readVectorPointer1DFlexi

    subroutine readArrayVectorFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      type(VECTOR) :: value(:)
      logical :: fileFormatted
      integer :: n

      call testDataType("v1darray", fileFormatted)
      if (fileFormatted) then
         read(lUnit,*) n
         if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
         read(lUnit,*) value(1:n)
      else
         read(lUnit) n
         if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
         read(lUnit) value(1:n)
      endif
    end subroutine readArrayVectorFlexi

    subroutine readArrayLogicalFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      logical :: value(:)
      logical :: fileFormatted
      integer :: n

      call testDataType("l1darray", fileFormatted)
      if (fileFormatted) then
         read(lUnit,*) n
         if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
         read(lUnit,*) value(1:n)
      else
         read(lUnit) n
         if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
         read(lUnit) value(1:n)
      endif
    end subroutine readArrayLogicalFlexi

    subroutine readArrayIntegerFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      integer :: value(:)
      logical :: fileFormatted
      integer :: n

      call testDataType("i1darray", fileFormatted)
      if (fileFormatted) then
         read(lUnit,*) n
         if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
         read(lUnit,*) value(1:n)
      else
         read(lUnit) n
         if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
         read(lUnit) value(1:n)
      endif
    end subroutine readArrayIntegerFlexi

    subroutine readArrayRealFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      real :: value(:)
      logical :: fileFormatted
      integer :: n

      call testDataType("r1darray", fileFormatted)
      if (fileFormatted) then
         read(lUnit,*) n
         if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
         read(lUnit,*) value(1:n)
      else
         read(lUnit) n
         if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
         read(lUnit) value(1:n)
      endif
    end subroutine readArrayRealFlexi

    subroutine readArrayDoubleFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      real(double) :: value(:)
      logical :: fileFormatted
      integer :: n

      call testDataType("d1darray", fileFormatted)
      if (fileFormatted) then
         read(lUnit,*) n
         if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
         read(lUnit,*) value(1:n)
      else
         read(lUnit) n
         if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
         read(lUnit) value(1:n)
      endif
    end subroutine readArrayDoubleFlexi

    subroutine testDataType(thisType, fileFormatted)
      character(len=*) :: thisType
      logical :: fileFormatted
      character(len=10) :: dataType
      character(len=80) :: message
      if (fileFormatted) then
         read(20, *) dataType
      else
         read(20) dataType
      endif
      dataType = ADJUSTL(dataType)
      if (dataType.ne.thisType) then
         write(message,*) &
              "Data type for read statement "//trim(thisType)//" does not agree with file type of "//trim(dataType)
         call writeFatal(message)
      endif
    end subroutine testDataType

    subroutine readDummyData(lUnit, fileFormatted)
      integer :: lUnit
      logical :: fileFormatted
      character(len=10) :: dataType
      type(VECTOR) :: vDummy
      real :: rDummy
      real, pointer :: rArray(:)
      integer, pointer :: iArray(:)
      logical, pointer :: lArray(:)
      type(VECTOR), pointer :: vArray(:)
      real(double), pointer :: dArray(:)
      real(double), pointer :: d2Array(:,:)
      real(double), pointer :: d3Array(:,:,:)
      real(double) :: dDummy
      integer :: iDummy
      logical :: lDummy
      character(len=20) cDummy
      integer :: n, m, j

      if (fileFormatted) then
         read(lUnit,*) dataType
      else
         read(lUnit) dataType
      endif
      dataType = ADJUSTL(dataType)
      select case (dataType)
         case("rsingle")
            if (fileFormatted) then
               read(lUnit,*) rDummy
            else
               read(lUnit) rDummy
            endif
         case("dsingle")
            if (fileFormatted) then
               read(lUnit,*) dDummy
            else
               read(lUnit) dDummy
            endif
         case("isingle")
            if (fileFormatted) then
               read(lUnit,*) iDummy
            else
               read(lUnit) iDummy
            endif
         case("vsingle")
            if (fileFormatted) then
               read(lUnit,*) vDummy
            else
               read(lUnit) vDummy
            endif
         case("lsingle")
            if (fileFormatted) then
               read(lUnit,*) lDummy
            else
               read(lUnit) lDummy
            endif
         case("csingle")
            if (fileFormatted) then
               read(lUnit,*) n
               read(lUnit,*) cDummy(1:n)
            else
               read(lUnit) n
               read(lUnit) cDummy(1:n)
            endif
         case("r1darray")
            if (fileFormatted) then
               read(lUnit,*) n
               allocate(rArray(1:n))
               read(lUnit,*) rArray(1:n)
            else
               read(lUnit) n
               allocate(rArray(1:n))
               read(lUnit) rArray(1:n)
            endif
            deallocate(rArray)
         case("v1darray")
            if (fileFormatted) then
               read(lUnit,*) n
               allocate(vArray(1:n))
               read(lUnit,*) vArray(1:n)
            else
               read(lUnit) n
               allocate(vArray(1:n))
               read(lUnit) vArray(1:n)
            endif
            deallocate(vArray)
         case("d1darray")
            if (fileFormatted) then
               read(lUnit,*) n
               allocate(dArray(1:n))
               read(lUnit,*) dArray(1:n)
            else
               read(lUnit) n
               allocate(dArray(1:n))
               read(lUnit) dArray(1:n)
            endif
            deallocate(dArray)
         case("i1darray")
            if (fileFormatted) then
               read(lUnit,*) n
               allocate(iArray(1:n))
               read(lUnit,*) iArray(1:n)
            else
               read(lUnit) n
               allocate(iArray(1:n))
               read(lUnit) iArray(1:n)
            endif
            deallocate(iArray)
         case("l1darray")
            if (fileFormatted) then
               read(lUnit,*) n
               allocate(lArray(1:n))
               read(lUnit,*) lArray(1:n)
            else
               read(lUnit) n
               allocate(lArray(1:n))
               read(lUnit) lArray(1:n)
            endif
            deallocate(iArray)
         case("d2darray")
            if (fileFormatted) then
               read(lUnit,*) n,m
               allocate(d2Array(1:n,1:m))
               read(lUnit,*) d2Array(1:n,1:m)
            else
               read(lUnit) n,m
               allocate(d2Array(1:n,1:m))
               read(lUnit) d2Array(1:n,1:m)
            endif
            deallocate(dArray)
         case("d3darray")
            if (fileFormatted) then
               read(lUnit,*) n,m,j
               allocate(d3Array(1:n,1:m,1:j))
               read(lUnit,*) d3Array(1:n,1:m,1:j)
            else
               read(lUnit) n,m,j
               allocate(d3Array(1:n,1:m,1:j))
               read(lUnit) d3Array(1:n,1:m,1:j)
            endif
            deallocate(d3Array)
         case DEFAULT
            write(*,*) "Data type not recognised: ",trim(dataType)
            stop
         end select
       end subroutine readDummyData
            


end module gridio_mod
