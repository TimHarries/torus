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
  use cmfgen_class
  use messages_mod
  use mpi_global_mod
#ifdef MPI
  use mpi_amr_mod
#endif
  use parallel_mod, only : torus_mpi_barrier
  use zlib_mod
 implicit none

  public


  interface writeAttributePointerFlexi
     module procedure writeAttributePointerInteger1dFlexi
     module procedure writeAttributePointerLogical1dFlexi
     module procedure writeAttributePointerReal1dFlexi
     module procedure writeAttributePointerVector1dFlexi
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
     module procedure readLogicalPointer1DFlexi
     module procedure readIntegerPointer1DFlexi
     module procedure readVectorPointer1DFlexi
     module procedure readDoublePointer2DFlexi
     module procedure readDoublePointer3DFlexi
     module procedure readRealPointer2DFlexi
  end interface

#ifdef MPI

  interface receiveSingleFlexi
     module procedure receiveSingleIntegerFlexi
     module procedure receiveSingleRealFlexi
     module procedure receiveSingleDoubleFlexi
     module procedure receiveSingleCharacterFlexi
     module procedure receiveSingleLogicalFlexi
     module procedure receiveSingleVectorFlexi
  end interface

  interface receivePointerFlexi
     module procedure receiveDoublePointer1DFlexi
     module procedure receiveRealPointer1DFlexi
     module procedure receiveLogicalPointer1DFlexi
     module procedure receiveIntegerPointer1DFlexi
     module procedure receiveVectorPointer1DFlexi
     module procedure receiveDoublePointer2DFlexi
     module procedure receiveDoublePointer3DFlexi
     module procedure receiveRealPointer2DFlexi
  end interface


  interface receiveArrayFlexi
     module procedure receiveArrayVectorFlexi
     module procedure receiveArrayRealFlexi
     module procedure receiveArrayDoubleFlexi
     module procedure receiveArrayIntegerFlexi
     module procedure receiveArrayLogicalFlexi
  end interface

  interface sendAttributePointerFlexi
     module procedure sendAttributePointerInteger1dFlexi
     module procedure sendAttributePointerLogical1dFlexi
     module procedure sendAttributePointerReal1dFlexi
     module procedure sendAttributePointerVector1dFlexi
     module procedure sendAttributePointerDouble1dFlexi
     module procedure sendAttributePointerDouble2dFlexi
     module procedure sendAttributePointerDouble3dFlexi
     module procedure sendAttributePointerReal2dFlexi
  end interface

  interface sendAttributeStaticFlexi
     module procedure sendAttributeStaticIntegerSingleFlexi
     module procedure sendAttributeStaticRealSingleFlexi
     module procedure sendAttributeStaticCharacterSingleFlexi
     module procedure sendAttributeStaticVectorSingleFlexi
     module procedure sendAttributeStaticInteger1DFlexi
     module procedure sendAttributeStaticLogical1DFlexi
     module procedure sendAttributeStaticVector1DFlexi
     module procedure sendAttributeStaticDouble1dFlexi
     module procedure sendAttributeStaticReal1dFlexi
     module procedure sendAttributeStaticDoubleSingleFlexi
     module procedure sendAttributeStaticLogicalSingleFlexi
  end interface

#endif

contains


  subroutine writeAMRgrid(rootfilename,fileFormatted,grid)
#ifdef MPI
    use mpi
#endif
    use inputs_mod, only : iModel
    use utils_mod, only : findMultiFilename
    character(len=*) :: rootfilename
    character(len=80) :: filename
    logical :: fileFormatted
    type(GRIDTYPE) :: grid
    logical :: writeFile
#ifdef MPI
     integer :: status(MPI_STATUS_SIZE)
    integer :: iThread, ierr
    integer, parameter :: tag = 33
#endif


    call findMultiFilename(rootFilename, iModel, filename)

    writeFile = .true.

#ifdef MPI
    if (.not.grid%splitOverMPI) then
       if (myrankGlobal == 0) then
          writeFile = .true.
       else
          writeFile = .false.
       endif
    endif
#endif


    if (.not.grid%splitOverMPI) then
       if (writeFile) call writeAMRgridSingle(filename, fileFormatted, grid)
       goto 666
    endif
    
#ifdef MPI
    if (grid%splitOverMPI) then

       if (myHydroSetGlobal /= 0) goto 666

       do iThread = 0, nHydroThreadsGlobal

         if (iThread == myRankGlobal) then
 
          if (.not.uncompressedDumpFiles) then
             if (iThread > 0) then
                call MPI_RECV(nBuffer, 1, MPI_INTEGER, ithread-1, tag, localWorldCommunicator, status, ierr)
                call MPI_RECV(buffer, maxBuffer, MPI_BYTE, ithread-1, tag, localWorldCommunicator, status, ierr)
             endif
          endif

          call writeAmrGridSingle(filename, fileFormatted, grid)

          if (.not.uncompressedDumpFiles) then
             if (iThread < nHydroThreadsGlobal) then
                call MPI_SEND(nBuffer, 1, MPI_INTEGER, iThread+1, tag, localWorldCommunicator, ierr)
                call MPI_SEND(buffer, maxBuffer, MPI_BYTE, iThread+1, tag, localWorldCommunicator, ierr)
             endif
          endif

          endif
          
          if (uncompressedDumpFiles) call MPI_BARRIER(localWorldCommunicator, ierr)
       enddo
       if (.not.uncompressedDumpFiles) call zeroZlibBuffer()

    endif
#endif


#ifdef MPI
    call MPI_BARRIER(localWorldCommunicator, ierr)
#endif

666 continue


  end subroutine writeAMRgrid

  subroutine writeAMRgridSingle(filename,fileFormatted,grid)
    ! writes out the 'grid' for an adaptive mesh geometry  

    implicit none
  
    character(len=*)           :: filename
    character(len=80)          :: updatedFilename
    logical, intent(in)        :: fileFormatted
    type(GRIDTYPE), intent(in) :: grid
    
    integer, dimension(8) :: timeValues ! system date and time
    integer               :: error      ! error code
    logical :: writeHeader
    character(len=20) :: positionStatus

    updatedFilename = filename


    writeHeader = .true.
    positionStatus = "rewind"

#ifdef MPI
    if (grid%splitOverMPI) then
       if (myrankGlobal == 0) then 
          writeHeader = .true.
          positionStatus = "newfile"
       else
          writeHeader = .false.
          positionStatus = "append"
       endif
    endif
#endif
          

    if (fileFormatted) then 
       call openUncompressedFile(20, updatedFilename)
    else 
       if (uncompressedDumpFiles) then
          call openUncompressedFile(20, updatedFilename)
       else
#ifdef USEZLIB
!          write(*,*) myrankGlobal, " attempting to open ",trim(updatedFilename)
          call openCompressedFile(20, updatedFilename, positionStatus=positionStatus)
!          write(*,*) myrankGlobal, " has opened ",trim(updatedFilename)
#else
          call writeFatal("zlib is needed to read compressed files")
#endif
       endif
    end if
    if (uncompressedDumpFiles) then
       call writeInfo("Writing AMR grid file to: "//trim(filename),TRIVIAL)
    else
       call writeInfo("Writing compressed AMR grid file to: "//trim(filename),TRIVIAL)
    endif
    
    if (writeHeader) then
       call date_and_time(values=timeValues)
       if (fileFormatted) then
          write(unit=20,fmt=*,iostat=error) timeValues(:)
       else
          write(unit=20,iostat=error) timeValues(:)
       endif

       call writeFileTag(20, "GRIDBEGINS", fileFormatted)
       call writeAttributeStaticFlexi(20,"version", grid%version, fileFormatted)
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
    endif
    call writeOctreePrivateFlexi(grid%octreeRoot,fileFormatted, grid)
    if (uncompressedDumpFiles) then
       close(unit=20)
    else
#ifdef USEZLIB
       if (.not.grid%splitOverMpi) then
          call closeCompressedFile(20,flushBuffer=.true.)
       else
          if (myrankGlobal == nHydroThreadsGlobal) then
             call closeCompressedFile(20,flushBuffer=.true.)
          else
             call closeCompressedFile(20,flushBuffer=.false.)
          endif
       endif
#else
       call writeFatal("zlib is needed for writing compressed files")
#endif
    endif
    call writeInfo("File written and closed",TRIVIAL)
       
    
  contains
  
! Conterpart to openCompressedFile for uncompressed files
    subroutine openUncompressedFile(lunit, thisFilename)

      integer, intent(in) :: lunit
      character(len=*), intent(in) :: thisFilename
      logical :: fileExists
      logical :: luOpened
      
      if (trim(positionStatus)=="append") then
         open(unit=lunit, file=thisFilename, form="unformatted", status="unknown",position=positionStatus)
      else if (trim(positionStatus)=="newfile") then
         open(unit=lunit, file=thisFilename, form="unformatted", status="replace")
      endif

    end subroutine openUncompressedFile

    recursive subroutine writeOctreePrivateFlexi(thisOctal,fileFormatted, grid)
       ! writes out an octal from the grid octree

       type(octal), intent(in), target :: thisOctal
       logical, intent(in)             :: fileFormatted
       type(gridtype) :: grid
       type(octal), pointer :: thisChild
       integer              :: iChild
       logical :: writeThisOctal
       integer :: tempNChildren
       integer :: tempIndexChild(8)
       logical :: tempHasChild(8)

#ifdef MPI
       integer :: i
       type(OCTAL), pointer :: thisOctalPointer
#endif

       writeThisOctal = .true.

       tempNChildren = thisOctal%nChildren
       tempIndexChild = thisOctal%indexChild
       tempHasChild = thisOctal%hasChild


#ifdef MPI


       if (grid%splitOverMPI) then ! this section chooses whether we need to write the octal

          if (nHydroThreadsGlobal == 4) then
             if (myrankGlobal == 0) then 
                if (thisOctal%nDepth == 1) then
                   writeThisOctal = .true.
                else
                   writeThisOctal = .false.
                endif
             else
                if (thisOctal%nDepth == 1) then
                   writeThisOctal = .false.
                else
                   thisOctalPointer => thisOctal
                   writeThisOctal = octalOnThread(thisOctalPointer, 1, myRankGlobal)
                endif
             endif
          endif

          if (nHydroThreadsGlobal == 8) then
             if (myrankGlobal == 0) then 
                if (thisOctal%nDepth == 1) then
                   writeThisOctal = .true.
                else
                   writeThisOctal = .false.
                endif
             else
                if (thisOctal%nDepth == 1) then
                   writeThisOctal = .false.
                else
                   thisOctalPointer => thisOctal
                   writeThisOctal = octalOnThread(thisOctalPointer, 1, myRankGlobal)
                endif
             endif
          endif

          if (nHydroThreadsGlobal == 16) then
             writeThisOctal = .false.
             if (myrankGlobal == 0) then 
                writeThisOctal = .false.
             else
                if ((thisOctal%nDepth == 1).and.(myRankGlobal == 1)) then
                   writeThisOctal = .true.
                endif

                if (thisOctal%nDepth == 2) then
                   thisOctalPointer => thisOctal
                   writeThisOctal = octalOnThread(thisOctalPointer, 1, myRankGlobal)
                endif


                if (thisOctal%nDepth >= 3) then
                   thisOctalPointer => thisOctal
                   writeThisOctal = octalOnThread(thisOctalPointer, 1, myRankGlobal)
                endif
             endif
          endif

          if (nHydroThreadsGlobal == 64) then
             writeThisOctal = .false.
             if (myrankGlobal == 0) then 
                writeThisOctal = .false.
             else
                if ((thisOctal%nDepth == 1).and.(myRankGlobal == 1)) then
                   writeThisOctal = .true.
                endif

                if (thisOctal%nDepth == 2) then
                   thisOctalPointer => thisOctal
                   writeThisOctal = octalOnThread(thisOctalPointer, 1, myRankGlobal)
                endif


                if (thisOctal%nDepth >= 3) then
                   thisOctalPointer => thisOctal
                   writeThisOctal = octalOnThread(thisOctalPointer, 1, myRankGlobal)
                endif
             endif
          endif

          if (nHydroThreadsGlobal == 512) then
             writeThisOctal = .false.
             if (myrankGlobal == 0) then 
                writeThisOctal = .false.
             else
                if ((thisOctal%nDepth == 1).and.(myRankGlobal == 1)) then
                   writeThisOctal = .true.
                endif

                if ((thisOctal%nDepth == 2).and.(myrankGlobal==1).and.(thisOctal%parentSubcell==1)) then
                   writeThisOctal = .true.
                endif

                if ((thisOctal%nDepth == 2).and.(myrankGlobal==65).and.(thisOctal%parentSubcell==2)) then
                   writeThisOctal = .true.
                endif

                if ((thisOctal%nDepth == 2).and.(myrankGlobal==129).and.(thisOctal%parentSubcell==3)) then
                   writeThisOctal = .true.
                endif

                if ((thisOctal%nDepth == 2).and.(myrankGlobal==193).and.(thisOctal%parentSubcell==4)) then
                   writeThisOctal = .true.
                endif

                if ((thisOctal%nDepth == 2).and.(myrankGlobal==257).and.(thisOctal%parentSubcell==5)) then
                   writeThisOctal = .true.
                endif

                if ((thisOctal%nDepth == 2).and.(myrankGlobal==321).and.(thisOctal%parentSubcell==6)) then
                   writeThisOctal = .true.
                endif

                if ((thisOctal%nDepth == 2).and.(myrankGlobal==385).and.(thisOctal%parentSubcell==7)) then
                   writeThisOctal = .true.
                endif

                if ((thisOctal%nDepth == 2).and.(myrankGlobal==449).and.(thisOctal%parentSubcell==8)) then
                   writeThisOctal = .true.
                endif



                if (thisOctal%nDepth >= 3) then
                   thisOctalPointer => thisOctal
                   writeThisOctal = octalOnThread(thisOctalPointer, 1, myRankGlobal)
                endif
             endif
          endif

       endif

       tempNChildren = thisOctal%nChildren
       tempIndexChild = thisOctal%indexChild
       tempHasChild = thisOctal%hasChild

       if (grid%splitOverMpi) then ! this sets up the number of children for the zeroth thread

          if (nHydroThreadsGlobal == 4) then
             if ((thisOctal%nDepth == 1) .and. (myRankGlobal==0)) then
                tempNChildren = 4
                do i = 1, 4
                   tempIndexChild(i) = i
                enddo
                tempHasChild = .true.
             endif
          endif


          if (nHydroThreadsGlobal == 8) then
             if ((thisOctal%nDepth == 1) .and. (myRankGlobal==0)) then
                tempNChildren = 8
                do i = 1, 8
                   tempIndexChild(i) = i
                enddo
                tempHasChild = .true.
             endif
          endif

          if (nHydroThreadsGlobal == 16) then
             if (thisOctal%nDepth == 2) then
                tempNChildren = 4
                do i = 1, 4
                   tempIndexChild(i) = i
                enddo
                tempHasChild = .true.
             endif
          endif


          if (nHydroThreadsGlobal == 64) then
             if (thisOctal%nDepth == 2) then
                tempNChildren = 8
                do i = 1, 8
                   tempIndexChild(i) = i
                enddo
                tempHasChild = .true.
             endif
          endif

          if (nHydroThreadsGlobal == 512) then
             if (thisOctal%nDepth == 2) then
                tempNChildren = 8
                do i = 1, 8
                   tempIndexChild(i) = i
                enddo
                tempHasChild = .true.
             endif

             if (thisOctal%nDepth == 3) then
                tempNChildren = 8
                do i = 1, 8
                   tempIndexChild(i) = i
                enddo
                tempHasChild = .true.
             endif
          endif


       endif
#endif


       if (writeThisOctal) then

          call writeFileTag(20, "OCTALBEGINS", fileFormatted)
          call writeAttributeStaticFlexi(20, "nDepth", thisOctal%nDepth, fileFormatted)
          call writeAttributeStaticFlexi(20, "nChildren", tempNchildren, fileFormatted)
          call writeAttributeStaticFlexi(20, "indexChild", tempIndexChild, fileFormatted)
          call writeAttributeStaticFlexi(20, "hasChild", tempHasChild, fileFormatted)
          call writeAttributeStaticFlexi(20, "centre", thisOctal%centre, fileFormatted)
          call writeAttributeStaticFlexi(20, "rho", thisOctal%rho, fileFormatted)
          call writeAttributeStaticFlexi(20, "temperature", thisOctal%temperature, fileFormatted)
          call writeAttributeStaticFlexi(20, "label", thisOctal%label, fileFormatted)
          call writeAttributeStaticFlexi(20, "subcellSize", thisOctal%subcellSize, fileFormatted)
          call writeAttributeStaticFlexi(20, "threeD", thisOctal%threeD, fileFormatted)
          call writeAttributeStaticFlexi(20, "twoD", thisOctal%twoD, fileFormatted)
          call writeAttributeStaticFlexi(20, "oneD", thisOctal%oneD, fileFormatted)
          call writeAttributeStaticFlexi(20, "maxChildren", thisOctal%maxChildren, fileFormatted)
          call writeAttributeStaticFlexi(20, "cylindrical", thisOctal%cylindrical, fileFormatted)
          call writeAttributeStaticFlexi(20, "splitAzimuthally", thisOctal%splitAzimuthally, fileFormatted)
          call writeAttributeStaticFlexi(20, "phi", thisOctal%phi, fileFormatted)
          call writeAttributeStaticFlexi(20, "dphi", thisOctal%dphi, fileFormatted)
          call writeAttributeStaticFlexi(20, "phimin", thisOctal%phiMin, fileFormatted)
          call writeAttributeStaticFlexi(20, "phimax", thisOctal%phiMax, fileFormatted)


          call writeAttributeStaticFlexi(20, "r", thisOctal%r, fileFormatted)
          call writeAttributeStaticFlexi(20, "parentSubcell", thisOctal%parentSubcell, fileFormatted)
          call writeAttributeStaticFlexi(20, "inStar", thisOctal%inStar, fileFormatted)
          call writeAttributeStaticFlexi(20, "inFlow", thisOctal%inFlow, fileFormatted)
          call writeAttributeStaticFlexi(20, "velocity", thisOctal%velocity, fileFormatted)
          
          call writeAttributeStaticFlexi(20, "xMax", thisOctal%xMax, fileFormatted)
          call writeAttributeStaticFlexi(20, "yMax", thisOctal%yMax, fileFormatted)
          call writeAttributeStaticFlexi(20, "zMax", thisOctal%zMax, fileFormatted)
          call writeAttributeStaticFlexi(20, "xMin", thisOctal%xMin, fileFormatted)
          call writeAttributeStaticFlexi(20, "yMin", thisOctal%yMin, fileFormatted)
          call writeAttributeStaticFlexi(20, "zMin", thisOctal%zMin, fileFormatted)
          
          call writeAttributePointerFlexi(20, "cornervelocity", thisOctal%cornervelocity, fileFormatted)
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
          
          call writeAttributePointerFlexi(20, "distanceGrid", thisOctal%distanceGrid, fileFormatted)
          
          call writeAttributePointerFlexi(20, "nCrossings", thisOctal%nCrossings, fileFormatted)
          call writeAttributePointerFlexi(20, "hHeating", thisOctal%hHeating, fileFormatted)
          call writeAttributePointerFlexi(20, "heHeating", thisOctal%heHeating, fileFormatted)
          call writeAttributePointerFlexi(20, "undersampled", thisOctal%undersampled, fileFormatted)
          call writeAttributePointerFlexi(20, "nDiffusion", thisOctal%nDiffusion, fileFormatted)
          call writeAttributePointerFlexi(20, "diffusionApprox", thisOctal%diffusionApprox, fileFormatted)
          
          call writeAttributePointerFlexi(20, "molecularLevel", thisOctal%molecularLevel, fileFormatted)
          call writeAttributePointerFlexi(20, "jnu", thisOctal%jnu, fileFormatted)
          call writeAttributePointerFlexi(20, "bnu", thisOctal%bnu, fileFormatted)
          call writeAttributePointerFlexi(20, "molAbundance", thisOctal%molAbundance, fileFormatted)
          call writeAttributePointerFlexi(20, "nh2", thisOctal%nh2, fileFormatted)
          call writeAttributePointerFlexi(20, "microTurb", thisOctal%microTurb, fileFormatted)
          call writeAttributePointerFlexi(20, "cornerrho", thisOctal%cornerrho, fileFormatted)          

          call writeAttributePointerFlexi(20, "N", thisOctal%n, fileFormatted)
          call writeAttributePointerFlexi(20, "departCoeff", thisOctal%departCoeff, fileFormatted)
          call writeAttributePointerFlexi(20, "dustTypeFraction", thisOctal%dustTypeFraction, fileFormatted)
          call writeAttributePointerFlexi(20, "oldFrac", thisOctal%oldFrac, fileFormatted)
          call writeAttributePointerFlexi(20, "scatteredIntensity", thisOctal%scatteredIntensity, fileFormatted)

          call writeAttributePointerFlexi(20, "meanIntensity", thisOctal%meanIntensity, fileFormatted)
          
          
          call writeAttributePointerFlexi(20, "atomAbundance", thisOctal%atomAbundance, fileFormatted)
          call writeAttributePointerFlexi(20, "atomLevel", thisOctal%atomLevel, fileFormatted)
          call writeAttributePointerFlexi(20, "jnuCont", thisOctal%jnuCont, fileFormatted)
          call writeAttributePointerFlexi(20, "jnuLine", thisOctal%jnuLine, fileFormatted)
          
          call writeAttributePointerFlexi(20, "q_i", thisOctal%q_i, fileFormatted)
          call writeAttributePointerFlexi(20, "q_i_plus_1", thisOctal%q_i_plus_1, fileFormatted)
          call writeAttributePointerFlexi(20, "q_i_minus_1", thisOctal%q_i_minus_1, fileFormatted)
          call writeAttributePointerFlexi(20, "q_i_minus_2", thisOctal%q_i_minus_2, fileFormatted)
          
          call writeAttributePointerFlexi(20, "x_i", thisOctal%x_i, fileFormatted)
          call writeAttributePointerFlexi(20, "x_i_plus_1", thisOctal%x_i_plus_1, fileFormatted)
          call writeAttributePointerFlexi(20, "x_i_minus_1", thisOctal%x_i_minus_1, fileFormatted)
          call writeAttributePointerFlexi(20, "x_i_minus_2", thisOctal%x_i_minus_2, fileFormatted)
          
          call writeAttributePointerFlexi(20, "u_interface", thisOctal%u_interface, fileFormatted)
          call writeAttributePointerFlexi(20, "u_i", thisOctal%u_i, fileFormatted)
          call writeAttributePointerFlexi(20, "u_i_plus_1", thisOctal%u_i_plus_1, fileFormatted)
          call writeAttributePointerFlexi(20, "u_i_minus_1", thisOctal%u_i_minus_1, fileFormatted)

          call writeAttributePointerFlexi(20, "flux_i", thisOctal%flux_i, fileFormatted)
          call writeAttributePointerFlexi(20, "flux_i_plus_1", thisOctal%flux_i_plus_1, fileFormatted)
          call writeAttributePointerFlexi(20, "flux_i_minus_1", thisOctal%flux_i_minus_1, fileFormatted)


          call writeAttributePointerFlexi(20, "phiLimit", thisOctal%phiLimit, fileFormatted)
       

          call writeAttributePointerFlexi(20, "ghostCell", thisOctal%ghostCell, fileFormatted)
          call writeAttributePointerFlexi(20, "corner", thisOctal%corner, fileFormatted)
          call writeAttributePointerFlexi(20, "feederCell", thisOctal%feederCell, fileFormatted)
          call writeAttributePointerFlexi(20, "edgeCell", thisOctal%edgeCell, fileFormatted)
          call writeAttributePointerFlexi(20, "refinedLastTime", thisOctal%refinedLastTime, fileFormatted)
          
          call writeAttributePointerFlexi(20, "pressure_i", thisOctal%pressure_i, fileFormatted)
          call writeAttributePointerFlexi(20, "pressure_i_plus_1", thisOctal%pressure_i_plus_1, fileFormatted)
          call writeAttributePointerFlexi(20, "pressure_i_minus_1", thisOctal%pressure_i_minus_1, fileFormatted)
          
          call writeAttributePointerFlexi(20, "rhou", thisOctal%rhou, fileFormatted)
          call writeAttributePointerFlexi(20, "rhov", thisOctal%rhov, fileFormatted)
          call writeAttributePointerFlexi(20, "rhow", thisOctal%rhow, fileFormatted)

          call writeAttributePointerFlexi(20, "rhorv_i_plus_1", thisOctal%rhorv_i_plus_1, fileFormatted)
          call writeAttributePointerFlexi(20, "rhorv_i_minus_1", thisOctal%rhorv_i_minus_1, fileFormatted)

          
          call writeAttributePointerFlexi(20, "rhoe", thisOctal%rhoe, fileFormatted)
          call writeAttributePointerFlexi(20, "rhoeLast", thisOctal%rhoeLastTime, fileFormatted)
          call writeAttributePointerFlexi(20, "energy", thisOctal%energy, fileFormatted)

          call writeAttributePointerFlexi(20, "qViscosity", thisOctal%qViscosity, fileFormatted)
       


          call writeAttributePointerFlexi(20, "phi_i", thisOctal%phi_i, fileFormatted)
          call writeAttributePointerFlexi(20, "phi_gas", thisOctal%phi_gas, fileFormatted)
          call writeAttributePointerFlexi(20, "phi_stars", thisOctal%phi_stars, fileFormatted)
          call writeAttributePointerFlexi(20, "phi_i_plus_1", thisOctal%phi_i_plus_1, fileFormatted)
          call writeAttributePointerFlexi(20, "phi_i_minus_1", thisOctal%phi_i_minus_1, fileFormatted)
          
          call writeAttributePointerFlexi(20, "rho_i_plus_1", thisOctal%rho_i_plus_1, fileFormatted)
          call writeAttributePointerFlexi(20, "rho_i_minus_1", thisOctal%rho_i_minus_1, fileFormatted)

!          print *, " write : thisOctal%boundaryCondition ", thisOctal%boundaryCondition
!          print *, " write : ghostcell ", thisOctal%ghostCell
          call writeAttributePointerFlexi(20, "boundaryCondition", thisOctal%boundaryCondition, fileFormatted)
          call writeAttributePointerFlexi(20, "boundaryPartner", thisOctal%boundaryPartner, fileFormatted)
          call writeAttributePointerFlexi(20, "boundaryCell", thisOctal%boundaryCell, fileFormatted)
          call writeAttributePointerFlexi(20, "radiationMomentum", thisOctal%radiationMomentum, fileFormatted)
          call writeAttributePointerFlexi(20, "kappaTimesFlux", thisOctal%kappaTimesFlux, fileFormatted)
          call writeAttributePointerFlexi(20, "UVvector", thisOctal%UVvector, fileFormatted)

          call writeAttributePointerFlexi(20, "gravboundaryPartner", thisOctal%GravboundaryPartner, fileFormatted)
          call writeAttributePointerFlexi(20, "changed", thisOctal%changed, fileFormatted)
          call writeAttributePointerFlexi(20, "rLimit", thisOctal%rLimit, fileFormatted)
          
          call writeAttributePointerFlexi(20, "iEquationOfState", thisOctal%iEquationOfState, fileFormatted)

          call writeAttributePointerFlexi(20, "iAnalyticalVelocity", thisOctal%iAnalyticalVelocity, fileFormatted)

          call writeAttributePointerFlexi(20, "gamma", thisOctal%gamma, fileFormatted)


          call writeAttributePointerFlexi(20, "udens", thisOctal%udens, fileFormatted)
          call writeAttributePointerFlexi(20, "adot", thisOctal%adot, fileFormatted)
          call writeAttributePointerFlexi(20, "adotdist", thisOctal%distanceGridAdot, fileFormatted)
          call writeAttributePointerFlexi(20, "dfromgas", thisOctal%distanceGridPhotonFromGas, fileFormatted)
          call writeAttributePointerFlexi(20, "dfromsource", thisOctal%distanceGridPhotonFromSource, fileFormatted)
          call writeAttributePointerFlexi(20, "ufromgas", thisOctal%photonEnergyDensityFromGas, fileFormatted)
          call writeAttributePointerFlexi(20, "ufromsource", thisOctal%photonEnergyDensityFromSource, fileFormatted)
          call writeAttributePointerFlexi(20, "utotal", thisOctal%photonEnergyDensity, fileFormatted)
          call writeAttributePointerFlexi(20, "oldutotal", thisOctal%oldphotonEnergyDensity, fileFormatted)

          call writeAttributePointerFlexi(20, "fixedtemperature", thisOctal%fixedTemperature, fileFormatted)

          call writeAttributeStaticFlexi(20, "mpiThread", thisOctal%mpiThread, fileFormatted)
          call writeFileTag(20, "OCTALENDS", fileFormatted)
       endif
       
       if (thisOctal%nChildren > 0) then 
          do iChild = 1, thisOctal%nChildren, 1
             thisChild => thisOctal%child(iChild)
             call writeOctreePrivateFlexi(thisChild,fileFormatted,grid)
          end do
       end if

     end subroutine writeOctreePrivateFlexi
    
    
   end subroutine writeAMRgridSingle


  subroutine readAMRgrid(rootfilename,fileFormatted,grid)
    use inputs_mod, only: gridUsesAMR, iModel, modelwashydro
    use utils_mod, only : findMultiFilename
    type(romanova) :: romdata
    character(len=*) :: rootfilename
    character(len=80) :: filename
    logical :: fileFormatted
    type(GRIDTYPE) :: grid
    logical :: readFile
    logical, save :: firstTime=.true.
    character(len=80) :: message
    integer ::  nOctals, nVoxels
#ifdef MPI
!    integer :: iThread
#endif

    readFile = .true.
    gridUsesAMR = .true. 

    call findMultiFilename(rootFilename, iModel, filename)

    if (associated(grid%octreeRoot)) then
       call deleteOctreeBranch(grid%octreeRoot,onlyChildren=.false., adjustParent=.false.)
       grid%octreeRoot => null()
    endif

    if (.not.grid%splitOverMPI) then


#ifdef MPI
       call readGridMPI(grid, filename, fileFormatted)
#else
       call readAmrGridSingle(filename, fileFormatted, grid)
#endif
             call updateMaxDepth(grid)
             call setSmallestSubcell(grid)
!             call checkAMRgrid(grid, .false.)
             call countVoxels(grid%octreeRoot,nOctals,nVoxels)
             grid%nOctals = nOctals
             
!#ifdef MPI
!       do iThread = 0, nThreadsGlobal-1
!          if (iThread == myRankGlobal) then
!             write(*,'(a,i3)') "Reading grid for thread: ",ithread
!#endif
!             call readAmrGridSingle(filename, fileFormatted, grid)
!             call updateMaxDepth(grid)
!             call setSmallestSubcell(grid)
!!             call checkAMRgrid(grid, .false.)
!             call countVoxels(grid%octreeRoot,nOctals,nVoxels)
!             grid%nOctals = nOctals
!!             CALL checkAMRgrid(grid,checkNoctals=.FALSE.)
!#ifdef MPI
!          endif
!          call torus_mpi_barrier
!       Enddo
!#endif
    endif
    
#ifdef MPI
    if (grid%splitOverMPI) then
       call readGridSplitOverMPI(grid, filename, fileFormatted)
!       do iThread = 0, nThreadsGlobal-1
!          if (iThread == myRankGlobal) then
!             call readAmrGridSingle(filename, fileFormatted, grid)
!             call updateMaxDepth(grid)
!             call setSmallestSubcell(grid)
!             call checkAMRgrid(grid, .false.)
!             call countVoxels(grid%octreeRoot,nOctals,nVoxels)
!             grid%nOctals = nOctals
!          endif
!       enddo
       call grid_info_mpi(grid, "info_grid.dat")
       call torus_mpi_barrier
    endif
#endif

    if(modelwashydro) then
       if(firstTime) then
          write(message,'(a,i3,a)') "Adding cell corner values"
          
          call writeInfo(message,TRIVIAL)
          firsttime = .false.
       end if
       call finishgrid(grid%octreeRoot, grid, romData=romData)                
    end if

  end subroutine readAMRgrid

  subroutine readAMRgridSingle(filename,fileFormatted,grid)
    ! reads in a previously saved 'grid' for an adaptive mesh geometry  
    use unix_mod, only: unixGetenv
    implicit none

    character(len=*)            :: filename
    logical, intent(in)         :: fileFormatted
    type(GRIDTYPE), intent(inout) :: grid
    
    integer               :: error         ! status code
    integer :: nOctal
    character(len=80) :: absolutePath, inFile, updatedFilename
    integer :: ithread

    type(OCTAL), pointer :: mynull => null()

    absolutePath = " "
    error = 0

  call unixGetEnv("TORUS_JOB_DIR",absolutePath)
  inFile = trim(absolutePath)//trim(filename)


  updatedFilename = inFile

  call openGridFile(updatedFilename, fileformatted)


    call readStaticComponents(grid, fileformatted)
    call writeInfo("Static components read")

       allocate(grid%octreeRoot)
       grid%octreeRoot%nDepth = 1
       nOctal = 0
       mynull => grid%octreeRoot
       call readOctreePrivateFlexi(grid%octreeRoot,mynull,fileFormatted, nOctal, grid)

    call writeInfo("Grid read done")
    
    close(unit=20)

  contains
   
    recursive subroutine readOctreePrivateFlexi(thisOctal,parent,fileFormatted, noctal, grid)
      ! read in an octal to the grid octree

      implicit none
      type(octal), target :: thisOctal
      type(octal), target :: parent
      type(gridtype) :: grid

      logical, intent(in)  :: fileFormatted
      integer :: nOctal

      type(octal), pointer :: childPointer
      type(octal), pointer :: thisChild => null()
      type(octal), pointer :: topOctal => null()
      type(octal), pointer :: tempChildPointer => null()
      type(octal), pointer :: tempChildPointer2 => null()
      integer              :: iChild
      logical :: foundBranch

      

      nOctal = nOctal+1
      thisOctal%parent => parent

      call readOctalViaTags(thisOctal, fileFormatted)

#ifdef MPI

      if (grid%splitOverMPI) then ! this sets the nummber of children correctly for the threads

         if (nHydroThreadsGlobal == 4) then
            if (myrankGlobal == 0) then
               thisOctal%nChildren = 0
            else
               if (thisOctal%nDepth == 1) then
                  thisOctal%nChildren = 1
                  thisOctal%hasChild = .false.
                  thisOctal%indexChild(1) = myRankGlobal
                  thisOctal%hasChild(myRankGlobal) = .true.
               endif
            endif
         endif

         if (nHydroThreadsGlobal == 8) then
            if (myrankGlobal == 0) then
               thisOctal%nChildren = 0
            else
               if (thisOctal%nDepth == 1) then
                  thisOctal%nChildren = 1
                  thisOctal%hasChild = .false.
                  thisOctal%indexChild(1) = myRankGlobal
                  thisOctal%hasChild(myRankGlobal) = .true.
               endif
            endif
         endif

         if (nHydroThreadsGlobal == 16) then
            if (myrankGlobal == 0) then
               if (thisOctal%nDepth == 2) then
                  thisOctal%nChildren = 0
                  thisOctal%hasChild = .false.
               endif
            else
               if (thisOctal%nDepth == 2) then
                  if (.not.octalOnThread(thisOctal%parent, thisOctal%parentSubcell, myrankGlobal)) then
                     thisOctal%nChildren = 0
                     thisOctal%hasChild = .false.
                  else
                     thisOctal%nChildren = 1
                     thisOctal%hasChild = .false.
                     iChild = thisOctal%parentSubcell
                     thisOctal%indexChild(1) = iChild
                     thisOctal%hasChild(iChild) = .true.
                  endif
               endif
            endif
         endif

         if (nHydroThreadsGlobal == 64) then
            if (myrankGlobal == 0) then
               if (thisOctal%nDepth == 2) then
                  thisOctal%nChildren = 0
                  thisOctal%hasChild = .false.
               endif
            else
               if (thisOctal%nDepth == 2) then
                  if (.not.octalOnThread(thisOctal%parent, thisOctal%parentSubcell, myrankGlobal)) then
                     thisOctal%nChildren = 0
                     thisOctal%hasChild = .false.
                  else
                     thisOctal%nChildren = 1
                     thisOctal%hasChild = .false.
                     iChild = thisOctal%parentSubcell
                     thisOctal%indexChild(1) = iChild
                     thisOctal%hasChild(iChild) = .true.
                  endif
               endif
            endif
         endif

         if (nHydroThreadsGlobal == 512) then
            if (myrankGlobal == 0) then
               if (thisOctal%nDepth == 2) then
                  thisOctal%nChildren = 0
                  thisOctal%hasChild = .false.
               endif
            else
               if (thisOctal%nDepth == 3) then
                  if (.not.octalOnThread(thisOctal%parent, thisOctal%parentSubcell, myrankGlobal)) then
                     thisOctal%nChildren = 0
                     thisOctal%hasChild = .false.
                  else
                     thisOctal%nChildren = 1
                     thisOctal%hasChild = .false.
                     iChild = thisOctal%parentSubcell
                     thisOctal%indexChild(1) = iChild
                     thisOctal%hasChild(iChild) = .true.
                  endif
               endif
            endif
         endif


      endif
#endif


      if (.not.grid%splitOVerMPI) then
         if (thisOctal%nChildren > 0) then 
            allocate(thisOctal%child(1:thisOctal%nChildren)) 
            do iChild = 1, thisOctal%nChildren, 1
               thisChild => thisOctal%child(iChild)
               call readOctreePrivateFlexi(thisChild,thisOctal,fileFormatted, nOctal, grid)               
            end do
         end if
      endif

      if (grid%splitOVerMPI) then

         if (nHydroThreadsGlobal == 4) then
            if (thisOctal%nDepth > 1) then
               if (thisOctal%nChildren > 0) then 
                  allocate(thisOctal%child(1:thisOctal%nChildren)) 
                  do iChild = 1, thisOctal%nChildren, 1
                     thisChild => thisOctal%child(iChild)
                     call readOctreePrivateFlexi(thisChild,thisOctal,fileFormatted, nOctal, grid)               
                  end do
               end if
            else
               if (myrankGlobal == 0) then
                  thisOctal%nChildren = 0
                  thisOctal%hasChild = .false.
               else
                  if (thisOctal%nChildren > 0) then 
                     allocate(thisOctal%child(1:thisOctal%nChildren)) 
                     do iChild = 1, thisOctal%nChildren, 1
                        do iThread = 1, myRankGlobal
                           thisChild => thisOctal%child(iChild)
                           topOctal => thisChild
                           call readOctreePrivateFlexi(thisChild,thisOctal,fileFormatted, nOctal, grid)               
                           if (iThread /= myRankGlobal) then
                              call deleteOctreeBranch(topOctal,onlyChildren=.false., adjustParent=.false.)
                           else
                              exit
                           endif
                        enddo
                     end do
                  end if
               endif
            endif
         endif

         if (nHydroThreadsGlobal == 8) then
            if (thisOctal%nDepth > 1) then
               if (thisOctal%nChildren > 0) then 
                  allocate(thisOctal%child(1:thisOctal%nChildren)) 
                  do iChild = 1, thisOctal%nChildren, 1
                     write(*,*) "myrank ",myrankglobal, " ichild ",ichild, size(thisOctal%child),thisOctal%nchildren
                     thisChild => thisOctal%child(iChild)
                     call readOctreePrivateFlexi(thisChild,thisOctal,fileFormatted, nOctal, grid)               
                  end do
               end if
            else
               if (myrankGlobal == 0) then
                  thisOctal%nChildren = 0
                  thisOctal%hasChild = .false.
               else
                  if (thisOctal%nChildren > 0) then 
                     allocate(thisOctal%child(1:thisOctal%nChildren)) 
                     do iChild = 1, thisOctal%nChildren, 1
                        do iThread = 1, myRankGlobal
                           thisChild => thisOctal%child(iChild)
                           topOctal => thisChild
                           call readOctreePrivateFlexi(thisChild,thisOctal,fileFormatted, nOctal, grid)               
                           if (iThread /= myRankGlobal) then
                              call deleteOctreeBranch(topOctal,onlyChildren=.false., adjustParent=.false.)
                           else
                              exit
                           endif
                        enddo
                     end do
                  end if
               endif
            endif
         endif

         if (nHydroThreadsGlobal == 16) then
            if (thisOctal%nDepth > 2) then
               if (thisOctal%nChildren > 0) then 
                  allocate(thisOctal%child(1:thisOctal%nChildren)) 
                  do iChild = 1, thisOctal%nChildren, 1
                     tempChildPointer2 => thisOctal%child(iChild)
                     call readOctreePrivateFlexi(tempChildPointer2,thisOctal,fileFormatted, nOctal, grid)               
                  end do
               end if
            else if (thisOctal%nDepth ==   1) then 
               allocate(thisOctal%child(1:thisOctal%nChildren)) 
               do iChild = 1, thisOctal%nChildren, 1
                  thisChild => thisOctal%child(iChild)
                  call readOctreePrivateFlexi(thisChild,thisOctal,fileFormatted, nOctal, grid)               
               end do
            else  if (thisOctal%nDepth ==  2) then 
               if (thisOctal%nChildren == 1) then
                  allocate(thisOctal%child(1))
               endif
               foundBranch = .false.
               do iChild = 1, 4


                  allocate(grid%tempBranch)
                  tempChildPointer => grid%tempBranch

                  call readOctreePrivateFlexi(tempChildPointer,thisOctal,fileFormatted, nOctal, grid)               

                  if (thisOctal%mpiThread(iChild) /= myRankGlobal) then
                     call deleteOctreeBranch(tempChildPointer,onlyChildren=.true., adjustParent=.false.)
                     deallocate(grid%tempBranch)
                     grid%tempBranch => null()
!                     call skipOctalsToDepth(fileformatted, 2)

                  else
                     childPointer => thisOctal%child(1) ! TJH 9 JULY
                     call insertOctreeBranch(childPointer , grid%tempBranch, onlyChildren = .false.)
!                     call insertOctreeBranch(thisOctal%child(1), grid%tempBranch, onlyChildren = .false.)

!                     call readOctreePrivateFlexi(thisOctal%child(1),thisOctal,fileFormatted, nOctal, grid)               
                     thisOctal%hasChild = .false.
                     thisOctal%hasChild(iChild) = .true.
                     thisOctal%indexChild(1) = iChild
                     thisOctal%child(1)%parent => thisOctal
                     thisOctal%child(1)%parentSubcell = iChild
                  endif

               end do
            end if
         endif

         if (nHydroThreadsGlobal == 64) then
            if (thisOctal%nDepth > 2) then
               if (thisOctal%nChildren > 0) then 
                  allocate(thisOctal%child(1:thisOctal%nChildren)) 
                  do iChild = 1, thisOctal%nChildren, 1
                     tempChildPointer2 => thisOctal%child(iChild)
                     call readOctreePrivateFlexi(tempChildPointer2,thisOctal,fileFormatted, nOctal, grid)               
                  end do
               end if
            else if (thisOctal%nDepth ==   1) then 
               allocate(thisOctal%child(1:thisOctal%nChildren)) 
               do iChild = 1, thisOctal%nChildren, 1
                  thisChild => thisOctal%child(iChild)
                  call readOctreePrivateFlexi(thisChild,thisOctal,fileFormatted, nOctal, grid)               
               end do
            else  if (thisOctal%nDepth ==  2) then 
               if (thisOctal%nChildren == 1) then
                  allocate(thisOctal%child(1))
               endif
               foundBranch = .false.
               do iChild = 1, 8


                  allocate(grid%tempBranch)
                  tempChildPointer => grid%tempBranch

                  call readOctreePrivateFlexi(tempChildPointer,thisOctal,fileFormatted, nOctal, grid)               

                  if (thisOctal%mpiThread(iChild) /= myRankGlobal) then
                     call deleteOctreeBranch(tempChildPointer,onlyChildren=.true., adjustParent=.false.)
                     deallocate(grid%tempBranch)
                     grid%tempBranch => null()
!                     call skipOctalsToDepth(fileformatted, 2)

                  else
                     childPointer => thisOctal%child(1) ! TJH 9 JULY
                     call insertOctreeBranch(childPointer , grid%tempBranch, onlyChildren = .false.)
!                     call insertOctreeBranch(thisOctal%child(1), grid%tempBranch, onlyChildren = .false.)

!                     call readOctreePrivateFlexi(thisOctal%child(1),thisOctal,fileFormatted, nOctal, grid)               
                     thisOctal%hasChild = .false.
                     thisOctal%hasChild(iChild) = .true.
                     thisOctal%indexChild(1) = iChild
                     thisOctal%child(1)%parent => thisOctal
                     thisOctal%child(1)%parentSubcell = iChild
                  endif

               end do
            end if
         endif

         if (nHydroThreadsGlobal == 512) then
            if (thisOctal%nDepth > 3) then
               if (thisOctal%nChildren > 0) then 
                  allocate(thisOctal%child(1:thisOctal%nChildren)) 
                  do iChild = 1, thisOctal%nChildren, 1
                     tempChildPointer2 => thisOctal%child(iChild)
                     call readOctreePrivateFlexi(tempChildPointer2,thisOctal,fileFormatted, nOctal, grid)               
                  end do
               end if
            else if (thisOctal%nDepth ==   1) then 
               allocate(thisOctal%child(1:thisOctal%nChildren)) 
               do iChild = 1, thisOctal%nChildren, 1
                  thisChild => thisOctal%child(iChild)
                  call readOctreePrivateFlexi(thisChild,thisOctal,fileFormatted, nOctal, grid)               
               end do
            else if (thisOctal%nDepth ==   2) then 
               allocate(thisOctal%child(1:thisOctal%nChildren)) 
               do iChild = 1, thisOctal%nChildren, 1
                  thisChild => thisOctal%child(iChild)
                  call readOctreePrivateFlexi(thisChild,thisOctal,fileFormatted, nOctal, grid)               
               end do
            else  if (thisOctal%nDepth ==  3) then 
               if (thisOctal%nChildren == 1) then
                  allocate(thisOctal%child(1))
               endif
               foundBranch = .false.
               do iChild = 1, 8


                  allocate(grid%tempBranch)
                  tempChildPointer => grid%tempBranch

                  call readOctreePrivateFlexi(tempChildPointer,thisOctal,fileFormatted, nOctal, grid)               

                  if (thisOctal%mpiThread(iChild) /= myRankGlobal) then
                     call deleteOctreeBranch(tempChildPointer,onlyChildren=.true., adjustParent=.false.)
                     deallocate(grid%tempBranch)
                     grid%tempBranch => null()
!                     call skipOctalsToDepth(fileformatted, 2)

                  else
                     childPointer => thisOctal%child(1) ! TJH 9 JULY
                     call insertOctreeBranch(childPointer , grid%tempBranch, onlyChildren = .false.)
!                     call insertOctreeBranch(thisOctal%child(1), grid%tempBranch, onlyChildren = .false.)

!                     call readOctreePrivateFlexi(thisOctal%child(1),thisOctal,fileFormatted, nOctal, grid)               
                     thisOctal%hasChild = .false.
                     thisOctal%hasChild(iChild) = .true.
                     thisOctal%indexChild(1) = iChild
                     thisOctal%child(1)%parent => thisOctal
                     thisOctal%child(1)%parentSubcell = iChild
                  endif

               end do
            end if
         endif



      endif


 end subroutine readOctreePrivateFlexi
 
  end subroutine readAMRgridSingle




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
       if (uncompressedDumpFiles) then
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) value
       else
          call writeCompressedFile(lunit, attributeName)
          call writeCompressedFile(lunit, dataType)
          call writeCompressedFile(lunit, value)
       endif
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
       if (uncompressedDumpFiles) then
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) value
       else
          call writeCompressedFile(lunit, attributeName)
          call writeCompressedFile(lunit, dataType)
          call writeCompressedFile(lunit, value)
       endif
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
       if (uncompressedDumpFiles) then
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) len(trim(value))
          write(lUnit) trim(value)
       else
          call writeCompressedFile(lunit, attributeName)
          call writeCompressedFile(lunit, dataType)
          call writeCompressedFile(lunit, len(trim(value)))
          call writeCompressedFile(lunit, trim(value))
       endif
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
       if (uncompressedDumpFiles) then
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) size(value)
          write(lUnit) value(1:size(value))
       else
          call writeCompressedFile(lunit, attributeName)
          call writeCompressedFile(lunit, dataType)
          call writeCompressedFile(lunit, size(value))
          call writeCompressedFile(lunit, value(1:size(value)))
       endif
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
       if (uncompressedDumpFiles) then
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) size(value)
          write(lUnit) value(1:size(value))
       else
          call writeCompressedFile(lunit, attributeName)
          call writeCompressedFile(lunit, dataType)
          call writeCompressedFile(lunit, size(value))
          call writeCompressedFile(lunit, value(1:size(value)))
       endif
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
       if (uncompressedDumpFiles) then
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) size(value)
          write(lUnit) value(1:size(value))
       else
          call writeCompressedFile(lunit, attributeName)
          call writeCompressedFile(lunit, dataType)
          call writeCompressedFile(lunit, size(value))
          call writeCompressedFile(lunit, value(1:size(value)))
       endif
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
       if (uncompressedDumpFiles) then
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) value
       else
          call writeCompressedFile(lunit, attributeName)
          call writeCompressedFile(lunit, dataType)
          call writeCompressedFile(lunit, value)
       endif
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
       if (uncompressedDumpFiles) then
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) value
       else
          call writeCompressedFile(lunit, attributeName)
          call writeCompressedFile(lunit, dataType)
          call writeCompressedFile(lunit, value)
       endif
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
       if (uncompressedDumpFiles) then
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) size(value)
          write(lUnit) value(1:size(value))
       else
          call writeCompressedFile(lunit, attributeName)
          call writeCompressedFile(lunit, dataType)
          call writeCompressedFile(lunit, size(value))
          call writeCompressedFile(lunit, value(1:size(value)))
       endif
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
       if (uncompressedDumpFiles) then
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) size(value)
          write(lUnit) value(1:size(value))
       else
          call writeCompressedFile(lunit, attributeName)
          call writeCompressedFile(lunit, dataType)
          call writeCompressedFile(lunit, size(value))
          call writeCompressedFile(lunit, value(1:size(value)))
       endif
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
       if (uncompressedDumpFiles) then
          write(lUnit) attributeName
          write(lUnit) dataType
          write(lUnit) value
       else
          call writeCompressedFile(lunit, attributeName)
          call writeCompressedFile(lunit, dataType)
          call writeCompressedFile(lunit, value)
       endif
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
          if (uncompressedDumpFiles) then
             write(lUnit) attributeName
             write(lUnit) dataType
             write(lUnit) size(value)
             write(lUnit) value(1:SIZE(value))
          else
             call writeCompressedFile(lunit, attributeName)
             call writeCompressedFile(lunit, dataType)
             call writeCompressedFile(lunit, size(value))
             call writeCompressedFile(lunit, value(1:SIZE(value)))
          endif
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
          if (uncompressedDumpFiles) then
             write(lUnit) attributeName
             write(lUnit) dataType
             write(lUnit) SIZE(value)
             write(lUnit) value(1:SIZE(value))
          else
             call writeCompressedFile(lunit, attributeName)
             call writeCompressedFile(lunit, dataType)
             call writeCompressedFile(lunit, SIZE(value))
             call writeCompressedFile(lunit, value(1:SIZE(value)))
          endif
       endif
    endif
  end subroutine writeAttributePointerInteger1DFlexi


  subroutine writeAttributePointerLogical1DFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    logical, pointer :: value(:)
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "l1darray"

    attributeName = name
    if (associated(value)) then
       if (fileFormatted) then
          write(lUnit,*) attributeName
          write(lUnit,*) dataType
          write(lUnit,*) SIZE(value)
          write(lUnit,*) value(1:SIZE(value))
       else
          if (uncompressedDumpFiles) then
             write(lUnit) attributeName
             write(lUnit) dataType
             write(lUnit) SIZE(value)
             write(lUnit) value(1:SIZE(value))
          else
             call writeCompressedFile(lunit, attributeName)
             call writeCompressedFile(lunit, dataType)
             call writeCompressedFile(lunit, SIZE(value))
             call writeCompressedFile(lunit, value(1:SIZE(value)))
          endif
       endif
    endif
  end subroutine writeAttributePointerLogical1DFlexi

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
          if (uncompressedDumpFiles) then
             write(lUnit) attributeName
             write(lUnit) dataType
             write(lUnit) SIZE(value)
             write(lUnit) value(1:SIZE(value))
          else
             call writeCompressedFile(lunit, attributeName)
             call writeCompressedFile(lunit, dataType)
             call writeCompressedFile(lunit, SIZE(value))
             call writeCompressedFile(lunit, value(1:SIZE(value)))
          endif
       endif
    endif
  end subroutine writeAttributePointerReal1DFlexi

  subroutine writeAttributePointerVector1DFlexi(lUnit, name, value, fileFormatted)
    integer :: lUnit
    character(len=*) :: name
    character(len=20) :: attributeName
    type(VECTOR), pointer :: value(:)
    logical :: fileFormatted
    character(len=10) :: dataType

    dataType = "v1darray"

    attributeName = name
    if (associated(value)) then
       if (fileFormatted) then
          write(lUnit,*) attributeName
          write(lUnit,*) dataType
          write(lUnit,*) SIZE(value)
          write(lUnit,*) value(1:SIZE(value))
       else
          if (uncompressedDumpFiles) then
             write(lUnit) attributeName
             write(lUnit) dataType
             write(lUnit) SIZE(value)
             write(lUnit) value(1:SIZE(value))
          else
             call writeCompressedFile(lunit, attributeName)
             call writeCompressedFile(lunit, dataType)
             call writeCompressedFile(lunit, SIZE(value))
             call writeCompressedFile(lunit, value(1:SIZE(value)))
          endif
       endif
    endif
  end subroutine writeAttributePointerVector1DFlexi

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
          if (uncompressedDumpFiles) then
             write(lUnit) attributeName
             write(lUnit) dataType
             write(lUnit) SIZE(value,1),SIZE(value,2)
             write(lUnit) value(1:SIZE(value,1),1:SIZE(value,2))
          else
             call writeCompressedFile(lunit, attributeName)
             call writeCompressedFile(lunit, dataType)
             call writeCompressedFile(lunit, SIZE(value,1))
             call writeCompressedFile(lunit, SIZE(value,2))
             call writeCompressedFile(lunit, value(1:SIZE(value,1),1:SIZE(value,2)))
          endif
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
          if (uncompressedDumpFiles) then
             write(lUnit) attributeName
             write(lUnit) dataType
             write(lUnit) SIZE(value,1),SIZE(value,2), SIZE(value,3)
             write(lUnit) value(1:SIZE(value,1),1:SIZE(value,2),1:SIZE(value,3))
          else
             call writeCompressedFile(lunit, attributeName)
             call writeCompressedFile(lunit, dataType)
             call writeCompressedFile(lunit, SIZE(value,1))
             call writeCompressedFile(lunit, SIZE(value,2))
             call writeCompressedFile(lunit, SIZE(value,3))
             call writeCompressedFile(lunit, value(1:SIZE(value,1),1:SIZE(value,2),1:SIZE(value,3)))
          endif
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
          if (uncompressedDumpFiles) then
             write(lUnit) attributeName
             write(lUnit) dataType
             write(lUnit) SIZE(value,1),SIZE(value,2)
             write(lUnit) value(1:SIZE(value,1),1:SIZE(value,2))
          else
             call writeCompressedFile(lunit, attributeName)
             call writeCompressedFile(lunit, dataType)
             call writeCompressedFile(lunit, SIZE(value,1))
             call writeCompressedFile(lunit, SIZE(value,2))
             call writeCompressedFile(lunit, value(1:SIZE(value,1),1:SIZE(value,2)))
          endif
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
         if (uncompressedDumpFiles) then
            write(lUnit) wtag
         else
            call writeCompressedFile(lUnit, wTag)
         endif
      endif
    end subroutine writeFileTag

#ifdef MPI

    subroutine receiveSingleIntegerFlexi(value)
      use mpi
      integer, intent(out) :: value
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      call MPI_RECV(value, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)

    end subroutine receiveSingleIntegerFlexi

    subroutine receiveSingleRealFlexi(value)
      use mpi
      real, intent(out) :: value
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      call MPI_RECV(value, 1, MPI_REAL, 0, tag, localWorldCommunicator, status, ierr)

    end subroutine receiveSingleRealFlexi

    subroutine receiveSingleDoubleFlexi(value)
      use mpi
      real(double), intent(out) :: value
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      call MPI_RECV(value, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

    end subroutine receiveSingleDoubleFlexi

    subroutine receiveSingleLogicalFlexi(value)
      use mpi
      logical, intent(out) :: value
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      call MPI_RECV(value, 1, MPI_LOGICAL, 0, tag, localWorldCommunicator, status, ierr)

    end subroutine receiveSingleLogicalFlexi

    subroutine receiveSingleVectorFlexi(value)
      use mpi
      type(VECTOR), intent(out) :: value
      real(double) :: v(3)
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      call MPI_RECV(v, 3, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
      
      value = VECTOR(v(1), v(2), v(3))

    end subroutine receiveSingleVectorFlexi

    subroutine receiveSingleCharacterFlexi(value)
      use mpi
      character(len=*), intent(out) :: value
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr, i

      do i = 1, len(value)
         value(i:i) = " "
      enddo
      call MPI_RECV(i, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      call MPI_RECV(value, i, MPI_CHARACTER, 0, tag, localWorldCommunicator, status, ierr)

    end subroutine receiveSingleCharacterFlexi
#endif    

    subroutine readSingleIntegerFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      integer, intent(out) :: value
      logical :: fileFormatted

      call testDataType("isingle", fileFormatted)
      if (fileFormatted) then
         read(lUnit,*) value
      else
         if (uncompressedDumpFiles) then
            read(lUnit) value
         else
            call readCompressedFile(lUnit, value)
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) value
         else
            call readCompressedFile(lUnit, value)
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) value
         else
            call readCompressedFile(lUnit, value)
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) value
         else
            call readCompressedFile(lUnit, value)
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) value
         else
            call readCompressedFile(lUnit, value)
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) i
            read(lUnit) value(1:i)
         else
            call readCompressedFile(lUnit, i)
            call readCompressedFile(lUnit, value(1:i))
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) n
            allocate(value(1:n))
            read(lUnit) value(1:n)
         else
            call readCompressedFile(lunit, n)
            allocate(value(1:n))
            call readCompressedFile(lunit,value)
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) n, m
            allocate(value(1:n,1:m))
            read(lUnit) value(1:n,1:m)
         else
            call readCompressedFile(lunit, n)
            call readCompressedFile(lunit, m)
            allocate(value(1:n,1:m))
            call readCompressedFile(lunit,value)
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) n, m, j
            allocate(value(1:n,1:m,1:j))
            read(lUnit) value(1:n,1:m,1:j)
         else
            call readCompressedFile(lunit, n)
            call readCompressedFile(lunit, m)
            call readCompressedFile(lunit, j)
            allocate(value(1:n,1:m,1:j))
            call readCompressedFile(lunit,value)
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) n, m
            allocate(value(1:n,1:m))
            read(lUnit) value(1:n,1:m)
         else
            call readCompressedFile(lunit, n)
            call readCompressedFile(lunit, m)
            allocate(value(1:n,1:m))
!            write(*,*) "allocated array ",n,m
            call readCompressedFile(lunit,value)
!            write(*,*) "found ",value(1:n,1:m)
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) n
            allocate(value(1:n))
            read(lUnit) value(1:n)
         else
            call readCompressedFile(lunit, n)
            allocate(value(1:n))
            call readCompressedFile(lunit,value)
         endif
      endif
    end subroutine readRealPointer1DFlexi

    subroutine readLogicalPointer1DFlexi(lUnit, value, fileFormatted)
      integer :: lUnit
      logical, pointer :: value(:)
      logical :: fileFormatted
      integer :: n

      call testDataType("l1darray", fileFormatted)
      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      if (fileFormatted) then
         read(lUnit,*) n
         allocate(value(1:n))
         read(lUnit,*) value(1:n)
      else
         if (uncompressedDumpFiles) then
            read(lUnit) n
            allocate(value(1:n))
            read(lUnit) value(1:n)
         else
            call readCompressedFile(lunit, n)
            allocate(value(1:n))
            call readCompressedFile(lunit,value)
         endif
      endif
    end subroutine readLogicalPointer1DFlexi

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
         if (uncompressedDumpFiles) then
            read(lUnit) n
            allocate(value(1:n))
            read(lUnit) value(1:n)
         else
            call readCompressedFile(lunit, n)
            allocate(value(1:n))
            call readCompressedFile(lunit,value)
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) n
            allocate(value(1:n))
            read(lUnit) value(1:n)
         else
            call readCompressedFile(lunit, n)
            allocate(value(1:n))
            call readCompressedFile(lunit,value)
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) n
            read(lUnit) value(1:n)
         else
            call readCompressedFile(lunit, n)
            call readCompressedFile(lunit,value(1:n))
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) n
            if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
            read(lUnit) value(1:n)
         else
            call readCompressedFile(lunit, n)
            if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
            call readCompressedFile(lunit,value(1:n))
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) n
            if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
            read(lUnit) value(1:n)
         else
            call readCompressedFile(lunit, n)
            if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
            call readCompressedFile(lunit,value(1:n))
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) n
            if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
            read(lUnit) value(1:n)
         else
            call readCompressedFile(lunit, n)
            if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
            call readCompressedFile(lunit,value(1:n))
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) n
            if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
            read(lUnit) value(1:n)
         else
            call readCompressedFile(lunit, n)
            if (n > size(value)) call writeFatal("Error with array sizes in readVectorArrayFlexi")
            call readCompressedFile(lunit,value(1:n))
         endif
      endif
    end subroutine readArrayDoubleFlexi

#ifdef MPI

  subroutine sendAttributePointerInteger1DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    integer, pointer :: value(:)
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name
    if (associated(value)) then

       call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(SIZE(value), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(value, SIZE(value), MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)

    endif
  end subroutine sendAttributePointerInteger1DFlexi

  subroutine sendAttributePointerReal1DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    real, pointer :: value(:)
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name
    if (associated(value)) then

       call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(SIZE(value), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(value, SIZE(value), MPI_REAL, iThread, tag, localWorldCommunicator, ierr)

    endif
  end subroutine sendAttributePointerReal1DFlexi

  subroutine sendAttributePointerDouble1DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    real(double), pointer :: value(:)
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name
    if (associated(value)) then

       call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(SIZE(value), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(value, SIZE(value), MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)

    endif
  end subroutine sendAttributePointerDouble1DFlexi

  subroutine sendAttributePointerDouble2DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    real(double), pointer :: value(:,:)
    real(double), allocatable :: temp(:)
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name
    if (associated(value)) then

       call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(SIZE(value,1), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(SIZE(value,2), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
       allocate(temp(SIZE(value,1)*SIZE(value,2)))
       temp = RESHAPE(value,SHAPE(temp))
       call MPI_SEND(temp, SIZE(temp), MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       deallocate(temp)

    endif
  end subroutine sendAttributePointerDouble2DFlexi

  subroutine sendAttributePointerReal2DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    real, pointer :: value(:,:)
    real, allocatable :: temp(:)
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name
    if (associated(value)) then

       call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(SIZE(value,1), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(SIZE(value,2), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
       allocate(temp(SIZE(value,1)*SIZE(value,2)))
       temp = RESHAPE(value,SHAPE(temp))
       call MPI_SEND(temp, SIZE(temp), MPI_REAL, iThread, tag, localWorldCommunicator, ierr)
       deallocate(temp)

    endif
  end subroutine sendAttributePointerREAL2DFlexi

  subroutine sendAttributePointerDouble3DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    real(double), pointer :: value(:,:,:)
    real(double), allocatable :: temp(:)
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name
    if (associated(value)) then

       call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(SIZE(value,1), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(SIZE(value,2), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(SIZE(value,3), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
       allocate(temp(SIZE(value,1)*SIZE(value,2)*SIZE(value,3)))
       temp = RESHAPE(value,SHAPE(temp))
       call MPI_SEND(temp, SIZE(temp), MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       deallocate(temp)

    endif
  end subroutine sendAttributePointerDouble3DFlexi

  subroutine sendAttributeStaticLogical1DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    logical :: value(:)
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name
    call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(SIZE(value), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(value, SIZE(value), MPI_LOGICAL, iThread, tag, localWorldCommunicator, ierr)
       
  end subroutine sendAttributeStaticLogical1DFlexi

  subroutine sendAttributeStaticLogicalSingleFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    logical :: value
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name
    call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(value, 1, MPI_LOGICAL, iThread, tag, localWorldCommunicator, ierr)
       
  end subroutine sendAttributeStaticLogicalSingleFlexi

  subroutine sendAttributeStaticIntegerSingleFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    integer :: value
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name
    call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(value, 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
       
  end subroutine sendAttributeStaticIntegerSingleFlexi

  subroutine sendAttributeStaticRealSingleFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    real :: value
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name
    call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(value, 1, MPI_REAL, iThread, tag, localWorldCommunicator, ierr)
       
  end subroutine sendAttributeStaticRealSingleFlexi

  subroutine sendAttributeStaticDoubleSingleFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    real(double) :: value
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name
    call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(value, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       
  end subroutine sendAttributeStaticDoubleSingleFlexi

  subroutine sendAttributeStaticVectorSingleFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    type(VECTOR) :: value
    real(double) :: loc(3)
    integer, parameter :: tag = 50
    integer :: ierr

    loc(1) = value%x
    loc(2) = value%y
    loc(3) = value%z
    attributeName = name
    call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       
  end subroutine sendAttributeStaticVECTORSingleFlexi

  subroutine sendAttributeStaticCharacterSingleFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    character(len=*) :: value
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name
    call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(len(value), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(value, len(value), MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
       
  end subroutine sendAttributeStaticCharacterSingleFlexi

  subroutine sendAttributeStaticVector1DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    TYPE(VECTOR) :: value(:)
    real(double), allocatable :: temp(:)
    integer, parameter :: tag = 50
    integer :: ierr
    integer :: n, n3
    integer :: i,j

    attributeName = name

    call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
    n = SIZE(value)
    n3 = 3 * n
    allocate(temp(1:n3))
    do i = 1, n
       j = (i-1)*3
       temp(j+1) = value(i)%x
       temp(j+2) = value(i)%y
       temp(j+3) = value(i)%z
    enddo
    call MPI_SEND(n, 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(temp, n3, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
    deallocate(temp)
  end subroutine sendAttributeStaticVector1DFlexi
  

  subroutine sendAttributeStaticInteger1DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    integer :: value(:)
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name

    call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(SIZE(value), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(value, SIZE(value), MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)

  end subroutine sendAttributeStaticInteger1DFlexi

  subroutine sendAttributeStaticReal1DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    real :: value(:)
    integer, parameter :: tag = 50
    integer :: ierr

    attributeName = name

    call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(SIZE(value), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(value, SIZE(value), MPI_REAL, iThread, tag, localWorldCommunicator, ierr)

  end subroutine sendAttributeStaticReal1DFlexi

  subroutine sendAttributeStaticDouble1DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    real(double) :: value(:)
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name

    call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(SIZE(value), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(value, SIZE(value), MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)

  end subroutine sendAttributeStaticDouble1DFlexi

  subroutine sendAttributeStaticDouble2DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    real(double) :: value(:,:)
    real(double), allocatable :: temp(:)
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name

    call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(SIZE(value,1), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(SIZE(value,2), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
    allocate(temp(SIZE(value,1)*SIZE(value,2)))
    temp = RESHAPE(value,SHAPE(temp))
    call MPI_SEND(temp, SIZE(temp), MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
    deallocate(temp)

     end subroutine sendAttributeStaticDouble2DFlexi

  subroutine sendAttributeStaticReal2DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    real :: value(:,:)
    real, allocatable :: temp(:)
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name

    call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(SIZE(value,1), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(SIZE(value,2), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
    allocate(temp(SIZE(value,1)*SIZE(value,2)))
    temp = RESHAPE(value,SHAPE(temp))
    call MPI_SEND(temp, SIZE(temp), MPI_REAL, iThread, tag, localWorldCommunicator, ierr)
    deallocate(temp)

  end subroutine sendAttributeStaticReal2DFlexi

  subroutine sendAttributeStaticDouble3DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    real(double) :: value(:,:,:)
    real(double), allocatable :: temp(:)
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name

    call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(SIZE(value,1), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(SIZE(value,2), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
    call MPI_SEND(SIZE(value,3), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
    allocate(temp(SIZE(value,1)*SIZE(value,2)*SIZE(value,3)))
    temp = RESHAPE(value,SHAPE(temp))
    call MPI_SEND(temp, SIZE(temp), MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
    deallocate(temp)
    
  end subroutine sendAttributeStaticDouble3DFlexi

  subroutine sendAttributePointerLogical1DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    logical, pointer :: value(:)
    integer, parameter :: tag = 50
    integer :: ierr


    attributeName = name
    if (associated(value)) then

       call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(SIZE(value), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(value, SIZE(value), MPI_LOGICAL, iThread, tag, localWorldCommunicator, ierr)

    endif
  end subroutine sendAttributePointerLogical1DFlexi

  subroutine sendAttributePointerVector1DFlexi(iThread, name, value)
    use mpi
    integer :: iThread
    character(len=*) :: name
    character(len=20) :: attributeName
    TYPE(VECTOR), pointer :: value(:)
    real(double), allocatable :: temp(:)
    integer, parameter :: tag = 50
    integer :: ierr
    integer :: n, n3
    integer :: i,j

    attributeName = name
    if (associated(value)) then

       call MPI_SEND(attributeName, 20, MPI_CHARACTER, iThread, tag, localWorldCommunicator, ierr)
       n = SIZE(value)
       n3 = 3 * n
       allocate(temp(1:n3))
       do i = 1, n
          j = (i-1)*3
          temp(j+1) = value(i)%x
          temp(j+2) = value(i)%y
          temp(j+3) = value(i)%z
       enddo
       call MPI_SEND(n, 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(temp, n3, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       deallocate(temp)
    endif
  end subroutine sendAttributePointerVector1DFlexi



    subroutine receiveArrayIntegerFlexi(value)
      use mpi
      integer, intent(out) :: value(:)
      integer :: n 
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      call MPI_RECV(n, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      call MPI_RECV(value, n, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)

    end subroutine receiveArrayIntegerFlexi

    subroutine receiveArrayRealFlexi(value)
      use mpi
      real, intent(out) :: value(:)
      integer :: n 
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      call MPI_RECV(n, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      call MPI_RECV(value, n, MPI_REAL, 0, tag, localWorldCommunicator, status, ierr)

    end subroutine receiveArrayRealFlexi

    subroutine receiveArrayLogicalFlexi(value)
      use mpi
      logical, intent(out) :: value(:)
      integer :: n 
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      call MPI_RECV(n, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      call MPI_RECV(value, n, MPI_LOGICAL, 0, tag, localWorldCommunicator, status, ierr)

    end subroutine receiveArrayLogicalFlexi

    subroutine receiveArrayDoubleFlexi(value)
      use mpi
      real(double), intent(out) :: value(:)
      integer :: n 
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      call MPI_RECV(n, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      call MPI_RECV(value, n, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

    end subroutine receiveArrayDoubleFlexi

    subroutine receiveArrayVectorFlexi(value)
      use mpi
      type(VECTOR), intent(out) :: value(:)
      integer :: n , n3, i, j
      real(double), allocatable :: temp(:)
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      call MPI_RECV(n, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      n3 = 3 * n
      allocate(temp(1:n3))
      call MPI_RECV(temp, n3, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
      do i = 1, n
         j = (i-1)*3
         value(i) = VECTOR(temp(j+1), temp(j+2), temp(j+3))
      enddo
      deallocate(temp)

    end subroutine receiveArrayVectorFlexi



    subroutine receiveIntegerPointer1dFlexi(value)
      use mpi
      integer, pointer :: value(:)
      integer :: n 
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      call MPI_RECV(n, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      allocate(value(1:n))
      call MPI_RECV(value, n, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)

    end subroutine receiveIntegerPointer1dFlexi

    subroutine receiveRealPointer1dFlexi(value)
      use mpi
      real, pointer :: value(:)
      integer :: n 
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      call MPI_RECV(n, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      allocate(value(1:n))
      call MPI_RECV(value, n, MPI_REAL, 0, tag, localWorldCommunicator, status, ierr)

    end subroutine receiveRealPointer1dFlexi

    subroutine receiveDoublePointer1dFlexi(value)
      use mpi
      real(double), pointer :: value(:)
      integer :: n 
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      call MPI_RECV(n, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      allocate(value(1:n))
      call MPI_RECV(value, n, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)

    end subroutine receiveDoublePointer1dFlexi

    subroutine receiveDoublePointer2dFlexi(value)
      use mpi
      real(double), pointer :: value(:,:)
      integer :: n, m
      real(double), allocatable :: temp(:)
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      call MPI_RECV(n, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      call MPI_RECV(m, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      allocate(value(1:n, 1:m))
      allocate(temp(1:(n*m)))
      call MPI_RECV(temp, n*m, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
      value(1:n,1:m) = RESHAPE(temp, SHAPE(value))
      deallocate(temp)
    end subroutine receiveDoublePointer2dFlexi

    subroutine receiveDoublePointer3dFlexi(value)
      use mpi
      real(double), pointer :: value(:,:,:)
      integer :: n, m, l
      real(double), allocatable :: temp(:)
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      call MPI_RECV(n, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      call MPI_RECV(m, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      call MPI_RECV(l, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      allocate(value(1:n, 1:m, 1:l))
      allocate(temp(1:(n*m*l)))
      call MPI_RECV(temp, n*m*l, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
      value(1:n,1:m,1:l) = RESHAPE(temp, SHAPE(value))
      deallocate(temp)
    end subroutine receiveDoublePointer3dFlexi

    subroutine receiveRealPointer2dFlexi(value)
      use mpi
      real, pointer :: value(:,:)
      integer :: n, m
      real, allocatable :: temp(:)
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      call MPI_RECV(n, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      call MPI_RECV(m, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      allocate(value(1:n, 1:m))
      allocate(temp(1:(n*m)))
      call MPI_RECV(temp, n*m, MPI_REAL, 0, tag, localWorldCommunicator, status, ierr)
      value(1:n,1:m) = RESHAPE(temp, SHAPE(value))
      deallocate(temp)
    end subroutine receiveRealPointer2dFlexi


    subroutine receiveLogicalPointer1dFlexi(value)
      use mpi
      logical, pointer :: value(:)
      integer :: n 
      integer :: status(MPI_STATUS_SIZE)
      integer, parameter :: tag = 50
      integer :: ierr

      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      call MPI_RECV(n, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      allocate(value(1:n))
      call MPI_RECV(value, n, MPI_LOGICAL, 0, tag, localWorldCommunicator, status, ierr)

    end subroutine receiveLogicalPointer1dFlexi



    subroutine receiveVectorPointer1dFlexi(value)
      use mpi
      type(VECTOR), pointer :: value(:)
      integer :: n, n3, i, j
      integer :: status(MPI_STATUS_SIZE)
      real(double), allocatable :: temp(:)
      integer, parameter :: tag = 50
      integer :: ierr

      if (associated(value)) then
         deallocate(value)
         nullify(value)
      endif
      call MPI_RECV(n, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
      allocate(value(1:n))
      n3 = 3 * n
      allocate(temp(1:n3))
      call MPI_RECV(temp, n3, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
      do i = 1, n
         j = (i-1)*3
         value(i) = VECTOR(temp(j+1), temp(j+2), temp(j+3))
      enddo
      deallocate(temp)
    end subroutine receiveVectorPointer1dFlexi

#endif

    subroutine testDataType(thisType, fileFormatted)
      character(len=*) :: thisType
      logical :: fileFormatted
      character(len=10) :: dataType
      character(len=80) :: message
      if (fileFormatted) then
         read(20, *) dataType
      else
         if (uncompressedDumpFiles) then
            read(20) dataType
         else
            call readCompressedFile(20, dataType)
         endif
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
         if (uncompressedDumpFiles) then
            read(lUnit) dataType
         else
            call readCompressedFile(lUnit, dataType)
         endif
      endif
      dataType = ADJUSTL(dataType)
      select case (dataType)
         case("rsingle")
            if (fileFormatted) then
               read(lUnit,*) rDummy
            else
               if (uncompressedDumpFiles) then
                  read(lUnit) rDummy
               else
                  call readCompressedFile(lUnit, rDummy)
               endif
            endif
         case("dsingle")
            if (fileFormatted) then
               read(lUnit,*) dDummy
            else
               if (uncompressedDumpFiles) then
                  read(lUnit) dDummy
               else
                  call readCompressedFile(lUnit, dDummy)
               endif
            endif
         case("isingle")
            if (fileFormatted) then
               read(lUnit,*) iDummy
            else
               if (uncompressedDumpFiles) then
                  read(lUnit) iDummy
               else
                  call readCompressedFile(lUnit, iDummy)
               endif
            endif
         case("vsingle")
            if (fileFormatted) then
               read(lUnit,*) vDummy
            else
               if (uncompressedDumpFiles) then
                  read(lUnit) vDummy
               else
                  call readCompressedFile(lUnit, vDummy)
               endif
            endif
         case("lsingle")
            if (fileFormatted) then
               read(lUnit,*) lDummy
            else
               if (uncompressedDumpFiles) then
                  read(lUnit) lDummy
               else
                  call readCompressedFile(lUnit, lDummy)
               endif
            endif
         case("csingle")
            if (fileFormatted) then
               read(lUnit,*) n
               read(lUnit,*) cDummy(1:n)
            else
               if (uncompressedDumpFiles) then
                  read(lUnit) n
                  read(lUnit) cDummy(1:n)
               else
                  call readCompressedFile(lUnit, n)
                  call readCompressedFile(lUnit, cDummy(1:n))
               endif
            endif
         case("r1darray")
            if (fileFormatted) then
               read(lUnit,*) n
               allocate(rArray(1:n))
               read(lUnit,*) rArray(1:n)
            else
               if (uncompressedDumpFiles) then
                  read(lUnit) n
                  allocate(rArray(1:n))
                  read(lUnit) rArray(1:n)
               else
                  call readCompressedFile(lUnit, n)
                  allocate(rArray(1:n))
                  call readCompressedFile(lUnit, rArray(1:n))
               endif
            endif
            deallocate(rArray)
         case("v1darray")
            if (fileFormatted) then
               read(lUnit,*) n
               allocate(vArray(1:n))
               read(lUnit,*) vArray(1:n)
            else
               if (uncompressedDumpFiles) then
                  read(lUnit) n
                  allocate(vArray(1:n))
                  read(lUnit) vArray(1:n)
               else
                  call readCompressedFile(lUnit, n)
                  allocate(vArray(1:n))
                  call readCompressedFile(lUnit, vArray(1:n))
               endif
            endif
            deallocate(vArray)
         case("d1darray")
            if (fileFormatted) then
               read(lUnit,*) n
               allocate(dArray(1:n))
               read(lUnit,*) dArray(1:n)
            else
               if (uncompressedDumpFiles) then
                  read(lUnit) n
                  allocate(dArray(1:n))
                  read(lUnit) dArray(1:n)
               else
                  call readCompressedFile(lUnit, n)
                  allocate(dArray(1:n))
                  call readCompressedFile(lUnit, dArray(1:n))
               endif
            endif
            deallocate(dArray)
         case("i1darray")
            if (fileFormatted) then
               read(lUnit,*) n
               allocate(iArray(1:n))
               read(lUnit,*) iArray(1:n)
            else
               if (uncompressedDumpFiles) then
                  read(lUnit) n
                  allocate(iArray(1:n))
                  read(lUnit) iArray(1:n)
               else
                  call readCompressedFile(lUnit, n)
                  allocate(iArray(1:n))
                  call readCompressedFile(lUnit, iArray(1:n))
               endif
            endif
            deallocate(iArray)
         case("l1darray")
            if (fileFormatted) then
               read(lUnit,*) n
               allocate(lArray(1:n))
               read(lUnit,*) lArray(1:n)
            else
               if (uncompressedDumpFiles) then
                  read(lUnit) n
                  allocate(lArray(1:n))
                  read(lUnit) lArray(1:n)
               else
                  call readCompressedFile(lUnit, n)
                  allocate(lArray(1:n))
                  call readCompressedFile(lUnit, lArray(1:n))
               endif
            endif
            deallocate(lArray)
         case("d2darray")
            if (fileFormatted) then
               read(lUnit,*) n,m
               allocate(d2Array(1:n,1:m))
               read(lUnit,*) d2Array(1:n,1:m)
            else
               if (uncompressedDumpFiles) then
                  read(lUnit) n,m
                  allocate(d2Array(1:n,1:m))
                  read(lUnit) d2Array(1:n,1:m)
               else
                  call readCompressedFile(lUnit, n)
                  call readCompressedFile(lUnit, m)
                  allocate(d2Array(1:n,1:m))
                  call readCompressedFile(lUnit, d2Array(1:n,1:m))
               endif
            endif
            deallocate(d2Array)
         case("d3darray")
            if (fileFormatted) then
               read(lUnit,*) n,m,j
               allocate(d3Array(1:n,1:m,1:j))
               read(lUnit,*) d3Array(1:n,1:m,1:j)
            else
               if (uncompressedDumpFiles) then
                  read(lUnit) n,m,j
                  allocate(d3Array(1:n,1:m,1:j))
                  read(lUnit) d3Array(1:n,1:m,1:j)
               else
                  call readCompressedFile(lUnit, n)
                  call readCompressedFile(lUnit, m)
                  call readCompressedFile(lUnit, j)
                  allocate(d3Array(1:n,1:m,1:j))
                  call readCompressedFile(lUnit, d3Array(1:n,1:m,1:j))
               endif
            endif
            deallocate(d3Array)
         case DEFAULT
            write(*,*) myrankGlobal, " Data type not recognised: ",trim(dataType)
            stop
         end select
       end subroutine readDummyData
            

       subroutine skipOctalsToDepth(fileFormatted, nDepth)
         logical :: fileFormatted
         integer :: nDepth, thisNdepth
         character(len=20) :: tag
         logical :: runToEndofOctal
         integer :: iChildren


         runToEndofOctal = .false.
         iChildren = 0
         do while (.true.)

            if (fileFormatted) then
               read(20, *, end = 666) tag
            else
               if (uncompressedDumpFiles) then
                  read(20, end = 666) tag
               else
                  call readCompressedFile(20, tag)
               endif
            endif
            tag = ADJUSTL(tag)
!            write(*,*) myrankglobal, " tag ",tag
            if (tag == "OCTALBEGINS") cycle
            if ((tag == "OCTALENDS").and.(runToEndofOctal)) exit
            if (tag == "OCTALENDS") cycle

            select case (tag)
            case ("nDepth") 
               call readSingleFlexi(20, thisNDepth, fileFormatted)
               if (thisNDepth == nDepth) then
                  runToEndOfOctal = .true.
               endif
            case DEFAULT
               call readDummyData(20, fileFormatted)
            end select



         enddo
666      continue

       end subroutine skipOctalsToDepth

#ifdef MPI

       subroutine readGridSplitOverMPI(grid, gridFilename, fileFormatted)
         type(GRIDTYPE) :: grid
         character(len=*) :: gridFilename
         
         character(len=80) :: message
         logical :: fileFormatted
         integer :: iThread, nOctal, nOctals, nVoxels, ierr


!         do iSet = 0, nHydroSetsGlobal - 1
!            if (myHydroSetGlobal == iSet) then
               
               if (associated(grid%octreeRoot)) then
                  call deleteOctreeBranch(grid%octreeRoot,onlyChildren=.false., adjustParent=.false.)
                  grid%octreeRoot => null()
               endif
               
               do iThread = 1, nHydroThreadsGlobal
                  if (myrankGlobal == iThread) then
                     call openGridFile(gridFilename, fileformatted)
                     call readStaticComponents(grid, fileFormatted)
                     close(20)
                  endif
                  call mpi_barrier(localWorldCommunicator, ierr)
               enddo
               
               if (myrankGlobal /= 0) then
                  allocate(grid%octreeRoot)
                  grid%octreeRoot%nDepth = 1
                  nOctal = 0
                  call getBranchOverMPI(grid%octreeRoot, null())
                  write(message, '(a,i3.3,a,i3.3)') "AMR grid read for thread: ",myrankGlobal, " set ",myHydroSetGlobal
                  write(*,*) trim(message)
               else
                  call openGridFile(gridFilename, fileformatted)
                  call readStaticComponents(grid, fileFormatted)
!                  write(*,*) "static components for zero thread read in readgridsplitovermpi"
                  allocate(grid%octreeRoot)
                  grid%octreeRoot%nDepth = 1
                  nOctal = 0
                  call readZerothThread(grid%octreeRoot, null(), fileFormatted)
               endif
               call updateMaxDepth(grid)
               call setSmallestSubcell(grid)
               call countVoxels(grid%octreeRoot,nOctals,nVoxels)
               grid%nOctals = nOctals
               !         call checkAMRgrid(grid, .false.)
               close(20)
!            endif
!            call torus_mpi_barrier
!         enddo

       end subroutine readGridSplitOverMPI

       subroutine readZerothThread(thisOctal, parent, fileFormatted)
         logical :: fileFormatted
         type(OCTAL), pointer :: thisOctal, parent, child, depth3, nextOctal
!         character(len=80) :: message
         integer :: i, iChild, ithread, ichildbig


         allocate(depth3)
         thisOctal%parent => parent
         call readOctalViaTags(thisOctal, fileFormatted)
         
         if (nHydroThreadsGlobal == 8) then
            thisOctal%nChildren = 1
            do i = 1, nHydroThreadsGlobal
               thisOctal%hasChild = .false.
               thisOctal%hasChild(i) = .true.
               thisOctal%indexChild(1) = i
               call sendOctalviaMPI(thisOctal,i)
            enddo

            do i = 1, nHydroThreadsGlobal
               call readBranchFromFile(i, fileFormatted)
            enddo

            thisOctal%nChildren = 0 
            thisOctal%hasChild = .false.
         endif

         if (nHydroThreadsGlobal == 4) then
            thisOctal%nChildren = 1
            do i = 1, nHydroThreadsGlobal
               thisOctal%hasChild = .false.
               thisOctal%hasChild(i) = .true.
               thisOctal%indexChild(1) = i
               call sendOctalviaMPI(thisOctal,i)
            enddo

            do i = 1, nHydroThreadsGlobal
               call readBranchFromFile(i, fileFormatted)
            enddo

!            write(*,*) myrankglobal, "setting nchildren to zero"
            thisOctal%nChildren = 0 
            thisOctal%hasChild = .false.
         endif


         if (nHydroThreadsGlobal == 16) then

!            allocate(child)
!            do while (.true.)
!               call readOctalViaTags(child, fileFormatted)
!               write(*,*) "from tags depth: ",child%ndepth, child%mpiThread, child%nchildren
!            enddo

            do i = 1, nHydroThreadsGlobal
               call sendOctalviaMPI(thisOctal,i)
            enddo


            allocate(thisOctal%child(1:thisOctal%nChildren))
            do iChild = 1, thisOctal%nChildren
               child => thisOctal%child(iChild)
               call readOctalViaTags(child, fileFormatted)
               do i = 1, nHydroThreadsGlobal
                  if (.not.octalonThread(thisOctal, iChild, i)) then 
                     child%hasChild = .false.
                     child%nChildren = 0
                     call sendOctalviaMPI(child,i)
                  endif
               enddo

               do i = 1, 4
                  ithread = (ichild-1)*4 + i
                  child%hasChild = .false.
                  child%nChildren = 1
                  child%indexChild(1) = i
                  child%hasChild(i) = .true.
                  call sendOctalviaMPI(child,ithread)
                  call readBranchFromFile(ithread, fileFormatted)
               enddo
               child%hasChild = .false.
               child%nChildren = 0 
               child%parent => thisOctal
            enddo

         endif

         if (nHydroThreadsGlobal == 64) then

!            allocate(child)
!            do while (.true.)
!               call readOctalViaTags(child, fileFormatted)
!               write(*,*) "from tags depth: ",child%ndepth, child%mpiThread
!            enddo

            do i = 1, nHydroThreadsGlobal
               call sendOctalviaMPI(thisOctal,i)
            enddo


            allocate(thisOctal%child(1:thisOctal%nChildren))
            do iChild = 1, thisOctal%nChildren
               child => thisOctal%child(iChild)
               call readOctalViaTags(child, fileFormatted)
               do i = 1, nHydroThreadsGlobal
                  if (.not.octalonThread(thisOctal, iChild, i)) then 
                     child%hasChild = .false.
                     child%nChildren = 0
                     call sendOctalviaMPI(child,i)
                  endif
               enddo

               do i = 1, 8
                  ithread = (ichild-1)*8 + i
                  child%hasChild = .false.
                  child%nChildren = 1
                  child%indexChild(1) = i
                  child%hasChild(i) = .true.
                  call sendOctalviaMPI(child,ithread)
                  call readBranchFromFile(ithread, fileFormatted)
               enddo
               child%hasChild = .false.
               child%nChildren = 0 
               child%parent => thisOctal
            enddo

         endif

         if (nHydroThreadsGlobal == 512) then

!            allocate(child)
!            do while (.true.)
!               call readOctalViaTags(child, fileFormatted)
!               write(*,*) "from tags depth: ",child%ndepth, child%mpiThread
!            enddo

! level 1
            do i = 1, nHydroThreadsGlobal
               call sendOctalviaMPI(thisOctal,i)
            enddo
! level 2
            allocate(thisOctal%child(1:thisOctal%nChildren))

            do iChildBig = 1, thisOctal%nChildren
               child => thisOctal%child(ichildBig)
               call readOctalViaTags(child, fileFormatted)
               do i = 1, nHydroThreadsGlobal
                  call sendOctalviaMPI(child,i)
               enddo
               nextOctal => child
               allocate(nextOctal%child(1:nextOctal%nChildren))
               
               do iChild = 1, nextOctal%nChildren
                  child => nextOctal%child(iChild)
                  call readOctalViaTags(child, fileFormatted)
                  do i = 1, nHydroThreadsGlobal
                     if (.not.octalonThread(nextOctal, iChild, i)) then 
                        child%hasChild = .false.
                        child%nChildren = 0
                        call sendOctalviaMPI(child,i)
                     endif
                  enddo
                  
                  do i = 1, 8
                     ithread = (iChildBig-1)*64+(ichild-1)*8 + i
                     child%hasChild = .false.
                     child%nChildren = 1
                     child%indexChild(1) = i
                     child%hasChild(i) = .true.
                     call sendOctalviaMPI(child,ithread)
                     call readBranchFromFile(ithread, fileFormatted)
                  enddo
                  child%hasChild = .false.
                  child%nChildren = 0 
                  child%parent => nextOctal
               enddo
            enddo

         endif

       end subroutine readZerothThread


       subroutine readGridMPI(grid, gridfilename, fileFormatted)
         use inputs_mod, only : splitOverMPI, modelwashydro
         type(GRIDTYPE) :: grid
         character(len=*) :: gridFilename
         logical :: fileFormatted
         integer :: iThread
         integer :: nOctals, nVoxels
         character(len=80) :: message
         logical, save :: firstTime = .true.
         type(romanova) :: romdata

         if (associated(grid%octreeRoot)) then
            call deleteOctreeBranch(grid%octreeRoot,onlyChildren=.false., adjustParent=.false.)
            grid%octreeRoot => null()
         endif
         
         if (splitOverMPI) then
            do iThread = 1, nThreadsGlobal - 1
               if (myrankGlobal == iThread) then
                  call openGridFile(gridFilename, fileformatted)
                  call readStaticComponents(grid, fileFormatted)
                  close(20)
               endif
               call torus_mpi_barrier
            enddo
            
            if (myrankGlobal == 0) then
               call readAmrGridSingle(gridfilename, fileFormatted, grid)
               call sendGridToAllThreads(grid%octreeRoot)
            else
               allocate(grid%octreeRoot)
               grid%octreeRoot%nDepth = 1
               nOctals = 0
               call getBranchOverMPI(grid%octreeRoot, null())
            endif
            write(message,'(a,i3,a)') "Thread ",ithread, " read"
         else
!            do iThread = 0, nThreadsGlobal - 1
!               if (myrankGlobal == iThread) call readAmrGridSingle(gridfilename, fileFormatted, grid)
!               write(message,'(a,i3,a)') "Thread ",ithread, " read"
!               call writeInfo(message,TRIVIAL)
!               call torus_mpi_barrier
!            enddo
            
            call readAmrGridSingle(gridfilename, fileFormatted, grid)
            call writeInfo(message,TRIVIAL)

         endif


         if(modelwashydro) then
            if(firstTime) then
               write(message,'(a,i3,a)') "Adding cell corner values"               
               call writeInfo(message,TRIVIAL)
               firsttime = .false.
            end if
            call finishgrid(grid%octreeRoot, grid, romData=romData)                
         end if


          call torus_mpi_barrier
          call updateMaxDepth(grid)
          call setSmallestSubcell(grid)
          call countVoxels(grid%octreeRoot,nOctals,nVoxels)
          grid%nOctals = nOctals
        end subroutine readGridMPI

        recursive subroutine sendGridToAllThreads(thisOctal)
          type(OCTAL), pointer :: thisOctal, child
          integer :: iThread, n, i
          do iThread = 1, nThreadsGlobal-1
             call sendOctalViaMPI(thisOctal, iThread)
          enddo
          n = thisOctal%nChildren
          if (n > 0) then
             do i = 1, n
                child => thisOctal%child(i)
                call sendGridToAllThreads(child)
             enddo
          endif
        end subroutine sendGridToAllThreads

       recursive subroutine readBranchFromFile(iThread, fileFormatted)
         integer :: ithread
         logical :: fileFormatted
         type(OCTAL), pointer :: thisOctal => null()
         integer :: i, n

         allocate(thisOctal)
         call readOctalViaTags(thisOctal, fileFormatted)
         call sendOctalViaMPI(thisOctal, iThread)
         n = thisOctal%nChildren
         call deallocateOctalDynamicAttributes(thisOctal)
         deallocate(thisOctal)
         if (n > 0) then
            do i = 1, n
               call readBranchFromFile(iThread, fileFormatted)
            enddo
         endif
       end subroutine readBranchFromFile


       recursive subroutine getBranchOverMPI(thisOctal, parent)
         type(OCTAL), pointer :: thisOctal
         type(OCTAL), pointer :: parent
         type(OCTAL), pointer :: child
         integer :: i

         thisOctal%parent => parent
!         write(*,*) myrankGlobal, " receiving octal from 0"
         call receiveOctalViaMPI(thisOctal)
!         write(*,'(i4,a,i4,i4,8i4)') myrankGlobal, " received successfully. depth ",thisOctal%nDepth, thisOctal%nChildren, thisOctal%mpiThread
         if (thisOctal%nChildren > 0) then
            allocate(thisOctal%child(1:thisOctal%nChildren)) 
            do i = 1, thisOctal%nChildren
               child => thisOctal%child(i)
!               write(*,*) myrankGlobal," getting branch for child ",i, " of ",thisOctal%nchildren, " at depth ",thisOctal%nDepth
               call getBranchOverMpi(child, thisOctal)
            enddo
         endif
       end subroutine getBranchOverMPI

#endif

       subroutine openGridFile(gridFilename, fileFormatted)
         character(len=*) :: gridFilename
         logical :: fileFormatted
         integer, dimension(8) :: timeValues    ! system date and time
         integer               :: error         ! status code
         character(len=80) :: message
         if (fileFormatted) then
            open(unit=20, iostat=error, file=gridFilename, form="formatted", status="old")
            if (error /=0) then
               print *, 'Panic: file open error in readAMRgrid, file:',trim(gridFilename) ; stop               
            end if
            ! read the file's time stamp
            read(unit=20,fmt=*,iostat=error) timeValues 
            if (error /=0) then
               print *, 'Panic: read error in readAMRgrid (formatted timeValues)' ; stop
            end if
         else
            if (uncompressedDumpFiles) then
               open(unit=20, iostat=error, file=gridFilename, form="unformatted", status="old")
            else
#ifdef USEZLIB
               error = 0
               call openCompressedFile(20, gridFilename)
#else
               call writeFatal("zlib is needed to read compressed files")
#endif
            endif
            if (error /=0) then
               print *, 'Panic: file open error in readAMRgrid, file:',trim(gridFilename) ; stop
            end if
            ! read the file's time stamp
            read(unit=20,iostat=error) timeValues
            if (error /=0) then
               print *, 'Panic: read error in readAMRgrid (unformatted timeValues)' ; stop
            end if
         end if

         write(message,'(a,a)') "Reading flexible AMR file from: ",trim(gridfilename)
         call writeInfo(message,TRIVIAL)

         write(message,'(a,i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') ' - data file written at: ', &
              timeValues(1),'/',timeValues(2),'/',&
              timeValues(3),'  ',timeValues(5),':',timeValues(6)
         call writeInfo(message,TRIVIAL)


       end subroutine openGridFile


       subroutine readStaticComponents(grid, fileformatted)
         type(GRIDTYPE) :: grid
         logical :: fileFormatted
         character(len=20) :: tag
         character(len=80) :: message

         do while(.true.)
            if (fileFormatted) then 
               read(20, '(a)') tag
            else
               if (uncompressedDumpFiles) then
                  read(20) tag
               else
                  call readCompressedFile(20, tag)
               endif
            endif
            tag = ADJUSTL(tag)
!            write(*,*) myrankglobal," tag: ",trim(tag)

            if (tag == "GRIDBEGINS") cycle
            if (tag == "GRIDENDS") exit

            select case (trim(tag))
            case("version")
               call readSingleFlexi(20, grid%version, fileFormatted)
               write(message,'(a,a,a,a)') "Dump file written by TORUS version  ", trim(torusVersion)
               call writeInfo(message, IMPORTANT)
               if (grid%version /= torusVersion) then
                  write(message,'(a,a,a,a)') "This dump file written with ", trim(grid%version), &
                       " and read with ", trim(torusVersion)
                  call writeWarning(message)
                  grid%version = torusVersion
               endif
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
   end subroutine readStaticComponents


   subroutine readOctalViaTags(thisOctal, fileformatted)
     type(OCTAL)  :: thisOctal
     logical :: fileFormatted
     character(len=20) :: tag
     character(len=80) :: message
!     integer :: iThread

      do while (.true.)

         if (fileFormatted) then
            read(20, *) tag
         else
            if (uncompressedDumpFiles) then
               read(20) tag
            else
               call readCompressedfile(20, tag)
            endif
         endif
         tag = ADJUSTL(tag)
!         write(*,*) myrankglobal, " tag ",tag
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
         case("oneD")
            call readSingleFlexi(20, thisOctal%oneD, fileFormatted)
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
         case("phimin")
            call readSingleFlexi(20, thisOctal%phimin, fileFormatted)
         case("phimax")
            call readSingleFlexi(20, thisOctal%phimax, fileFormatted)
         case("r")
            call readSingleFlexi(20, thisOctal%r, fileFormatted)

         case("xMax")
            call readSingleFlexi(20, thisOctal%xMax, fileFormatted)
         case("yMax")
            call readSingleFlexi(20, thisOctal%yMax, fileFormatted)
         case("zMax")
            call readSingleFlexi(20, thisOctal%zMax, fileFormatted)
         case("xMin")
            call readSingleFlexi(20, thisOctal%xMin, fileFormatted)
         case("yMin")
            call readSingleFlexi(20, thisOctal%yMin, fileFormatted)
         case("zMin")
            call readSingleFlexi(20, thisOctal%zMin, fileFormatted)

         case("parentSubcell")
            call readSingleFlexi(20, thisOctal%parentSubcell, fileFormatted)
         case("inStar")
            call readArrayFlexi(20, thisOctal%inStar, fileFormatted)
         case("inFlow")
            call readArrayFlexi(20, thisOctal%inFlow, fileFormatted)
         case("velocity")
            call readArrayFlexi(20, thisOctal%velocity, fileFormatted)
         case("cornervelocity")
            call readPointerFlexi(20, thisOctal%cornervelocity, fileFormatted)
         case("cornerrho")
            call readPointerFlexi(20, thisOctal%cornerrho, fileFormatted)
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
         case("oldFrac")
            call readPointerFlexi(20, thisOctal%oldFrac, fileFormatted)
         case("scatteredIntensity")
            call readPointerFlexi(20, thisOctal%scatteredIntensity, fileFormatted)
         case("meanIntensity")
            call readPointerFlexi(20, thisOctal%meanIntensity, fileFormatted)
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
         case("distanceGrid")
            call readPointerFlexi(20, thisOctal%distanceGrid, fileFormatted)

         case("nCrossings")
            call readPointerFlexi(20, thisOctal%nCrossings, fileFormatted)
         case("hHeating")
            call readPointerFlexi(20, thisOctal%hHeating, fileFormatted)
         case("heHeating")
            call readPointerFlexi(20, thisOctal%heHeating, fileFormatted)

         case("undersampled")
            call readPointerFlexi(20, thisOctal%undersampled, fileFormatted)
         case("nDiffusion")
            call readPointerFlexi(20, thisOctal%nDiffusion, fileFormatted)
         case("diffusionApprox")
            call readPointerFlexi(20, thisOctal%diffusionApprox, fileFormatted)


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

         case("q_i")
            call readPointerFlexi(20, thisOctal%q_i, fileFormatted)
         case("q_i_plus_1")
            call readPointerFlexi(20, thisOctal%q_i_plus_1, fileFormatted)
         case("q_i_minus_1")
            call readPointerFlexi(20, thisOctal%q_i_minus_1, fileFormatted)
         case("q_i_minus_2")
            call readPointerFlexi(20, thisOctal%q_i_minus_2, fileFormatted)

         case("x_i")
            call readPointerFlexi(20, thisOctal%x_i, fileFormatted)
         case("x_i_plus_1")
            call readPointerFlexi(20, thisOctal%x_i_plus_1, fileFormatted)
         case("x_i_minus_1")
            call readPointerFlexi(20, thisOctal%x_i_minus_1, fileFormatted)
         case("x_i_minus_2")
            call readPointerFlexi(20, thisOctal%x_i_minus_2, fileFormatted)

         case("u_interface")
            call readPointerFlexi(20, thisOctal%u_interface, fileFormatted)
         case("u_i")
            call readPointerFlexi(20, thisOctal%u_i, fileFormatted)
         case("u_i_plus_1")
            call readPointerFlexi(20, thisOctal%u_i_plus_1, fileFormatted)
         case("u_i_minus_1")
            call readPointerFlexi(20, thisOctal%u_i_minus_1, fileFormatted)

         case("flux_i")
            call readPointerFlexi(20, thisOctal%flux_i, fileFormatted)
         case("flux_i_plus_1")
            call readPointerFlexi(20, thisOctal%flux_i_plus_1, fileFormatted)
         case("flux_i_minus_1")
            call readPointerFlexi(20, thisOctal%flux_i_minus_1, fileFormatted)


         case("phiLimit")
            call readPointerFlexi(20, thisOctal%phiLimit, fileFormatted)

         case("ghostCell")
            call readPointerFlexi(20, thisOctal%ghostCell, fileFormatted)
            !print *, "read ghostcell ", thisOctal%ghostCell 

         case("corner")
            call readPointerFlexi(20, thisOctal%corner, fileFormatted)


         case("feederCell")
            call readPointerFlexi(20, thisOctal%feederCell, fileFormatted)
         case("edgeCell")
            call readPointerFlexi(20, thisOctal%edgeCell, fileFormatted)
         case("refinedLastTime")
            call readPointerFlexi(20, thisOctal%refinedLastTime, fileFormatted)

         case("pressure_i")
            call readPointerFlexi(20, thisOctal%pressure_i, fileFormatted)
         case("pressure_i_plus_1")
            call readPointerFlexi(20, thisOctal%pressure_i_plus_1, fileFormatted)
         case("pressure_i_minus_1")
            call readPointerFlexi(20, thisOctal%pressure_i_minus_1, fileFormatted)

         case("rhou")
            call readPointerFlexi(20, thisOctal%rhou, fileFormatted)
         case("rhov")
            call readPointerFlexi(20, thisOctal%rhov, fileFormatted)
         case("rhow")
            call readPointerFlexi(20, thisOctal%rhow, fileFormatted)

         case("rhorv_i_plus_1")
            call readPointerFlexi(20, thisOctal%rhorv_i_plus_1, fileFormatted)

         case("rhorv_i_minus_1")
            call readPointerFlexi(20, thisOctal%rhorv_i_minus_1, fileFormatted)


         case("rhoe")
            call readPointerFlexi(20, thisOctal%rhoe, fileFormatted)
         case("rhoeLast")
            call readPointerFlexi(20, thisOctal%rhoeLastTime, fileFormatted)
         case("energy")
            call readPointerFlexi(20, thisOctal%energy, fileFormatted)

         case("qViscosity")
            call readPointerFlexi(20, thisOctal%qViscosity, fileFormatted)


         case("phi_i")
            call readPointerFlexi(20, thisOctal%phi_i, fileFormatted)
         case("phi_stars")
            call readPointerFlexi(20, thisOctal%phi_stars, fileFormatted)
         case("phi_gas")
            call readPointerFlexi(20, thisOctal%phi_gas, fileFormatted)
         case("phi_i_plus_1")
            call readPointerFlexi(20, thisOctal%phi_i_plus_1, fileFormatted)
         case("phi_i_minus_1")
            call readPointerFlexi(20, thisOctal%phi_i_minus_1, fileFormatted)

         case("rho_i_plus_1")
            call readPointerFlexi(20, thisOctal%rho_i_plus_1, fileFormatted)
         case("rho_i_minus_1")
            call readPointerFlexi(20, thisOctal%rho_i_minus_1, fileFormatted)

         case("boundaryCondition")
            call readPointerFlexi(20, thisOctal%boundaryCondition, fileFormatted)

         case("boundaryCell")
            call readPointerFlexi(20, thisOctal%boundaryCell, fileFormatted)

           

!#ifdef MPI 
!            if(myRankGlobal == 0) then
!               if(thisOctal%twoD) then
!                  do iThread = 1, 4
!                     if(thisOctal%boundaryCondition(iThread) == 0) then
!                        thisOctal%boundaryCondition(iThread) = 1
!                     end if
!                  enddo
!               end if
!            end if
!#endif
!            print *, "myRankGlobal ", myRankGlobal
!            print *, " read : thisOctal%boundaryCondition ", thisOctal%boundaryCondition
!            print *, " read : ghostcell ", thisOctal%ghostCell

         case("boundaryPartner")
            call readPointerFlexi(20, thisOctal%boundaryPartner, fileFormatted)

         case("radiationMomentum")
            call readPointerFlexi(20, thisOctal%radiationMomentum, fileFormatted)
         case("kappaTimesFlux")
            call readPointerFlexi(20, thisOctal%kappaTimesFlux, fileFormatted)
         case("UVvector")
            call readPointerFlexi(20, thisOctal%UVvector, fileFormatted)
         case("gravboundaryPartner")
            call readPointerFlexi(20, thisOctal%gravboundaryPartner, fileFormatted)
         case("changed")
            call readPointerFlexi(20, thisOctal%changed, fileFormatted)
         case("rLimit")
            call readPointerFlexi(20, thisOctal%rLimit, fileFormatted)

         case("iEquationOfState")
            call readPointerFlexi(20, thisOctal%iEquationOfState, fileFormatted)


         case("iAnalyticalVelocity")
            call readPointerFlexi(20, thisOctal%iAnalyticalVelocity, fileFormatted)

         case("gamma")
            call readPointerFlexi(20, thisOctal%gamma, fileFormatted)

         case("udens")
            call readPointerFlexi(20, thisOctal%udens, fileFormatted)
         case("adot")
            call readPointerFlexi(20, thisOctal%aDot, fileFormatted)
         case("adotdist")
            call readPointerFlexi(20, thisOctal%distanceGridaDot, fileFormatted)
         case("dfromgas")
            call readPointerFlexi(20, thisOctal%distanceGridPhotonFromGas, fileFormatted)
         case("dfromsource")
            call readPointerFlexi(20, thisOctal%distanceGridPhotonFromSource, fileFormatted)
         case("ufromgas")
            call readPointerFlexi(20, thisOctal%photonEnergyDensityFromGas, fileFormatted)
         case("ufromsource")
            call readPointerFlexi(20, thisOctal%photonEnergyDensityFromSource, fileFormatted)
         case("utotal")
            call readPointerFlexi(20, thisOctal%photonEnergyDensityFromSource, fileFormatted)
         case("oldutotal")
            call readPointerFlexi(20, thisOctal%oldphotonEnergyDensity, fileFormatted)

         case("fixedtemperature")
            call readPointerFlexi(20, thisOctal%fixedtemperature, fileFormatted)


         case DEFAULT
            write(message,*) "Unrecognised tag on read: "//trim(tag)
            call writeWarning(message)
            call readDummyData(20, fileFormatted)
         end select

      end do
    end subroutine readOctalViaTags

#ifdef MPI

   subroutine receiveOctalViaMPI(thisOctal)
     use mpi
     type(OCTAL), pointer :: thisOctal
     integer :: status(MPI_STATUS_SIZE)
     integer, parameter :: mpitag = 50
     character(len=20) :: tag
     character(len=80) :: message
     integer :: ierr

      do while (.true.)

         call MPI_RECV(tag, 20, MPI_CHARACTER, 0, mpitag, localWorldCommunicator, status, ierr)
         tag = ADJUSTL(tag)
         if (tag == "OCTALBEGINS") cycle
         if (tag == "OCTALENDS") exit
         select case (tag)
         case("nDepth")
            call receiveSingleFlexi(thisOctal%nDepth)
         case("nChildren")
            call receiveSingleFlexi(thisOctal%nChildren)
         case("indexChild")
            call receiveArrayFlexi(thisOctal%indexChild)
         case("hasChild")
            call receiveArrayFlexi(thisOctal%hasChild)
         case("centre")
            call receiveSingleFlexi(thisOctal%centre)
         case("rho")
            call receiveArrayFlexi(thisOctal%rho)
         case("temperature")
            call receiveArrayFlexi(thisOctal%temperature)
         case("label")
            call receiveArrayFlexi(thisOctal%label)
         case("subcellSize")
            call receiveSingleFlexi(thisOctal%subcellSize)
         case("threeD")
            call receiveSingleFlexi(thisOctal%threeD)
         case("twoD")
            call receiveSingleFlexi(thisOctal%twoD)
         case("oneD")
            call receiveSingleFlexi(thisOctal%oneD)
         case("maxChildren")
            call receiveSingleFlexi(thisOctal%maxChildren)
         case("cylindrical")
            call receiveSingleFlexi(thisOctal%cylindrical)
         case("splitAzimuthally")
            call receiveSingleFlexi(thisOctal%splitAzimuthally)
         case("phi")
            call receiveSingleFlexi(thisOctal%phi)
         case("dphi")
            call receiveSingleFlexi(thisOctal%dphi)
         case("phimin")
            call receiveSingleFlexi(thisOctal%phimin)
         case("phimax")
            call receiveSingleFlexi(thisOctal%phimax)
         case("r")
            call receiveSingleFlexi(thisOctal%r)

         case("xMax")
            call receiveSingleFlexi(thisOctal%xMax)
         case("yMax")
            call receiveSingleFlexi(thisOctal%yMax)
         case("zMax")
            call receiveSingleFlexi(thisOctal%zMax)
         case("xMin")
            call receiveSingleFlexi(thisOctal%xMin)
         case("yMin")
            call receiveSingleFlexi(thisOctal%yMin)
         case("zMin")
            call receiveSingleFlexi(thisOctal%zMin)

         case("parentSubcell")
            call receiveSingleFlexi(thisOctal%parentSubcell)
         case("inStar")
            call receiveArrayFlexi(thisOctal%inStar)
         case("inFlow")
            call receiveArrayFlexi(thisOctal%inFlow)
         case("velocity")
            call receiveArrayFlexi(thisOctal%velocity)
         case("cornervelocity")
            call receivePointerFlexi(thisOctal%cornervelocity)
         case("cornerrho")
            call receivePointerFlexi(thisOctal%cornerrho)
         case("chiLine")
            call receivePointerFlexi(thisOctal%chiLine)
         case("etaLine")
            call receivePointerFlexi(thisOctal%etaLine)
         case("etaCont")
            call receivePointerFlexi(thisOctal%etaCont)
         case("biasLine3D")
            call receivePointerFlexi(thisOctal%biasLine3D)
         case("biasCont3D")
            call receivePointerFlexi(thisOctal%biasCont3D)
         case("probDistLine")
            call receivePointerFlexi(thisOctal%probDistLine)
         case("probDistCont")
            call receivePointerFlexi(thisOctal%probDistCont)
         case("ne")
            call receivePointerFlexi(thisOctal%ne)
         case("nH")
            call receivePointerFlexi(thisOctal%nH)
         case("nTot")
            call receivePointerFlexi(thisOctal%nTot)
         case("dustType")
            call receivePointerFlexi(thisOctal%dustType)
         case("dustTypeFraction")
            call receivePointerFlexi(thisOctal%dustTypeFraction)
         case("oldFrac")
            call receivePointerFlexi(thisOctal%oldFrac)
         case("scatteredIntensity")
            call receivePointerFlexi(thisOctal%scatteredIntensity)
         case("meanIntensity")
            call receivePointerFlexi(thisOctal%meanIntensity)
         case("mpiThread")
            call receiveArrayFlexi(thisOctal%mpiThread)
         case("kappaAbs")
            call receivePointerFlexi(thisOctal%kappaAbs)
         case("kappaSca")
            call receivePointerFlexi(thisOctal%kappaSca)
         case("ionFrac")
            call receivePointerFlexi(thisOctal%ionFrac)
         case("photoIonCoeff")
            call receivePointerFlexi(thisOctal%photoIonCoeff)
         case("distanceGrid")
            call receivePointerFlexi(thisOctal%distanceGrid)

         case("nCrossings")
            call receivePointerFlexi(thisOctal%nCrossings)
         case("hHeating")
            call receivePointerFlexi(thisOctal%hHeating)
         case("heHeating")
            call receivePointerFlexi(thisOctal%heHeating)

         case("undersampled")
            call receivePointerFlexi(thisOctal%undersampled)
         case("nDiffusion")
            call receivePointerFlexi(thisOctal%nDiffusion)
         case("diffusionApprox")
            call receivePointerFlexi(thisOctal%diffusionApprox)


         case("molecularLevel")
            call receivePointerFlexi(thisOctal%molecularLevel)
         case("jnu")
            call receivePointerFlexi(thisOctal%jnu)
         case("bnu")
            call receivePointerFlexi(thisOctal%bnu)
         case("molAbundance")
            call receivePointerFlexi(thisOctal%molAbundance)
         case("nh2")
            call receivePointerFlexi(thisOctal%nh2)
         case("microTurb")
            call receivePointerFlexi(thisOctal%microTurb)
         case("N")
            call receivePointerFlexi(thisOctal%n)
         case("departCoeff")
            call receivePointerFlexi(thisOctal%departCoeff)
         case("atomAbundance")
            call receivePointerFlexi(thisOctal%atomAbundance)
         case("atomLevel")
            call receivePointerFlexi(thisOctal%atomLevel)
         case("jnuCont")
            call receivePointerFlexi(thisOctal%jnuCont)
         case("jnuLine")
            call receivePointerFlexi(thisOctal%jnuLine)

         case("q_i")
            call receivePointerFlexi(thisOctal%q_i)
         case("q_i_plus_1")
            call receivePointerFlexi(thisOctal%q_i_plus_1)
         case("q_i_minus_1")
            call receivePointerFlexi(thisOctal%q_i_minus_1)
         case("q_i_minus_2")
            call receivePointerFlexi(thisOctal%q_i_minus_2)

         case("x_i")
            call receivePointerFlexi(thisOctal%x_i)
         case("x_i_plus_1")
            call receivePointerFlexi(thisOctal%x_i_plus_1)
         case("x_i_minus_1")
            call receivePointerFlexi(thisOctal%x_i_minus_1)

         case("x_i_minus_2")
            call receivePointerFlexi(thisOctal%x_i_minus_2)

         case("u_interface")
            call receivePointerFlexi(thisOctal%u_interface)
         case("u_i")
            call receivePointerFlexi(thisOctal%u_i)
         case("u_i_plus_1")
            call receivePointerFlexi(thisOctal%u_i_plus_1)
         case("u_i_minus_1")
            call receivePointerFlexi(thisOctal%u_i_minus_1)

         case("flux_i")
            call receivePointerFlexi(thisOctal%flux_i)
         case("flux_i_plus_1")
            call receivePointerFlexi(thisOctal%flux_i_plus_1)
         case("flux_i_minus_1")
            call receivePointerFlexi(thisOctal%flux_i_minus_1)


         case("phiLimit")
            call receivePointerFlexi(thisOctal%phiLimit)

         case("ghostCell")
            call receivePointerFlexi(thisOctal%ghostCell)
            
         case("corner")
            call receivePointerFlexi(thisOctal%corner)


         case("feederCell")
            call receivePointerFlexi(thisOctal%feederCell)
         case("edgeCell")
            call receivePointerFlexi(thisOctal%edgeCell)
         case("refinedLastTime")
            call receivePointerFlexi(thisOctal%refinedLastTime)

         case("pressure_i")
            call receivePointerFlexi(thisOctal%pressure_i)
         case("pressure_i_plus_1")
            call receivePointerFlexi(thisOctal%pressure_i_plus_1)
         case("pressure_i_minus_1")
            call receivePointerFlexi(thisOctal%pressure_i_minus_1)

         case("rhou")
            call receivePointerFlexi(thisOctal%rhou)
         case("rhov")
            call receivePointerFlexi(thisOctal%rhov)
         case("rhow")
            call receivePointerFlexi(thisOctal%rhow)

         case("rhorv_i_plus_1")
            call receivePointerFlexi(thisOctal%rhorv_i_plus_1)
         case("rhorv_i_minus_1")
            call receivePointerFlexi(thisOctal%rhorv_i_minus_1)



         case("rhoe")
            call receivePointerFlexi(thisOctal%rhoe)
         case("rhoeLast")
            call receivePointerFlexi(thisOctal%rhoeLastTime)
         case("energy")
            call receivePointerFlexi(thisOctal%energy)

         case("qViscosity")
            call receivePointerFlexi(thisOctal%qViscosity)


         case("phi_i")
            call receivePointerFlexi(thisOctal%phi_i)
         case("phi_stars")
            call receivePointerFlexi(thisOctal%phi_stars)
         case("phi_gas")
            call receivePointerFlexi(thisOctal%phi_gas)
         case("phi_i_plus_1")
            call receivePointerFlexi(thisOctal%phi_i_plus_1)
         case("phi_i_minus_1")
            call receivePointerFlexi(thisOctal%phi_i_minus_1)

         case("rho_i_plus_1")
            call receivePointerFlexi(thisOctal%rho_i_plus_1)
         case("rho_i_minus_1")
            call receivePointerFlexi(thisOctal%rho_i_minus_1)

         case("boundaryCondition")
            call receivePointerFlexi(thisOctal%boundaryCondition)

         case("boundaryCell")
            call receivePointerFlexi(thisOctal%boundaryCell)

         case("boundaryPartner")
            call receivePointerFlexi(thisOctal%boundaryPartner)
         case("gravboundaryPartner")
            call receivePointerFlexi(thisOctal%gravboundaryPartner)
         case("radiationMomentum")
            call receivePointerFlexi(thisOctal%radiationMomentum)
         case("kappaTimesFlux")
            call receivePointerFlexi(thisOctal%kappaTimesFlux)
         case("UVvector")
            call receivePointerFlexi(thisOctal%UVvector)
         case("changed")
            call receivePointerFlexi(thisOctal%changed)
         case("rLimit")
            call receivePointerFlexi(thisOctal%rLimit)

         case("iEquationOfState")
            call receivePointerFlexi(thisOctal%iEquationOfState)

         case("iAnalyticalVelocity")
            call receivePointerFlexi(thisOctal%iAnalyticalVelocity)

         case("gamma")
            call receivePointerFlexi(thisOctal%gamma)

         case("udens")
            call receivePointerFlexi(thisOctal%udens)
         case("adot")
            call receivePointerFlexi(thisOctal%aDot)
         case("adotdist")
            call receivePointerFlexi(thisOctal%distanceGridaDot)
         case("dfromgas")
            call receivePointerFlexi(thisOctal%distanceGridPhotonFromGas)
         case("dfromsource")
            call receivePointerFlexi(thisOctal%distanceGridPhotonFromSource)
         case("ufromgas")
            call receivePointerFlexi(thisOctal%photonEnergyDensityFromGas)
         case("ufromsource")
            call receivePointerFlexi(thisOctal%photonEnergyDensityFromSource)
         case("utotal")
            call receivePointerFlexi(thisOctal%photonEnergyDensityFromSource)
         case("oldutotal")
            call receivePointerFlexi(thisOctal%oldphotonEnergyDensity)
         case("fixedtemperature")
            call receivePointerFlexi(thisOctal%fixedtemperature)


         case DEFAULT
            write(message,*) "Unrecognised tag on receive: "//trim(tag)
            call writeFatal(message)
            stop
         end select

      end do
    end subroutine receiveOctalViaMPI

    subroutine sendOctalViaMPI(thisOctal, ithread)
     use mpi
     type(OCTAL), pointer :: thisOctal
     integer, parameter :: mpitag = 50
     character(len=20) :: tmp
     integer :: ierr
     integer :: iThread

      tmp = "OCTALBEGINS"
      call MPI_SEND(tmp, 20, MPI_CHARACTER, iThread, mpitag, localWorldCommunicator, ierr)
      call sendAttributeStaticFlexi(iThread, "nDepth", thisOctal%nDepth)
      call sendAttributeStaticFlexi(iThread, "nChildren", thisOctal%nChildren)
      call sendAttributeStaticFlexi(iThread, "indexChild", thisOCtal%IndexChild)
      call sendAttributeStaticFlexi(iThread, "hasChild", thisOctal%HasChild)
      call sendAttributeStaticFlexi(iThread, "centre", thisOctal%centre)
      call sendAttributeStaticFlexi(iThread, "rho", thisOctal%rho)
      call sendAttributeStaticFlexi(iThread, "temperature", thisOctal%temperature)
      call sendAttributeStaticFlexi(iThread, "label", thisOctal%label)
      call sendAttributeStaticFlexi(iThread, "subcellSize", thisOctal%subcellSize)
      call sendAttributeStaticFlexi(iThread, "threeD", thisOctal%threeD)
      call sendAttributeStaticFlexi(iThread, "twoD", thisOctal%twoD)
      call sendAttributeStaticFlexi(iThread, "oneD", thisOctal%oneD)
      call sendAttributeStaticFlexi(iThread, "maxChildren", thisOctal%maxChildren)
      call sendAttributeStaticFlexi(iThread, "cylindrical", thisOctal%cylindrical)
      call sendAttributeStaticFlexi(iThread, "splitAzimuthally", thisOctal%splitAzimuthally)
      call sendAttributeStaticFlexi(iThread, "phi", thisOctal%phi)
      call sendAttributeStaticFlexi(iThread, "dphi", thisOctal%dphi)
      call sendAttributeStaticFlexi(iThread, "phimin", thisOctal%phimin)
      call sendAttributeStaticFlexi(iThread, "phimax", thisOctal%phimax)
      call sendAttributeStaticFlexi(iThread, "r", thisOctal%r)
      call sendAttributeStaticFlexi(iThread, "parentSubcell", thisOctal%parentSubcell)
      call sendAttributeStaticFlexi(iThread, "inStar", thisOctal%inStar)
      call sendAttributeStaticFlexi(iThread, "inFlow", thisOctal%inFlow)
      call sendAttributeStaticFlexi(iThread, "velocity", thisOctal%velocity)

      call sendAttributeStaticFlexi(iThread, "xMax", thisOctal%xMax)
      call sendAttributeStaticFlexi(iThread, "yMax", thisOctal%yMax)
      call sendAttributeStaticFlexi(iThread, "zMax", thisOctal%zMax)
      call sendAttributeStaticFlexi(iThread, "xMin", thisOctal%xMin)
      call sendAttributeStaticFlexi(iThread, "yMin", thisOctal%yMin)
      call sendAttributeStaticFlexi(iThread, "zMin", thisOctal%zMin)

      call sendAttributePointerFlexi(iThread, "cornervelocity", thisOctal%cornervelocity)
      call sendAttributePointerFlexi(iThread, "chiLine", thisOctal%chiLine)
      call sendAttributePointerFlexi(iThread, "etaLine", thisOctal%etaLine)
      call sendAttributePointerFlexi(iThread, "etaCont", thisOctal%etaCont)
      call sendAttributePointerFlexi(iThread, "biasLine3D", thisOctal%biasLine3D)
      call sendAttributePointerFlexi(iThread, "biasCont3D", thisOctal%biasCont3D)
      call sendAttributePointerFlexi(iThread, "probDistLine", thisOctal%probDistLine)
      call sendAttributePointerFlexi(iThread, "probDistCont", thisOctal%probDistCont)
      call sendAttributePointerFlexi(iThread, "ne", thisOctal%ne)
      call sendAttributePointerFlexi(iThread, "nH", thisOctal%nH)
      call sendAttributePointerFlexi(iThread, "nTot", thisOctal%nTot)
      call sendAttributePointerFlexi(iThread, "dustType", thisOctal%dustType)


      call sendAttributePointerFlexi(iThread, "kappaAbs", thisOctal%kappaAbs)
      call sendAttributePointerFlexi(iThread, "kappaSca", thisOctal%kappaSca)

      call sendAttributePointerFlexi(iThread, "ionFrac", thisOctal%ionFrac)
      call sendAttributePointerFlexi(iThread, "photoIonCoeff", thisOctal%photoIonCoeff)

      call sendAttributePointerFlexi(iThread, "distanceGrid", thisOctal%distanceGrid)

      call sendAttributePointerFlexi(iThread, "nCrossings", thisOctal%nCrossings)
      call sendAttributePointerFlexi(iThread, "hHeating", thisOctal%hHeating)
      call sendAttributePointerFlexi(iThread, "heHeating", thisOctal%heHeating)
      call sendAttributePointerFlexi(iThread, "undersampled", thisOctal%undersampled)
      call sendAttributePointerFlexi(iThread, "nDiffusion", thisOctal%nDiffusion)
      call sendAttributePointerFlexi(iThread, "diffusionApprox", thisOctal%diffusionApprox)

      call sendAttributePointerFlexi(iThread, "molecularLevel", thisOctal%molecularLevel)
      call sendAttributePointerFlexi(iThread, "jnu", thisOctal%jnu)
      call sendAttributePointerFlexi(iThread, "bnu", thisOctal%bnu)
      call sendAttributePointerFlexi(iThread, "molAbundance", thisOctal%molAbundance)
      call sendAttributePointerFlexi(iThread, "nh2", thisOctal%nh2)
      call sendAttributePointerFlexi(iThread, "microTurb", thisOctal%microTurb)
      call sendAttributePointerFlexi(iThread, "cornerrho", thisOctal%cornerrho)

      call sendAttributePointerFlexi(iThread, "N", thisOctal%n)
      call sendAttributePointerFlexi(iThread, "departCoeff", thisOctal%departCoeff)
      call sendAttributePointerFlexi(iThread, "dustTypeFraction", thisOctal%dustTypeFraction)
      call sendAttributePointerFlexi(iThread, "oldFrac", thisOctal%oldFrac)
      call sendAttributePointerFlexi(iThread, "scatteredIntensity", thisOctal%scatteredIntensity)

      call sendAttributePointerFlexi(iThread, "meanIntensity", thisOctal%meanIntensity)


      call sendAttributePointerFlexi(iThread, "atomAbundance", thisOctal%atomAbundance)
      call sendAttributePointerFlexi(iThread, "atomLevel", thisOctal%atomLevel)
      call sendAttributePointerFlexi(iThread, "jnuCont", thisOctal%jnuCont)
      call sendAttributePointerFlexi(iThread, "jnuLine", thisOctal%jnuLine)

      call sendAttributePointerFlexi(iThread, "q_i", thisOctal%q_i)
      call sendAttributePointerFlexi(iThread, "q_i_plus_1", thisOctal%q_i_plus_1)
      call sendAttributePointerFlexi(iThread, "q_i_minus_1", thisOctal%q_i_minus_1)
      call sendAttributePointerFlexi(iThread, "q_i_minus_2", thisOctal%q_i_minus_2)

      call sendAttributePointerFlexi(iThread, "x_i", thisOctal%x_i)
      call sendAttributePointerFlexi(iThread, "x_i_plus_1", thisOctal%x_i_plus_1)
      call sendAttributePointerFlexi(iThread, "x_i_minus_1", thisOctal%x_i_minus_1)
      call sendAttributePointerFlexi(iThread, "x_i_minus_2", thisOctal%x_i_minus_2)

      call sendAttributePointerFlexi(iThread, "u_interface", thisOctal%u_interface)
      call sendAttributePointerFlexi(iThread, "u_i", thisOctal%u_i)
      call sendAttributePointerFlexi(iThread, "u_i_plus_1", thisOctal%u_i_plus_1)
      call sendAttributePointerFlexi(iThread, "u_i_minus_1", thisOctal%u_i_minus_1)

      call sendAttributePointerFlexi(iThread, "flux_i", thisOctal%flux_i)
      call sendAttributePointerFlexi(iThread, "flux_i_plus_1", thisOctal%flux_i_plus_1)
      call sendAttributePointerFlexi(iThread, "flux_i_minus_1", thisOctal%flux_i_minus_1)


      call sendAttributePointerFlexi(iThread, "phiLimit", thisOctal%phiLimit)

      call sendAttributePointerFlexi(iThread, "ghostCell", thisOctal%ghostCell)
      call sendAttributePointerFlexi(iThread, "corner", thisOctal%corner)
      call sendAttributePointerFlexi(iThread, "feederCell", thisOctal%feederCell)
      call sendAttributePointerFlexi(iThread, "edgeCell", thisOctal%edgeCell)
      call sendAttributePointerFlexi(iThread, "refinedLastTime", thisOctal%refinedLastTime)

      call sendAttributePointerFlexi(iThread, "pressure_i", thisOctal%pressure_i)
      call sendAttributePointerFlexi(iThread, "pressure_i_plus_1", thisOctal%pressure_i_plus_1)
      call sendAttributePointerFlexi(iThread, "pressure_i_minus_1", thisOctal%pressure_i_minus_1)

      call sendAttributePointerFlexi(iThread, "rhou", thisOctal%rhou)
      call sendAttributePointerFlexi(iThread, "rhov", thisOctal%rhov)
      call sendAttributePointerFlexi(iThread, "rhow", thisOctal%rhow)

      call sendAttributePointerFlexi(iThread, "rhorv_i_minus_1", thisOctal%rhorv_i_minus_1)
      call sendAttributePointerFlexi(iThread, "rhorv_i_plus_1", thisOctal%rhorv_i_plus_1)


      call sendAttributePointerFlexi(iThread, "rhoe", thisOctal%rhoe)
      call sendAttributePointerFlexi(iThread, "rhoeLast", thisOctal%rhoeLastTime)
      call sendAttributePointerFlexi(iThread, "energy", thisOctal%energy)

      call sendAttributePointerFlexi(iThread, "qViscosity", thisOctal%qViscosity)


      call sendAttributePointerFlexi(iThread, "phi_i", thisOctal%phi_i)
      call sendAttributePointerFlexi(iThread, "phi_gas", thisOctal%phi_gas)
      call sendAttributePointerFlexi(iThread, "phi_stars", thisOctal%phi_stars)
      call sendAttributePointerFlexi(iThread, "phi_i_plus_1", thisOctal%phi_i_plus_1)
      call sendAttributePointerFlexi(iThread, "phi_i_minus_1", thisOctal%phi_i_minus_1)

      call sendAttributePointerFlexi(iThread, "rho_i_plus_1", thisOctal%rho_i_plus_1)
      call sendAttributePointerFlexi(iThread, "rho_i_minus_1", thisOctal%rho_i_minus_1)


      call sendAttributePointerFlexi(iThread, "boundaryCondition", thisOctal%boundaryCondition)


      call sendAttributePointerFlexi(iThread, "boundaryPartner", thisOctal%boundaryPartner)

      call sendAttributePointerFlexi(iThread, "boundaryCell", thisOctal%boundaryCell)

      call sendAttributePointerFlexi(iThread, "gravboundaryPartner", thisOctal%GravboundaryPartner)
      call sendAttributePointerFlexi(iThread, "radiationMomentum", thisOctal%radiationMomentum)
      call sendAttributePointerFlexi(iThread, "kappaTimesFlux", thisOctal%kappaTimesFlux)
      call sendAttributePointerFlexi(iThread, "UVvector", thisOctal%UVvector)
      call sendAttributePointerFlexi(iThread, "changed", thisOctal%changed)
      call sendAttributePointerFlexi(iThread, "rLimit", thisOctal%rLimit)

      call sendAttributePointerFlexi(iThread, "iEquationOfState", thisOctal%iEquationOfState)

      call sendAttributePointerFlexi(iThread, "iAnalyticalVelocity", thisOctal%iAnalyticalVelocity)

      call sendAttributePointerFlexi(iThread, "gamma", thisOctal%gamma)


      call sendAttributePointerFlexi(iThread, "udens", thisOctal%udens)
      call sendAttributePointerFlexi(iThread, "adot", thisOctal%adot)
      call sendAttributePointerFlexi(iThread, "adotdist", thisOctal%distanceGridAdot)
      call sendAttributePointerFlexi(iThread, "dfromgas", thisOctal%distanceGridPhotonFromGas)
      call sendAttributePointerFlexi(iThread, "dfromsource", thisOctal%distanceGridPhotonFromSource)
      call sendAttributePointerFlexi(iThread, "ufromgas", thisOctal%photonEnergyDensityFromGas)
      call sendAttributePointerFlexi(iThread, "ufromsource", thisOctal%photonEnergyDensityFromSource)
      call sendAttributePointerFlexi(iThread, "utotal", thisOctal%photonEnergyDensity)
      call sendAttributePointerFlexi(iThread, "oldutotal", thisOctal%oldphotonEnergyDensity)

      call sendAttributePointerFlexi(iThread, "fixedtemperature", thisOctal%fixedTemperature)

      call sendAttributeStaticFlexi(iThread, "mpiThread", thisOctal%mpiThread)
      tmp = "OCTALENDS"
      call MPI_SEND(tmp, 20, MPI_CHARACTER, iThread, mpitag, localWorldCommunicator, ierr)


    end subroutine sendOctalViaMPI

#endif

 end module gridio_mod
