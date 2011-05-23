  subroutine fillGridDisk2(grid, rhoNought, rCore, rInner, rOuter, height, &
       mCore, diskTemp)


    implicit none
    type(GRIDTYPE) :: grid
    real :: rOuter, rInner, rCore
    type(VECTOR) :: rVec, vVec
    real :: r
    real :: rho, rhoNought
    real :: diskTemp
    real :: height
    real :: vel, mCore
    type(VECTOR) :: spinAxis
    integer :: i, j, k 

    grid%geometry = "disk"
    grid%lineEmission = .true.
    grid%rCore = rCore


    spinAxis = VECTOR(0.,0.,1.)

    grid%rho = 0.
    grid%etaLine = 1.e-30
    grid%etaCont = 1.e-30
    grid%chiLine = 1.e-30
    grid%kappaSca = 1.e-30
    grid%kappaAbs = 1.e-30
    grid%velocity = VECTOR(1.e-30,1.e-30,1.e-30)
    grid%temperature = diskTemp

    do i = 1, grid%nx
       grid%xAxis(i) = 2.*real(i-1)/real(grid%nx-1)-1.
    enddo
    do i = 1, grid%ny
       grid%yAxis(i) = 2.*real(i-1)/real(grid%ny-1)-1.
    enddo
    do i = 1, grid%nz
       grid%zAxis(i) = 2.*real(i-1)/real(grid%nz-1)-1.
    enddo

    grid%xAxis(1:grid%nx) = grid%xAxis(1:grid%nx) * rOuter
    grid%yAxis(1:grid%ny) = grid%yAxis(1:grid%ny) * rOuter
    grid%zAxis(1:grid%nz) = grid%zAxis(1:grid%nz) * 5.*height
    do i = 1, grid%nx
       write(*,*) i ,grid%xAxis(i), grid%yaxis(i), grid%zaxis(i)
    enddo

    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             rVec = VECTOR(grid%xAxis(i),grid%yAxis(j),grid%zAxis(k))
             r = modulus(rVec)
             if ((r > rInner).and. (r < rOuter)) then
                rho = rhoNought * (rInner/r)**2
                grid%rho(i,j,k) = rho * exp(-abs(grid%zAxis(k)/height))
                rVec = VECTOR(grid%xAxis(i),grid%yAxis(j),grid%zAxis(k))
                rVec = rVec / dble(r)
                vVec = rVec  .cross. spinAxis
                if (modulus(vVec) /= 0.) then
                   call normalize(vVec)
                   vel = sqrt(bigG * mCore / r) / cSpeed
                   vVec = dble(vel)  * vVec 
                   grid%velocity(i,j,k) = vVec
                endif
             endif
          enddo
       enddo
    enddo

    
    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             if (grid%rho(i,j,k) /= 0.) then
                grid%kappaSca(i,j,k,1) = max(1.e-30,grid%rho(i,j,k) * real(sigmaE))
             endif
          enddo
       enddo
    enddo

  end subroutine fillGridDisk2

  subroutine writeAMRgridOld(filename,fileFormatted,grid)
    use input_variables, only : molecular, cmf
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
            
       ! write a time stamp to the file
       write(unit=20,fmt=*,iostat=error) timeValues(:)
       if (error /=0) then
         print *, 'Panic: write error in writeAMRgrid (formatted timeValues)' ; stop
       end if
           
       ! write the variables that are stored in the top-level 'grid' structure
       write(unit=20,fmt=*,iostat=error) grid%nLambda, grid%flatSpec, grid%adaptive,& 
               grid%cartesian, grid%isotropic, grid%hitCore, grid%diskRadius, &
               grid%diskNormal, grid%DipoleOffset, grid%geometry, grid%rCore, &
               grid%lCore, grid%chanceWindOverTotalContinuum,                 &
               grid%lineEmission, grid%contEmission, grid%doRaman,            &
               grid%resonanceLine, grid%rStar1, grid%rStar2, grid%lumRatio,   &
               grid%tempSource, grid%starPos1, grid%starPos2, grid%lambda2,   &
               grid%maxLevels, grid%maxDepth, grid%halfSmallestSubcell,       &
               grid%nOctals, grid%smoothingFactor, grid%oneKappa, grid%rInner,&
               grid%rOuter, grid%amr2dOnly, grid%photoionization, grid%iDump, grid%currentTime
       if (error /=0) then
         print *, 'Panic: write error in writeAMRgrid (formatted variables)' ; stop
       end if

!       if (cmf) then
!          write(unit=20,fmt=*) grid%nFreqArray, grid%freqArray(1:2000)
!       endif
               
    else
            
       ! write a time stamp to the file
       write(unit=20,iostat=error) timeValues(:)
       if (error /=0) then
         print *, 'Panic: write error in writeAMRgrid (unformatted timeValues)' ; stop
       end if
    
       ! write the variables that are stored in the top-level 'grid' structure
       write(unit=20,iostat=error) grid%nLambda, grid%flatSpec, grid%adaptive,& 
               grid%cartesian, grid%isotropic, grid%hitCore, grid%diskRadius, &
               grid%diskNormal, grid%DipoleOffset, grid%geometry, grid%rCore, &
               grid%lCore, grid%chanceWindOverTotalContinuum,                 &
               grid%lineEmission, grid%contEmission, grid%doRaman,            &
               grid%resonanceLine, grid%rStar1, grid%rStar2, grid%lumRatio,   &
               grid%tempSource, grid%starPos1, grid%starPos2, grid%lambda2,   &
               grid%maxLevels, grid%maxDepth, grid%halfSmallestSubcell,       &
               grid%nOctals, grid%smoothingFactor, grid%oneKappa, grid%rInner,& 
               grid%rOuter, grid%amr2dOnly, grid%photoionization, grid%iDump, grid%currentTime
       if (error /=0) then
         print *, 'Panic: write error in writeAMRgrid (unformatted variables)' ; stop
       end if

!       if (cmf) then
!          write(unit=20) grid%nFreqArray, grid%freqArray(1:2000)
!       endif

               
    end if 

    
    call writeReal1D(grid%lamarray,fileFormatted)
    call writeReal2D(grid%oneKappaAbs,fileFormatted)
    call writeReal2D(grid%oneKappaSca,fileFormatted)
    call writeClumps(fileFormatted)

    ! now we call the recursive subroutine to store the tree structure 
    if (associated(grid%octreeRoot)) then
       if (fileFormatted) then 
          write(unit=20,fmt=*) .true.
       else
          write(unit=20) .true.
       end if
       call writeOctreePrivate(grid%octreeRoot,fileFormatted, grid)
    else 
       if (fileFormatted) then 
          write(unit=20,fmt=*) .false.
       else
          write(unit=20) .false.
       end if
    end if

    endfile 20
    close(unit=20)
    
  contains
  
    recursive subroutine writeOctreePrivate(thisOctal,fileFormatted, grid)
      use octal_mod, only: writeattributepointer, writeattributestatic
       ! writes out an octal from the grid octree

       type(octal), intent(in), target :: thisOctal
       logical, intent(in)             :: fileFormatted
       type(gridtype) :: grid
       type(octal), pointer :: thisChild
       integer              :: iChild

       call writeAttributeStatic(20, thisOctal%nDepth, fileFormatted)
       call writeAttributeStatic(20, thisOctal%nChildren, fileFormatted)
       call writeAttributeStatic(20, thisOctal%indexChild, fileFormatted)
       call writeAttributeStatic(20, thisOctal%hasChild, fileFormatted)
       call writeAttributeStatic(20, thisOctal%centre, fileFormatted)
       call writeAttributeStatic(20, thisOctal%rho, fileFormatted)
       call writeAttributeStatic(20, thisOctal%temperature, fileFormatted)
       call writeAttributeStatic(20, thisOctal%label, fileFormatted)
       call writeAttributeStatic(20, thisOctal%subcellSize, fileFormatted)
       call writeAttributeStatic(20, thisOctal%threeD, fileFormatted)
       call writeAttributeStatic(20, thisOctal%twoD, fileFormatted)
       call writeAttributeStatic(20, thisOctal%maxChildren, fileFormatted)
       call writeAttributeStatic(20, thisOctal%cylindrical, fileFormatted)
       call writeAttributeStatic(20, thisOctal%splitAzimuthally, fileFormatted)
       call writeAttributeStatic(20, thisOctal%phi, fileFormatted)
       call writeAttributeStatic(20, thisOctal%dphi, fileFormatted)
       call writeAttributeStatic(20, thisOctal%r, fileFormatted)
       call writeAttributeStatic(20, thisOctal%parentSubcell, fileFormatted)
       call writeAttributeStatic(20, thisOctal%inStar, fileFormatted)
       call writeAttributeStatic(20, thisOctal%inFlow, fileFormatted)
       call writeAttributeStatic(20, thisOctal%velocity, fileFormatted)
       call writeAttributeStatic(20, thisOctal%cornervelocity, fileFormatted)

       call writeAttributePointer(20, thisOctal%chiLine, fileFormatted)
       call writeAttributePointer(20, thisOctal%etaLine, fileFormatted)
       call writeAttributePointer(20, thisOctal%etaCont, fileFormatted)
       call writeAttributePointer(20, thisOctal%biasLine3D, fileFormatted)
       call writeAttributePointer(20, thisOctal%biasCont3D, fileFormatted)
       call writeAttributePointer(20, thisOctal%probDistLine, fileFormatted)
       call writeAttributePointer(20, thisOctal%probDistCont, fileFormatted)
       call writeAttributePointer(20, thisOctal%ne, fileFormatted)
       call writeAttributePointer(20, thisOctal%nH, fileFormatted)
       call writeAttributePointer(20, thisOctal%nTot, fileFormatted)
       call writeAttributePointer(20, thisOctal%dustType, fileFormatted)


       
!       if (fileFormatted) then 
!          write(iostat=error,fmt=*,unit=20) thisOctal%nDepth, thisOctal%nChildren,&
!                  thisOctal%indexChild, thisOctal%hasChild, thisOctal%centre,& 
!                  thisOctal%rho, thisOctal%velocity,thisOctal%cornerVelocity,&
!                  thisOctal%temperature, thisOctal%chiLine,thisOctal%etaLine,&
!                  thisOctal%etaCont, thisOctal%biasLine3D,                   &
!                  thisOctal%biasCont3D, thisOctal%probDistLine,              &
!                  thisOctal%probDistCont, thisOctal%Ne,  thisOctal%nh,thisOctal%nTot,      &
!                  thisOctal%inStar, thisOctal%inFlow, thisOctal%label,       &
!                  thisOctal%subcellSize,thisOctal%threed, thisOctal%twoD,    &
!                  thisOctal%maxChildren, thisOctal%dustType,  &
!                  thisOctal%cylindrical, thisOctal%splitAzimuthally, thisOctal%phi, &
!                  thisOctal%dPhi, thisOctal%r, thisOctal%parentSubcell
!       else
!          write(iostat=error,unit=20) thisOctal%nDepth, thisOctal%nChildren, &
!                  thisOctal%indexChild, thisOctal%hasChild, thisOctal%centre,& 
!                  thisOctal%rho, thisOctal%velocity,thisOctal%cornerVelocity,&
!                  thisOctal%temperature, thisOctal%chiLine,thisOctal%etaLine,&
!                  thisOctal%etaCont, thisOctal%biasLine3D,                   &
!                  thisOctal%biasCont3D, thisOctal%probDistLine,              &
!                  thisOctal%probDistCont, thisOctal%Ne, thisOctal%nh, thisOctal%nTot,      &
!                  thisOctal%inStar, thisOctal%inFlow, thisOctal%label,       &
!                  thisOctal%subcellSize, thisOctal%threeD, thisOctal%twoD,   &
!                  thisOctal%maxChildren, thisOctal%dustType, &
!                  thisOctal%cylindrical, thisOctal%splitAzimuthally, thisOctal%phi, &
!                  thisOctal%dPhi, thisOctal%r, thisOctal%parentSubcell
!       end if 
       if (.not.grid%oneKappa) then
!          call writeReal2D(thisOctal%kappaAbs,fileFormatted)
!          call writeReal2D(thisOctal%kappaSca,fileFormatted)
          call writeDouble2D(thisOctal%kappaAbs,fileFormatted)
          call writeDouble2D(thisOctal%kappaSca,fileFormatted)
       endif
       if (grid%photoionization) then
          call writeDouble2D(thisOctal%ionFrac,fileFormatted)
          call writeDouble2D(thisOctal%photoIonCoeff,fileFormatted)
       endif
       if (molecular) then
          call writeDouble2D(thisOctal%molecularLevel,fileFormatted)
          call writeDouble2D(thisOctal%jnu,fileFormatted)
          call writeDouble2D(thisOctal%bnu,fileFormatted)
          call writeReal1D(thisOctal%molAbundance,fileFormatted)
          call writeDouble1D(thisOctal%nh2,fileFormatted)
          call writeDouble1D(thisOctal%microturb,fileFormatted)
       endif

       call writeDouble2D(thisOctal%N,fileFormatted)
       call writeReal2D(thisOctal%departCoeff,fileFormatted)
       call writeDouble2D(thisOctal%dustTypeFraction, fileFormatted)

       if (cmf) then
          call writeDouble2D(thisOctal%atomAbundance, fileFormatted)
          call writeDouble3D(thisOctal%atomLevel,fileFormatted)
          call writeDouble2D(thisOctal%jnuCont,fileFormatted)
          call writeDouble2D(thisOctal%jnuLine,fileFormatted)
          if (fileformatted) then
             write(unit=20,iostat=error,fmt=*) thisOctal%microturb
          else
             write(unit=20,iostat=error) thisOctal%microturb
          endif
       endif

       if (grid%splitOverMpi) then
          if (fileformatted) then
             write(unit=20,iostat=error,fmt=*) thisOctal%mpiThread
          else
             write(unit=20,iostat=error) thisOctal%mpiThread
          endif
       endif
       
       if (thisOctal%nChildren > 0) then 
          do iChild = 1, thisOctal%nChildren, 1
             thisChild => thisOctal%child(iChild)
             call writeOctreePrivate(thisChild,fileFormatted,grid)
          end do
       end if

    end subroutine writeOctreePrivate
    
    
  end subroutine writeAMRgridOld

  subroutine readAMRgridOld(filename,fileFormatted,grid)
    ! reads in a previously saved 'grid' for an adaptive mesh geometry  

    use input_variables, only: geometry,dipoleOffset,amr2dOnly,statEq2d, molecular, cmf, hydrodynamics
    use unix_mod, only: unixGetenv
    implicit none

    character(len=*)            :: filename
    character(len=80)            :: message
    logical, intent(in)         :: fileFormatted
    type(GRIDTYPE), intent(inout) :: grid
    
    integer, dimension(8) :: timeValues    ! system date and time
    integer               :: dummy         
    integer               :: error         ! status code
    logical               :: octreePresent ! true if grid has an octree
    integer :: nOctal
    character(len=80) :: absolutePath, inFile
    integer :: iJunk
    real, pointer :: junk(:), junk2(:,:), junk3(:,:)
    absolutePath = " "
    error = 0
    junk => null()
    junk2 => null()
    junk3 => null()

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

    write(message,'(a,a)') "Reading populations file from: ",trim(filename)
    call writeInfo(message,TRIVIAL)
    write(message,'(a,i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') ' - data file written at: ', &
                          timeValues(1),'/',timeValues(2),'/',&
                          timeValues(3),'  ',timeValues(5),':',timeValues(6)
    call writeInfo(message,TRIVIAL)
                          
    ! read the variables to be stored in the top-level 'grid' structure
    if (fileFormatted) then
       read(unit=20,fmt=*,iostat=error) ijunk, grid%flatSpec, grid%adaptive,& 
               grid%cartesian, grid%isotropic, grid%hitCore, grid%diskRadius, &
               grid%diskNormal, grid%DipoleOffset, grid%geometry, grid%rCore, &
               grid%lCore, grid%chanceWindOverTotalContinuum,                 &
               grid%lineEmission, grid%contEmission, grid%doRaman,            &
               grid%resonanceLine, grid%rStar1, grid%rStar2, grid%lumRatio,   &
               grid%tempSource, grid%starPos1, grid%starPos2, grid%lambda2,   &
               grid%maxLevels, grid%maxDepth, grid%halfSmallestSubcell,       &
               grid%nOctals, grid%smoothingFactor, grid%oneKappa, grid%rInner,&
               grid%rOuter, grid%amr2dOnly, grid%photoionization, grid%iDump, grid%currentTime
    else
       read(unit=20,iostat=error) ijunk, grid%flatSpec, grid%adaptive, & 
               grid%cartesian, grid%isotropic, grid%hitCore, grid%diskRadius, &
               grid%diskNormal, grid%DipoleOffset, grid%geometry, grid%rCore, &
               grid%lCore, grid%chanceWindOverTotalContinuum,                 &
               grid%lineEmission, grid%contEmission, grid%doRaman,            &
               grid%resonanceLine, grid%rStar1, grid%rStar2, grid%lumRatio,   &
               grid%tempSource, grid%starPos1, grid%starPos2, grid%lambda2,   &
               grid%maxLevels, grid%maxDepth, grid%halfSmallestSubcell,       &
               grid%nOctals, grid%smoothingFactor, grid%oneKappa, grid%rInner,& 
               grid%rOuter, grid%amr2dOnly, grid%photoionization, grid%iDump, grid%currentTime
    end if    

    if (error /=0) then
       print *, 'Panic: read error in readAMRgrid 1'
       stop
    end if
!    call readReal1D(grid%lamarray,fileFormatted)
!    call readReal2D(grid%oneKappaAbs,fileFormatted)
!    call readReal2D(grid%oneKappaSca,fileFormatted)


    call readReal1D(junk,fileFormatted)
    call readReal2D(junk2,fileFormatted)
    call readReal2D(junk3,fileFormatted)

!    deallocate(junk, junk2, junk3)
   
    Call readClumps(fileFormatted)
    
    ! now we call the recursive subroutine to read the tree structure 
    if (fileFormatted) then
       read(unit=20,fmt=*) octreePresent
    else
       read(unit=20) octreePresent
    end if
     
    if (octreePresent) then
       allocate(grid%octreeRoot)
       nOctal = 0
       call readOctreePrivate(grid%octreeRoot,null(),fileFormatted, nOctal, grid)
!       write(*,*) noctal,"octals read"
    end if

    ! check that we are at the end of the file
    if (fileFormatted) then
       read(unit=20,fmt=*, iostat=error) dummy
    else
       read(unit=20, iostat=error) dummy
    end if
    if (error == 0) then
       print *, 'Panic: read error (expected end of file) in readAMRgrid' ; stop
    end if
    
    close(unit=20)    grid%geometry = trim(geometry)
    grid%dipoleOffset = dipoleOffset
    grid%amr2donly = amr2donly
    grid%statEq2d = statEq2d

  contains
   
    recursive subroutine readOctreePrivate(thisOctal,parent,fileFormatted, noctal, grid)
       ! read in an octal to the grid octree
      use octal_mod, only: readattributepointer, readattributestatic

       implicit none
       type(octal), pointer :: thisOctal
       type(octal), pointer :: parent
       type(gridtype) :: grid

       logical, intent(in)  :: fileFormatted
       integer :: nOctal

       type(octal), pointer :: thisChild
       integer              :: iChild

       nOctal = nOctal+1
       thisOctal%parent => parent


       call readAttributeStatic(20, thisOctal%nDepth, fileFormatted)
       call readAttributeStatic(20, thisOctal%nChildren, fileFormatted)
       call readAttributeStatic(20, thisOctal%indexChild, fileFormatted)
       call readAttributeStatic(20, thisOctal%hasChild, fileFormatted)
       call readAttributeStatic(20, thisOctal%centre, fileFormatted)
       call readAttributeStatic(20, thisOctal%rho, fileFormatted)
       call readAttributeStatic(20, thisOctal%temperature, fileFormatted)
       call readAttributeStatic(20, thisOctal%label, fileFormatted)
       call readAttributeStatic(20, thisOctal%subcellSize, fileFormatted)
       call readAttributeStatic(20, thisOctal%threeD, fileFormatted)
       call readAttributeStatic(20, thisOctal%twoD, fileFormatted)
       call readAttributeStatic(20, thisOctal%maxChildren, fileFormatted)
       call readAttributeStatic(20, thisOctal%cylindrical, fileFormatted)
       call readAttributeStatic(20, thisOctal%splitAzimuthally, fileFormatted)
       call readAttributeStatic(20, thisOctal%phi, fileFormatted)
       call readAttributeStatic(20, thisOctal%dphi, fileFormatted)
       call readAttributeStatic(20, thisOctal%r, fileFormatted)
       call readAttributeStatic(20, thisOctal%parentSubcell, fileFormatted)
       call readAttributeStatic(20, thisOctal%inStar, fileFormatted)
       call readAttributeStatic(20, thisOctal%inFlow, fileFormatted)
       call readAttributeStatic(20, thisOctal%velocity, fileFormatted)
       call readAttributeStatic(20, thisOctal%cornervelocity, fileFormatted)

       call readAttributePointer(20, thisOctal%chiLine, fileFormatted)
       call readAttributePointer(20, thisOctal%etaLine, fileFormatted)
       call readAttributePointer(20, thisOctal%etaCont, fileFormatted)
       call readAttributePointer(20, thisOctal%biasLine3D, fileFormatted)
       call readAttributePointer(20, thisOctal%biasCont3D, fileFormatted)
       call readAttributePointer(20, thisOctal%probDistLine, fileFormatted)
       call readAttributePointer(20, thisOctal%probDistCont, fileFormatted)
       call readAttributePointer(20, thisOctal%ne, fileFormatted)
       call readAttributePointer(20, thisOctal%nH, fileFormatted)
       call readAttributePointer(20, thisOctal%nTot, fileFormatted)
       call readAttributePointer(20, thisOctal%dustType, fileFormatted)

       
!       if (fileFormatted) then
!          read(unit=20,fmt=*,iostat=error) thisOctal%nDepth, thisOctal%nChildren,  &
!                  thisOctal%indexChild, thisOctal%hasChild, thisOctal%centre,& 
!                  thisOctal%rho, thisOctal%velocity,thisOctal%cornerVelocity,&
!                  thisOctal%temperature, thisOctal%chiLine,thisOctal%etaLine,&
!                  thisOctal%etaCont, thisOctal%biasLine3D,                   &
!                  thisOctal%biasCont3D, thisOctal%probDistLine,              &
!                  thisOctal%probDistCont, thisOctal%Ne,  thisOctal%nh,thisOctal%nTot,      &
!                  thisOctal%inStar, thisOctal%inFlow, thisOctal%label,       &
!                  thisOctal%subcellSize, thisOctal%threeD, thisOctal%twoD,   &
!                  thisOctal%maxChildren, thisOctal%dustType, &
!                  thisOctal%cylindrical, thisOctal%splitAzimuthally, thisOctal%phi, &
!                  thisOctal%dPhi, thisOctal%r, thisOctal%parentSubcell
!       else
!          read(unit=20,iostat=error) thisOctal%nDepth, thisOctal%nChildren,  &
!                  thisOctal%indexChild, thisOctal%hasChild, thisOctal%centre,& 
!                  thisOctal%rho, thisOctal%velocity,thisOctal%cornerVelocity,&
!                  thisOctal%temperature, thisOctal%chiLine,thisOctal%etaLine,&
!                  thisOctal%etaCont, thisOctal%biasLine3D,                   &
!                  thisOctal%biasCont3D, thisOctal%probDistLine,              &
!                  thisOctal%probDistCont, thisOctal%Ne,  thisOctal%nh,thisOctal%nTot,      &
!                  thisOctal%inStar, thisOctal%inFlow, thisOctal%label,       &
!                  thisOctal%subcellSize, thisOctal%threeD,  thisOctal%twoD,  &
!                  thisOctal%maxChildren, thisOctal%dustType, &
!                  thisOctal%cylindrical, thisOctal%splitAzimuthally, thisOctal%phi, &
!                  thisOctal%dPhi, thisOctal%r, thisOctal%parentSubcell
!
!       end if
!       
       thisOctal%oneD = .not.(thisOctal%twoD.or.thisOctal%threeD)

       if (.not.grid%oneKappa) then
!          call readReal2D(thisOctal%kappaAbs,fileFormatted)
!          call readReal2D(thisOctal%kappaSca,fileFormatted)
          call readDouble2D(thisOctal%kappaAbs,fileFormatted)
          call readDouble2D(thisOctal%kappaSca,fileFormatted)
       endif
       if (grid%photoionization) then
          call readDouble2D(thisOctal%ionFrac,fileFormatted)
          call readDouble2D(thisOctal%photoIonCoeff,fileFormatted)
       endif
       if (molecular) then
          call readDouble2D(thisOctal%molecularLevel,fileFormatted)
          call readDouble2D(thisOctal%jnu,fileFormatted)
          call readDouble2D(thisOctal%bnu,fileFormatted)
          call readReal1D(thisOctal%molAbundance,fileFormatted)

          if (fileformatted) then
             read(unit=20,iostat=error,fmt=*) thisOctal%nh2
             read(unit=20,iostat=error,fmt=*) thisOctal%microturb
          else
             read(unit=20,iostat=error) thisOctal%nh2
             read(unit=20,iostat=error) thisOctal%microturb
          endif

       endif
       call readDouble2D(thisOctal%N,fileFormatted)
       call readReal2D(thisOctal%departCoeff,fileFormatted)
       call readDouble2D(thisOctal%dustTypeFraction, fileFormatted)

       if (cmf) then
          call readDouble2D(thisOctal%atomAbundance, fileFormatted)
          call readDouble3D(thisOctal%atomLevel,fileFormatted)
          call readDouble2D(thisOctal%jnuCont,fileFormatted)
          call readDouble2D(thisOctal%jnuLine,fileFormatted)

          if (fileformatted) then
             read(unit=20,iostat=error,fmt=*) thisOctal%microturb
          else
             read(unit=20,iostat=error) thisOctal%microturb
          endif

       endif

       if (hydrodynamics) then
          if (fileformatted) then
             write(unit=20,iostat=error,fmt=*) thisOctal%boundaryCondition
          else
             write(unit=20,iostat=error) thisOctal%boundaryCondition
          endif
       endif


       if (grid%splitOverMpi) then
          if (fileformatted) then
             read(unit=20,iostat=error,fmt=*) thisOctal%mpiThread
          else
             read(unit=20,iostat=error) thisOctal%mpiThread

          endif
       endif

       if (thisOctal%nChildren > 0) then 
          allocate(thisOctal%child(1:thisOctal%nChildren)) 
          do iChild = 1, thisOctal%nChildren, 1
             thisChild => thisOctal%child(iChild)
             call readOctreePrivate(thisChild,thisOctal,fileFormatted, nOctal, grid)
          end do
       end if

    end subroutine readOctreePrivate
 
  end subroutine readAMRgridOld

  logical function oddBall(grid,i1,i2,i3)
    type(GRIDTYPE) :: grid
    integer :: i1, i2, i3
    integer :: j1, j2, j3
    integer :: k1, k2, k3
    integer :: i, j, k
    oddBall = .true.
    
    j1 = max(i1-1,1)
    j2 = max(i2-1,1)
    j3 = max(i3-1,1)

    k1 = min(i1+1,grid%na1)
    k2 = min(i2+1,grid%na2)
    k3 = min(i3+1,grid%na3)

    do i = j1, k1
       do j = j2, k2
          do k = j3, k3
             if (grid%inUse(i,j,k)) then
                oddball = .false.
                EXIT
             endif
          enddo
          if (.not.oddBall) EXIT
       enddo
       if (.not.oddBall) EXIT
    enddo
  end function oddBall



  subroutine fillGridDonati2(grid, resonanceLine, misc)
    use utils_mod, only: hunt
    type(GRIDTYPE) :: grid
    integer :: i, j
    integer :: i1, i2, i3, j1, j2
    real :: zPrime, r, fac, ang
    logical :: resonanceLine
    character(len=*) :: misc
    type(VECTOR) :: rHat, rVec, vHat
    real :: surfaceVel, vel
    integer :: nz, ny
    real :: t1, t2, t3
    real :: thickness
    real, allocatable :: z(:), y(:), rho(:,:), t(:,:), phi(:,:), v(:,:)
    character(len=80) :: junkline
    integer, parameter :: ndisk = 8
    real :: rStar
    real :: diskDensity
    real :: xDisk(ndisk) = (/1.13, 1.22, 1.33, 1.49, 1.70, 2.00, 2.50, 3.0/)
    real :: yDisk(ndisk) = (/0.021, 0.055, 0.052, 0.031, 0.018, 0.010, 0.0046, &
                                                                   0.0025 /)
    real :: rhoDisk(ndisk) = (/5., 7.6, 9.5, 10.1, 10.0, 9.4, 8.5, 7.6 /)

    real :: yDisk2(nDisk) = (/0.004, 0.011, 0.025, 0.060, 0.090, 0.050, 0.023, 0.0125 /)
    real :: rhoDisk2(nDisk) = (/ 1.0, 1.5, 1.9, 2.0, 2.0, 1.9, 1.7, 1.5 /)


    grid%geometry = "donati"
    grid%lineEmission = .true.


    diskDensity = 1.5e12

    rStar = 8.2 * rSol

    surfaceVel = 250.*1e5

    write(*,'(a)') "Reading Donati maps..."

    write(*,'(a)') "Density..."
    select case (misc)
       case("norm")
          open(20, file = "rho35_closed.map", status = "old", form="formatted")
       case("low")
          open(20, file = "rho35_m-0.7_closed.map", status = "old", form="formatted")
       case DEFAULT
          write(*,'(a)') "! unrecognised misc"
    end select
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,*) nz, ny
    allocate(z(1:nz))
    allocate(y(1:ny))
    allocate(rho(ny, nz))
    allocate(t(ny, nz))
    allocate(v(ny, nz))
    allocate(phi(ny, nz))
    read(20,*) (z(i),i=1,nz)
    read(20,*) (y(i),i=1,ny)
    read(20,*) ((rho(i,j),i=1,nz),j=1,ny)
    close(20)

    write(*,'(a)') "Velocity..."
    select case (misc)
       case("norm")
          open(20, file = "vel35_closed.map", status = "old", form="formatted")
       case("low")
          open(20, file = "vel35_m-0.7_closed.map", status = "old", form="formatted")
       case DEFAULT
          write(*,'(a)') "! unrecognised misc"
    end select
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,*) nz, ny
    read(20,*) (z(i),i=1,nz)
    read(20,*) (y(i),i=1,ny)
    read(20,*) ((v(i,j),i=1,nz),j=1,ny)
    close(20)

    write(*,'(a)') "Phi..."
    select case (misc)
       case("norm")
          open(20, file = "phi35_closed.map", status = "old", form="formatted")
       case("low")
          open(20, file = "phi35_m-0.7_closed.map", status = "old", form="formatted")
       case DEFAULT
          write(*,'(a)') "! unrecognised misc"
    end select

    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,*) nz, ny
    read(20,*) (z(i),i=1,nz)
    read(20,*) (y(i),i=1,ny)
    read(20,*) ((phi(i,j),i=1,nz),j=1,ny)
    close(20)

    write(*,'(a)') "Temperature..."
    select case (misc)
       case("norm")
          open(20, file = "temp35_closed.map", status = "old", form="formatted")
       case("low")
          open(20, file = "temp35_m-0.7_closed.map", status = "old", form="formatted")
       case DEFAULT
          write(*,'(a)') "! unrecognised misc"
    end select
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,*) nz, ny
    read(20,*) (z(i),i=1,nz)
    read(20,*) (y(i),i=1,ny)
    read(20,*) ((t(i,j),i=1,nz),j=1,ny)
    close(20)
    write(*,'(a)') "Done."

    z = 10.**z
    v = 10.**v
    rho = 10.**rho
    t = 10.**t




    y(1) = 0.
    z(1) = 0.


    do i = 1, grid%nx
       grid%xAxis(i) = 2.*real(i-1)/real(grid%nx-1)-1.
    enddo
    do i = 1, grid%ny
       grid%yAxis(i) = 2.*real(i-1)/real(grid%ny-1)-1.
    enddo
    grid%xAxis = grid%xAxis * 3.
    grid%yAxis = grid%yAxis * 3.

    i1 = grid%nz/2 + 1



    grid%zAxis(i1) = 0.
    grid%zAxis(i1+1) = 0.001
    fac = 1.1
    do i = i1+2,grid%nz
       grid%zAxis(i) = grid%zAxis(i-1) * fac
    enddo


    do i = 1,i1-1
       grid%zAxis(i) = -grid%zAxis(grid%nz-i+1)
    enddo

    grid%zAxis = 3.*grid%zAxis/ abs(grid%zAxis(1))

    do i = 1, grid%nz
       write(*,*) i, grid%zAxis(i)
    enddo
       



    grid%rho = 1.e-30
    grid%temperature = 1.e-30

    grid%inUse = .false.

    write(*,'(a)') "Re-mapping grids..."
    do i1 = 1, grid%nx
       do i2 = 1, grid%ny
          do i3 = 1, grid%nz
             r = sqrt(grid%xAxis(i1)**2 + grid%yAxis(i2)**2 + grid%zAxis(i3)**2)
             if ((r > xDisk(1)).and.(r < 3.)) then
                grid%inUse(i1,i2,i3) = .true.
                r = sqrt(grid%xAxis(i1)**2 + grid%yAxis(i2)**2)
                call hunt(y, ny, r, j)
                call hunt(z, nz, abs(grid%zAxis(i3)), i)
                grid%rho(i1, i2, i3) = rho(i,j)
                grid%temperature(i1, i2, i3) = t(i,j)
                ang = atan2(grid%yAxis(i2),grid%xAxis(i1))
                rHat%z = cos(phi(i,j))
                if (grid%zAxis(i3) < 0.) rHat%z = -rHat%z
                fac = sqrt(1.d0 - rHat%z**2)
                rHat%x = fac * cos(ang)
                rHat%y = fac * sin(ang)
                grid%velocity(i1, i2, i3) = (v(i,j)/cSpeed) * rHat
             endif
          enddo
       enddo
    enddo

    where((grid%temperature < 1.).or.(grid%rho < 1.e-19))
       grid%rho = 1.e-30
       grid%inUse = .false.
    endwhere
    write(*,'(a)') "Done."

 
    write(*,'(a)') "Adding disk..."
    do i1 = 1, grid%nx
       do i2 = 1, grid%ny
          r = sqrt(grid%xAxis(i1)**2 + grid%yAxis(i2)**2)
          if ((r >= xDisk(1)).and.(r <= xDisk(nDisk))) then

             select case(misc)
                case("norm")
                   call locate(xDisk, ndisk, r, j)
                   thickness = yDisk(j) + (r - xDisk(j))*(yDisk(j+1)-yDisk(j))/ &
                        (xDisk(j+1)-xDisk(j))
                   diskDensity = rhoDisk(j) + (r - xDisk(j))*(rhoDisk(j+1)-rhoDisk(j))/ &
                        (xDisk(j+1)-xDisk(j))
                case("low")
                   call locate(xDisk, ndisk, r, j)
                   thickness = yDisk2(j) + (r - xDisk(j))*(yDisk2(j+1)-yDisk2(j))/ &
                        (xDisk(j+1)-xDisk(j))
                   diskDensity = rhoDisk2(j) + (r - xDisk(j))*(rhoDisk2(j+1)-rhoDisk2(j))/ &
                        (xDisk(j+1)-xDisk(j))
             end select

             rVec = VECTOR(grid%xAxis(i1),grid%yAxis(i2),-thickness/2.)
             call getIndices(grid, rVec, i, j, j1, t1, t2, t3)
             rVec = VECTOR(grid%xAxis(i1),grid%yAxis(i2), thickness/2.)
             call getIndices(grid, rVec, i, j, j2, t1, t2, t3)
             do j = j1, j2
                grid%inUse(i1,i2,j) = .true.
                grid%rho(i1,i2,j) = (diskDensity*1.e12) * mHydrogen 
                grid%temperature(i1,i2,j) = 40000.
                zPrime = grid%zAxis(j)
                if (zPrime > 0.) then
                   vHat = VECTOR(0.,0.,-1.)
                else
                   vHat = VECTOR(0.,0.,+1.)
                endif
                vel = abs(zPrime) * (2.e3 * 1.e5) ! 2000km/s per rstar
                grid%velocity(i1,i2,j) = (vel/cSpeed) * vHat
             enddo
          endif
       enddo
    enddo
    
    write(*,'(a)') "Done."


    where (grid%temperature > 1.e6)
       grid%temperature = 1.e6
    end where


    grid%xAxis = grid%xAxis * rStar / 1.e10
    grid%yAxis = grid%yAxis * rStar / 1.e10
    grid%zAxis = grid%zAxis * rStar / 1.e10



    grid%rCore = rStar/1.e10

    write(*,'(a,1p,2e12.5)') "Maximum density: ",&
         MAXVAL(grid%rho, mask=grid%inUse)/mHydrogen,MAXVAL(grid%rho, mask=grid%inUse)
    write(*,'(a,1p,2e12.5)') "Minimum density: ", &
         MINVAL(grid%rho, mask=grid%inUse)/mHydrogen,MINVAL(grid%rho, mask=grid%inUse)
    write(*,'(a,1p,2e12.5)') "Maximum/min temp: ",&
         MAXVAL(grid%temperature, mask=grid%inUse),MINVAL(grid%temperature, mask=grid%inUse)

    grid%rstar1 = grid%rCore
    grid%rstar2 = 0.
    grid%starpos1 = VECTOR(0.,0.,0.)
    grid%starpos2 = VECTOR(1.e20,0.,0.)


    grid%etaLine = 1.e-30
    grid%etaCont = 1.e-30
    grid%chiLine = 1.e-30
    grid%kappaSca = 1.e-30
    grid%kappaAbs = 1.e-30

    write(*,'(a)') "Done."


    if (resonanceLine) then

       grid%etaLine = 1.e-30
       grid%etaCont = 1.e-30
       grid%chiLine = 1.e-30
       grid%kappaSca = 1.e-30
       grid%kappaAbs = 1.e-30
       grid%resonanceLine = .true.
       grid%lambda2 = 1548.202

       grid%chiLine = grid%rho*1.e26
       write(*,'(a)') "This is a resonance line model."

    endif



  end subroutine fillGridDonati2


  !  
  !  This subrouine recursively writes out the distance from the center of the
  !  given plane,  with a give value of the 3rd dimension, and the correcsponding 
  !  value stored in a octal to a file specified by the logical unit number. 
  !  
  recursive subroutine radial_profile(thisOctal, name, plane, val_3rd_dim, luout, center, grid)
    use octal_mod, only: subcellCentre
    implicit none
    !
    type(octal), pointer   :: thisOctal
    type(gridtype) :: grid
    character(len=*), intent(in)  :: name     ! "rho", "temperature", chiLine", "etaLine",  
    !                                         ! "etaCont", "Vx", "Vy" or "Vz"
    character(len=*), intent(in)  :: plane    ! must be 'x-y', 'y-z', 'z-x', 'x-z'
    ! The value of the third dimension.
    ! For example, if plane = "x-y", the third dimension is the value of z.
    ! Then, when  val_3rd_dim = 0.0, this will plot the density (and the grid)
    ! on the z=0 plane, and so on...
    real, intent(in)              :: val_3rd_dim 
    integer, intent(in) :: luout  ! file unit number for the output
    type(VECTOR), intent(in) :: center   ! position of the center of the root node
    !
    !
    type(octal), pointer  :: child 
    type(VECTOR) :: rvec
    real :: value

    integer :: subcell, i
    real :: xp, yp, xm, ym, zp, zm
!    real(double) :: kabs, ksca
!    integer :: ilam
    real(double) :: d, L, eps, distance
    logical :: use_this_subcell

  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call radial_profile(child, name, plane, val_3rd_dim, &
                     luout, center, grid)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal,subcell)
          L = thisOctal%subcellSize
          d = L/2.0d0
          eps = d/1000.0d0
          xp = REAL(rVec%x + d)
          xm = REAL(rVec%x - d)
          yp = REAL(rVec%y + d)
          ym = REAL(rVec%y - d)
          zp = REAL(rVec%z + d)
          zm = REAL(rVec%z - d)     
          
          use_this_subcell = .false.
          if ( plane(1:3) == "x-y" .and. (ABS(rVec%z-val_3rd_dim)-d) < eps) then
             use_this_subcell = .true.      
          elseif ( plane(1:3) == "y-z" .and. (ABS(rVec%x-val_3rd_dim)-d) < eps) then 
             use_this_subcell = .true.
          elseif ( plane(1:3) == "z-x" .and. (ABS(rVec%y-val_3rd_dim)-d) < eps) then 
             use_this_subcell = .true.
          elseif ( plane(1:3) == "x-z" .and. (ABS(rVec%y-val_3rd_dim)-d) < eps) then 
             use_this_subcell = .true.
          else
             use_this_subcell = .false.
          end if
          
          if ((plane(1:3)=="x-y").and.(thisOctal%twoD)) then
             use_this_subcell = .false.
             if  ( ((rVec%z-L/2.d0) < 0.d0).and.((rvec%z+L/2.d0) >= 0.)) use_this_subcell = .true.
          endif

          
          if (use_this_subcell) then
             select case (name)
             case("rho")
                value = thisOctal%rho(subcell)
             case("temperature")
                value = thisOctal%temperature(subcell)
             case("chiLine")
                value = thisOctal%chiLine(subcell)
             case("etaLine")
                value = thisOctal%etaLine(subcell)
             case("etaCont")
                value = thisOctal%etaCont(subcell)
             case("Vx")
                value = thisOctal%velocity(subcell)%x * cSpeed/1.0d5 ![km/s]
             case("Vy")
                value = thisOctal%velocity(subcell)%y * cSpeed/1.0d5 ![km/s]
             case("Vz")
                value = thisOctal%velocity(subcell)%z * cSpeed/1.0d5 ![km/s]
!             case("tau")
!                call returnKappa(grid, thisOctal, subcell, ilam, grid%lamArray(ilam), kappaSca=ksca, kappaAbs=kabs)
!                value = thisOctal%subcellsize * (kSca+kAbs)
             case default
                value = 666.
             end select

!             distance = modulus(rvec-center)   ! length of the vector
             distance = modulus(rVec)
             write(luout, '(2(2x, 1PE18.4))')  distance, value


          end if
       end if

    end do

  end subroutine radial_profile

