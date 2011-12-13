module datacube_mod

  use kind_mod
  use interferometry_mod
  use messages_mod

  implicit none

  type TELESCOPE
     
     character(len=10) :: label
     real :: diameter
     real :: beamsize ! in arcsecs
     
  end type TELESCOPE

  type DATACUBE
     TYPE(TELESCOPE) :: telescope

     character(len=80) :: label
     character(len=10) :: vUnit ! units for velocity
     character(len=10) :: xUnit ! units for space
     character(len=10) :: IntensityUnit ! units for intensity

     integer, pointer :: nsubpixels(:,:,:) => null() ! contains resolution information 
     integer, pointer :: converged(:,:,:) => null() ! contains convergence information (should take 1 or 0)
     real(double),pointer :: weight(:,:) => null()    ! Weighting for integration (used to find spectra)
     integer :: nx 
     integer :: ny
     integer :: nv

     real(double) :: obsdistance ! observation distance for use with instrument functions
     real(double), pointer :: xAxis(:)
     real(double), pointer :: yAxis(:)
     real(double), pointer :: vAxis(:)
     character(len=8) :: xaxisType ! To be written as a FITS keyword
     character(len=8) :: yaxisType ! To be written as a FITS keyword
     character(len=8) :: vaxisType ! To be written as a FITS keyword
     real, pointer :: intensity(:,:,:) => null()
     real, pointer :: tau(:,:,:) => null()
     real, pointer :: nCol(:,:) => null()
     real, pointer :: i0_pos(:,:,:) => null() ! Positive contribution to flux
     real, pointer :: i0_neg(:,:,:) => null() ! Negative contribution to flux

  end type DATACUBE

  real, save :: cubePositionAngle

contains

#ifdef USECFITSIO
  subroutine writeDataCube(thisCube, filename, write_Intensity, write_ipos, write_ineg, write_Tau, write_nCol, write_axes)

    use fits_utils_mod

    implicit none
    
    type(DATACUBE), intent(in) :: thisCube
    character(len=*) :: filename

! Optional arguments which allow allocated components of the data cube to not be written.
    logical, optional, intent(in) :: write_Intensity
    logical, optional, intent(in) :: write_ipos
    logical, optional, intent(in) :: write_ineg
    logical, optional, intent(in) :: write_Tau
    logical, optional, intent(in) :: write_nCol
    logical, optional, intent(in) :: write_axes

    integer :: status,unit,blocksize,bitpix,naxis
    integer, dimension(5) :: naxes
    integer :: group,fpixel,nelements
    logical :: simple, extend
    
    logical :: do_write_Intensity
    logical :: do_write_ipos
    logical :: do_write_ineg
    logical :: do_write_Tau
    logical :: do_write_nCol
    logical :: do_write_xaxis
    logical :: do_write_yaxis
    logical :: do_write_vaxis


! Decide which arrays to write. Default is to write all which are allocated.
! Optional arguments override the defaults. 
    if ( .not. associated(thiscube%intensity) ) then 
       do_write_Intensity = .false. 
    else if ( present(write_Intensity) ) then 
       do_write_Intensity = write_Intensity
    else
       do_write_Intensity = .true. 
    end if

    if ( .not. associated(thisCube%tau) ) then 
       do_write_Tau = .false.
    else if ( present(write_Tau) ) then 
       do_write_Tau = write_Tau
    else
       do_write_Tau = .true. 
    end if

    if ( .not. associated(thisCube%nCol) ) then 
       do_write_nCol = .false.
    else if ( present(write_nCol) ) then 
       do_write_nCol = write_nCol
    else
       do_write_nCol = .true. 
    end if

    if ( .not. associated(thisCube%i0_pos) ) then 
       do_write_ipos = .false.
    else if ( present(write_ipos) ) then 
       do_write_ipos = write_ipos
    else
       do_write_ipos = .true. 
    end if

    if ( .not. associated(thisCube%i0_neg) ) then 
       do_write_ineg = .false.
    else if ( present(write_ineg) ) then 
       do_write_ineg = write_ineg
    else
       do_write_ineg = .true. 
    end if

    if ( .not. associated(thisCube%xAxis) ) then 
       do_write_xaxis = .false.
    else
       do_write_xaxis = .true.
    end if

    if ( .not. associated(thisCube%yAxis) ) then 
       do_write_yaxis = .false.
    else
       do_write_yaxis = .true.
    end if

    if ( .not. associated(thisCube%vAxis) ) then 
       do_write_vaxis = .false.
    else
       do_write_vaxis = .true.
    end if

    if ( present(write_axes) ) then 
       if (.not. write_axes ) then 
          do_write_xaxis = .false.
          do_write_yaxis = .false.
          do_write_vaxis = .false.
       end if
    end if

    call checkBitpix(FitsBitpix)

    status=0
    !
    !  Delete the file if it already exists, so we can then recreate it.
    !
    call deleteFitsFile ( trim(filename), status )
    
    !  Get an unused Logical Unit Number to use to open the FITS file.
    call ftgiou ( unit, status )
    !  Create the new empty FITS file.
    blocksize=1
    call ftinit(unit,trim(filename),blocksize,status)

    !  Initialize parameters about the FITS image
    simple=.true.
    extend=.true.
    group=1
    fpixel=1


    ! 1st HDU : flux
    if( do_write_Intensity ) then

       bitpix=FitsBitpix
       naxis=3
       naxes(1)=thisCube%nx
       naxes(2)=thisCube%ny
       naxes(3)=thisCube%nv
       nelements=thisCube%nx*thisCube%ny*thisCube%nv

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

       !  Write keywords to the header.
       call ftpkys(unit,'LABEL',thisCube%label,"Telescope label",status) 
       call ftpkys(unit,'VUNIT',thisCube%vUnit,"Velocity unit",status)
       call ftpkys(unit,'XUNIT',thisCube%xUnit,"Spatial unit",status)
       call ftpkys(unit,'BUNIT',thisCube%IntensityUnit,"Intensity unit",status)
       call ftpkyd(unit,'DISTANCE',thisCube%obsdistance,-3,'observation distance',status)
       call addScalingKeywords(maxval(thisCube%intensity), minval(thisCube%intensity), unit, FitsBitpix)
       call addWCSinfo

       !  Write the array to the FITS file.
       call ftppre(unit,group,fpixel,nelements,thisCube%intensity,status)
    endif

    if( do_write_Tau ) then

       ! 2nd HDU : tau
       call FTCRHD(unit, status)
       bitpix=FitsBitpix
       naxis=3
       naxes(1)=thisCube%nx
       naxes(2)=thisCube%ny
       naxes(3)=thisCube%nv
       nelements=naxes(1)*naxes(2)*naxes(3)

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       call ftpkys(unit,'BUNIT',"Optical depth","Optical depth unit",status)
       call addScalingKeywords(maxval(thisCube%tau), minval(thisCube%tau), unit, FitsBitpix)
       call addWCSinfo

       !  Write the array to the FITS file.
       call ftppre(unit,group,fpixel,nelements,thisCube%tau,status)
    endif

    if(associated(thisCube%weight)) then
       
       ! 3rd HDU : weight
       call FTCRHD(unit, status)
       bitpix=-64
       naxis=2
       naxes(1)=thisCube%nx
       naxes(2)=thisCube%ny
       nelements=naxes(1)*naxes(2)

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       
       !  Write the array to the FITS file.
       call ftpprd(unit,group,fpixel,nelements,thisCube%weight,status)
    endif

    if( do_write_xaxis ) then
       ! 4th HDU : xAxis
       call FTCRHD(unit, status)
       bitpix=-64
       naxis=1
       naxes(1)=thisCube%nx     
       nelements=naxes(1)

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       
       !  Write the array to the FITS file.
       call ftpprd(unit,group,fpixel,nelements,thisCube%xAxis,status)
    endif

    if( do_write_yaxis ) then
       ! 5th HDU : yAxis
       call FTCRHD(unit, status)
       bitpix=-64
       naxis=1
       naxes(1)=thisCube%ny    
       nelements=naxes(1)

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       
       !  Write the array to the FITS file.
       call ftpprd(unit,group,fpixel,nelements,thisCube%yAxis,status)
    endif

    if( do_write_vaxis) then
       ! 6th HDU : vAxis
       call FTCRHD(unit, status)
       bitpix=-64 
       naxis=1
       naxes(1)=thisCube%nv        
       nelements=naxes(1)
       
       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       
       !  Write the array to the FITS file.
       call ftpprd(unit,group,fpixel,nelements,thisCube%vAxis,status)
    endif
       
    if(associated(thisCube%nSubpixels)) then
       ! 7th HDU : nsubpixels
       call FTCRHD(unit, status)
       bitpix=32
       naxis=3
       naxes(1)=thisCube%nx
       naxes(2)=thisCube%ny
       naxes(3)=thisCube%nv
       nelements=naxes(1)*naxes(2)*naxes(3)
       
       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       
       !  Write the array to the FITS file.
       call ftpprj(unit,group,fpixel,nelements,thisCube%nsubpixels,status)
    endif

    if( do_write_nCol ) then
       ! 8th HDU : nCol
       call FTCRHD(unit, status)
       bitpix=FitsBitpix
       naxis=2
       naxes(1)=thisCube%nx
       naxes(2)=thisCube%ny
       nelements=naxes(1)*naxes(2)
       
       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       call ftpkys(unit,'BUNIT',"cm^-2","Column density unit",status)
       call addScalingKeywords(maxval(thisCube%nCol), minval(thisCube%nCol), unit, FitsBitpix)
       call addWCSinfo

       !  Write the array to the FITS file.
       call ftppre(unit,group,fpixel,nelements,thisCube%nCol,status)
       call printFitsError(status)
    endif

    if( do_write_ipos) then

       ! 9th HDU : positive intensity contribution
       call FTCRHD(unit, status)
       bitpix=FitsBitpix
       naxis=3
       naxes(1)=thisCube%nx
       naxes(2)=thisCube%ny
       naxes(3)=thisCube%nv
       nelements=naxes(1)*naxes(2)*naxes(3)

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       call ftpkys(unit,'BUNIT',thisCube%IntensityUnit,"Intensity unit",status)
       call addScalingKeywords(maxval(thisCube%i0_pos), minval(thisCube%i0_pos), unit, FitsBitpix)
       call addWCSinfo

       !  Write the array to the FITS file.
       call ftppre(unit,group,fpixel,nelements,thisCube%i0_pos,status)
    endif


    if( do_write_ineg) then

       ! 10th HDU : negative intensity contribution
       call FTCRHD(unit, status)
       bitpix=FitsBitpix
       naxis=3
       naxes(1)=thisCube%nx
       naxes(2)=thisCube%ny
       naxes(3)=thisCube%nv
       nelements=naxes(1)*naxes(2)*naxes(3)

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       call ftpkys(unit,'BUNIT',thisCube%IntensityUnit,"Intensity unit",status)
       call addScalingKeywords(maxval(thisCube%i0_neg), minval(thisCube%i0_neg), unit, FitsBitpix)
       call addWCSinfo

       !  Write the array to the FITS file.
       call ftppre(unit,group,fpixel,nelements,thisCube%i0_neg,status)
    endif

    !  Close the file and free the unit number.
    call ftclos(unit, status)
    call ftfiou(unit, status)

     !  Check for any error, and if so print out error messages
    if (status > 0) then
       call printFitsError(status)
    end if
    

    contains

      subroutine addWCSinfo
        implicit none

        ! Write WCS keywords to the header
        call ftpkyd(unit,'CRPIX1',0.5_db,-3,'reference pixel',status)
        call ftpkyd(unit,'CDELT1',thisCube%xAxis(2)-thisCube%xAxis(1),-3,'coordinate increment at reference point',status)
        call ftpkys(unit,'CTYPE1',thisCube%xAxisType,"x axis", status)
        call ftpkyd(unit,'CRVAL1',thisCube%xAxis(1),-3,'coordinate value at reference point',status)
        call ftpkys(unit,'CUNIT1', thisCube%xUnit, "x axis unit", status)

        call ftpkyd(unit,'CRPIX2',0.5_db,-3,'reference pixel',status)
        call ftpkyd(unit,'CDELT2',thisCube%yAxis(2)-thisCube%yAxis(1),-3,'coordinate increment at reference point',status)
        call ftpkys(unit,'CTYPE2',thisCube%yAxisType, "y axis", status)
        call ftpkyd(unit,'CRVAL2',thisCube%yAxis(1),-3,'coordinate value at reference point',status)
        call ftpkys(unit,'CUNIT2', thisCube%xUnit, "y axis unit", status)

        call ftpkyd(unit,'CRPIX3',0.5_db,-3,'reference pixel',status)
        if (SIZE(thisCube%vAxis)  > 1) then
           call ftpkyd(unit,'CDELT3',thisCube%vAxis(2)-thisCube%vAxis(1),-3,'coordinate increment at reference point',status)
        endif
        call ftpkys(unit,'CTYPE3',thisCube%vAxisType, "velocity axis", status)
        call ftpkyd(unit,'CRVAL3',thisCube%vAxis(1),-3,'coordinate value at reference point',status)

      end subroutine addWCSinfo
    
    end subroutine writeDataCube

#endif

  !***********************************************************

! Initialises cube - sets intensity for cube to 0 
  subroutine initCube(thisCube, nx, ny, nv, mytelescope, splitCubes, wantTau, galacticPlaneSurvey)

    type(DATACUBE) :: thisCube
    type(TELESCOPE), optional :: mytelescope
    integer :: nx, ny, nv
    logical, optional, intent(in) :: splitCubes, wantTau, galacticPlaneSurvey

    if(present(mytelescope)) then
       
       thisCube%telescope = mytelescope

    else

       thisCube%telescope%label=" "
       thisCube%telescope%diameter = 1.
       thisCube%telescope%beamsize = 1d10 !arcsecs - set large to give unconvolved image

    endif

    thisCube%label        = " "
    thisCube%vaxisType    = "velocity"
    thisCube%vUnit        = "km/s"
    thisCube%xAxisType    = "X"
    thisCube%yAxisType    = "Y"
    thisCube%xUnit        = "10^10cm "
    thisCube%IntensityUnit= " "

! Set up cube for a Galactic plane survey, as generated by angularImage_mod.
    if ( present(galacticPlaneSurvey) ) then 
       if (galacticPlaneSurvey) then 
          thisCube%xUnit     = "degrees"
          thisCube%xAxisType = "GLON-CAR"
          thisCube%yAxisType = "GLAT-CAR"
          thisCube%vAxisType = "VELO-LSR"
          thisCube%intensityUnit = "K (Tb)  "
       end if
    end if

    thisCube%nx = nx
    thisCube%ny = ny
    thisCube%nv = nv
    allocate(thisCube%xAxis(1:nx))
    allocate(thisCube%yAxis(1:ny))
    allocate(thisCube%vAxis(1:nv))
    allocate(thisCube%intensity(1:nx,1:ny,1:nv))

    allocate(thisCube%nCol(1:nx,1:ny))
!    allocate(thisCube%nsubpixels(1:nx,1:ny,1:nv))
!    allocate(thisCube%converged(1:nx,1:ny,1:nv))
!    allocate(thisCube%weight(1:nx,1:ny))

    thisCube%intensity = 0.d0
    thisCube%nCol = 0.d0

!    thisCube%nsubpixels = 0.d0
!    thisCube%converged = 0
!    thisCube%weight = 1.d0

    if ( present(wantTau) ) then 
       if (wantTau) then
          allocate(thisCube%tau(1:nx,1:ny,1:nv))
          thisCube%tau = 0.d0
       end if
    end if

    if ( present(splitCubes) ) then 
       if (splitCubes) then 
          allocate(thisCube%i0_pos(1:nx,1:ny,1:nv))
          allocate(thisCube%i0_neg(1:nx,1:ny,1:nv))
          thisCube%i0_pos(:,:,:) = 0.0
          thisCube%i0_neg(:,:,:) = 0.0
       end if
    end if

  end subroutine initCube

! Set spatial axes for datacube - Equally spaced (linearly) between min and max
  subroutine addSpatialAxes(cube, xMin, xMax, yMin, yMax, griddistance)
    use constants_mod, only: auTocm, pi, rSol
    type(DATACUBE) :: cube
    real(double) :: xMin, xMax, yMax, yMin, dx, dy
    integer :: i
    character(len=100) :: message
    real, intent(in) :: gridDistance

    dx = (xMax - xMin)/dble(cube%nx)
    dy = (yMax - yMin)/dble(cube%ny)

    write(message,'(a,1pe12.3,a)') "Linear pixel resolution  : ", dx*1e10/autocm, " AU"
    call writeinfo(message,TRIVIAL)
    write(message,'(a,1pe12.3,a)') "Linear pixel resolution  : ", dx*1e10/rSol, " Rsol"
    call writeinfo(message,TRIVIAL)
    write(message,*) "Angular pixel resolution : ", (dx*1e10/griddistance)*(180./pi)*60.*60., " arcseconds"
    call writeinfo(message,TRIVIAL)

    do i = 1, cube%nx
       cube%xAxis(i) = xmin + dble(i-0.5d0)*dx 
    enddo

    do i = 1, cube%ny
       cube%yAxis(i) = ymin + dble(i-0.5d0)*dx 
    enddo

  end subroutine addSpatialAxes

  subroutine convertSpatialAxes(cube, newUnit)

! Purpose: Convert spatial axis units from torus length units to something else
! Author: D. Acreman

    use constants_mod, only: autocm, pctocm
    type(DATACUBE) :: cube
    character(len=*), intent(in) :: newUnit

    select case (newUnit) 

    case ('pc')
       cube%xAxis(:) = cube%xAxis(:) * 1.0e10_db / pctocm
       cube%yAxis(:) = cube%yAxis(:) * 1.0e10_db / pctocm
       cube%xUnit    = "pc"

    case ('kpc')
       cube%xAxis(:) = cube%xAxis(:) * 1.0e10_db / (1000.0_db * pctocm)
       cube%yAxis(:) = cube%yAxis(:) * 1.0e10_db / (1000.0_db * pctocm)
       cube%xUnit    = "kpc"

    case('au')
       cube%xAxis(:) = cube%xAxis(:) * 1.0e10_db / autocm
       cube%yAxis(:) = cube%yAxis(:) * 1.0e10_db / autocm
       cube%xUnit    = "AU"

    case('cm')
       cube%xAxis(:) = cube%xAxis(:) * 1.0e10_db
       cube%yAxis(:) = cube%yAxis(:) * 1.0e10_db
       cube%xUnit    = "cm"

    case('m')
       cube%xAxis(:) = cube%xAxis(:) * 1.0e10_db / 100.0
       cube%yAxis(:) = cube%yAxis(:) * 1.0e10_db / 100.0 
       cube%xUnit    = "m"

    case DEFAULT 
       call writewarning('Unrecognised unit in convertSpatialAxes')

    end select

  end subroutine convertSpatialAxes

  subroutine convertVelocityAxis(cube, newUnit)

! Purpose: Convert velocty axis units from km/s to something else
! Author: D. Acreman

    type(DATACUBE) :: cube
    character(len=*), intent(in) :: newUnit

    select case (newUnit) 

    case ('m/s')
       cube%vAxis(:) = cube%vAxis(:) * 1000.0 
       cube%vUnit    = "m/s      "

    case DEFAULT 
       call writewarning('Unrecognised unit in convertVelocity')

    end select

  end subroutine convertVelocityAxis

! Set velocity axis for datacube - Equally spaced (linearly) between min and max
  subroutine addVelocityAxis(cube, vMin, vMax)
    type(DATACUBE) :: cube
    real(double):: vMin, vMax, dv
    integer :: i
    character(len=100) :: message

    if (cube%nv > 1) then
       dv = (vMax - vMin) / dble(cube%nv-1)

       write(message, '(a,f9.4,a)') "Velocity pixel resolution: ", dv, " km/s"
       call writeinfo(message,TRIVIAL)
       
       do i = 1, cube%nv
          cube%vAxis(i) = vmin + dv * real(i-1)
       enddo
    else
       cube%vAxis(1) = vmin
    endif
  end subroutine addVelocityAxis

subroutine TranslateCubeIntensity(cube,constant)

  type(DATACUBE) :: cube
  real(double) :: constant
  
  cube%intensity = cube%intensity + real(constant)
    
end subroutine TranslateCubeIntensity

subroutine freeDataCube(thiscube)
  type(DATACUBE) :: thiscube

    if (associated(thisCube%nsubpixels)) then
       deallocate(thiscube%nSubpixels)
       nullify(thiscube%nSubpixels)
    end if

    if (associated(thisCube%converged)) then 
       deallocate(thiscube%converged)
       nullify(thiscube%converged)
    end if

    if (associated(thisCube%weight)) then 
       deallocate(thiscube%weight)
       nullify(thiscube%weight)
    end if

    if (associated(thisCube%xAxis)) then 
       deallocate(thiscube%xAxis)
       nullify(thiscube%xAxis)
    end if

    if (associated(thisCube%yAxis)) then 
       deallocate(thiscube%yAxis)
       nullify(thiscube%yAxis)
    end if

    if (associated(thisCube%vAxis)) then 
       deallocate(thiscube%vAxis)
       nullify(thiscube%vAxis)
    end if

    if (associated(thisCube%intensity)) then 
       deallocate(thiscube%intensity)
       nullify(thiscube%intensity)
    end if

    if (associated(thisCube%tau)) then 
       deallocate(thiscube%tau)
       nullify(thiscube%tau)
    end if

    if (associated(thisCube%nCol)) then 
       deallocate(thiscube%nCol)
       nullify(thiscube%nCol)
    end if

    if (associated(thisCube%i0_pos)) then 
       deallocate(thiscube%i0_pos)
       nullify(thiscube%i0_pos)
    end if

    if (associated(thisCube%i0_neg)) then 
       deallocate(thiscube%i0_neg)
       nullify(thiscube%i0_neg)
    end if

  end subroutine freeDataCube
  
#ifdef USECFITSIO
  subroutine writeCollapsedDataCube(thisCube, filename)
    use fits_utils_mod

    type(DATACUBE) :: thisCube
    character(len=*) :: filename

    integer :: status,unit,blocksize,bitpix,naxis,naxes(2)
    integer :: group,fpixel,nelements
    real, allocatable :: array(:,:)
    logical :: simple,extend
    integer :: i, j

    allocate(array(1:thisCube%nx, 1:thisCube%ny))
    call writeInfo("Writing fits image",TRIVIAL)
    status=0

    do i = 1, thisCube%nx
       do j = 1, thisCube%ny
          array(i,j) = SUM(thisCube%intensity(i,j,:))
       enddo
    enddo

    !
    !  Delete the file if it already exists, so we can then recreate it.
    !
    call deleteFitsFile ( filename, status )
    !
    !  Get an unused Logical Unit Number to use to open the FITS file.
    !
    call ftgiou ( unit, status )
    !
    !  Create the new empty FITS file.
    !
    blocksize=1
    call ftinit(unit,filename,blocksize,status)
    !
    !  Initialize parameters about the FITS image (300 x 200 16-bit integers).
    !
    simple=.true.
    bitpix=FitsBitpix
    naxis=2
    naxes(1)=thisCube%nx
    naxes(2)=thisCube%ny
    extend=.true.
    !
    !  Write the required header keywords.
    !
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
    !
    !  Write the array to the FITS file.
    !
    group=1
    fpixel=1
    nelements=naxes(1)*naxes(2)

    call ftppre(unit,group,fpixel,nelements,array,status)
    ! Write WCS keywords to the header
    call ftpkyd(unit,'CRPIX1',0.5_db,-3,'reference pixel',status)
    call ftpkyd(unit,'CDELT1',thisCube%xAxis(2)-thisCube%xAxis(1),-3,'coordinate increment at reference point',status)
    call ftpkys(unit,'CTYPE1',thisCube%xAxisType,"x axis", status)
    call ftpkyd(unit,'CRVAL1',thisCube%xAxis(1),-3,'coordinate value at reference point',status)
    call ftpkys(unit,'CUNIT1', thisCube%xUnit, "x axis unit", status)
    
    call ftpkyd(unit,'CRPIX2',0.5_db,-3,'reference pixel',status)
    call ftpkyd(unit,'CDELT2',thisCube%yAxis(2)-thisCube%yAxis(1),-3,'coordinate increment at reference point',status)
    call ftpkys(unit,'CTYPE2',thisCube%yAxisType, "y axis", status)
    call ftpkyd(unit,'CRVAL2',thisCube%yAxis(1),-3,'coordinate value at reference point',status)
    call ftpkys(unit,'CUNIT2', thisCube%xUnit, "y axis unit", status)

    !
    !  Close the file and free the unit number.
    !
    call ftclos(unit, status)
    call ftfiou(unit, status)
    !
    !  Check for any error, and if so print out error messages
    !
    if (status > 0) then
       call printFitserror(status)
    end if
    deallocate(array)

  end subroutine writeCollapsedDataCube
#endif


  subroutine dumpCubeToVisibilityCurves(thisCube, visFile, wavelength, gridDistance)

    type(DATACUBE) :: thisCube
    character(len=*) :: visFile
    character(len=80) :: tempFile
    real(double) :: wavelength
    real(double), allocatable :: dimage(:,:), xAxis(:), yAxis(:)
    integer :: iv
    real, intent(in) :: gridDistance 

    allocate(dImage(1:thisCube%nx,1:thisCube%ny))
    allocate(xAxis(1:thisCube%nx))
    allocate(yAxis(1:thisCube%ny))
    do iv = 1, thiscube%nv
       write(tempFile,'(a, a, SP, i5.4, a)') trim(visFile),"_v=",nint(thiscube%vAxis(iv)), ".dat"
       dImage(1:thisCube%nx,1:thisCube%ny) = &
            reshape(thisCube%intensity(1:thisCube%nx,1:thisCube%ny,iv:iv),(/thisCube%nx,thisCube%ny/))
       xAxis = thisCube%xAxis * 1.d10 / gridDistance
       yAxis = thisCube%yAxis * 1.d10 / gridDistance
       call visibilityCurve(tempFile, dImage, &
            thisCube%nx, thisCube%ny, xAxis, yAxis, wavelength, 300.d2)
    enddo
    deallocate(dimage,xaxis, yaxis)
  end subroutine dumpCubeToVisibilityCurves

  subroutine dumpCubeToSpectrum(thisCube, specFile)
    type(DATACUBE) :: thisCube
    character(len=*) :: specFile

    real(double), allocatable :: v(:),y(:)
    integer :: i


    allocate(v(1:thisCube%nv))
    allocate(y(1:thisCube%nv))
    do i = 1, thisCube%nv
       v(i) = thisCube%vAxis(i)
       y(i) = SUM(thisCube%intensity(:,:,i))
    enddo
    y = y / y (1)
    open(38, file=specfile, status="unknown",form="formatted")
    do i = 1, thisCube%nv
       write(38,*) v(i), y(i)
    enddo
    deallocate(v, y)
    close(38)
  end subroutine dumpCubeToSpectrum
       

end module datacube_mod

