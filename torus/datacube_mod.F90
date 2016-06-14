#ifdef FITSCUBE
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
! Weighting for integration (used to find spectra) (molecular_mod::GaussinWeighting) or
! Density weighted temperature (angularImage_mod)
     real(double),pointer :: weight(:,:) => null()    
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
  character(len=80), save :: datacubeFilename
  integer, save          :: npixels  ! 
  integer, save, private :: npixelsX ! number of spatial pixels in x-axis
  integer, save, private :: npixelsY ! number of spatial pixels in y-axis
  real, save, private    :: axisRatio ! ratio of npixelsX to npixelsY
  real(double), save :: RA, DEC
  logical, save, private :: useFixedBg
  real, save, private    :: fixedBg

contains

  subroutine setCubeParams(npix_in, aspectratio, WV_background)
    implicit none 
    integer, intent(in) :: npix_in
    real, intent(in)    :: aspectRatio 
    real, intent(in)    :: WV_background

    npixelsX  = int(npix_in * aspectratio)
    npixelsY  = npixels
    axisRatio = aspectratio

! Set background to subtract when generating moment maps
! Negative values of WV_background set the bg to the minimum value in the data cube
    if (WV_background < 0.0 ) then 
       useFixedBg=.false.
    else
       useFixedBg=.true.
       fixedBg=WV_background
    endif

  end subroutine setCubeParams

#ifdef USECFITSIO
  subroutine writeDataCube(thisCube, filename, write_Intensity, write_ipos, write_ineg, write_Tau, &
       write_Weight, write_nCol, write_axes, write_WV, frequency)
    use inputs_mod, only : ALMA
    use fits_utils_mod
    implicit none
    
    type(DATACUBE), intent(in) :: thisCube
    character(len=*) :: filename

! Optional arguments which allow allocated components of the data cube to not be written.
    logical, optional, intent(in) :: write_Intensity
    logical, optional, intent(in) :: write_ipos
    logical, optional, intent(in) :: write_ineg
    logical, optional, intent(in) :: write_Tau
    logical, optional, intent(in) :: write_Weight
    logical, optional, intent(in) :: write_nCol
    logical, optional, intent(in) :: write_axes
    logical, optional, intent(in) :: write_WV
    real(double), optional, intent(in) :: frequency

    integer :: status,unit,blocksize,bitpix,naxis
    integer, dimension(5) :: naxes
    integer :: group,fpixel,nelements
    logical :: simple, extend
    
    logical :: do_write_Intensity
    logical :: do_write_ipos
    logical :: do_write_ineg
    logical :: do_write_Tau
    logical :: do_write_Weight
    logical :: do_write_nCol
    logical :: do_write_xaxis
    logical :: do_write_yaxis
    logical :: do_write_vaxis
    logical :: do_write_WV


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

    if ( .not. associated(thisCube%weight) ) then 
       do_write_Weight = .false.
    else if ( present(write_Weight) ) then 
       do_write_Weight = write_Weight
    else
       do_write_Weight = .true. 
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

    if ( present(write_WV) ) then 
       do_write_WV = write_WV
    else
       do_write_WV = .false.
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


    if(present(frequency) .and. ALMA) then
       call convertVelocityToHz(thiscube, frequency)
    endif
    

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

    if( do_write_Weight ) then
       
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
       
    if( do_write_nCol ) then
       ! 7th HDU : nCol
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

       ! 8th HDU : positive intensity contribution
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

       ! 9th HDU : negative intensity contribution
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

    if( do_write_WV ) call writeWeightedVelocity

    !  Close the file and free the unit number.
    call ftclos(unit, status)
    call ftfiou(unit, status)

     !  Check for any error, and if so print out error messages
    if (status > 0) then
       call printFitsError(status)
    end if
    

    contains

      subroutine addWCSinfo
        use inputs_mod, only : ALMA
        implicit none
        real(double) :: refPix, refVal, deltaPix, startVal
        real(double) :: refValX, refValY, dx, dy
! 
! Axis 1:
!
        ! Write WCS keywords to the header
        if(ALMA) then
           dx = thisCube%xAxis(2)-thisCube%xAxis(1)
           dy = thisCube%yAxis(2)-thisCube%yAxis(1)
           dx = ((dx * 1.d20)/thisCube%obsdistance)*radtodeg
           dy = ((dy * 1.d20)/thisCube%obsdistance)*radtodeg
           refValX = DEC
           refValY = RA
           call ftpkyd(unit,'CRPIX1',0.5_db,-3,'reference pixel',status)
           call ftpkyd(unit,'CDELT1',dx,10,' ',status)
           call ftpkys(unit,'CTYPE1',thisCube%xAxisType,"x axis", status)
           call ftpkyd(unit,'CRVAL1',refValX,-5,'coordinate value at reference point',status)
           call ftpkys(unit,'CUNIT1', "deg", "x axis unit", status)
           
        else
           call ftpkyd(unit,'CRPIX1',0.5_db,-3,'reference pixel',status)
           call ftpkyd(unit,'CDELT1',thisCube%xAxis(2)-thisCube%xAxis(1),-3,'coordinate increment at reference point',status)
           call ftpkys(unit,'CTYPE1',thisCube%xAxisType,"x axis", status)
           call ftpkyd(unit,'CRVAL1',thisCube%xAxis(1),-3,'coordinate value at reference point',status)
           call ftpkys(unit,'CUNIT1', thisCube%xUnit, "x axis unit", status)
        endif

! 
! Axis 2:
!

        if(ALMA) then
           dx = thisCube%xAxis(2)-thisCube%xAxis(1)
           dy = thisCube%yAxis(2)-thisCube%yAxis(1)
           dx = ((dx * 1.d20)/thisCube%obsdistance)*radtodeg
           dy = ((dy * 1.d20)/thisCube%obsdistance)*radtodeg
           refValX = DEC
           refValY = RA
           call ftpkyd(unit,'CRPIX2',0.5_db,-3,'reference pixel',status)
           call ftpkyd(unit,'CDELT2',dy,10 ,' ',status)
           call ftpkys(unit,'CTYPE2',"DEC--SIN","y axis", status)
           call ftpkyd(unit,'CRVAL2',refValY,-5,'coordinate value at reference point',status)
           call ftpkys(unit,'CUNIT2', "deg", "y axis unit", status)

           call ftpkyd(unit,'CD1_1',dx,10,' ',status)
           call ftpkyd(unit,'CD2_2',dy,10,' ',status)

        else
           ! When dealing with Galactic co-ordinates some software assumes that the reference pixel is at b=0 so we'll
           ! handle Galactic co-ordinates as a special case
           if (thisCube%yAxisType(1:8)=="GLAT-CAR") then
              ! Pixel size assuming uniform pixels
              deltaPix  = thisCube%yAxis(2)-thisCube%yAxis(1)
              ! startVal is the lower edge of the grid which is offset by 1/2 pixel from centre of first pixel
              startVal  = thisCube%yAxis(1) - deltaPix*0.5
              ! -1 x distance from start of the grid to zero terms of in number of pixels
              refPix    = -1.0*startVal / (deltaPix)
              ! Reference latitude is zero by definition
              refVal    = 0.0_db                              
           else
              deltaPix  = thisCube%yAxis(2)-thisCube%yAxis(1)
              refPix    = 0.5_db
              refVal    = thisCube%yAxis(1)
           endif

           call ftpkyd(unit,'CRPIX2',refPix,-3,'reference pixel',status)
           call ftpkyd(unit,'CDELT2',deltaPix,-3,'coordinate increment at reference point',status)
           call ftpkys(unit,'CTYPE2',thisCube%yAxisType, "y axis", status)
           call ftpkyd(unit,'CRVAL2',refVal,-3,'coordinate value at reference point',status)
           call ftpkys(unit,'CUNIT2',thisCube%xUnit, "y axis unit", status)
        endif

! 
! Axis 3:
!
        if(ALMA) then

           call ftpkyd(unit,'CRPIX3',0.5_db,-3,'reference pixel',status)
           if (SIZE(thisCube%vAxis)  > 1) then
              call ftpkyd(unit,'CDELT3',thisCube%vAxis(2)-thisCube%vAxis(1),-3,'coordinate increment at reference point',status)
              call ftpkyd(unit,'CD3_3',thisCube%vAxis(2)-thisCube%vAxis(1),-3,'coordinate increment at reference point',status)
           endif
           
           
           call ftpkys(unit,'CTYPE3',thisCube%vAxisType, "velocity axis", status)
           call ftpkyd(unit,'CRVAL3',thisCube%vAxis(1),-9,'coordinate value at reference point',status)
           call ftpkys(unit,'CUNIT3', "Hz", "vel axis unit", status)           

        else
           call ftpkyd(unit,'CRPIX3',0.5_db,-3,'reference pixel',status)
           if (SIZE(thisCube%vAxis)  > 1) then
              call ftpkyd(unit,'CDELT3',thisCube%vAxis(2)-thisCube%vAxis(1),-3,'coordinate increment at reference point',status)
           endif
           
           call ftpkys(unit,'CTYPE3',thisCube%vAxisType, "velocity axis", status)
           call ftpkyd(unit,'CRVAL3',thisCube%vAxis(1),-3,'coordinate value at reference point',status)
        endif
      end subroutine addWCSinfo

      subroutine writeWeightedVelocity

        real, allocatable :: zeroMoment(:,:), firstMoment(:,:), secondMoment(:,:)
        real, allocatable :: S(:) ! background subtracted intensity
        real, allocatable :: vAxis_sp(:)
        real :: intensitySum, background, deltaV
        integer :: i, j
        character(len=80) :: message

        real, parameter :: blankVal=-1.0e33

        allocate ( zeroMoment(thisCube%nx, thisCube%ny) )
        allocate ( firstMoment(thisCube%nx, thisCube%ny) )
        allocate ( secondMoment(thisCube%nx, thisCube%ny) )
        allocate ( S(thisCube%nv) )
        allocate ( vAxis_sp(thisCube%nv) )
        vAxis_sp(:) = real(thisCube%vAxis(:))
        deltaV      = vAxis_sp(2) - vAxis_sp(1)

        if (useFixedBg) then 
           background = fixedBg
           write(message,*) "Using fixed background of: ", background
           call writeInfo(message,FORINFO)
        else
           background = 0.0
           write(message,*) "Taking background as zero"
           call writeInfo(message,FORINFO)
        endif

        ! Calculate emission weighted velocity
        do j=1, thisCube%ny
           do i=1, thisCube%nx

              S = thisCube%intensity(i,j,:) - background
              where (S<0.0) S=0.0
              intensitySum = sum( S(:) )

              if (intensitySum /= 0.0 ) then 
                 ! Multiply intensity sum by channel width to get zero moment in K.km/s
                 zeroMoment(i,j)   = intensitySum * deltaV
                 firstMoment(i,j)  = sum( S(:)*vAxis_sp(:) ) / intensitySum
                 secondMoment(i,j) = sum (S(:) * ( (vAxis_sp(:)-firstMoment(i,j))**2 ) )
                 secondMoment(i,j) = sqrt(secondMoment(i,j) / intensitySum)
              else
                 zeroMoment(i,j)   = blankVal
                 firstMoment(i,j)  = blankVal
                 secondMoment(i,j) = blankVal
              endif

           end do
        end do


        bitpix=FitsBitpix
        naxis=2
        naxes(1)=thisCube%nx
        naxes(2)=thisCube%ny
        nelements=naxes(1)*naxes(2)

! Zero moment 
        call FTCRHD(unit, status)
        
        ! Write the required header keywords.
        call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
        call ftpkys(unit,'BUNIT',"K km/s","Velocity unit",status)

        call addScalingKeywords(maxval(zeroMoment), minval(zeroMoment), unit, FitsBitpix)
        call addWCSinfo

       !  Write the array to the FITS file.
        call ftppne(unit,group,fpixel,nelements,zeroMoment,blankVal,status)
        call printFitsError(status)

! First moment 
        call FTCRHD(unit, status)
        
        ! Write the required header keywords.
        call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
        call ftpkys(unit,'BUNIT',"km/s","Velocity unit",status)

        call addScalingKeywords(maxval(firstMoment), minval(firstMoment), unit, FitsBitpix)
        call addWCSinfo

       !  Write the array to the FITS file.
        call ftppne(unit,group,fpixel,nelements,firstMoment,blankVal,status)
        call printFitsError(status)

! Second moment
        call FTCRHD(unit, status)
        
        ! Write the required header keywords.
        call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
        call ftpkys(unit,'BUNIT',"km/s","Velocity unit",status)

        call addScalingKeywords(maxval(secondMoment), minval(secondMoment), unit, FitsBitpix)
        call addWCSinfo

       !  Write the array to the FITS file.
        call ftppne(unit,group,fpixel,nelements,secondMoment,blankVal,status)
        call printFitsError(status)

        deallocate (zeroMoment)
        deallocate (firstMoment)
        deallocate (secondMoment)
        deallocate (S)
        deallocate (vAxis_sp)

      end subroutine writeWeightedVelocity
    
    end subroutine writeDataCube

#endif

  !***********************************************************

! Initialises cube - sets intensity for cube to 0 
!
! thisCube%weight is only allocated when using angularImage_mod. If you want to re-instate 
! the GaussianWeighting subroutine in molecular_mod then you'll need to allocate 
! thisCube%weight here. D. Acreman, March 2014.
!
  subroutine initCube(thisCube, nv, mytelescope, splitCubes, wantTau, galacticPlaneSurvey)
    use inputs_mod, only : ALMA
    type(DATACUBE) :: thisCube
    type(TELESCOPE), optional :: mytelescope
    integer :: nx, ny, nv
    logical, optional, intent(in) :: splitCubes, wantTau, galacticPlaneSurvey

    nx = npixelsX
    ny = npixelsY

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
! Weighted temperature is only required when nv=1 i.e. column density mode
          if ( nv==1 ) then 
             allocate(thisCube%weight(1:nx,1:ny))
             thisCube%weight(:,:) = 0.0 
          end if
       endif
    end if



    if(ALMA) then
       thisCube%xUnit     = "degrees "
       thisCube%xAxisType = "RA---SIN"
       thisCube%yAxisType = "DEC--SIN"
       thisCube%vAxisType = "FREQ    "
       thisCube%vUnit =     "Hz      "
       thisCube%intensityUnit =  "Jy/pixel"
       if ( nv==1 ) then          
          allocate(thisCube%weight(1:nx,1:ny))
          thisCube%weight(:,:) = 0.0
       end if
    endif

    thisCube%nx = nx
    thisCube%ny = ny
    thisCube%nv = nv
    allocate(thisCube%xAxis(1:nx))
    allocate(thisCube%yAxis(1:ny))
    allocate(thisCube%vAxis(1:nv))
    allocate(thisCube%intensity(1:nx,1:ny,1:nv))
    allocate(thisCube%nCol(1:nx,1:ny))

    thisCube%intensity = 0.d0
    thisCube%nCol = 0.d0

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
! The axis ratio is taken into account so only need to pass in the y-axis min/max
  subroutine addSpatialAxes(cube, yMin, yMax, griddistance, smallestCell)
    use constants_mod
    type(DATACUBE) :: cube
    real(double), intent(in) :: yMax, yMin
    real(double) :: xMin, xMax, dx, dy, dxInCm
    integer :: i
    character(len=100) :: message
    real, intent(in) :: gridDistance
    real(double), intent(in), optional :: smallestCell ! In Torus units

! Calculate xMin and xMax from the axis ratio 
    xMin = yMin * axisRatio
    xMax = yMax * axisRatio

    dx = (xMax - xMin)/dble(cube%nx)
    dy = (yMax - yMin)/dble(cube%ny)
    dxInCm = dx*1e10

! Smallest grid cell (in AU) is reported by createimage, so report pixel size in AU for comparison.
    write(message,'(a,1pe12.3,a)') "Linear pixel resolution  : ", dxInCm/autocm, " AU"
    call writeinfo(message,TRIVIAL)

! If the smallestCell argument is present then compare with the pixel size to check we aren't doing something daft
    if (present(smallestCell)) then 
       if (dx < 0.1*smallestCell) call writeWarning("X-axis pixel size is more than 10x smaller than the smallest grid cell")
       if (dy < 0.1*smallestCell) call writeWarning("Y-axis pixel size is more than 10x smaller than the smallest grid cell")
    endif

! Report the linear pixel size in suitable units
    if (dxInCm < pcToCm) then 
       write(message,'(a,1pe12.3,a)') "Linear pixel resolution  : ", dxInCm/rSol, " Rsol"
       call writeinfo(message,TRIVIAL)
    else
       write(message,'(a,1pe12.3,a)') "Linear pixel resolution  : ", dxInCm/pcToCm, " pc"
       call writeinfo(message,TRIVIAL)
    endif

    if (gridDistance > 0.0 ) then 
       write(message,'(a,1pe12.3,a)') "Angular pixel resolution : ", (dxInCm/griddistance)*radiansToArcSec, " arcseconds"
       call writeinfo(message,TRIVIAL)
    else
       write(message,*) "Distance to grid should be > 0 but it is ", gridDistance
       call writeWarning(message)
    endif

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
       cube%xAxis(:) = cube%xAxis(:) * 1.0e10_db / kpctocm
       cube%yAxis(:) = cube%yAxis(:) * 1.0e10_db / kpctocm
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

    case('torus')
       continue

    case DEFAULT 
       call writewarning('Unrecognised unit in convertSpatialAxes')

    end select

  end subroutine convertSpatialAxes

  subroutine convertVelocityToHz(cube, freq)

    type(DATACUBE) :: cube
    real(double) :: freq

    cube%vAxis(:) = freq*(1.d0 + ((cube%vAxis(:)*1.d5)/cspeed))

  end subroutine convertVelocityToHz


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

  subroutine addWavelengthAxis(cube, lambdaArray)
    type(DATACUBE) :: cube
    real(double):: lambdaArray(:)
    integer :: i
    do i = 1, cube%nv
       cube%vAxis(i) = lambdaArray(i)/1.e4
    enddo
    cube%vUnit    = "microns"

  end subroutine addWavelengthAxis

  subroutine convertIntensityToBrightnessTemperature(cube, thisWavelength)
    use constants_mod, only: kErg
    
    type(DATACUBE) :: cube
    real(double), intent(in) :: thisWavelength

    cube%intensity(:,:,:) = real(cube%intensity(:,:,:) * (thisWavelength**2) / (2.0 * kErg))

  end subroutine convertIntensityToBrightnessTemperature

subroutine TranslateCubeIntensity(cube,constant)

  type(DATACUBE) :: cube
  real(double) :: constant
  
  cube%intensity = cube%intensity + real(constant)
    
end subroutine TranslateCubeIntensity

subroutine freeDataCube(thiscube)
  type(DATACUBE) :: thiscube

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
    real(double) :: f
    integer :: i

    allocate(v(1:thisCube%nv))
    allocate(y(1:thisCube%nv))
    do i = 1, thisCube%nv
       v(i) = thisCube%vAxis(i)
       y(i) = SUM(thisCube%intensity(:,:,i))
    enddo
    f = y(1)
    open(38, file=specfile, status="unknown",form="formatted")
    do i = 1, thisCube%nv
       write(38,*) v(i), y(i)/f, y(i)
    enddo
    deallocate(v, y)
    close(38)
  end subroutine dumpCubeToSpectrum
       
  subroutine writeSpectroastrometry(cube, saFilename)

    type(DATACUBE) :: cube
    character(len=*) :: saFilename
    real(double) :: meanPosX, meanPosY, tot
    integer :: ix, iy, iv

    open(34, file=saFilename)
    do iv = 1, cube%nv
       meanPosX = 0.d0
       meanPosY = 0.d0
       tot = 0.d0
       do ix = 1, cube%nx
          do iy = 1, cube%ny
             meanPosX = meanPosX + cube%intensity(ix,iy,iv)*cube%xAxis(ix)
             meanPosY = meanPosY + cube%intensity(ix,iy,iv)*cube%yAxis(iy)
             tot = tot + cube%intensity(ix,iy,iv)
          enddo
       enddo
       meanPosX = meanPosX / tot
       meanPosY = meanPosY / tot
       write(34, '(1p,3e12.3)') cube%vAxis(iv), meanPosX, meanPosY
    enddo
    close(34)
  end subroutine writeSpectroastrometry

end module datacube_mod
#endif
