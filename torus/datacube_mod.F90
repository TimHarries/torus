module datacube_mod

  use kind_mod
  use vector_mod
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
     character(len=10) :: FluxUnit ! units for flux

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
     real, pointer :: flux(:,:,:) => null()
     real, pointer :: tau(:,:,:) => null()
     real, pointer :: nCol(:,:) => null()
     real, pointer :: i0_pos(:,:,:) => null() ! Positive contribution to flux
     real, pointer :: i0_neg(:,:,:) => null() ! Negative contribution to flux

  end type DATACUBE

contains

  subroutine writeDataCube(thisCube, filename, write_Intensity, write_ipos, write_ineg, write_Tau, write_nCol, write_axes)

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

#ifdef USECFITSIO
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

    status=0
    
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

       bitpix=-32
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
       
       call addWCSinfo

       !  Write the array to the FITS file.
       call ftppre(unit,group,fpixel,nelements,thisCube%intensity,status)
              
    endif

    if( do_write_Tau ) then

       ! 2nd HDU : tau
       call FTCRHD(unit, status)
       bitpix=-32
       naxis=3
       naxes(1)=thisCube%nx
       naxes(2)=thisCube%ny
       naxes(3)=thisCube%nv
       nelements=naxes(1)*naxes(2)*naxes(3)

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       call ftpkys(unit,'BUNIT',"Optical depth","Optical depth unit",status)
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
       bitpix=-32
       naxis=2
       naxes(1)=thisCube%nx
       naxes(2)=thisCube%ny
       nelements=naxes(1)*naxes(2)
       
       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       call ftpkys(unit,'BUNIT',"cm^-2","Column density unit",status)
       call addWCSinfo

       !  Write the array to the FITS file.
       call ftppre(unit,group,fpixel,nelements,thisCube%nCol,status)
       call print_error(status)
    endif

    if( do_write_ipos) then

       ! 9th HDU : positive intensity contribution
       call FTCRHD(unit, status)
       bitpix=-32
       naxis=3
       naxes(1)=thisCube%nx
       naxes(2)=thisCube%ny
       naxes(3)=thisCube%nv
       nelements=naxes(1)*naxes(2)*naxes(3)

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       call ftpkys(unit,'BUNIT',thisCube%IntensityUnit,"Intensity unit",status)
       call addWCSinfo

       !  Write the array to the FITS file.
       call ftppre(unit,group,fpixel,nelements,thisCube%i0_pos,status)
    endif


    if( do_write_ineg) then

       ! 10th HDU : negative intensity contribution
       call FTCRHD(unit, status)
       bitpix=-32
       naxis=3
       naxes(1)=thisCube%nx
       naxes(2)=thisCube%ny
       naxes(3)=thisCube%nv
       nelements=naxes(1)*naxes(2)*naxes(3)

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       call ftpkys(unit,'BUNIT',thisCube%IntensityUnit,"Intensity unit",status)
       call addWCSinfo

       !  Write the array to the FITS file.
       call ftppre(unit,group,fpixel,nelements,thisCube%i0_neg,status)
    endif

    !  Close the file and free the unit number.
    call ftclos(unit, status)
    call ftfiou(unit, status)

     !  Check for any error, and if so print out error messages
    if (status > 0) then
       call print_error(status)
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
        call ftpkyd(unit,'CDELT3',thisCube%vAxis(2)-thisCube%vAxis(1),-3,'coordinate increment at reference point',status)
        call ftpkys(unit,'CTYPE3',thisCube%vAxisType, "velocity axis", status)
        call ftpkyd(unit,'CRVAL3',thisCube%vAxis(1),-3,'coordinate value at reference point',status)

      end subroutine addWCSinfo

#endif
    
    
  end subroutine writeDataCube

  !**********************************************************************

  subroutine readDataCube(thisCube)

    implicit none
    
#ifdef USECFITSIO
    type(DATACUBE), intent(out) :: thisCube
    integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,nbuffer,npixels, hdutype, hdu, junk
    character(len=512) :: filename
    character(len=80) :: comment
    real :: nullval
    integer, dimension(4) :: naxes
    logical :: anynull
    character(len=80) :: keyword

    status=0
    !  Get an unused Logical Unit Number to use to open the FITS file.
    call ftgiou(unit,status)

    filename="datacube.fits.gz"
    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)

    group=1
    firstpix=1
    nullval=-999
     

    ! 1st HDU : intensity
    hdu=1

    call ftgcrd(unit,keyword,status)

    write(*,*) "Keyword", keyword

    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    if (nfound /= 3) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of datacube.fits.gz file HDU',hdu,'. Exiting.'
       stop
    endif
    npixels=naxes(1)*naxes(2)*naxes(3)
    nbuffer=npixels
    ! read_image
    call ftgpve(unit,group,firstpix,nbuffer,nullval,thisCube%intensity,anynull,status)


    !  Read keywords from the header.
    call FTGKYJ(unit,"LABEL", junk,comment,status)
    thisCube%label = comment
    call FTGKYJ(unit,"VUNIT", junk,comment,status)
    thisCube%vUnit = comment
    call FTGKYJ(unit,"XUNIT", junk,comment,status)
     thisCube%xUnit = comment
    call FTGKYJ(unit,"IUNIT", junk,comment,status)
    thisCube%IntensityUnit = comment
    call FTGKYD(unit,"DISTANCE", thisCube%obsdistance,comment,status)

    ! 2nd HDU : converged
    hdu=2
    call ftmahd(unit,hdu,hdutype,status)
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    if (nfound /= 3) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of datacube.fits.gz file HDU',hdu,'. Exiting.'
       stop
    endif
    npixels=naxes(1)*naxes(2)*naxes(3)
    nbuffer=npixels
    ! read_image
    call ftgpvj(unit,group,firstpix,nbuffer,nullval,thisCube%converged,anynull,status)

    ! 3rd HDU : weight
    hdu=3
    call ftmahd(unit,hdu,hdutype,status)
    call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
    if (nfound /= 2) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of datacube.fits.gz file HDU',hdu,'. Exiting.' 
       stop
    endif
    npixels=naxes(1)*naxes(2)
    nbuffer=npixels
    ! read_image
    call ftgpvd(unit,group,firstpix,nbuffer,nullval,thisCube%weight,anynull,status)

    ! 4th HDU : xAxis
    hdu=4
    call ftmahd(unit,hdu,hdutype,status)
    call ftgknj(unit,'NAXIS',1,1,naxes,nfound,status)
    if (nfound /= 1) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of datacube.fits.gz file HDU',hdu,'. Exiting.'     
       stop
    endif
    npixels=naxes(1)
    nbuffer=npixels
    ! read_image
    call ftgpvd(unit,group,firstpix,nbuffer,nullval,thisCube%xAxis,anynull,status)
    
    ! 5th HDU : yAxis
    hdu=5
    call ftmahd(unit,hdu,hdutype,status)
    call ftgknj(unit,'NAXIS',1,1,naxes,nfound,status)
    if (nfound /= 1) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of datacube.fits.gz file HDU',hdu,'. Exiting.'   
       stop
    endif
    npixels=naxes(1)
    nbuffer=npixels
    ! read_image
    call ftgpvd(unit,group,firstpix,nbuffer,nullval,thisCube%yAxis,anynull,status)
    
    ! 6th HDU : vAxis
    hdu=6
    call ftmahd(unit,hdu,hdutype,status)
    call ftgknj(unit,'NAXIS',1,1,naxes,nfound,status)
    if (nfound /= 1) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of datacube.fits.gz file HDU',hdu,'. Exiting.' 
       stop
    endif
    npixels=naxes(1)
    nbuffer=npixels
    ! read_image
    call ftgpvd(unit,group,firstpix,nbuffer,nullval,thisCube%vAxis,anynull,status)
    
    ! 7th HDU : intensity
    hdu=7
    call ftmahd(unit,hdu,hdutype,status)
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    if (nfound /= 3) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of datacube.fits.gz file HDU', hdu,'. Exiting.'
       stop
    endif
    npixels=naxes(1)*naxes(2)*naxes(3)
    nbuffer=npixels
    ! read_image
    call ftgpvj(unit,group,firstpix,nbuffer,nullval,thisCube%nsubpixels,anynull,status)   

#else
    type(DATACUBE) :: thisCube

#endif
    return

  end subroutine readDataCube

  !**********************************************************************


  subroutine print_error(status)
    ! PRINT_ERROR prints out the FITSIO error messages to the user.
    
    integer status
#ifdef USECFITSIO
    character ( len = 30 ) errtext
    character ( len = 80 ) errmessage

    !  Check if status is OK (no error); if so, simply return.
    if (status <= 0) then
       return
    end if

    !  Get the text string which describes the error
    call ftgerr(status,errtext)
    print *,'FITSIO Error Status =',status,': ',errtext

    !  Read and print out all the error messages on the FITSIO stack
    call ftgmsg(errmessage)
    do while (errmessage .ne. ' ')
       print *,errmessage
       call ftgmsg(errmessage)
    end do
#endif    

    return
  end subroutine print_error

  !***********************************************************

! Initialises cube - sets intensity for cube to 0 
  subroutine initCube(thisCube, nx, ny, nv, mytelescope)
    use input_variables, only: splitCubes
    type(DATACUBE) :: thisCube
    type(TELESCOPE), optional :: mytelescope
    integer :: nx, ny, nv
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
    thisCube%FluxUnit     = " "

    thisCube%nx = nx
    thisCube%ny = ny
    thisCube%nv = nv
    allocate(thisCube%xAxis(1:nx))
    allocate(thisCube%yAxis(1:ny))
    allocate(thisCube%vAxis(1:nv))
    allocate(thisCube%intensity(1:nx,1:ny,1:nv))
    allocate(thisCube%flux(1:nx,1:ny,1:nv))
    allocate(thisCube%tau(1:nx,1:ny,1:nv))
    allocate(thisCube%nCol(1:nx,1:ny))
!    allocate(thisCube%nsubpixels(1:nx,1:ny,1:nv))
!    allocate(thisCube%converged(1:nx,1:ny,1:nv))
!    allocate(thisCube%weight(1:nx,1:ny))

    thisCube%intensity = 0.d0
    thisCube%tau =  0.d0
    thisCube%nCol = 0.d0

!    thisCube%flux = 0.d0
!    thisCube%nsubpixels = 0.d0
!    thisCube%converged = 0
!    thisCube%weight = 1.d0

    if (splitCubes) then 
       ! flux component not required in configurations which use split cubes 
       deallocate (thisCube%flux)
       nullify (thisCube%flux)
       allocate(thisCube%i0_pos(1:nx,1:ny,1:nv))
       allocate(thisCube%i0_neg(1:nx,1:ny,1:nv))
       thisCube%i0_pos(:,:,:) = 0.0
       thisCube%i0_neg(:,:,:) = 0.0
    end if

  end subroutine initCube

! Set spatial axes for datacube - Equally spaced (linearly) between min and max
  subroutine addSpatialAxes(cube, xMin, xMax, yMin, yMax)
    use input_variables , only : gridDistance
    use constants_mod, only: auTocm, pi, rSol
    type(DATACUBE) :: cube
    real(double) :: xMin, xMax, yMax, yMin, dx, dy
    integer :: i
    character(len=100) :: message

    dx = (xMax - xMin)/dble(cube%nx)
    dy = (yMax - yMin)/dble(cube%ny)

    write(message,'(a,f7.3,a)') "Linear pixel resolution  : ", dx*1e10/autocm, " AU"
    call writeinfo(message,TRIVIAL)
    write(message,'(a,f7.3,a)') "Linear pixel resolution  : ", dx*1e10/rSol, " Rsol"
    call writeinfo(message,TRIVIAL)
    write(message,'(a,f10.4,a)') "Angular pixel resolution : ", (dx*1e10/griddistance)*(180./pi)*60.*60., " arcseconds"
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

! Set velocity axis for datacube - Equally spaced (linearly) between min and max
  subroutine addVelocityAxis(cube, vMin, vMax)
    type(DATACUBE) :: cube
    real(double):: vMin, vMax, dv
    integer :: i
    character(len=100) :: message

    dv = (vMax - vMin) / dble(cube%nv)

    write(message, '(a,f7.4,a)') "Velocity pixel resolution: ", dv, " km/s"
    call writeinfo(message,TRIVIAL)

    do i = 1, cube%nv
       cube%vAxis(i) = vmin + dv * real(i-1)
    enddo
  end subroutine addVelocityAxis

  subroutine getSpectrum(cube, ix1, ix2, iy1, iy2, spec)
    type(DATACUBE) :: cube
    integer :: ix1, ix2, iy1, iy2
    real(double) ::  spec(:)
    integer :: i,j, n
    n = 0
    spec = 0.d0
    do i = ix1, ix2
       do j = iy1, iy2
          n = n + 1
          spec(1:cube%nv) = spec(1:cube%nv) + cube%flux(i,j,1:cube%nv)
       enddo
    enddo
!    spec = spec / dble(n)
  end subroutine getSpectrum

  subroutine getWeightedSpectrum(cube, ix1, ix2, iy1, iy2, spec)
    type(DATACUBE) :: cube
    integer :: ix1, ix2, iy1, iy2
    real(double) ::  spec(:) , tot
    integer :: i,j
    tot = 0.d0
    spec = 0.d0
    do i = ix1, ix2
       do j = iy1, iy2
          tot =  tot + cube%weight(i,j)
          spec(1:cube%nv) = spec(1:cube%nv) + cube%intensity(i,j,1:cube%nv) * cube%weight(i,j)
       enddo
    enddo
!    spec = spec / dble(tot)
  end subroutine getWeightedSpectrum

  subroutine convolveCube(cube, beamSize)
    use constants_mod, only: pcTocm, pi, twoPi
    type(DATACUBE) :: cube
    real(double) :: beamSize ! beamsize in arcsec
    real(double), allocatable :: newArray(:,:)
    real(double) :: rrinArcSec, fac 
    integer :: ix, iy, iv, i, j
    real(double) :: sigma, sigma2, dx, dy, tot, flux, background
    real(double) :: deltaX, deltaY

    sigma = beamsize/2.35d0 ! changed from 2.35 (FWHM) as a result of reading IRAM-30m paper
    sigma2 = sigma**2
    dx = 3600.d0*((cube%xAxis(2) - cube%xAxis(1))/(cube%obsDistance/1.d10))*180.d0/pi
    dy = 3600.d0*((cube%yAxis(2) - cube%yAxis(1))/(cube%obsDistance/1.d10))*180.d0/pi

    allocate(newArray(1:cube%nx, 1:cube%ny))
    call writeInfo("Convolving data cube with beam size", TRIVIAL)
    write(*,*) "cube%obsdist",cube%obsDistance/pctocm, "dx",dx,"dy",dy

    do iv = 1, cube%nv

       background = cube%flux(cube%nx,cube%ny,iv)
       newArray = 0.d0

       do ix = 1, cube%nx
          do iy = 1, cube%ny
             tot = 0.d0
             do i = ix - cube%nx, ix + cube%nx
                do j = iy - cube%ny, iy + cube%ny
                   if ((i > 0).and.(i<=cube%nx).and.(j > 0).and.(j <= cube%ny)) then
                      flux = cube%flux(i,j,iv)
                   else
                      flux = background
                   endif

                   deltaX = dble(ix-i)*dx
                   deltaY = dble(iy-j)*dy
                   rrInArcSec = deltaX**2 + deltaY**2
             
                   fac = (1.d0/(twoPi*sigma2))*exp(-0.5d0*(rrInArcSec/sigma2))*dx*dy

                   newArray(ix,iy) = newArray(ix, iy) + flux*fac
                   tot = tot + fac
                enddo
             enddo
          enddo
          !write(*,*) "Weight",tot
       enddo

       cube%flux(1:cube%nx, 1:cube%ny, iv) = newArray(1:cube%nx, 1:cube%ny)!/dble(cube%nx*cube%ny)
    enddo
    deallocate(newArray)
    call writeInfo("Done.",TRIVIAL)
  end subroutine convolveCube

subroutine TranslateCubeIntensity(cube,constant)

  type(DATACUBE) :: cube
  real(double) :: constant
  
  cube%intensity = cube%intensity + constant
    
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

    if (associated(thisCube%flux)) then 
       deallocate(thiscube%flux)
       nullify(thiscube%flux)
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

end module datacube_mod

