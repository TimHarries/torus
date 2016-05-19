module detector_mod
  use constants_mod
  use vector_mod
  use spectrum_mod
  implicit none

  type DETECTORTYPE

     character(len=80) :: type
     type(VECTOR) :: centre
     type(VECTOR) :: normal

     type(VECTOR) :: xAxis
     type(VECTOR) :: yAxis


     real(double) :: xSize
     real(double) :: ySize
     integer :: nx
     integer :: ny

     real(double),pointer :: intensity(:,:)

     real(double),pointer :: rIntensity(:,:)
     real(double),pointer :: gIntensity(:,:)
     real(double),pointer :: bIntensity(:,:)

     type(SPECTRUMTYPE) :: spectrum
     real(double) :: detRadius
     real(double) :: detectorNA


  end type DETECTORTYPE

contains

  subroutine createNewDetector(thisDetector)

    use inputs_mod, only : detectorPosition, detectorTheta, detectorPhi, &
         detectorXsize, detectorYsize, detectorNx, detectorNy, detType, detectorNA, &
         detRadius
    type(DETECTORTYPE) :: thisDetector
    type(VECTOR), parameter :: yAxis = VECTOR(0.d0, 1.d0, 0.d0)
    
    thisDetector%type = detType
    thisDetector%centre = detectorPosition
    thisDetector%normal = VECTOR(cos(detectorPhi)*sin(detectorTheta), &
                                 sin(detectorPhi)*sin(detectorTheta), &
                                 cos(detectorTheta))

    select case(detType)
       case("ccd")
          thisDetector%xAxis = thisDetector%normal .cross. yAxis
          if (modulus(thisDetector%xAxis) == 0.d0) then
             thisDetector%xAxis = VECTOR(1.d0, 0.d0, 0.d0)
          endif
          call normalize(thisDetector%xAxis)
          thisDetector%yAxis = thisDetector%normal .cross. thisDetector%xAxis
          call normalize(thisDetector%yAxis)
          thisDetector%nx = detectornx
          thisDetector%ny = detectorny
          thisDetector%xSize = detectorxSize
          thisDetector%ySize = detectorYsize
          allocate(thisDetector%intensity(thisDetector%nx, thisDetector%ny))
          allocate(thisDetector%rIntensity(thisDetector%nx, thisDetector%ny))
          allocate(thisDetector%gIntensity(thisDetector%nx, thisDetector%ny))
          allocate(thisDetector%bIntensity(thisDetector%nx, thisDetector%ny))
          thisDetector%intensity = 0.d0
          thisDetector%rintensity = 0.d0
          thisDetector%gintensity = 0.d0
          thisDetector%bintensity = 0.d0
       case("fibre")
          thisDetector%detRadius = detRadius
          thisDetector%detectorNA = detectorNA
          call newSpectrum(thisDetector%spectrum, 400.d0, 800.d0, 1000)
       case DEFAULT
          call writeFatal("Detector type not recognised "//trim(detType))
       end select
  end subroutine createNewDetector
    
    
  logical function hitDetector(thisDetector, rVec, uHat)
    type(DETECTORTYPE) :: thisDetector
    type(VECTOR) :: rVec, uHat, towardsDetector, hitPosition
    real(double) :: distanceToPlane, xCoord, yCoord

    hitDetector = .false.

    towardsDetector = thisDetector%centre - rVec
    call normalize(towardsDetector)

    if ( ((uHat.dot.towardsDetector) > 0.d0) .and. &
         ((uHat.dot.thisDetector%normal) < 0.d0) ) then ! photon is moving towards detector (from the correct side!)
       select case(thisDetector%type)
          case("ccd")

             distanceToPlane = ((thisDetector%centre - rVec).dot.thisDetector%normal)/(uHat.dot.thisDetector%normal)
             hitPosition = rVec + distanceToPlane * uHat
             
             xCoord = (hitPosition - thisDetector%centre).dot.thisDetector%xAxis
             yCoord = (hitPosition - thisDetector%centre).dot.thisDetector%yAxis
             
             if ( (abs(xCoord) < thisDetector%xSize/2.d0) .and. (abs(yCoord) < thisDetector%ySize/2.d0) ) then
                hitDetector = .true.
             endif

          case("fibre")
             distanceToPlane = ((thisDetector%centre - rVec).dot.thisDetector%normal)/(uHat.dot.thisDetector%normal)
             hitPosition = rVec + distanceToPlane * uHat
             if (modulus(hitPosition-thisDetector%centre) < thisDetector%detRadius) then
                if (acos(-(uhat.dot.thisdetector%normal)) < asin(thisDetector%detectorna)) then
                   hitDetector = .true.
                endif
             endif
          end select
       endif
  end function hitDetector


  subroutine storePhotonOnDetector(thisDetector, photonPos, photonDirection, photonPower, photonWavelength, stored)
    type(DETECTORTYPE) :: thisDetector
    type(VECTOR) :: photonPos, photonDirection, hitPosition
    real(double) :: photonPower, distanceToPlane, xCoord, yCoord, photonWavelength
    integer :: ix, iy
    logical :: stored

    stored = .false.


    if (hitDetector(thisDetector, photonPos, photonDirection)) then
       stored = .true.
       select case(thisDetector%type)
          case("ccd") 
             distanceToPlane = ((thisDetector%centre - photonPos).dot.thisDetector%normal)/(photonDirection.dot.thisDetector%normal)
             hitPosition = photonPos + distanceToPlane * photonDirection
             
             xCoord = (hitPosition - thisDetector%centre).dot.thisDetector%xAxis
             yCoord = (hitPosition - thisDetector%centre).dot.thisDetector%yAxis
             
             ix = int(real(thisDetector%nx)*(xCoord + thisDetector%xSize/2.d0)/thisDetector%xSize)+1
             iy = int(real(thisDetector%ny)*(yCoord + thisDetector%ySize/2.d0)/thisDetector%ySize)+1
             thisDetector%intensity(ix, iy) = thisDetector%intensity(ix, iy) + photonPower
             thisDetector%rIntensity(ix, iy) = thisDetector%rIntensity(ix, iy) + photonPower*responseFunction("R", photonWavelength)
             thisDetector%gIntensity(ix, iy) = thisDetector%gIntensity(ix, iy) + photonPower*responseFunction("G", photonWavelength)
             thisDetector%bIntensity(ix, iy) = thisDetector%bIntensity(ix, iy) + photonPower*responseFunction("B", photonWavelength)
          case("fibre")
             call addPhotonToSpectrum(thisDetector%spectrum, photonWavelength, photonPower)
       end select
    endif

  end subroutine storePhotonOnDetector

  real(double) function ResponseFunction(band, wavelength)
    character(len=1) :: band
    real(double) :: wavelength

    select case(band)
       case("R")
          responseFunction = gaussianResponse(600.d0, 50.d0, wavelength)
       case("G")
          responseFunction = gaussianResponse(550.d0, 50.d0, wavelength)
       case("B")
          responseFunction = gaussianResponse(450.d0, 50.d0, wavelength)
       case DEFAULT
          call writeFatal("Band not recognised: "//band)
       end select
   end function ResponseFunction

real(double) function gaussianResponse(centre, sigma, x)
  real(double) :: centre, sigma, x
  gaussianResponse = 1.d0/(sigma*sqrt(twoPi)) * exp(-(centre-x)**2/(2.d0*sigma**2))
end function gaussianResponse

subroutine writeDetectorFile(thisDetector, detFilename)
  type(DETECTORTYPE) :: thisDetector
  character(len=*) :: detFilename
  select case(thisdetector%type)
     case("ccd")
        call writePPMfile(thisDetector, detFilename)
     case("fibre")
        call writeSpectrumFile(thisDetector, detFilename)
  end select 
end subroutine writeDetectorFile

subroutine writeSpectrumFile(thisDetector, detFilename)
  type(DETECTORTYPE) :: thisDetector
  character(len=*) :: detFilename
  character(len=80) :: tmpfile
  integer :: i
  tmpfile = trim(detFilename)//".dat"
  open(22, file=tmpFile, status="unknown", form="formatted")
  do i = 1, thisDetector%spectrum%nlambda
     write(22,*) thisDetector%spectrum%lambda(i),thisDetector%spectrum%flux(i)
  enddo
  close(22)
end subroutine writeSpectrumFile

subroutine writePPMfile(thisDetector, filename)
  type(DETECTORTYPE) :: thisDetector
  character(len=*) :: filename
  integer, allocatable :: ir(:,:), ig(:,:), ib(:,:)
  real(double) :: norm
  integer :: nx, ny

  nx = size(thisDetector%rintensity,1)
  ny = size(thisDetector%rintensity,2)

  allocate(ir(1:nx,1:ny))
  allocate(ig(1:nx,1:ny))
  allocate(ib(1:nx,1:ny))
  norm = MAX(MAXVAL(thisDetector%rIntensity),MAXVAL(thisDetector%gIntensity),MAXVAL(thisDetector%bIntensity))
  if (norm == 0.d0) norm = 1.d0
  ir = nint((thisDetector%rIntensity/norm)*dble(2**16-1))
  ig = nint((thisDetector%gIntensity/norm)*dble(2**16-1))
  ib = nint((thisDetector%bIntensity/norm)*dble(2**16-1))

  call writeppm3(ir, ig, ib, filename)
  deallocate(ir, ig, ib)
end subroutine writePPMfile


subroutine writeppm3(R,G,B,text)
  integer :: R(:,:),G(:,:),B(:,:)
  character(len=*) :: text
  integer :: cols,rows
  integer :: i,j
  integer :: maxvalue

  ! Open File   
  open(unit=100, file=trim(text)//".ppm", status='unknown')

  ! Write Header and ppm file type
  write(100,'( A )') "P3"
  write(100,'( A )') "# PPM Type 2 File (generated with fortran)"

  ! Write Image Size
  cols = size(R,2)
  rows = size(R,1)
  write(100,'( i5, 1x, i5 )') cols, rows

  ! Write Maximum Value
  maxvalue = max( maxval(maxval(R,dim=1),dim=1)&
       ,maxval(maxval(G,dim=1),dim=1)&
       ,maxval(maxval(B,dim=1),dim=1))
  write(100,'( i5 )') maxvalue

  ! Write Image
  do i=1,rows
     do j=1,cols
        write(100,'( 3(i5,1x) )') R(i,j),G(i,j),B(i,j)
     enddo
  enddo
  close(100)
end subroutine writeppm3

end module detector_mod
