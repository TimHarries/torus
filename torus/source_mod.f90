module source_mod
  
  use spectrum_mod
  use constants_mod
  use vector_mod
  use utils_mod

  implicit none

  public

  type SOURCETYPE
     type(VECTOR) :: position
     real(kind=doubleKind) :: radius
     real(kind=doubleKind) :: luminosity
     real(kind=doubleKind) :: teff
     type(SPECTRUMTYPE) :: spectrum
  end type SOURCETYPE

  contains

    function newSource(position, luminosity, radius, teff, spectrum)
      type(SOURCETYPE) :: newSource
      type(VECTOR), intent(in) :: position
      real(kind=doubleKind), intent(in) :: luminosity, radius, teff
      type(SPECTRUMTYPE) :: spectrum
      newSource%position = position
      newSource%radius = radius
      newSource%luminosity = luminosity
      newSource%teff = teff
      newSource%spectrum = spectrum
    end function newSource

    subroutine randomSource(source, nSource, iSource)
      integer :: nSource
      type(SOURCETYPE) :: source(1:nSource)
      integer :: iSource
      real, allocatable :: prob(:)
      real :: r, t
      real :: lRatio
      integer :: i

      if (nSource == 1) then
         iSource = 1
      else if (nSource == 2) then
         lRatio = source(1)%luminosity / (source(1)%luminosity + source(2)%luminosity)
         call random_number(r)
         if ( r  < lRatio ) then
            iSource = 1
         else
            iSource = 2
         endif
      else if (nSource > 2) then
         prob(1:nSource) = source(1:nSource)%luminosity
         do i = 2, nSource
            prob(i) = prob(i) + prob(i-1)
         enddo
         prob(1:nSource) = prob(1:nSource) - prob(1)
         prob(1:nSource) = prob(1:nSource) / prob(nSource)
         call random_number(r)
         call locate(prob, nSource, r, iSource)
         t = (r - prob(iSource))/(prob(iSource+1) - prob(iSource))
         if (t > 0.5) iSource = iSource + 1
      endif

    end subroutine randomSource

    subroutine getPhotonPositionDirection(source, position, direction, rHat)
      type(SOURCETYPE) :: source
      type(VECTOR) :: position, direction, rVec, rHat

      rHat = randomUnitVector()
      position = source%position + source%radius*rHat
      direction = fromPhotosphereVector(rHat)
    end subroutine getPhotonPositionDirection

  end module source_mod
