module source_mod
  
  use spectrum_mod
  use constants_mod
  use vector_mod
  use utils_mod
  use octal_mod
  use surface_mod

  implicit none

  public 


  type SOURCETYPE
     type(VECTOR)    :: position   ! [10^10cm]
     real(double) :: radius     ! [10^10cm]
     real(double) :: luminosity ! [erg/s]
     real(double) :: teff       ! [K]
     real(double) :: initialMass ! [solar]
     real(double) :: mass       ! [solar]
     real(double) :: age        ! [years]
     type(SPECTRUMTYPE)    :: spectrum   ! [???]
     type(SURFACETYPE) :: surface
     logical :: outsideGrid
  end type SOURCETYPE



  contains

    function newSource(position, luminosity, radius, teff, spectrum)
      type(SOURCETYPE) :: newSource
      type(VECTOR), intent(in) :: position
      real(double), intent(in) :: luminosity, radius, teff
      type(SPECTRUMTYPE) :: spectrum
      newSource%position = position
      newSource%radius = radius
      newSource%luminosity = luminosity
      newSource%teff = teff
      newSource%spectrum = spectrum
    end function newSource






!    subroutine randomSource(source, nSource, iSource)
!      integer :: nSource
!      type(SOURCETYPE) :: source(1:nSource)
!      integer :: iSource
!      real, save, allocatable :: prob(:)
!      real :: r, t
!      real :: lRatio
!      integer :: i
!      logical, save :: first_time = .true.
!
!      if (nSource == 1) then
!         iSource = 1
!      else if (nSource == 2) then
!         lRatio = source(1)%luminosity / (source(1)%luminosity + source(2)%luminosity)
!         call random_number(r)
!         if ( r  < lRatio ) then
!            iSource = 1
!         else
!            iSource = 2
!         endif
!      else if (nSource > 2) then
!	 if (first_time) then
!	    ! allocate array
!	    ALLOCATE(prob(1:nSource))
!	    ! Create the prob. dist. function.
!	    prob(1:nSource) = source(1:nSource)%luminosity/lsol
!	    do i = 2, nSource
!	       prob(i) = prob(i) + prob(i-1)
!	    enddo
!	    prob(1:nSource) = prob(1:nSource) - prob(1)
!	    prob(1:nSource) = prob(1:nSource) / prob(nSource)
!
!	    first_time = .false.
!	    
!	 end if
!	 
!         call random_number(r)
!         call locate(prob, nSource, r, iSource)
!         if (iSource < nSource) then
!            t = (r - prob(iSource))/(prob(iSource+1) - prob(iSource))
!            if (t > 0.5) iSource = iSource + 1
!         endif
!      endif
!
!    end subroutine randomSource


    subroutine randomSource(source, nSource, iSource, lamArray, nLambda, initialize)
      integer :: nSource
      type(SOURCETYPE) :: source(1:nSource)
      integer :: iSource
      real, save, allocatable :: prob(:)
      real, optional :: lamArray(:)
      integer,optional :: nlambda
      real :: r, t
      integer :: i
      logical, optional :: initialize

      if (nSource == 1) then
         iSource = 1
      else
	 if (PRESENT(initialize)) then
	    ! allocate array
	    if (allocated(prob)) then
               deallocate(prob)
            endif
            ALLOCATE(prob(1:nSource))
	    ! Create the prob. dist. function.
            do i = 1, nSource
               prob(i) = integrateSpectrumOverBand(source(i)%spectrum, dble(lamArray(1)) , &
                    dble(lamArray(nLambda))) /lsol
            enddo
	    do i = 2, nSource
	       prob(i) = prob(i) + prob(i-1)
	    enddo
	    prob(1:nSource) = prob(1:nSource) - prob(1)
	    prob(1:nSource) = prob(1:nSource) / prob(nSource)
!            write(*,*) "Source prob: ",prob(1:nSource)
         end if
	 
         call random_number(r)
         call locate(prob, nSource, r, iSource)
         if (iSource < nSource) then
            t = (r - prob(iSource))/(prob(iSource+1) - prob(iSource))
            if (t > 0.5) iSource = iSource + 1
         endif
      endif

    end subroutine randomSource

    real(double) function sumSourceLuminosity(source, nsource, lam1, lam2) result (tot)
      integer :: nSource
      type(SOURCETYPE) :: source(:)
      real :: lam1, lam2
      integer :: i

      tot = 0
      do i = 1, nSource
         tot = tot + integrateSpectrumOverBand(source(i)%spectrum, dble(lam1) , &
                    dble(lam2)) * (fourPi*(1.d10*source(i)%radius)**2)
      enddo
    end function sumSourceLuminosity

    real(double) function sumSourceLuminosityMonochromatic(source, nsource, lam) result (tot)
      integer :: nSource
      type(SOURCETYPE) :: source(:)
      real(double) :: lam, flux
      integer :: i !, j

      tot = 0
      do i = 1, nSource
         if (lam > source(i)%spectrum%lambda(source(i)%spectrum%nLambda)) then
           tot = tot + 1.d-200
         else
!           call locate(source(i)%spectrum%lambda, source(i)%spectrum%nLambda, lam, j)

           flux =  loginterp_dble(source(i)%spectrum%flux(1:source(i)%spectrum%nLambda), &
                source(i)%spectrum%nLambda, source(i)%spectrum%lambda(1:source(i)%spectrum%nLambda), lam)
           tot = tot + flux * fourPi * (source(i)%radius*1.d10)**2
         endif
      enddo
    end function sumSourceLuminosityMonochromatic




    subroutine distanceToSource(source, nSource, rVec, uHat, hitSource, distance, sourcenumber)
      integer :: nSource
      type(SOURCETYPE) :: source(nSource)
      type(VECTOR) :: rVec, uHat
      integer, optional :: sourceNumber
      real(double) :: distance, cosTheta, sintheta
      logical :: hitSource
      integer :: i

      hitSource = .false.
      if (present(sourceNumber)) sourceNumber = 0 
      distance = 1.d30
      mainloop: do i = 1, nSource
         cosTheta = (uHat .dot. (source(i)%position - rVec))/modulus(source(i)%position - rVec)
         if (cosTheta < 0.) then ! source is behind you
            cycle mainloop
         endif
         sinTheta = sqrt(max(0.d0,1.d0 - cosTheta*cosTheta))
         if (source(i)%radius > (modulus(source(i)%position - rVec)*sinTheta)) then
            distance = min(distanceToSphere(rVec, uHat, source(i)%position, source(i)%radius),distance)
            hitSource = .true. ! bang
            if (present(sourcenumber)) sourceNumber = i
         endif
      enddo mainloop
    end subroutine distanceToSource

    subroutine getPhotonPositionDirection(source, position, direction, rHat, grid)
      use input_variables, only : pointSource
      type(GRIDTYPE) :: grid
      real(double) :: r 
      type(SOURCETYPE) :: source
      type(VECTOR) :: position, direction, rHat

      if (.not.source%outsideGrid) then
         if (pointSource) then
            !      ! simply treating as a point source
            position = source%position
            direction = randomUnitVector()
            
         else
            rHat = randomUnitVector()
            position = source%position + source%radius*rHat
            ! A limb darkening law should be applied here for 
            ! for general case here.....
            direction = fromPhotoSphereVector(rHat)
         endif
      else
         position%x = -grid%octreeRoot%subcellSize+1.d-10*grid%octreeRoot%subcellSize
         call random_number(r)
         r = 2.d0 * r - 1.d0
         position%y = r * grid%octreeRoot%subcellSize
         call random_number(r)
         r = 2.d0 * r - 1.d0
         position%z = r * grid%octreeRoot%subcellSize
         direction = VECTOR(1.d0, 0.d0, 0.d0)
      endif
    end subroutine getPhotonPositionDirection




  function distanceToSphere(rVec, uHat, centre, radius) result (distance)
    real(double) :: distance, radius
    real(double) :: x1, x2,a,b,c,d,costheta
    logical :: ok
    type(VECTOR) :: rVec, uHat, centre

    d = modulus(centre - rVec)
    cosTheta = (uHat .dot. (centre-rVec)) / d

    a = 1.d0
    b = -2.d0*d*cosTheta
    c = d**2 - radius**2
    call solveQuadDble(a, b, c, x1, x2, ok)
    if (.not.ok) then 
       write(*,*) "! Distance to sphere called with no intersection"
    endif
    if ((x1 > 0.).and.(x2 > 0.)) then
       distance = min(x1, x2)
    else
       if (x1 > 0.) then
          distance = x1
       else
          distance = x2
       endif
    endif
    if (distance < 0.) then
       write(*,*) "! screw up in distance to sphere"
    endif
  end function distanceToSphere




    !
    ! For a given octal and sourcetyupe objects, this routine checks 
    ! if the source is with in the range this octal.
    
    function source_within_octal(this, an_octal) RESULT(out)
      implicit none
      logical :: out 
      type(sourcetype), intent(in) :: this
      type(octal), intent(in) :: an_octal 
      
      real(oct) :: xc, yc, zc ! position of the octal center.
      real(oct) :: d          ! size of the subcell
      real(oct) :: x, y, z    ! position of the source.
      
      xc = an_octal%centre%x
      yc = an_octal%centre%y
      zc = an_octal%centre%z
      d  = an_octal%subcellsize

      x = this%position%x
      y = this%position%y
      z = this%position%z
      
    if ( x > (xc+d) ) then
       out = .false.
    else if ( x < (xc-d)) then
       out = .false.      
    elseif ( y > (yc+d) ) then
       out = .false.
    elseif ( y < (yc-d)) then
       out = .false.
    elseif ( z > (zc+d) ) then
       out = .false.
    elseif ( z < (zc-d)) then
       out = .false.
    else
       out = .true.
    end if

  end function source_within_octal




  !
  ! Pick a random direction (with length 1) from the surface of a sphere. 
  ! If the direction points inward of the sphere, it will be rejected, and a new one will 
  ! be calculated. 
  !
  ! note: rvec is a vector which connects a point on the surface of the sphere and its origin.
  !       The length of the rvec does not matther when passed to this routine as long as 
  !        it has a correct direction, this routine should work.
  function random_direction_from_sphere(rvec) RESULT(out)
    implicit none
    type(VECTOR) :: out 
    !
    type(VECTOR), intent(in) :: rvec
    !
    type(VECTOR) :: dir 
    real(oct) :: inner_product
    
    inner_product = -1.0d0
    do while (inner_product<0.0d0) 
       dir = randomUnitVECTOR()
       inner_product = dir .dot. rvec
    end do
    
    out = dir

  end function random_direction_from_sphere

  real(double) function I_nu(source, nu, iElement) 
    type(SOURCETYPE) :: source
    real(double) :: nu, fnu !, flambda, lam
    integer :: i, iElement
    real(double) :: tAccretion, ic_hot

!    lam = 1.d8 * cSpeed/ nu ! angs
!    if (lam < source%spectrum%lambda(1)) then
!       I_nu = tiny(i_nu)
!    endif
!    if (lam > source%spectrum%lambda(source%spectrum%nlambda)) then
!       I_nu = tiny(i_nu)
!    endif


!    call locate(source%spectrum%lambda, source%spectrum%nLambda, lam, i)
!    fLambda = source%spectrum%flux(i)
!    fnu = flambda * (cSpeed*1.d8) /nu**2 ! 1.d8 to go from cm/s to angs/c
!    i_nu = 1.5d0 * fnu / pi

!    write(*,*) "i_nu ",nu, source%surface%nuArray(1), &
!         source%surface%nuArray(source%surface%nNuHotFlux)

    if (nu < source%surface%nuArray(1)) then
       i_nu = tiny(i_nu)
!       write(*,*) nu,source%surface%nuArray(1)
    else if (nu > source%surface%nuArray(source%surface%nNuHotFlux)) then
       i_nu = tiny(i_nu)
!       write(*,*) nu,source%surface%nuArray(source%surface%nNuHotFlux)
    else
       call locate(source%surface%nuArray, source%surface%nNuHotFlux, real(nu), i)
       fnu = source%surface%hnuArray(i) 
       i_nu = 1.5d0 * fnu / pi
    endif

    if (isHot(source%surface,ielement)) then
       tAccretion = source%surface%element(ielement)%temperature
       !================CHECK UNITS HERE!! ===========================
       IC_hot = blackbody(REAL(tAccretion), 1.e8*REAL(cSpeed/nu)) ! [B_nu]
!       write(*,*) "hit hotspot",ic_hot,i_nu
       i_nu  = I_nu + ic_hot
    endif

  end function I_nu

  logical function insideSource(thisOctal, subcell, nsource, source)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    type(VECTOR) :: rVec, corner
    integer :: nsource
    type(SOURCETYPE) :: source(:)
    integer :: i
    real(double) :: r

    insideSource = .false.
    rvec = subcellCentre(thisOctal, subcell)

    do i = 1, nSource
       
       if (thisOctal%twoD) then
          corner = rVec + thisOctal%subcellSize/rootTwo * VECTOR(oneOverrootTwo,0.d0,oneOverrootTwo)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/rootTwo * VECTOR(oneOverrootTwo,0.d0,-oneOverrootTwo)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/rootTwo * VECTOR(-oneOverrootTwo,0.d0,oneOverrootTwo)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/rootTwo * VECTOR(-oneOverrootTwo, 0.d0, -oneOverrootTwo)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
       else if (thisOctal%oneD) then
          corner = rVec - VECTOR(thisOctal%subcellSize/2.d0,0.d0,0.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
       else
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(1.d0, 1.d0, 1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(1.d0, -1.d0, 1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(1.d0, 1.d0, -1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(1.d0, -1.d0, -1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(-1.d0, -1.d0, -1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(-1.d0, 1.d0, 1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(-1.d0, 1.d0, -1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif
          corner = rVec + thisOctal%subcellSize/2.d0 * VECTOR(-1.d0, -1.d0, 1.d0)
          r = modulus(corner - source(i)%position)
          if (r < source(i)%radius) then
             insideSource = .true.
          endif

       endif
    enddo
  end function insideSource



  end module source_mod
