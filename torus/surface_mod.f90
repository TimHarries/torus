module surface_mod

  use kind_mod
  use constants_mod
  use vector_mod
  use utils_mod
  use gridtype_mod
  use amr_mod

  implicit none

  public

  type ELEMENTTYPE
     type(VECTOR) :: norm
     real :: area ! 1.e20 cm^2
     type(VECTOR) :: position
     real, dimension(:), pointer :: hotFlux => null()
     real :: temperature
     real :: prob
     logical :: hot
  end type ELEMENTTYPE

  type SURFACETYPE
     integer :: nElements
     real :: radius ! 1.e10 cm
     type(ELEMENTTYPE),pointer :: element(:) => null()
     type(VECTOR) :: centre ! 1.e10 cm
     integer :: nNuHotFlux
     real, dimension(:), pointer :: nuArray  => null()
     real, dimension(:), pointer :: hnuArray => null()
     real(double), dimension(:), pointer :: totalPhotosphere => null()
     real(double), dimension(:), pointer :: totalAccretion => null()
  end type SURFACETYPE

contains

  subroutine buildSphere(centre, radius, surface, nTheta, contFile)
    type(VECTOR),intent(in) :: centre
    real,intent(in) :: radius ! 1.e10 cm
    real :: area ! 1.e20 cm^2
    type(SURFACETYPE),intent(out) :: surface
    integer,intent(in) :: nTheta
    character(len=*),intent(in) :: contfile
    integer :: nPhi, i, j, n
    real :: theta, phi, dTheta, dPhi
    real :: dPhase
    type(VECTOR) :: rVec
    integer :: nNuHotFlux
    integer :: returnVal
    real :: dummyReal

    surface%centre = centre
    surface%radius = radius
    surface%nElements = 0

    ! open the continuum flux file to get the number of points
    open(20,file=contfile,status="old",form="formatted")
    nNuHotFlux = 0
    do
       read(20,*,iostat=returnVal) dummyReal, dummyReal
       if (returnVal /= 0) exit
       nNuHotFlux = nNuHotFlux + 1
    end do
    close(unit=20)
    ! store the value 
    surface%nNuHotFlux = nNuHotFlux 
    allocate(surface%nuArray(nNuHotFlux))
    allocate(surface%hnuArray(nNuHotFlux))
    ! now read in the array again, and store it.
    open(20,file=contfile,status="old",form="formatted")
    nNuHotFlux = 1
    do
       read(20,*,iostat=returnVal) surface%nuArray(nNuHotFlux), surface%hnuArray(nNuHotFlux)
       if (returnVal /= 0) exit
       nNuHotFlux = nNuHotFlux + 1
    end do
    close(unit=20)
    
    n = 0
    do i = 1, nTheta
       theta = pi*real(i-1)/real(nTheta-1)
       nphi = max(1,nint(real(nTheta)*sin(theta)))
       do j = 1, nPhi
          n = n + 1
       enddo
    enddo
    write(*,*) "Creating a sphere of ",n," elements..."
    if (associated(surface%element)) then
      print *, 'Trying to ALLOCATE surface%element in buildSphere, but it''s'
      print *, '  already ASSOCIATED! (size is ',SIZE(surface%element),')'         
      stop
    else
      allocate(surface%element(1:n))
    end if

    do i = 1, nTheta
       theta = pi*real(i-1)/real(nTheta-1)
       nphi = max(1,nint(real(nTheta)*sin(theta)))
       call random_number(dPhase)
       dPhase = dPhase * twoPi
       do j = 1, nPhi
          n = n + 1
          if (nPhi > 1) then
             phi = twoPi * real(j-1)/real(nPhi-1) + dPhase
          else
             phi = 0.
          endif
          rVec = VECTOR(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta))

          dTheta = pi / real(nTheta)
          dPhi = twoPi / real(nPhi)
          area = radius * dTheta * radius * sin(theta) * dPhi
          call addElement(surface, radius, rVec, area)
       enddo
    enddo
    write(*,*) "done."
  end subroutine buildSphere

  subroutine emptySurface(surface)
    type(SURFACETYPE),intent(inout) :: surface
    integer :: i 

    if (associated(surface%element)) then
       do i = 1, surface%nelements
          if (associated(surface%element(i)%hotFlux)) then
             deallocate(surface%element(i)%hotFlux)
             nullify(surface%element(i)%hotFlux)
          end if
       end do
       deallocate(surface%element)
    end if
    if (associated(surface%nuarray)) deallocate(surface%nuarray)
    if (associated(surface%hnuarray)) deallocate(surface%hnuarray)

    surface%nElements = 0 
    nullify(surface%element) 
    nullify(surface%nuarray) 
    nullify(surface%hnuarray) 
    
    print *, 'Deallocated memory used for surface'

    
  end subroutine emptySurface
    
  subroutine addElement(surface, radius, rVec, area)
    real,intent(in) :: radius, area
    type(SURFACETYPE),intent(inout) :: surface
    type(VECTOR),intent(in) :: rVec

    surface%nElements = surface%nElements + 1
    surface%element(surface%nElements)%norm = rVec
    surface%element(surface%nElements)%position = radius * rVec
    surface%element(surface%nElements)%area = area
  end subroutine addElement

  subroutine testSurface(surface)
    type(SURFACETYPE),intent(in) :: surface
    
    real :: fillFactor
    real, dimension(:), allocatable :: meanTemp, area
    
    print *, ''
    print *, 'Stellar surface tests:'
    if (surface%nElements /= SIZE(surface%element)) then
      print *, 'In testSurface, surface%nElements /= SIZE(surface%element)'
      stop
    end if

    fillfactor = SUM(surface%element%area,MASK=surface%element%hot) / &
                 SUM(surface%element%area)
    
    write(*,*) "Surface area test (should be ~1): ", &
               SUM(surface%element(1:surface%nElements)%area)/(fourPi * surface%radius**2)

               
    write(*,'(a,e14.3)') "Photospheric luminosity (ergs/s): ",&
               SUM(surface%totalPhotosphere) 
    write(*,'(a,f14.3)') "Photospheric luminosity / solar luminosity: ",&
               SUM(surface%totalPhotosphere) / lSol

    write(*,'(a,f14.1)') "Photospheric black-body temperature estimated to be about: ",&
               (real(SUM(surface%totalPhotosphere),kind=db) / &
               (fourPi*(surface%radius*1.e10)**2.*stefanBoltz))**0.25
    
    if (any(surface%element%hot)) then
      print *, 'Surface contains accretion hotspots'
      allocate(meanTemp(SIZE(surface%element)))
      allocate(area(SIZE(surface%element)))
      
      where(surface%element%hot)
        meanTemp = surface%element%temperature
        area = surface%element%area
        meanTemp = meanTemp * area
      elsewhere
        meanTemp = 0.0
        area = 0.0
      end where

      write (*,'(a,f12.3)') 'Hotspot filling factor: ',fillFactor
                            
      write (*,'(a,f12.3)') 'Min temperature of accretion region: ',&
                            MINVAL(surface%element%temperature,MASK=surface%element%hot)
      write (*,'(a,f12.3)') 'Max temperature of accretion region: ',&
                            MAXVAL(surface%element%temperature,MASK=surface%element%hot)
      
      write (*,'(a,f12.3)') 'Mean temperature of accretion region: ',SUM(meanTemp) / SUM(area)
      
      write(*,'(a,e14.3)') "Accretion luminosity (ergs/s): ",&
                 SUM(surface%totalAccretion)
      write(*,'(a,f14.3)') "Accretion luminosity / solar luminosity: ",&
                 SUM(surface%totalAccretion) / lSol
      write(*,'(a,f14.3)') "Accretion luminosity / total photospheric luminosity: ",&
                 SUM(surface%totalAccretion) / SUM(surface%totalPhotosphere)
      write(*,'(a,f14.3)') "Accretion luminosity  / local photospheric luminosity (approximate!): ",&
                 SUM(surface%totalAccretion) / (SUM(surface%totalPhotosphere) * fillfactor)
      write(*,'(a,f14.1)') "Accretion black-body temperature estimated to be about: ",&
                 (real(SUM(surface%totalAccretion),kind=db) / &
                 ((fourPi*(surface%radius*1.e10)**2.*stefanBoltz)*fillFactor) )**0.25
    else
      print *, 'Surface does not contain accretion hotspots'
    end if
    print *, ''
    
  end subroutine testSurface

  subroutine createTTauriSurface(surface, grid, lineFreq, &
                                 coreContFlux,fAccretion)

    use input_variables, only: TTauriRinner, TTauriRouter,&
                               TTauriRstar, TTauriMstar
   
    type(SURFACETYPE),intent(inout) :: surface
    type(gridType), intent(in) :: grid
    real(double), intent(in) :: coreContFlux
    real, intent(in) :: lineFreq
    real, intent(out) :: fAccretion ! erg s^-1 Hz^-1
    integer :: iElement
    type(octalVector) :: aboveSurface
    real(double) :: Laccretion, Taccretion
    real :: ttauriMdotLocal
    real :: theta1, theta2
    
    print *, 'Creating T Tauri stellar surface'
    
    do iElement = 1, SIZE(surface%element)
      aboveSurface = surface%element(iElement)%position - surface%centre
      aboveSurface = aboveSurface * 1.001_oc 

      TTauriMdotLocal = TTauriVariableMdot(aboveSurface,grid) ! g s^-1
      if (TTauriMdotLocal > 1.e15) then
        surface%element(iElement)%hot = .true.
        allocate(surface%element(iElement)%hotFlux(surface%nNuHotFlux))
        
        Laccretion = (REAL(bigG,KIND=db)* REAL(TTauriMstar,KIND=db)* &
          REAL(TTauriMdotLocal,KIND=db)/ REAL(TTauriRstar,kind=db))* &
          REAL((1.0_db-(2.0_db*TTauriRstar/(TTauriRouter+TTauriRinner))),KIND=db)
     
        ! to calculate the accretion temperature, we assume that the MdotLocal
        !   applies *everywhere* in the accretion impact rings.

        theta1 = asin(sqrt(TTauriRstar/TTauriRouter))
        theta2 = asin(sqrt(TTauriRstar/TTauriRinner))
        
        Taccretion = Laccretion / (fourPi * TTauriRstar**2 * stefanBoltz) 
        Taccretion = Taccretion / (abs(cos(theta1)-cos(theta2)))
        Taccretion = Taccretion**0.25
!! FOR DEBUG==================================================
!        Taccretion = 8000.0
!! FOR DEBUG==================================================
        surface%element(iElement)%hotFlux(:) = &
           pi*blackbody(REAL(tAccretion), 1.e8*cSpeed/surface%nuArray(:)) !* &
!                  ((1.e20*surface%element(iElement)%area)/(fourPi*TTauriRstar**2))
        surface%element(iElement)%temperature = Taccretion
      else 
        surface%element(iElement)%hot = .false.
      end if
    end do
    call createProbs(surface,lineFreq,coreContFlux,fAccretion)
    call sumSurface(surface)
    
  end subroutine createTTauriSurface
  
  subroutine createHotRing(surface, phiStart1, phiEnd1, phiStart2, phiEnd2, theta1, theta2, fluxRing, fluxPhoto)
    type(SURFACETYPE),intent(inout) :: surface
    real,intent(in) :: phiStart1, phiStart2, phiEnd1, phiEnd2
    real,intent(in) :: theta1, theta2
    real :: cosTheta1, cosTheta2
    real,intent(in) :: fluxRing, fluxPhoto
    logical :: inRing
    integer :: i
    real :: cosTheta, phi

    print *, 'createHotRing subroutine needs updated before you can use it!'
    stop
    
    cosTheta1 = cos(theta1)
    cosTheta2 = cos(theta2)

    do i = 1, SIZE(surface%element)
       inRing = .false.
       cosTheta = abs((surface%element(i)%position%z)/(modulus(surface%element(i)%position)))
       phi = atan2(surface%element(i)%position%y, surface%element(i)%position%x)
       if (phi < 0.) phi = phi+twoPi
       if ((cosTheta > cosTheta2).and.(cosTheta < cosTheta1)) then
          if ( ((phi > phiStart1).and.(phi < phiEnd1)).or.((phi > phiStart2).and.(phi < phiEnd2))) then
             inRing = .true.
          endif
       endif
       if (inRing) then
          !surface%element(i)%flux = fluxRing
          surface%element(i)%hot = .true.
       else
          !surface%element(i)%flux = fluxPhoto
          surface%element(i)%hot = .false.
       endif

    enddo
    !call createProbs(surface,lineFreq,coreContFlux)
  end subroutine createHotRing

  !
  !
  ! For a uniform suraface (for a general use)
  !
  subroutine createSurface(surface, grid, lineFreq, coreContFlux, fAccretion)

    type(SURFACETYPE),intent(inout) :: surface
    type(gridType), intent(in) :: grid
    real(double), intent(in) :: coreContFlux
    real, intent(in) :: lineFreq
    real, intent(out) :: fAccretion ! erg s^-1 Hz^-1
    integer :: iElement
    type(octalVector) :: aboveSurface
    
    print *, 'Creating a unirform stellar surface'
    
    do iElement = 1, SIZE(surface%element)
       surface%element(iElement)%hot = .false.
    end do
    call createProbs(surface,lineFreq,coreContFlux,fAccretion)
    call sumSurface(surface)
    
  end subroutine createSurface



  subroutine sumSurface(surface)
    type(SURFACETYPE),intent(inout) :: surface
    integer :: iElement
    real :: surfaceArea
    integer :: iNu

    surfaceArea = fourPi * surface%radius**2. ! 1.e20 cm^2

    if (.not. associated(surface%totalPhotosphere)) &
      allocate(surface%totalPhotosphere(surface%nNuHotFlux))    
    if (.not. associated(surface%totalAccretion)) &
      allocate(surface%totalAccretion(surface%nNuHotFlux))    
      
    do iNu = 1, SIZE(surface%hnuArray)-1, 1
       surface%totalPhotosphere(iNu) = 0.5*ABS((surface%hnuArray(iNu+1)+surface%hnuArray(iNu)) &
                                       * (surface%nuarray(iNu+1)-surface%nuarray(iNu))) &
                                       * SUM(surface%element%area)*1.e20
    enddo
    surface%totalPhotosphere(SIZE(surface%hnuArray)) = 0.0
    ! just checking ..
    write(*,*) "SUM(surface%totalPhotosphere(:)", SUM(surface%totalPhotosphere(:))

    surface%totalAccretion = 0.0
    do iElement = 1, SIZE(surface%element), 1
      if (surface%element(iElement)%hot) then
        do iNu = 1, SIZE(surface%hnuArray)-1, 1
          surface%totalAccretion(iNu) = surface%totalAccretion(iNu) + ABS( &
            0.5*SUM(surface%element(iElement)%hotFlux(iNu:iNu+1)) &
            * (surface%nuarray(iNu+1)-surface%nuarray(iNu)) &
            * surface%element(iElement)%area*1.e20)
        enddo
      end if 
    end do

  end subroutine sumSurface
  
  subroutine createProbs(surface,lineFreq,coreContFlux,fAccretion)
    type(SURFACETYPE),intent(inout) :: surface
    real(double), intent(in) :: coreContFlux ! erg s^-1 cm^-2 Hz^-1
    real, intent(in) :: lineFreq
    real, intent(out) :: fAccretion ! 1.e-20 erg s^-1 Hz^-1
    integer :: i, lineIndex
    real :: accretionFlux ! erg s^-1 cm^-2 Hz^-1

    fAccretion = 0.0
    
    ! find the line frequency in the tabulated fluxes
    call locate(surface%nuArray,SIZE(surface%nuArray),lineFreq,lineIndex)
    
    do i = 1, SIZE(surface%element)
      if (surface%element(i)%hot) then
        !accretionFlux = logint(lineFreq,surface%nuArray(lineIndex),&
        !                       surface%nuArray(lineIndex+1),&
        !                       surface%element(i)%hotFlux(lineIndex),&
        !                       surface%element(i)%hotFlux(lineIndex+1))
        accretionFlux = logint(lineFreq,surface%nuArray(lineIndex),&
                               surface%nuArray(lineIndex+1),&
                               surface%element(i)%hotFlux(lineIndex),&
                               surface%element(i)%hotFlux(lineIndex+1))
        fAccretion = fAccretion + (accretionFlux * surface%element(i)%area)                 
        surface%element(i)%prob = (accretionFlux + coreContFlux) * surface%element(i)%area
      else
        surface%element(i)%prob = coreContFlux * surface%element(i)%area
      end if
    end do
    !! have to average fAccretion over stellar surface to get it in same units as coreContFlux 
    !fAccretion = fAccretion / SUM(surface%element(1:surface%nElements)%area)
    
    do i = 2, SIZE(surface%element)
       surface%element(i)%prob = surface%element(i)%prob + surface%element(i-1)%prob
    enddo
    surface%element(:)%prob = surface%element(:)%prob - surface%element(1)%prob
    surface%element(:)%prob = surface%element(:)%prob / surface%element(SIZE(surface%element))%prob
  end subroutine createProbs
    
  subroutine getPhotoVec(surface, position, direction)
    type(SURFACETYPE),intent(in) :: surface
    type(OCTALVECTOR),intent(out) :: position
    type(OCTALVECTOR),intent(out) :: direction
    real :: r
    integer :: j

    call random_number(r)
    call locate(surface%element(:)%prob, SIZE(surface%element), r, j)
    position = surface%centre + 1.0001*surface%element(j)%position
    direction = fromPhotosphereVector(surface%element(j)%norm)
  end subroutine getPhotoVec
  
  pure function photoFluxIntegral(position, surface, nFreqs) result(photoFlux)
    type(vector),intent(in) :: position
    type(surfaceType),intent(in) :: surface
    integer, intent(in) :: nFreqs
    real, dimension(nFreqs) :: photoFlux

    integer :: i
    type(vector) :: rVec
    real :: dist, photoOmega
    real :: cosTheta, geometricalFactor
    
    photoFlux = 0.0
    photoOmega = 0.0
    
    do i = 1, SIZE(surface%element)
       rVec = position - (surface%centre+surface%element(i)%position)
       dist = modulus(rVec)
       rVec = (-1.) * rVec / dist
       cosTheta = rVec .dot. surface%element(i)%norm
       if (cosTheta < 0.) then
          geometricalFactor = abs(cosTheta)*surface%element(i)%area / dist**2
          photoOmega = photoOmega + geometricalFactor
          if (surface%element(i)%hot) then
             photoFlux = photoFlux + &
                         geometricalFactor * &
                         (surface%element(i)%hotFlux + surface%hNuArray)
          else
             photoFlux = photoFlux + geometricalFactor * surface%hNuArray
          endif
       endif
    enddo
    if (photoOmega>0) then
       photoFlux = photoFlux / photoOmega
    else
       ! the point must be below the surafce!!
       photoFlux = 1.0e-20
    end if
  
  end function photoFluxIntegral

  pure subroutine photoSolidAngle(position, surface, hotOmega, photoOmega)
    type(SURFACETYPE),intent(in) :: surface
    type(VECTOR),intent(in) :: position
    type(VECTOR) :: rVec
    real,intent(out) :: hotOmega, photoOmega
    integer :: i
    real :: cosTheta
    real :: dist
    hotOmega = 0.
    photoOmega = 0.

    do i = 1, SIZE(surface%element)
       rVec = position - (surface%centre+surface%element(i)%position)
       dist = modulus(rVec)
       rVec = (-1.) * rVec / dist
       cosTheta = rVec .dot. surface%element(i)%norm
       if (cosTheta < 0.) then
          if (surface%element(i)%hot) then
             hotOmega = hotOmega + abs(cosTheta)*surface%element(i)%area / dist**2
          else
             photoOmega = photoOmega + abs(cosTheta)*surface%element(i)%area / dist**2
          endif
       endif
    enddo
  end subroutine photoSolidAngle


end module surface_mod
