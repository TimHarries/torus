module surface_mod

  use kind_mod
  use constants_mod
  use vector_mod
  use utils_mod
  use gridtype_mod
  use density_mod
  use romanova_class
  use messages_mod
  use spectrum_mod

  implicit none

  public :: &
       buildSphere, &
       emptySurface, &
       addElement, &
       testSurface, &
       createTTauriSurface, &
       createHotRing, &
       createSurface, &
       createTTauriSurface2, &
       sumSurface, &
       createProbs, &
       getPhotoVec, &
       photoFluxIntegral, &
       photoSolidAngle, &
       whichElement, &
       isHot

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
     integer :: ntheta
     real(double) :: surfaceMdot, trueMdot
     real :: radius ! 1.e10 cm
     integer, pointer :: angleArray(:,:)
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
    use input_variables, only : teff, lamstart, lamEnd, nLambda
    type(VECTOR),intent(in) :: centre
    real(double),intent(in) :: radius ! 1.e10 cm
    real(double) :: area ! 1.e20 cm^2
    type(SURFACETYPE),intent(out) :: surface
    integer,intent(in) :: nTheta
    character(len=*),intent(in) :: contfile
    integer :: nPhi, i, j, n
    real(double) :: theta, phi, dTheta, dPhi
    real(double) :: dPhase
    type(VECTOR) :: rVec
    integer :: nNuHotFlux
    logical :: ok
    type(SPECTRUMTYPE) :: hotspec

    ok = .true.

    surface%centre = centre
    surface%radius = radius
    surface%nElements = 0
    surface%nTheta = ntheta

    allocate(surface%angleArray(1:nTheta,1:nTheta))
    ! open the continuum flux file to get the number of points

    if (contfile /= "blackbody") then
       call readSpectrum(hotSpec, contfile, ok)
    else
       call fillSpectrumBB(hotSpec, dble(teff),  dble(lamstart), dble(lamEnd), nLambda)
    endif


    nNuHotFlux = hotSpec%nLambda
    allocate(surface%nuArray(nNuHotFlux))
    allocate(surface%hnuArray(nNuHotFlux))
    surface%nNuHotFlux = hotSpec%nLambda
    surface%nuArray = hotSpec%lambda
    surface%hnuArray = hotSpec%flux

!    open(20,file=contfile,status="old",form="formatted")
!    nNuHotFlux = 0
!    do
!       read(20,*,iostat=returnVal) dummyReal, dummyReal
!       if (returnVal /= 0) exit
!       nNuHotFlux = nNuHotFlux + 1
!    end do
!    close(unit=20)
!    ! store the value 
!    surface%nNuHotFlux = nNuHotFlux 
!    allocate(surface%nuArray(nNuHotFlux))
!    allocate(surface%hnuArray(nNuHotFlux))
!    ! now read in the array again, and store it.
!    open(20,file=contfile,status="old",form="formatted")
!    nNuHotFlux = 1
!    do
!       read(20,*,iostat=returnVal) nu, hnu
!       if (returnVal /= 0) exit
!       surface%nuArray(nNuHotFlux)=nu;  surface%hnuArray(nNuHotFlux)=hnu
!       nNuHotFlux = nNuHotFlux + 1
!    end do
!    close(unit=20)


    call convertToFnu(surface%nuarray, surface%hnuArray, surface%nNuHotFlux)

    n = 0
    do i = 1, nTheta
       theta = pi*real(i-1)/real(nTheta-1)
!       nphi = max(1,nint(real(nTheta)*sin(theta)))
       nPhi = ntheta
       do j = 1, nPhi
          n = n + 1
       enddo
    enddo
    if (Writeoutput) write(*,*) "Creating a sphere of ",n," elements..."
    if (associated(surface%element)) then
      print *, 'Trying to ALLOCATE surface%element in buildSphere, but it''s'
      print *, '  already ASSOCIATED! (size is ',SIZE(surface%element),')'         
      stop
    else
      allocate(surface%element(1:n))
    end if

    n = 0
    do i = 1, nTheta
       theta = pi*real(i-1)/real(nTheta-1)
!       nphi = max(1,nint(real(nTheta)*sin(theta)))
       nPhi = nTheta ! makes it easier to find element
       call random_number(dPhase)
       dPhase = dphase * twoPi
       do j = 1, nPhi
          n = n + 1
          if (nPhi > 1) then
             phi = twoPi * real(j-1)/real(nPhi) !+ dPhase
          else
             phi = 0.
          endif
          if (phi > twoPi) phi = phi - twoPi
          rVec = VECTOR(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta))
          dTheta = pi / real(nTheta)
          dPhi = twoPi / real(nPhi)
          area = radius * dTheta * radius * sin(theta) * dPhi
          call addElement(surface, radius, rVec, area)
          surface%angleArray(i,j) = n
       enddo
    enddo

    if (writeoutput) write(*,*) "done."
  end subroutine buildSphere


  integer function getElement(surface, direction)
    type(VECTOR) :: direction , td
    type(SURFACETYPE) :: surface
    real :: theta, phi
    integer :: itheta, iphi
!    td = (-1.d0) * direction

    td = direction
    theta = acos(td%z/modulus(td))
    phi = atan2(td%y,td%x)
    if (phi < 0.d0) phi = phi + twoPi

    itheta = real(surface%ntheta)*theta/pi + 1
    iphi =  real(surface%ntheta)*phi/twopi + 1

    iPhi = min(surface%nTheta, iphi)
    iTheta = min(surface%nTheta, iTheta)

    getElement = surface%angleArray(iTheta, iPhi)


  end function getElement

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
    real(double),intent(in) :: radius, area
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

      if (Writeoutput) then
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
      endif
    else
      if (Writeoutput) print *, 'Surface does not contain accretion hotspots'
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
    type(VECTOR) :: aboveSurface
    real(double) :: Laccretion, Taccretion
    real :: ttauriMdotLocal
    real :: theta1, theta2
    
    call writeInfo("Creating T Tauri stellar surface",TRIVIAL)
    
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
!        Taccretion = 7500.0
!! FOR DEBUG==================================================
        surface%element(iElement)%hotFlux(:) = &
           pi*blackbody(REAL(tAccretion), 1.e8*real(cSpeed)/surface%nuArray(:)) !* &
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
    
    print *, 'Creating a unirform stellar surface'
    
    do iElement = 1, SIZE(surface%element)
       surface%element(iElement)%hot = .false.
    end do
    call createProbs(surface,lineFreq,coreContFlux,fAccretion)
    call sumSurface(surface)
    
  end subroutine createSurface



  !
  !
  ! Creating the surface with hot spots...
  ! Similar to createTTauriSurface, but this routine 
  ! finds the mass-flux hence the kinetic energy flux (F_k) across 
  ! of the near the surface and its assumes its balanced with 
  ! the blackbody radiation (F_b).  Ignorging the internal energy 
  ! of the plasma.
  !        
  !  F_k = Vn (1/2 rho V^2)  and F_b = sigma*T^4
  !
  !             -->      ^       |-->|
  !   where Vn = V .dot. n,  V = | V |.  Sigma is Stefan-Boltzmann constant.  
  !
  !-----------------------------------------------------------------------
  ! This is for "romanova" geometry now, but can be generalized later. 
  ! 
  !
  !
  subroutine createTTauriSurface2(surface, grid, romData, lineFreq, &
       coreContFlux, fAccretion)
   
    type(SURFACETYPE),intent(inout) :: surface
    type(gridType), intent(in) :: grid
    type(romanova), intent(in) :: romData
    real(double), intent(in) :: coreContFlux
    real, intent(in) :: lineFreq
    real, intent(out) :: fAccretion ! erg s^-1 Hz^-1
    !
    integer :: iElement
    type(VECTOR) :: aboveSurface
    real(double) :: rho, Vn
    type(vector) :: V  ! velocity 
    type(vector) :: n  ! normal vector
    real(double) :: F_k, Vabs, T, Tmax, Tmin, Tsum, T_gas !, T_eff
    real(double) :: T_ref, V_ref
    real(double) :: units  ! unit conversion factor
    integer :: nhot
    ! Adiabatic exponent (It cannot be 1.0).
!    real(double), parameter :: gam = 3.0d0/4.0d0 ! ideal gas monoatomic)
    real(double), parameter :: gam = 1.1d0 ! Romanova's choice
    real(double) :: tmp, fac, w
    
    print *, ' '
    print *, 'Creating stellar surface elements by createTTauriSurface2.'

    ! First check if the data is alive
    if (.not. romanova_data_alive(romData)) then
       print*, "Error:: romData (romanova's data) is not alive anymore. [createTTauriSurface2]."
       stop
    end if

!    ! finding the effective temperature of the photosphere without hot spots first
!    T_eff = (real(SUM(surface%totalPhotosphere),kind=db) / &
!         (fourPi*(surface%radius*1.e10)**2.*stefanBoltz))**0.25
!    print *, "  Effective temperature of the photosphere is : ", T_eff, " [K]."
!    print *, "  The surface is will be flaged as hot surface T > T_eff ...."
    
    
    !
    nhot = 0       ! number of hot surface elemets
    Tmin = 1.0d10  ! 
    Tmax = 1.0d-10 !
    Tsum = 0.0d0
    do iElement = 1, SIZE(surface%element)
       aboveSurface = surface%element(iElement)%position - surface%centre
       aboveSurface = aboveSurface * 1.2_oc  ! location just above the surface element
       
       ! Finding the density and velocity at this position.
       ! -- using the routines in romanova_class
       rho = romanova_density(romData, aboveSurface)        ! [g/cm^3]
       V = romanova_velocity(romData, aboveSurface)         ! [c]
       T_gas = romanova_temperature(romData, aboveSurface)  ! [K]
       n = surface%element(iElement)%norm
       
       Vabs = MODULUS(V)
       Vn = DBLE( v .dot. n  )  
       ! Unit vonversions [c] --> [cm/s]
       Vn = Vn*cSpeed_dbl
       Vabs = Vabs*cSpeed_dbl

       ! finding the specific enthalpy (w) of gas note T=(P/rho) in Romanova's code
       ! units and w = gam*(gam-1)^-1 * (P/rho), so we can convert T to w here.
       ! Just need to be careful about the unit vonversion.
       T_ref = get_dble_parameter(romData, "T_ref")
       V_ref = get_dble_parameter(romData, "v_ref")
       units = V_ref*V_ref/T_ref
       fac = gam/(gam-1.0d0)    
       
       w = fac * T_gas * units   ! should be in cgs now.

       tmp  = 0.5d0*Vabs*Vabs
       tmp =  tmp + w

       F_k = ABS(rho*Vn*tmp)     ! matter energy flux [erg/s/cm^2]


       T = (F_k/DBLE(stefanBoltz))**0.25d0   ! effective temperature of the [K]

!       !
!       ! For debug ---------------
!       !  Should be eliminated later. 
!       T = 1.36d0*T
       

       ! For now, we define the surface to be "HOT" when T > 4000K. 
       ! This shold be changed to be read from parameter file later. 

       if (T > 4000.0d0) then 
          ! in futire this should be somehow compared to the effective
          ! temperature of photosphere. 
          !       if (T > T_eff) then
          nhot = nhot + 1
          surface%element(iElement)%hot = .true.
          allocate(surface%element(iElement)%hotFlux(surface%nNuHotFlux))
          
          surface%element(iElement)%hotFlux(:) = &
               pi*blackbody(REAL(T), 1.e8*REAL(cSpeed)/surface%nuArray(:)) !* &
          !                  ((1.e20*surface%element(iElement)%area)/(fourPi*TTauriRstar**2))
         surface%element(iElement)%temperature = T
         Tmin = MIN(Tmin, T)
         Tmax = MAX(Tmax, T)
         Tsum = Tsum + T
      else 
         surface%element(iElement)%hot = .false.
      end if
   end do

   if (nhot > 0) then
      write(*,*) "The minimum temperature of hot surface elements is ", Tmin
      write(*,*) "The maximum temperature of hot surface elements is ", Tmax
      write(*,*) "The average temperature of hot surface elements is ", Tsum/DBLE(nhot)
   end if

   call createProbs(surface,lineFreq,coreContFlux,fAccretion)
   call sumSurface(surface)
    
 end subroutine createTTauriSurface2





  subroutine sumSurface(surface, luminosity)
    type(SURFACETYPE),intent(inout) :: surface
    real(double), optional :: luminosity
    integer :: iElement
    integer :: iNu


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
    if (writeoutput) write(*,*) "SUM(surface%totalPhotosphere(:)", SUM(surface%totalPhotosphere(:))
    if (PRESENT(luminosity)) luminosity =  SUM(surface%totalPhotosphere(:))
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
    type(VECTOR),intent(out) :: position
    type(VECTOR),intent(out) :: direction
    real :: r
    integer :: j

    call random_number(r)
    call locate(surface%element(:)%prob, SIZE(surface%element), r, j)
    position = surface%centre + 1.0001d0*surface%element(j)%position
    direction = fromPhotosphereVector(surface%element(j)%norm)
  end subroutine getPhotoVec
  
  pure function photoFluxIntegral(position, surface, nFreqs) result(photoFlux)
    type(vector),intent(in) :: position
    type(surfaceType),intent(in) :: surface
    integer, intent(in) :: nFreqs
    real, dimension(nFreqs) :: photoFlux

    integer :: i
    type(vector) :: rVec
    real(double) :: dist, photoOmega
    real(double) :: cosTheta, geometricalFactor
    
    photoFlux = 0.0
    photoOmega = 0.0
    
    do i = 1, SIZE(surface%element)
       rVec = position - (surface%centre+surface%element(i)%position)
       dist = modulus(rVec)
       rVec = (-1.d0) * rVec / dist
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
    real(double) :: cosTheta
    real(double) :: dist
    hotOmega = 0.
    photoOmega = 0.

    do i = 1, SIZE(surface%element)
       rVec = position - (surface%centre+surface%element(i)%position)
       dist = modulus(rVec)
       rVec = (-1.d0) * rVec / dist
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



  !
  ! Given a surface object and a position vector,  
  ! this routine find the element which has
  ! the smallest angular distance from the position vector. 
  ! The function returns the index of the elements in the surface
  ! rather than the element itself.
  function whichElement(thisSurface, position) RESULT(out)
    implicit none
    integer :: out
    type(surfacetype), intent(in) :: thisSurface
    type(vector) :: position
    !
    integer :: i 
    real :: ip_max, ip
    real :: mp
    real :: imax
    
    
    ! now we go through all the surface elements.
    ! this is not a very smart thing to do....
    ! What you really need is a better data structure.
    ip_max = -1.0e20
    imax =1 
    mp = modulus(position)
    do i = 1, thisSurface%nElements
       ! finding the inner product of the surface position vector 
       ! and the input position vector.
       ip = position .dot. thisSurface%element(i)%position
       ! normalize it just in case..
       ip = ip / mp / MODULUS(thisSurface%element(i)%position)

       if (ip > ip_max) then
          ip_max = ip
          imax = i
       end if
    end do

     out = imax

  end function whichElement
  

  !
  ! Given a surface, and an index of the surface elemnt
  ! this function checks if the element is "hot" .
  logical function isHot(thisSurface, i)
    implicit none
    type(surfacetype), intent(in) :: thisSurface
    integer, intent(in) :: i
    isHot = thisSurface%element(i)%hot
  end function isHot
    


  SUBROUTINE createMagStreamSurface(surface, grid, lineFreq, &
                                 coreContFlux,fAccretion)

    USE magField                             
   
    TYPE(surfacetype),INTENT(INOUT) :: surface
    TYPE(gridType), INTENT(IN) :: grid
    REAL(double), INTENT(IN) :: coreContFlux
    REAL, INTENT(IN) :: lineFreq
    REAL, INTENT(OUT) :: fAccretion ! erg s^-1 Hz^-1
    INTEGER :: iElement, iSpot
    TYPE(VECTOR) :: surfacePoint
    
    PRINT *, 'Creating stellar surface for magStream system'
    write(*,*) "number of hotspots",size(magfieldhotspots), magfieldhotspots(1)%temperature
    
    DO iElement = 1, SIZE(surface%element)

      surface%element(iElement)%hot = .FALSE.
      
      ! find the distance between each element's centre and
      !   the "subsurface" position 
      
      surfacePoint = surface%element(iElement)%position - surface%centre

      DO iSpot = 1, SIZE(magFieldHotspots)
        IF ( modulus(magFieldHotspots(iSpot)%subSurfacePosn - surfacePoint)   &
             < magFieldHotspots(iSpot)%radius_1e10 ) THEN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF ( surface%element(iElement)%hot ) THEN
            ! if this element has already been set as "hot" because of another 
            !   hotspot, we use the hotest of the spots to compute the properties
            IF ( surface%element(iElement)%temperature < magFieldHotspots(iSpot)%temperature ) THEN
              surface%element(iElement)%temperature = magFieldHotspots(iSpot)%temperature

            !  ! add mass accretion rate from this surface element
            !  surface%surfaceMdot =  surface%surfaceMdot + &
            !     surface%element(iElement)*area * &
            !     magFieldHotspots(iSpot)%rho * &
            !     magFieldHotspots(iSpot)%speed / mSol / secsToYears
              
              surface%element(iElement)%hotFlux(:) = &
                pi*blackbody(surface%element(iElement)%temperature, 1.e8*real(cSpeed)/surface%nuArray(:))
            ELSE
              ! new spot is cooler than last one. do nothing.
            END IF

          ELSE

            surface%element(iElement)%hot = .TRUE.
            ALLOCATE(surface%element(iElement)%hotFlux(surface%nNuHotFlux))
            surface%element(iElement)%temperature = magFieldHotspots(iSpot)%temperature
            
            ! add mass accretion rate from this surface element
             surface%surfaceMdot =  surface%surfaceMdot + &
                 surface%element(iElement)%area * &
                 magFieldHotspots(iSpot)%rho * &
                 magFieldHotspots(iSpot)%speed / mSol / secsToYears

            surface%element(iElement)%hotFlux(:) = &
              pi*blackbody(magFieldHotspots(iSpot)%temperature, 1.e8*real(cSpeed)/surface%nuArray(:))
          END IF
        ELSE
          ! not near hot spot
        END IF
      END DO
          
    END DO
    PRINT *, '         stellar surface created, finishing...'
    
    ! calculate Mdot directly from hotspot variables
    surface%trueMdot = 0.0_db
    DO iSpot = 1, SIZE(magFieldHotspots), 1
      surface%trueMdot = surface%trueMdot + &
         (magFieldHotspots(iSpot)%radius_1e10 * 1.e10)**2 * pi * &
          magFieldHotspots(iSpot)%rho * &
          magFieldHotspots(iSpot)%speed
    END DO
    ! convert from g per second to Msol per year
    surface%trueMdot = surface%trueMdot / mSol 
    surface%trueMdot = surface%trueMdot / secsToYears

    CALL createProbs(surface,lineFreq,coreContFlux,fAccretion)
    CALL sumSurface(surface)
    PRINT *, '         stellar surface finished'
    
  END SUBROUTINE createMagStreamSurface






end module surface_mod
