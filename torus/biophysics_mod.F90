module biophysics_mod

  use vector_mod
  use kind_mod
  use unix_mod
  use utils_mod
  use gridtype_mod
  use spectrum_mod

  implicit none

  type SOURCETYPE
     type(VECTOR) :: position
     type(VECTOR) :: direction
     real(double) :: power
     type(SPECTRUMTYPE) :: spectrum
     integer :: iType ! 1 fibre, 2 laser
     real(double) :: radius
     real(double) :: wavelength
     real(double) :: numericalAperture
     real(double) :: xSize, ySize ! cm
  end type SOURCETYPE

  type PHOTONSOURCE
     type(VECTOR) :: position
     type(VECTOR) :: direction
     integer :: nComponent
     type(SOURCETYPE) :: component(10)
     real(double) :: totalPower
  end type PHOTONSOURCE

     

  type TISSUETYPE
     character(len=80) :: name
     integer :: nlambda
     real(double), pointer :: lamArray(:) => null()
     real(double), pointer :: muAbs(:) => null()
     real(double), pointer :: muSca(:) => null()
     real(double), pointer :: gFac(:) => null()
     real(double), pointer :: n(:) => null()
  end type TISSUETYPE

  integer :: nTissue
  type(TISSUETYPE), pointer :: tissueArray(:)
  integer :: nh2o, nHb
  real(double), pointer :: lamH2O(:), muH2O(:)
  real(double), pointer :: lamHb(:), muHb(:), muHbo2(:)

contains

  subroutine createSource(source)
    use inputs_mod, only : sourcePosition, sourceTheta, sourcePhi, nComponent
    use inputs_mod, only : componentPower, componentPosition, componentRadius
    use inputs_mod, only : componentNa, componentWavelength, componentType
    type(VECTOR) :: newPos
    type(PHOTONSOURCE) :: source
    integer :: i

    source%position = sourcePosition
    source%direction = VECTOR(sin(sourceTheta)*cos(sourcePhi),sin(sourceTheta)*sin(sourcePhi), cos(sourceTheta))
    source%nComponent = nComponent
    do i = 1, nComponent

       source%component(i)%direction = source%direction
       newPos = componentPosition(i)
       newPos = rotateY(newPos, sourceTheta)
       newPos = rotateZ(newPos, sourcePhi)
       newPos = newPos
       select case(componentType(i))
          case("fibre")
             source%component(i)%iType = 1
             source%component(i)%power = componentPower(i)
             source%component(i)%position = newPos
             source%component(i)%radius = componentRadius(i)
             source%component(i)%numericalAperture = componentNA(i)
             source%component(i)%wavelength = componentwavelength(i)
          case("laser")
             source%component(i)%iType = 2
             source%component(i)%power = componentPower(i)
             source%component(i)%radius = componentRadius(i)
             source%component(i)%wavelength = componentwavelength(i)
          case DEFAULT
             call writeFatal("Cant understand component type: "//trim(componentType(i)))
             stop
       end select

    enddo
    source%totalPower = SUM(source%component(1:nComponent)%power)
  end subroutine createSource
  
  subroutine emitPhoton(source, photonPosition, photonDirection, photonWavelength, photonWeight, photonNormal)
    type(PHOTONSOURCE) :: source
    type(VECTOR) :: photonPosition, photonDirection, photonNormal
    real(double) :: photonWavelength, photonWeight
    real(double), allocatable, save :: prob(:)
    integer :: i, iComponent
    real(double) :: phi, r, thisR
    type(VECTOR) :: normVec
    type(VECTOR) :: xAxis, vec
    xAxis = VECTOR(1.d0, 0.d0, 0.d0)

    if (source%nComponent > 1) then
       if (.not.allocated(prob)) then
          
          allocate(prob(1:source%nComponent))
          prob(1:source%ncomponent) = source%component(1:source%ncomponent)%power
          prob(1:source%ncomponent) = prob(1:source%ncomponent) / source%totalPower
          if (source%nComponent > 2) then
             do i = 2, source%nComponent
                prob(i) = prob(i) + prob(i-1)
             enddo
             prob(1:source%nComponent) = prob(1:source%nComponent) / prob(source%nComponent)
          endif
       endif
       
       if (source%nComponent > 2) then
          call randomNumberGenerator(getDouble=r)
          if (r < prob(1)) then
             iComponent = 1
          else
             call locate(prob, source%nComponent, r, iComponent)
             iComponent = iComponent + 1
          endif
       else
          call randomNumberGenerator(getDouble=r)
          if (r < prob(1)) then
             iComponent = 1
          else 
             iComponent = 2
          endif
       endif
    else
       iComponent = 1
    endif

    photonPosition = source%position + source%component(iComponent)%position
    call randomNumberGenerator(getDouble=r)
    thisR = source%component(iComponent)%radius * sqrt(r)
    call randomNumberGenerator(getDouble=phi)
    phi = phi * twoPi

    normVec = source%component(iComponent)%direction .cross. xAxis
    photonNormal = normVec
    normVec = normVec * thisR
    vec = source%component(iComponent)%direction
    normVec = arbitraryRotate(normVec, phi, vec)
    photonPosition = photonPosition + normVec
    photonDirection = source%component(iComponent)%direction
    photonWavelength = source%component(iComponent)%wavelength
    photonWeight  = source%component(iComponent)%power/source%totalPower

    if (source%component(iComponent)%iType == 1) then
       call randomNumberGenerator(getDouble=r)
       phi = asin(source%component(iComponent)%numericalAperture)
       phi = acos((1.d0-cos(phi))*r+cos(phi))
       photonDirection = arbitraryRotate(photonDirection,phi,photonNormal)
       call randomNumberGenerator(getDouble=r)
       phi = twoPi * r
       photonDirection = arbitraryRotate(photonDirection,phi,source%component(iComponent)%direction)
    endif

  end subroutine emitPhoton





  subroutine readWater(nlam, lam, ua)
    integer :: nLam
    real(double), pointer :: lam(:),ua(:)
    character(len=120) :: dataDirectory, charLine, filename
    integer :: i
    call unixGetenv("TORUS_DATA", dataDirectory, i)
    filename = trim(dataDirectory)//"/"//"buiteveld1994.dat"
    
    
    open(21, file=filename,status="old", form="formatted")
    nLam = 0
10  continue
    read(21,*,end=20) charLine
    if (charLine(1:1)=="#") goto 10
    nLam = nLam + 1
    goto 10
20 continue
    close(21)
    allocate(lam(nLam), ua(nLam))
    open(21, file=filename,status="old", form="formatted")
    nLam = 0
30  continue
    read(21,'(a)',end=40) charLine
    if (charLine(1:1)=="#") goto 30
    nLam = nLam + 1
    read(charLine,*) lam(nLam), ua(nLam)
    goto 30
40 continue
    close(21)
  end subroutine readWater

  subroutine readHaemoglobin(nlam, lam, mua_hbo2, mua_hb)
    integer :: nLam
    real(double), pointer :: lam(:),mua_hbo2(:),mua_hb(:)
    character(len=120) :: dataDirectory, charLine, filename
    integer :: i

    call unixGetenv("TORUS_DATA", dataDirectory, i)
    filename = trim(dataDirectory)//"/"//"prahl1999.dat"
    
    
    open(21, file=filename,status="old", form="formatted")
    nLam = 0
10  continue
    read(21,*,end=20) charLine
    if (charLine(1:1)=="#") goto 10
    nLam = nLam + 1
    goto 10
20 continue
    close(21)
    allocate(lam(nLam), mua_hbo2(nLam), mua_hb(nLam))
    open(21, file=filename,status="old", form="formatted")
    nLam = 0
30  continue
    read(21,'(a)',end=40) charLine
    if (charLine(1:1)=="#") goto 30
    nLam = nLam + 1
    read(charLine,*) lam(nLam), mua_hbo2(nLam), mua_hb(nLam)
    goto 30
40 continue
    close(21)

    mua_hbo2 = mua_hbo2 * 2.303d0 * 150.d0 / 64500.d0
    mua_hb   = mua_hbo2 * 2.303d0 * 150.d0 / 64500.d0

  end subroutine readHaemoglobin
    


  subroutine setupMatcherSkinTissue()
    
    real(double) ::  CH2O
    real(double), allocatable :: lamArray(:)
    real(double) :: cMelanin
    integer :: nLam, i
    
    call readWater(nH2O, lamH2O, muH2O)
    call readHaemoglobin(nHb, lamHb, muhbo2, muhb)
    if (associated(tissueArray)) deallocate(tissueArray)
    nTissue = 7
    allocate(tissueArray(0:nTissue))

    nLam = 100
    allocate(lamArray(nLam))
    do i = 1, nLam
       lamArray(i) = 400.d0 + 400.d0 * dble(i-1)/dble(nLam-1)
    enddo

! air

    call createNewTissue(tissueArray(0), "Air", nLam, lamArray)
    tissueArray(0)%muAbs = 0.d0
    tissueArray(0)%muSca = 0.d0
    tissueArray(0)%gFac = 0.d0
    tissueArray(0)%n = 1.d0


! stratum corneum

    call createNewTissue(tissueArray(1), "Stratum corneum", nLam, lamArray)
    CH2O = 0.05d0
    do i = 1, tissueArray(1)%nLambda
       tissueArray(1)%muAbs(i) = ( (1.d0 - 0.3d-3 * tissueArray(1)%lamArray(i)) + &
            0.125d0 * mu_a_0(tissueArray(1)%lamArray(i))) & 
            * (1.d0 - CH2O) + CH2O * interpArray(tissueArray(1)%lamArray(i), nh2o, lamh2o, muh2o)
    enddo
    tissueArray(1)%muSca = 1000.d0
    tissueArray(1)%gFac = 0.86
    tissueArray(1)%n = 1.5d0

! living epidermis

    call createNewTissue(tissueArray(2), "Living epidermis", nLam, lamArray)
    do i = 1, tissueArray(2)%nLambda

       cMelanin = 0.5!!!!!!!!!!

       tissueArray(2)%muAbs(i) = (Cmelanin * muMelanin(tissueArray(2)%lamArray(i)) + &
            (1.d0 - Cmelanin)*mu_a_0(tissueArray(2)%lamArray(i)))*(1.d0-CH2O) + &
            CH2O * interpArray(tissueArray(1)%lamArray(i), nh2o, lamh2o, muh2o)
    enddo
    tissueArray(2)%muSca = 450.d0
    tissueArray(2)%gFac = 0.8
    tissueArray(2)%n = 1.34d0

! Papillary dermis

    call createNewTissue(tissueArray(3), "Papillary dermis", nLam, lamArray)
    call fillTissueMatcher(tissueArray(3), 0.04d0, 0.5d0, 300.d0, 0.9d0, 1.4d0)

! Upper blood net dermis

    call createNewTissue(tissueArray(4), "Upper blood net dermis", nLam, lamArray)
    call fillTissueMatcher(tissueArray(4), 0.3d0, 0.6d0, 350.d0, 0.95d0, 1.39d0)

! Reticular dermis

    call createNewTissue(tissueArray(5), "Reticular dermis", nLam, lamArray)
    call fillTissueMatcher(tissueArray(5), 0.04d0, 0.7d0, 250.d0, 0.8d0, 1.4d0)

! Deep blood net dermis

    call createNewTissue(tissueArray(6), "Deep blood net dermis", nLam, lamArray)
    call fillTissueMatcher(tissueArray(6), 0.1d0, 0.7d0, 300.d0, 0.95d0, 1.38d0)

! Subcutaneous fat

    call createNewTissue(tissueArray(7), "Subcutaneous fat", nLam, lamArray)
    call fillTissueMatcher(tissueArray(7), 0.05d0, 0.7d0, 50.d0, 0.75d0, 1.44d0)


  end subroutine setupMatcherSkinTissue

  subroutine fillTissueMatcher(tissue, cBlood, cH2O, muSca, g, n)
    type(TISSUETYPE) :: tissue
    real(double) :: cBlood, cH2O, muSca, g, n
    real(double), parameter :: FRBC = 0.99d0
    real(double), parameter :: FHb = 0.25d0
    real(double), parameter :: S = 0.6d0
    real(double), parameter :: Ht = 0.45d0
    real(double) :: gamma
    integer :: i
    gamma = FHb * FRBC * Ht

    do i = 1, tissue%nLambda

       tissue%muAbs(i) = (1.d0-S) * gamma * cBlood * interpArray(tissue%lamArray(i), nHb, lamHb, muHb) + &
            S * gamma * interpArray(tissue%lamArray(i), nHb, lamHb, muHbO2) + &
            (1.d0 - gamma*cBlood) * CH2O * interpArray(tissue%lamArray(i), nH2o, lamH2o, muH2o) + &
            (1.d0 - gamma*cBlood) * (1- CH2O) * mu_a_0(tissue%lamArray(i))
    enddo
    tissue%muSca = muSca
    tissue%gFac = g
    tissue%n = n
  end subroutine fillTissueMatcher

  subroutine createNewTissue(tissue, nameString, nLam, lamArray)
    type(TISSUETYPE) :: tissue
    character(len=*) :: nameString
    integer :: nLam
    real(double) :: lamArray(:)

    tissue%name = trim(nameString)
    tissue%nLambda = nLam
    allocate(tissue%lamArray(1:nLam))
    allocate(tissue%muSca(1:nLam))
    allocate(tissue%muAbs(1:nLam))
    allocate(tissue%gFac(1:nLam))
    allocate(tissue%n(1:nLam))
    tissue%lamArray = lamArray
  end subroutine createNewTissue

  real(double) function mu_a_0 (lambda)
    real(double) :: lambda
    mu_a_0 = 7.84d8 * lambda**(-3.255d0) ! equation 10
  end function mu_a_0

  real(double) function muMelanin (lambda)
    real(double) :: lambda
    muMelanin = 5.d10 * lambda**(-3.33d0)
  end function muMelanin

  real(double)  function interpArray(x, nx, xArray, yArray) result (out)
    real(double) :: x, xArray(:), yArray(:)
    integer :: nx, i

    call locate(xArray, nx, x, i)
    out = yArray(i) + (yArray(i+1)-yArray(i)) * (x - xArray(i))/(xArray(i+1)-xArray(i))
  end function interpArray

  real(double) function skinInterface(x, y, Zk, Akx, Aky, omegakx, omegaky, phikx, phiky) result (z)
    real(double) :: x, y, Akx, Aky, omegakx, omegaky, phikx, phiky, zk
    z = zk + (Akx * sin (omegakx * x + phikx)) + (Aky * sin (omegaky * y + phiky))
!    z = zk
  end function skinInterface

  logical function splitSkin(thisOctal, subcell)
    type(OCTAL) :: thisOctal
    integer :: subcell
    type(VECTOR) :: cen, probe(8)
    real(double) :: x, y, z, d
    real(double) :: zk(8)
    real(double) :: Akx(8)
    real(double) :: Aky(8)
    real(double) :: omegakx(8)
    real(double) :: omegaky(8)
    real(double) :: phikx(8)
    real(double) :: phiky(8)
    real(double) :: zLayer
    integer :: iLayer, iCorner
    
    
    splitSkin = .false.

    Akx = (/2.d0, 2.5d0, 20d0, 2.d0, 2.d0, 2.d0, 5.d0, 5.d0/)
    Aky = Akx

    omegakx = (/100.d0, 80.d0, 50.d0, 20.d0, 20.d0, 20.d0, 20.d0, 25.d0/)
    omegaky = (/150.d0, 80.d0, 45.d0, 40.d0, 50.d0, 50.d0, 50.d0, 30.d0/)
    omegakx = omegakx/1.d4
    omegaky = omegaky/1.d4

    omegakx = pi/omegakx
    omegaky = pi/omegaky
    
    zk = (/0.d0, 20.d0, 100.d0, 200.d0, 280.d0, 1900.d0, 2100.d0, 8000.d0/)

    Akx = Akx / 1.d4
    Aky = Aky / 1.d4

    zk = zk / 1.d4

    phikx = 0.d0
    phiky = 0.d0

    d = thisOctal%subcellSize / 2.d0

    probe(1) = vector(+1.d0, +1.d0, 1.0)
    probe(2) = vector(+1.d0, -1.d0, 1.0)
    probe(3) = vector(-1.d0, +1.d0, 1.0)
    probe(4) = vector(-1.d0, -1.d0, 1.0)
    probe(5) = vector(+1.d0, +1.d0,-1.0)
    probe(6) = vector(+1.d0, -1.d0,-1.0)
    probe(7) = vector(-1.d0, +1.d0,-1.0)
    probe(8) = vector(-1.d0, -1.d0,-1.0)
    
    do iCorner = 1, 8
       cen = subcellCentre(thisOctal,subcell) + probe(iCorner)*d
       x = cen%x
       z = cen%z
       y = cen%y
    

       do iLayer = 1, 7
          zLayer = skinInterface(x, y, Zk(iLayer), Akx(iLayer), Aky(iLayer), omegakx(iLayer), &
               omegaky(iLayer), phikx(iLayer), phiky(iLayer))
          if (abs(zLayer-z)< d) splitSkin =.true.
       enddo
    enddo
  end function splitSkin


  recursive subroutine fillMatcherSkinLayers(thisOctal)
    use amr_utils_mod, only : inSubcell
    use inputs_mod, only : smallestCellsize
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, iLayer
    real(double) :: x, y, z, d
    real(double) :: zk(8)
    real(double) :: Akx(8)
    real(double) :: Aky(8)
    real(double) :: omegakx(8)
    real(double) :: omegaky(8)
    real(double) :: phikx(8)
    real(double) :: phiky(8)
    real(double) :: zUpper, zLower, zLayer
    type(VECTOR) :: cen

    Akx = (/2.d0, 2.5d0, 20d0, 2.d0, 2.d0, 2.d0, 5.d0, 5.d0/)
    Aky = Akx

    omegakx = (/100.d0, 80.d0, 50.d0, 20.d0, 20.d0, 20.d0, 20.d0, 25.d0/)
    omegaky = (/150.d0, 80.d0, 45.d0, 40.d0, 50.d0, 50.d0, 50.d0, 30.d0/)
    omegakx = omegakx/1.d4
    omegaky = omegaky/1.d4

    omegakx = pi/omegakx
    omegaky = pi/omegaky
    
    zk = (/0.d0, 20.d0, 100.d0, 200.d0, 280.d0, 1900.d0, 2100.d0, 8000.d0/)

    Akx = Akx / 1.d4
    Aky = Aky / 1.d4

    zk = zk / 1.d4

    phikx = 0.d0
    phiky = 0.d0
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fillMatcherSkinLayers(child)
                exit
             end if
          end do
       else
          cen = subcellCentre(thisOctal, subcell)
          thisOctal%surfaceNormal(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
          x = cen%x
          y = cen%y
          z = cen%z
          d = thisOctal%subcellSize/2.d0 + 0.01d0*smallestCellSize
          if (.not.associated(thisOctal%iTissue)) allocate(thisOctal%iTissue(1:thisOctal%maxChildren))

          do iLayer = 1, 7
             zLower = skinInterface(x, y, Zk(iLayer), Akx(iLayer), Aky(iLayer), omegakx(iLayer), &
                  omegaky(iLayer), phikx(iLayer), phiky(iLayer))
             zUpper = skinInterface(x, y, Zk(iLayer+1), Akx(iLayer+1), Aky(iLayer+1), omegakx(iLayer+1), &
                  omegaky(iLayer+1), phikx(iLayer+1), phiky(iLayer+1))


             if ((z < zUpper).and.(z > zLower)) then
                thisOctal%rho(subcell) = real(iLayer)
                thisOctal%iTissue(subcell) = iLayer
             endif
             zLayer = skinInterface(x, y, Zk(iLayer), Akx(iLayer), Aky(iLayer), omegakx(iLayer), &
                  omegaky(iLayer), phikx(iLayer), phiky(iLayer))

             if (inSubcell(thisOctal,subcell,VECTOR(x,y,zlayer))) then
                thisOctal%iTissue(subcell) = iLayer
             endif

          enddo

      endif
    enddo
  end subroutine fillMatcherSkinLayers

  recursive subroutine setSurfaceNormals(grid,thisOctal)
    use inputs_mod, only : smallestCellsize
    use amr_utils_mod
    type(octal), pointer   :: thisOctal, neighbourOctal
    type(octal), pointer  :: child 
    type(GRIDTYPE) :: grid
    integer :: subcell, i, iLayer, neighbourSubcell
    real(double) :: x, y, z, d
    real(double) :: zk(8)
    real(double) :: Akx(8)
    real(double) :: Aky(8)
    real(double) :: omegakx(8)
    real(double) :: omegaky(8)
    real(double) :: phikx(8)
    real(double) :: phiky(8)
    real(double) :: zLayer, xp, yp, zp
    
    type(VECTOR) :: cen, vec

    Akx = (/2.d0, 2.5d0, 20d0, 2.d0, 2.d0, 2.d0, 5.d0, 5.d0/)
    Aky = Akx

    omegakx = (/100.d0, 80.d0, 50.d0, 20.d0, 20.d0, 20.d0, 20.d0, 25.d0/)
    omegaky = (/150.d0, 80.d0, 45.d0, 40.d0, 50.d0, 50.d0, 50.d0, 30.d0/)
    omegakx = omegakx/1.d4
    omegaky = omegaky/1.d4

    omegakx = pi/omegakx
    omegaky = pi/omegaky
    
    zk = (/0.d0, 20.d0, 100.d0, 200.d0, 280.d0, 1900.d0, 2100.d0, 8000.d0/)

    Akx = Akx / 1.d4
    Aky = Aky / 1.d4

    zk = zk / 1.d4

    phikx = 0.d0
    phiky = 0.d0
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setSurfaceNormals(grid,child)
                exit
             end if
          end do
       else
          cen = subcellCentre(thisOctal, subcell)
          x = cen%x
          y = cen%y
          z = cen%z
          d = thisOctal%subcellSize/2.d0 + 0.01d0*smallestCellSize

          do iLayer = 1, 7
             zLayer = skinInterface(x, y, Zk(iLayer), Akx(iLayer), Aky(iLayer), omegakx(iLayer), &
                  omegaky(iLayer), phikx(iLayer), phiky(iLayer))

             if (inSubcell(thisOctal,subcell,VECTOR(x,y,zlayer))) then
                thisOctal%surfaceNormal(subcell) = VECTOR(0.d0, 0.d0, 1.d0)
                xp = x
                yp = y
                zp = z - d
                vec = VECTOR(xp, yp, zp)
                if (inOctal(grid%octreeRoot, vec)) then
                   neighbourOctal => thisOctal
                   call findSubcellLocal(vec, neighbourOctal, neighbourSubcell)
                   neighbourOctal%surfaceNormal(neighbourSubcell) = VECTOR(0.d0, 0.d0, -1.d0)
                endif
             endif
          enddo

      endif
    enddo
  end subroutine setSurfaceNormals

  subroutine scatter(inComingDirection, normal, gFac, outGoingDirection)
    type(VECTOR) :: inComingDirection, normal, outgoingDirection
    integer, parameter :: nTheta = 1000, ng = 200
    real(double), save :: costhetaArray(nTheta), gArray(ng)
    real(double), save :: prob(ng, ntheta)
    logical, save :: firstTime = .true.
    real(double) :: gFac, theta, r
    integer :: i, j

    if (firstTime) then
       firstTime = .false.
       prob = 0.d0
       do i = 1, nTheta
          costhetaArray(i) = 2.d0 * dble(i-1)/dble(nTheta-1) - 1.d0
       enddo
       do i = 1, ng
          gArray(i) = dble(i-1)/dble(ng)
       enddo
       do i = 1, ng
          do j = 2, ntheta
             prob(i,j) = prob(i,j-1) + (1.d0 - gArray(i)**2)/(1.d0 + gArray(i)**2 - 2.d0 * gArray(i) * costhetaArray(j))**1.5d0
          enddo
          prob(i,1:nTheta) = prob(i,1:nTheta)/prob(i,nTheta)
       enddo
    endif
    call locate(gArray, ng, gfac, i)
    call randomNumberGenerator(getdouble=r)
    call locate(prob(i,1:nTheta), nTheta, r, j)
    theta = acos(min(1.d0,costhetaArray(j) + (costhetaArray(j+1)-costhetaArray(j))*(r - prob(i,j))/(prob(i,j+1)-prob(i,j))))
    outgoingDirection = arbitraryRotate(incomingDirection, theta,normal)
    call normalize(outgoingDirection)
    call randomNumberGenerator(getdouble=theta)
    theta = theta * twoPi
    outgoingDirection = arbitraryRotate(outgoingDirection, theta, incomingDirection)
    call normalize(outgoingDirection)
    normal = outgoingDirection .cross. incomingDirection
    call normalize(normal)
    if (modulus(outgoingDirection) < 0.9d0) write(*,*) "bug ",outGoingDirection,incomingDirection
  end subroutine scatter


  
  subroutine photonLoop(grid)
    use random_mod
    use inputs_mod, only : smallestCellSize
    use amr_utils_mod
    use messages_mod

    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, nextOctal
    type(VECTOR) :: photonPosition
    type(VECTOR) :: photonDirection
    type(VECTOR) :: photonNormal
    real(double) :: photonWeight, photonWavelength
    real(double) :: r
    real(double) :: albedo
    integer :: i, iTissue, subcell, nextSubcell, nPhoton
    real(double) :: muAbs, muSca, nRef, tau, tauToBoundary, tVal, gfac
    real(double) :: distToGrid, alpha, alphat
    logical :: hitGrid, absorbed
    real(double) :: detectorRadius
    real(double) :: nk, nkp1, criticalAngle, Rs, Rp
    type(VECTOR) :: detectorPosition, vec, normvec, probe
    type(PHOTONSOURCE) :: source


    call writeInfo("Creating photon source...",TRIVIAL)
    call createSource(source)
    call writeInfo("Done.",TRIVIAL)
    
    r = 0.d0
    do i = 1, 100000
       photonDirection = randomUnitVector()
       photonNormal = photonDirection.cross.VECTOR(1.d0,0.d0,0.d0)
       call normalize(photonNormal)
       call scatter(photonDirection, photonNormal, 0.5d0, vec)
       r = r + (photonDirection.dot.vec)
    enddo
    r = r/ 100000
    write(*,*) "average cos theta ", r
    detectorRadius = 500.d0*microntocm
    detectorPosition = VECTOR(0.02d0, 0.d0, 0.d0)
    call zeroDistanceGrid(grid%octreeRoot)
    nphoton = 0
    call writeInfo("Running photon loop...",TRIVIAL)
    do i = 1, nPhoton


       call  emitPhoton(source, photonPosition, photonDirection, photonWavelength, photonWeight, photonNormal)
       if (.not.inOctal(grid%octreeRoot, photonPosition)) then
          distToGrid = distanceToGridFromOutside(grid, photonPosition, photonDirection, hitGrid)
          if (.not.hitgrid) then
             write(*,*) "photon missed grid ",photonPosition,photonDirection
             cycle
          endif
          photonPosition = photonPosition + (distToGrid+1.d-6*smallestCellSize) * photonDirection
       endif

       thisOctal => grid%octreeRoot
       call findSubcellLocal(photonPosition, thisOctal, subcell)

       do while(inOctal(grid%octreeRoot, photonPosition))
          call findSubcellLocal(photonPosition, thisOctal, subcell)
          call randomNumberGenerator(getDouble=r)
          tau = -log(r)
          call distanceToCellBoundary(grid, photonPosition, photonDirection, tVal, thisOctal)
          iTissue = thisOctal%iTissue(subcell)
          muAbs = getMuAbs(tissueArray(iTissue), photonWavelength)
          muSca = getMuSca(tissueArray(iTissue), photonWavelength)
          nRef = getN(tissueArray(iTissue), photonWavelength)
          gfac = getgfac(tissueArray(iTissue), photonWavelength)
          tauToBoundary = tVal * (muAbs + muSca)
          if (tau > tauToBoundary) then
             thisOctal%distanceGrid(subcell) = thisOctal%distanceGrid(subcell) + tVal/(cSpeed/nRef) * photonWeight
             thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1

             probe = photonPosition + photonDirection*(tVal + 1.d-6*smallestCellSize)
             if (inOctal(grid%octreeRoot, probe)) then
                nextOctal => thisOctal
                call findSubcellLocal(probe, nextOctal, nextsubcell)
                if (thisOctal%iTissue(subcell) /= nextOctal%iTissue(nextSubcell)) then
!                   write(*,*) "crossed boundary ",photonDirection
                   alpha = acos(((-1.d0)*photonDirection).dot.thisOctal%surfaceNormal(Subcell))
!                   write(*,*) "alpha ",alpha*radtodeg
                   nk = getN(tissueArray(thisOctal%iTissue(subcell)), photonWavelength)
                   nkp1 = getN(tissueArray(nextOctal%iTissue(nextsubcell)), photonWavelength)
!                   write(*,*) "nk nkp1 ",nk,nkp1

                   criticalAngle = pi/2.d0
                   if (nkp1 < nk) then
                      criticalAngle = asin(nkp1/nk)
                   endif
                   if (alpha > criticalAngle) then
                      albedo = 1.d0
                   else
                      alphat  = asin((nk/nkp1)*sin(alpha))
!                      write(*,*) "alphat ",alphat*radtodeg
                      ! fresnel coefficent part
                      Rs = ((nk*cos(alpha) - nkp1*cos(alphat))/(nk*cos(alpha)+nkp1*cos(alphat)))**2
                      Rp = ((nk*cos(alphat) - nkp1*cos(alpha))/(nk*cos(alphat)+nkp1*cos(alpha)))**2
                      albedo = 0.5d0*(Rs + Rp)
                   endif
!                   write(*,*) "albedo ",albedo, thisOctal%surfaceNormal(subcell)
                   normVec = photonDirection .cross. thisOctal%surfaceNormal(subcell)
                   call normalize(normVec)
                   
                   call randomNumberGenerator(getDouble=r)
                   if (r < albedo) then ! reflection
                      photonPosition = photonPosition + photonDirection*(tVal - 1.d-6*smallestCellSize)
                      vec = photonDirection
                      photonDirection = (-1.d0)*photonDirection
                      photonDirection = arbitraryRotate(photonDirection, -2.d0*alpha, normVec)
                      call normalize(photonDirection)
!                      write(*,*) "reflection test ",radtodeg*alpha,&
!                           radtodeg*acos(photonDirection.dot.thisOctal%surfaceNormal(subcell))
                   else ! refraction
                      photonPosition = photonPosition + photonDirection*(tVal + 1.d-6*smallestCellSize)
                      vec = photonDirection
                      photonDirection = arbitraryRotate(photonDirection, alphat-alpha, normVec)
                      call normalize(photonDirection)
!                      write(*,*) "refraction test ", alphat, acos(photonDirection.dot.((-1.d0)*thisOctal%surfaceNormal(subcell)))
                   endif
                else
                   photonPosition = photonPosition + photonDirection*(tVal + 1.d-6*smallestCellSize)
                endif
             else
                photonPosition = photonPosition + photonDirection*(tVal + 1.d-6*smallestCellSize)
             endif
          else
             tVal = tval * tau / tauToBoundary
             photonPosition = photonPosition + photonDirection*tval
             thisOctal%distanceGrid(subcell) = thisOctal%distanceGrid(subcell) + tVal/(cSpeed/nRef) * photonWeight
             thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1
             call randomNumberGenerator(getDouble=r)
             absorbed = .false.
             albedo = muSca / (muAbs + muSca)
             if (r < albedo) then
                call scatter(photonDirection, photonNormal, gfac, vec)
                photonNormal = photonDirection.cross.vec
                call normalize(photonNormal)
                photonDirection = Vec

             else
                absorbed = .true.
                exit
             endif
          endif
       enddo
       if (.not.absorbed.and.(photonDirection%z < 0.d0)) then
          if (modulus(photonPosition - detectorPosition) < detectorRadius) then
!             write(*,*) "photon captured ", photonPosition
          endif
       endif

    end do
    call writeInfo("Done.",TRIVIAL)
    call calculateEnergyDensity(grid%octreeRoot)


  end subroutine photonLoop
  
  real(double) function getMuAbs(tissue, wavelength)
    type(TISSUETYPE) :: tissue
    real(double) :: wavelength
    getMuAbs = interpArray(wavelength, tissue%nLambda, tissue%lamArray, tissue%muAbs)
  end function getMuAbs

  real(double) function getMuSca(tissue, wavelength)
    type(TISSUETYPE) :: tissue
    real(double) :: wavelength
    getMuSca = interpArray(wavelength, tissue%nLambda, tissue%lamArray, tissue%muSca)
  end function getMuSca

  real(double) function getN(tissue, wavelength)
    type(TISSUETYPE) :: tissue
    real(double) :: wavelength
    getN = interpArray(wavelength, tissue%nLambda, tissue%lamArray, tissue%n)
  end function getN

  real(double) function getGfac(tissue, wavelength)
    type(TISSUETYPE) :: tissue
    real(double) :: wavelength
    getGfac = interpArray(wavelength, tissue%nLambda, tissue%lamArray, tissue%gfac)
  end function getGfac

  recursive subroutine zeroDistanceGrid(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call zeroDistanceGrid(child)
                exit
             end if
          end do
       else
          thisOctal%distanceGrid(subcell) = 0.d0
          thisOctal%nCrossings(subcell) = 0
          thisOctal%undersampled(subcell) = .false.
       endif
    enddo
  end subroutine zeroDistanceGrid

  recursive subroutine calculateEnergyDensity(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateEnergyDensity(child)
                exit
             end if
          end do
       else
          thisOctal%uDens(subcell) = thisOctal%distanceGrid(subcell)/cellVolume(thisOctal,Subcell)
       endif
    enddo
  end subroutine calculateEnergyDensity

end module biophysics_mod
