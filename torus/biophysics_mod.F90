module biophysics_mod

    use vector_mod
    use kind_mod
    use utils_mod
    use gridtype_mod
    use spectrum_mod
    use detector_mod

    use octal_mod

    use triangle_mesh_mod, only: triangleMesh
    use object_media_mod

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


    type pMedium
        type(medium), pointer :: p =>  null()
    end type pMedium

    type(triangleMesh), dimension(:), allocatable :: objectMeshList
    type(pMedium), dimension(:), allocatable :: objectMediaList
    type(medium),  pointer, dimension(:) :: tissueList


contains


    subroutine setupMatcherSkinTissue()

        integer, parameter :: nTissue = 7

        real(double) ::  lambda, cH2O, cMelanin
        real(double), allocatable :: lambdaList(:)
        integer :: nLambda, i

        type(medium) :: water, Hb, HbO2
        type(mediumProps) :: props

        water = defineWater(35.5d0)
        Hb = defineHemoglobin(.false.)
        HbO2 = defineHemoglobin(.true.)

        if(associated(tissueList)) return ! if list already exists
        allocate(tissueList(0:nTissue))

        nLambda = 100
        allocate(LambdaList(nLambda))
        do i = 1, nLambda
            LambdaList(i) = 400.d0 + 400.d0 * dble(i-1)/dble(nLambda-1)
        enddo

        ! air

        tissueList(0) =  createNewMedium("Air", lambdaList)

        ! stratum corneum

        tissueList(1) = createNewMedium("Stratum corneum", lambdaList)
        cH2O = 0.05d0
        do i = 1, nLambda
            lambda = lambdaList(i)
            props = getMediumProps(water, lambda)
            tissueList(1)%propList(i)%muAbsorb = ( (1.d0-0.3d-3*lambda)+0.125d0*mu_a_0(lambda) ) &
                * (1.d0 - cH2O) + cH2O * props%muAbsorb
        enddo
        tissueList(1)%propList%muScatter = 1000.d0
        tissueList(1)%propList%gFactor = 0.86
        tissueList(1)%propList%nRefract = 1.5d0

        ! living epidermis

        tissueList(2) = createNewMedium("Living epidermis", lambdaList)
        do i = 1, nLambda

            lambda = lambdaList(i)
            cMelanin = 0.5!!!!!!!!!!
            props = getMediumProps(water, lambda)

            tissueList(2)%propList(i)%muAbsorb = (1.d0-cH2O)*(cMelanin * muMelanin(lambda) + &
                (1.d0 - cMelanin)*mu_a_0(lambda) ) + cH2O * props%muAbsorb
        enddo
        tissueList(2)%propList%muScatter = 450.d0
        tissueList(2)%propList%gFactor = 0.8
        tissueList(2)%propList%nRefract = 1.34d0

        ! Papillary dermis

        tissueList(3) = createNewMedium("Papillary dermis", lambdaList)
        call fillTissueMatcher(tissueList(3), 0.04d0, 0.5d0, 300.d0, 0.9d0, 1.4d0)

        ! Upper blood net dermis

        tissueList(4) = createNewMedium("Upper blood net dermis", lambdaList)
        call fillTissueMatcher(tissueList(4), 0.9d0, 0.1d0, 350.d0, 0.95d0, 1.39d0)


        ! Reticular dermis

        tissueList(5) =  createNewMedium("Reticular dermis", lambdaList)
        call fillTissueMatcher(tissueList(5), 0.04d0, 0.7d0, 250.d0, 0.8d0, 1.4d0)

        ! Deep blood net dermis

        tissueList(6) = createNewMedium("Deep blood net dermis", lambdaList)
        call fillTissueMatcher(tissueList(6), 0.1d0, 0.7d0, 300.d0, 0.95d0, 1.38d0)

        ! Subcutaneous fat

        tissueList(7) = createNewMedium("Subcutaneous fat", lambdaList)
        call fillTissueMatcher(tissueList(7), 0.05d0, 0.7d0, 50.d0, 0.75d0, 1.44d0)

    contains

        subroutine fillTissueMatcher(tissue, cBlood, cH2O, muSca, g, n)
            ! Mixes material properties based on given
            !   blood concentration cBlood
            !   water concentration cH20
            !   blood oxygen saturation S
            !use
            type(medium), intent(inout) :: tissue
            real(double), intent(in) :: cBlood, cH2O, muSca, g, n
            real(double), parameter :: FRBC = 0.99d0
            real(double), parameter :: FHb = 0.25d0
            real(double), parameter :: S = 0.6d0
            real(double), parameter :: Ht = 0.45d0
            real(double) :: gamma, lambda
            integer :: i

            type(mediumprops) :: h2o_props, hb_props, hbo2_props

            gamma = FHb * FRBC * Ht

            do i = 1, size(tissue%lambdaList)
                lambda = tissue%lambdaList(i)
                h2o_props = getMediumProps(water, lambda)
                hb_props = getMediumProps(Hb, lambda)
                hbo2_props = getMediumProps(HbO2, lambda)
                !ho_props = getMediumProps(water, lambda)
                !tissue%propList(i)%muAbsorb = gamma*cBlood*((1.d0-S)*getMediumProps(Hb,lambda)%muAbsorb &
                tissue%propList(i)%muAbsorb = gamma*cBlood*((1.d0-S)*hb_props%muAbsorb &
                                                               + S * hbo2_props%muAbsorb) &
                                      + (1.d0-gamma*cBlood)*( (1.d0- cH2O) * mu_a_0(lambda) &
                                                            + cH2O * h2o_props%muAbsorb)
            enddo
            tissue%propList%muScatter = muSca
            tissue%propList%gFactor   = g
            tissue%propList%nRefract  = n
        end subroutine fillTissueMatcher

        real(double) function mu_a_0 (lambda)
            real(double) :: lambda
            mu_a_0 = 7.84d8 * lambda**(-3.255d0) ! equation 10
            end function mu_a_0

        real(double) function muMelanin (lambda)
            real(double) :: lambda
            muMelanin = 5.d10 * lambda**(-3.33d0)
        end function muMelanin

    end subroutine setupMatcherSkinTissue


    subroutine sortTriangleMeshList(meshList)
      use triangle_mesh_mod
      use utils_mod, only : indexx
        type(triangleMesh), dimension(:) :: meshList
        type(triangleMesh), allocatable, dimension(:) :: tmpMeshList
        integer :: n, i
        real(double), allocatable :: zList(:)
        integer, allocatable :: iArray(:)
        n = size(meshList)
        allocate(zList(1:n))
        allocate(iArray(1:n))
        allocate(tmpMeshList(1:n))

        do i = 1, n
           zList(i) = averageZ(meshList(i))
           if (meshList(i)%isClosed)  zList(i) = -1.d10
        enddo

        call indexx(n, zList, iarray)

        do i = 1, n
           tmpMeshList(i) = meshList(iarray(i))
        enddo
        
        do i = 1, n
           meshList(i) = tmpmeshList(i)
        enddo
        
    end subroutine sortTriangleMeshList



    subroutine initialize3dObjects()
        use inputs_mod, only : wavefrontFile
        use triangle_mesh_mod, only: readTriangleMesh, printTriangleMeshSummary
        integer :: i, n
        type(medium), pointer :: vacuum, water, crownglass, flintglass, air, blood, rubber, &
             intralipid10

        allocate(vacuum, water, crownglass, flintglass, air, blood, rubber, intralipid10)
        !vacuum = createNewMedium('vacuum', (/0.d0, 1.d4/) )
        water = defineWater(25.d0)
        flintglass = defineGlass("denseflint.dat")
        crownglass = defineGlass("crownglass.dat")
        air = defineAir()
        blood = defineHemoglobin(.true.)
        rubber = defineRubber()
        intralipid10 = defineFatSolution("intralipid10")
        call setupMatcherSkinTissue()

        objectmeshlist = readTriangleMesh(wavefrontFile)
        n = size(objectMeshList)
        if (n > 1) call sortTriangleMeshList(objectMeshList)

        allocate(objectMediaList(0:n))

        objectMediaList(0)%p => air
        do i = 1, n
           select case(objectMeshList(i)%material)
              case("air")
                 objectMediaList(i)%p => air
              case("crownglass")
                 objectMediaList(i)%p => crownglass
              case("flintglass")
                 objectMediaList(i)%p => flintglass
              case("water")
                 objectMediaList(i)%p => water
              case("blood")
                 objectMediaList(i)%p => blood
              case("rubber")
                 objectMediaList(i)%p => rubber
              case("keratin")
                 objectMediaList(i)%p => rubber
              case("intralipid10")
                 objectMediaList(i)%p => intralipid10
              case("stratum_corneum")
                 objectMediaList(i)%p => tissueList(1)
              case("living_epidermis")
                 objectMediaList(i)%p => tissueList(2)
              case("papillary_dermis")
                 objectMediaList(i)%p => tissueList(3)
              case("upper_blood_net_dermis")
                 objectMediaList(i)%p => tissueList(4)
              case("dermis")
                 objectMediaList(i)%p => tissueList(5)
              case("deep_blood_net_dermis")
                 objectMediaList(i)%p => tissueList(6)
              case DEFAULT
                 call writeFatal("Do not recognise material: "//trim(objectMeshList(i)%material))
                 stop
            end select
         enddo

    end subroutine initialize3dObjects

    subroutine summarize3dObjects()
        use triangle_mesh_mod, only: printTriangleMeshSummary
        integer :: i, n

        if(.not. allocated(objectMeshList) .or. .not.allocated(objectMediaList)) then
            print *, 'objectMeshList or objectMedium list not yet defined. Initialize first.'
            return
        end if

        if (writeoutput) print *, 'summary of objects and their materials:'

        n = size(objectMeshList)
        do i = 1,n
            call printTriangleMeshSummary(objectMeshList(i))
            if (writeoutput) then
               if(associated(objectMediaList(i)%p)) then
                  print *,'medium: ', objectMediaList(i)%p%label
               else
                  print *, 'medium undefined'
               end if
            endif
        end do

    end subroutine summarize3dObjects

    type(medium) function findMedium(j)
        integer :: j

        findMedium = objectMediaList(j)%p

    end function findMedium

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
                    call setupSpectrum(source%component(i)%spectrum, 5000.d0)
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
        real(double) :: phi, r, thisR, weight
!$OMP THREADPRIVATE (prob)
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

        !normVec = source%component(iComponent)%direction .cross. xAxis
        !photonNormal = normVec
        !normVec = normVec * thisR
        !vec = source%component(iComponent)%direction
        !normVec = arbitraryRotate(normVec, phi, vec)
        photonNormal = source%component(iComponent)%direction .cross. randomUnitVector()
        call normalize(photonNormal)
        photonPosition = photonPosition + thisR * photonNormal
        photonDirection = source%component(iComponent)%direction
        photonWavelength = source%component(iComponent)%wavelength

        call getWavelength(source%component(iComponent)%spectrum, photonWavelength, weight)
        photonWavelength = photonWavelength / 10.d0 ! angstrom to nm

        photonWeight  = source%component(iComponent)%power/source%totalPower

        if (source%component(iComponent)%iType == 1) then
            call randomNumberGenerator(getDouble=r)
            phi = asin(source%component(iComponent)%numericalAperture/1.5d0)
            phi = acos((1.d0-cos(phi))*r+cos(phi))
            photonDirection = arbitraryRotate(photonDirection,phi,photonNormal)
            call randomNumberGenerator(getDouble=r)
            phi = twoPi * r
            photonDirection = arbitraryRotate(photonDirection,phi,source%component(iComponent)%direction)
        endif
    end subroutine emitPhoton

    subroutine scatter(inComingDirection, normal, gFac, outGoingDirection)
        type(VECTOR) :: inComingDirection, normal, outgoingDirection
        integer, parameter :: nTheta = 1000, ng = 200
        real(double), save :: costhetaArray(nTheta), gArray(ng)
        real(double), save :: prob(ng, ntheta)
        logical, save :: firstTime = .true.
        real(double) :: gFac, theta, r
        integer :: i, j

!$OMP THREADPRIVATE(firstTime, costhetaarray, garray,prob)

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
                    prob(i,j) = prob(i,j-1) + (1.d0 - gArray(i)**2)/(1.d0 + gArray(i)**2 &
                                - 2.d0 * gArray(i) * costhetaArray(j))**1.5d0
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



    function photonInterfaceInteraction(directionIn, surfNormal, n1, n2, enter) result(directionOut)
        use random_mod

        type(vector), intent(in) :: directionIn
        type(vector), intent(in) :: surfNormal
        real(double), intent(in) :: n1, n2
        type(vector) :: directionOut
        logical, intent(out), optional :: enter ! true if photon enters other material (refract), false for reflect
        type(vector) :: normal, reflectDirection, refractDirection
        real(double) :: Rs, Rp, albedo, nRatio, criticalVal
        real(double) :: cosIn, cosOut
        real(double) :: r

        ! use Snell's law in vector format (equations from Wikipedia)
        normal = surfNormal
        cosIn = -(directionIn .dot. normal)
        if (cosIn < 0.d0) then ! orient surfNormal if necessary
            normal = (-1.d0)*normal
            cosIn = -cosIn
        endif
        reflectDirection = directionIn + 2.d0*cosIn*normal

        nRatio = n1/n2

        criticalVal = 1.d0 - nRatio*nRatio*(1.d0-cosIn*cosIn)
        if (criticalVal <= 0.d0) then ! full internal reflection
            albedo = 1.d0
        else ! determine refraction
            cosOut = sqrt(criticalVal)
            refractDirection = nRatio*directionIn + (nRatio*cosIn - cosOut)*normal

            ! Fresnel coefficients (equations from Wikipedia)
            Rs = ( (n1*cosIn - n2*cosOut)/(n1*cosIn + n2*cosOut) )**2
            Rp = ( (n1*cosOut - n2*cosIn)/(n1*cosOut + n2*cosIn) )**2
            albedo = 0.5d0*(Rs + Rp)
        endif

        call randomNumberGenerator(getDouble=r)
        if (r < albedo) then ! reflection
            directionOut = reflectDirection
            if(present(enter)) enter = .false.
        else ! refraction
            directionOut = refractDirection
            if(present(enter)) enter = .true.
        endif
        call normalize(directionOut)
        return
    end function photonInterfaceInteraction

    logical function splitSkin(thisOctal, subcell)
      use interface_hierarchy_mod, only: makeLocalInterfaceHierarchy, refineLocalInterfaceHierarchy, &
                                         countInterfaceTriangles, printLocalInterfaceHierarchy, findGlobalObjectID, &
                                         returnGlobalIdFromfullMesh, combineAllIds
      use triangle_mesh_mod, only : isPointInTriangleMesh
      use amr_utils_mod
      type(OCTAL) :: thisOctal
      integer :: subcell
      type(VECTOR) :: cellCentre, boxMin, boxMax
      real(double) :: d
      type(localInterfaceHierarchy) :: parentRoot
      !type(localInterfaceHierarchy) :: test
      type(localInterfaceHierarchy), pointer :: thisNode, thisControl

      splitSkin = .false.

      if(.not.allocated(objectMeshList)) then
         call initialize3dObjects()
         call summarize3dObjects()
      end if


      cellCentre = subcellCentre(thisOctal,subcell)
      d = thisOctal%subcellSize/2.d0
      boxMin = VECTOR(cellcentre%x - d, cellCentre%y - d, cellCentre%z - d)
      boxMax = VECTOR(cellcentre%x + d, cellCentre%y + d, cellCentre%z + d)


      !if (.not.associated(thisOctal%localInterfaces)) allocate(thisOctal%localInterfaces(1:thisOctal%maxChildren))
      !if (associated(thisOctal%parent) .and. associated(thisOctal%parent%localInterfaces) ) then
      !   parentRoot = thisOctal%parent%localInterfaces(thisOctal%parentSubcell)
      !   thisOctal%localInterfaces(subcell) = refineLocalInterfaceHierarchy(parentRoot, boxMin, boxMax)
      !else
      !    thisOctal%localInterfaces(subcell) = makeLocalInterfaceHierarchy(objectMeshList, boxMin, boxMax)
      !end if

      ! === TEST CODE START
      if (.not.associated(thisOctal%localInterfaces)) allocate(thisOctal%localInterfaces(1:thisOctal%maxChildren))

!      if (associated(thisOctal%parent) .and. associated(thisOctal%parent%localInterfaces) ) then
!         parentRoot = thisOctal%parent%localInterfaces(thisOctal%parentSubcell)
!         debug = .false.
!         if (inSubcell(thisOctal, subcell, VECTOR(-2d-3, 1.e-5,1.e-4))) debug = .true.
!         thisOctal%localInterfaces(subcell) = refineLocalInterfaceHierarchy(parentRoot, boxMin, boxMax, debug=debug)
!      else

      if (associated(thisOctal%parent) .and. associated(thisOctal%parent%localInterfaces) ) then
         parentRoot = thisOctal%parent%localInterfaces(thisOctal%parentSubcell)
         thisOctal%localInterfaces(subcell) = refineLocalInterfaceHierarchy(parentRoot, boxMin, boxMax)
      else
          thisOctal%localInterfaces(subcell) = makeLocalInterfaceHierarchy(objectMeshList, boxMin, boxMax)
      end if

      if (.not.associated(thisOctal%localInterfacesCONTROL)) allocate(thisOctal%localInterfacesCONTROL(1:thisOctal%maxChildren))
      thisOctal%localInterfacesCONTROL(subcell) = makeLocalInterfaceHierarchy(objectMeshList, boxMin, boxMax, &
                                                                    thisOctal%localInterfaces(subcell)%reference)
      thisNode => thisOctal%localInterfaces(subcell)
      thisControl => thisOctal%localInterfacesCONTROL(subcell)

!        print *,'  start test'
        do while(associated(thisNode) .or. associated(thisControl) )

!            print *, '    check node'
            if(.not. (associated(thisNode) .and. associated(thisControl) ) ) then
                print *, 'node mismatch'
                exit
            end if

            !print *, '    check id'
!            print *, '    id ref', thisNode%globalObjectId
!            print *, '    id con', thisControl%globalObjectId
            if(thisNode%globalObjectId /= thisControl%globalObjectId) then
                print *, 'ID mismatch'
                exit
            end if

!            print *, '    hasTri ref', associated(thisNode%triangleList)
!            print *, '    hasTri con', associated(thisControl%triangleList)
            if(associated(thisNode%triangleList) .NEQV. associated(thisControl%triangleList)) then
                print *, 'triangle mismatch'
                exit
            end if

            if (associated(thisNode%triangleList) .and. associated(thisControl%triangleList)) then
!                print *, '    nTri ref', size(thisNode%triangleList)
!                print *, '    nTri con', size(thisControl%triangleList)
                if(size(thisNode%triangleList) /= size(thisControl%triangleList)) then
                    print *, 'triangle count mismatch'
                    exit
                end if
            end if

!            print *, '  next'
            thisNode => thisNode%nextNode
            thisControl => thisControl%nextNode
        end do

      ! === TEST CODE END

      if ( countInterfaceTriangles(thisOctal%localInterfaces(subcell)) > 0) then
         splitSkin = .true.
      endif

      if (.not.associated(thisOctal%iTissue)) allocate(thisOctal%itissue(1:thisOctal%maxChildren))
      thisOctal%iTissue(Subcell) = combineAllIds(thisOctal%localInterfaces(subcell))

    end function splitSkin


    recursive subroutine compareMesh(cell)
        use interface_hierarchy_mod

        type(octal),  intent(in) :: cell

        type(localInterfaceHierarchy), pointer :: thisNode, thisControl
        integer :: i

        print *, 'depth:', cell%nDepth

        if( associated(cell%localInterfaces) .neqv. associated(cell%localInterfacesCONTROL) ) then
            print *, 'interface mismatch'
            return
        end if

        if( .not. associated(cell%localInterfaces) ) return

        do i = 1,cell%nChildren
            print *,'  i = ', i
            thisNode => cell%localInterfaces(i)
            thisControl => cell%localInterfacesCONTROL(i)

            print *,'  start loop'
            do while(associated(thisNode) .or. associated(thisControl) )

                print *, '    check node'
                if(.not. (associated(thisNode) .and. associated(thisControl) ) ) then
                    print *, 'node mismatch'
                    return
                end if

                !print *, '    check id'
                print *, '    id ref', thisNode%globalObjectId
                print *, '    id con', thisControl%globalObjectId
                if(thisNode%globalObjectId /= thisControl%globalObjectId) then
                    print *, 'ID mismatch'
                    return
                end if

                print *, '    hasTri ref', associated(thisNode%triangleList)
                print *, '    hasTri con', associated(thisControl%triangleList)
                if(associated(thisNode%triangleList) .NEQV. associated(thisControl%triangleList)) then
                    print *, 'triangle mismatch'
                    return
                end if

                if (associated(thisNode%triangleList) .and. associated(thisControl%triangleList)) then
                    print *, '    nTri ref', size(thisNode%triangleList)
                    print *, '    nTri con', size(thisControl%triangleList)
                    if(size(thisNode%triangleList) /= size(thisControl%triangleList)) then
                        print *, 'triangle count mismatch'
                        return
                    end if
                end if

                print *, '  next'
                thisNode => thisNode%nextNode
                thisControl => thisControl%nextNode
            end do
            print *, '  end'

            print *, '  children'
            !if(associated(cell%child)) then
!            if(cell%hasChild(i)) then
                call compareMesh(cell%child(i))
!                if(associated(child)) call compareMesh(child)
!            end if
            print *, '  end children'

        end do

    end subroutine compareMesh


    subroutine photonLoop(grid)
#ifdef MPI
      use mpi_global_mod
#endif
        use random_mod
        use inputs_mod, only : smallestCellSize, nPhotons, detFilename
        use amr_utils_mod
        use messages_mod

        use interface_hierarchy_mod, only: objectIntersection, findGlobalObjectId, findFirstIntersection

        type(GRIDTYPE) :: grid
        type(OCTAL), pointer :: thisOctal
        type(VECTOR) :: photonPosition, photonDirection, photonNormal
        real(double) :: photonWeight, photonWavelength
        real(double) :: r, albedo
        integer(bigInt) :: i
        integer :: thisObjectId, subcell, countEvent
        real(double) :: tau, tauToBoundary, dist
        logical :: hitGrid, absorbed
        type(VECTOR) :: vec
        type(PHOTONSOURCE) :: source

        logical :: subcellHasInterface, hitsInterface, doesEnter, firstwrite, stored
        type(mediumProps) :: thisMedium, nextMedium
        type(objectIntersection) :: intersect
        type(DETECTORTYPE) :: thisDetector
        integer, parameter :: maxEvents = 100000
        type(VECTOR) :: cellEvent(maxEvents)
        real(double) :: lengthEvent(maxEvents)
        integer :: nEvent
        integer(bigint) :: iMonte_beg, iMonte_end
#ifdef MPI
        integer(bigint) :: m, n_rmdr, np
#endif
        call createNewDetector(thisDetector)

        call writeInfo("Creating photon source...",TRIVIAL)
        call createSource(source)
        call writeInfo("Done.",TRIVIAL)

        call zeroDistanceGrid(grid%octreeRoot)

        call writeInfo("Running photon loop...",TRIVIAL)



          imonte_beg=1; imonte_end=nPhotons  ! default value

#ifdef MPI
                 ! Set the range of index for a photon loop used later.     
                 np = nThreadsGlobal
                 n_rmdr = MOD(nPhotons,np)
                 m = nPhotons/np
          
                 if (myRankGlobal .lt. n_rmdr ) then
                    imonte_beg = (m+1)*myRankGlobal + 1
                    imonte_end = imonte_beg + m
                 else
                    imonte_beg = m*myRankGlobal + 1 + n_rmdr
                    imonte_end = imonte_beg + m -1
                 end if
             !    print *, ' '
             !    print *, 'imonte_beg = ', imonte_beg
             !    print *, 'imonte_end = ', imonte_end
            
                
                 !  Just for safety.
                 if (imonte_end .gt. nPhotons .or. imonte_beg < 1) then
                    print *, 'Index out of range: i_beg and i_end must be ' 
                    print *, ' 0< index < ', nPhotons , '    ... [lucy_mod::lucyRadiativeEquilibriumAMR]'
                    print *, 'imonte_beg = ', imonte_beg
                    print *, 'imonte_end = ', imonte_end
                    stop
                 end if
#endif    


!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(i, firstWrite, photonNormal, photonPosition, photonWeight, photonDirection, photonWavelength, nEvent) &
!$OMP PRIVATE(dist, hitgrid, absorbed, countEvent, thisOctal, subcell, thisObjectID, thisMedium, subcellHasInterface) &
!$OMP PRIVATE(intersect, tautoboundary, r, tau, cellEvent, lengthEvent, albedo, hitsInterface, vec, nextMedium, doesenter, stored) &
!$OMP SHARED(nPhotons, thisDetector, source, grid, smallestCellSize, iMonte_beg, iMonte_end, writeoutput) 

!$OMP DO SCHEDULE(DYNAMIC,2)
        do i = iMonte_beg, iMonte_end

           if (mod(i-iMonte_beg+1,(iMonte_end-iMonte_beg)/10)==0) then
              if (writeoutput) write(*,*) nint(100.*real(i-iMonte_beg+1)/real(iMonte_end-iMonte_beg)), " % complete"
           endif


           firstwrite=.true.

            call  emitPhoton(source, photonPosition, photonDirection, photonWavelength, photonWeight, photonNormal)
            nEvent = 0
            if (.not.inOctal(grid%octreeRoot, photonPosition)) then
                dist = distanceToGridFromOutside(grid, photonPosition, photonDirection, hitGrid)
                if (.not.hitgrid) then
                    write(*,*) "photon missed grid ",photonPosition,photonDirection
                    cycle
                endif
                photonPosition = photonPosition + (dist+1.d-6*smallestCellSize) * photonDirection
            endif
            absorbed = .false.
            countEvent = 0

            thisOctal => grid%octreeRoot
            call findSubcellLocal(photonPosition, thisOctal, subcell)
            if (.not.associated(thisOctal%localInterfaces)) then
               write(*,*) "local interfaces not assocaited for octal of depth ",thisoctal%nDepth
               write(*,*) "check parent ",associated(thisOctal%parent%localinterfaces)
               thisObjectID = 0 ! background
            else
               thisObjectId = findGlobalObjectId( thisOctal%localInterfaces(subcell), photonPosition)
            endif

            thisMedium = getMediumProps(findMedium(thisObjectId)  , photonWaveLength)
            subcellHasInterface = .false.
            if (associated(thisOctal%localInterfaces)) then
               subcellHasInterface = associated(thisOctal%localInterfaces(subcell)%triangleList)
            endif


            do ! infinite loop
!               write(*,*) "pos ",photonPosition
                !call findSubcellLocal(photonPosition, thisOctal, subcell)
                call distanceToCellBoundary(grid, photonPosition, photonDirection, dist, thisOctal)

                if(subcellHasInterface) then
                    intersect = findFirstIntersection(thisOctal%localInterfaces(subcell), &
                                                      photonPosition, photonDirection, thisObjectId)


                    if(intersect%globalObjectId > -1 .and. intersect%distance < dist) then ! if interface is closer than boudary
                        hitsInterface = .true.
                        dist = intersect%distance



                    if (intersect%globalObjectID==1.and.firstwrite) then
                       firstwrite=.false.
                    endif


                    else
                        hitsInterface = .false.
                    end if
                else ! no interface to hit
                    hitsInterface = .false.
                end if
                tauToBoundary = dist * (thisMedium%muAbsorb + thisMedium%muScatter)
                call randomNumberGenerator(getDouble=r)
                absorbed = .false.
                tau = -log(r)
!                call getCompositeBiasedTau(tau, weight)
!                write(*,*) "tau weight " ,tau, weight
!                photonWeight = photonWeight * weight
                if (tau < tauToBoundary ) then ! photon has event before hitting any boundary/interface
                    dist = dist * tau / tauToBoundary

                    ! some admin
!$OMP ATOMIC
                    thisOctal%distanceGrid(subcell) = thisOctal%distanceGrid(subcell) + &
                                                      dist/(cSpeed/thisMedium%nRefract) * photonWeight

!$OMP ATOMIC
                    thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1

                    nEvent = nEvent + 1
                    cellEvent(nEvent) = subcellCentre(thisOctal, subcell)
                    lengthEvent(nEvent) = dist/(cSpeed/thisMedium%nRefract) * photonWeight

                    photonPosition = photonPosition + photonDirection*dist
                    ! determine event: scatter or absorption
                    call randomNumberGenerator(getDouble=r)
                    absorbed = .false.
                    albedo = thisMedium%muScatter / (thisMedium%muAbsorb + thisMedium%muScatter)
                    if (r < albedo) then ! scatter
                        call scatter(photonDirection, photonNormal, thisMedium%gFactor, vec)
                        photonNormal = photonDirection.cross.vec
                        call normalize(photonNormal)
                        photonDirection = vec
                    else ! absorbed
                        absorbed = .true.
                        exit
                    endif

                else if (hitsInterface) then ! photon hits an interface inside current subcell
                    ! some admin

!$OMP ATOMIC
                    thisOctal%distanceGrid(subcell) = thisOctal%distanceGrid(subcell) + &
                                                      dist/(cSpeed/thisMedium%nRefract) * photonWeight
!$OMP ATOMIC
                    thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1

                    nEvent = nEvent + 1
                    cellEvent(nEvent) = subcellCentre(thisOctal, subcell)
                    lengthEvent(nEvent) = dist/(cSpeed/thisMedium%nRefract) * photonWeight


                    photonPosition = photonPosition + photonDirection*dist
                    nextMedium = getMediumProps(findMedium(intersect%globalObjectId), photonWavelength)

                    photonDirection = photonInterfaceInteraction(photonDirection, intersect%normal, &
                                                                 thisMedium%nRefract, nextMedium%nRefract, &
                                                                 doesEnter) ! <= optional output!
                    photonPosition = photonPosition + 1d-6*smallestCellSize*photonDirection ! tweak away from intersection
!                    if (doesEnter) then
!                       write(*,*) " photon enters next medium ",thisObjectId,intersect%globalObjectId
!                    endif
                    if(doesEnter) then ! assign next objectid and medium
                        thisObjectId = intersect%globalObjectId
                        thisMedium = nextMedium
                    end if

                else ! photon crosses an (sub)cell boundary
                    ! some admin
!$OMP ATOMIC
                    thisOctal%distanceGrid(subcell) = thisOctal%distanceGrid(subcell) + &
                                                      dist/(cSpeed/thisMedium%nRefract) * photonWeight
!$OMP ATOMIC
                    thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1


                    nEvent = nEvent + 1
                    cellEvent(nEvent) = subcellCentre(thisOctal, subcell)
                    lengthEvent(nEvent) = dist/(cSpeed/thisMedium%nRefract) * photonWeight

                    photonPosition = photonPosition + photonDirection*(dist+1d-6*smallestCellSize)

                    ! find next subcell
                    if(inOctal(grid%octreeRoot, photonPosition)) then
                        call findSubcellLocal(photonPosition, thisOctal, subcell)
                        subcellHasInterface = .false.
                        if (associated(thisOctal%localInterfaces)) then
                           subcellHasInterface = associated(thisOctal%localInterfaces(subcell)%triangleList)
                        endif
                    else
                        exit
                    end if
                end if

                countEvent = countEvent + 1

            end do ! end while: photon left the grid

            if (.not.absorbed) then
               call storePhotonOnDetector(thisDetector, photonPosition, photonDirection, &
                    photonWeight, photonWavelength, stored)
               if (stored) then
                  call addEventsToGrid(nEvent, cellEvent, lengthEvent, grid)
               endif
            else
!               write(*,*) "photon absorbed after ", countEvent, "events"
            endif
        end do ! all photons
!$OMP END DO

!$OMP END PARALLEL

#ifdef MPI
        call updateGridMPI(grid)
        call gatherDetector(thisDetector)
#endif

        if (writeoutput) call writeDetectorfile(thisDetector, detFilename)
        call writeInfo("Done.",TRIVIAL)
        call calculateEnergyDensity(grid%octreeRoot)


    end subroutine photonLoop

    subroutine getCompositeBiasedTau(tau, weight)
      real(double) :: tau, weight, r, alpha, fac
      real(double), parameter :: tauPath = 20.d0
      integer, parameter :: ntau = 1000
      real(double), save :: prob(ntau),tauArray(ntau)
      real(double), parameter :: epsilon = 0.1d0
      integer :: i
      logical, save :: firstTime = .true.

 !$OMP THREADPRIVATE(prob,tauArray,firstTime)

      alpha = 1.d0 / (1.d0 + tauPath)
      if (firstTime) then
         firstTime = .false.
         prob(1) = 0.d0
         do i = 2, nTau
            tauArray(i) = tauPath * dble(i-1)/dble(ntau-1)
            prob(i) = prob(i-1) + (1.d0-epsilon) * exp(-tauArray(i)) + epsilon*alpha*exp(-alpha*tauArray(i))
         enddo
         prob(1:nTau) = prob(1:nTau)/prob(nTau)
      endif

      call randomNumberGenerator(getDouble=r)
      call locate(prob, ntau, r, i)
      fac = (r-prob(i))/(prob(i+1)-prob(i))
      tau = tauArray(i) + (tauArray(i+1)-tauArray(i))*fac
      weight = 1.d0/( (1.d0-epsilon) + epsilon*alpha*exp( (1.d0-alpha)*tau))
    end subroutine getCompositeBiasedTau


    subroutine addEventsToGrid(nEvent, cellEvent, lengthEvent, grid)
      use amr_utils_mod
      integer :: nEvent, i
      type(VECTOR) :: cellEvent(:)
      real(double) :: lengthEvent(:)
      type(GRIDTYPE) :: grid
      type(OCTAL), pointer :: thisOctal
      integer :: subcell
      thisOctal => grid%octreeRoot
      do i = 1, nEvent
         call findsubcellLocal(cellEvent(i), thisOctal, subcell)
         thisOctal%etaLine(subcell) = thisOctal%etaLine(subcell) + lengthEvent(subcell)
      enddo
    end subroutine addEventsToGrid
         
      


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
                thisOctal%etaline(subcell) = 0.d0
                thisOctal%nCrossings(subcell) = 0
                thisOctal%undersampled(subcell) = .false.
            endif
        enddo
    end subroutine zeroDistanceGrid

    recursive subroutine setTissues(thisOctal)
        use interface_hierarchy_mod, only: findGlobalObjectId
        type(octal), pointer   :: thisOctal
        type(octal), pointer  :: child
        integer :: subcell, i
        
        do subcell = 1, thisOctal%maxChildren
            if (thisOctal%hasChild(subcell)) then
                ! find the child
                do i = 1, thisOctal%nChildren, 1
                    if (thisOctal%indexChild(i) == subcell) then
                        child => thisOctal%child(i)
                        call setTissues(child)
                        exit
                    end if
                end do
            else

               if (.not.associated(thisOctal%localInterfaces)) then
                  write(*,*) "local interfaces not assocaited for octal of depth ",thisoctal%nDepth
                  write(*,*) "check parent ",associated(thisOctal%parent%localinterfaces)
               endif
               if (.not.associated(thisOctal%iTissue)) allocate(thisOctal%itissue(1:thisOctal%maxChildren))
               thisOctal%iTissue(Subcell) = &
                    findGlobalObjectId( thisOctal%localInterfaces(subcell), subcellCentre(thisOctal,subcell))
            endif
        enddo
    end subroutine setTissues

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
                thisOctal%etaline(subcell) = thisOctal%etaline(subcell)/cellVolume(thisOctal,Subcell)
            endif
        enddo
    end subroutine calculateEnergyDensity


    subroutine setupSpectrum(spectrum, teff)
      use inputs_mod, only : fibreSpectrum
      real(double) :: teff
      type(SPECTRUMTYPE) :: spectrum
      integer :: nLambda
      real(double) :: lamStart, lamEnd
      logical :: ok

      lamStart = 4000.
      lamEnd = 8000.
      nLambda = 1000

      select case (fibreSpectrum)
         case("blackbody")
            call fillSpectrumBB(spectrum, teff, lamStart, lamEnd, nLambda)
         case("led")
            call readSpectrum(spectrum, "led.dat", ok)
      end select
    end subroutine setupSpectrum



#ifdef MPI

  subroutine updateGridMPI(grid)

    use mpi
    use amr_utils_mod
    implicit none

    type(gridtype) :: grid
    integer :: nOctals, nVoxels
    real, allocatable :: tempRealArray(:)
    real, allocatable :: nCrossings(:)
    real(double), allocatable :: distanceGrid(:), etaLine(:), tempDoubleArray(:)
    integer :: ierr, nIndex

    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
    nOctals = 0
    nVoxels = 0
    call countVoxels(grid%octreeRoot,nOctals,nVoxels)
    allocate(nCrossings(1:nVoxels))
    allocate(distanceGrid(1:nVoxels))
    allocate(etaLine(1:nVoxels))

    nIndex = 0
    call packValues(grid%octreeRoot,nIndex, &
         distanceGrid,nCrossings, etaLine)

    allocate(tempRealArray(nVoxels))
    allocate(tempDoubleArray(nVoxels))


    tempDoubleArray = 0.d0
    call MPI_ALLREDUCE(distanceGrid,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    distanceGrid = tempDoubleArray 

    tempDoubleArray = 0.d0
    call MPI_ALLREDUCE(etaLine,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    etaLine = tempDoubleArray 

    tempRealArray = 0.0
    call MPI_ALLREDUCE(nCrossings,tempRealArray,nVoxels,MPI_REAL,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    nCrossings = tempRealArray 

    
    deallocate(tempRealArray, tempDoubleArray)
     
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
    
    nIndex = 0
    call unpackValues(grid%octreeRoot,nIndex, &
         distanceGrid,nCrossings, etaLine)

    deallocate(nCrossings, distanceGrid, etaLine)

  end subroutine updateGridMPI

  recursive subroutine packvalues(thisOctal,nIndex, &
       distanceGrid, nCrossings, etaLine)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: distanceGrid(:)
  real(double) :: etaLine(:)
  real :: nCrossings(:)
  integer :: nIndex
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call packvalues(child,nIndex,distanceGrid,nCrossings,etaLine)
                exit
             end if
          end do
       else
          nIndex = nIndex + 1
          distanceGrid(nIndex) = thisOctal%distanceGrid(subcell)
          etaLine(nIndex) = thisOctal%etaLine(subcell)
          nCrossings(nIndex) = real(thisOctal%nCrossings(subcell))
       endif
    enddo
  end subroutine packvalues

  recursive subroutine unpackvalues(thisOctal,nIndex,distanceGrid,nCrossings, etaLine)
!    use inputs_mod, only : storeScattered
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: distanceGrid(:), etaLine(:)
  real :: ncrossings(:)
  integer :: nIndex
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unpackvalues(child,nIndex,distanceGrid,nCrossings, &
                     etaLine)
                exit
             end if
          end do
       else
          nIndex = nIndex + 1
          thisOctal%distanceGrid(subcell) = distanceGrid(nIndex)
          thisOctal%etaLine(subcell) = etaLine(nIndex)
          thisOctal%nCrossings(subcell) = int(nCrossings(nIndex))
       endif
    enddo
  end subroutine unpackvalues

  subroutine gatherDetector(thisDetector)
    use mpi
    type(DETECTORTYPE) :: thisDetector
    real(double), allocatable :: tArrayd(:), tempArrayd(:)
    integer :: n, ierr

    n = thisDetector%nx * thisDetector%ny
    allocate(tArrayd(1:n))
    allocate(tempArrayd(1:n))

    tArrayd = RESHAPE(thisDetector%intensity, SHAPE=SHAPE(tArrayd))
    call MPI_ALLREDUCE(tArrayd,tempArrayd,n,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    thisDetector%intensity = RESHAPE(tempArrayd,SHAPE=SHAPE(thisDetector%intensity))

    tArrayd = RESHAPE(thisDetector%rintensity, SHAPE=SHAPE(tArrayd))
    call MPI_ALLREDUCE(tArrayd,tempArrayd,n,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    thisDetector%rintensity = RESHAPE(tempArrayd,SHAPE=SHAPE(thisDetector%intensity))

    tArrayd = RESHAPE(thisDetector%gintensity, SHAPE=SHAPE(tArrayd))
    call MPI_ALLREDUCE(tArrayd,tempArrayd,n,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    thisDetector%gintensity = RESHAPE(tempArrayd,SHAPE=SHAPE(thisDetector%intensity))

    tArrayd = RESHAPE(thisDetector%bintensity, SHAPE=SHAPE(tArrayd))
    call MPI_ALLREDUCE(tArrayd,tempArrayd,n,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    thisDetector%bintensity = RESHAPE(tempArrayd,SHAPE=SHAPE(thisDetector%intensity))

    deallocate(tArrayd, tempArrayd)

  end subroutine gatherDetector


#endif

end module biophysics_mod


! =============================================================
! =============================================================
! =======  CODES FOR BUILDING INTERFACES FROM FUNCTIONS   =====
! =============================================================
! =============================================================
!
!
!  real(double) function skinInterface(x, y, Zk, Akx, Aky, omegakx, omegaky, phikx, phiky) result (z)
!    real(double) :: x, y, Akx, Aky, omegakx, omegaky, phikx, phiky, zk
!    z = zk + (Akx * sin (omegakx * x + phikx)) + (Aky * sin (omegaky * y + phiky))
!!    z = zk
!  end function skinInterface
!
!  logical function splitSkin(thisOctal, subcell)
!    type(OCTAL) :: thisOctal
!    integer :: subcell
!    type(VECTOR) :: cen, probe(8)
!    real(double) :: x, y, z, d
!    real(double) :: zk(8)
!    real(double) :: Akx(8)
!    real(double) :: Aky(8)
!    real(double) :: omegakx(8)
!    real(double) :: omegaky(8)
!    real(double) :: phikx(8)
!    real(double) :: phiky(8)
!    real(double) :: zLayer
!    integer :: iLayer, iCorner
!
!
!    splitSkin = .false.
!
!    Akx = (/2.d0, 2.5d0, 20d0, 2.d0, 2.d0, 2.d0, 5.d0, 5.d0/)
!    Aky = Akx
!
!    omegakx = (/100.d0, 80.d0, 50.d0, 20.d0, 20.d0, 20.d0, 20.d0, 25.d0/)
!    omegaky = (/150.d0, 80.d0, 45.d0, 40.d0, 50.d0, 50.d0, 50.d0, 30.d0/)
!    omegakx = omegakx/1.d4
!    omegaky = omegaky/1.d4
!
!    omegakx = pi/omegakx
!    omegaky = pi/omegaky
!
!    zk = (/0.d0, 20.d0, 100.d0, 200.d0, 280.d0, 1900.d0, 2100.d0, 8000.d0/)
!
!    Akx = Akx / 1.d4
!    Aky = Aky / 1.d4
!
!    zk = zk / 1.d4
!
!    phikx = 0.d0
!    phiky = 0.d0
!
!    d = thisOctal%subcellSize / 2.d0
!
!    probe(1) = vector(+1.d0, +1.d0, 1.0)
!    probe(2) = vector(+1.d0, -1.d0, 1.0)
!    probe(3) = vector(-1.d0, +1.d0, 1.0)
!    probe(4) = vector(-1.d0, -1.d0, 1.0)
!    probe(5) = vector(+1.d0, +1.d0,-1.0)
!    probe(6) = vector(+1.d0, -1.d0,-1.0)
!    probe(7) = vector(-1.d0, +1.d0,-1.0)
!    probe(8) = vector(-1.d0, -1.d0,-1.0)
!
!    do iCorner = 1, 8
!       cen = subcellCentre(thisOctal,subcell) + probe(iCorner)*d
!       x = cen%x
!       z = cen%z
!       y = cen%y
!
!
!       do iLayer = 1, 7
!          zLayer = skinInterface(x, y, Zk(iLayer), Akx(iLayer), Aky(iLayer), omegakx(iLayer), &
!               omegaky(iLayer), phikx(iLayer), phiky(iLayer))
!          if (abs(zLayer-z)< d) splitSkin =.true.
!       enddo
!    enddo
!  end function splitSkin
!
!
!  recursive subroutine fillMatcherSkinLayers(thisOctal)
!    use amr_utils_mod, only : inSubcell
!    use inputs_mod, only : smallestCellsize
!    type(octal), pointer   :: thisOctal
!    type(octal), pointer  :: child
!    integer :: subcell, i, iLayer
!    real(double) :: x, y, z, d
!    real(double) :: zk(8)
!    real(double) :: Akx(8)
!    real(double) :: Aky(8)
!    real(double) :: omegakx(8)
!    real(double) :: omegaky(8)
!    real(double) :: phikx(8)
!    real(double) :: phiky(8)
!    real(double) :: zUpper, zLower, zLayer
!    type(VECTOR) :: cen
!
!    Akx = (/2.d0, 2.5d0, 20d0, 2.d0, 2.d0, 2.d0, 5.d0, 5.d0/)
!    Aky = Akx
!
!    omegakx = (/100.d0, 80.d0, 50.d0, 20.d0, 20.d0, 20.d0, 20.d0, 25.d0/)
!    omegaky = (/150.d0, 80.d0, 45.d0, 40.d0, 50.d0, 50.d0, 50.d0, 30.d0/)
!    omegakx = omegakx/1.d4
!    omegaky = omegaky/1.d4
!
!    omegakx = pi/omegakx
!    omegaky = pi/omegaky
!
!    zk = (/0.d0, 20.d0, 100.d0, 200.d0, 280.d0, 1900.d0, 2100.d0, 8000.d0/)
!
!    Akx = Akx / 1.d4
!    Aky = Aky / 1.d4
!
!    zk = zk / 1.d4
!
!    phikx = 0.d0
!    phiky = 0.d0
!    do subcell = 1, thisOctal%maxChildren
!       if (thisOctal%hasChild(subcell)) then
!          ! find the child
!          do i = 1, thisOctal%nChildren, 1
!             if (thisOctal%indexChild(i) == subcell) then
!                child => thisOctal%child(i)
!                call fillMatcherSkinLayers(child)
!                exit
!             end if
!          end do
!       else
!          cen = subcellCentre(thisOctal, subcell)
!          thisOctal%surfaceNormal(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
!          x = cen%x
!          y = cen%y
!          z = cen%z
!          d = thisOctal%subcellSize/2.d0 + 0.01d0*smallestCellSize
!          if (.not.associated(thisOctal%iTissue)) allocate(thisOctal%iTissue(1:thisOctal%maxChildren))
!
!          do iLayer = 1, 7
!             zLower = skinInterface(x, y, Zk(iLayer), Akx(iLayer), Aky(iLayer), omegakx(iLayer), &
!                  omegaky(iLayer), phikx(iLayer), phiky(iLayer))
!             zUpper = skinInterface(x, y, Zk(iLayer+1), Akx(iLayer+1), Aky(iLayer+1), omegakx(iLayer+1), &
!                  omegaky(iLayer+1), phikx(iLayer+1), phiky(iLayer+1))
!
!
!             if ((z < zUpper).and.(z > zLower)) then
!                thisOctal%rho(subcell) = real(iLayer)
!                thisOctal%iTissue(subcell) = iLayer
!             endif
!             zLayer = skinInterface(x, y, Zk(iLayer), Akx(iLayer), Aky(iLayer), omegakx(iLayer), &
!                  omegaky(iLayer), phikx(iLayer), phiky(iLayer))
!
!             if (inSubcell(thisOctal,subcell,VECTOR(x,y,zlayer))) then
!                thisOctal%iTissue(subcell) = iLayer
!             endif
!
!          enddo
!
!      endif
!    enddo
!  end subroutine fillMatcherSkinLayers
!
!  recursive subroutine setSurfaceNormals(grid,thisOctal)
!    use inputs_mod, only : smallestCellsize
!    use amr_utils_mod
!    type(octal), pointer   :: thisOctal, neighbourOctal
!    type(octal), pointer  :: child
!    type(GRIDTYPE) :: grid
!    integer :: subcell, i, iLayer, neighbourSubcell
!    real(double) :: x, y, z, d
!    real(double) :: zk(8)
!    real(double) :: Akx(8)
!    real(double) :: Aky(8)
!    real(double) :: omegakx(8)
!    real(double) :: omegaky(8)
!    real(double) :: phikx(8)
!    real(double) :: phiky(8)
!    real(double) :: zLayer, xp, yp, zp
!
!    type(VECTOR) :: cen, vec
!
!    Akx = (/2.d0, 2.5d0, 20d0, 2.d0, 2.d0, 2.d0, 5.d0, 5.d0/)
!    Aky = Akx
!
!    omegakx = (/100.d0, 80.d0, 50.d0, 20.d0, 20.d0, 20.d0, 20.d0, 25.d0/)
!    omegaky = (/150.d0, 80.d0, 45.d0, 40.d0, 50.d0, 50.d0, 50.d0, 30.d0/)
!    omegakx = omegakx/1.d4
!    omegaky = omegaky/1.d4
!
!    omegakx = pi/omegakx
!    omegaky = pi/omegaky
!
!    zk = (/0.d0, 20.d0, 100.d0, 200.d0, 280.d0, 1900.d0, 2100.d0, 8000.d0/)
!
!    Akx = Akx / 1.d4
!    Aky = Aky / 1.d4
!
!    zk = zk / 1.d4
!
!    phikx = 0.d0
!    phiky = 0.d0
!    do subcell = 1, thisOctal%maxChildren
!       if (thisOctal%hasChild(subcell)) then
!          ! find the child
!          do i = 1, thisOctal%nChildren, 1
!             if (thisOctal%indexChild(i) == subcell) then
!                child => thisOctal%child(i)
!                call setSurfaceNormals(grid,child)
!                exit
!             end if
!          end do
!       else
!          cen = subcellCentre(thisOctal, subcell)
!          x = cen%x
!          y = cen%y
!          z = cen%z
!          d = thisOctal%subcellSize/2.d0 + 0.01d0*smallestCellSize
!
!          do iLayer = 1, 7
!             zLayer = skinInterface(x, y, Zk(iLayer), Akx(iLayer), Aky(iLayer), omegakx(iLayer), &
!                  omegaky(iLayer), phikx(iLayer), phiky(iLayer))
!
!             if (inSubcell(thisOctal,subcell,VECTOR(x,y,zlayer))) then
!                thisOctal%surfaceNormal(subcell) = VECTOR(0.d0, 0.d0, 1.d0)
!                xp = x
!                yp = y
!                zp = z - d
!                vec = VECTOR(xp, yp, zp)
!                if (inOctal(grid%octreeRoot, vec)) then
!                   neighbourOctal => thisOctal
!                   call findSubcellLocal(vec, neighbourOctal, neighbourSubcell)
!                   neighbourOctal%surfaceNormal(neighbourSubcell) = VECTOR(0.d0, 0.d0, -1.d0)
!                endif
!             endif
!          enddo
!
!      endif
!    enddo
!  end subroutine setSurfaceNormals
!
! =============================================================

