  subroutine setKappaTest(grid, scale, aMin, aMax, a0, qDist, pDist, grainType, &
       ngrain, abundance, grainname, lambdaTau)
    use mieDistCrossSection_mod, only: mieDistCrossSection

    implicit none
    type(GRIDTYPE) :: grid
    real :: aMin, aMax, a0, qDist, pDist
    real :: sigmaAbs, sigmaSca, sigmaExt
    real :: scale
    real, allocatable :: mReal(:), mImg(:)
    real, allocatable :: mReal2D(:,:), mImg2D(:,:)  ! size = ngrain x nlambda
    character(len=*) :: grainType
    integer, intent(in) :: ngrain  ! number of grain types
    real, intent(in) :: abundance(ngrain)   ! relative abundance of grains
    character(len=*) :: grainname(ngrain)   ! names of grains available
    real :: sig_ext, sig_scat, sig_abs
    real :: total_abundance, gsca
    real :: meanParticleMass
    real :: lambdaTau
    real :: getMeanMass2

    integer :: i, j 

    scale = 1.
    allocate(mReal(1:grid%nLambda))
    allocate(mImg(1:grid%nLambda))
    call locate(grid%lamArray, grid%nLambda, lambdaTau, i)


    if (writeoutput) write(*,*) "kappa test set for: ",grid%lamarray(i)


    scale = 1.
    if (writeoutput) write(*,'(a)') "NEW: Filling grid with mie cross-sections..."

    if (graintype(1:5) == "mixed") then
       ! Synthetic grains

       ! quick test for zero total abundance.
       total_abundance = SUM(abundance)
       if ( total_abundance <= 0.0 ) then
          write(*,*) "Error:: total_abundance <= 0.0 in  grain_mod::setKappaTest."
          write(*,*) "  ==> You probably forgot to assign abundance in your parameter file!"
          write(*,*) "  ==> Exiting the prograim ... "
          stop 
       end if

       ! allocate mem for temp arrays
       allocate(mReal2D(1:ngrain, 1:grid%nLambda))
       allocate(mImg2D(1:ngrain, 1:grid%nLambda))
       ! initializing the values
       mReal2D(:,:) = 0.0; mImg2D(:,:) = 0.0

       ! Find the index of refractions for all types of grains available
       do j = 1, ngrain
          call getRefractiveIndex(grid%lamArray, grid%nLambda, grainname(j), mReal, mImg)
          mReal2D(j,:) = mReal(:)  ! copying the values to a 2D maxtrix
          mImg2D(j,:)  = mImg(:)   ! copying the values to a 2D maxtrix            
       end do

       ! finding the cross sections
       sigmaExt = 0.0; sigmaAbs=0.0; sigmaSca=0.0 ! initializing the values

       do j = 1, ngrain
          call mieDistCrossSection(aMin, aMax, a0, qDist, pDist, grid%lamArray(i), &
               mReal2D(j,i), mImg2D(j,i), sig_ext, sig_scat, sig_abs ,gsca)

          ! Weighting the cross section according to their abundance...            
          sigmaExt = sig_ext*abundance(j)+ sigmaExt
          sigmaAbs = sig_abs*abundance(j)+ sigmaAbs
          sigmaSca = sig_scat*abundance(j)+ sigmaSca
       end do
       sigmaExt =    sigmaExt/total_abundance 
       sigmaAbs =    sigmaAbs/total_abundance 
       sigmaSca =    sigmaSca/total_abundance 

    else 
       ! Do a single grain calculations...       
       call getRefractiveIndex(grid%lamArray, grid%nLambda, graintype, mReal, mImg)         
       call mieDistCrossSection(aMin, aMax, a0, qDist, pDist,grid%lamArray(i),  &
            mReal(i), mImg(i), sigmaExt, sigmaSca, sigmaAbs, gsca)
    end if
    meanParticleMass = getMeanMass2(aMin, aMax, a0, qDist, pDist, graintype)

    grid%kappaTest = sigmaExt * 1.e10 / meanParticleMass

    if (writeoutput) write(*,*) "kappa test is: ",grid%kappatest
  end subroutine setKappaTest

  subroutine setKappa(kappaAbs, kappaSca, lambda, nLambda, aMin, aMax, a0, qDist, pDist, grainType)

    implicit none
    real :: aMin, aMax, qDist, a0, pDist
    real :: sigmaAbs, sigmaSca, sigmaExt
    real, allocatable :: mReal(:), mImg(:)
    character(len=*) :: grainType
    real :: lambda(:), kappaAbs(:), kappaSca(:), gSca
    integer :: nLambda
    integer :: i

    if (writeoutput) write(*,'(a)') "Setting kappas for: ",trim(grainType)

    allocate(mReal(1:nLambda))
    allocate(mImg(1:nLambda))
    call getRefractiveIndex(lambda, nLambda, graintype, mReal, mImg)
    do i = 1, nLambda
       call mieDistCrossSection(aMin, aMax, a0, qDist, pDist, lambda(i),  mReal(i), mImg(i), sigmaExt, &
            sigmaSca, sigmaAbs, gSca)
       kappaAbs(i) = sigmaAbs * 1.e10
       kappaSca(i) = sigmaSca * 1.e10
    enddo
  end subroutine setKappa

  !
  ! Computes Kappa at a single fgrequency/wavelength.
  !
  subroutine MieCrossSection(sigmaExt, sigmaAbs, sigmaSca, &
       aMin, aMax, a0, qDist, pDist, grainType, &
       ngrain, abundance, grainname, lambda)
    use mieDistCrossSection_mod, only: mieDistCrossSection

    implicit none
    real, intent(out) :: sigmaExt, sigmaAbs, sigmaSca ! total, absorption and scattering  x-sections
    real, intent(in) :: aMin, aMax, a0, qDist, pDist
    character(len=*), intent(in) :: grainType
    integer, intent(in) :: ngrain  ! number of grain types
    real, intent(in) :: abundance(ngrain)   ! relative abundance of grains
    character(len=*) :: grainname(ngrain)   ! names of grains available
    real, intent(in) :: lambda  ! wavelength at which sigma are cmoputed
    !
    real :: mReal(ngrain), mImg(ngrain)  ! automatic arrays
    !
    real :: sig_ext, sig_scat, sig_abs
    real :: total_abundance, gsca
!    real :: meanParticleMass
!    real :: getMeanMass2
    integer, parameter :: nlambda = 1
    integer :: j 
    ! dummy array for interface with subroutines
    real :: lamArray(nlambda)   
    real :: mRealArray(nlambda), mImgArray(nlambda)  ! automatic arrays


    lamArray(:) = lambda  ! same for all elemets
    if (graintype(1:5) == "mixed") then
       ! Synthetic grains

       ! quick test for zero total abundance.
       total_abundance = SUM(abundance)
       if ( total_abundance <= 0.0 ) then
          write(*,*) "Error:: total_abundance <= 0.0 in  dust_mod::getKappa."
          write(*,*) "  ==> You probably forgot to assign abundance in your parameter file!"
          write(*,*) "  ==> Exiting the prograim ... "
          stop 
       end if

       ! initializing the values
       mReal(:) = 0.0; mImg(:) = 0.0

       ! Find the index of refractions for all types of grains available
       do j = 1, ngrain
          call getRefractiveIndex(lamArray, nLambda, grainname(j), mRealArray, mImgArray)
          mReal(j) = mRealArray(1) 
          mImg(j)  = mImgArray(1)  
       end do

       ! finding the cross sections
       sigmaExt = 0.0; sigmaAbs=0.0; sigmaSca=0.0 ! initializing the values

       do j = 1, ngrain
          call mieDistCrossSection(aMin, aMax, a0, qDist, pDist, lamArray(1), &
               mReal(j), mImg(j), sig_ext, sig_scat, sig_abs ,gsca)

          ! Weighting the cross section according to their abundance...            
          sigmaExt = sig_ext*abundance(j)+ sigmaExt
          sigmaAbs = sig_abs*abundance(j)+ sigmaAbs
          sigmaSca = sig_scat*abundance(j)+ sigmaSca
       end do
       sigmaExt =    sigmaExt/total_abundance 
       sigmaAbs =    sigmaAbs/total_abundance 
       sigmaSca =    sigmaSca/total_abundance 

    else 
       ! Do a single grain calculations... 
       call getRefractiveIndex(lamArray, nLambda, graintype, mRealArray, mImgArray)         
       call mieDistCrossSection(aMin, aMax, a0, qDist, pDist,lamArray(1),  &
            mRealArray(1), mImgArray(1), sigmaExt, sigmaSca, sigmaAbs, gsca)
    end if

    !      meanParticleMass = getMeanMass2(aMin, aMax, a0, qDist, pDist, graintype)
  end subroutine MieCrossSection

  recursive subroutine fillDustABAur(grid, thisOctal)

    use input_variables, only : rInner, rOuter
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    type(VECTOR) :: rVec
    real :: x, z
    real :: height
    real(double) :: fac
    integer :: subcell, i
    height = 0.d0

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fillDustABAur(grid, child)
                exit
             end if
          end do
       else

          thisOctal%DustTypeFraction(subcell,1:2) = 0.d0
          thisOctal%DustTypeFraction(subcell,1) = 1.d0
          rVec = subcellCentre(thisOctal, subcell)
          x = rVec%x
          z = rVec%z
          if ( (x > rInner).and.(x < rOuter)) then
             call returnScaleHeight(grid, x, height)
             fac = exp(-abs(z/height))
             thisOctal%dustTypeFraction(subcell,1) = 1.d0 - fac
             thisOctal%dustTypeFraction(subcell,2) = fac
          endif

       end if
    end do

  end subroutine fillDustABAur

  recursive subroutine fillDustWhitney(grid, thisOctal)
    use input_variables
    TYPE(VECTOR) :: point
    TYPE(gridtype), INTENT(IN)    :: grid
    real :: r, mu, mu_0, muCavity, rhoEnv, r_c
    real :: h, rhoDisc, alpha
    real(double) :: fac
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fillDustWhitney(grid, child)
                exit
             end if
          end do
       else

          point = subcellCentre(thisOctal, subcell)
          muCavity = cos(cavAngle/2.)
          r = modulus(point)*1.e10

          mu = (point%z*1.e10) /r

          r_c = erInner
          alpha = 2.25
          beta = 1.25

          rhoEnv = cavdens * mHydrogen

          ! by default use envelope grains

          thisOctal%dustTypeFraction(subcell,1:4) = 0.d0
          thisOctal%dustTypeFraction(subcell,3) = 1.d0


          if ((r > erInner).and.(r < erOuter)) then
             mu_0 = rtnewtdust(-0.2 , 1.5 , 1.e-4, r/r_c, abs(mu))
             ! equation 1 for Whitney 2003 ApJ 591 1049 has a mistake in it
             ! this is from Momose et al. 1998 ApJ 504 314

             rhoEnv = (mdotenv / fourPi) * (bigG * mCore)**(-0.5) * r**(-1.5) * &
                  (1. + abs(mu)/mu_0)**(-0.5) * &
                  (abs(mu)/mu_0 + (2.*mu_0**2 * r_c/r))**(-1.)

             fac =  1.d0-min(dble(r - erInner)/(0.02d0*erinner),1.d0)
             fac = exp(-fac*10.d0)
             rhoEnv = rhoEnv * fac
             rhoEnv = max(rhoEnv, tiny(rhoEnv))
          endif

          if (mu_0 > muCavity) then
             !       fac =  1.d0-min(dble(r - drInner)/(0.02d0*erinner),1.d0)
             !       fac = exp(-fac*10.d0)
             rhoEnv = cavdens * 2.*mHydrogen

             ! outflow cavity grains

             thisOctal%dustTypeFraction(subcell,1:4) = 0.d0
             thisOctal%dustTypeFraction(subcell,4) = 1.d0
          endif

          if (r < erInner) then
             thisOctal%dustTypeFraction(subcell,1:4) = 0.d0
             thisOctal%dustTypeFraction(subcell,4) = 1.d0
          endif

          rho0  = mDisc *(beta-alpha+2.) / ( twoPi**1.5 * 0.01*rStellar * rStellar**(alpha-beta) * ( &
               (drouter**(beta-alpha+2.)-drInner**(beta-alpha+2.))) )

          r = sqrt(point%x**2 + point%y**2)*1.e10
          h = 0.01 * rStellar * (r/rStellar)**beta
          rhoDisc = 1.e-30
          if ((r > drInner).and.(r < drOuter)) then
             rhoDisc = rho0 * (rStellar/r)**alpha  * exp(-0.5*((point%z*1.e10)/h)**2)
             fac =  1.d0-min(dble(r - drInner)/(0.02d0*drinner),1.d0)
             fac = exp(-fac*10.d0)
             rhoDisc = rhoDisc * fac
             rhoDisc = max(rhoDisc, tiny(rhoDisc))

             if (rhoDisc > 1.e-9) then ! large midplane dust grains
                thisOctal%dusttypeFraction(subcell,1:4) = 0.d0
                thisOctal%dusttypeFraction(subcell,1) = 1.d0
             else if (rhoDisc > rhoEnv) then  ! normal disc dust
                thisOctal%dusttypeFraction(subcell,1:4) = 0.d0
                thisOctal%dusttypeFraction(subcell,2) = 1.d0
             endif
          endif


       endif
    enddo

  end subroutine fillDustWhitney

