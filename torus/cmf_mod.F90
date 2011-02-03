module cmf_mod

  ! written by tjh


  use kind_mod
  use constants_mod
  use messages_mod
  use vector_mod
  use octal_mod, only: octal, octalWrapper, subcellCentre
  use amr_mod, only: inOctal, distanceToCellBoundary, findsubcellLocal, amrGridVelocity
  use gridtype_mod, only: GRIDTYPE
  use utils_mod, only: gaussj, locate, toPerAngstrom
  use random_mod
  use modelatom_mod, only: MODELATOM, BoltzSahaGeneral, bfOpacity, bfEmissivity, photoCrossSection, collisionRate, &
        addcrosssectionstoatom, createcontfreqarray, createrbbarrays, returneinsteincoeffs
  use source_mod, only: SOURCETYPE, I_nu, insideSource, distanceToSource
  use surface_mod, only: getElement
  use timing, only: tune
  use parallel_mod, only: torus_mpi_barrier
  use vtk_mod, only: writeVtkfile
  use octal_mod, only : allocateAttribute
  implicit none

  private
  public:: atomLoop, calculateAtomSpectrum

contains

  subroutine solveLevels(thisOctal, subcell, nPops, jnuLine,  &
       temperature, nAtom, thisAtom, ne, rho, jnuCont, freq, dfreq, nfreq)
    use input_variables, only : debug
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer :: nFreq
    integer :: num(100)
    real(double) :: freq(:), dfreq(:), jnuCont(:)
    real(double), allocatable :: hCgsFreq(:)
    real(double) :: nPops(:,:)
    real(double) :: temperature, ne
    real(double) :: jnuLine(:)
    type(MODELATOM) :: thisAtom(:)
    real(double), allocatable :: matrixA(:,:), matrixB(:)
!    integer, allocatable :: indx(:)
!    real(double), allocatable :: vMatrix(:,:), wMatrix(:,:), xMatrix(:)
!    real(double) :: wMax, wMin
!    real(double) :: d 
    real(double) :: boltzFac
    integer :: nLevels
    integer :: i
    integer :: itrans, l, k
    real(double) :: collEx, colldeEx
    real(double) :: a, Bul, Blu
    real(double) :: photoRatelk, recombRatekl, xSection, xSection2
    real(double) :: rho
    real(double) :: NstarRatio, totRecomb, totPHotoIon, totcion
    real(double) :: tot1, tot2
    real(double) :: nuStart !,thisFreq,test
    integer :: iStart
    integer :: iJnu
    real(double) :: jnu
    real(double) :: radtot, colltot
    integer :: nAtom, iAtom, nMatrix
    integer, allocatable :: nOffset(:)
    logical, allocatable :: continuumGround(:)
    integer, allocatable :: nCons(:)
    logical :: ok
    integer :: nDifferentElements
    
    blu = 0.d0; bul = 0.d0; a = 0.d0; ok = .true.
    nDifferentElements = 0
    do iAtom = 1, nAtom
       if (iAtom /= nAtom) then
          if (thisAtom(iAtom)%nz /= thisAtom(iAtom+1)%nz) then ! if next atom is different element
             nDifferentElements = nDifferentElements + 1
          else
             nDifferentElements = nDifferentElements + 1
          endif
       endif
    enddo

    allocate(nOffset(1:nAtom))
    allocate(nCons(1:nAtom))
    allocate(continuumGround(1:nAtom))
    allocate(hCgsFreq(1:nFreq))
    hCgsFreq(1:nFreq) = hCgs * freq(1:nFreq)
    continuumGround = .false.

    nMatrix = 0
    nOffset(1) = 0
    nCons = 0
    do iAtom = 1, nAtom
       if (iAtom /= nAtom) then
          if (thisAtom(iAtom)%nz /= thisAtom(iAtom+1)%nz) then ! if next atom is different element
             nMatrix = nMatrix + thisAtom(iAtom)%nLevels 
             nMatrix = nMatrix + 1  
             nCons(iAtom) = nMatrix
             nOffset(iAtom+1) = nMatrix
          else
             nMatrix = nMatrix + thisAtom(iAtom)%nLevels-1  ! continuum state of this atom is ground state of next
             continuumGround(iAtom) = .true.
             nOffset(iAtom+1) = nMatrix
             nCons(iAtom) = -1
          endif
       else
          nMatrix = nMatrix + thisAtom(iAtom)%nLevels  ! last atom
          nMatrix = nMatrix + 1  
          nCons(iAtom) = nMatrix
       endif
    enddo


    allocate(matrixA(1:nMatrix,1:nMatrix))
    allocate(matrixB(1:nMatrix))


    matrixA = 0.d0
    matrixB = 0.d0



    do i = 1, nMatrix
       num(i) = i
    enddo


    do iAtom = 1, nAtom
       if (nCons(iAtom) /= 0) then
          if (.not.continuumGround(iAtom)) then
             matrixA(nCons(iAtom),1+nOffset(iAtom):nOffset(iAtom)+thisAtom(iAtom)%nLevels) = 1.d0
             matrixB(nCons(iAtom)) = thisOctal%atomAbundance(subcell,iAtom) * rho
          else
             do i = iAtom, nAtom
                if (nCons(i) /= -1) exit
             enddo
             matrixA(nCons(i),1+nOffset(iAtom): nOffset(iAtom)+thisAtom(iAtom)%nLevels-1) = 1.d0
          endif
       endif
    enddo


    do iAtom = 1, nAtom

       colltot = 0.d0
       radtot = 0.d0

       nLevels = thisAtom(iAtom)%nLevels

       ! recombination rates

       do iTrans = 1, thisAtom(iAtom)%nTrans
          totRecomb = 0.d0
          tot1 = 0.d0
          tot2 = 0.d0

          k = thisAtom(iAtom)%iUpper(iTrans)
          l = thisAtom(iAtom)%iLower(iTrans)

          NstarRatio = boltzSahaGeneral(thisAtom(iAtom),  l, ne, temperature)

          if ((thisAtom(iAtom)%transType(iTrans) == "RBF").and.&
               (.not.thisAtom(iAtom)%inDetailedBalance(iTrans)))  then ! radiative recomb

             recombratekl = 0.d0
             nuStart = (thisAtom(iAtom)%iPot - thisAtom(iAtom)%energy(l))*evtoerg/hcgs

             call locate(freq, nfreq, nuStart, istart)

             do i = istart+1, nFreq
                xSection = photoCrossSection(thisAtom(iAtom), l, freq(i-1))
                xSection2 = photoCrossSection(thisAtom(iAtom), l, freq(i))

                recombRatekl = recombRatekl + &
                     0.5 * ( (fourPi/(hCgs*freq(i-1)))*xSection*((2.d0*hCgs*freq(i-1)**3)/cSpeed**2 + jnuCont(i-1)) * &
                     exp(-(hCgs*freq(i-1))/(kErg*temperature)) + &
                             (fourPi/(hCgs*freq(i)))*xSection2*((2.d0*hCgs*freq(i)**3)/cSpeed**2 + jnuCont(i)) * &
                     exp(-(hCgs*freq(i))/(kErg*temperature)) )*dfreq(i)
             enddo

             matrixA(l + noffset(iatom), k+nOffset(iAtom)) = matrixA(l + noffset(iatom), k+nOffset(iAtom)) + &
                  NstarRatio * recombratekl
             matrixA(k+nOffset(iAtom),k+nOffset(iAtom)) = matrixA(k+nOffset(iAtom),k+nOffset(iAtom)) &
                  - NstarRatio * recombratekl 


             tot1 = tot1  + NstarRatio * recombratekl

!             if ((l == 1).and.iatom==2) then
!                write(*,*) "rad recomb to ground state ",nStarRatio * recombratekl * npops(iatom, k+nOffset(iatom))
!             endif


             radtot = radtot + tot1 * ne 

          endif


          if (thisAtom(iAtom)%transType(iTrans) == "CBF") then ! collisional recomb
             collEx = collisionRate(thisAtom(iAtom), iTrans, temperature) 

             collex = tiny(collex)


             matrixA(l + noffset(iatom), k+nOffset(iAtom)) = matrixA(l + noffset(iatom), k+nOffset(iAtom)) + &
                  NstarRatio * collex * ne

             matrixA(k+nOffset(iAtom),k+nOffset(iAtom)) = matrixA(k+nOffset(iAtom),k+nOffset(iAtom)) - NstarRatio*collEx*ne


             tot2 = tot2 + NstarRatio*collEx*ne


             colltot = colltot + tot2 * ne

          
          endif


!          totRecomb =  tot1 + tot2
!          if (thisAtom(iAtom)%transType(iTrans)(2:3) == "BF") then
!             matrixA(k+nOffset(iAtom), k+nOffset(iAtom)) = matrixA(k+nOffset(iAtom), k+nOffset(iAtom))-totRecomb
!          endif


       enddo



       do iTrans = 1, thisAtom(iAtom)%nTrans
          totPhotoIon = 0.d0
          totcion = 0.d0

          k = thisAtom(iAtom)%iUpper(iTrans)
          l = thisAtom(iAtom)%iLower(iTrans)



          if (thisAtom(iAtom)%transType(iTrans) == "RBB") then ! radiative BB rates
             iJnu = thisAtom(iAtom)%indexRBBTrans(iTrans)
             call returnEinsteinCoeffs(thisAtom(iAtom), iTrans, a, Bul, Blu)



             call locate(freq, nfreq, thisAtom(iatom)%transFreq(itrans), istart)

             jnu = jNuLine(iJnu)
!             if ((thisAtom(iAtom)%name == "HeII").and.(l==1).and.(k==2)) write(*,*) "jnu ",jnu,ijnu


             if (.not.thisAtom(iAtom)%inDetailedBalance(iTrans)) then
                

                matrixA(k+nOffset(iAtom),k+nOffset(iAtom)) = matrixA(k+nOffset(iAtom),k+nOffset(iAtom)) - Bul * jnu - a
                matrixA(l+nOffset(iAtom),l+nOffset(iAtom)) = matrixA(l+nOffset(iAtom),l+nOffset(iAtom)) - Blu * jnu
                matrixA(k+nOffset(iAtom),l+nOffset(iAtom)) = matrixA(k+nOffset(iAtom),l+nOffset(iAtom)) + Blu * jnu
                matrixA(l+nOffset(iAtom),k+nOffset(iAtom)) = matrixA(l+nOffset(iAtom),k+nOffset(iAtom)) + Bul * jnu + a


             endif

          endif

          if (thisAtom(iAtom)%transType(iTrans) == "CBB") then ! collision BB  rates


             collEx = collisionRate(thisAtom(iAtom), iTrans, temperature)
             boltzFac =  exp(-abs(thisAtom(iAtom)%energy(k)-thisAtom(iAtom)%energy(l)) / (kev*temperature))
             colldeEx = collEx / (boltzFac * thisAtom(iAtom)%g(k)/thisAtom(iAtom)%g(l))
	     
             matrixA(k+nOffset(iAtom),k+nOffset(iAtom)) = matrixA(k+nOffset(iAtom),k+nOffset(iAtom)) - colldeex * ne
             matrixA(l+nOffset(iAtom),l+nOffset(iAtom)) = matrixA(l+nOffset(iAtom),l+nOffset(iAtom)) - collex * ne
             matrixA(k+nOffset(iAtom),l+nOffset(iAtom)) = matrixA(k+nOffset(iAtom),l+nOffset(iAtom)) + collex * ne
             matrixA(l+nOffset(iAtom),k+nOffset(iAtom)) = matrixA(l+nOffset(iAtom),k+nOffset(iAtom)) + colldeex * ne

          endif

          ! now do bound-free rates



          if (thisAtom(iAtom)%transType(iTrans) == "RBF") then ! photoionization

             if (.not.thisAtom(iAtom)%inDetailedBalance(iTrans)) then


                nuStart = (thisAtom(iAtom)%iPot - thisAtom(iAtom)%energy(l))*evtoerg/hcgs

                call locate(freq, nfreq, nuStart, istart)
                
                photoRatelk = 0.d0
                do i = istart+1, nFreq
                   xSection = photoCrossSection(thisAtom(iAtom), l, freq(i-1))
                   xSection2 = photoCrossSection(thisAtom(iAtom), l, freq(i))
                   photoRatelk = photoRatelk + 0.5d0 * ((jnuCont(i-1)/(hCgsfreq(i-1)))*xSection &
                        + (jnuCont(i)/(hCgsfreq(i)))*xSection2) * dfreq(i)
                enddo
                photoRatelk = photoRatelk * fourPi


                matrixA(l+nOffset(iAtom),l+nOffset(iAtom)) = matrixA(l+nOffset(iAtom),l+nOffset(iAtom)) - photoRatelk
                matrixA(k+nOffset(iAtom),l+nOffset(iAtom)) = matrixA(k+nOffset(iAtom),l+nOffset(iAtom)) + photoRatelk
                totphotoion=totphotoion + photoRatelk

             endif




          endif

          if (thisAtom(iAtom)%transType(iTrans) == "CBF") then ! collisional ionization

             k = thisAtom(iAtom)%iUpper(iTrans)
             l = thisAtom(iAtom)%iLower(iTrans)


             collEx = collisionRate(thisAtom(iAtom), iTrans, temperature)

             matrixA(l+nOffset(iAtom),l+nOffset(iAtom)) = matrixA(l+nOffset(iAtom),l+nOffset(iAtom)) - collex * ne
             matrixA(k+nOffset(iAtom),l+nOffset(iAtom)) = matrixA(k+nOffset(iAtom),l+nOffset(iAtom)) + collex * ne
             totcion = totcion + collex*ne

             
          endif


       enddo


!       if (iatom == 2) write(*,*) "total photoionization: " , totphotoion, totcion

    enddo

!
    if (debug) then
       write(*,'(4x,100i9)') num(1:nMatrix)
       do i = 1, nMatrix
          write(*,'(i4,1p,100e9.1)') i,matrixA(i,1:nMatrix),matrixB(i)
       enddo
    endif

    !    call luSlv(matrixA, matrixB, nMatrix)
    call GAUSSJ(matrixA, nMatrix, nMatrix, matrixB, 1, nMatrix, ok)

    !    allocate(indx(1:nMatrix))
    !    call ludcmp(matrixA, nMatrix, nMatrix, indx, d)
    !    call lubksb(matrixA, nMatrix, nMatrix, indx, matrixB)
    !    deallocate(indx)

    !    allocate(wMatrix(1:nMatrix,1:nMatrix))
    !    allocate(vMatrix(1:nMatrix,1:nMatrix))
    !    allocate(xMatrix(1:nMatrix))
    !    call svdcmp(matrixA, nMatrix, nMatrix, nMatrix, nMatrix, wMatrix, vMatrix)
    !    wMax = maxval(wMatrix)
    !    wMin = 1.d-13 * wMax
    !    where (wMatrix < wMin)
    !       wMatrix = 0.d0
    !    end where
    !    call svbksb(matrixA, Wmatrix, vMatrix, nMatrix, nMatrix, nMatrix, nMatrix, matrixB, xmatrix)
    !    matrixB = xmatrix
    !    deallocate(wmatrix,vmatrix, xmatrix)


    do iAtom = 1, nAtom
       nPops(iAtom,1:thisAtom(iAtom)%nLevels) = matrixB(1+nOffset(iAtom):thisAtom(iatom)%nLevels+nOffset(iAtom))
       if (continuumGround(iatom)) npops(iAtom,thisAtom(iatom)%nlevels) = matrixB(1+nOffset(iatom+1))
!       write(*,*) "iatom ",iatom, continuumGround(iatom)
    enddo
    deallocate(matrixA, matrixB)


!    if (debug) then
!       do iAtom = 1, nAtom
!          write(*,*) trim(thisAtom(iAtom)%name),SUM(nPops(iAtom,1:thisAtom(iAtom)%nLevels-1))
!          write(*,*) trim(thisAtom(iAtom)%name)//"I",nPops(iAtom,thisAtom(iAtom)%nLevels)
!       enddo
!       write(*,*) "Ne",ne,rho/mhydrogen
!    endif
    deallocate(hCgsFreq)
  end subroutine solveLevels


  subroutine getRay(grid, fromOctal, fromSubcell, position, direction, rayDeltaV, ds, phi, i0, &
       hCol, HeIcol, HeIIcol, nAtom, thisAtom, source, nSource, &
       hitPhotosphere, sourceNumber, cosTheta, weightFreq, weightOmega, &
       nRBBTrans, indexRBBTrans, indexAtom, &
       nFreq, freq, iCont)
    use input_variables, only : opticallyThickContinuum, nLambda, mie, onTheSpot
    use amr_mod, only: randomPositionInCell, returnKappa
    use utils_mod, only : findIlambda
    use atom_mod, only : bnu
    type(SOURCETYPE) :: source(:)
    real(double) :: rayDeltaV
    integer :: nfreq
    real(double) :: freq(:)
    real(double) :: iCont(:)
    integer :: nAtom
    integer :: nSource
    integer :: nRBBTrans, indexRBBTrans(:), indexAtom(:)
    type(MODELATOM) :: thisAtom(:)
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, startOctal, fromOctal
    integer :: fromSubcell
    integer :: subcell
    real(double) :: ds, phi, i0(:), r, cosTheta
    real(double) :: Hcol, HeICol, HeIICol
    integer :: iTrans
    type(VECTOR) :: position, direction, currentPosition, thisPosition, thisVel
    type(VECTOR) :: rayVel, startVel, endVel, endPosition, pvec, photoDirection
    real(double) :: alphanu, snu, jnu
    integer :: iLower , iUpper
    real(double) :: dv, deltaV
    integer :: i
    real(double) :: distArray(200), tval
    integer :: nTau
    real(double) :: nLower, nUpper
    real(double) :: dTau, etaline
    real(double), allocatable :: tau(:)
    real(double) :: intensityIntegral
    real(double) :: dvAcrossCell, projVel
    integer :: iCount
    real(double) :: a, blu, bul
    real(double) :: distTosource, totDist
    integer :: sourceNumber
    logical :: hitPhotosphere, hitSource
    real(double) :: dv1, dv2
    real(double) :: weightFreq
    real(double) :: weightOmega
    integer :: iAtom, iRBB
    real(double) :: phiAv, phiNorm
    logical :: firstSubcell
    integer :: nBug
    integer :: iFreq
    integer :: iElement
    integer :: j
    real(double), allocatable :: tauCont(:), jnuCont(:), alphanuCont(:), snuCont(:)
    logical, save :: firstWarning = .true.
    logical, save :: firstTime = .true.
    real(double) :: nStar(5,40)
    logical :: passThroughResonance, velocityIncreasing, ok
    real(double) :: x1, x2, fac, deltaDist, lambda, kappaSca, kappaExt, kappaAbs
    real(double) :: dustOpac, dustEmiss
    Integer :: ilambda
    integer, allocatable,save :: iFreqRBB(:)
    integer, save :: iflagRBB

    !$OMP THREADPRIVATE (firstWarning, firstTime, ifreqRBB, iflagRBB)

    if (firstTime) then
       allocate(iFreqRBB(1:nRBBTrans))
       do iRBB = 1, nRBBTRans
          iAtom = indexAtom(iRBB)
          iTrans = indexRBBTrans(iRBB)
          
          if ((thisAtom(iatom)%name == "HeII").and.(thisAtom(iatom)%ilower(itrans) == 1).and. &
               (thisAtom(iatom)%iupper(itrans) == 2)) iflagRBB = iRbb

          call locate(freq, nfreq, thisAtom(iAtom)%transFreq(iTrans), iFreqRBB(iRBB))
       enddo
       firstTime = .false.
    endif


    distToSource = 0.d0; hitSource = .false.
    a = 0.d0; blu = 0.d0; bul =0.d0
    allocate(tauCont(1:nFreq))
    allocate(jnuCont(1:nFreq))
    allocate(alphanuCont(1:nFreq))
    allocate(snuCont(1:nFreq))

    allocate(tau(1:nRBBTrans))

    position = randomPositionInCell(fromOctal, fromsubcell)

    icount = 1
    thisOctal => grid%octreeRoot
    call findSubcellLocal(position, thisOctal, subcell)

    call randomRayDirection(0.9d0, position, source, nSource, direction, weightOmega)


    call randomNumberGenerator(getDouble=r)

    totDist = 0.d0
    call distanceToSource(source, nSource, position, direction, hitSource, disttoSource, sourcenumber)

    if (hitSource) then
       pVec = (position + (direction * distToSource) - source(sourceNumber)%position)
       call normalize(pVec)
       cosTheta = -1.d0*(pVec.dot.direction)
       photoDirection = pVec
       call normalize(PhotoDirection)
    endif


    rayVel = amrGridVelocity(grid%octreeRoot, position, startOctal = thisOctal, actualSubcell = subcell)


    call randomNumberGenerator(getDouble=r)


    deltaV = 4.3 * thisOctal%microturb(subcell) * (r-0.5d0) ! random frequency near line spectrum peak. 

    rayDeltaV = deltaV
    ! 4.3 corresponds to the width where the peak of the line profile has dropped to 1% of its peak
    ! microturulence is assumed gaussian - b is FULL WIDTH

    weightFreq = phiProf(deltaV, thisOctal%microturb(subcell))

    deltaV = deltaV +  (rayVel .dot. direction)

    projVel = deltaV - (rayVel .dot. direction)

    Hcol = 0.d0
    HeICol = 0.d0
    HeIICol = 0.d0
    i0 = 0.d0
    intensityIntegral = 0.0
    tau = 0.d0
    tauCont = 0.d0
    iCont = 0.d0

    hitPhotosphere = .false.
    firstSubcell = .true.

    currentPosition = position

    nBug = 0
    do while(inOctal(grid%octreeRoot, currentPosition).and.(.not.hitPhotosphere))

       nBug = 0
       phiAv = 0.d0
       phiNorm = 0.d0

       call findSubcellLocal(currentPosition, thisOctal, subcell)
       call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal)

       if ((totDist + tval) > distTosource) then
          tVal = distToSource - totDist
          hitPhotosphere = .true.
       endif

       if (firstSubcell) then
          ds = tval * 1.d10
       endif


       startVel = amrGridVelocity(grid%octreeRoot, currentPosition, startOctal = thisOctal, actualSubcell = subcell) 
       !       startVel = getVel(grid, thisOctal, subcell, currentposition, direction)

       endPosition = currentPosition + tval * direction
       endVel = amrGridVelocity(grid%octreeRoot, endPosition)
       !       endVel = getVel(grid, thisOctal, subcell, endposition, direction)


       dvAcrossCell = ((startVel - rayVel).dot.direction) - ((endVel - rayVel).dot.direction)
       dvAcrossCell = abs(dvAcrossCell / thisOctal%microturb(subcell))
       dv1 = abs(deltaV - (startVel .dot. direction))
       dv2 = abs(deltaV - (endVel .dot. direction))

       distArray(1) = 0.d0
       distArray(2) = tVal
       nTau = 2

       passThroughResonance =.false.


       if (dv1*dv2 < 0.d0) passThroughResonance = .true.

       if (modulus(endVel)==0.d0) passThroughResonance = .false.

       if (passthroughresonance.or.(min(abs(dv1),abs(dv2)) < 4.d0*thisOctal%microturb(subcell))) then

          if (dv1 <= dv2) then
             velocityIncreasing = .true.
          else
             velocityIncreasing = .false.
          endif

          if ( (dv2 - dv1) /= 0.d0) then
             if (velocityIncreasing) then
                x1 = tval*(-4.d0*thisOctal%microturb(subcell) - dv1)/(dv2 - dv1)
                x2 = tval*(+4.d0*thisOctal%microturb(subcell) - dv1)/(dv2 - dv1)
             else
                x1 = tval*(+4.d0*thisOctal%microturb(subcell) - dv1)/(dv2 - dv1)
                x2 = tval*(-4.d0*thisOctal%microturb(subcell) - dv1)/(dv2 - dv1)
             endif
          else
             x1 = 0.d0
             x2 = tVal
          endif

          if (x1 > x2) then
             fac = x1
             x1 = x2
             x2 = fac
          endif
          x1 = max(0.d0, x1)
          x2 = min(x2, tVal)
          nTau = 8
          distArray(1) = 0.d0
          do i = 1, nTau-1
             distArray(i+1) = x1 + (x2 - x1)*dble(i-1)/dble(ntau-2)
          enddo
          distArray(nTau) = tVal
       endif

       if (.not.thisOctal%inflow(subcell)) then
          distArray(1) = 0.d0
          distArray(2) = tVal
          nTau = 2
       endif

       do i = 2, nTau

          distArray(i) = tval * dble(i-1)/dble(nTau-1)

          startOctal => thisOctal
          thisPosition = currentPosition + distArray(i)*direction

          thisVel = amrGridVelocity(grid%octreeRoot, thisPosition, startOctal = startOctal, actualSubcell = subcell) 
          !          thisVel = getVel(grid, thisOctal, subcell, thisposition, direction)


          ! only need to calc continuum opacities once per cell

          if (opticallyThickContinuum.and.(i==2).and.(.not.onTheSpot)) then

             do iAtom = 1, nAtom
                do j = 1, thisAtom(iAtom)%nLevels - 1
                   nStar(iAtom,j) = BoltzSahaGeneral(thisAtom(iAtom),  j, thisOctal%ne(subcell), &
                        dble(thisOctal%temperature(subcell))) * &
                        Thisoctal%atomlevel(subcell, iAtom,thisAtom(iatom)%nLevels)
                enddo
             enddo


             do iFreq = 1, nFreq

                alphanuCont(ifreq) = bfOpacity(freq(ifreq), nAtom, thisAtom, thisOctal%atomLevel(subcell,:,:), &
                     thisOctal%ne(subcell),  dble(thisOctal%temperature(subcell)), ifreq=ifreq)

                jnuCont(iFreq) = bfEmissivity(freq(ifreq), nAtom, thisAtom, &
                     thisOctal%atomLevel(subcell, :, :), nStar, dble(thisOctal%temperature(subcell)), thisOctal%ne(subcell), &
                     ifreq=ifreq)

                dustOpac = 0.d0
                dustEmiss = 0.d0
                if (mie) then
                   lambda = (cSpeed/freq(ifreq)) /angstromTocm
                   iLambda = findIlambda(real(lambda), grid%lamArray, nLambda, ok)
                   call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, kappaSca=kappaSca, kappaAbs=kappaAbs)
                   kappaExt = kappaAbs + kappaSca
                   dustOpac = kappaExt/1.d10
                   dustEmiss = kappaAbs * bnu(freq(ifreq), dble(thisOctal%temperature(subcell)))/1.d10
                endif
                
                jnuCont(ifreq) = jnuCont(ifreq) + dustEmiss
                alphaNuCont(iFreq) = alphaNuCont(ifreq) + dustOpac

                if (alphanuCont(ifreq) /= 0.d0) then
                   snuCont(iFreq) = jnuCont(iFreq)/alphanuCont(iFreq)
                else
                   snuCont(iFreq) = tiny(snuCont(iFreq))
                endif
             enddo
          endif

          if (opticallyThickContinuum.and.(.not.onTheSpot)) then
             deltaDist = (distArray(i)-distArray(i-1)) * 1.d10
             do iFreq = 1, nFreq
                dTau = alphaNuCont(iFreq) *  deltaDist
                iCont(iFreq) = iCont(ifreq) + exp(-tauCont(iFreq)) * (1.d0-exp(-dtau))*snuCont(iFreq)
                if (icont(ifreq) < 0.d0) then
                   write(*,*) "icont negative ",icont(ifreq)
                   write(*,*) "tau ",tauCont(ifreq), " dtau ",dtau, " snu ",snuCont(ifreq)
                   write(*,*) "expression ",exp(-tauCont(iFreq)) * (1.d0-exp(-dtau))*snuCont(iFreq)
                   stop
                endif
                tauCont(iFreq) = tauCont(iFreq) + dtau
             enddo
          endif
          !          write(*,*) tauCont(48),iCont(48),snuCont(48)
          dv = deltaV - (thisVel .dot. direction)

          icount = icount + 1

          do iRBB = 1, nRBBTRans
             iAtom = indexAtom(iRBB)
             iTrans = indexRBBTrans(iRBB)

            

!             call locate(freq, nfreq, thisAtom(iAtom)%transFreq(iTrans), iFreq)
             call returnEinsteinCoeffs(thisAtom(iAtom), iTrans, a, Bul, Blu)


             alphanu = hCgsOverFourPi * phiProf(dv, thisOctal%microturb(subcell)) 


             if (firstSubcell) then
                phiAv = phiAv + phiProf(dv, thisOctal%microturb(subcell)) * &
                     (distArray(i)-distArray(i-1)) * 1.d10
                phiNorm = phiNorm + (distArray(i)-distArray(i-1)) * 1.d10
             endif


             iUpper = thisAtom(iAtom)%iUpper(iTrans)
             iLower = thisAtom(iAtom)%iLower(iTrans)

             nLower = thisOctal%atomLevel(subcell,iAtom, iLower)
             nUpper = thisOctal%atomLevel(subcell,iAtom, iUpper)

             alphanu = alphanu * (nLower * Blu - nUpper * Bul)

             if (opticallyThickContinuum.and.(.not.onTheSpot)) alphanu = alphanu + alphaNuCont(iFreqRBB(iRBB))


             if (alphanu < 0.d0) then
                alphanu = 0.d0
                if (firstWarning) then
                   write(*,*) "negative opacity warning in getray",iUpper,iLower,nLower,nUpper,thisAtom%name
                   firstWarning = .false.
                endif
             endif



             dTau = max(0.d0, alphaNu *  (distArray(i)-distArray(i-1)) * 1.d10)


             etaLine = hCgs * a * thisAtom(iAtom)%transFreq(iTrans)
             etaLine = etaLine * thisOctal%atomLevel(subcell, iAtom, iUpper)
             jnu = (etaLine/fourPi) * phiProf(dv, thisOctal%microturb(subcell)) /thisAtom(iAtom)%transFreq(iTrans)


             if (opticallyThickContinuum.and.(.not.onTheSpot)) jnu = jnu + jnuCont(iFreqRBB(iRBB)) 


             if (alphanu > 1.d-30) then
                snu = jnu/alphanu
             else
                snu = tiny(snu)
             endif


             if (tau(irbb) < 0.d0) then
                write(*,*) "Negative tau(irbb) problem ",tau(irbb), irbb
             endif
             if (dtau < 0.d0) then
                write(*,*) "negative dtau ",dtau
             endif
             if (i0(irbb) > 1.e20) then
                write(*,*) "intensity bug ",i0(irbb), irbb
             endif
             if (snu > 1.d20) then
                write(*,*) "snu bug ",snu, jnu, alphanu
             endif

             i0(iRBB) = i0(iRBB) +  exp(-tau(irBB)) * (1.d0-exp(-dtau))*snu

             tau(iRBB) = tau(iRBB) + dtau

          enddo
       enddo
       currentPosition = currentPosition + (distArray(ntau)+1.d-3*grid%halfSmallestSubcell) * direction
       totdist = totdist + (distArray(ntau)+1.d-3*grid%halfSmallestSubcell)

       if (nBug > 10000) then
          write(*,*) "bug",currentPosition,nTau,distArray(nTau),modulus(currentPosition)
       endif

       if (firstSubcell) then
          phi = phiAv / phiNorm
          firstSubcell = .false.
       endif

    enddo

    if (onTheSpot) then
       iCont = 0.d0
       i0 = 0.d0
    endif

    if (hitPhotosphere) then ! don't include weight below - that's done when jbar is calculated
       iElement = getElement(source(sourcenumber)%surface, photoDirection)
       do iFreq = 1, nFreq
          !          write(*,*) freq(ifreq),iCont(ifreq),i_nu(source(sourceNumber), &
          !               freq(iFreq), iElement)*cosTheta*exp(-tauCont(iFreq)),taucont(ifreq)
          !          write(*,*) freq(ifreq),i_nu(source(sourcenumber),freq(ifreq),ielement)
          !          iCont(iFreq) = iCont(iFreq) + i_nu(source(sourceNumber), freq(iFreq), iElement)*cosTheta*exp(-tauCont(iFreq))
          if (.not.onTheSpot) then
             iCont(iFreq) = iCont(iFreq) + i_nu(source(sourceNumber), freq(iFreq), iElement)*exp(-tauCont(iFreq))
          else
             iCont(iFreq) = iCont(iFreq) + i_nu(source(sourceNumber), freq(iFreq), iElement)
          endif
       enddo
    endif


    if (hitPhotosphere) then ! don't include weight below - that's done when jbar is calculated
       iElement = getElement(source(sourcenumber)%surface, photoDirection)
       do iRBB = 1, nRBBTrans
          iAtom = indexAtom(iRBB)
          iTrans = indexRBBTrans(iRBB)
          if (.not.onTheSpot) then
             i0(iRBB) = i0(iRBB) + i_nu(source(sourceNumber), thisAtom(iATom)%transFreq(iTrans), iElement)*exp(-tau(iRBB))
          else
             i0(iRBB) = i0(iRBB) + i_nu(source(sourceNumber), thisAtom(iATom)%transFreq(iTrans), iElement)
          endif
       enddo
    endif

    deallocate(tau)
  end subroutine getRay

  function phiProf(dv, b) result (phi)
    real(double) :: dv, b
    real(double) :: fac, phi

    phi = 1.d0 / (b * sqrtPi)
    fac = (dv/b)**2
    phi = phi * exp(-fac)

  end function phiProf

  subroutine calculateJbar(grid, thisOctal, subcell, thisAtom, nRay, position, direction, rayDeltaV, &
       ds,  i0, iTrans, jbar, nPops, &
       weight,  tauAv)
    type(GRIDTYPE) :: grid
    type(VECTOR) :: position(:), direction(:), thisPosition, startVel, thisVel
    real(double) :: dds, dv, rayDeltaV(:)
    integer, parameter :: ns = 2
    real(double) :: inu
    integer :: i
    real(double) :: weight(:)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    type(MODELATOM) :: thisAtom
    integer :: nRay
    Real(double) :: ds(:), i0(:), nPops(:)
    integer :: iTrans
    real(double) :: jbar
    integer :: iRay
    real(double) :: nLower, nUpper
    real(double) :: jBarInternal, jBarExternal
    real(double) :: alphanu, jnu, etaline
    integer :: iUpper, iLower
    real(double) :: tau, snu, sumWeight,inuAv
    real(double) :: a, bul, blu
    logical,save :: first = .true.
    real(double) :: dtau, tauav

    !$OMP THREADPRIVATE (first)

    jBarExternal = 0.d0
    jBarInternal = 0.d0
    a = 0.d0; bul = 0.d0; blu = 0.d0

    if (thisAtom%transType(iTrans) == "RBB") then

       iUpper = thisAtom%iUpper(iTrans)
       iLower = thisAtom%iLower(iTrans)


       sumWeight = 0.d0
       tauAv = 0.d0
       inuAv = 0.d0
       do iRay = 1, nRay
          nLower = nPops(iLower)
          nUpper = nPops(iUpper)

          call returnEinsteinCoeffs(thisAtom, iTrans, a, Bul, Blu)

          etaLine = hCgs * a * thisAtom%transFreq(iTrans)
          etaLine = etaLine *  nPops(iUpper)

          thisPosition = position(iray)
          startVel = amrGridVelocity(grid%octreeRoot, thisPosition, startOctal = thisOctal, actualSubcell = subcell) 

          dds = (ds(iray)/1.d10)/dble(ns)
          inu = 0.d0
          tau = 0.d0
          do i = 1, ns
             thisPosition = thisPosition + dds * direction(iray)
             thisVel = amrGridVelocity(grid%octreeRoot, thisPosition, startOctal = thisOctal, actualSubcell = subcell) 
             thisVel= thisVel - startVel

             dv = rayDeltaV(iray) - (thisVel .dot. direction(iray))

             alphanu = (hCgs*thisAtom%transFreq(iTrans)/fourPi)
             alphanu = alphanu * (nLower * Blu - nUpper * Bul) * &
                  phiProf(dv, thisOctal%microturb(subcell)) /thisAtom%transFreq(iTrans)

             if (alphanu < 0.d0) then
                alphanu = 0.d0
                if (first) then
                   write(*,*) "negative opacity warning in calcjbar",iUpper,iLower,nLower,nUpper,thisAtom%name
                   first = .false.
                endif
             endif

             dtau = alphaNu * dds * 1.d10
 
             jnu = (etaLine/fourPi) * phiProf(dv, thisOctal%microturb(subcell)) /thisAtom%transFreq(iTrans)

             if (alphanu /= 0.d0) then
                snu = jnu/alphanu
             else
                snu = tiny(snu)
             endif


             inu = inu + exp(-tau) * snu * (1.d0 - exp(-dtau))


             tau = tau + dtau
             
          enddo

!          if ((thisAtom%name == "HeII").and.(ilower==1).and.(iupper==2)) &
!               write(*,*) "delta v ",raydeltaV(iray)*cspeed/1e5,tau,inu

          tauAv = tauAv + tau * weight(iray)
          jBarInternal = jbarInternal + inu  * weight(iray)
          jBarExternal = jbarExternal + i0(iray) * exp(-tau) * weight(iray)

          inuAv = inuAv + inu * weight(iray)

!          if ((thisAtom%name == "HeII").and.(ilower==1).and.(iupper==2).and.(tau < 1.d0)) &
!               write(*,*) "jext ",jbarExternal,i0(iray),tau,weight(iray),dv*cspeed/1.e5
          sumWeight = sumWeight + weight(iRay) 

       enddo
       tauAv = tauAv / sumWeight
       inuAv = inuAv / sumWeight

 
       jBarExternal = jBarExternal / sumWeight
       jBarInternal = jBarInternal / sumWeight

       

          
       jbar = (jBarExternal + jBarInternal)


       if ((thisAtom%name == "HeII").and.(ilower==1))&
          write(*,*) ilower," line data ",jbarExternal,jbarInternal,jbar,tauav

    endif
          
  end subroutine calculateJbar

  subroutine calculateJbarCont(thisOctal, subcell, nAtom, thisAtom, ne, nray, ds, freq, nfreq, &
       iCont, jBarCont, weight)
    use input_variables, only : opticallyThickContinuum, onTheSpot
    real(double) :: iCont(:,:), jBarCont(:), ds(:), ne
    integer :: nAtom
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    type(MODELATOM) :: thisAtom(:)
    integer :: nfreq
    real(double) :: freq(:), weight(:)
    integer :: iray, nRay
    real :: tau
    real(double), allocatable :: jBarContExternal(:), jBarContInternal(:)
    real(double) :: alphanu
    integer :: iFreq
    real(double) :: jnu, snu
    real(double) :: nstar(10,50), expMinusTau
    integer :: iAtom, j
    allocate(jBarContExternal(1:nFreq))
    allocate(jBarContInternal(1:nFreq))

    jBarCont = 0.d0
    jBarContInternal = 0.d0
    jBarContExternal = 0.d0

    do iAtom = 1, nAtom
       do j = 1, thisAtom(iAtom)%nLevels - 1
          nStar(iAtom,j) = BoltzSahaGeneral(thisAtom(iAtom), j, ne, &
               dble(thisOctal%temperature(subcell))) * &
               Thisoctal%atomlevel(subcell, iAtom,thisAtom(iatom)%nLevels)
!          write(*,*) "iatom ", iatom, " j ",nstar(iatom,j), thisOctal%atomlevel(subcell, iatom,thisAtom(iatom)%nLevels)

       enddo
       if (any(nstar(iatom,1:thisAtom(iatom)%nlevels-1) < 0.d0)) then
          write(*,*) "atom: ",thisAtom(iatom)%name
          write(*,*) "nstar bug ",nstar(iAtom,1:thisAtom(iatom)%nlevels-1)
          stop
       endif
    enddo


    do iRay = 1, nRay
       do iFreq = 1, nFreq

          tau = 0.d0
          snu = 0.d0
          if (opticallyThickContinuum) then
             alphanu = bfOpacity(freq(ifreq), nAtom, thisAtom, thisOctal%atomLevel(subcell,:,:), ne, &
                  dble(thisOctal%temperature(subcell)),ifreq=ifreq)
             jnu =  bfEmissivity(freq(ifreq), nAtom, thisAtom, &
                  thisOctal%atomLevel(subcell,:,:), nstar, &
                  dble(thisOctal%temperature(subcell)), thisOctal%ne(subcell), &
                  ifreq=ifreq)
             tau = alphaNu * ds(iray)
             if (alphanu /= 0.d0) then
                snu = jnu/alphanu
             else
                snu = tiny(snu)
             endif
         endif
         if (tau < 0.d0) then
            write(*,*) "tau bug ", tau
            stop
         endif
         if (tau < 60.d0) then
            expMinusTau = exp(-tau)
         else
            expMinusTau = 0.d0
         endif

         if (.not.onTheSpot) then
            jBarContExternal(iFreq) = jBarContExternal(iFreq) + iCont(iray, iFreq) * expMinusTau * weight(iRay)
            jBarContInternal(iFreq) = jBarContInternal(iFreq) + snu * (1.d0 - expMinusTau) * weight(iRay)
         else
            jBarContExternal(iFreq) = jBarContExternal(iFreq) + iCont(iray, iFreq)  * weight(iRay)
            jBarContInternal(iFreq) = jBarContInternal(iFreq) + 0.d0
         endif

       enddo
       if (any(jbarcontexternal(1:nfreq) < 0.d0)) then
          write(*,*) "negative jbarcontexternal"
          write(*,*) icont(iray,1:nfreq)
          stop
       endif
    enddo
    jBarContExternal = jBarContExternal / SUM(weight(1:nRay))
    jBarContInternal = jBarContInternal / SUM(weight(1:nRay))
    if (onthespot) jbarcontinternal = 0.d0
    jbarCont = jBarContExternal + jBarContInternal

    if (any(jbarCont < 0.d0)) then
       write(*,*) "fatal jbar cont bug"
       write(*,*) "ext ",jbarContExternal(1:nfreq)
       write(*,*) "int ",jbarContInternal(1:nfreq)
       write(*,*) "weight ",weight(1:nray)
       write(*,*) "icont ",icont(1:nray,1:nfreq)
       stop
    endif

  end subroutine calculateJbarCont


  subroutine atomLoop(grid, nAtom, thisAtom, nSource, source)

    use input_variables, only : debug, rcore, lte, vturb
    use messages_mod, only : myRankIsZero
    use gridio_mod, only : writeAmrGrid
    use utils_mod, only : ngstep
    use random_mod
    use amr_mod, only: getOctalArray, sortOctalArray
#ifdef MEMCHECK
    use memory_mod, only : resetGlobalMemory
#endif
#ifdef MPI
    use amr_mod, only : countVoxels
    include 'mpif.h'
#endif

    type(SOURCETYPE) :: source(:)
    integer :: nSource
    integer :: nAtom
    type(GRIDTYPE) :: grid
    type(MODELATOM) :: thisAtom(:)
    type(VECTOR), allocatable :: position(:), direction(:)
    real(double), allocatable :: rayDeltaV(:)
    integer :: nOctal, iOctal, subcell
    real(double), allocatable :: ds(:), phi(:), i0(:,:)
    real(double), allocatable :: Hcol(:), HeICol(:), HeIICol(:)
    integer :: nRay
    type(OCTAL), pointer :: thisOctal
    integer, parameter :: maxIter = 1000, maxRay = 200000
    logical :: popsConverged, gridConverged 
    integer :: iRay, iTrans, iter,i 
    integer :: iStage
    real(double), allocatable :: oldpops(:,:), newPops(:,:), dPops(:,:), mainoldpops(:,:)
    real(double) :: newNe
    real(double), parameter :: underCorrect = 1.d0
    real(double), parameter :: underCorrectne = 1.d0
    real(double) :: dne

    real(double) :: fac
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: ioctal_beg, ioctal_end
    real(double) :: maxFracChange
    logical :: fixedRays
    integer(bigint) :: iseed
    integer, allocatable :: sourceNumber(:)
!    type(VECTOR) :: posVec
!    real(double) :: r
    real(double), allocatable :: cosTheta(:)
    real(double), allocatable :: weightFreq(:), weightOmega(:)
    real(double), allocatable :: iCont(:,:)
    logical, allocatable :: hitPhotosphere(:)
    integer, parameter :: maxFreq = 2000
    real(double) :: freq(maxFreq), dFreq(maxFreq)
    real(double), allocatable :: jnuCont(:)
    integer :: nFreq, nhit, iRBB
    integer :: nRBBTrans
    integer :: indexRBBTrans(1000), indexAtom(1000)
    real(double) :: ne, tauav
    integer :: iAtom
    integer :: nHAtom, nHeIAtom, nHeIIatom !, ir, ifreq
    real(double) :: nstar, ratio, ntot
    real(double), parameter :: convergeTol = 1.d-4, gridtolerance = 1.d-2
    integer :: neIter, itmp
    logical :: recalcJbar,  firstCheckonTau
    character(len=80) :: message, ifilename
    real :: r
    logical :: ionized
    integer :: nIter, idump, nt, nInuse, nConverged
    real(double) :: percentageConverged
    real(double), save, allocatable :: oldpops1(:,:), oldpops2(:,:), oldpops3(:,:), oldpops4(:,:)
    integer, parameter :: iNgStep = 5

#ifdef MPI
    ! For MPI implementations
    integer       ::   my_rank        ! my processor rank
    integer       ::   np             ! The number of processes
    integer       ::   n_rmdr, m
    integer       ::   ierr           ! error flag
    integer       ::   nVoxels
    real(double), allocatable :: tArrayd(:),tempArrayd(:)
#endif

#ifdef _OPENMP
    integer :: omp_get_thread_num
    integer(bigInt) ::  iCellSeed
#endif
!$OMP THREADPRIVATE (oldpops1, oldpops2, oldpops3, oldpops4)



!    blockHandout = .false. 

#ifdef MPI
    ! FOR MPI IMPLEMENTATION=======================================================
    !  Get my process rank # 
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  
    ! Find the total # of precessor being used in this run
    call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
#endif


    freq = 0.d0
    indexRBBTrans = 0
    nfreq = 0; nRBBTrans = 0; tauAv = 0.d0
    call createRBBarrays(nAtom, thisAtom, nRBBtrans, indexAtom, indexRBBTrans)

     call createContFreqArray(nFreq, freq, nAtom, thisAtom, nsource, source)

     dfreq(1) = freq(2)-freq(1)
     do i = 2, nFreq
        dfreq(i) = freq(i) - freq(i-1)
     enddo


!     open(69,file="photoionization.dat",status="unknown",form="formatted")
!     do i = 1, nfreq
!        do itrans = 1, thisAtom(1)%nTrans
!          if (thisAtom(1)%name == "HI".and.thisAtom(1)%transType(iTrans) == "RBF") then ! photoionization
!             if (thisAtom(1)%iLower(iTrans) == 1) then
!                xH = photoCrossSection(thisAtom(1), iTrans, 1, freq(i))
!             endif
!          endif
!       enddo
!        do itrans = 1, thisAtom(2)%nTrans
!          if (thisAtom(2)%name == "HeI".and.thisAtom(2)%transType(iTrans) == "RBF") then ! photoionization
!             if (thisAtom(2)%iLower(iTrans) == 2) then
!                xHeI = photoCrossSection(thisAtom(2), iTrans, 2, freq(i))
!             endif
!          endif
!       enddo
!        do itrans = 1, thisAtom(3)%nTrans
!          if (thisAtom(3)%name == "HeII".and.thisAtom(3)%transType(iTrans) == "RBF") then ! photoionization
!             if (thisAtom(3)%iLower(iTrans) == 2) then
!                xHeII = photoCrossSection(thisAtom(3), iTrans, 2, freq(i))
!             endif
!          endif
!       enddo
!
!       write(69,*) freq(i),  xHeI, xHeII
!    end do
!    close(69)


     ionized = .false.
     if (grid%geometry == "gammavel") ionized = .true.
     if (grid%geometry == "wrshell") ionized = .true.
     if (grid%geometry == "wind") ionized = .true.

     ionized = .true.
     call writeInfo("Allocating levels and applying LTE...")
     call allocateLevels(grid, grid%octreeRoot, nAtom, thisAtom, nRBBTrans, nFreq, ionized)
#ifdef MEMCHECK
     call resetGlobalMemory(grid)
#endif
     call setMicroturb(grid%octreeRoot, dble(vTurb))

     call writeInfo("Done.")

    do i = 1, nAtom
       call addCrossSectionstoAtom(thisAtom(i), nFreq, freq)
    enddo



    if (LTE) goto 666

    nHAtom = 0
    nHeIAtom = 0
    nHeIIAtom = 0
    do iAtom = 1, nAtom
       if (thisAtom(iAtom)%name == "HI") nHatom = iAtom
       if (thisAtom(iAtom)%name == "HeI") nHeIatom = iAtom
       if (thisAtom(iAtom)%name == "HeII") nHeIIatom = iAtom
    enddo

! ly-alpha in detailed balance

!    do iAtom = 1, nAtom
!       do iTrans = 1, thisAtom(iAtom)%nTrans
!          if (thisAtom(iAtom)%transType(iTrans) == "RBB") then 
!             if ((thisAtom(iAtom)%iLower(iTrans) == 1).and.(thisAtom(iatom)%iUpper(itrans)==2)) then
!                thisAtom(iAtom)%inDetailedBalance(iTrans) = .true.
!             endif
!          endif
!       enddo
!    enddo


! Lyman continuum in detailed balance

    do iAtom = 1, nAtom
       do iTrans = 1, thisAtom(iAtom)%nTrans
          if (thisAtom(iAtom)%name == "HI".and.thisAtom(iAtom)%transType(iTrans) == "RBF") then ! photoionization
             if (thisAtom(iAtom)%iLower(iTrans) == 1) thisAtom(iAtom)%inDetailedBalance(iTrans) = .true.
          endif
       enddo
    enddo
!!!

! Lyman continuum in detailed balance

!    do iAtom = 1, nAtom
!       do iTrans = 1, thisAtom(iAtom)%nTrans
!          if (thisAtom(iAtom)%name == "HeI".and.thisAtom(iAtom)%transType(iTrans) == "RBF") then ! photoionization
!             if (thisAtom(iAtom)%iLower(iTrans) == 1) thisAtom(iAtom)%inDetailedBalance(iTrans) = .true.
!          endif
!       enddo
!    enddo
!
!! Lyman continuum in detailed balance

!    do iAtom = 1, nAtom
!       do iTrans = 1, thisAtom(iAtom)%nTrans
!          if (thisAtom(iAtom)%name == "HeII".and.thisAtom(iAtom)%transType(iTrans) == "RBF") then ! photoionization
!             if (thisAtom(iAtom)%iLower(iTrans) == 1) thisAtom(iAtom)%inDetailedBalance(iTrans) = .true.
!             write(*,*) "lyman cont detailed balance ",itrans
!          endif
!       enddo
!    enddo

!! Lyman lines in detailed balance
!
!    do iAtom = 1, nAtom
!       do iTrans = 1, thisAtom(iAtom)%nTrans
!          if (thisAtom(iAtom)%name == "HeII".and.thisAtom(iAtom)%transType(iTrans) == "RBB") then ! photoionization
!             if (thisAtom(iAtom)%iLower(iTrans) == 1) thisAtom(iAtom)%inDetailedBalance(iTrans) = .true.
!             write(*,*) "lyman cont detailed balance ",itrans
!          endif
!       enddo
!    enddo

    allocate(octalArray(grid%nOctals))
    nOctal = 0
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)
    call  sortOctalArray(octalArray,grid)

    nInUse = 0 
    do iOctal = 1, SIZE(octalArray)
       thisOctal => octalArray(iOctal)%content
       do subcell = 1, thisOctal%maxChildren
          if (.not.thisOctal%hasChild(subcell)) then
             
             thisOctal%inFlow(subcell) = &
                  thisOctal%inFlow(subcell).and.(.not.insideSource(thisOctal, subcell, nsource, Source))
             
             if (thisOctal%inflow(subcell).and.(thisOctal%temperature(subcell) > 3000.)) nInuse = nInUse + 1
          endif
       enddo
    enddo
    if (writeoutput) write(*,*) "Total number of cells in use ",nInuse


    call randomNumberGenerator(randomSeed=.true.)
!    if (myrankiszero) then
!       call testRays(grid, nSource, source)
!    endif
    call torus_mpi_barrier

    call randomNumberGenerator(getiseed = iseed)
    call randomNumberGenerator(syncIseed = .true.)

     if (myRankisZero) then
        open(69, file="cmf_convergence.dat", status="unknown", form="formatted")
        write(69,'(a)') &
!             01234567890123456789012345678901234567890123456789
             "       Iter    Tol  %Converged   nRays     Fix ray"
        close(69)
     endif

     nRay = 100
    do iStage = 1, 2


       if (iStage == 1) then
          fixedRays = .true.
       else
          fixedRays = .false.
       endif

       gridConverged = .false.
       nIter = 0
       idump = 0
       itmp = 0
       do while (.not.gridConverged)

          nIter = nIter + 1
          idump = idump + 1
          write(ifilename,'(a,i2.2,a)') "ionization",idump,".dat"
          if (myRankisZero) then
             open(69, file=ifilename, status="unknown", form="formatted")
             write(69,'(a)') &
                  "# log(r/R) , log N(HI), log N(HeI), log N(HeII), log N(HeIII)"
             close(69)
          endif

          if (doTuning) call tune(6, "One cmf iteration")  ! start a stopwatch

          
          ! default loop indecies
          ioctal_beg = 1
          ioctal_end = SIZE(octalArray)       

          

#ifdef MPI
    
            ! Set the range of index for octal loop used later.     
            np = nThreadsGlobal
            n_rmdr = MOD(SIZE(octalArray),np)
            m = SIZE(octalArray)/np
            
            if (myRankGlobal .lt. n_rmdr ) then
               ioctal_beg = (m+1)*myRankGlobal + 1
               ioctal_end = ioctal_beg + m
            else
               ioctal_beg = m*myRankGlobal + 1 + n_rmdr
               ioctal_end = ioctal_beg + m - 1
            end if
            
#endif
            if (fixedRays) then
               call randomNumberGenerator(synciSeed = .true.)
            else
               call randomNumberGenerator(randomSeed=.true.)
            endif


            !$OMP PARALLEL DEFAULT (NONE) &
            !$OMP PRIVATE (iOctal, thisOctal, subcell, i0, position, direction, nt, iCellSeed) &
	    !$OMP PRIVATE(rayDeltaV, ds, phi, hcol, heicol, heiicol, hitphotosphere, sourcenumber, costheta,weightfreq) &
	    !$OMP PRIVATE(weightOmega, icont, neiter, iter,popsConverged, oldpops, mainoldpops, firstCheckonTau) &
	    !$OMP PRIVATE(fac,dne,message,ifilename,itmp,ne,recalcjbar,ratio,nstar,dpops,newne) &
	    !$OMP PRIVATE(nhit, jnucont,tauav,newpops,ntot,r,iatom,itrans) &
            !$OMP SHARED(octalArray, grid, ioctal_beg, ioctal_end, nsource, nray, nrbbtrans, indexRbbtrans, indexatom) &
	    !$OMP SHARED(freq,dfreq,nfreq, natom,myrankiszero,debug,rcore, iseed, fixedRays, source, thisAtom, myrankGlobal, writeoutput)

            !$OMP MASTER
            call randomNumberGenerator(getIseed=iseed)
            !$OMP END MASTER
            allocate(oldPops(1:nAtom,1:maxval(thisAtom(1:nAtom)%nLevels)))
            allocate(mainoldPops(1:nAtom,1:maxval(thisAtom(1:nAtom)%nLevels)))
            allocate(newPops(1:nAtom,1:maxval(thisAtom(1:nAtom)%nLevels)))
            allocate(dPops(1:nAtom,1:maxval(thisAtom(1:nAtom)%nLevels)))
            allocate(jnuCont(1:nFreq)) 
            allocate(ds(1:nRay))
            allocate(phi(1:nRay))
            allocate(rayDeltaV(1:nRay))
            allocate(i0(1:nRBBTrans, 1:nRay))
            allocate(Hcol(1:nRay))
            allocate(HeIcol(1:nRay))
            allocate(HeIIcol(1:nRay))
            allocate(sourceNumber(1:nRay))
            allocate(cosTheta(1:nRay))
            allocate(weightFreq(1:nRay))
            allocate(weightOmega(1:nRay))
            allocate(hitPhotosphere(1:nRay))
            allocate(iCont(1:nRay,1:nFreq))
            allocate(position(1:nray))
            allocate(direction(1:nray))

            allocate(oldpops1(nAtom, maxval(thisAtom(1:nAtom)%nLevels)))
            allocate(oldpops2(nAtom, maxval(thisAtom(1:nAtom)%nLevels)))
            allocate(oldpops3(nAtom, maxval(thisAtom(1:nAtom)%nLevels)))
            allocate(oldpops4(nAtom, maxval(thisAtom(1:nAtom)%nLevels)))
            
            if (fixedRays) then
               call randomNumberGenerator(putIseed = iseed)
               call randomNumberGenerator(reset = .true.)
            else
               call randomNumberGenerator(randomSeed=.true.)
            endif
            call randomNumberGenerator(getReal=r)
            write(*,*) myrankGlobal, r




            !$OMP DO SCHEDULE(DYNAMIC,1)
          do iOctal = ioctal_beg, ioctal_end

             thisOctal => octalArray(iOctal)%content

             nt = 0
#ifdef _OPENMP
             nt = omp_get_thread_num()
             iCellSeed = iseed + thisOctal%label(1)*10000
            if (fixedRays) then
               call randomNumberGenerator(putIseed = iCellSeed)
               call randomNumberGenerator(reset = .true.)
            endif
#endif

!          if (doTuning) call tune(6, "One octal iteration")  ! start a stopwatch
             
!             do subcell = thisOctal%maxChildren,1,-1
             do subcell = 1, thisOctal%maxChildren

!          open(31, file="r.dat",status="unknown",form="formatted")
!          do ir = 1, 50
!             r = log10(20.d0*rsol/1.d10) + (log10(2000.d0*rSol/1.d10)-log10(20.d0*rsol/1.d10))*real(ir-1)/49.d0
!             r = 10.d0**r
!             posVec = VECTOR(r, 0.d0, 0.d0)
!             thisOctal => grid%octreeRoot
!             call findSubcellLocal(posVec, thisOctal, subcell)

                if (.not.thisOctal%hasChild(subcell)) then

                   thisOctal%inFlow(subcell) = &
                        thisOctal%inFlow(subcell).and.(.not.insideSource(thisOctal, subcell, nsource, Source))

                   if (thisOctal%inflow(subcell).and.(thisOctal%temperature(subcell) > 3000.)) then 

                      write(*,*) myrankGlobal, nt, subcell, " doing octal ",iOctal, " of ",ioctal_end

                        !.and.(.not.thisOctal%fixedTemperature(subcell))) then
                      nHit = 0
                      do iRay = 1, nRay
                         call getRay(grid, thisOCtal, subcell, position(iray), direction(iray), rayDeltaV(iray),&
                              ds(iRay), phi(iRay), i0(1:nRBBTrans,iRay), Hcol(iray), HeICol(iRay), HeIICol(iRay),&
                              nAtom, thisAtom, source, nSource, hitPhotosphere(iRay), sourceNumber(iray), &
                              cosTheta(iRay), weightFreq(iRay), weightOmega(iray), &
                              nRBBTrans, indexRBBTrans, indexAtom,  &
                              nfreq, freq, iCont(iray,1:nFreq))
                         if (hitPhotosphere(iray)) nHit = nHit + 1
                      enddo
                      iter = 0
                      neiter = 0
                      popsConverged = .false.
                      thisOctal%newatomLevel(subcell,:, :) = thisOctal%atomLevel(subcell,:, :)
                      ne = thisOctal%ne(subcell)

                      recalcJbar = .true.
                      firstCheckOnTau = .true.

                      do while (.not.popsConverged)
                         iter = iter + 1
                         mainoldpops = thisOctal%newatomLevel(subcell,1:nAtom,1:)

                         oldpops = thisOctal%newatomLevel(subcell,1:nAtom,1:)

                         if (recalcJbar) then
                            call calculateJbarCont(thisOctal, subcell, nAtom, thisAtom, ne, nray, ds, &
                                 freq, nfreq, &
                                 iCont, jnuCont, weightOmega)

                            thisOctal%JnuCont(subcell, :) = jnuCont(:)

                            do iRBB = 1, nRBBTrans
                               iTrans = indexRBBTrans(iRBB)
                               iAtom =  indexAtom(iRBB)
                               call calculateJbar(grid, thisOctal, subcell, thisAtom(iAtom), nRay, position(1:nray), &
                                    direction(1:nray), rayDeltaV(1:nray), ds(1:nRay), &
                                    i0(iRBB,1:nRay), iTrans, thisOctal%jnuLine(subcell,iRBB), &
                                    thisOctal%newAtomLevel(subcell,iAtom,1:), &
                                    weightFreq(1:nRay) * weightOmega(1:nRay), tauAv)
!                               if  ( (thisAtom(iatom)%name == "HeII").and.(thisAtom(iatom)%ilower(itrans)==1)&
!                                    .and.(thisAtom(iatom)%iupper(itrans)==2) ) then
!                                  write(*,*) "after calc jbar ",thisOctal%jnuline(subcell,irBB),irbb
!                               endif
                            enddo
                         endif



                            neIter = neIter + 1
                            oldpops = thisOctal%newatomLevel(subcell,1:nAtom,1:)

                            newpops = thisOctal%newatomLevel(subcell, 1:nAtom, 1:)
                            call solveLevels(thisOctal, subcell, newPops, &
                                 thisOctal%jnuLine(subcell,1:nRBBTrans),  &
                                 dble(thisOctal%temperature(subcell)), nAtom, thisAtom, ne, &
                                 thisOctal%rho(subcell),&
                                 jnuCont, freq, dfreq, nfreq)
                            
                            dPops = newPops - thisOctal%newAtomLevel(subcell,1:nAtom,:) 
                            

                            thisOctal%newAtomLevel(subcell,1:nAtom,:)  = &
                                 thisOctal%newAtomLevel(subcell,1:nAtom,:)  + underCorrect * dPops


!                            ! improve convergence...
!                            thisOctal%atomLevel(subcell,:,:) = thisOctal%newatomLevel(subcell,:,:)!!!!!!!!!!!!!!!!!!!!!!!!!!


                            where (abs(thisOctal%newAtomLevel(subcell,1:nAtom,:)) < 1.d-30)
                               thisOctal%newAtomLevel(subcell,1:nAtom,:) = 1.d-30
                            end where

                            oldpops1 = oldpops2
                            oldpops2 = oldpops3
                            oldpops3 = oldpops4
                            oldpops4(1:nAtom,:) = thisOctal%newAtomLevel(subcell,1:nAtom,:)

                            if (mod(iter, iNgStep) == 0) then
!                               if (writeoutput) write(*,*) "Doing Ng acceleration step"
                               do iAtom = 1, nAtom
!                                  if (writeoutput) then
!                                     write(*,*) "oldpops 1",real(oldpops1(iatom, 1:6))
!                                     write(*,*) "oldpops 2",real(oldpops2(iatom, 1:6))
!                                     write(*,*) "oldpops 3",real(oldpops3(iatom, 1:6))
!                                     write(*,*) "oldpops 4",real(oldpops4(iatom, 1:6))
!                                  endif
                                  call ngStep(thisOctal%newAtomLevel(subcell, iAtom, :), &
                                       oldpops1(iAtom, :), oldpops2(iAtom, :), &
                                       oldpops3(iAtom, :), oldpops4(iAtom, :), length=thisAtom(iAtom)%nLevels)
!                               if (writeoutput) write(*,*) "newpops ",real(thisOctal%newAtomLevel(subcell, iAtom, 1:6))
                               enddo
                            endif
                            newNe = tiny(newNe)
                            do iAtom = 1, nAtom
                               if (iAtom /= nAtom) then
                                  if (thisAtom(iAtom)%nz /= thisAtom(iAtom+1)%nz) then
                                     newne =  newNe +  &
                                          SUM(thisOctal%newAtomLevel(subcell, iAtom, 1:thisAtom(iAtom)%nLevels-1)) &
                                          * dble(thisAtom(iAtom)%charge)
                                     newne = newne + thisOctal%newAtomLevel(subcell, iAtom, thisAtom(iAtom)%nLevels) &
                                          * dble((thisAtom(iAtom)%charge+1))
                                  else
                                     newne =  newNe +  &
                                          SUM(thisOctal%newAtomLevel(subcell, iAtom, 1:thisAtom(iAtom)%nLevels-1)) &
                                          * dble(thisAtom(iAtom)%charge)
                                  endif
                               else
                                  newne =  newNe +  &
                                       SUM(thisOctal%newAtomLevel(subcell, iAtom, 1:thisAtom(iAtom)%nLevels-1)) &
                                       * dble(thisAtom(iAtom)%charge)
                                  newne = newne + thisOctal%newAtomLevel(subcell, iAtom, thisAtom(iAtom)%nLevels) &
                                       * dble((thisAtom(iAtom)%charge+1))
                               endif
                            enddo
                            if (newNe < 0.d0) then
                               write(*,*) "newNe is negative"
                               newNE = TINY(newNE)
                            endif

                           ntot = 0.d0
                           ntot = ntot + SUM(thisOctal%newAtomLevel(subcell,1, &
                                1:thisAtom(1)%nLevels))
                           if (natom > 1) then
                              ntot = ntot + SUM(thisOctal%newAtomLevel(subcell,2, &
                                   1:thisAtom(2)%nLevels-1))
                           endif
                           if (natom  > 2) then
                              ntot = ntot + SUM(thisOctal%newAtomLevel(subcell,3, &
                                   1:thisAtom(3)%nLevels))
                           endif

                            if (debug) then
                               write(*,*) "Iteration: ",iter
                               do iatom = 1, size(thisAtom)
                                  write(*,*) "Atom: ",thisAtom(iatom)%name
                                  write(*,*) "Number density: ",thisOctal%rho(subcell)*thisOctal%atomAbundance(subcell, iatom)
                                  do i = 1, thisAtom(iatom)%nLevels-1


                                     nStar = boltzSahaGeneral(thisAtom(iAtom), i, ne, &
                                          dble(thisOctal%temperature(subcell))) * &
                                          thisOctal%newAtomLevel(subcell, iatom, thisAtom(iAtom)%nLevels)
                                     ratio = boltzSahaGeneral(thisAtom(iAtom), i, ne, &
                                          dble(thisOctal%temperature(subcell)))
                        
                                 write(*,'(i2,1x,1p,7e12.3)') i,thisOctal%newAtomLevel(subcell,iatom,i),mainoldpops(iatom,i), &
                                      abs((thisOctal%newAtomLevel(subcell,iAtom,i)-max(mainoldpops(iAtom,i),1.d-20)) &
                                      / max(mainoldpops(iAtom,i),1.d-20)), thisOctal%newAtomLevel(subcell,iatom,i)/nstar, &
                                      thisOctal%newAtomLevel(subcell,iAtom,i)-mainoldpops(iAtom,i),nstar,ratio
                                 
                              enddo
                           enddo
                           !                           if (thisatom(1)%indetailedbalance(4)) then
!                              write(*,*) "Halpha is in detailed balance",tauHalpha
                           !                           else
                           !                              write(*,*) "Halpha is NOT in detailed balance",tauHalpha
                           !                           endif
                           write(*,*) "Ne: ",newne
                           write(*,*) "T: ",thisOctal%temperature(subcell)
                           write(*,*) "R: ",modulus(subcellCentre(thisOctal,subcell))/rCore
                           write(*,*) "log(H I/ Ntot): ", &
                                log10(SUM(thisOctal%newAtomLevel(subcell,1,1:thisAtom(1)%nlevels-1)) /ntot)
                           if (nAtom > 1) &
                           write(*,*) "log(He I/ Ntot): ", &
                                log10(SUM(thisOctal%newAtomLevel(subcell,2,1:thisAtom(2)%nlevels-1)) / ntot)
                           if (nAtom > 2) &
                           write(*,*) "log(He II/ Ntot): ", &
                                log10(SUM(thisOctal%newAtomLevel(subcell,3,1:thisAtom(3)%nlevels-1)) / ntot)
                           if (nAtom > 2) &
                           write(*,*) "log(He III/ Ntot): ", &
                                log10(thisOctal%newAtomLevel(subcell,3,thisAtom(3)%nlevels) / ntot)
                        endif

                         fac = -1.d30
                         where (mainoldPops == 0.d0)
                            mainoldPops = 1.d-20
                         end where


                         fac = -1.d30
                         do iAtom = 1, nAtom
                            fac = max(fac,maxval(abs((thisOctal%newAtomLevel(subcell,iatom,1:thisAtom(iAtom)%nLevels-1) &
                                 - mainoldpops(iAtom,1:thisAtom(iAtom)%nLevels-1))/mainoldpops(iAtom,1:thisAtom(iAtom)%nLevels-1))))
                         enddo
                         dne = newne - ne
                         ne = ne + undercorrectne * dne
                         if ((fac < convergeTol).and.(iter > 1)) popsConverged = .true.
                         if (iter == maxIter) then
                            popsConverged = .true.
                            write(message,'(a,e12.5)') "Maximum number of iterations reached in pop solver. fac = ",fac
                            write(*,*) "fac, converged tol", fac,convergetol
                            write(*,*) trim(message)
                            call writeWarning(message)
                         endif
!                         write(*,*) "iter",iter,fac
                      enddo
!                      write(*,*) "Main iteration route converged after ", iter, " iterations"
                      thisOctal%ne(subcell) = ne
                      
!                        if (myRankisZero) then
!                           open(69, file=ifilename, status="old", position = "append", form="formatted")
!                           if (nAtom == 1) write(69,'(6f10.4)') log10(modulus(subcellCentre(thisOctal,subcell))/rCore), &
!                                log10(SUM(thisOctal%newAtomLevel(subcell,1,1:thisAtom(1)%nlevels-1)) /ntot)
!
!                           if (nAtom == 2) write(69,'(6f10.4)') log10(modulus(subcellCentre(thisOctal,subcell))/rCore), &
!                                log10(SUM(thisOctal%newAtomLevel(subcell,1,1:thisAtom(1)%nlevels-1)) /ntot), &
!                                log10(SUM(thisOctal%newAtomLevel(subcell,2,1:thisAtom(2)%nlevels-1)) / ntot)
!
!                           if (nAtom == 3) write(69,'(6f10.4)') log10(modulus(subcellCentre(thisOctal,subcell))/rCore), &
!                                log10(SUM(thisOctal%newAtomLevel(subcell,1,1:thisAtom(1)%nlevels-1)) /ntot), &
!                                log10(SUM(thisOctal%newAtomLevel(subcell,2,1:thisAtom(2)%nlevels-1)) / ntot), &
!                                log10(SUM(thisOctal%newAtomLevel(subcell,3,1:thisAtom(3)%nlevels-1)) / ntot), &
!                                log10(thisOctal%newAtomLevel(subcell,3,thisAtom(3)%nlevels) / ntot), &
!                                log10(thisOctal%rho(subcell))
!                           close(69) 
!                           itmp = itmp + 1
!                           write(tfilename,'(a,i2.2,a)') "jbar",itmp,".dat"
!                           open(69, file=tfilename, status="unknown",form="formatted")
!                           x1 = sqrt(max(0.d0,(1.d0 - source(1)%radius**2 / modulus(subcellCentre(thisOctal,subcell))**2)))
!                           w = 0.5d0*(1.d0 - x1)
!                           do i = 1, nFreq
!                              write(69,*) freq(i),thisOctal%jnucont(subcell,i), w*i_nu(source(1), freq(i), 1)
!                           enddo
!                           close(69)
!
!                        endif

                   endif

!                   write(*,*) "IMMEDIATELY replaceing populations...."
!                   thisOctal%atomLevel(subcell,:,:) = thisOctal%newatomLevel(subcell,:,:)


                endif

!                call locate(freq,nfreq,cspeed/6562.8d-8,ifreq)
!                nStar = boltzSahaGeneral(thisAtom(1), 6, thisOctal%ne(subcell), &
!                     dble(thisOctal%temperature(subcell))) * &
!                     thisOctal%newAtomLevel(subcell, 1, thisAtom(1)%nLevels)
!                write(31,*) real(r*1.e10), real(thisOctal%ne(subcell)), real(thisOctal%newAtomLevel(subcell,1,1:6)),&
!                     real(jnuCont(ifreq)),real(nstar)

!             enddo
!             close(31)
!             stop

             enddo
          end do
          !$OMP END DO
          !$OMP BARRIER

            deallocate(oldPops)
            deallocate(mainoldPops)
            deallocate(newPops)
            deallocate(dPops)
            deallocate(jnuCont)
            deallocate(ds)
            deallocate(phi)
            deallocate(rayDeltaV)
            deallocate(i0)
            deallocate(Hcol)
            deallocate(HeIcol)
            deallocate(HeIIcol)
            deallocate(sourceNumber)
            deallocate(cosTheta)
            deallocate(weightFreq)
            deallocate(weightOmega)
            deallocate(hitPhotosphere)
            deallocate(iCont)
            deallocate(position)
            deallocate(direction)
            deallocate(oldpops1, oldpops2, oldpops3, oldpops4)

          !$OMP END PARALLEL

!          if (doTuning) call tune(6, "One octal iteration")  ! start a stopwatch

#ifdef MPI

     write(*,*) myrankGlobal, " now waiting at barrier"
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
       if(my_rank == 0) write(*,*) "Updating MPI grids"

     call countVoxels(grid%octreeRoot,nOctal,nVoxels)
     allocate(tArrayd(1:nVoxels))
     allocate(tempArrayd(1:nVoxels))
     tArrayd = 0.d0
     tempArrayd = 0.d0
     do iAtom = 1, nAtom
       do i = 1, thisAtom(iAtom)%nLevels
         tArrayd = 0.d0
          call packAtomLevel(octalArray, nVoxels, tArrayd,  iAtom, i, ioctal_beg, ioctal_end)
          call MPI_ALLREDUCE(tArrayd,tempArrayd,nVoxels,MPI_DOUBLE_PRECISION,&
            MPI_SUM,MPI_COMM_WORLD,ierr)
            tArrayd = tempArrayd
          call unpackAtomLevel(octalArray, nVoxels, tArrayd, iAtom, i)
       enddo
     enddo
     do i = 1, nFreq
       tArrayd = 0.d0
       call packJnu(octalArray, nVoxels, tArrayd, i, ioctal_beg, ioctal_end)
       call MPI_ALLREDUCE(tArrayd,tempArrayd,nVoxels,MPI_DOUBLE_PRECISION,&
            MPI_SUM,MPI_COMM_WORLD,ierr)
       tArrayd = tempArrayd
       call unpackJnu(octalArray, nVoxels, tArrayd, i)
     enddo
     deallocate(tArrayd, tempArrayd)

       if(my_rank == 0) write(*,*) "Done updating"
#endif

          maxFracChange = -1.d30
          nInUse = 0 ; nConverged = 0
          call swapPops(grid%octreeRoot, gridtolerance, nInUse, nConverged)
          percentageConverged = 100.d0 * dble(nConverged)/dble(nInUse)
          write(ifilename,'(a,i2.2,a)') "fracchange",idump,".vtk"
          call writeVTKfile(grid,ifilename, valueTypeString = (/"adot"/))
          if (writeoutput) write(*,*) "Percentage converged: ",percentageConverged

          if (myRankIsZero) &
               call writeAmrGrid("atom_tmp.grid",.false.,grid)


          if (myRankisZero) then
             open(69, file="cmf_convergence.dat", status="old", position = "append", form="formatted")
             write(69,'(i10, 1p, e10.2, 0p,  f10.2, i10, l10)') &
                  nIter,  gridtolerance, real(percentageConverged), nRay, fixedRays
             close(69)
          endif

          if (percentageConverged > 99.d0) then
             gridConverged = .true.
          endif
!          gridconverged = .true.
!          write(*,*) "forcing convergence !!!!!!!!!!"


!         deallocate(ds, phi, i0, sourceNumber, cosTheta, hitPhotosphere, &
!               weightFreq, weightOmega, hcol, heicol, heiicol, iCont)
!         deallocate(position, direction, rayDeltaV)

          if (.not.gridConverged) then
             if (.not.fixedRays) nRay = nRay * 2
          endif
          if (nRay > maxRay) then
             nRay = maxRay
             call writeWarning("Maximum number of rays exceeded - capping")
          endif

!          gridconverged = .true.
!          write(*,*) "!!!!!!!!!!!!! Forcing exit after only one iteration"


          if (doTuning)  call tune(6, "One cmf iteration")  ! stop a stopwatch
       enddo
    enddo
666 continue

    call writeInfo( "ATOM loop done.")
  end subroutine atomLoop

  recursive  subroutine  swapPops(thisOctal, tolerance, nInuse, nConverged)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, j, iAtom
    real(double) :: maxFrac, temp, tolerance
    integer :: nInUse, nConverged
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call swapPops(child, tolerance, nInUse, nConverged)
                exit
             end if
          end do
       else

          if (.not.associated(thisOctal%adot)) then
             allocate(thisOctal%adot(1:thisOctal%maxChildren))
             thisOctal%adot = 0.d0
          endif
          
          if (thisOctal%inflow(subcell).and.(thisOctal%temperature(subcell) > 3000.)) then
             nInuse = nInuse + 1
             do iAtom = 1, size(thisOctal%newAtomLevel,2)
                maxFrac = -1.d30
                do j = 1 , 6
                   if (thisOctal%atomLevel(subcell,iAtom,j) /= 0.d0) then
                      temp = abs((thisOctal%newatomLevel(subcell,iAtom,j) - &
                           thisOctal%atomLevel(subcell,iAtom,j)) / &
                           thisOctal%atomLevel(subcell,iAtom,j))
                      thisOctal%adot(subcell) = temp
                      maxFrac = max(temp, maxFrac)
                   endif
                enddo
             enddo
             if (maxFrac < tolerance) then
                nConverged = nConverged + 1
             endif
          endif
          thisOctal%atomLevel(subcell,:,:) = &
               thisOctal%newatomLevel(subcell,:,:)
       endif
    enddo
  end subroutine swapPops


  recursive  subroutine  calcEtaLine(thisOctal, thisAtom, nAtom, iAtom, iTrans)
    type(MODELATOM) :: thisAtom(:)
    integer :: nAtom
    integer :: iTrans
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, iUpper, iAtom
    real(double) :: a, bul, blu
    real(double) :: etaLine

    a = 0.d0; blu = 0.d0; bul = 0.d0
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calcEtaLine(child, thisAtom, nAtom, iAtom, iTrans)
                exit
             end if
          end do
       else

          call returnEinsteinCoeffs(thisAtom(iatom), iTrans, a, Bul, Blu)
          iUpper = thisAtom(iatom)%iUpper(iTrans)
          etaLine = (hCgs/fourPi) * a * thisAtom(iatom)%transFreq(iTrans)
          etaLine = etaLine * thisOctal%atomLevel(subcell, iAtom,iUpper)
	  
          thisOctal%etaLine(subcell) = etaLine * 1.d10
          
       endif
    enddo
  end subroutine calcEtaLine

  recursive  subroutine  calcChiLine(thisOctal, thisAtom, nAtom, iAtom, iTrans)
    type(MODELATOM) :: thisAtom(:)
    integer :: nAtom
    integer :: iTrans
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, iUpper, iAtom
    real(double) :: a, bul, blu
    real(double) :: alphaNu
    integer :: ilower
    real(double) :: nlower, nupper
    a = 0.d0; blu = 0.d0; bul = 0.d0
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calcChiLine(child, thisAtom, nAtom, iAtom, iTrans)
                exit
             end if
          end do
       else

       iUpper = thisAtom(iatom)%iUpper(iTrans)
       iLower = thisAtom(iatom)%iLower(iTrans)
          nLower = thisOctal%atomLevel(subcell,iAtom, iLower)
          nUpper = thisOctal%atomLevel(subcell,iAtom, iUpper)
          alphanu = (hCgs*thisAtom(iatom)%transFreq(iTrans)/fourPi)

          call returnEinsteinCoeffs(thisAtom(iatom), iTrans, a, Bul, Blu)

          alphanu = alphanu * (nLower * Blu - nUpper * Bul) !/thisAtom(iatom)%transFreq(iTrans)
          thisOctal%chiLine(subcell) = alphanu * 1.d10
          
       endif
    enddo
  end subroutine calcChiLine


  recursive  subroutine calcContinuumOpacities(thisOctal, thisAtom, nAtom, freq)
    type(MODELATOM) :: thisAtom(:)
    real(double) :: freq
    real(double) :: nstar(10,50)
    integer :: nAtom
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, iAtom, iLevel

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calcContinuumOpacities(child, thisAtom, nAtom, freq)
                exit
             end if
          end do
       else

          do iatom = 1, size(thisAtom)
             do iLevel = 1, thisAtom(iatom)%nLevels-1
                nStar = boltzSahaGeneral(thisAtom(iAtom), iLevel, thisOctal%ne(subcell), &
                     dble(thisOctal%temperature(subcell))) * &
                     thisOctal%atomLevel(subcell, iatom, thisAtom(iAtom)%nLevels)
             enddo
          enddo

          thisOctal%kappaAbs(subcell, 1) = bfOpacity(freq, nAtom, thisAtom, thisOctal%atomLevel(subcell,:,:), &
                     thisOctal%ne(subcell), dble(thisOctal%temperature(subcell))) * 1.d10

          thisOctal%kappaSca(subcell, 1) = thisOctal%ne(subcell) * sigmaE * 1.d10
          thisOctal%etaCont(subcell)  = bfEmissivity(freq, nAtom, thisAtom,  thisOctal%atomLevel(subcell,:,:), nstar, &
                     dble(thisOctal%temperature(subcell)), thisOctal%ne(subcell)) * 1.d10
          
       endif
    enddo
  end subroutine calcContinuumOpacities

  recursive subroutine  allocateLevels(grid, thisOctal, nAtom, thisAtom, nRBBTrans, nFreq, ionized)
    use stateq_mod, only : z_hi
    type(GRIDTYPE) :: grid
    type(MODELATOM) :: thisAtom(:)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    typE(VECTOR) :: rvec
    real(double) :: ne, n_h, ntot, phit, t
    real(double), parameter  :: CI = 2.07d-16   ! in cgs units
    integer :: i, subcell
    integer :: nAtom
    integer :: nRBBTrans
    integer :: iatom
    logical :: ionized
    integer :: nFreq
    

    if ( thisOctal%nChildren > 0) then
       do i = 1, thisOctal%nChildren
          child => thisOctal%child(i)
          call  allocateLevels(grid, child, nAtom, thisAtom, nRBBTrans, nfreq, ionized)
       end do
    end if
    call allocateAttribute(thisOctal%atomAbundance,thisOctal%maxchildren, nAtom)
    thisOctal%atomAbundance(:, 1) = 1.d0 / mHydrogen
    if (nAtom > 1) then
       thisOctal%atomAbundance(:, 2:nAtom) =  0.2d0 / mHydrogen !assume higher atoms are helium
    endif
    call allocateAttribute(thisOctal%microturb, thisOctal%maxChildren)
    call allocateAttribute(thisOctal%etaLine, thisOctal%maxChildren)
    call allocateAttribute(thisOctal%chiLine, thisOctal%maxChildren)
    call allocateAttribute(thisOctal%ne, thisOctal%maxChildren)


    call allocateAttribute(thisOctal%atomLevel, thisOctal%maxChildren, nAtom, maxval(thisAtom(1:nAtom)%nLevels))
    call allocateAttribute(thisOctal%newatomLevel, thisOctal%maxChildren, nAtom, maxval(thisAtom(1:nAtom)%nLevels))


    call allocateAttribute(thisOctal%jnuLine, thisOctal%maxChildren, nRBBTrans)
    call allocateAttribute(thisOctal%jnuCont, thisOctal%maxChildren, nFreq)

    call allocateAttribute(thisOctal%microturb, thisOctal%maxChildren)

    if (ionized) then
       thisOctal%ne = thisOctal%rho(1:SIZE(thisOctal%ne))/mHydrogen
    else
       thisOctal%ne = tiny(1.d-30 * thisOctal%rho(1:SIZE(thisOctal%ne))/mHydrogen)
    endif


    do subcell = 1, thisOctal%maxChildren
       t = thisOctal%temperature(subcell)
       rVec = subcellCentre(thisOctal, subcell)
!       write(*,*) " t ",t , " r ",modulus(rVec)/grid%rcore
       N_H = thisOctal%rho(subcell)/mHydrogen  ! number density of HI plus number density of HII
       if (real(hydE0eV,kind=double)/(kev*T) < 60.d0) then 
          phiT = CI*Z_HI(10,T)*(T**(-1.5))*EXP(real(hydE0eV,kind=double)/(kev*T))

          ! Solving for phi(T)*ne^2 + 2ne -nTot =0 and ne+N_H = nTot for ne where
          ! nTot is the number density of particles includeing all species.
          ! ==> phi(T)*ne^2 + ne - N_H =0
          ! Th physical solution  is chosen out of two ...  
          !    Ne = (sqrt(nTot*phiT+1.0_db) -1.0_db)/phiT
          Ne = (sqrt(4.0_db*N_H*phiT+1.0_db) -1.0_db)/(2.0_db*phiT)
       else 
          ne = tiny(ne)
       endif
       nTot = Ne + N_H
       !    Ne = min(Ne, nTot)     ! to avoid unphysical solution.
       !    if (Ne<=0) Ne =nTot   ! to avoid unphysical solution.
       if (Ne<=0) Ne =1.0d-40 ! to avoid unphysical solution.
       thisOctal%ne(subcell) = ne
       !          write(*,*) "ne ",thisOctal%ne(subcell),4.0_db*N_H*phiT+1.0_db, t,phit,n_h
       do iAtom = 1, nAtom
          thisOctal%atomLevel(subcell,iAtom,thisAtom(iAtom)%nLevels) = thisOctal%rho(subcell) &
               * thisOctal%atomAbundance(subcell, iAtom)/ thisAtom(iatom)%mass

          do i = 1, thisAtom(iatom)%nLevels-1
             thisOctal%atomLevel(subcell, iAtom,i) = boltzSahaGeneral(thisAtom(iAtom), i, &
                  thisOctal%ne(subcell), &
                  dble(thisOctal%temperature(subcell))) * &
                  thisOctal%AtomLevel(subcell, iatom, thisAtom(iAtom)%nLevels)
          enddo
       enddo
    enddo
    thisOctal%newAtomLevel = thisOctal%atomLevel

  end subroutine allocateLevels
  

  subroutine randomRayDirection(probTowardsSource, point, source, nSource, direction, weight)
    type(VECTOR) :: point, direction, toStar
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    real(double) :: probTowardsSource, weight
    integer :: i
    real(double) :: chanceSource
    real(double) :: omegaSubtendedBySource
    real(double) :: theta, dOmega, r, cosTheta, ang
    logical :: hitCore
    integer :: nBug

    omegaSubtendedBySource = 0.d0
    do i = 1 , nSource
       theta = asin(source(i)%radius/modulus(point-source(i)%position))
       domega = twoPi * (1.d0 - cos(theta))
       omegaSubtendedBySource = omegaSubtendedbySource + dOmega
    enddo

    chanceSource = omegaSubtendedBySource / fourPi

    weight = 1.d0
 
    call randomNumberGenerator(getDouble=r)

    if (r < probTowardsSource) then
       weight = chanceSource/probTowardsSource
       call randomNumberGenerator(getDouble=r)
       i = int(r*real(nSource))+1
       toStar = source(i)%position - point
       call normalize(toStar)
       theta = asin(min(1.d0, max(-1.d0,source(i)%radius/modulus(point-source(i)%position))))
       cosTheta = cos(theta)
       hitCore = .false.
       nbug = 0
       do while(.not.hitCore)
          direction = randomUnitVector()
          ang = acos(min(1.d0, max(-1.d0, direction.dot.toStar)))
          if (ang < theta) hitCore = .true.
          nbug = nbug  + 1
          if (nbug > 10000000) then
             write(*,*) "bug in direction to star"
             direction = toStar
             exit
          endif
       enddo
    else
       weight = (1.d0-chanceSource)/(1.d0-probTowardsSource)
       hitCore = .true.
       do while(hitCore)
          direction = randomUnitVector()
          do i = 1, nSource
             call distanceToSource(source, i, point, direction, hitCore, r)
          enddo
       end do
    endif
  end subroutine randomRayDirection



  function intensityAlongRay(position, direction, grid, thisAtom, nAtom, iAtom, iTrans, deltaV, source, nSource, &
       nFreq, freqArray, forceFreq, occultingDisc) result (i0)
    use input_variables, only : lineOff,  mie
    use amr_mod, only: distanceToGridFromOutside, returnKappa
    use utils_mod, only : findIlambda
    use atom_mod, only : bnu
    logical ::     justPhotosphere
    type(VECTOR) :: position, direction, pvec, photoDirection
    type(GRIDTYPE) :: grid
    logical, optional :: occultingDisc
    integer :: nSource
    real(double) :: freqArray(:)
    real(double), optional :: forceFreq
    integer :: nFreq
    type(SOURCETYPE) :: source(:)
    real(double) :: transitionFreq
    integer :: iAtom, nAtom
    type(MODELATOM) :: thisAtom(:)
    real(double) :: disttoGrid
    integer :: itrans, iFreq
    real(double) :: totDist
    logical :: hitSource
    real(double) :: i0
    type(OCTAL), pointer :: thisOctal, startOctal !, endOctal
    !    integer :: endSubcell
    integer :: subcell
    real(double) :: costheta
    type(VECTOR) :: currentPosition, thisPosition, thisVel, oldposition
    type(VECTOR) :: rayVel, startVel, endVel, endPosition !, rvec
    real(double) :: alphanu, snu, jnu
    integer :: iLower , iUpper
    real(double) :: dv, deltaV
    integer :: i, icount
    real(double) :: distArray(1000), tval
    integer :: nTau
    real(double) :: nLower, nUpper
    real(double) :: dTau, etaline, tau
    real(double) :: intensityIntegral
    real(double) :: dvAcrossCell
    real(double) :: dv1, dv2
    real(double) :: a, bul, blu
    integer :: nHatom,nHeIAtom, nHeIIAtom
    real(double) :: distToSource,disttoDisc
    integer :: sourcenumber
    integer :: iElement
    logical :: endLoopAtPhotosphere
    real(double) :: nstar(10,50), rhoCol
    real(double) :: bfOpac, bfEmiss, x1, x2, fac
    real(double) :: dustOpac, dustEmiss
    integer :: ilambda
    real(double) :: transitionLambda, kappaSca, kappaAbs, kappaExt
    integer :: j, k
    logical :: passThroughResonance, velocityIncreasing, ok

    i0 = tiny(i0)
    justPhotosphere = .false.

    hitsource = .false.; disttosource = 0.d0; sourceNumber = 0
    a = 0.d0; blu = 0.d0; bul = 0.d0
    nHAtom = 0
    nHeIAtom = 0
    nHeIIAtom = 0
    do i = 1, nAtom
       if (thisAtom(i)%name == "HI") nHatom = i
       if (thisAtom(i)%name == "HeI") nHeIatom = i
       if (thisAtom(i)%name == "HeII") nHeIIatom = i
    enddo

    transitionFreq = thisAtom(iAtom)%transFreq(iTrans)

    if (PRESENT(forceFreq)) then
       transitionFreq = forceFreq
    endif

    call locate(freqArray, nFreq, transitionFreq, iFreq)

    transitionLambda = (cSpeed/transitionFreq) /angstromTocm
    iLambda = findIlambda(real(transitionLambda), grid%lamArray, grid%nLambda, ok)


    distToGrid = distanceToGridFromOutside(grid, position, direction)

    currentposition = position
    distToDisc = 1.d30
    if ((currentposition.dot.direction) < 0.d0) then
       if (direction%z /= 0.d0) then
          distToDisc = abs(currentPosition%z/direction%z)
          oldPosition = currentPosition + distToDisc * direction
          if (sqrt(oldPosition%x**2 + oldPosition%y**2) < (2.d0*grid%octreeRoot%subcellSize)) then
             distToDisc = 1.d30
          endif
       else
          distToDisc = 1.d30
       endif
    endif
    if (present(OccultingDisc)) then
       if (occultingDisc) then
          if (distToDisc < distToGrid) then
             i0 = tiny(i0)
             goto 666
          endif
       endif
    endif

    if (distToGrid > 1.e29) then
!              write(*,*) "ray does not intersect grid",position,direction
       i0 = tiny(i0)
       goto 666
    endif

    iUpper = thisAtom(iAtom)%iUpper(iTrans)
    iLower = thisAtom(iAtom)%iLower(iTrans)

    currentPosition = position + (distToGrid + 1.d-3*grid%halfSmallestSubcell) * direction
 
 

    if (.not.inOctal(grid%octreeRoot, currentPosition)) then
       write(*,*) "initial position not in grid"
       write(*,*) "curre pos",currentPosition
       write(*,*) "dir",direction
       write(*,*) "pos",position
       write(*,*) "modulsu",modulus(currentPosition - position)
       stop
    endif

    totDist = 0.d0
    call distanceToSource(source, nSource, currentposition, direction, hitSource, disttoSource, sourcenumber)

    if ((.not.hitsource).and.justPhotosphere) goto 666
    if (hitSource) then
       pVec = (currentposition + (direction * distToSource) - source(sourceNumber)%position)
       call normalize(pVec)
       cosTheta = -1.d0*(pVec.dot.direction)
       photoDirection = pVec
       call normalize(photoDirection)
    endif


    !    write(*,*) "currentposition",sqrt(currentPosition%x**2+currentPosition%y**2),currentPosition%z, &
    !         inOctal(grid%octreeRoot, currentPosition),distTogrid
    i0 = tiny(i0)!0.d0
    intensityIntegral = 0.0
    tau = 0.d0
    rayVel = VECTOR(0.d0, 0.d0, 0.d0)

    thisOctal => grid%octreeRoot
    icount = 0
    rhoCol = 0.d0
    endLoopAtPhotosphere = .false.

    !    if (hitSource) endLoopAtphotosphere = .true.

    !    write(*,*) lineoff,hitsource,endloopatphotosphere



       do while(inOctal(grid%octreeRoot, currentPosition).and.(.not.endloopAtPhotosphere))
          icount = icount + 1 
          call findSubcellLocal(currentPosition, thisOctal, subcell)

          !       rVec = subcellCentre(thisOctal,subcell)

          call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal)


          if ((totDist + tval) > distTosource) then
             tVal = distToSource - totDist
             endLoopAtPhotosphere = .true.
          endif


          if (.not.lineOff) then
             startVel = amrGridVelocity(grid%octreeRoot, currentPosition, startOctal = thisOctal, actualSubcell = subcell) 

             endPosition = currentPosition + tval * direction
             endVel = amrGridVelocity(grid%octreeRoot, endPosition)

             dv1 = deltaV + (startVel .dot. direction)
             dv2 = deltaV + (endVel .dot. direction)

             dvAcrossCell = abs((dv2-dv1) / thisOctal%microturb(subcell))

             distArray(1) = 0.d0
             distArray(2) = tVal
             nTau = 2

             passThroughResonance =.false.


             if (dv1*dv2 < 0.d0) passThroughResonance = .true.

             if (modulus(endVel)==0.d0) passThroughResonance = .false.

             if (passthroughresonance.or.(min(abs(dv1),abs(dv2)) < 4.d0*thisOctal%microturb(subcell))) then

                if (dv1 <= dv2) then
                   velocityIncreasing = .true.
                else
                   velocityIncreasing = .false.
                endif

                if ( (dv2 - dv1) /= 0.d0) then
                   if (velocityIncreasing) then
                      x1 = tval*(-4.d0*thisOctal%microturb(subcell) - dv1)/(dv2 - dv1)
                      x2 = tval*(+4.d0*thisOctal%microturb(subcell) - dv1)/(dv2 - dv1)
                   else
                      x1 = tval*(+4.d0*thisOctal%microturb(subcell) - dv1)/(dv2 - dv1)
                      x2 = tval*(-4.d0*thisOctal%microturb(subcell) - dv1)/(dv2 - dv1)
                   endif
                else
                   x1 = 0.d0
                   x2 = tVal
                endif

                if (x1 > x2) then
                   fac = x1
                   x1 = x2
                   x2 = fac
                endif
                x1 = max(0.d0, x1)
                x2 = min(x2, tVal)
                nTau = 100
                distArray(1) = 0.d0
                do i = 1, nTau-1
                   distArray(i+1) = x1 + (x2 - x1)*dble(i-1)/dble(ntau-2)
                enddo
                distArray(nTau) = tVal
             endif

             if (.not.thisOctal%inflow(subcell)) then
                distArray(1) = 0.d0
                distArray(2) = tVal
                nTau = 2
             endif
          else
             distArray(1) = 0.d0
             distArray(2) = tVal
             nTau = 2
          endif


          bfOpac = 0.d0
          bfEmiss = 0.d0

          do i = 2, nTau

             startOctal => thisOctal
             thisPosition = currentPosition + distArray(i)*direction




             if (.not.lineoff) then
                thisVel = amrGridVelocity(grid%octreeRoot, thisPosition, startOctal = startOctal, actualSubcell = subcell) 
                thisVel= thisVel - rayVel
                dv = (thisVel .dot. direction) + deltaV
                call returnEinsteinCoeffs(thisAtom(iAtom), iTrans, a, Bul, Blu)
                alphanu = (hCgs*thisAtom(iAtom)%transFreq(iTrans)/fourPi) * &
                     phiProf(dv, thisOctal%microturb(subcell)) /thisAtom(iAtom)%transFreq(iTrans)          
                iUpper = thisAtom(iAtom)%iUpper(iTrans)
                iLower = thisAtom(iAtom)%iLower(iTrans)
                nLower = thisOctal%atomLevel(subcell,iAtom, iLower)
                nUpper = thisOctal%atomLevel(subcell,iAtom, iUpper)
                alphanu = alphanu * (nLower * Blu - nUpper * Bul)
             else
                alphanu = 0.d0
             endif


             if (i == 2) then
                do k = 1, nAtom
                   do j = 1, thisAtom(k)%nLevels - 1
                      nStar(k,j) = BoltzSahaGeneral(thisAtom(k), j, thisOctal%ne(subcell), &
                           dble(thisOctal%temperature(subcell))) * &
                           Thisoctal%atomlevel(subcell, k,thisAtom(k)%nLevels)
                   enddo
                enddo
                bfOpac = bfOpacity(transitionFreq, nAtom, thisAtom, thisOctal%atomLevel(subcell,:,:), &
                     thisOctal%ne(subcell), dble(thisOctal%temperature(subcell)))
                bfEmiss = bfEmissivity(transitionFreq, nAtom, thisAtom,  thisOctal%atomLevel(subcell,:,:), nstar, &
                     dble(thisOctal%temperature(subcell)), thisOctal%ne(subcell))
                dustOpac = 0.d0
                dustEmiss = 0.d0
                if (mie) then
                   call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, kappaSca=kappaSca, kappaAbs=kappaAbs)
                   kappaExt = kappaAbs + kappaSca
                   dustOpac = kappaExt/1.d10
                   dustEmiss = kappaAbs * bnu(transitionFreq, dble(thisOctal%temperature(subcell)))/1.d10
                   dustEmiss = dustEmiss + kappaSca * thisOctal%meanIntensity(subcell)/1.d10
                endif
             endif


             if (.not.lineoff) then
                etaLine = hCgs * a * transitionFreq
                etaLine = etaLine * thisOctal%atomLevel(subcell, iAtom, iUpper)
                jnu = (etaLine/fourPi) * phiProf(dv, thisOctal%microturb(subcell)) /transitionFreq
             else
                jnu = 0.d0
                etaline = 0.d0
             endif

             if (associated(thisOctal%fixedTemperature)) then
                if (.not.thisOctal%fixedTemperature(subcell)) then
                   jnu = 0.d0
                   etaline = 0.d0
                endif
             endif

             if (thisOctal%rho(subcell) > 0.1d0) then ! opaque disc
                bfOpac = 1.d30
                bfEmiss = 0.d0
                alphanu = 1.d30
                etaLine = 0.d0
                jnu = 0.d0
                snu = 0.d0
             endif

             alphanu = alphanu + bfOpac + dustOpac

             ! add continuous bf and ff emissivity of hydrogen

             jnu = jnu + bfEmiss + dustEmiss

             if (alphanu /= 0.d0) then
                snu = jnu/alphanu
             else
                snu = tiny(snu)
             endif


             dTau = alphaNu *  (distArray(i)-distArray(i-1)) * 1.d10

             !          write(*,*) iCount,ntau,totdist, &
             !               modulus(currentPosition)/(grid%rCore),tau,modulus(thisOctal%velocity(subcell))*cspeed/1.d5,i0,&
             !               passthroughresonance, dv*cspeed/1.e5, phiProf(dv, thisOctal%microturb(subcell)), &
             !               snu,dtau,tau

!                       write(*,'(i4,i4,l3,10(1pe12.3))') iCount, ntau, passThroughResonance, &
!                            dv1*cspeed/1.e5,dv2*cspeed/1.e5,dv*cspeed/1.e5, &
!                          i0,tau, jnu,alphanu,snu, nlower,nupper

             if (thisOctal%inflow(subcell)) then
                i0 = i0 +  exp(-tau) * (1.d0-exp(-dtau))*snu
             endif
             tau = tau + dtau
          enddo
          rhoCol = rhoCol + distArray(ntau)*thisOctal%rho(subcell)*1.d10
          oldPosition = currentPosition
          currentPosition = currentPosition + (distArray(ntau)+1.d-3*grid%halfSmallestSubcell) * direction
          totdist = totdist + (distArray(ntau)+1.d-3*grid%halfSmallestSubcell)

          if (PRESENT(occultingDisc)) then
             if (occultingDisc) then
                if ((oldPosition%z >= 0.d0).and.(currentPosition%z < 0.d0)) then
                   if (sqrt(currentPosition%x**2 + currentPosition%y**2) > (2.d0*grid%octreeRoot%subcellSize)) then
                      goto 666
                   endif
                endif
             endif
          endif
       enddo

    if (endLoopAtPhotosphere) then

       iElement = getElement(source(sourcenumber)%surface, photoDirection)

       i0 = i0 + i_nu(source(sourceNumber), transitionFreq, iElement)*exp(-tau)

    endif
666 continue 

  end function intensityAlongRay


  function intensityAlongRayGeneric(position, direction, grid,deltaV, source, nSource, &
      forceFreq, occultingDisc) result (i0)
    use input_variables, only : lineOff,  mie, lamLine
    use amr_mod, only: distanceToGridFromOutside, returnKappa
    use utils_mod, only : findIlambda
    use atom_mod, only : bnu
    logical ::     justPhotosphere
    type(VECTOR) :: position, direction, pvec, photoDirection
    type(GRIDTYPE) :: grid
    logical, optional :: occultingDisc
    integer :: nSource
    real(double), optional :: forceFreq
    type(SOURCETYPE) :: source(:)
    real(double) :: transitionFreq
    real(double) :: disttoGrid
    real(double) :: totDist
    logical :: hitSource
    real(double) :: i0
    type(OCTAL), pointer :: thisOctal, startOctal !, endOctal
    !    integer :: endSubcell
    integer :: subcell
    real(double) :: costheta
    type(VECTOR) :: currentPosition, thisPosition, thisVel, oldposition
    type(VECTOR) :: rayVel, startVel, endVel, endPosition !, rvec
    real(double) :: alphanu, snu, jnu
    real(double) :: dv, deltaV
    integer :: i, icount
    real(double) :: distArray(1000), tval
    integer :: nTau
    real(double) :: dTau, etaline, tau
    real(double) :: intensityIntegral
    real(double) :: dvAcrossCell
    real(double) :: dv1, dv2
    real(double) :: a, bul, blu
    real(double) :: distToSource,disttoDisc
    integer :: sourcenumber
    integer :: iElement
    logical :: endLoopAtPhotosphere
    real(double) ::  rhoCol
    real(double) :: bfOpac, bfEmiss, x1, x2, fac
    real(double) :: dustOpac, dustEmiss
    integer :: ilambda
    real(double) :: transitionLambda, kappaSca, kappaAbs, kappaExt
    logical :: passThroughResonance, velocityIncreasing, ok

    i0 = tiny(i0)
    justPhotosphere = .false.

    hitsource = .false.; disttosource = 0.d0; sourceNumber = 0
    a = 0.d0; blu = 0.d0; bul = 0.d0


    if (PRESENT(forceFreq)) then
       transitionFreq = forceFreq
    endif


    transitionLambda = lamLine 

    transitionFreq = cSpeed / (lamLine * angstromtocm)
    iLambda = findIlambda(real(transitionLambda), grid%lamArray, grid%nLambda, ok)


    distToGrid = distanceToGridFromOutside(grid, position, direction)

    currentposition = position
    distToDisc = 1.d30
    if ((currentposition.dot.direction) < 0.d0) then
       if (direction%z /= 0.d0) then
          distToDisc = abs(currentPosition%z/direction%z)
          oldPosition = currentPosition + distToDisc * direction
          if (sqrt(oldPosition%x**2 + oldPosition%y**2) < (2.d0*grid%octreeRoot%subcellSize)) then
             distToDisc = 1.d30
          endif
       else
          distToDisc = 1.d30
       endif
    endif
    if (present(OccultingDisc)) then
       if (occultingDisc) then
          if (distToDisc < distToGrid) then
             i0 = tiny(i0)
             goto 666
          endif
       endif
    endif

    if (distToGrid > 1.e29) then
!              write(*,*) "ray does not intersect grid",position,direction
       i0 = tiny(i0)
       goto 666
    endif

    currentPosition = position + (distToGrid + 1.d-3*grid%halfSmallestSubcell) * direction
 
 

    if (.not.inOctal(grid%octreeRoot, currentPosition)) then
       write(*,*) "initial position not in grid"
       write(*,*) "curre pos",currentPosition
       write(*,*) "dir",direction
       write(*,*) "pos",position
       write(*,*) "modulsu",modulus(currentPosition - position)
       stop
    endif

    totDist = 0.d0
    call distanceToSource(source, nSource, currentposition, direction, hitSource, disttoSource, sourcenumber)

    if ((.not.hitsource).and.justPhotosphere) goto 666
    if (hitSource) then
       pVec = (currentposition + (direction * distToSource) - source(sourceNumber)%position)
       call normalize(pVec)
       cosTheta = -1.d0*(pVec.dot.direction)
       photoDirection = pVec
       call normalize(photoDirection)
    endif


    !    write(*,*) "currentposition",sqrt(currentPosition%x**2+currentPosition%y**2),currentPosition%z, &
    !         inOctal(grid%octreeRoot, currentPosition),distTogrid
    i0 = tiny(i0)!0.d0
    intensityIntegral = 0.0
    tau = 0.d0
    rayVel = VECTOR(0.d0, 0.d0, 0.d0)

    thisOctal => grid%octreeRoot
    icount = 0
    rhoCol = 0.d0
    endLoopAtPhotosphere = .false.

    !    if (hitSource) endLoopAtphotosphere = .true.

    !    write(*,*) lineoff,hitsource,endloopatphotosphere



       do while(inOctal(grid%octreeRoot, currentPosition).and.(.not.endloopAtPhotosphere))
          icount = icount + 1 
          call findSubcellLocal(currentPosition, thisOctal, subcell)

          !       rVec = subcellCentre(thisOctal,subcell)

          call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal)


          if ((totDist + tval) > distTosource) then
             tVal = distToSource - totDist
             endLoopAtPhotosphere = .true.
          endif


          if (.not.lineOff) then
             startVel = amrGridVelocity(grid%octreeRoot, currentPosition, startOctal = thisOctal, actualSubcell = subcell, linearInterp=.false.) 

             endPosition = currentPosition + tval * direction
             endVel = amrGridVelocity(grid%octreeRoot, endPosition, linearInterp=.false.)

             dv1 = deltaV + (startVel .dot. direction)
             dv2 = deltaV + (endVel .dot. direction)

             dvAcrossCell = abs((dv2-dv1) / thisOctal%microturb(subcell))

             distArray(1) = 0.d0
             distArray(2) = tVal
             nTau = 2

             passThroughResonance =.false.


             if (dv1*dv2 < 0.d0) passThroughResonance = .true.

             if (modulus(endVel)==0.d0) passThroughResonance = .false.

             if (passthroughresonance.or.(min(abs(dv1),abs(dv2)) < 4.d0*thisOctal%microturb(subcell))) then

                if (dv1 <= dv2) then
                   velocityIncreasing = .true.
                else
                   velocityIncreasing = .false.
                endif

                if ( (dv2 - dv1) /= 0.d0) then
                   if (velocityIncreasing) then
                      x1 = tval*(-4.d0*thisOctal%microturb(subcell) - dv1)/(dv2 - dv1)
                      x2 = tval*(+4.d0*thisOctal%microturb(subcell) - dv1)/(dv2 - dv1)
                   else
                      x1 = tval*(+4.d0*thisOctal%microturb(subcell) - dv1)/(dv2 - dv1)
                      x2 = tval*(-4.d0*thisOctal%microturb(subcell) - dv1)/(dv2 - dv1)
                   endif
                else
                   x1 = 0.d0
                   x2 = tVal
                endif

                if (x1 > x2) then
                   fac = x1
                   x1 = x2
                   x2 = fac
                endif
                x1 = max(0.d0, x1)
                x2 = min(x2, tVal)
                nTau = 100
                distArray(1) = 0.d0
                do i = 1, nTau-1
                   distArray(i+1) = x1 + (x2 - x1)*dble(i-1)/dble(ntau-2)
                enddo
                distArray(nTau) = tVal
             endif

             if (.not.thisOctal%inflow(subcell)) then
                distArray(1) = 0.d0
                distArray(2) = tVal
                nTau = 2
             endif
          else
             distArray(1) = 0.d0
             distArray(2) = tVal
             nTau = 2
          endif


          bfOpac = 0.d0
          bfEmiss = 0.d0

          do i = 2, nTau

             startOctal => thisOctal
             thisPosition = currentPosition + distArray(i)*direction




             if (.not.lineoff) then
                thisVel = amrGridVelocity(grid%octreeRoot, thisPosition, startOctal = startOctal, actualSubcell = subcell, linearInterp=.false.) 
                thisVel= thisVel - rayVel
                dv = (thisVel .dot. direction) + deltaV
                alphanu = thisOctal%chiLine(subcell) * phiProf(dv, thisOctal%microturb(subcell))
             else
                alphanu = 0.d0
             endif


             if (i == 2) then
                bfOpac = thisOctal%kappaAbs(subcell,1)
                bfEmiss = thisOctal%etaCont(subcell)
                dustOpac = 0.d0
                dustEmiss = 0.d0
                if (mie) then
                   call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, kappaSca=kappaSca, kappaAbs=kappaAbs)
                   kappaExt = kappaAbs + kappaSca
                   dustOpac = kappaExt
                   dustEmiss = kappaAbs * bnu(transitionFreq, dble(thisOctal%temperature(subcell)))
                   dustEmiss = dustEmiss + kappaSca * thisOctal%meanIntensity(subcell)
                endif
             endif


             if (.not.lineoff) then
                etaLine = thisOctal%etaLine(subcell)
                jnu = etaLine * phiProf(dv, thisOctal%microturb(subcell))
             else
                jnu = 0.d0
                etaline = 0.d0
             endif

             if (associated(thisOctal%fixedTemperature)) then
                if (.not.thisOctal%fixedTemperature(subcell)) then
                   jnu = 0.d0
                   etaline = 0.d0
                endif
             endif

             if (thisOctal%rho(subcell) > 0.1d0) then ! opaque disc
                bfOpac = 1.d30
                bfEmiss = 0.d0
                alphanu = 1.d30
                etaLine = 0.d0
                jnu = 0.d0
                snu = 0.d0
             endif

             alphanu = alphanu + bfOpac + dustOpac

             ! add continuous bf and ff emissivity of hydrogen

             jnu = jnu + bfEmiss + dustEmiss

             if (alphanu /= 0.d0) then
                snu = jnu/alphanu
             else
                snu = tiny(snu)
             endif


             dTau = alphaNu *  (distArray(i)-distArray(i-1))


             if (thisOctal%inflow(subcell)) then
                i0 = i0 +  exp(-tau) * (1.d0-exp(-dtau))*snu
             endif
             tau = tau + dtau
          enddo
          rhoCol = rhoCol + distArray(ntau)*thisOctal%rho(subcell)*1.d10
          oldPosition = currentPosition
          currentPosition = currentPosition + (distArray(ntau)+1.d-3*grid%halfSmallestSubcell) * direction
          totdist = totdist + (distArray(ntau)+1.d-3*grid%halfSmallestSubcell)

          if (PRESENT(occultingDisc)) then
             if (occultingDisc) then
                if ((oldPosition%z >= 0.d0).and.(currentPosition%z < 0.d0)) then
                   if (sqrt(currentPosition%x**2 + currentPosition%y**2) > (2.d0*grid%octreeRoot%subcellSize)) then
                      goto 666
                   endif
                endif
             endif
          endif
       enddo

    if (endLoopAtPhotosphere) then

       iElement = getElement(source(sourcenumber)%surface, photoDirection)

       i0 = i0 + i_nu(source(sourceNumber), transitionFreq, iElement)*exp(-tau)
    endif
666 continue 

  end function intensityAlongRayGeneric

  
  subroutine calculateAtomSpectrum(grid, thisAtom, nAtom, iAtom, iTrans, viewVec, distance, source, nsource, nfile, &
       totalFlux, forceLambda, occultingDisc)
    use input_variables, only : vturb, lineoff, nv, calcDataCube, lamLine, cmf
    use messages_mod, only : myRankIsZero
    use datacube_mod, only: DATACUBE, freedatacube
    use modelatom_mod, only : identifyTransitionCmf
#ifdef USECFITSIO
    use input_variables, only : dataCubeFilename
    use datacube_mod, only : writedataCube
#endif
#ifdef MPI
    include 'mpif.h'
#endif
    logical, optional :: occultingDisc
    type(GRIDTYPE) :: grid
    type(MODELATOM) :: thisAtom(:)
    integer :: nSource
    real(double), optional :: forceLambda
    integer :: nFile
    type(SOURCETYPE) :: source(:)
    integer :: nAtom, iAtom
    real(double) :: distance, totalFlux
    integer :: itrans
    integer :: nRay,iray1,iray2
    integer, parameter :: maxRay  = 1000000
    type(VECTOR),allocatable :: rayPosition(:)
    real(double),allocatable :: da(:), dOmega(:)
    type(VECTOR) :: viewVec
    real(double) :: deltaV
    integer :: iv, iray
    real(double) :: i0
    real(double), allocatable :: vArray(:), spec(:)
    integer :: iv1, iv2, i
    character(len=30) :: plotfile,message
    type(DATACUBE) :: cube
    integer :: nFreqArray
    real(double) :: totalOmega
    integer, parameter :: maxFreq = 2000
    real(double) :: freqArray(maxFreq), broadBandFreq, transitionFreq
    logical :: doCube, doSpec

#ifdef MPI
    ! For MPI implementations
    integer       ::   my_rank        ! my processor rank
    integer       ::   np             ! The number of processes
    integer       ::   ierr           ! error flag
    real(double), allocatable :: tempArray(:)



    ! FOR MPI IMPLEMENTATION=======================================================
    !  Get my process rank # 
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  
    ! Find the total # of precessor being used in this run
    call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
#endif

    if (cmf) then
       call identifyTransitionCmf(dble(lamLine), thisAtom, iAtom, iTrans)
       transitionFreq = thisAtom(iAtom)%transfreq(iTrans)
       call calcChiLine(grid%octreeRoot, thisAtom, nAtom, iAtom, iTrans)
       call calcEtaLine(grid%octreeRoot, thisAtom, nAtom, iAtom, iTrans)
       call calcContinuumOpacities(grid%octreeRoot, thisAtom, nAtom, transitionfreq)
    endif

    call setMicroturb(grid%octreeRoot, dble(vTurb))
    call writeVTKfile(grid,"eta.vtk", valueTypeString = (/"etaline   ","chiline   ","sourceline"/))


    doCube = calcDataCube
    doSpec = .false.

    broadBandFreq = 1.d15
    if (PRESENT(forceLambda)) then
       broadBandFreq = cSpeed/(forceLambda * angstromToCm)
       write(message,'(a,f7.1)') "Calculating flux at ",forceLambda
       call writeInfo(message)
       lineoff = .true.
       nv = 1
    endif

    if (lineoff) call writeWarning("Line transfer switched off")

    cube%label = " "
    freqArray = 0.d0; nFreqArray = 0
    call createContFreqArray(nFreqArray, freqArray, nAtom, thisAtom, nsource, source)


    if (myRankIsZero.and.(.not.PRESENT(forcelambda))) &
         write(*,*) "Calculating spectrum for: ",lamLine

    if (doCube) then
       call createDataCube(cube, grid, viewVec, nAtom, thisAtom, iAtom, iTrans, nSource, source, nFreqArray, freqArray, &
            occultingDisc)
       
#ifdef MPI
       write(*,*) "Process ",my_rank, " create data cube done"
#endif
       if (myrankiszero) then

#ifdef USECFITSIO
          call writeDataCube(cube,datacubefilename)
#endif
!          write(plotfile,'(a,i3.3,a)') "flatimage",nfile,".fits.gz"
!          call writeCollapsedDataCube(cube,plotfile)
       endif
       call torus_mpi_barrier
       call freeDataCube(cube)
    endif
    if (.not.doSpec) goto 666
#ifdef MPI
  call randomNumberGenerator(syncIseed=.true.)
#endif
  allocate(da(1:maxray), domega(1:maxRay),  rayPosition(1:maxray))
    da = 0.d0; dOmega = 0.d0
    nray = 0; rayPosition = VECTOR(0.d0, 0.d0, 0.d0)

    call createRayGrid(nRay, rayPosition, da, dOmega, viewVec, distance, grid)

    iv1 = 1
    iv2 = nv
 

    allocate(spec(1:nv), vArray(1:nv))
    spec = 0.d0
    if (PRESENT(forceLambda)) then
       varray(1) = 0.d0
    else
       if (nv == 1) then
          varray(1) = 0.d0
       else
          do iv = 1, nv
             vArray(iv) = 2100.e5/cspeed * (2.d0*dble(iv-1)/dble(nv-1)-1.d0)
          enddo
       endif
    endif
    
    do iv = iv1, iv2
       write(*,*) iv
       deltaV  = vArray(iv)

       iray1 = 1
       iray2 = nray
#ifdef MPI
    iray1 = (my_rank) * (nray / (np)) + 1
    iray2 = (my_rank+1) * (nray / (np))
    if (my_rank == (np-1)) iv2 = nray
#endif
    totalOmega = 0.d0
       do iRay = iray1, iray2
          if (PRESENT(occultingDisc)) then
             i0 = intensityAlongRay(rayposition(iRay), viewvec, grid, thisAtom, nAtom, iAtom, iTrans, -deltaV, source, nSource, &
                  nFreqArray, freqArray, occultingDisc=.true., forceFreq=broadBandFreq)
          else
             i0 = intensityAlongRay(rayposition(iRay), viewvec, grid, thisAtom, nAtom, iAtom, iTrans, -deltaV, source, nSource, &
                  nFreqArray, freqArray)
          endif
          spec(iv) = spec(iv) + i0 * domega(iRay) 
          if (i0 > 1.d-20) then
             totalOmega = totalOmega + domega(iray)
          endif
       enddo
    enddo
#ifdef MPI
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     allocate(tempArray(1:nv))
     call MPI_ALLREDUCE(spec,tempArray,nv,MPI_DOUBLE_PRECISION,&
          MPI_SUM,MPI_COMM_WORLD,ierr)
     spec(1:nv) = tempArray(1:nv)
     deallocate(tempArray)

     allocate(tempArray(1:1))
     call MPI_ALLREDUCE(totalOmega,tempArray,1,MPI_DOUBLE_PRECISION,&
          MPI_SUM,MPI_COMM_WORLD,ierr)
     totalOmega=temparray(1)
     deallocate(tempArray)
     

#endif
     if (myrankiszero) then
        write(*,*) "Solid angle check: ", totalOmega, pi* grid%rCore**2/distance**2
     endif

    if (myRankIsZero) then
       write(plotfile,'(a,i3.3,a)') "spec",nfile,".dat"
       open(42, file=plotfile,status="unknown",form="formatted")
       do i = 1, nv
          transitionFreq = thisAtom(iAtom)%transFreq(iTrans)
          write(42, *) vArray(i)*cspeed/1.d5, toPerAngstrom(spec(i), transitionFreq)
       enddo
       close(42)
    endif
    totalFlux = toPerAngstrom(spec(1), broadBandFreq)
    deallocate(vArray, spec)
    deallocate(da, domega, rayPosition)
666 continue
  end subroutine calculateAtomSpectrum


  subroutine createRayGrid(nRay, rayPosition, da, dOmega, viewVec, distance, grid)
    type(GRIDTYPE) :: grid
    integer :: nRay
    type(VECTOR) :: rayPosition(:), viewVec, xProj,yProj
    real(double) :: da(:), dOmega(:), distance
    real(double), allocatable :: rGrid(:), dr(:), phigrid(:), dphi(:)
    real(double) :: rMax, rMin
    integer :: nr, nphi, ir, iphi
    real(double) :: r1 , r2, phi1, phi2, phiOffset
    real(double) :: xPos, yPos, zPos
    integer :: nr1, nr2, i

    nr1 = 20
    nr2 = 20
    nr = nr1 + nr2
    nphi = 20
    nray = 0
    i = 0

    allocate(rGrid(1:nr), dr(1:nr), phiGrid(1:nPhi), dphi(1:nPhi))
    rmin = grid%rCore
!    rMax = 2.d0*grid%octreeRoot%subcellSize
    rMax = 5.d0*grid%rCore

    do ir = 1, nr1
       r1 = rMin * dble(ir-1)/dble(nr1)
       r2 = rMin * dble(ir)/dble(nr1)
       i = i + 1
       rgrid(i) = 0.5d0 * (r1 + r2)
       dr(i) = r2 - r1
!       write(*,*) 1.d10*rGrid(i)/(20.d0*rSol)
    enddo
    

    do ir = 1, nr2
       r1 = rMin + (rmax-rMin) * (dble(ir-1)/dble(nr2))**3
       r2 = rMin + (rmax-rMin) * (dble(ir)/dble(nr2))**3
       i = i + 1
       rgrid(i) = 0.5d0 * (r1 + r2)
       dr(i) = r2 - r1
!       write(*,*) 1.d10*rGrid(i)/(20.d0*rSol)
    enddo
    do iphi = 1, nPhi
       phi1 = twoPi * dble(iphi-1)/dble(nPhi)
       phi2 = twoPi * dble(iphi)/dble(nPhi)
       phiGrid(iPhi) = 0.5d0 * (phi1 + phi2)
       dphi(iPhi) = phi2 - phi1
    enddo

    do ir = 1, nr
       r1 = rGrid(ir)

       call randomNumberGenerator(getDouble=phiOffset)
       phiOffset = phiOffset * dphi(1)
       do iPhi = 1, nPhi
          phi1 = phiGrid(iPhi) + phiOffset
          if (phi1 > twoPi) phi1 = phi1 - twoPi

          xPos = r1 * sin(phi1)
          yPos = 0.d0
          zPos = r1 * cos(phi1)

          xProj =  VECTOR(0.d0, 0.d0, 1.d0)  .cross. viewVec
          call normalize(xProj)
          yProj = viewVec .cross. xProj
         call normalize(yProj)

          nRay = nRay + 1
         rayPosition(nray) =  (xPos * xProj) + (zPos * yProj)
         rayposition(nray) = rayPosition(nRay) + ((-1.d0*grid%octreeRoot%subcellSize*10.d0) * viewVec)

          da(nRay) = pi*( (r1 + dr(ir)/2.d0)**2 - (r1 - dr(ir)/2.d0)**2) * dphi(iPhi)/twoPi
          dOmega(nRay) = da(nRay) / distance**2
       enddo
    enddo
  end subroutine createRayGrid


  subroutine createDataCube(cube, grid, viewVec, nAtom, thisAtom, iAtom, iTrans, nSource, source, &
       nFreqArray, freqArray, occultingDisc)
    use mpi_global_mod
    use input_variables, only : npixels, nv, imageSide, maxVel, &
         positionAngle
    use datacube_mod, only: DATACUBE, initCube, addspatialaxes, addvelocityAxis
    use amr_mod, only : countVoxels
#ifdef MPI
    include 'mpif.h'
#endif
    logical, optional :: occultingDisc
    integer :: nSource
    type(SOURCETYPE) :: source(:)
    type(MODELATOM) :: thisAtom(:)
    real(double) :: freqArray(:)
    integer :: nFreqArray
    integer :: iAtom
    integer :: nAtom
    type(GRIDTYPE) :: grid
    type(DATACUBE) :: cube
    type(VECTOR) :: viewvec, rayPos, xProj, yProj, northVec
    real(double) :: deltaV
    integer :: iTrans
    integer :: ix, iy, iv
    real(double) :: r, xval, yval, vstart,vend
    real(double), allocatable :: vArray(:)
    integer ::  i
    integer :: nMonte, imonte
    integer :: iv1, iv2, nx, ny
    integer :: nPoints, nVoxels, nOctals
    real(double), allocatable :: xPoints(:), yPoints(:)
    real(double) :: dx, dy
    integer :: nRay
    integer, parameter :: maxRay = 10000
    real(double) :: xRay(maxray), yRay(maxray)
    real :: area(maxray)
    real(double) :: totArea
    integer :: iRay

    ! For MPI implementations
    integer       ::   my_rank        ! my processor rank
#ifdef MPI
    integer :: j
    integer       ::   np             ! The number of processes
    integer       ::   ierr           ! error flag
    integer       ::   n
    integer       ::   tag, tag2, tag3
    integer       :: iThread, status(MPI_STATUS_SIZE)

    real(double), allocatable :: tempArray(:)

    ! FOR MPI IMPLEMENTATION=======================================================
    !  Get my process rank # 
    status = 0
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

    ! Find the total # of precessor being used in this run
    call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
    tag = 77
    tag2 = 88
    tag3 = 99
#else
! Set rank to zero for non-MPI cases
    my_rank = 0
#endif

    nMonte = 1

    vStart = -maxVel
    vEnd = maxVel
    nx = npixels
    ny = npixels



    iv1 = 1
    iv2 = nv


#ifdef MPI
    iv1 = int(real(my_rank) * (real(nv) / real(np))) + 1
    iv2 = int(real(my_rank+1) * (real(nv) / real(np)))
    if (my_rank == (np-1)) iv2 = nv
#endif

    if (myRankGlobal == 0) then
       call initCube(cube, nx, ny, nv)
    else
       call initCube(cube, nx, ny, iv2-iv1+1)
    endif
    allocate(vArray(1:nv))
    if (nv > 1) then
       do i = 1, nv
          vArray(i) = vStart + (vEnd-vStart)*dble(i-1)/dble(nv-1)
       enddo
    else
       vArray(1) = 0.d0
    endif


    call addSpatialAxes(cube, -dble(imageSide/2.), dble(imageSide/2.), &
         -dble(imageSide/2.), dble(imageSide/2.))


    if (myRankGlobal == 0) then
       call addVelocityAxis(cube, vStart, vEnd)
    else
       call addvelocityAxis(cube, vArray(iv1), vArray(iv2))
    endif

    northVec = VECTOR(0.d0, 0.d0, 1.d0)
    northVec = rotateY(northVec, dble(positionAngle(1)))

    xProj =   viewVec .cross. northVec
    call normalize(xProj)
    yProj =  xProj .cross.viewVec
    call normalize(yProj)

    call countVoxels(grid%octreeRoot, nOctals, nVoxels)
    nPoints = nVoxels + 1000 * nSource + cube%nx*cube%ny
    allocate(xPoints(1:nPoints),yPoints(1:nPoints))
    call createRayGridGeneric(grid, nSource, source, viewVec, xProj, yProj, xPoints, yPoints, nPoints)
    do ix = 1, cube%nx
       do iy = 1, cube%ny
          nPoints = nPoints + 1
          xPoints(nPoints) = cube%xAxis(ix)
          yPoints(nPoints) = cube%yAxis(iy)
       enddo
    enddo

!    do i = 1, nPoints
!       write(76,*) xpoints(i),ypoints(i)
!    enddo
    
    dx = cube%xAxis(2)-cube%xAxis(1)
    dy = cube%yAxis(2)-cube%yAxis(1)

    do iv = iv1, iv2
       deltaV = cube%vAxis(iv-iv1+1)*1.d5/cSpeed
       !$OMP PARALLEL DEFAULT (NONE) &
       !$OMP PRIVATE (ix, iy, iMonte, r, xval, yval, rayPos, nRay, xRay, yRay, area,totArea) &
       !$OMP SHARED (cube, viewVec, grid, thisAtom, nAtom, iAtom, iTrans) &
       !$OMP SHARED (deltaV, source, nSource, nFreqArray, freqArray, occultingDisc) &
       !$OMP SHARED (iv, iv1, xproj, yproj, nMonte, dx, dy, xPoints, yPoints, nPoints)

       !$OMP DO SCHEDULE(DYNAMIC,2)
       do ix = 1, cube%nx
          do iy = 1, cube%ny
             call findRaysInPixel(cube%xAxis(ix),cube%yAxis(iy),dx,dy,xPoints, yPoints, &
                 nPoints,  nRay, xRay, yRay, area)
             
             totArea = 0.d0
             cube%intensity(ix,iy,iv-iv1+1) = 0.d0
             do iRay = 1, nRay
                
                rayPos =  (xRay(iRay) * xProj) + (yRay(iRay) * yProj)
                raypos = rayPos + ((-1.d0*grid%octreeRoot%subcellsize*30.d0) * Viewvec)
                
                cube%intensity(ix,iy,iv-iv1+1) = cube%intensity(ix,iy,iv-iv1+1) &
                     + intensityAlongRayGeneric(rayPos, viewVec, grid,  &
                     -deltaV, source, nSource) * area(iRay)
                totArea = totArea + Area(iray)
             enddo
!             write(*,*) "Pixel done with ",nRay, " rays. check on area ",totArea/(dx**2)
             cube%intensity(ix,iy,iv-iv1+1) = cube%intensity(ix,iy,iv-iv1+1) / SUM(area(1:nRay))
             
          enddo
       enddo
       !$OMP END DO
       !$OMP BARRIER
       !$OMP END PARALLEL
       write(*,*) "Velocity bin ",iv, " done."

    enddo

    deallocate(xPoints, yPoints)


#ifdef MPI
    write(*,*) "Process ",my_rank, " done. awaiting reduce"
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 


    n = (cube%nx * cube%ny)

    do iThread = 1, nThreadsGlobal-1
       if (my_rank == iThread) then
          call MPI_SEND(iv2-iv1+1, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD,  ierr)
          do iv = 1, nv                
             if ((iv >= iv1).and.(iv <= iv2)) then
                allocate(tempArray(1:n))
                tempArray = reshape(cube%intensity(:,:,iv-iv1+1), (/  n /))
                call MPI_SEND(iv, 1, MPI_INTEGER, 0, tag2, MPI_COMM_WORLD,  ierr)
                call MPI_SEND(tempArray, n, MPI_DOUBLE_PRECISION, 0, tag3, MPI_COMM_WORLD,  ierr)
                deallocate(tempArray)
             endif
          enddo
       endif
       if (my_rank == 0) then
          call MPI_RECV(j, 1, MPI_INTEGER, iThread, tag, MPI_COMM_WORLD, status, ierr)
          do i = 1, j
             call MPI_RECV(iv, 1, MPI_INTEGER, iThread, tag2, MPI_COMM_WORLD, status, ierr)
             allocate(tempArray(1:n))
             tempArray = reshape(cube%intensity(:,:,iv-iv1+1), (/  n /))
             call MPI_RECV(tempArray, n, MPI_DOUBLE_PRECISION, iThread, tag3, MPI_COMM_WORLD, status, ierr)
             cube%intensity(:, :, iv) = reshape(tempArray, (/ cube%nx, cube%ny /))
             deallocate(tempArray)
          enddo
       endif
       call torus_mpi_barrier
    enddo



    write(*,*) "Process ",my_rank, " reduce done."
#endif

!    cube%flux = cube%intensity

  end subroutine createDataCube


  subroutine findRaysInPixel(xcen, yCen, dx, dy, xPoints, yPoints, nPoints, &
                  nRay, xRay, yRay, area)
    use utils_mod, only : voron2

    real(double) :: xCen, yCen, dx, dy
    real(double) :: xPoints(:), yPoints(:)
    integer :: nPoints
    integer :: nRay, i
    real(double) :: xRay(:), yRay(:)
    real :: area(:)
    real, allocatable :: xTmp(:), yTmp(:)
    nRay = 0
    do i = 1, nPoints
       if ((abs(xCen-xPoints(i)) < dx/2.d0).and. &
            (abs(yCen-yPoints(i)) < dy/2.d0)) then
          nRay = nRay + 1
          xRay(nRay) = xPoints(i)
          yRay(nRay) = yPoints(i)
          if (nRay == SIZE(xRay)) then
             write(*,*) "max array sized reached for nray"
             exit
          endif
       endif
    enddo
    if (nRay == 1) then
       area(1) = dx*dy
    else if (nRay == 2) then
       area(1:2) = 0.5*dx*dy
    else if (nRay == 3) then
       area(1:2) = 0.33333*dx*dy
    else
       allocate(xtmp(1:nray),ytmp(1:nRay))
       xtmp = xRay - (xcen-dx/2.d0)
       ytmp = yRay - (yCen-dy/2.d0)
       call voron2(nRay, xTmp, yTmp, real(dx), area(1:nRay))
       deallocate(xTmp, yTmp)
    endif
  end subroutine findRaysInPixel
          
          


#ifdef MPI
      subroutine packAtomLevel(octalArray, nTemps, tArray, iAtom, iLevel, ioctal_beg, ioctal_end)
    include 'mpif.h'
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: nTemps
        integer :: ioctal_beg, ioctal_end
        real(double) :: tArray(:)
        integer :: iOctal, iSubcell, my_rank, ierr, iAtom
        integer :: iLevel
        type(OCTAL), pointer :: thisOctal

       call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
       !
       ! Update the edens values of grid computed by all processors.
       !
       nTemps = 0
       do iOctal = 1, SIZE(octalArray)

          thisOctal => octalArray(iOctal)%content
          
          do iSubcell = 1, thisOctal%maxChildren
              if (.not.thisOctal%hasChild(iSubcell)) then
                 nTemps = nTemps + 1
                 if ((ioctal >= ioctal_beg).and.(iOctal<=ioctal_end)) then
                   tArray(nTemps) = thisOctal%newAtomLevel(isubcell, iAtom, iLevel)
                 else 
                   tArray(nTemps) = 0.d0
                 endif
              endif
          end do
       end do
     end subroutine packAtomLevel

      subroutine unpackAtomLevel(octalArray, nTemps, tArray, iAtom, iLevel)
    include 'mpif.h'
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: nTemps
        real(double) :: tArray(:)
        integer :: iOctal, iSubcell, my_rank, ierr, iAtom
        integer :: iLevel
        type(OCTAL), pointer :: thisOctal

       call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
       !
       ! Update the edens values of grid computed by all processors.
       !
       nTemps = 0
       do iOctal = 1, SIZE(octalArray)

          thisOctal => octalArray(iOctal)%content
          
          do iSubcell = 1, thisOctal%maxChildren
              if (.not.thisOctal%hasChild(iSubcell)) then
                 nTemps = nTemps + 1
                 thisOctal%newAtomLevel(isubcell, iAtom, iLevel) = tArray(nTemps) 
              endif
          end do
       end do
     end subroutine unpackAtomLevel

      subroutine packjnu(octalArray, nTemps, tArray, iFreq, ioctal_beg, ioctal_end)
    include 'mpif.h'
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: ioctal_beg, ioctal_end
        integer :: nTemps
        real(double) :: tArray(:)
        integer :: iOctal, iSubcell, my_rank, ierr, iFreq
        type(OCTAL), pointer :: thisOctal

       call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
       !
       ! Update the edens values of grid computed by all processors.
       !
       nTemps = 0
       do iOctal = 1, SIZE(octalArray)

          thisOctal => octalArray(iOctal)%content
          
          do iSubcell = 1, thisOctal%maxChildren
              if (.not.thisOctal%hasChild(iSubcell)) then
                 nTemps = nTemps + 1
                 if ((ioctal >= ioctal_beg).and.(iOctal<=ioctal_end)) then
                   tArray(nTemps) = thisOctal%jnuCont(isubcell, ifreq)
                 else 
                   tArray(nTemps) = 0.d0
                 endif
              endif
          end do
       end do
     end subroutine packJnu

      subroutine unpackJnu(octalArray, nTemps, tArray, iFreq)
    include 'mpif.h'
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: nTemps
        real(double) :: tArray(:)
        integer :: iOctal, iSubcell, my_rank, ierr, iFreq
        type(OCTAL), pointer :: thisOctal

       call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
       !
       ! Update the edens values of grid computed by all processors.
       !
       nTemps = 0
       do iOctal = 1, SIZE(octalArray)

          thisOctal => octalArray(iOctal)%content
          
          do iSubcell = 1, thisOctal%maxChildren
              if (.not.thisOctal%hasChild(iSubcell)) then
                 nTemps = nTemps + 1
                 thisOctal%jnuCont(isubcell, iFreq) = tArray(nTemps) 
              endif
          end do
       end do
     end subroutine unpackJnu
#endif

  recursive  subroutine  setMicroturb(thisOctal, microTurb)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: microTurb
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setMicroturb(child, microTurb)
                exit
             end if
          end do
       else
          thisOctal%microTurb(subcell) = microturb
       endif
    enddo
  end subroutine setMicroturb


  subroutine testRays(grid, nSource, source)
    type(GRIDTYPE) :: grid
    integer :: nSource
    type(SOURCETYPE) :: source(:)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: weight
    type(VECTOR) :: position, direction, pvec, photoDirection
    integer :: sourceNumber, i, j, nray, iElement
    real(double) :: i0, distTosource, cosTheta, freq
    logical :: hitSource
    
    freq  = cSpeed / (6562.8d0 * 1.d-8)
    position = VECTOR(0.d0 , 0.d0 , source(1)%radius*5.d0)

    thisOctal => grid%octreeRoot
    call findSubcellLocal(position, thisOctal, subcell)

    do i = 1, 10
       nray = 100 * 2**(i-1)
       i0 = 0.d0
       do j = 1, nray
          call randomRayDirection(0.8d0, position, source, nSource, direction, weight)
          
          call distanceToSource(source, nSource, position, direction, hitSource, disttoSource, sourcenumber)
          if (hitSource) then
             pVec = (position + (direction * distToSource) - source(sourceNumber)%position)
             call normalize(pVec)
             cosTheta = -1.d0*(pVec.dot.direction)
             photoDirection = pVec
             call normalize(photoDirection)
             iElement = getElement(source(sourcenumber)%surface, photoDirection)
             i0 = i0 + i_nu(source(sourceNumber), freq, iElement)*weight
          endif
       enddo
       write(*,*) "nray ",nray, "j_nu ",i0/dble(nray)
    enddo
  end subroutine testRays


  subroutine createRayGridGeneric(grid, nSource, SourceArray, viewVec, xProj, yProj, xPoints, yPoints, nPoints)
    type(GRIDTYPE) :: grid
    integer :: nSource
    type(SOURCETYPE) :: sourceArray(:)
    type(VECTOR) :: viewVec, xProj, yProj
    integer :: nPoints, i, j
    real(double) :: xPoints(:), yPoints(:)
    integer :: nr, nphi
    real(double) :: r, phi, rStar, rMin, rMax, dphi
    nr = 10
    nphi = 20

    nPoints = 0
    rStar = 1.1*sourceArray(1)%radius
    rMin = rStar/20.
    do i = 1, nr
       r = rMin + (rStar-rMin)*real(i-1)/real(nr-1)
       call randomNumberGenerator(getDouble=dphi)
       dphi = dphi * twoPi
       do j = 1, nPhi
          phi = dphi + twoPi * real(j-1)/real(nPhi)

          nPoints = nPoints + 1
          xPoints(nPoints) = r * cos(phi)
          yPoints(nPoints) = r * sin(phi)
       enddo
    enddo

    nr = 100
    nphi = 50
    rMin = 1.2 * rStar
    rMax = grid%octreeRoot%subcellSize*2.d0

    do i = 1, nr
       r = log10(rMin) + log10(rMax/rmin)*real(i-1)/real(nr-1)
       r = 10.d0**r
       call randomNumberGenerator(getDouble=dphi)
       dphi = dphi * twoPi
       do j = 1, nPhi
          phi = dphi + twoPi * real(j-1)/real(nPhi)

          nPoints = nPoints + 1
          xPoints(nPoints) = r * cos(phi)
          yPoints(nPoints) = r * sin(phi)
       enddo
    enddo
    
  end subroutine createRayGridGeneric

  recursive  subroutine  getProjectedPoints(thisOctal, viewVec, xProj, yProj, xPoints, yPoints, nPoints)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    integer :: nPoints
    type(VECTOR) :: xProj, yProj, viewVec, rVec
    real(double) :: xPoints(:), yPoints(:)
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call getProjectedPoints(child, viewVec, xProj, yProj, xPoints, yPoints, nPoints)
                exit
             end if
          end do
       else
          nPoints = nPoints + 1
          rVec = subcellCentre(thisOctal, subcell)
          xPoints(nPoints) = xProj.dot.rVec
          yPoints(nPoints) = yProj.dot.rVec
       endif
    enddo
  end subroutine getProjectedPoints

  subroutine getSurfacePoints(source, viewVec, xProj, yProj, xPoints, yPoints, nPoints)
    type(SOURCETYPE) :: source
    type(VECTOR) :: viewVec, xProj, yProj
    real(double) :: xPoints(:), yPoints(:)
    integer :: nPoints, i

  end subroutine getSurfacePoints

end module cmf_mod
