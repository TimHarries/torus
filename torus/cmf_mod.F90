module cmf_mod

  ! written by tjh


  use kind_mod
  use constants_mod
  use messages_mod
  use vector_mod
  use octal_mod, only: octal, octalWrapper, subcellCentre
  use amr_mod, only: inOctal, distanceToCellBoundary, findsubcellLocal, amrGridVelocity
  use gridtype_mod, only: GRIDTYPE
  use utils_mod, only: myexp, init_random_seed, gaussj, locate
  use modelatom_mod, only: MODELATOM, BoltzSahaGeneral, bfOpacity, bfEmissivity, photoCrossSection, collisionRate, &
        addcrosssectionstoatom, createcontfreqarray, createrbbarrays, returneinsteincoeffs
  use source_mod, only: SOURCETYPE, I_nu, insideSource, distanceToSource
  use surface_mod, only: getElement
  use timing, only: tune
  use parallel_mod, only: torus_mpi_barrier
  use vtk_mod, only: writeVtkfile


  implicit none

  private
  public:: atomLoop, calculateAtomSpectrum

contains

  subroutine solveLevels(thisOctal, subcell, nPops, jnuLine,  &
       temperature, nAtom, thisAtom, ne, rho, jnuCont, freq, dfreq, nfreq)
!    use input_variables, only : debug
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

          NstarRatio = boltzSahaGeneral(thisAtom(iAtom), 1, l, ne, temperature)



          if ((thisAtom(iAtom)%transType(iTrans) == "RBF").and.(.not.thisAtom(iAtom)%inDetailedBalance(iTrans)))  then ! radiative recomb

             recombratekl = 0.d0
             nuStart = (thisAtom(iAtom)%iPot - thisAtom(iAtom)%energy(l))*evtoerg/hcgs


             call locate(freq, nfreq, nuStart, istart)

             !             istart = 1
             do i = istart+1, nFreq
                xSection = photoCrossSection(thisAtom(iAtom), iTrans, l, freq(i))
                recombRatekl = recombRatekl + &
                     (fourPi/(hCgs*freq(i)))*xSection*((2.d0*hCgs*freq(i)**3)/cSpeed**2 + jnuCont(i)) * &
                     exp(-(hCgs*freq(i))/(kErg*temperature))*(freq(i)-freq(i-1))
             enddo


             !             thisFreq = (thisAtom(iatom)%iPot-thisAtom(iAtom)%energy(l))*evtoerg/hcgs
             !             write(*,*) "thisFreq",thisFreq,thisFreq*hcgs*ergtoev
             !             xSection = photoCrossSection(thisAtom(iAtom), l, thisFreq*1.01d0)
             !             test = (8.d0*pi/cspeed**2) * xSection *thisFreq**3*expint(1,hcgs*thisFreq/(kerg*temperature))
             !             write(*,*) "test",l,test,recombratekl



             matrixA(l+nOffset(iAtom),k+nOffset(iAtom)) = matrixA(l+nOffset(iAtom),k+nOffset(iAtom)) + NstarRatio * recombratekl
             tot1 = tot1  + NstarRatio * recombratekl


             radtot = radtot + tot1 * ne 

          endif


          if (thisAtom(iAtom)%transType(iTrans) == "CBF") then ! collisional recomb

             collEx = collisionRate(thisAtom(iAtom), iTrans, temperature) 

             matrixA(l+nOffset(iAtom),k+nOffset(iAtom)) = matrixA(l+nOffset(iAtom),k+nOffset(iAtom)) + NstarRatio*collEx*ne
             tot2 = tot2 + NstarRatio*collEx*ne


             colltot = colltot + tot2 * ne
          
          endif


          totRecomb =  tot1 + tot2
          if (thisAtom(iAtom)%transType(iTrans)(2:3) == "BF") then
             matrixA(k+nOffset(iAtom), k+nOffset(iAtom)) = matrixA(k+nOffset(iAtom), k+nOffset(iAtom))-totRecomb
          endif


       enddo

!      if (iatom == 1) write(*,*) "total recomb: ", totrecomb, radtot,colltot


       totPhotoIon = 0.d0
       totcion = 0.d0
       do iTrans = 1, thisAtom(iAtom)%nTrans

          k = thisAtom(iAtom)%iUpper(iTrans)
          l = thisAtom(iAtom)%iLower(iTrans)



          if (thisAtom(iAtom)%transType(iTrans) == "RBB") then ! radiative BB rates
             iJnu = thisAtom(iAtom)%indexRBBTrans(iTrans)
             call returnEinsteinCoeffs(thisAtom(iAtom), iTrans, a, Bul, Blu)



             call locate(freq, nfreq, thisAtom(iatom)%transFreq(itrans), istart)

             jnu = jNuLine(iJnu)


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
                   xSection = photoCrossSection(thisAtom(iAtom), iTrans, l, freq(i-1))
                   xSection2 = photoCrossSection(thisAtom(iAtom), iTrans, l, freq(i))
                   photoRatelk = photoRatelk + 0.5d0 * ((jnuCont(i-1)/(hCgsfreq(i-1)))*xSection &
                        + (jnuCont(i)/(hCgsfreq(i)))*xSection2) * dfreq(i)
!                   if (iAtom == 1) write(*,*) i,l,freq(i),photoratelk,xsection,jnucont(i)
                enddo
                photoRatelk = photoRatelk * fourPi
                

                matrixA(l+nOffset(iAtom),l+nOffset(iAtom)) = matrixA(l+nOffset(iAtom),l+nOffset(iAtom)) - photoRatelk
                matrixA(k+nOffset(iAtom),l+nOffset(iAtom)) = matrixA(k+nOffset(iAtom),l+nOffset(iAtom)) + photoRatelk
                totphotoion=totphotoion + photoRatelk*npops(iatom,l)

             endif



!             if ((l == 1).and.iatom==1) then
!                write(*,*) "total photoionization from ground state ",photoRatelk*npops(iatom,l)
!             endif
!             if ((l == 2).and.iatom==1) then
!                write(*,*) "total photoionization from l=2 ",photoRatelk*npops(iatom,l)
!             endif
!             if ((l == 3).and.iatom==1) then
!                write(*,*) "total photoionization from l=3 ",photoRatelk*npops(iatom,l)
!             endif
!             if ((l == 4).and.iatom==1) then
!                write(*,*) "total photoionization from l=4 ",photoRatelk*npops(iatom,l)
!             endif

          endif


          if (thisAtom(iAtom)%transType(iTrans) == "CBF") then ! collisional ionization

             k = thisAtom(iAtom)%iUpper(iTrans)
             l = thisAtom(iAtom)%iLower(iTrans)


             collEx = collisionRate(thisAtom(iAtom), iTrans, temperature)



             matrixA(l+nOffset(iAtom),l+nOffset(iAtom)) = matrixA(l+nOffset(iAtom),l+nOffset(iAtom)) - collex * ne
             matrixA(k+nOffset(iAtom),l+nOffset(iAtom)) = matrixA(k+nOffset(iAtom),l+nOffset(iAtom)) + collex * ne
             totcion = totcion + collex*ne*npops(iatom,l)
          endif


       enddo

!       if (iatom == 1) write(*,*) "total photoionization: " , totphotoion, totcion

    enddo

!
!    if (debug) then
!       write(*,'(4x,100i8)') num(1:nMatrix)
!       do i = 1, nMatrix
!          write(*,'(i4,1p,100e9.1)') i,matrixA(i,1:nMatrix),matrixB(i)
!       enddo
!    endif

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
       nPops(iAtom,1:thisAtom(iAtom)%nLevels) = matrixB(1+nOffset(iAtom):thisAtom(1)%nLevels+nOffset(iAtom))
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


  subroutine getRay(grid, fromOctal, fromSubcell, position, direction, ds, phi, i0, &
       hCol, HeIcol, HeIIcol, nAtom, thisAtom, source, nSource, &
       hitPhotosphere, sourceNumber, cosTheta, weight, nRBBTrans, indexRBBTrans, indexAtom, nHAtom, nHeIAtom, nHeIIAtom, &
       nFreq, freq, iCont)
    use input_variables, only : opticallyThickContinuum
    use amr_mod, only: randomPositionInCell
    type(SOURCETYPE) :: source(:)
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
    integer :: nHatom, nHeIAtom, nHeIIAtom
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
    real(double) :: weight
    integer :: iAtom, iRBB
    real(double) :: phiAv, phiNorm
    logical :: firstSubcell
    integer :: nBug
    integer :: iFreq
    integer :: iElement
    integer :: j
    real(double), allocatable :: tauCont(:), jnuCont(:), alphanuCont(:), snuCont(:)
    logical, save :: first = .true.
    real(double) :: nStar(5,40)
    logical :: passThroughResonance, velocityIncreasing
    real(double) :: x1, x2, fac, deltaDist

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

    call randomRayDirection(0.9d0, position, source, nSource, direction, weight)



    call random_number(r)

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

    !   rayVel = getVel(grid, thisOctal, subcell, position, direction)


    call random_number(r)

    !    deltaV = 4.3 * thisOctal%microturb(subcell) * (r-0.5d0)

    ! assume microturb is gaussian
    !    deltaV = gasdev()*thisOctal%microturb(subcell)


    deltaV = 4.3 * thisOctal%microturb(subcell) * (r-0.5d0) ! random frequency near line spectrum peak. 
    ! 4.3 corresponds to the width where the peak of the line profile has dropped to 1% of its peak
    ! microturulence is assumed gaussian - b is FULL WIDTH

    weight = weight * phiProf(deltaV, thisOctal%microturb(subcell))

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
          nTau = 80
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

          if (opticallyThickContinuum.and.(i==2)) then

             do iAtom = 1, nAtom
                do j = 1, thisAtom(iAtom)%nLevels - 1
                   nStar(iAtom,j) = BoltzSahaGeneral(thisAtom(iAtom), 1, j, thisOctal%ne(subcell), &
                        dble(thisOctal%temperature(subcell))) * &
                        Thisoctal%atomlevel(subcell, iAtom,thisAtom(iatom)%nLevels)
                enddo
             enddo


             do iFreq = 1, nFreq

                alphanuCont(ifreq) = bfOpacity(freq(ifreq), nAtom, thisAtom, thisOctal%atomLevel(subcell,:,:), &
                     thisOctal%ne(subcell), nstar, dble(thisOctal%temperature(subcell)), ifreq=ifreq)

                jnuCont(iFreq) = bfEmissivity(freq(ifreq), nAtom, thisAtom, &
                     nStar, dble(thisOctal%temperature(subcell)), thisOctal%ne(subcell), &
                     thisOctal%jnuCont(subcell,ifreq), ifreq=ifreq)/fourpi


                if (alphanuCont(ifreq) /= 0.d0) then
                   snuCont(iFreq) = jnuCont(iFreq)/alphanuCont(iFreq)
                else
                   snuCont(iFreq) = tiny(snuCont(iFreq))
                endif
             enddo
          endif

          if (opticallyThickContinuum) then
             deltaDist = (distArray(i)-distArray(i-1)) * 1.d10
             do iFreq = 1, nFreq
                dTau = alphaNuCont(iFreq) *  deltaDist
                iCont(iFreq) = iCont(ifreq) +  myexp(-tauCont(iFreq)) * (1.d0-myexp(-dtau))*snuCont(iFreq)
                tauCont(iFreq) = tauCont(iFreq) + dtau
             enddo
          endif
          !          write(*,*) tauCont(48),iCont(48),snuCont(48)

          dv = deltaV - (thisVel .dot. direction)

          icount = icount + 1

          do iRBB = 1, nRBBTRans
             iAtom = indexAtom(iRBB)
             iTrans = indexRBBTrans(iRBB)


             call locate(freq, nfreq, thisAtom(iAtom)%transFreq(iTrans), iFreq)

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

             !                if (isnan(alphanu)) then
             !                   write(*,*) "alphanu isnan",nlower,nupper
             !                   stop
             !                endif

             if (opticallyThickContinuum) alphanu = alphanu + alphaNuCont(iFreq)


             if (alphanu < 0.d0) then
                alphanu = 0.d0
                if (first) then
                   write(*,*) "negative opacity warning in getray",iUpper,iLower,nLower,nUpper,thisAtom%name
                   first = .false.
                endif
             endif



             dTau = alphaNu *  (distArray(i)-distArray(i-1)) * 1.d10


             etaLine = hCgs * a * thisAtom(iAtom)%transFreq(iTrans)
             etaLine = etaLine * thisOctal%atomLevel(subcell, iAtom, iUpper)
             jnu = (etaLine/fourPi) * phiProf(dv, thisOctal%microturb(subcell))/thisAtom(iAtom)%transFreq(iTrans)

             if (opticallyThickContinuum) jnu = jnu + jnuCont(iFreq) 


             if (alphanu > 1.d-30) then
                snu = jnu/alphanu
             else
                snu = tiny(snu)
             endif


             i0(iRBB) = i0(iRBB) +  exp(-tau(irBB)) * (1.d0-exp(-dtau))*snu
             !                if (dtau > 1.d10) write(*,*) dtau

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
    if (hitPhotosphere) then ! don't include weight below - that's done when jbar is calculated
       iElement = getElement(source(sourcenumber)%surface, photoDirection)
       do iFreq = 1, nFreq
          !          write(*,*) freq(ifreq),iCont(ifreq),i_nu(source(sourceNumber), &
          !               freq(iFreq), iElement)*cosTheta*exp(-tauCont(iFreq)),taucont(ifreq)
          !          write(*,*) freq(ifreq),i_nu(source(sourcenumber),freq(ifreq),ielement)
          !          iCont(iFreq) = iCont(iFreq) + i_nu(source(sourceNumber), freq(iFreq), iElement)*cosTheta*exp(-tauCont(iFreq))
          iCont(iFreq) = iCont(iFreq) + i_nu(source(sourceNumber), freq(iFreq), iElement)*exp(-tauCont(iFreq))
       enddo
    endif

    !    call locate(freq, nfreq, cspeed/6562.8d-8,iFreq)
    !    write(*,*) "icont ",icont(ifreq)



    if (hitPhotosphere) then ! don't include weight below - that's done when jbar is calculated
       iElement = getElement(source(sourcenumber)%surface, photoDirection)
       do iRBB = 1, nRBBTrans
          iAtom = indexAtom(iRBB)
          iTrans = indexRBBTrans(iRBB)
          !          i0(iRBB) = i0(iRBB) + i_nu(source(sourceNumber), thisAtom(iATom)%transFreq(iTrans), iElement)*cosTheta*exp(-tau(iRBB))
          i0(iRBB) = i0(iRBB) + i_nu(source(sourceNumber), thisAtom(iATom)%transFreq(iTrans), iElement)*exp(-tau(iRBB))
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

  subroutine calculateJbar(thisOctal, subcell, thisAtom, nRay, ds, phi, i0, iTrans, jbar, nPops, &
       freq, nfreq, weight, iRBB, tauAv)
    real(double) :: freq(:), weight(:)
    integer :: nfreq
    integer :: iRBB
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    type(MODELATOM) :: thisAtom
    integer :: nRay
    Real(double) :: ds(:), phi(:), i0(:), nPops(:)
    integer :: iTrans
    real(double) :: jbar
    integer :: iRay
    real(double) :: nLower, nUpper
    real(double) :: jBarInternal, jBarExternal
    real(double) :: alphanu, jnu, etaline
    integer :: iUpper, iLower
    real(double) :: tau, snu, sumWeight
    real(double) :: a, bul, blu
    logical,save :: first = .true.
    real(double) :: tauAv
    jBarExternal = 0.d0
    jBarInternal = 0.d0
    a = 0.d0; bul = 0.d0; blu = 0.d0

    if (thisAtom%transType(iTrans) == "RBB") then

       iUpper = thisAtom%iUpper(iTrans)
       iLower = thisAtom%iLower(iTrans)

       sumWeight = 0.d0
       tauAv = 0.d0
       do iRay = 1, nRay
          alphanu = (hCgs*thisAtom%transFreq(iTrans)/fourPi)
          nLower = nPops(iLower)
          nUpper = nPops(iUpper)

          call returnEinsteinCoeffs(thisAtom, iTrans, a, Bul, Blu)

          alphanu = alphanu * (nLower * Blu - nUpper * Bul) * phi(iray)/thisAtom%transFreq(iTrans)

          if (alphanu < 0.d0) then
             alphanu = 0.d0
             if (first) then
                write(*,*) "negative opacity warning in calcjbar",iUpper,iLower,nLower,nUpper,thisAtom%name
                first = .false.
             endif
          endif

          tau = alphaNu * ds(iray)
          tauAv = tauAv + tau * weight(iray)

          etaLine = hCgs * a * thisAtom%transFreq(iTrans)
          etaLine = etaLine *  nPops(iUpper)
          jnu = (etaLine/fourPi) * phi(iRay)/thisAtom%transFreq(iTrans)

          if (alphanu /= 0.d0) then
             snu = jnu/alphanu
          else
             snu = tiny(snu)
          endif

! need to include continuum here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!          jBarExternal = jBarExternal + i0(iray) * exp(-tau) * phi(iRay) * weight(iRay)
!          jBarInternal = jBarInternal + snu * (1.d0 - exp(-tau)) * phi(iRay) * weight(iRay)
          
! if we have choosen the photon velocity according to a gaussian we 
! don't need to weight jbar by phi


          jBarExternal = jBarExternal + i0(iray) * exp(-tau) * weight(iRay)
          jBarInternal = jBarInternal + snu * (1.d0 - exp(-tau))  * weight(iRay)

 
!          sumWeight = sumWeight + weight(iRay) * phi (iRay)
          sumWeight = sumWeight + weight(iRay) 

       enddo

       jbar = (jBarExternal + jBarInternal)/sumWeight
       tauAv = tauAv / sumWeight

    endif
          
  end subroutine calculateJbar

  subroutine calculateJbarCont(thisOctal, subcell, nAtom, thisAtom, ne, nray, ds, source, nSource, &
       hitPhotosphere, sourceNumber, freq, nfreq, &
       iCont, jBarCont, cosTheta, weight)
    use input_variables, only : opticallyThickContinuum
    real(double) :: iCont(:,:), jBarCont(:), ds(:), ne
    integer :: nAtom
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    type(MODELATOM) :: thisAtom(:)
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    logical :: hitPhotosphere(:)
    integer :: sourceNumber(:)
    integer :: nfreq
    real(double) :: freq(:),  cosTheta(:), weight(:)
    integer :: iray, nRay
    real :: tau
    real(double), allocatable :: jBarContExternal(:), jBarContInternal(:)
    real(double) :: alphanu
    integer :: iFreq
    real(double) :: jnu, snu
    real(double) :: nstar(10,50)
    integer :: iAtom, j
    allocate(jBarContExternal(1:nFreq))
    allocate(jBarContInternal(1:nFreq))

    jBarCont = 0.d0
    jBarContInternal = 0.d0
    jBarContExternal = 0.d0

    do iAtom = 1, nAtom
       do j = 1, thisAtom(iAtom)%nLevels - 1
          nStar(iAtom,j) = BoltzSahaGeneral(thisAtom(iAtom), 1, j, thisOctal%ne(subcell), &
               dble(thisOctal%temperature(subcell))) * &
               Thisoctal%atomlevel(subcell, iAtom,thisAtom(iatom)%nLevels)
       enddo
    enddo


    do iRay = 1, nRay
       do iFreq = 1, nFreq

          tau = 0.d0
          snu = 0.d0
          if (opticallyThickContinuum) then
             alphanu = bfOpacity(freq(ifreq), nAtom, thisAtom, thisOctal%atomLevel(subcell,:,:), ne, &
                  nstar, dble(thisOctal%temperature(subcell)),ifreq=ifreq)
             jnu =  bfEmissivity(freq(ifreq), nAtom, thisAtom, &
                  nstar, &
                  dble(thisOctal%temperature(subcell)), thisOctal%ne(subcell), thisOctal%jnuCont(subcell,ifreq), &
                  ifreq=ifreq)/fourPi
             tau = alphaNu * ds(iray)
             if (alphanu /= 0.d0) then
                snu = jnu/alphanu
             else
                snu = tiny(snu)
             endif
          endif
          jBarContExternal(iFreq) = jBarContExternal(iFreq) + iCont(iray, iFreq) * exp(-tau) * weight(iRay)
!          write(*,*) 1.d8*cspeed/freq(iFreq)," icont(iray,ifreq)",iCont(iRay,iFreq),tau
          jBarContInternal(iFreq) = jBarContInternal(iFreq) + snu * (1.d0 - exp(-tau)) * weight(iRay)
       enddo
    enddo

    jbarCont = (jBarContExternal + jBarContInternal)/SUM(weight(1:nRay))


  end subroutine calculateJbarCont


  subroutine atomLoop(grid, nAtom, thisAtom, nSource, source)

    use input_variables, only : debug, rcore, lte
    use messages_mod, only : myRankIsZero
    use gridio_mod, only : writeAmrGrid
    use mpi_global_mod, only: myRankGlobal
    use amr_mod, only: countVoxels, getOctalArray
#ifdef MPI
    use input_variables, only : blockhandout
    use parallel_mod, only: mpiBlockHandout, mpiGetBlock
    include 'mpif.h'
#endif

    type(SOURCETYPE) :: source(:)
    integer :: nSource
    integer :: nAtom
    type(GRIDTYPE) :: grid
    type(MODELATOM) :: thisAtom(:)
    type(VECTOR) :: position, direction
    integer :: nOctal, iOctal, subcell
    real(double), allocatable :: ds(:), phi(:), i0(:,:)
    real(double), allocatable :: Hcol(:), HeICol(:), HeIICol(:)
    integer :: nRay
    type(OCTAL), pointer :: thisOctal
    integer, parameter :: maxIter = 10000, maxRay = 100000
    logical :: popsConverged, gridConverged 
    integer :: iRay, iTrans, iter,i 
    integer :: iStage
    real(double), allocatable :: oldpops(:,:), newPops(:,:), dPops(:,:), mainoldpops(:,:)
    real(double) :: newNe, dNe
    real(double), parameter :: underCorrect = 0.9d0
    real(double) :: fac
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: ioctal_beg, ioctal_end
    real(double) :: maxFracChange
    logical :: fixedRays
    integer :: isize
    integer, allocatable :: iseed(:)
    real(double) :: tolerance
    integer, allocatable :: sourceNumber(:)
!    type(VECTOR) :: posVec
!    real(double) :: r
    real(double), allocatable :: cosTheta(:)
    real(double), allocatable :: weight(:)
    real(double), allocatable :: iCont(:,:)
    logical, allocatable :: hitPhotosphere(:)
    integer, parameter :: maxFreq = 2000
    real(double) :: freq(maxFreq), dFreq(maxFreq)
    real(double), allocatable :: jnuCont(:)
    integer :: nFreq, nhit, iRBB
    integer :: nRBBTrans
    integer :: indexRBBTrans(1000), indexAtom(1000)
    real(double) :: ne
    integer :: iAtom
    integer :: nHAtom, nHeIAtom, nHeIIatom !, ir, ifreq
    real(double) :: nstar, tauAv, tauHalpha, ratio
    real(double), parameter :: convergeTol = 1.d-4, neTolerance = 1.d-5
    integer :: neIter
    logical :: recalcJbar, neConverged, firstCheckonTau
    character(len=80) :: message
    logical :: ionized
    integer :: iLev, nIter, iLab
    real(double) :: lev1, lev2

#ifdef MPI
    ! For MPI implementations
    integer       ::   my_rank        ! my processor rank
    integer       ::   np             ! The number of processes
    integer       ::   ierr           ! error flag
    integer       ::   nVoxels
    integer       ::   tag=0
    integer, dimension(:), allocatable :: octalsBelongRank
    logical       ::   rankComplete
    real(double), allocatable :: tArrayd(:),tempArrayd(:)
#endif



!    blockHandout = .false. 

#ifdef MPI
    ! FOR MPI IMPLEMENTATION=======================================================
    !  Get my process rank # 
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  
    ! Find the total # of precessor being used in this run
    call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
#endif

    position = VECTOR(0.d0, 0.d0, 0.d0)
    direction = VECTOR(0.d0, 0.d0, 0.d0)
    freq = 0.d0
    indexRBBTrans = 0
    ilev = 0; ilab = 0; indexAtom = 0; lev1 = 0; lev2 = 0
    nfreq = 0; nRBBTrans = 0; tauAv = 0.d0
    call createRBBarrays(nAtom, thisAtom, nRBBtrans, indexAtom, indexRBBTrans)

     call createContFreqArray(nFreq, freq, nAtom, thisAtom, nsource, source, maxFreq)

     dfreq(1) = freq(2)-freq(1)
     do i = 2, nFreq
        dfreq(i) = freq(i) - freq(i-1)
     enddo

     allocate(jnuCont(1:nFreq)) 
     ionized = .false.
     if (grid%geometry == "gammavel") ionized = .true.
     if (grid%geometry == "wrshell") ionized = .true.
    call allocateLevels(grid, grid%octreeRoot, nAtom, thisAtom, nRBBTrans, nFreq, ionized)

    do i = 1, nAtom
       call addCrossSectionstoAtom(thisAtom(i), nFreq, freq)
    enddo


    if (myRankIsZero) &
         call writeAmrGrid("atom_tmp.grid",.false.,grid)

    if (LTE) goto 666

    nHAtom = 0
    nHeIAtom = 0
    nHeIIAtom = 0
    do iAtom = 1, nAtom
       if (thisAtom(iAtom)%name == "HI") nHatom = iAtom
       if (thisAtom(iAtom)%name == "HeI") nHeIatom = iAtom
       if (thisAtom(iAtom)%name == "HeII") nHeIIatom = iAtom
    enddo

! Lyman continuum in detailed balance

    do iAtom = 1, nAtom
       do iTrans = 1, thisAtom(iAtom)%nTrans
          if (thisAtom(iAtom)%name == "HI".and.thisAtom(iAtom)%transType(iTrans) == "RBF") then ! photoionization
             if (thisAtom(iAtom)%iLower(iTrans) == 1) thisAtom(iAtom)%inDetailedBalance(iTrans) = .true.
          endif
       enddo
    enddo

! Lyman continuum in detailed balance

    do iAtom = 1, nAtom
       do iTrans = 1, thisAtom(iAtom)%nTrans
          if (thisAtom(iAtom)%name == "HeI".and.thisAtom(iAtom)%transType(iTrans) == "RBF") then ! photoionization
             if (thisAtom(iAtom)%iLower(iTrans) == 1) thisAtom(iAtom)%inDetailedBalance(iTrans) = .true.
          endif
       enddo
    enddo

! Lyman continuum in detailed balance

    do iAtom = 1, nAtom
       do iTrans = 1, thisAtom(iAtom)%nTrans
          if (thisAtom(iAtom)%name == "HeII".and.thisAtom(iAtom)%transType(iTrans) == "RBF") then ! photoionization
             if (thisAtom(iAtom)%iLower(iTrans) == 1) thisAtom(iAtom)%inDetailedBalance(iTrans) = .true.
             write(*,*) "lyman cont detailed balance ",itrans
          endif
       enddo
    enddo

       

    allocate(oldPops(1:nAtom,1:maxval(thisAtom(1:nAtom)%nLevels)))
    allocate(mainoldPops(1:nAtom,1:maxval(thisAtom(1:nAtom)%nLevels)))
    allocate(newPops(1:nAtom,1:maxval(thisAtom(1:nAtom)%nLevels)))
    allocate(dPops(1:nAtom,1:maxval(thisAtom(1:nAtom)%nLevels)))

    allocate(octalArray(grid%nOctals))
    nOctal = 0
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)

    call init_random_seed()

    call random_seed(size=iSize)
    allocate(iSeed(1:iSize))
    call random_seed(get=iSeed)
#ifdef MPI
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     call MPI_BCAST(iSeed, iSize, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
     if (myRankisZero) then
        open(69, file="cmf_convergence.dat", status="unknown", form="formatted")
        write(69,'(a)') &
!             012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
             "      Iter Max Frac      Tol      Subcell    Level  New pop   Old pop      nRays   Fix ray"
        close(69)
     endif


    nRay = 100
    do iStage = 1, 2


       if (iStage == 1) then
          fixedRays = .true.
          tolerance = 1.d-3
       else
          fixedRays = .false.
          tolerance = 1.d-5
       endif

       gridConverged = .false.
       nIter = 0


       do while (.not.gridConverged)

          nIter = nIter + 1
          if (doTuning) call tune(6, "One cmf iteration")  ! start a stopwatch

          allocate(ds(1:nRay))
          allocate(phi(1:nRay))
          allocate(i0(1:nRBBTrans, 1:nRay))
          allocate(Hcol(1:nRay))
          allocate(HeIcol(1:nRay))
          allocate(HeIIcol(1:nRay))
          allocate(sourceNumber(1:nRay))
          allocate(cosTheta(1:nRay))
          allocate(weight(1:nRay))
          allocate(hitPhotosphere(1:nRay))
          allocate(iCont(1:nRay,1:nFreq))

          if (fixedRays) then
             call random_seed(put=iseed)   ! same seed for fixed rays
          else
             call init_random_seed()
          endif


          ! default loop indecies
          ioctal_beg = 1
          ioctal_end = SIZE(octalArray)       

          

#ifdef MPI
    
    ! we will use an array to store the rank of the process
    !   which will calculate each octal's variables
    allocate(octalsBelongRank(size(octalArray)))
    
    if (my_rank == 0) then
       call mpiBlockHandout(np,octalsBelongRank,blockDivFactor=1,tag=tag,&
                            maxBlockSize=10,setDebug=.false.)
    
    endif
    ! ============================================================================


 if (my_rank /= 0) then
  blockLoop: do     
 call mpiGetBlock(my_rank,iOctal_beg,iOctal_end,rankComplete,tag,setDebug=.false.)
   if (rankComplete) exit blockLoop 
#endif

          do iOctal = ioctal_beg, ioctal_end
!          do iOctal =  ioctal_end, ioctal_beg, -1



!          if (doTuning) call tune(6, "One octal iteration")  ! start a stopwatch
             
             thisOctal => octalArray(iOctal)%content
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

                   if (thisOctal%inflow(subcell)) then
                      nHit = 0
                      do iRay = 1, nRay
                         call getRay(grid, thisOCtal, subcell, position, direction, &
                              ds(iRay), phi(iRay), i0(1:nRBBTrans,iRay), Hcol(iray), HeICol(iRay), HeIICol(iRay),&
                              nAtom, thisAtom, source, nSource, hitPhotosphere(iRay), sourceNumber(iray), &
                              cosTheta(iRay), weight(iRay), nRBBTrans, indexRBBTrans, indexAtom, nHatom, nHeIAtom, nHeIIatom, &
                              nfreq, freq, iCont(iray,1:nFreq))
                         if ((iray == 50).and.(iOctal==iOctal_beg).and.(myrankglobal==2)) &
                              write(*,*) myrankglobal, subcell, " direction ",direction
                         if (hitPhotosphere(iray)) nHit = nHit + 1
                      enddo
                      iter = 0
                      popsConverged = .false.
                      thisOctal%newatomLevel(subcell,:, :) = thisOctal%atomLevel(subcell,:, :)
                      ne = thisOctal%ne(subcell)

                      recalcJbar = .true.
                      firstCheckOnTau = .true.

                      do while (.not.popsConverged)
                         iter = iter + 1
                         mainoldpops = thisOctal%newatomLevel(subcell,1:nAtom,1:)

                         if (recalcJbar) then
                            call calculateJbarCont(thisOctal, subcell, nAtom, thisAtom, ne, nray, ds, source, nSource, &
                                 hitPhotosphere, sourceNumber, freq, nfreq, &
                                 iCont, jnuCont, cosTheta, weight)

                            thisOctal%JnuCont(subcell, :) = jnuCont(:)

                            do iRBB = 1, nRBBTrans
                               iTrans = indexRBBTrans(iRBB)
                               iAtom =  indexAtom(iRBB)
                               call calculateJbar(thisOctal, subcell, thisAtom(iAtom), nRay, ds(1:nRay), &
                                    phi(1:nRay), i0(iRBB,1:nRay), iTrans, thisOctal%jnuLine(subcell,iRBB), &
                                    thisOctal%newAtomLevel(subcell,iAtom,1:), &
                                    freq, nFreq, weight(1:nRay), iRBB,tauAv)


                               if (firstCheckonTau) then
                                  thisAtom(iAtom)%indetailedBalance(iTrans) = .false.
                               endif


                               if (iTrans == 4) tauHalpha = tauAv

                               if (iter < 10) then
                                 if (tauAv > 1.d2) then
                                    thisAtom(iAtom)%indetailedBalance(iTrans) = .true.
                                 else
                                    thisAtom(iAtom)%indetailedBalance(iTrans) = .false.
                                 endif
                                 firstCheckonTau = .false.
                              endif

                              if (debug) write(*,*) thisAtom(iAtom)%iLower(iTrans), " -> ", thisAtom(iAtom)%iUpper(iTrans), ": ", &
                                    tauAv,thisAtom(iAtom)%indetailedbalance(itrans)

                            enddo
                         endif



                         neConverged = .false.
                         neIter = 0
                         oldpops = thisOctal%newatomLevel(subcell,1:nAtom,1:)
                         do while (.not.neConverged)
                            neIter = neIter + 1
                            oldpops = thisOctal%newatomLevel(subcell,1:nAtom,1:)

                            newpops = thisOctal%newatomLevel(subcell, 1:nAtom, 1:)
                            call solveLevels(thisOctal, subcell, newPops, &
                                 thisOctal%jnuLine(subcell,1:nRBBTrans),  &
                                 dble(thisOctal%temperature(subcell)), nAtom, thisAtom, ne, &
                                 thisOctal%rho(subcell),&
                                 jnuCont, freq, dfreq, nfreq)
                            
                            dPops = newPops - thisOctal%newAtomLevel(subcell,1:nAtom,:) 
                            

!                            if (iter > 1)  then
!                               where (abs(dpops/newpops) > 0.1d0) dpops = sign(0.1d0 * newpops,dpops)
!                            endif

                            thisOctal%newAtomLevel(subcell,1:nAtom,:)  = &
                                 thisOctal%newAtomLevel(subcell,1:nAtom,:)  + underCorrect * dPops


                            where (abs(thisOctal%newAtomLevel(subcell,1:nAtom,:)) < 1.d-30)
                               thisOctal%newAtomLevel(subcell,1:nAtom,:) = 1.d-30
                            end where
                            
                            newNe = 0.d0
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
                            
                            dNe = newne - Ne
                            
                            ne = ne + undercorrect * dne

                            fac = -1.d30
                            where (oldPops == 0.d0)
                               oldPops = 1.d-20
                            end where
                            do iAtom = 1, nAtom
                               fac = max(fac,maxval(abs((thisOctal%newAtomLevel(subcell,iatom,1:thisAtom(iAtom)%nLevels-1) &
                                    - oldpops(iAtom,1:thisAtom(iAtom)%nLevels-1))/oldpops(iAtom,1:thisAtom(iAtom)%nLevels-1))))
                            enddo
                            fac = max(fac,abs(dne/ne))
                            if ((fac < neTolerance).and.(neIter > 1)) neConverged = .true.
!                            write(*,*) neiter, fac, dne/ne
                            if (neIter == maxIter) then
                               neConverged = .true.
!                               write(message,'(a,e12.5)') "Maximum number of iterations reached in ne solver. fac = ",fac
!                               call writeWarning(message)
                            endif
                         enddo


                        
                        if (debug) then
                           write(*,*) "Iteration: ",iter
                           do iatom = 1, size(thisAtom)
                              write(*,*) "Atom: ",thisAtom(iatom)%name
                              write(*,*) "Number density: ",thisOctal%rho(subcell)*thisOctal%atomAbundance(subcell, iatom)
                              do i = 1, thisAtom(iatom)%nLevels-1
                                 nStar = boltzSahaGeneral(thisAtom(iAtom), 1, i, ne, &
                                      dble(thisOctal%temperature(subcell))) * &
                                      thisOctal%newAtomLevel(subcell, iatom, thisAtom(iAtom)%nLevels)
                                 ratio = boltzSahaGeneral(thisAtom(iAtom), 1, i, ne, &
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
                           write(*,*) "Ne: ",ne,dne
                           write(*,*) "T: ",thisOctal%temperature(subcell)
                           write(*,*) "R: ",modulus(subcellCentre(thisOctal,subcell))/rCore
                           write(*,*) "log(H I/H II): ",log10(SUM(thisOctal%newAtomLevel(subcell,1,1:thisAtom(1)%nlevels-1)) / &
                                thisOctal%newAtomLevel(subcell, 1, thisAtom(1)%nlevels))
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
                         fac = max(fac,abs(dne/ne))

                         if ((fac < convergeTol).and.(iter > 1)) popsConverged = .true.
                         if (iter == maxIter) then
                            popsConverged = .true.
                            write(message,'(a,e12.5)') "Maximum number of iterations reached in pop solver. fac = ",fac
                            call writeWarning(message)
                         endif
!                         write(*,*) "iter",iter,fac
                      enddo
!                      write(*,*) "Main iteration route converged after ", iter, " iterations"
                      thisOctal%ne(subcell) = ne
                   endif


                endif

!                call locate(freq,nfreq,cspeed/6562.8d-8,ifreq)
!                nStar = boltzSahaGeneral(thisAtom(1), 1, 6, thisOctal%ne(subcell), &
!                     dble(thisOctal%temperature(subcell))) * &
!                     thisOctal%newAtomLevel(subcell, 1, thisAtom(1)%nLevels)
!                write(31,*) real(r*1.e10), real(thisOctal%ne(subcell)), real(thisOctal%newAtomLevel(subcell,1,1:6)),&
!                     real(jnuCont(ifreq)),real(nstar)

!             enddo
!             close(31)
!             stop

             enddo
          end do

!          if (doTuning) call tune(6, "One octal iteration")  ! start a stopwatch

#ifdef MPI
 if (.not.blockHandout) exit blockloop
 end do blockLoop        
 end if ! (my_rank /= 0)




     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
       if(my_rank == 0) write(*,*) "Updating MPI grids"


     ! have to send out the 'octalsBelongRank' array
     call MPI_BCAST(octalsBelongRank,SIZE(octalsBelongRank),  &
                    MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     call countVoxels(grid%octreeRoot,nOctal,nVoxels)
     allocate(tArrayd(1:nVoxels))
     allocate(tempArrayd(1:nVoxels))
     tArrayd = 0.d0
     tempArrayd = 0.d0
     do iAtom = 1, nAtom
       do i = 1, thisAtom(iAtom)%nLevels
         tArrayd = 0.d0
          call packAtomLevel(octalArray, nVoxels, tArrayd, octalsBelongRank, iAtom, i)
          call MPI_ALLREDUCE(tArrayd,tempArrayd,nVoxels,MPI_DOUBLE_PRECISION,&
            MPI_SUM,MPI_COMM_WORLD,ierr)
            tArrayd = tempArrayd
          call unpackAtomLevel(octalArray, nVoxels, tArrayd, octalsBelongRank, iAtom, i)
       enddo
     enddo
     do i = 1, nFreq
       tArrayd = 0.d0
       call packJnu(octalArray, nVoxels, tArrayd, octalsBelongRank, i)
       call MPI_ALLREDUCE(tArrayd,tempArrayd,nVoxels,MPI_DOUBLE_PRECISION,&
            MPI_SUM,MPI_COMM_WORLD,ierr)
       tArrayd = tempArrayd
       call unpackJnu(octalArray, nVoxels, tArrayd, octalsBelongRank, i)
     enddo
     deallocate(tArrayd, tempArrayd)

       if(my_rank == 0) write(*,*) "Done updating"
#endif

          maxFracChange = -1.d30
          call swapPops(grid%octreeRoot, maxFracChange, lev1, lev2, ilev, ilab)
          if (writeoutput) write(*,*) "Maximum fractional change this iteration",maxFracChange
          if (writeoutput) write(*,*) "Fractional change",maxFracChange,"tolerance",tolerance , &
               "fixed rays",fixedrays,"nray",nray
          if (writeoutput) write(*,*) "iLevel ", iLev, " new pops ", lev1, " old pops ", lev2

          if (myRankIsZero) &
               call writeAmrGrid("atom_tmp.grid",.false.,grid)


          if (myRankisZero) then
             open(69, file="cmf_convergence.dat", status="old", position = "append", form="formatted")
             write(69,'(i10, 1p, 2e10.2, 2i10, 2e10.2, i10, l10)') &
                  nIter, maxFracChange, tolerance, ilab, ilev, lev1, lev2, nRay, fixedRays
             close(69)
          endif

          if (maxFracChange < tolerance) then
             gridConverged = .true.
          endif


#ifdef MPI
       deallocate(octalsBelongRank)
#endif

          deallocate(ds, phi, i0, sourceNumber, cosTheta, hitPhotosphere, &
               weight, hcol, heicol, heiicol, iCont)

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

  recursive  subroutine  swapPops(thisOctal, maxFracChange, lev1, lev2, ilev,ilab)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, j, iAtom
    real(double) :: maxFracChange, temp, lev1, lev2
    integer :: ilev, ilab
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call swapPops(child, maxFracChange, lev1, lev2, ilev, ilab)
                exit
             end if
          end do
       else
          do iAtom = 1, size(thisOctal%newAtomLevel,2)
             do j = 1 , 6
                if (thisOctal%atomLevel(subcell,iAtom,j) /= 0.d0) then
                   temp = abs((thisOctal%newatomLevel(subcell,iAtom,j) - &
                        thisOctal%atomLevel(subcell,iAtom,j)) / &
                        thisOctal%atomLevel(subcell,iAtom,j))
                   if (temp > maxFracChange) then
                      maxFracChange = temp
                      ilev = j
                      lev1 = thisOctal%newatomLevel(subcell,iAtom,j)
                      lev2 = thisOctal%atomLevel(subcell,iAtom,j)
                      ilab = thisOctal%label(subcell)
                   endif
                endif
             enddo
          enddo
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
          etaLine = hCgs * a * thisAtom(iatom)%transFreq(iTrans)
          etaLine = etaLine * thisOctal%atomLevel(subcell, iAtom,iUpper)
	  
          thisOctal%etaLine(subcell) = etaLine 
          
       endif
    enddo
  end subroutine calcEtaLine

  recursive subroutine  allocateLevels(grid, thisOctal, nAtom, thisAtom, nRBBTrans, nFreq, ionized)
    type(GRIDTYPE) :: grid
    type(MODELATOM) :: thisAtom(:)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    integer :: nAtom
    integer :: nRBBTrans
    integer :: iatom
    logical :: ionized
    integer :: nFreq

  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call allocateLevels(grid, child, nAtom, thisAtom, nRBBTrans, nfreq, ionized)
                exit
             end if
          end do
       else

          if (.not.associated(thisOctal%atomLevel)) then
             allocate(thisOctal%atomLevel(1:thisOctal%maxChildren, 1:nAtom, 1:maxval(thisAtom(1:nAtom)%nLevels)))
          endif

          if (.not.associated(thisOctal%newatomLevel)) then
             allocate(thisOctal%newatomLevel(1:thisOctal%maxChildren, 1:nAtom, 1:maxval(thisAtom(1:nAtom)%nLevels)))
             thisOctal%newatomLevel = 1.d-10
             if (.not.ionized) then
                do iAtom = 1, nAtom
                   if (thisAtom(iAtom)%charge == 0) then
                      thisOctal%newAtomLevel(1:thisOctal%maxChildren, 1:nAtom, 1) = thisOctal%rho(subcell) &
                           * thisOctal%atomAbundance(subcell, iAtom)
                   endif
                enddo
             endif
          endif

          if (.not.associated(thisOctal%jnuLine)) then
             allocate(thisOctal%jnuLine(1:thisOctal%maxChildren, 1:nRBBTrans))
             thisOctal%jnuLine = 1.d-30
          endif

          if (.not.associated(thisOctal%jnuCont)) then
             allocate(thisOctal%jnuCont(1:thisOctal%maxChildren, 1:nFreq))
             thisOctal%jnuCont = 1.d-30
          endif

          if (ionized) then
             thisOctal%ne(subcell) = thisOctal%rho(subcell)/mHydrogen
          else
             thisOctal%ne(subcell) = 1.d-2 * thisOctal%rho(subcell)/mHydrogen
          endif

          do iAtom = 1, nAtom
             thisOctal%atomLevel(subcell,iAtom,thisAtom(iAtom)%nLevels) = thisOctal%rho(subcell) &
                  * thisOctal%atomAbundance(subcell, iAtom)

             do i = 1, thisAtom(iatom)%nLevels-1
                thisOctal%atomLevel(subcell, iAtom,i) = boltzSahaGeneral(thisAtom(iAtom), 1, i, &
                     thisOctal%ne(subcell), &
                     dble(thisOctal%temperature(subcell))) * &
                     thisOctal%newAtomLevel(subcell, iatom, thisAtom(iAtom)%nLevels)
             enddo
          enddo
          thisOctal%newAtomLevel(subcell,:,:) = thisOctal%atomLevel(subcell,:,:)


       endif
    enddo
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
 
    call random_number(r)

    if (r < probTowardsSource) then
       weight = chanceSource/probTowardsSource
       call random_number(r)
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
       nFreq, freqArray) result (i0)
    use amr_mod, only: distanceToGridFromOutside

    type(VECTOR) :: position, direction, pvec, photoDirection
    type(GRIDTYPE) :: grid
    integer :: nSource
    real(double) :: freqArray(:)
    integer :: nFreq
    type(SOURCETYPE) :: source(:)
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
    type(VECTOR) :: currentPosition, thisPosition, thisVel
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
    real(double) :: distToSource
    integer :: sourcenumber
    integer :: iElement
    logical :: endLoopAtPhotosphere
    real(double) :: nstar(10,50), rhoCol
    real(double) :: bfOpac, bfEmiss, x1, x2, fac
    integer :: j, k
    logical :: lineoff, passThroughResonance, velocityIncreasing


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


    call locate(freqArray, nFreq, thisAtom(iAtom)%transFreq(iTrans), iFreq)

    distToGrid = distanceToGridFromOutside(grid, position, direction)

    if (distToGrid > 1.e29) then
!       write(*,*) "ray does not intersect grid",position,direction
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
    lineOff = .false.

!    if (hitSource) endLoopAtphotosphere = .true.

!    write(*,*) lineoff,hitsource,endloopatphotosphere
   if (.not.lineOff) then



    do while(inOctal(grid%octreeRoot, currentPosition).and.(.not.endloopAtPhotosphere))
       icount = icount + 1 

       call findSubcellLocal(currentPosition, thisOctal, subcell)

!       rVec = subcellCentre(thisOctal,subcell)

       call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal)


       if ((totDist + tval) > distTosource) then
          tVal = distToSource - totDist
          endLoopAtPhotosphere = .true.
       endif



       startVel = amrGridVelocity(grid%octreeRoot, currentPosition, startOctal = thisOctal, actualSubcell = subcell) 
!       startVel = getVel(grid, thisOctal, subcell, currentposition, direction)

       endPosition = currentPosition + tval * direction
       endVel = amrGridVelocity(grid%octreeRoot, endPosition)
!       endVel = getVel(grid, thisOctal, subcell, endposition, direction)

!       endOctal => grid%octreeRoot
!       call findSubcellLocal(endPosition, endOctal, endSubcell)

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
          nTau = 80
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
       bfOpac = 0.d0
       bfEmiss = 0.d0

       do i = 2, nTau

          startOctal => thisOctal
          thisPosition = currentPosition + distArray(i)*direction

          thisVel = amrGridVelocity(grid%octreeRoot, thisPosition, startOctal = startOctal, actualSubcell = subcell) 
!          thisVel = getVel(grid, thisOctal, subcell, thisposition, direction)

          thisVel= thisVel - rayVel


          dv = (thisVel .dot. direction) + deltaV

!          if (abs(dv)*cspeed/1.d5 < 100.d0) then
!             write(*,*) i,distArray(i), dv1*cspeed/1.e5,dv2*cspeed/1.e5,x1,x2,tval,dv*cspeed/1.e5
!          endif
  
          call returnEinsteinCoeffs(thisAtom(iAtom), iTrans, a, Bul, Blu)


          alphanu = (hCgs*thisAtom(iAtom)%transFreq(iTrans)/fourPi) * &
               phiProf(dv, thisOctal%microturb(subcell))/thisAtom(iAtom)%transFreq(iTrans)

          

          iUpper = thisAtom(iAtom)%iUpper(iTrans)
          iLower = thisAtom(iAtom)%iLower(iTrans)


          nLower = thisOctal%atomLevel(subcell,iAtom, iLower)
          nUpper = thisOctal%atomLevel(subcell,iAtom, iUpper)
          alphanu = alphanu * (nLower * Blu - nUpper * Bul)


          if (i == 2) then
             do k = 1, nAtom
                do j = 1, thisAtom(k)%nLevels - 1
                   nStar(k,j) = BoltzSahaGeneral(thisAtom(k), 1, j, thisOctal%ne(subcell), &
                        dble(thisOctal%temperature(subcell))) * &
                        Thisoctal%atomlevel(subcell, k,thisAtom(k)%nLevels)
                enddo
             enddo
             bfOpac = bfOpacity(thisAtom(iAtom)%transFreq(iTrans), nAtom, thisAtom, thisOctal%atomLevel(subcell,:,:), &
                  thisOctal%ne(subcell), nstar, dble(thisOctal%temperature(subcell)))
             bfEmiss = bfEmissivity(thisAtom(iatom)%transFreq(iTrans), nAtom, thisAtom, nstar, &
               dble(thisOctal%temperature(subcell)), thisOctal%ne(subcell), thisOctal%jnuCont(subcell, iFreq))/fourPi
          endif



          etaLine = hCgs * a * thisAtom(iAtom)%transFreq(iTrans)
          etaLine = etaLine * thisOctal%atomLevel(subcell, iAtom, iUpper)
          jnu = (etaLine/fourPi) * phiProf(dv, thisOctal%microturb(subcell))/thisAtom(iAtom)%transFreq(iTrans)


          if (thisOctal%rho(subcell) > 0.1d0) then ! opaque disc
             bfOpac = 1.d30
             bfEmiss = 0.d0
             alphanu = 1.d30
             etaLine = 0.d0
             jnu = 0.d0
             snu = 0.d0
          endif
             
          alphanu = alphanu + bfOpac

          ! add continuous bf and ff emissivity of hydrogen
       
          jnu = jnu + bfEmiss

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

!          write(*,'(i4,i4,l3,10(1pe12.3))') iCount, ntau, passThroughResonance, dv1*cspeed/1.e5,dv2*cspeed/1.e5,dv*cspeed/1.e5, &
!             i0,tau, jnu,alphanu,snu, nlower,nupper

          i0 = i0 +  exp(-tau) * (1.d0-exp(-dtau))*snu
          tau = tau + dtau
       enddo
       rhoCol = rhoCol + distArray(ntau)*thisOctal%rho(subcell)*1.d10
       currentPosition = currentPosition + (distArray(ntau)+1.d-3*grid%halfSmallestSubcell) * direction
       totdist = totdist + (distArray(ntau)+1.d-3*grid%halfSmallestSubcell)
    enddo
 endif


    if (endLoopAtPhotosphere) then

       iElement = getElement(source(sourcenumber)%surface, photoDirection)
!       if (source(sourcenumber)%surface%element(iElement)%hot) then
!          write(*,*) i0,i_nu(source(sourceNumber), thisAtom(iATom)%transFreq(iTrans), iElement), &
!               i_nu(source(sourceNumber), thisAtom(iATom)%transFreq(iTrans), 1)
!       endif
!       write(*,*) i0, i_nu(source(sourceNumber), thisAtom(iATom)%transFreq(iTrans), iElement),tau
!       i0 = i0 + i_nu(source(sourceNumber), thisAtom(iATom)%transFreq(iTrans), iElement)*cosTheta*exp(-tau)
       i0 = i0 + i_nu(source(sourceNumber), thisAtom(iATom)%transFreq(iTrans), iElement)*exp(-tau)


!       write(*,*) "i0",i0,i_nu(source(sourceNumber), thisAtom(iATom)%transFreq(iTrans), iElement)*cosTheta*exp(-tau),tau
!       write(*,*) "pos",modulus(currentPosition)/source(1)%radius
    endif

!    write(*,*) "tau",tau,"deltav",deltaV*cspeed/1.d5

666 continue 

  end function intensityAlongRay

  
  subroutine calculateAtomSpectrum(grid, thisAtom, nAtom, iAtom, iTrans, viewVec, distance, source, nsource, nfile)
    use messages_mod, only : myRankIsZero
    use datacube_mod, only: DATACUBE, writeDataCube, freedatacube
#ifdef MPI
    include 'mpif.h'
#endif

    type(GRIDTYPE) :: grid
    type(MODELATOM) :: thisAtom(:)
    integer :: nSource
    integer :: nFile
    type(SOURCETYPE) :: source(:)
    integer :: nAtom, iAtom
    real(double) :: distance
    integer :: itrans
    integer :: nRay
    type(VECTOR) :: rayPosition(5000)
    real(double) :: da(5000), dOmega(5000)
    type(VECTOR) :: viewVec
    real(double) :: deltaV
    integer :: iv, iray
    integer :: nLambda
    real(double) :: i0
    real(double), allocatable :: vArray(:), spec(:)
    integer :: iv1, iv2, i
    character(len=30) :: plotfile
    type(DATACUBE) :: cube
    integer :: nFreqArray
    integer, parameter :: maxFreq = 2000
    real(double) :: freqArray(maxFreq)

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

    da = 0.d0; dOmega = 0.d0
    cube%label = " "
    freqArray = 0.d0; nFreqArray = 0
    nray = 0; rayPosition = VECTOR(0.d0, 0.d0, 0.d0)
    call createContFreqArray(nFreqArray, freqArray, nAtom, thisAtom, nsource, source, maxFreq)


    if (myRankIsZero) &
         write(*,*) "Calculating spectrum for: ",thisAtom(iatom)%name,(cspeed/thisAtom(iatom)%transFreq(iTrans))*1.d8

    call createDataCube(cube, grid, viewVec, nAtom, thisAtom, iAtom, iTrans, nSource, source, nFreqArray, freqArray)

#ifdef MPI
     write(*,*) "Process ",my_rank, " create data cube done"
#endif
    if (Writeoutput) then
       write(plotfile,'(a,i3.3,a)') "datacube",nfile,".fits.gz"
       call writeDataCube(cube,plotfile)
    endif
    call torus_mpi_barrier
    call freeDataCube(cube)
    return

!    call calcEtaLine(grid%octreeRoot, thisAtom, nAtom, iAtom, iTrans)
    

    call createRayGrid(nRay, rayPosition, da, dOmega, viewVec, distance, grid)

    nLambda = 100

    iv1 = 1
    iv2 = nlambda
 
#ifdef MPI
    iv1 = (my_rank) * (nLambda / (np)) + 1
    iv2 = (my_rank+1) * (nLambda / (np))
    if (my_rank == (np-1)) iv2 = nLambda
#endif

    allocate(spec(1:nLambda), vArray(1:nLambda))
    spec = 0.d0
    do iv = 1, nLambda
       vArray(iv) = 2200.e5/cspeed * (2.d0*dble(iv-1)/dble(nLambda-1)-1.d0)
    enddo
    do iv = iv1, iv2
       write(*,*) iv
       deltaV  = vArray(iv)
       do iRay = 1, nRay
          i0 = intensityAlongRay(rayposition(iRay), viewvec, grid, thisAtom, nAtom, iAtom, iTrans, -deltaV, source, nSource, &
               nFreqArray, freqArray) !minus v 
          spec(iv) = spec(iv) + i0 * domega(iRay) 
       enddo
    enddo
#ifdef MPI
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     allocate(tempArray(1:nLambda))
     call MPI_ALLREDUCE(spec,tempArray,nLambda,MPI_DOUBLE_PRECISION,&
          MPI_SUM,MPI_COMM_WORLD,ierr)
     spec(1:nLambda) = tempArray(1:nLambda)
     deallocate(tempArray)
#endif

    if (myRankIsZero) then
       open(42, file="spectrum.dat",status="unknown",form="formatted")
       do i = 1, nLambda
          write(42, *) vArray(i)*cspeed/1.d5, spec(i)
       enddo
       close(42)
    endif
    deallocate(vArray, spec)
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

    nr1 = 5
    nr2 = 20
    nr = nr1 + nr2
    nphi = 20
    nray = 0
    i = 0

    allocate(rGrid(1:nr), dr(1:nr), phiGrid(1:nPhi), dphi(1:nPhi))
    rmin = grid%rCore
    rMax = grid%rOuter

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
       call random_number(phiOffset)
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
          dOmega(nRay) = da(nRay) / (fourPi * distance**2)
       enddo
    enddo
  end subroutine createRayGrid


  subroutine createDataCube(cube, grid, viewVec, nAtom, thisAtom, iAtom, iTrans, nSource, source, &
       nFreqArray, freqArray)
    use input_variables, only : cylindrical, ttauriRouter
    use datacube_mod, only: DATACUBE, initCube, addspatialaxes, addvelocityAxis
#ifdef MPI
    include 'mpif.h'
#endif

    integer :: nSource
    type(SOURCETYPE) :: source(:)
    type(MODELATOM) :: thisAtom(:)
    real(double) :: freqArray(:)
    integer :: nFreqArray
    integer :: iAtom
    integer :: nAtom
    type(GRIDTYPE) :: grid
    type(DATACUBE) :: cube
    type(VECTOR) :: viewvec, rayPos, xProj, yProj
    real(double) :: deltaV
    integer :: iTrans
    integer :: ix, iy, iv
    real(double) :: r, xval, yval
    integer :: nMonte, imonte
    integer :: iv1, iv2

#ifdef MPI
    ! For MPI implementations
    integer       ::   my_rank        ! my processor rank
    integer       ::   np             ! The number of processes
    integer       ::   ierr           ! error flag
    integer       ::   n
    real(double), allocatable :: tempArray(:), tempArray2(:)

    ! FOR MPI IMPLEMENTATION=======================================================
    !  Get my process rank # 
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  
    ! Find the total # of precessor being used in this run
    call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
#endif

    nMonte = 1

    call initCube(cube, 200, 200, 2)
!    call addSpatialAxes(cube, -grid%octreeRoot%subcellSize*0.1d0, +grid%octreeRoot%subcellSize*0.1d0, &
!         -grid%octreeRoot%subcellSize*0.1d0, grid%octreeRoot%subcellSize*0.1d0)

!    call addSpatialAxes(cube, -grid%octreeRoot%subcellSize*1.9d0, +grid%octreeRoot%subcellSize*1.9d0, &
!         -grid%octreeRoot%subcellSize*1.9d0, grid%octreeRoot%subcellSize*1.9d0)

    if (grid%octreeRoot%threed) then
       if (.not.cylindrical) then
          call addSpatialAxes(cube, -grid%octreeRoot%subcellSize*0.9d0, +grid%octreeRoot%subcellSize*0.9d0, &
               -grid%octreeRoot%subcellSize*0.9d0, grid%octreeRoot%subcellSize*0.9d0)
       else
!          call addSpatialAxes(cube, -grid%octreeRoot%subcellSize*1.9d0, +grid%octreeRoot%subcellSize*1.9d0, &
!               -grid%octreeRoot%subcellSize*1.9d0, grid%octreeRoot%subcellSize*1.9d0)
          call addSpatialAxes(cube, -ttauriRouter*1.5d0/1.d10, ttauriRouter*1.5d0/1.d10, &
               -ttauriRouter*1.5d0/1.d10, ttauriRouter*1.5d0/1.d10)
!          call addSpatialAxes(cube, -2.d0*rsol/1.d10, 2.d0*rsol/1.d10, -2.d0*rsol/1.d10,  2.d0*rsol/1.d10)
       endif
    endif
!    call addSpatialAxes(cube, -grid%octreeRoot%subcellSize/1.d6, +grid%octreeRoot%subcellSize/1.d6, &
!         -grid%octreeRoot%subcellSize/1.d6, grid%octreeRoot%subcellSize/1.d6)
!    call addSpatialAxes(cube, -dble(3.1*grid%rInner), +dble(3.1*grid%rInner), -dble(3.1*grid%rInner), +dble(3.1*grid%rInner))
!    write(*,*) "rinner",grid%rinner/(rsol/1.e10)

    call addvelocityAxis(cube, -200.d0, 200.d0)

    xProj =   viewVec .cross. VECTOR(0.d0, 0.d0, 1.d0)
    call normalize(xProj)
    yProj =  xProj .cross.viewVec
    call normalize(yProj)
    iv1 = 1
    iv2 = cube%nv
 

#ifdef MPI
    iv1 = int(real(my_rank) * (real(cube%nv) / real(np))) + 1
    iv2 = int(real(my_rank+1) * (real(cube%nv) / real(np)))
    if (my_rank == (np-1)) iv2 = cube%nv
#endif

!    iv1 = 25

    do iv = iv1, iv2
       write(*,*) iv,iv1,iv2
       do ix = 1, cube%nx
          do iy = 1, cube%ny
             do iMonte = 1, nMonte
                if (nMonte > 1) then
                   call random_number(r)
                   xVal = cube%xAxis(ix) + (r-0.5d0)*(cube%xAxis(2)-cube%xAxis(1))
                   call random_number(r)
                   yVal = cube%yAxis(iy) + (r-0.5d0)*(cube%yAxis(2)-cube%yAxis(1))
                else
                   xVal = cube%xAxis(ix)
                   yVal = cube%yAxis(iy)
             endif
             rayPos =  (xval * xProj) + (yval * yProj)
             raypos = rayPos + ((-1.d0*grid%octreeRoot%subcellsize*3.d0) * Viewvec)
             deltaV = cube%vAxis(iv)*1.d5/cSpeed


                cube%intensity(ix,iy,iv) = intensityAlongRay(rayPos, viewVec, grid, thisAtom, nAtom, iAtom, &
                     iTrans, -deltaV, source, nSource, nFreqArray, freqArray)
!                write(*,*) ix,iy,iv,cube%intensity(ix,iy,iv)
             enddo
          enddo
       enddo
       cube%intensity(:,:,iv) = cube%intensity(:,:,iv) / dble(nMonte)
    enddo
#ifdef MPI
     write(*,*) "Process ",my_rank, " done. awaiting reduce"
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
    do iv = 1, cube%nv
      do ix = 1, cube%nx
        n = (cube%ny)
        allocate(tempArray(1:n), tempArray2(1:n))
        tempArray = reshape(cube%intensity(ix,:,iv), (/  n /))
         call MPI_ALLREDUCE(tempArray,tempArray2,n,MPI_DOUBLE_PRECISION,&
             MPI_SUM,MPI_COMM_WORLD,ierr)
         cube%intensity(ix,:,iv) = reshape(tempArray2, (/ cube%ny/))
         deallocate(tempArray, tempArray2)
       enddo
    enddo
    write(*,*) "Process ",my_rank, " reduce done."
#endif

 cube%flux = cube%intensity

  end subroutine createDataCube



#ifdef MPI
      subroutine packAtomLevel(octalArray, nTemps, tArray, octalsBelongRank, iAtom, iLevel)
    include 'mpif.h'
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: octalsBelongRank(:)
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
                 if (octalsBelongRank(iOctal) == my_rank) then
                   tArray(nTemps) = thisOctal%newAtomLevel(isubcell, iAtom, iLevel)
                 else 
                   tArray(nTemps) = 0.d0
                 endif
              endif
          end do
       end do
     end subroutine packAtomLevel

      subroutine unpackAtomLevel(octalArray, nTemps, tArray, octalsBelongRank, iAtom, iLevel)
    include 'mpif.h'
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: octalsBelongRank(:)
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

      subroutine packjnu(octalArray, nTemps, tArray, octalsBelongRank, iFreq)
    include 'mpif.h'
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: octalsBelongRank(:)
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
                 if (octalsBelongRank(iOctal) == my_rank) then
                   tArray(nTemps) = thisOctal%jnuCont(isubcell, ifreq)
                 else 
                   tArray(nTemps) = 0.d0
                 endif
              endif
          end do
       end do
     end subroutine packJnu

      subroutine unpackJnu(octalArray, nTemps, tArray, octalsBelongRank, iFreq)
    include 'mpif.h'
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: octalsBelongRank(:)
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

end module cmf_mod
