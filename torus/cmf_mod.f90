module cmf_mod

  ! written by tjh


  use kind_mod
  use constants_mod
  use utils_mod
  use messages_mod
  use grid_mod
  use math_mod
  use modelatom_mod
  use source_mod

  !$MPI    use parallel_mod

  implicit none

  private
  public:: atomLoop

contains

  subroutine solveLevels(nPops, jnuLine,  temperature, nAtom, thisAtom, ne, nh, jnuCont, freq, nfreq)
    integer :: nFreq
    real(double) :: freq(:), jnuCont(:)
    real(double) :: nPops(:,:)
    real(double) :: temperature, ne
    real(double) :: jnuLine(:)
    type(MODELATOM) :: thisAtom(:)
    real(double), allocatable :: matrixA(:,:), matrixB(:), collMatrix(:,:), cTot(:)
    real(double) :: arateji, boltzFac
    integer :: nLevels
    integer :: iLower, iUpper, iLevel, i, j
    integer :: itrans, l, k, iPart
    real(double) :: collEx, colldeEx
    real(double) :: a, Bul, Blu
    real(double) :: photoRatelk, recombRatekl, xSection
    real(double) :: partitionFunc
    real(double) :: fac, nh
    real(double) :: prate, rrate, Nstar
    integer :: iJnu

    integer :: nAtom, iAtom, nMatrix
    integer :: nOffset


    nMatrix = 0
    nMatrix = nMatrix + thisAtom(1)%nLevels-1
    nMatrix = nMatrix + thisAtom(2)%nLevels
    nMatrix = nMatrix + 1
       
    allocate(matrixA(1:nMatrix,1:nMatrix))
    allocate(matrixB(1:nMatrix))

    matrixA = 0.d0
    matrixB = 0.d0

    
       
    nOffset  = 0

    do iAtom = 1, nAtom

       nLevels = thisAtom(iAtom)%nLevels

       ! RHS of equation

       do iTrans = 1, thisAtom(iAtom)%nTrans

          k = thisAtom(iAtom)%iUpper(iTrans)
          l = thisAtom(iAtom)%iLower(iTrans)

          Nstar = boltzSahaGeneral(thisAtom(iAtom), 1, l, ne, temperature, nPops(iAtom,thisAtom(iATom)%nLevels))

          if (thisAtom(iAtom)%transType(iTrans) == "RBF") then ! radiative recomb
             recombratekl = 0.d0
             do i = 2, nFreq
                xSection = photoCrossSection(thisAtom(iAtom), l, freq(i))
                recombRatekl = recombRatekl + &
                     (fourPi/(hCgs*freq(i)))*xSection*((2.d0*hCgs*freq(i)**3)/cSpeed**2 + jnuCont(i)) * &
                     exp(-(hCgs*freq(i))/(kErg*temperature))*(freq(i)-freq(i-1))
             enddo
             matrixB(l+nOffset) = matrixB(l+nOffset) - Nstar * recombratekl
          endif
          if (thisAtom(iAtom)%transType(iTrans) == "CBF") then ! collisional recomb

             collEx = collisionRate(thisAtom(iAtom), iTrans, temperature)
             matrixB(l+nOffset) = matrixB(l+nOffset) - Nstar*collEx*ne
          endif
       enddo

       matrixB(nLevels+nOffset) = matrixB(nLevels+nOffset)-SUM(matrixB(1+nOffset:nLevels-1+nOffset))

       ! LHS

       do iTrans = 1, thisAtom(iAtom)%nTrans

          k = thisAtom(iAtom)%iUpper(iTrans)
          l = thisAtom(iAtom)%iLower(iTrans)
          if (thisAtom(iAtom)%transType(iTrans) == "RBB") then ! radiative BB rates
             iJnu = thisAtom(iAtom)%indexRBBTrans(iTrans)
             call returnEinsteinCoeffs(thisAtom(iAtom), iTrans, a, Bul, Blu)

             matrixA(k+nOffset,k+nOffset) = matrixA(k+nOffset,k+nOffset) - Bul * jnuLine(iJnu) - a
             matrixA(l+nOffset,l+nOffset) = matrixA(l+nOffset,l+nOffset) - Blu * jnuLine(iJnu)
             matrixA(k+nOffset,l+nOffset) = matrixA(k+nOffset,l+nOffset) + Blu * jnuLine(iJnu)
             matrixA(l+nOffset,k+nOffset) = matrixA(l+nOffset,k+nOffset) + Bul * jnuLine(iJnu) + a
          endif

          if (thisAtom(iAtom)%transType(iTrans) == "CBB") then ! collision BB  rates


             collEx = collisionRate(thisAtom(iAtom), iTrans, temperature)
             boltzFac =  exp(-abs(thisAtom(iAtom)%energy(k)-thisAtom(iAtom)%energy(l)) / (kev*temperature))
             colldeEx = collEx / (boltzFac * thisAtom(iAtom)%g(k)/thisAtom(iAtom)%g(l))

             matrixA(k+nOffset,k+nOffset) = matrixA(k+nOffset,k+nOffset) - colldeex * ne
             matrixA(l+nOffset,l+nOffset) = matrixA(l+nOffset,l+nOffset) - collex * ne
             matrixA(k+nOffset,l+nOffset) = matrixA(k+nOffset,l+nOffset) + collex * ne
             matrixA(l+nOffset,k+nOffset) = matrixA(l+nOffset,k+nOffset) + colldeex * ne

          endif

          ! now do bound-free rates


          if (thisAtom(iAtom)%transType(iTrans) == "RBF") then ! photoionization


             l = thisAtom(iAtom)%iLower(iTrans)
             k = thisAtom(iAtom)%iUpper(iTrans)

             photoRatelk = 0.d0
             do i = 2, nFreq
                xSection = photoCrossSection(thisAtom(iAtom), l, freq(i))
                photoRatelk = photoRatelk + (jnuCont(i)/(hCgs*freq(i)))*xSection*(freq(i)-freq(i-1))
             enddo
             photoRatelk = photoRatelk * fourPi

!             recombratekl = 0.d0
!             do i = 2, nFreq
!                xSection = photoCrossSection(thisAtom, l, freq(i))
!                recombRatekl = recombRatekl + &
!                     (fourPi/(hCgs*freq(i)))*xSection*((2.d0*hCgs*freq(i)**3)/cSpeed**2 + jnuCont(i)) * &
!                     exp(-(hCgs*freq(i))/(kErg*temperature))*(freq(i)-freq(i-1))
!             enddo

             matrixA(l+nOffset,l+nOffset) = matrixA(l+nOffset,l+nOffset) - photoRatelk
             matrixA(k+nOffset,l+nOffset) = matrixA(k+nOffset,l+nOffset) + photoRatelk


          endif


          if (thisAtom(iAtom)%transType(iTrans) == "CBF") then ! collisional ionization

             k = thisAtom(iAtom)%iUpper(iTrans)
             l = thisAtom(iAtom)%iLower(iTrans)


             collEx = collisionRate(thisAtom(iAtom), iTrans, temperature)

             matrixA(l+nOffset,l+nOffset) = matrixA(l+nOffset,l+nOffset) - collex * ne
             matrixA(k+nOffset,l+nOffset) = matrixA(k+nOffset,l+nOffset) + collex * ne

          endif


       enddo
       nOffset = nOffset + thisAtom(iAtom)%nLevels-1
    enddo

    matrixB(nMatrix) = nh ! total number of atoms
    
 
    nOffset = 0
    matrixA(nMatrix,1:thisAtom(1)%nLevels-1) = 1.d0
    nOffset = nOffset + thisAtom(1)%nLevels-1
    matrixA(nMatrix,1+nOffset:thisAtom(2)%nLevels+noffset) = 1.d0

    
    
    do i = 1, nMatrix
       write(*,'(1p,100e8.1)') matrixA(i,1:nMatrix),matrixB(i)
    enddo
    

    call luSlv(matrixA, matrixB, nMatrix)
    
    nOffset  = 0
    nPops(1,1:thisAtom(1)%nLevels-1) = matrixB(1:thisAtom(1)%nLevels-1)
    nOffset = nOffset + thisAtom(1)%nLevels-1
    nPops(2,1:thisAtom(2)%nLevels) = matrixB(1+nOffset:thisAtom(2)%nLevels+nOffset)

    deallocate(matrixA, matrixB)

!    write(*,*) "H",SUM(nPops(1,1:thisAtom(1)%nLevels-1))
!    write(*,*) "NP",npops(1,thisAtom(1)%nLevels)

    write(*,*) "HeI",SUM(nPops(1,1:thisAtom(1)%nLevels-1))
    write(*,*) "HeII",SUM(nPops(2,1:thisAtom(2)%nLevels-1))
    write(*,*) "HeIII",nPops(2,thisAtom(2)%nLevels)


  end subroutine solveLevels

  subroutine getRay(grid, fromOctal, fromSubcell, position, direction, ds, phi, i0, nAtom, thisAtom, source, nSource, &
       hitPhotosphere, sourceNumber, cosTheta, weight, nRBBTrans, indexRBBTrans, indexAtom)
    type(SOURCETYPE) :: source(:)
    integer :: nAtom
    integer :: nSource
    integer :: nRBBTrans, indexRBBTrans(:), indexAtom(:)
    type(MODELATOM) :: thisAtom(:)
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, startOctal, fromOctal
    integer :: fromSubcell
    integer :: subcell
    real(double) :: ds, phi, i0(:), r, phi1, phi2, cosTheta
    integer :: iTrans
    type(OCTALVECTOR) :: position, direction, currentPosition, thisPosition, thisVel
    type(OCTALVECTOR) :: rayVel, startVel, endVel, endPosition, pvec
    real(double) :: alphanu, snu, jnu
    integer :: iLower , iUpper
    real(double) :: dv, deltaV
    integer :: i
    real(double) :: distArray(200), tval
    integer :: nTau
    real(double) :: nLower, nUpper
    real(double) :: dTau, etaline, didtau
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
   endif


    rayVel = amrGridVelocity(grid%octreeRoot, position, startOctal = thisOctal, actualSubcell = subcell)

    call random_number(r)

    deltaV = 4.3 * thisOctal%microturb(subcell) * (r-0.5d0)

    deltaV = deltaV +  (rayVel .dot. direction)

    projVel = deltaV - (rayVel .dot. direction)

    phi1 = phiProf(projVel, thisOctal%microturb(subcell))

    call distanceToCellBoundary(grid, position, direction, ds, sOctal=thisOctal)

    if (ds > disttosource) ds = disttoSource

    currentPosition = position + ds * direction

    totDist = totDist  + ds

    phi2 = phi1 
    if (inOctal(grid%octreeRoot, currentPosition)) then
       thisVel = amrGridVelocity(grid%octreeRoot, currentPosition, startOctal = thisOctal) 
       projVel = deltaV - (thisVel.dot.direction)
       
       phi2 = phiProf(projVel, thisOctal%microturb(subcell))
    else
       phi2 = phi1
    endif

    phi  = (phi1 + phi2)/2.d0

    ds = ds * 1.d10



    i0 = 0.d0
    intensityIntegral = 0.0
    tau = 0.d0

    hitPhotosphere = .false.

    do while(inOctal(grid%octreeRoot, currentPosition).and.(.not.hitPhotosphere))

       call findSubcellLocal(currentPosition, thisOctal, subcell)
       call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal)

       if ((totDist + tval) > distTosource) then
          tVal = distToSource - totDist
          hitPhotosphere = .true.
       endif

       startVel = amrGridVelocity(grid%octreeRoot, currentPosition, startOctal = thisOctal, actualSubcell = subcell) 
       endPosition = currentPosition + tval * direction
       endVel = amrGridVelocity(grid%octreeRoot, endPosition)

       dvAcrossCell = ((startVel - rayVel).dot.direction) - ((endVel - rayVel).dot.direction)
       dvAcrossCell = abs(dvAcrossCell / thisOctal%microturb(subcell))
       dv1 = abs(deltaV - (startVel .dot. direction))
       dv2 = abs(deltaV - (endVel .dot. direction))

       nTau = 2
       
       if (max(dv1,dv2) < 4.d0*thisOctal%microturb(subcell)) then
          nTau = min(max(2, nint(dvAcrossCell * 5.d0)),200)
       endif

!       write(*,*) ntau,modulus(startVel)*cspeed/1.e5,thisOCtal%microturb(subcell)*cspeed/1.e5,tval, &
!            dv1*cspeed/1.e5,dv2*cspeed/1.e5
       distArray(1) = 0.d0
       do i = 2, nTau
          
          distArray(i) = tval * dble(i-1)/dble(nTau-1)

          startOctal => thisOctal
          thisPosition = currentPosition + distArray(i)*direction
          thisVel = amrGridVelocity(grid%octreeRoot, thisPosition, startOctal = startOctal, actualSubcell = subcell) 


          dv = deltaV - (thisVel .dot. direction)

          icount = icount + 1

          do iRBB = 1, nRBBTRans
             iAtom = indexAtom(iRBB)
             iTrans = indexRBBTrans(iRBB)

             if (thisAtom(iAtom)%transType(iTrans) == "RBB") then

                call returnEinsteinCoeffs(thisAtom(iAtom), iTrans, a, Bul, Blu)

                alphanu = (hCgs*thisAtom(iAtom)%transFreq(iTrans)/fourPi) * &
                     phiProf(dv, thisOctal%microturb(subcell))/thisAtom(iAtom)%transFreq(iTrans)
                
                iUpper = thisAtom(iAtom)%iUpper(iTrans)
                iLower = thisAtom(iAtom)%iLower(iTrans)
                
                nLower = thisOctal%atomLevel(subcell,iAtom, iLower)
                nUpper = thisOctal%atomLevel(subcell,iAtom, iUpper)
                
                alphanu = alphanu * (nLower * Blu - nUpper * Bul)
                
                dTau = alphaNu *  (distArray(i)-distArray(i-1)) * 1.d10
                
                etaLine = hCgs * a * thisAtom(iAtom)%transFreq(iTrans)
                etaLine = etaLine * thisOctal%atomLevel(subcell, iAtom, iUpper)
                jnu = (etaLine/fourPi) * phiProf(dv, thisOctal%microturb(subcell))/thisAtom(iAtom)%transFreq(iTrans)
                
                if (alphanu /= 0.d0) then
                   snu = jnu/alphanu
                else
                   snu = tiny(snu)
                endif
                i0(iRBB) = i0(iRBB) +  exp(-tau(irBB)) * (1.d0-exp(-dtau))*snu
                tau(iRBB) = tau(iRBB) + dtau
                

             endif
          enddo
       enddo
       currentPosition = currentPosition + (distArray(ntau)+1.d-3*grid%halfSmallestSubcell) * direction
       totdist = totdist + (distArray(ntau)+1.d-3*grid%halfSmallestSubcell)
    enddo

    deallocate(tau)
  end subroutine getRay

  function phiProf(dv, b) result (phi)
    real(double) :: dv, b
    real(double) :: fac, phi
    phi = 1.d0 / (b * sqrt(Pi))
    fac = dv**2 / b**2
    phi = phi * exp(-fac)
  end function phiProf

  subroutine calculateJbar(thisOctal, subcell, thisAtom, nRay, ds, phi, i0, iTrans, jbar, nPops, &
       freq, nfreq, jnuCont, weight)
    real(double) :: freq(:), jnuCont(:), weight(:)
    logical :: doCont
    integer :: nfreq
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    type(MODELATOM) :: thisAtom
    integer :: nRay
    real(double) :: ds(:), phi(:), i0(:), nPops(:)
    integer :: iTrans
    real(double) :: jbar
    integer :: iRay
    real(double) :: nLower, nUpper
    real(double) :: jBarInternal, jBarExternal
    real(double) :: alphanu, jnu, etaline
    integer :: iUpper, iLower
    real(double) :: tau, snu, sumPhi
    real(double) :: a, bul, blu
    integer :: i
    logical :: first

    first = .true.
    jBarExternal = 0.d0
    jBarInternal = 0.d0

    if (thisAtom%transType(iTrans) == "RBB") then

       iUpper = thisAtom%iUpper(iTrans)
       iLower = thisAtom%iLower(iTrans)

       sumPhi = 0.d0
       do iRay = 1, nRay
          alphanu = (hCgs*thisAtom%transFreq(iTrans)/fourPi)
          nLower = nPops(iLower)
          nUpper = nPops(iUpper)

          call returnEinsteinCoeffs(thisAtom, iTrans, a, Bul, Blu)

          alphanu = alphanu * (nLower * Blu - nUpper * Bul) * phi(iray)/thisAtom%transFreq(iTrans)

          if (alphanu < 0.d0) then
             alphanu = 0.d0
             if (first) then
                write(*,*) "negative opacity warning",iUpper,iLower
                first = .false.
             endif
          endif

          tau = alphaNu * ds(iray)

          etaLine = hCgs * a * thisAtom%transFreq(iTrans)
          etaLine = etaLine *  nPops(iUpper)
          jnu = (etaLine/fourPi) * phi(iRay)/thisAtom%transFreq(iTrans)

          if (alphanu /= 0.d0) then
             snu = jnu/alphanu
          else
             snu = tiny(snu)
          endif
          
          jBarExternal = jBarExternal + i0(iray) * exp(-tau) * phi(iRay) * weight(iRay)
          jBarInternal = jBarInternal + snu * (1.d0 - exp(-tau)) * phi(iRay) * weight(iRay)
 
          sumPhi = sumPhi + phi(iRay) * weight(iRay)
       enddo

       jbar = (jBarExternal + jBarInternal)/sumPhi
       call locate(freq, nfreq, thisAtom%transFreq(iTrans), i)
       jBar = jBar  + jnuCont(i)

    endif
          
  end subroutine calculateJbar

  subroutine calculateJbarCont(nray, source, nSource, hitPhotosphere, sourceNumber, freq, nfreq, jnuCont, cosTheta, weight)
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    logical :: hitPhotosphere(:)
    integer :: sourceNumber(:)
    integer :: nfreq
    real(double) :: freq(:), jnuCont(:), cosTheta(:), weight(:)
    integer :: iray, nRay, i

    jnuCont = 0.d0
    rayloop: do iRay = 1, nRay
       if (.not.hitPhotosphere(iRay)) then
          cycle rayloop
       endif
       do i = 1, nFreq
          jnuCont(i) = jnuCont(i) + weight(iRay)*i_nu(source(sourceNumber(iray)), freq(i))*cosTheta(iRay)
       enddo
    end do rayloop
    jnuCont = jnuCont / SUM(weight(1:nRay))
  end subroutine calculateJbarCont


  subroutine atomLoop(grid, nAtom, thisAtom, nSource, source)

    use input_variables, only : blockhandout

    type(SOURCETYPE) :: source(:)
    integer :: nSource
    integer :: nAtom
    type(GRIDTYPE) :: grid
    type(MODELATOM) :: thisAtom(:)
    type(OCTALVECTOR) :: position, direction
    integer :: nOctal, iOctal, subcell
    real(double), allocatable :: ds(:), phi(:), i0(:,:), tau(:)
    integer :: nRay
    type(OCTAL), pointer :: thisOctal
    integer, parameter :: maxIter = 100, maxRay = 100000
    logical :: popsConverged, gridConverged 
    integer :: iRay, iTrans, iter,i 
    integer :: iStage
    real(double), allocatable :: oldpops(:,:)
    real(double) :: fac
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    logical :: dcAllocated
    integer, dimension(:), allocatable :: octalsBelongRank
    real(double), allocatable :: tArrayd(:),tempArrayd(:)
    integer :: nVoxels
    integer :: ioctal_beg, ioctal_end, tag = 0
    real(double) :: maxFracChange
    logical :: fixedRays
    integer :: isize
    integer, allocatable :: iseed(:)
    real(double) :: tolerance
    integer, allocatable :: sourceNumber(:)
    real(double), allocatable :: cosTheta(:)
    real(double), allocatable :: weight(:)
    logical, allocatable :: hitPhotosphere(:)
    integer, parameter :: maxFreq = 2000
    real(double) :: freq(maxFreq), jnuCont(maxFreq)
    integer :: nFreq, nhit, iRBB
    integer :: nRBBTrans
    integer :: indexRBBTrans(1000), indexAtom(1000)
    real(double) :: nuStart, nuEnd
    integer :: iAtom

    nuStart = cSpeed / (100000.d0 * 1.d-8)
    nuEnd = cSpeed / (50.d0 * 1.d-8)

    nfreq = maxFreq
    do i = 1, nFreq
       freq(i) = log10(nuStart) + (dble(i-1)/dble(nFreq-1))*(log10(nuEnd)-log10(nuStart))
    enddo
    freq(1:nFreq) = 10.d0**freq(1:nFreq)

    blockHandout = .false. 

    call createRRBarrays(nAtom, thisAtom, nRBBtrans, indexAtom, indexRBBTrans)

    call allocateLevels(grid, grid%octreeRoot, nAtom, thisAtom, nRBBTrans)


    !    call solveAllPops(grid, grid%octreeRoot, thisAtom, jnuCont, freq, nfreq)


    allocate(oldPops(1:nAtom,1:maxval(thisAtom(1:nAtom)%nLevels)))

    allocate(octalArray(grid%nOctals))
    nOctal = 0
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)

    call random_seed

    call random_seed(size=iSize)
    allocate(iSeed(1:iSize))
    call random_seed(get=iSeed)

    nRay = 1000

    do iStage = 1, 2


       if (iStage == 1) then
          fixedRays = .true.
          tolerance = 1.d-2
       else
          fixedRays = .false.
          tolerance = 1.d-4
       endif

       gridConverged = .false.

       do while (.not.gridConverged)

          allocate(ds(1:nRay))
          allocate(phi(1:nRay))
          allocate(i0(1:nRBBTrans, 1:nRay))
          allocate(tau(1:nRay))
          allocate(sourceNumber(1:nRay))
          allocate(cosTheta(1:nRay))
          allocate(weight(1:nRay))
          allocate(hitPhotosphere(1:nRay))


          if (fixedRays) then
             call random_seed(put=iseed)   ! same seed for fixed rays
          else
             call random_seed
          endif

          ! default loop indecies
          ioctal_beg = 1
          ioctal_end = SIZE(octalArray)       



          do iOctal = ioctal_beg, ioctal_end
             thisOctal => octalArray(iOctal)%content
             do subcell = 1, thisOctal%maxChildren

                if (.not.thisOctal%hasChild(subcell)) then
                   if (modulus(subcellCentre(thisOctal,subcell)) > source(1)%radius) then
                      nHit = 0
                      do iRay = 1, nRay
                         call getRay(grid, thisOCtal, subcell, position, direction, &
                              ds(iRay), phi(iRay), i0(1:nRBBTrans,iRay), &
                              nAtom, thisAtom, source, nSource, hitPhotosphere(iRay), sourceNumber(iray), &
                              cosTheta(iRay), weight(iRay), nRBBTrans, indexRBBTrans, indexAtom)
                         if (hitPhotosphere(iray)) nHit = nHit + 1
                      enddo
                      write(*,*) "fraction of rays hitting core",real(nhit)/real(nray)
                      iter = 0
                      popsConverged = .false.
                      thisOctal%newatomLevel(subcell,:, :) = thisOctal%atomLevel(subcell,:, :)


                      do while (.not.popsConverged)
                         iter = iter + 1
                         oldpops = thisOctal%newatomLevel(subcell,1:nAtom,1:)

                         call calculateJbarCont(nray, source, nSource, hitPhotosphere, sourceNumber, &
                              freq, nfreq, jnuCont, cosTheta, weight)

                         do iRBB = 1, nRBBTrans
                            iTrans = indexRBBTrans(iRBB)
                            iAtom =  indexAtom(iRBB)
                            call calculateJbar(thisOctal, subcell, thisAtom(iAtom), nRay, ds(1:nRay), &
                                 phi(1:nRay), i0(iRBB,1:nRay), iTrans, thisOctal%jnu(subcell,iRBB), &
                                 thisOctal%newAtomLevel(subcell,iAtom,1:), &
                                 freq, nFreq, jnuCont, weight(1:nRay))

                         enddo

                         call solveLevels(thisOctal%newatomLevel(subcell, 1:nAtom, 1:), &
                              thisOctal%jnu(subcell,1:nRBBTrans),  &
                              dble(thisOctal%temperature(subcell)), nAtom, thisAtom, thisOctal%ne(subcell), &
                              thisOctal%rho(subcell)/mHydrogen,&
                              jnuCont, freq, nfreq)
                         write(*,*) "ne fixed!"
                         do iatom = 1, nAtom
                            write(*,*) "Atom: ",thisAtom(iatom)%name
                            do i = 1, thisAtom(iatom)%nLevels
                               write(*,*) i,thisOctal%newAtomLevel(subcell,iatom,i),oldpops(iatom,i)
                            enddo
                         enddo
                         fac = -1.d30
                         do iAtom = 1, nAtom
                            fac = max(fac,maxval(abs((thisOctal%newAtomLevel(subcell,1:nAtom,1:thisAtom(iAtom)%nLevels) &
                                 - oldpops(1:nAtom,1:thisAtom(iAtom)%nLevels))/oldpops(1:nAtom,1:thisAtom(iAtom)%nLevels))))
                         enddo
                         write(*,*) iter,fac
                         stop
                         if (fac < 1.d-4) popsConverged = .true.
                         if (iter == maxIter) then
                            popsConverged = .true.
                            call writeWarning("Maximum number of iterations reached in pop solver")
                         endif
                      enddo
                   endif
                endif
             enddo
          end do



          maxFracChange = -1.d30
          call swapPops(grid%octreeRoot, maxFracChange)
          write(*,*) "Maximum fractional change this iteration",maxFracChange
          write(*,*) "Fractional change",maxFracChange,"tolerance",tolerance , &
               "fixed rays",fixedrays,"nray",nray

          if (maxFracChange < tolerance) then
             gridConverged = .true.
          endif


          deallocate(ds, phi, i0, tau, sourceNumber, cosTheta, hitPhotosphere, weight)

          if (.not.gridConverged) then
             if (.not.fixedRays) nRay = nRay * 2
          endif
          if (nRay > maxRay) then
             nRay = maxRay
             call writeWarning("Maximum number of rays exceeded - capping")
          endif
       enddo
    enddo
    write(*,*) "ATOM loop done."
    stop
  end subroutine atomLoop

  recursive  subroutine  swapPops(thisOctal, maxFracChange)
    type(GRIDTYPE) :: grid
    type(MODELATOM) :: thisAtom
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, iUpper, iLower
    real(double) :: maxFracChange, temp
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call swapPops(child, maxFracChange)
                exit
             end if
          end do
       else
          where(thisOctal%atomLevel(subcell,:,:) /= 0.d0)
             temp = MAXVAL(abs((thisOctal%newatomLevel(subcell,:,1:6) - &
                  thisOctal%atomLevel(subcell,:,1:6)) / &
                  thisOctal%atomLevel(subcell,:,1:6)))
          end where
          if (temp > maxFracChange) then
             maxFracChange = temp
          endif
          thisOctal%atomLevel(subcell,:,:) = &
               thisOctal%newatomLevel(subcell,:,:)
       endif
    enddo
  end subroutine swapPops

  recursive subroutine  allocateLevels(grid, thisOctal, nAtom, thisAtom, nRBBTrans)
    type(GRIDTYPE) :: grid
    type(MODELATOM) :: thisAtom(:)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, iUpper, iLower
    real(double) :: etaLine
    integer :: nAtom
    integer :: nRBBTrans
    integer :: iatom
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call allocateLevels(grid, child, nAtom, thisAtom, nRBBTrans)
                exit
             end if
          end do
       else

          if (.not.associated(thisOctal%atomLevel)) then
             allocate(thisOctal%atomLevel(1:thisOctal%maxChildren, 1:nAtom, 1:maxval(thisAtom(1:nAtom)%nLevels)))
          endif
          thisOctal%atomLevel(subcell,:,:) = 1.d-30

          do iatom = 1, nAtom
             thisOctal%atomLevel(subcell, iAtom, thisAtom(iAtom)%nLevels) = thisOctal%rho(subcell)/mHydrogen
          enddo
          if (.not.associated(thisOctal%newatomLevel)) then
             allocate(thisOctal%newatomLevel(1:thisOctal%maxChildren, 1:nAtom, 1:maxval(thisAtom(1:nAtom)%nLevels)))
          endif
          thisOctal%newatomLevel = 1.d-30

          if (.not.associated(thisOctal%jnu)) then
             allocate(thisOctal%jnu(1:thisOctal%maxChildren, 1:nRBBTrans))
          endif
          thisOctal%jnu = 1.d-30
       endif
    enddo
  end subroutine allocateLevels
  

  subroutine randomRayDirection(probTowardsSource, point, source, nSource, direction, weight)
    type(OCTALVECTOR) :: point, direction, toStar
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    real(double) :: probTowardsSource, weight
    integer :: i
    real(double) :: chanceSource
    real(double) :: omegaSubtendedBySource
    real(double) :: theta, dOmega, r, cosTheta
    logical :: hitCore

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
       theta = asin(source(i)%radius/modulus(point-source(i)%position))
       cosTheta = cos(theta)
       hitCore = .false.
       do while(.not.hitCore)
          direction = randomUnitVector()
          if (acos(direction.dot.toStar) < theta) hitCore = .true.
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
end module cmf_mod
