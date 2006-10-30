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

  subroutine solveLevels(nPops, jnuLine,  temperature, thisAtom, ne, nh, jnuCont, freq, nfreq)
    integer :: nFreq
    real(double) :: freq(:), jnuCont(:)
    real(double) :: nPops(:)
    real(double) :: temperature, ne
    real(double) :: jnuLine(:)
    type(MODELATOM) :: thisAtom
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

    nLevels = thisAtom%nLevels

    allocate(matrixA(1:nLevels+2,1:nLevels+2))
    allocate(matrixB(1:nLevels+2))

    matrixA = 1.d-30
    matrixB = 0.d0


    matrixA(nLevels+2,1:nLevels+1) = 1.d0
    matrixB(nLevels+2) = nh

! RHS of equation

    do iTrans = 1, thisAtom%nTrans

       k = thisAtom%iUpper(iTrans)
       l = thisAtom%iLower(iTrans)

       Nstar = boltzSaha(l, ne, temperature)

       if (thisAtom%transType(iTrans) == "RBF") then ! radiative recomb
          recombratekl = 0.d0
          do i = 2, nFreq
             xSection = photoCrossSection(thisAtom, l, freq(i))
             recombRatekl = recombRatekl + &
                  (fourPi/(hCgs*freq(i)))*xSection*((2.d0*hCgs*freq(i)**3)/cSpeed**2 + jnuCont(i)) * &
                  exp(-(hCgs*freq(i))/(kErg*temperature))*(freq(i)-freq(i-1))
          enddo
          matrixB(l) = matrixB(l) - Nstar * recombratekl
       endif
       if (thisAtom%transType(iTrans) == "CBF") then ! collisional recomb

          collEx = collisionRate(thisAtom, iTrans, temperature)
          matrixB(l) = matrixB(l) - Nstar*collEx*ne
       endif
    enddo

    matrixB(nLevels+1) = -SUM(matrixB(1:nLevels))

! LHS

    do iTrans = 1, thisAtom%nTrans

       k = thisAtom%iUpper(iTrans)
       l = thisAtom%iLower(iTrans)
       if (thisAtom%transType(iTrans) == "RBB") then ! radiative BB rates
          iJnu = thisAtom%indexRBB(iTrans)
          call returnEinsteinCoeffs(thisAtom, iTrans, a, Bul, Blu)

          matrixA(k,k) = matrixA(k,k) - Bul * jnuLine(iJnu) - a
          matrixA(l,l) = matrixA(l,l) - Blu * jnuLine(iJnu)
          matrixA(k,l) = matrixA(k,l) + Blu * jnuLine(iJnu)
          matrixA(l,k) = matrixA(l,k) + Bul * jnuLine(iJnu) + a
       endif

       if (thisAtom%transType(iTrans) == "CBB") then ! collision BB  rates


          collEx = collisionRate(thisAtom, iTrans, temperature)
          boltzFac =  exp(-abs(thisAtom%energy(k)-thisAtom%energy(l)) / (kev*temperature))
          colldeEx = collEx / (boltzFac * thisAtom%g(k)/thisAtom%g(l))

          matrixA(k,k) = matrixA(k,k) - colldeex * ne
          matrixA(l,l) = matrixA(l,l) - collex * ne
          matrixA(k,l) = matrixA(k,l) + collex * ne
          matrixA(l,k) = matrixA(l,k) + colldeex * ne

       endif

! now do bound-free rates


       if (thisAtom%transType(iTrans) == "RBF") then ! photoionization


          l = thisAtom%iLower(iTrans)
          k = thisAtom%iUpper(iTrans)

          photoRatelk = 0.d0
          do i = 2, nFreq
             xSection = photoCrossSection(thisAtom, l, freq(i))
             photoRatelk = photoRatelk + (jnuCont(i)/(hCgs*freq(i)))*xSection*(freq(i)-freq(i-1))
!             write(*,*) i, jnucont(i),freq(i),xsection,freq(i)-freq(i-1)
          enddo
          photoRatelk = photoRatelk * fourPi
             
          recombratekl = 0.d0
          do i = 2, nFreq
             xSection = photoCrossSection(thisAtom, l, freq(i))
             recombRatekl = recombRatekl + &
                  (fourPi/(hCgs*freq(i)))*xSection*((2.d0*hCgs*freq(i)**3)/cSpeed**2 + jnuCont(i)) * &
                  exp(-(hCgs*freq(i))/(kErg*temperature))*(freq(i)-freq(i-1))
          enddo

          matrixA(l,l) = matrixA(l,l) - photoRatelk
          matrixA(k,l) = matrixA(k,l) + photoRatelk

       endif


       if (thisAtom%transType(iTrans) == "CBF") then ! collisional ionization

          k = thisAtom%iUpper(iTrans)
          l = thisAtom%iLower(iTrans)


          collEx = collisionRate(thisAtom, iTrans, temperature)

          matrixA(l,l) = matrixA(l,l) - collex * ne
          matrixA(k,l) = matrixA(k,l) + collex * ne
       endif


    enddo

!    do i = 1, nLevels+2
!       write(*,'(1p,18e10.1)') matrixA(i,1:17)
!    enddo

    call luSlv(matrixA, matrixB, nLevels+2)


    nPops(1:nLevels+1) = matrixB(1:nLevels+1)

    deallocate(matrixA, matrixB)


  end subroutine solveLevels

  subroutine getRay(grid, fromOctal, fromSubcell, position, direction, ds, phi, i0, thisAtom, source, nSource, &
       hitPhotosphere, sourceNumber, cosTheta, weight)
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    type(MODELATOM) :: thisAtom
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
    
    allocate(tau(1:thisAtom%nTrans))
    
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

          do iTrans = 1, thisAtom%nTrans

             if (thisAtom%transType(iTrans) == "RBB") then

                call returnEinsteinCoeffs(thisAtom, iTrans, a, Bul, Blu)

                alphanu = (hCgs*thisAtom%transFreq(iTrans)/fourPi) * &
                     phiProf(dv, thisOctal%microturb(subcell))/thisAtom%transFreq(iTrans)
                
                iUpper = thisAtom%iUpper(iTrans)
                iLower = thisAtom%iLower(iTrans)
                
                nLower = thisOctal%molecularLevel(subcell,iLower)
                nUpper = thisOctal%molecularLevel(subcell,iUpper)
                
                alphanu = alphanu * (nLower * Blu - nUpper * Bul)
                
                dTau = alphaNu *  (distArray(i)-distArray(i-1)) * 1.d10
                
                etaLine = hCgs * a * thisAtom%transFreq(iTrans)
                etaLine = etaLine * thisOctal%molecularLevel(subcell, iUpper)
                jnu = (etaLine/fourPi) * phiProf(dv, thisOctal%microturb(subcell))/thisAtom%transFreq(iTrans)
                
                if (alphanu /= 0.d0) then
                   snu = jnu/alphanu
                else
                   snu = tiny(snu)
                endif
                
                i0(iTrans) = i0(iTrans) +  exp(-tau(iTrans)) * (1.d0-exp(-dtau))*snu
                tau(iTrans) = tau(iTrans) + dtau
                

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


  subroutine atomLoop(grid, thisAtom, nSource, source)

    use input_variables, only : blockhandout
 
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    type(GRIDTYPE) :: grid
    type(MODELATOM) :: thisAtom
    type(OCTALVECTOR) :: position, direction
    integer :: nOctal, iOctal, subcell
    real(double), allocatable :: ds(:), phi(:), i0(:,:), tau(:)
    integer :: nRay
    type(OCTAL), pointer :: thisOctal
    integer, parameter :: maxIter = 100, maxRay = 100000
    logical :: popsConverged, gridConverged 
    integer :: iRay, iTrans, iter,i 
    integer :: iStage
    real(double), allocatable :: oldpops(:)
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
    real(double) :: nuStart, nuEnd

    nuStart = cSpeed / (100000.d0 * 1.d-8)
    nuEnd = cSpeed / (50.d0 * 1.d-8)

    nfreq = maxFreq
    do i = 1, nFreq
       freq(i) = log10(nuStart) + (dble(i-1)/dble(nFreq-1))*(log10(nuEnd)-log10(nuStart))
    enddo
    freq(1:nFreq) = 10.d0**freq(1:nFreq)

     blockHandout = .false. 

    call allocateMolecularLevels(grid, grid%octreeRoot, thisAtom)


!    call solveAllPops(grid, grid%octreeRoot, thisAtom, jnuCont, freq, nfreq)


    allocate(oldPops(1:thisAtom%nLevels))

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
          allocate(i0(1:thisAtom%nTrans, 1:nRay))
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
                              ds(iRay), phi(iRay), i0(1:thisAtom%nTrans,iRay), &
                              thisAtom, source, nSource, hitPhotosphere(iRay), sourceNumber(iray), &
                              cosTheta(iRay), weight(iRay))
                         if (hitPhotosphere(iray)) nHit = nHit + 1
                      enddo
                      write(*,*) "fraction of rays hitting core",real(nhit)/real(nray)
                      iter = 0
                      popsConverged = .false.
                      thisOctal%newMolecularLevel(subcell,:) = thisOctal%molecularLevel(subcell,:)
                      do while (.not.popsConverged)
                         iter = iter + 1
                         oldpops = thisOctal%newmolecularLevel(subcell,1:thisAtom%nLevels)
                         
                         call calculateJbarCont(nray, source, nSource, hitPhotosphere, sourceNumber, &
                              freq, nfreq, jnuCont, cosTheta, weight)
                         
                         do iRBB = 1, thisAtom%nRBBTrans
                            iTrans = thisAtom%indexRBB(iRBB)
                            call calculateJbar(thisOctal, subcell, thisAtom, nRay, ds(1:nRay), &
                                 phi(1:nRay), i0(iTrans,1:nRay), iTrans, thisOctal%jnu(subcell,iTrans), &
                                 thisOctal%newMolecularLevel(subcell,1:thisAtom%nLevels), freq, nFreq, jnuCont, weight(1:nRay))
                            
                         enddo
                         call solveLevels(thisOctal%newMolecularLevel(subcell,1:thisAtom%nLevels+1), &
                              thisOctal%jnu(subcell,1:thisAtom%nRBBTrans),  &
                              dble(thisOctal%temperature(subcell)), thisAtom, thisOctal%ne(subcell), &
                              thisOctal%rho(subcell)/mHydrogen,&
                              jnuCont, freq, nfreq)
                         thisOctal%ne(subcell) = thisOctal%newmolecularLevel(subcell,thisAtom%nLevels+1)
!                         write(*,*) "iter",iter
!                         do i = 1, thisAtom%nLevels
!                            write(*,*) i, thisOctal%newMolecularLevel(subcell,i), &
!                                 thisOctal%newMolecularLevel(subcell,i ) / &
!                                 boltzSaha(i,thisOctal%ne(subcell), dble(thisOctal%temperature(subcell)))
!                         enddo
!                         write(*,*) "nk",thisOctal%newmolecularLevel(subcell,thisAtom%nlevels+1)
!                         write(*,*) "ne",thisOctal%ne(subcell)
!
                         fac = maxval(abs((thisOctal%newMolecularLevel(subcell,1:thisAtom%nLevels) - oldpops)/oldpops))
                         write(*,*) iter,fac
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

  recursive subroutine  solveAllPops(grid, thisOctal, thisAtom, jnu, freq, nfreq)
    type(GRIDTYPE) :: grid
    type(MODELATOM) :: thisAtom
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real(double) :: jnu(:), freq(:)
    integer :: nfreq
    integer :: subcell, i, iUpper, iLower
    real(double) :: etaLine
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call solveAllpops(grid, child, thisAtom, jnu, freq, nfreq)
                exit
             end if
          end do
       else
          call solveLevels(thisOctal%molecularLevel(subcell,1:thisAtom%nLevels), &
               thisOctal%jnu(subcell,1:thisAtom%nTrans),  &
               dble(thisOctal%temperature(subcell)), thisAtom, thisOctal%ne(subcell), &
               thisOctal%rho(subcell)/mHydrogen, jnu, freq, &
               nFreq)
       endif
    enddo
  end subroutine solveAllPops

  recursive subroutine  swapPops(thisOctal, maxFracChange)
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
          where(thisOctal%molecularLevel(subcell,:) /= 0.d0)
             temp = MAXVAL(abs((thisOctal%newMolecularLevel(subcell,1:6) - &
                  thisOctal%molecularLevel(subcell,1:6)) / &
                  thisOctal%molecularLevel(subcell,1:6)))
          end where
          if (temp > maxFracChange) then
             maxFracChange = temp
          endif
          thisOctal%molecularLevel(subcell,:) = &
               thisOctal%newmolecularLevel(subcell,:)
       endif
    enddo
  end subroutine swapPops

  recursive subroutine  allocateMolecularLevels(grid, thisOctal, thisAtom)
    type(GRIDTYPE) :: grid
    type(MODELATOM) :: thisAtom
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, iUpper, iLower
    real(double) :: etaLine
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call allocateMolecularLevels(grid, child, thisAtom)
                exit
             end if
          end do
       else

          if (.not.associated(thisOctal%molecularLevel)) then
             allocate(thisOctal%molecularLevel(1:thisOctal%maxChildren, 1:thisAtom%nLevels+1))
          endif
          thisOctal%molecularLevel = 1.d-30
!          do i = 1, thisAtom%nLevels
!             thisOctal%molecularLevel(subcell, i) = boltzSaha(i, thisOctal%ne(subcell), &
!                  dble(thisOctal%temperature(subcell)))
!          enddo

          thisOctal%molecularLevel(subcell, thisAtom%nLevels+1) = thisOctal%rho(subcell)/mHydrogen

          if (.not.associated(thisOctal%newmolecularLevel)) then
             allocate(thisOctal%newmolecularLevel(1:thisOctal%maxChildren, 1:thisAtom%nLevels+1))
          endif
          thisOctal%newmolecularLevel = 1.d-30

          if (.not.associated(thisOctal%jnu)) then
             allocate(thisOctal%jnu(1:thisOctal%maxChildren, 1:thisAtom%nRBBTrans))
          endif
!          do i = 1, thisAtom%nTrans
!             thisOctal%jnu(subcell,i) = bnu(thisAtom%transFreq(i), dble(thisOctal%temperature(subcell)))
!         enddo
          
          thisOctal%jnu = 1.d-30
       endif
    enddo
  end subroutine allocateMolecularLevels


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
