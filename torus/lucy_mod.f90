module lucy_mod

  use constants_mod
  use grid_mod
  use vector_mod
  use phasematrix_mod
  use amr_mod
  use source_mod
  use spectrum_mod
  use timing
  implicit none


contains


  subroutine lucyRadiativeEquilibrium(grid, miePhase, nMuMie, nLambda, lamArray, temperature, nLucy)

    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec, uHat, uNew
    integer :: nlucy
    integer :: nLambda, nMuMie
    type(PHASEMATRIX):: miePhase(1:nLambda, 1:nMuMie)
    real :: lamArray(:)
    real(kind=doubleKind), allocatable :: distanceGrid(:, :, :)
    real(kind=doubleKind) :: r
    integer, parameter :: nFreq = 100
    integer :: i, j, k
    real(kind=doubleKind) :: freq(nFreq), dnu(nFreq), probDistPlanck(nFreq), probDistJnu(nFreq)
    real(kind=doubleKind) :: temperature
    integer :: nMonte, iMonte, nScat, nAbs
    real(kind=doubleKind) :: thisFreq,  logFreqStart, logFreqEnd
    real(kind=doubleKind) :: albedo
    logical :: escaped
    integer :: i1, i2, i3
    real :: t1, t2, t3
    real(kind=doubleKind) :: thisLam
    integer :: iLam
    integer :: nInf
    real(kind=doubleKind) :: t
    integer :: nt
    real(kind=doubleKind) :: ang
    integer :: iIter, nIter = 5
    real(kind=doubleKind) :: kabs
    real(kind=doubleKind) ::  kappaP, norm
    real(kind=doubleKind) :: adot, V, epsOverDeltaT
    real(kind=doubleKind) :: newT, deltaT
    real(kind=doubleKind) :: meanDeltaT
    real(kind=doubleKind) :: dx, dy, tr(6), fg, bg
    real(kind=doubleKind), allocatable :: tempImage(:,:)
    integer :: nDT
    real(kind=doubleKind) :: totalEmission

    allocate(tempImage(1:grid%nx,1:grid%ny))


    logFreqStart = log10((cSpeed / (lamArray(nLambda)*1.e-8)))
    logFreqEnd =  log10((cSpeed / (lamArray(1)*1.e-8)))
    write(*,*) logFreqStart, logFreqEnd
    do i = 1, nFreq
       freq(i) = logFreqStart + (dble(i-1)/dble(nFreq-1))*(logFreqEnd-logFreqStart)
       freq(i) = 10.**freq(i)
    enddo
    write(*,*) "Lam",(cSpeed/freq(1))*1.e8,(cSpeed/freq(nFreq))*1.e8

    do i = 2, nFreq-1
       dnu(i) = 0.5*((freq(i+1)+freq(i))-(freq(i)+freq(i-1)))
    enddo
    dnu(1) = freq(2)-freq(1)
    dnu(nFreq) = freq(nFreq)-freq(nFreq-1)



    call setupFreqProb(temperature, freq, dnu, nFreq, ProbDistPlanck)


    write(*,'(a)') "Computing lucy radiative equilibrium..."


    allocate(distanceGrid(1:grid%nx, 1:grid%ny, 1:grid%nz))


    do iIter = 1, nIter
       distanceGrid = 0.d0
       write(*,*) "Iteration",iIter
       nInf = 0
       nScat = 0
       nAbs = 0

       nMonte = nLucy

       do iMonte = 1, nMonte
          escaped = .false.
          call random_number(r)
          call locate(probDistPlanck, nFreq, r, j)
          thisFreq = freq(j) + (freq(j+1) - freq(j))* &
               (r - probDistPlanck(j))/(probDistPlanck(j+1)-probDistPlanck(j))


          rVec  = grid%rCore * randomUnitVector()
          uHat = fromPhotosphereVector(rVec)
          if ( (rVec .dot. uHat) < 0.) uHat = (-1.) * uHat

          do while(.not.escaped)
             call toNextEvent(grid, rVec, uHat, escaped, distanceGrid, thisFreq, nLambda, lamArray)

             if (escaped) nInf = nInf + 1

             if (.not. escaped) then

                thisLam = (cSpeed / thisFreq) * 1.e8
                call locate(lamArray, nLambda, real(thisLam), iLam)

                call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
                if (.not.grid%oneKappa) then
                   albedo = grid%kappaSca(i1,i2,i3,iLam)/(grid%kappaSca(i1,i2,i3,iLam)+grid%kappaAbs(i1,i2,i3,iLam))
                else
                   albedo = grid%oneKappaSca(iLam) / (grid%oneKappaSca(iLam)+grid%oneKappaAbs(iLam))
                endif
                call random_number(r)
                if (r < albedo) then

                   uNew = newDirectionMie(uhat, real(thisLam), lamArray, nLambda, miePhase, nMuMie)

                   nScat = nScat + 1
                   uHat = uNew

                else
                   nAbs = nAbs + 1
                   probDistJnu(1) = 0.
                   do i = 2, nFreq
                      thisLam = (cSpeed / freq(i)) * 1.e8
                      call hunt(lamArray, nLambda, real(thisLam), iLam)
                      if ((ilam >=1).and.(ilam <= nlambda)) then
                         if (.not.grid%oneKappa) then
                            kabs = grid%kappaAbs(i1,i2,i3,iLam)
                         else
                            kabs = grid%oneKappaAbs(iLam) * grid%rho(i1,i2,i3)
                         endif
                         probDistJnu(i) = probDistJnu(i-1) + &
                             (bnu(freq(i),dble(grid%temperature(i1,i2,i3)))) &
                                                     * kabs * dnu(i)
                      endif
                   enddo
                   probDistJnu(1:nFreq) = probDistJnu(1:nFreq) / probDistJnu(nFreq)
                   call random_number(r)
                   call locate(probDistJnu, nFreq, r, j)
                   if (j == nFreq) j = nFreq -1
                   thisFreq = freq(j) + (freq(j+1) - freq(j))* &
                        (r - probDistJnu(j))/(probDistJnu(j+1)-probDistJnu(j))
                   uHat = randomUnitVector()
                endif

             endif
          enddo
       enddo
       write(*,'(a,f7.2)') "Photons done.",real(ninf)/real(nmonte)
       write(*,'(a,f7.3)') "Mean number of scatters per photon: ",real(nScat)/real(nMonte)
       write(*,'(a,f7.3)') "Mean number of absorbs  per photon: ",real(nAbs)/real(nMonte)


       V = (dble(grid%xAxis(2)-grid%xAxis(1))*dble(grid%yAxis(2)-grid%yAxis(1))*dble(grid%zAxis(2)-grid%zAxis(1)))


       epsOverDeltaT = (dble(grid%lCore)) / dble(nInf)

       meanDeltaT = 0.
       nDT = 0


       totalEmission = 0.
       do i1 = 1, grid%nx
          do i2 = 1, grid%ny
             do i3 = 1, grid%nz
                if (grid%inUse(i1,i2,i3)) then
                   adot = epsoverDeltaT * (1.d0 / v) * distancegrid(i1,i2,i3) / 1.d30
                   kappaP = 0.d0
                   norm = 0.d0
                   do i = 1, nFreq
                      thisLam = (cSpeed / freq(i)) * 1.e8
                      call hunt(lamArray, nLambda, real(thisLam), iLam)
                      if ((iLam >=1) .and. (iLam <= nLambda)) then
                         if (.not.grid%oneKappa) then
                            kabs = grid%kappaAbs(i1,i2,i3,iLam)
                         else
                            kabs = grid%oneKappaAbs(iLam) * grid%rho(i1,i2,i3)
                         endif

                         kappaP = kappaP + kabs * &
                              bnu(freq(i),dble(grid%temperature(i1,i2,i3)))  * dnu(i)
                         norm = norm + bnu(freq(i),dble(grid%temperature(i1,i2,i3)))  * dnu(i)
                      endif
                   enddo
                   kappaP = kappaP / norm /1.d10


                   newT = (pi / stefanBoltz) * aDot / (fourPi * kappaP)
                   newT = newT**0.25

                   deltaT = newT - grid%temperature(i1,i2,i3)
                   grid%temperature(i1,i2,i3) = max(1.e-3,grid%temperature(i1,i2,i3) + 0.8 * real(deltaT))
                   nDT = nDT  + 1
                   meanDeltaT = meanDeltaT + deltaT
                   kappaP = 0.d0
                   norm = 0.d0
                   do i = 1, nFreq
                      thisLam = (cSpeed / freq(i)) * 1.e8
                      call hunt(lamArray, nLambda, real(thisLam), iLam)
                      if ((iLam >=1) .and. (iLam <= nLambda)) then
                         if (.not.grid%oneKappa) then
                            kabs = grid%kappaAbs(i1,i2,i3,iLam)
                         else
                            kabs = grid%oneKappaAbs(iLam) * grid%rho(i1,i2,i3)
                         endif
                         kappaP = kappaP + kabs * &
                              bnu(freq(i),dble(grid%temperature(i1,i2,i3))) * dnu(i)
                         norm = norm + bnu(freq(i),dble(grid%temperature(i1,i2,i3))) * dnu(i)
                      endif
                   enddo
                   kappaP = kappaP / norm /1.e10
                   grid%etaCont(i1,i2,i3) = fourPi * kappaP * (stefanBoltz/pi) * (grid%temperature(i1,i2,i3)**4)
                   totalEmission = totalEmission + grid%etaCont(i1,i2,i3) * V

                endif
             enddo
          enddo
       enddo

       write(*,*) meanDeltaT / real(nDT)

       write(*,*) "Emissivity of dust / core", totalEmission /grid%lCore * 1.e30

       do i = 1, grid%nx
          do j = 1, grid%ny
             tempImage(i,j) = grid%temperature(i,j,grid%nz/2)
          enddo
       enddo

       dx = grid%xAxis(2) - grid%xAxis(1)
       dy = grid%yAxis(2) - grid%yAxis(1) 
       tr(1) =  - dx + grid%xAxis(1)
       tr(2) = dx
       tr(3) = 0.
       tr(4) =  - dy + grid%yAxis(1)
       tr(5) = 0.
       tr(6) = dy

       bg = MAXVAL(tempImage, mask=grid%inUse(1:grid%nx,1:grid%ny,grid%nz/2))
       fg = MINVAL(tempImage, mask=grid%inUse(1:grid%nx,1:grid%ny,grid%nz/2))
       call pgbegin(0,"/xs",1,1)

       call pgvport(0.1,0.9,0.1,0.9)

       call pgwnad(grid%xAxis(1), grid%xAxis(grid%nx), grid%yAxis(1), grid%yAxis(grid%ny))
       call palette(3)
       call pgimag(tempImage, grid%nx, grid%ny, 1, grid%nx, 1, grid%ny, fg, bg, tr)
       call pgwedg("RI", 0.,10., fg,bg," ")
       call pgbox('bcnst',0,0,'bcnst',0,0)
       call pgend

    enddo
    open(21,file="r.dat",status="unknown",form="formatted")
    do i = 1, 100
       r = grid%xAxis(grid%nx) * real(i-1)/99.
       t = 0
       nT = 0
       do j = 1, 100
          ang = twoPi * real(j-1)/100.
          rVec = VECTOR(r*cos(ang), r*sin(ang),0.)
          call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
          if (grid%inUse(i1,i2,i3)) then
             nT = nT + 1
             t = t + grid%temperature(i1,i2,i3)
          endif
       enddo
       if (nt > 0) then
          t = t / real(nt)
          write(21,*) r / grid%rCore, t
       endif
    enddo
    close(21)

    open(22,file="temps.dat",status="unknown",form="unformatted")
    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             write(22) grid%temperature(i,j,k), grid%etaCont(i,j,k)
          enddo
       enddo
    enddo
    close(22)


  end subroutine lucyRadiativeEquilibrium

  subroutine lucyRadiativeEquilibriumAMR(grid, miePhase, nMuMie, nLambda, lamArray, temperature, &
       source, nSource, nLucy)

    type(GRIDTYPE) :: grid
    type(SOURCETYPE) :: source(*), thisSource
    integer :: nSource
    integer :: nLucy
    type(VECTOR) ::  uHat, uNew, rVec, rHat
    type(OCTALVECTOR) :: octVec
    type(OCTAL), pointer :: thisOctal
    integer :: nLambda, nMuMie
    type(PHASEMATRIX):: miePhase(1:nLambda, 1:nMuMie)
    real  :: lamArray(:)
    real(kind=doubleKind) :: r
    integer :: nFreq
    integer :: i, j
    real(kind=doubleKind) :: freq(2000), dnu(2000), probDistJnu(2000)
    real(kind=doubleKind) :: temperature
!    real(kind=doubleKind) :: probDistPlanck(nFreq)
    real :: kappaScaReal, kappaAbsReal
    integer :: nMonte, iMonte, nScat, nAbs
    real(kind=doubleKind) :: thisFreq,  logFreqStart, logFreqEnd
    real(kind=doubleKind) :: albedo
    logical :: escaped
    real(kind=doubleKind) :: t1
    real(kind=doubleKind) :: thisLam,wavelength
    integer :: iSource
    integer :: nRemoved
    real :: Tthresh
    integer :: iLam
    integer :: nInf
    real(kind=doubleKind) :: ang
    integer :: iIter, nIter = 5
    integer :: nt
    real(kind=doubleKind) ::   epsOverDeltaT
    real(kind=doubleKind) :: meanDeltaT, meant
    integer :: nDT, nUndersampled
    real(kind=doubleKind) :: totalEmission
    integer :: subcell
    real :: treal
    real(kind=doubleKind) :: lCore
    real(kind=doubleKind) :: kabs


    nMonte = nLucy

    nFreq = nLambda
    do i = 1, nFreq
       freq(nFreq-i+1) = cSpeed / (lamArray(i)*1.e-8)
    enddo


    do i = 2, nFreq-1
       dnu(i) = 0.5*((freq(i+1)+freq(i))-(freq(i)+freq(i-1)))
    enddo
    dnu(1) = freq(2)-freq(1)
    dnu(nFreq) = freq(nFreq)-freq(nFreq-1)



!    call setupFreqProb(temperature, freq, dnu, nFreq, ProbDistPlanck)


    write(*,'(a)') "Computing lucy radiative equilibrium in AMR..."


    lCore = 0.d0
    do i = 1, nSource
       lCore = lCore + source(i)%luminosity
    enddo
    write(*,'(a,1pe12.5)') "Total souce luminosity (lsol): ",lCore/lSol

    nRemoved = 1
    do while (nRemoved > 0)
    do iIter = 1, nIter

       call tune(6, "One Lucy Rad Eq Itr")  ! start a stopwatch
       
       call zeroDistanceGrid(grid%octreeRoot)
       write(*,*) "Iteration",iIter,",",nmonte," photons"
       nInf = 0
       nScat = 0
       nAbs = 0
       do iMonte = 1, nMonte

          call randomSource(source, nSource, iSource)
          thisSource = source(iSource)
           
          call getPhotonPositionDirection(thisSource, rVec, uHat, rHat)
          escaped = .false.
          call getWavelength(thisSource%spectrum, wavelength)
          thisFreq = cSpeed/(wavelength / 1.e8)

!          call random_number(r)
!          call locate(probDistPlanck, nFreq, r, j)
!          thisFreq = freq(j) + (freq(j+1) - freq(j))* &
!               (r - probDistPlanck(j))/(probDistPlanck(j+1)-probDistPlanck(j))
!          rVec  = grid%rCore * randomUnitVector()
!          uHat = fromPhotosphereVector(rVec)



          if ( (rVec .dot. uHat) < 0.) uHat = (-1.) * uHat

          do while(.not.escaped)
             call toNextEventAMR(grid, rVec, uHat, escaped, thisFreq, nLambda, lamArray)

             if (escaped) nInf = nInf + 1

             if (.not. escaped) then

                thisLam = (cSpeed / thisFreq) * 1.e8
                call locate(lamArray, nLambda, real(thisLam), iLam)
                octVec = rVec 
                call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, iLambda=iLam, &
                     kappaSca=kappaScaReal, kappaAbs=kappaAbsReal, grid=grid)
                albedo = kappaScaReal / (kappaScaReal + kappaAbsReal)

                call random_number(r)
                if (r < albedo) then

                   uNew = newDirectionMie(uhat, real(thisLam), lamArray, nLambda, miePhase, nMuMie)
                   nScat = nScat + 1
                   uHat = uNew

                else
                   call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, &
                        foundSubcell=subcell, temperature=treal, grid=grid)
                   t1 = dble(treal)
                   nAbs = nAbs + 1
                   probDistJnu(1) = 0.
                   do i = 2, nFreq
                      thisLam = (cSpeed / freq(i)) * 1.e8
                      call hunt(lamArray, nLambda, real(thisLam), iLam)
                      if ((ilam >=1).and.(ilam <= nlambda)) then
                         if (.not.grid%oneKappa) then
                            kabs = thisOctal%kappaAbs(subcell,ilam)
                         else
                            kabs = grid%oneKappaAbs(iLam) * thisOctal%rho(subcell)
                         endif
                         probDistJnu(i) = probDistJnu(i-1) + &
                              bnu(freq(i),t1) * kabs * dnu(i)
                      endif
                   enddo
                   if (probDistJnu(nFreq) /= 0.) then
                      probDistJnu(1:nFreq) = probDistJnu(1:nFreq) / probDistJnu(nFreq)
                      call random_number(r)
                      call locate(probDistJnu, nFreq, r, j)
                      if (j == nFreq) j = nFreq -1
                      thisFreq = freq(j) + (freq(j+1) - freq(j))* &
                           (r - probDistJnu(j))/(probDistJnu(j+1)-probDistJnu(j))
                   endif
                   uHat = randomUnitVector()
                endif

             endif
          enddo
       enddo
       write(*,'(a,f7.2)') "Photons done.",real(ninf)/real(nmonte)
       write(*,'(a,f7.3)') "Mean number of scatters per photon: ",real(nScat)/real(nMonte)
       write(*,'(a,f7.3)') "Mean number of absorbs  per photon: ",real(nAbs)/real(nMonte)



       epsOverDeltaT = (lCore) / dble(nInf)

       meanDeltaT = 0.
       nDT = 0
       nUndersampled = 0


       totalEmission = 0.
       call calculateTemperatureCorrections(grid%octreeRoot, totalEmission, epsOverDeltaT, &
       nFreq, freq, dnu, lamarray, nLambda, grid, nDt, nUndersampled)

       write(*,'(a,i5)') "Number of undersampled cells: ",nUndersampled
       write(*,'(a,f8.2)') "Percentage of undersampled cells: ",100.*real(nUndersampled)/real(nDt)
       write(*,*) "Emissivity of dust / core", totalEmission /lCore * 1.e30
       write(*,*) "Total emission",totalEmission

       call tune(6, "One Lucy Rad Eq Itr")  ! stop a stopwatch
       
    enddo

    Tthresh = 2000.

    write(*,*) "Removing dust with T > Tthresh: ",Tthresh

    nRemoved = 0
    call removeDust(grid%octreeRoot, Tthresh, nRemoved)

    write(*,*) "Number of cells removed: ",nRemoved
  enddo


  end subroutine lucyRadiativeEquilibriumAMR


  subroutine setupFreqProb(temperature, freq, dnu, nFreq, probDist)

    ! Lucy 1999, A&A, 344, 282 Equation 3

    real(kind=doubleKind) :: temperature
    integer :: nFreq
    real(kind=doubleKind) :: freq(*)
    real(kind=doubleKind) :: probDist(*)
    real(kind=doubleKind) :: dnu(*)

    integer :: i


    
    probDist(1:nFreq) = 0.
    do i = 2, nFreq
       probDist(i) = probDist(i-1) + bnu(dble(freq(i)),dble(temperature)) * dnu(i)

    enddo

    probDist(1:nFreq) = probDist(1:nFreq) / probDist(nFreq)

  end subroutine setupFreqProb


 subroutine toNextEvent(grid, rVec, uHat,  escaped, distanceGrid, thisFreq, nLambda, lamArray)


   type(GRIDTYPE) :: grid
   type(VECTOR) :: rVec, uHat
   integer :: i1, i2, i3
   real :: t1, t2, t3
   real(kind=doubleKind) :: tval, tau, r
   real :: lamArray(*)
   integer :: nLambda
   logical :: stillinGrid
   logical :: escaped
   real(kind=doubleKind) :: thisTau
   real(kind=doubleKind) :: kabs, ksca
   real(kind=doubleKind) :: thisFreq
   real(kind=doubleKind) :: distanceGrid(1:grid%nx, 1:grid%ny, 1:grid%nz)
   real(kind=doubleKind) :: thisLam
   integer :: iLam

    stillinGrid = .true.
    escaped = .false.


    thisLam = (cSpeed / thisFreq) * 1.e8
    call hunt(lamArray, nLambda, real(thisLam), iLam)
    if ((ilam < 1).or.(ilam > nlambda)) then
       write(*,*) "ilam errro",ilam
    endif

    call random_number(r)
    tau = -log(r)
    call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
    call intersectCube(grid, rVec, i1, i2, i3, uHat, tVal)

    if (.not.grid%oneKappa) then
       kabs = grid%kappaAbs(i1,i2,i3,iLam)
       ksca = grid%kappaSca(i1,i2,i3,iLam)
    else
       kabs = grid%oneKappaAbs(iLam) * grid%rho(i1,i2,i3)
       ksca = grid%oneKappaSca(iLam) * grid%rho(i1,i2,i3)
    endif

    thisTau = tVal * (kabs + ksca)

    do while(stillinGrid .and. (tau > thisTau)) 
       rVec = rVec + tVal * uHat


       if ((rVec%x > grid%xAxis(grid%nx)).or. &
           (rVec%y > grid%yAxis(grid%ny)).or. &
           (rVec%z > grid%zAxis(grid%nz)).or. &
           (rVec%x < grid%xAxis(1)).or. &
           (rVec%y < grid%yAxis(1)).or. &
           (rVec%z < grid%zAxis(1))) then
          stillinGrid = .false.
          escaped = .true.
          
       else
          if (.not.grid%oneKappa) then
             kabs = grid%kappaAbs(i1,i2,i3,iLam)
          else
             kabs = grid%oneKappaAbs(iLam) * grid%rho(i1,i2,i3)
          endif
          distanceGrid(i1,i2,i3) = distanceGrid(i1,i2,i3) + tVal * kabs
       endif


       if (stillinGrid) then
          call random_number(r)
          tau = -log(r)
          call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
          call intersectCube(grid, rVec, i1, i2, i3, uHat, tVal)
          if (.not.grid%oneKappa) then
             kabs = grid%kappaAbs(i1,i2,i3,iLam)
             ksca = grid%kappaSca(i1,i2,i3,iLam)
          else
             kabs = grid%oneKappaAbs(iLam) * grid%rho(i1,i2,i3)
             ksca = grid%oneKappaSca(iLam) * grid%rho(i1,i2,i3)
          endif
          thisTau = tVal * (kabs+ksca)
          if (tVal == 0.) then
             escaped = .true.
             stillingrid = .false.
          endif
       endif
    enddo

    if (.not.escaped) then
       call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
       call intersectCube(grid, rVec, i1, i2, i3, uHat, tVal)

       if ((tVal*tau/thisTau) > 2.*(grid%xAxis(2)-grid%xAxis(1))) then
          write(*,*) "tval*tau/thistau too big", tval*tau/thisTau
       else
          if (.not.grid%oneKappa) then
             kabs = grid%kappaAbs(i1,i2,i3,iLam)
          else
             kabs = grid%oneKappaAbs(iLam) * grid%rho(i1,i2,i3)
          endif
          distanceGrid(i1,i2,i3) = distanceGrid(i1,i2,i3) + (tVal*tau/thisTau) * kabs
       endif
       

       if (tau > thisTau) then
          write(*,*) "tau > thistau"
          stop
       endif
       
       rVec = rVec + (tVal*tau/thisTau) * uHat
    endif




 end subroutine toNextEvent


  subroutine intersectCube(grid, posVec, i1,i2,i3,direction, tval)
   use vector_mod
   use grid_mod
   implicit none
   type(GRIDTYPE) :: grid
   type(VECTOR) :: posVec, direction, norm(6), p3(6)
   real(kind=doubleKind) :: t(6),tval,denom(6)
   integer :: i,j
   logical :: ok, thisOk(6)
   integer :: i1, i2, i3

   ok = .true.

   norm(1) = VECTOR(1., 0., 0.)
   norm(2) = VECTOR(0., 1., 0.)
   norm(3) = VECTOR(0., 0., 1.)
   norm(4) = VECTOR(-1., 0., 0.)
   norm(5) = VECTOR(0., -1., 0.)
   norm(6) = VECTOR(0., 0., -1.)

   p3(1) = VECTOR(grid%xAxis(i1+1), 0., 0.)
   p3(2) = VECTOR(0.,grid%yAxis(i2+1),0.)
   p3(3) = VECTOR(0.,0.,grid%zAxis(i3+1))
   p3(4) = VECTOR(grid%xAxis(i1), 0., 0.)
   p3(5) = VECTOR(0.,grid%yAxis(i2),0.)
   p3(6) = VECTOR(0.,0.,grid%zAxis(i3))

   thisOk = .true.

   do i = 1, 6

      denom(i) = norm(i) .dot. direction
      if (denom(i) /= 0.) then
         t(i) = (norm(i) .dot. (p3(i)-posVec))/denom(i)
      else
         thisOk(i) = .false.
         t(i) = 0.
      endif
      if (t(i) < 0.) thisOk(i) = .false.
!      if (denom > 0.) thisOK(i) = .false.
 enddo




  
  j = 0
  do i = 1, 6
    if (thisOk(i)) j=j+1
  enddo

  if (j == 0) ok = .false.
   
  if (.not.ok) then
     write(*,*) i1, i2, i3
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  tval = minval(t, mask=thisOk)
  tval = max(tval * 1.001d0,dble((grid%xAxis(2)-grid%xAxis(1))/1000.))


  if (tval == 0.) then
     write(*,*) i1, i2, i3,tval
     write(*,*) posVec
     write(*,*) grid%xAxis(i1),grid%yAxis(i2),grid%zAxis(i3)
     write(*,*) grid%xAxis(i1+1),grid%yAxis(i2+1),grid%zAxis(i3+1)
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  if (tval > 2.*(grid%xAxis(2)-grid%xAxis(1))) then
!     write(*,*) "tval too big",tval,i1,i2,i3,posvec
!     write(*,*) "direction",direction
!     write(*,*) t(1:6)
!     write(*,*) denom(1:6)
  endif


  end subroutine intersectCube 


  recursive subroutine zeroDistanceGrid(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, 8
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
          thisOctal%distanceGrid(subcell) = 0.
          thisOctal%nCrossings(subcell) = 0
       endif
    enddo
  end subroutine zeroDistanceGrid


  recursive subroutine calculateTemperatureCorrections(thisOctal, totalEmission, epsOverDeltaT, &
       nFreq, freq, dnu, lamarray, nLambda, grid, nDt, nUndersampled)

    real(kind=doubleKind) :: totalEmission
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(gridtype) :: grid
    integer :: subcell, i
    real(kind=doubleKind) :: V, epsOverDeltaT
    real(kind=doubleKind) :: adot
    real(kind=doubleKind) :: kappaP, norm
    real(kind=doubleKind) :: thisLam
    integer :: ilam, nLambda
    real(kind=doubleKind) :: newT, deltaT
    integer :: ndt
    real(kind=doubleKind) :: meanDeltaT 
    real(kind=doubleKind) :: kabs
    real :: lamArray(*), dlam(2000)
    integer :: nFreq
    real(kind=doubleKind) :: freq(*)
    real(kind=doubleKind) :: dnu(*)
    integer :: nUndersampled
    do subcell = 1, 8
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateTemperatureCorrections(child, totalEmission, epsOverDeltaT, &
                     nFreq, freq, dnu, lamarray, nLambda, grid, nDt, nUndersampled)
                exit
             end if
          end do
       else
          if (thisOctal%inFlow(subcell)) then
             v = thisOctal%subcellSize**3
             adot = epsoverDeltaT * (1.d0 / v) * thisOctal%distancegrid(subcell) / 1.d30
             kappaP = 0.d0
             norm = 0.d0
             do i = 1, nFreq
                thisLam = (cSpeed / freq(i)) * 1.e8
                call hunt(lamArray, nLambda, real(thisLam), iLam)
                if ((iLam >=1) .and. (iLam <= nLambda)) then
                   if (.not.grid%oneKappa) then
                      kabs = thisOctal%kappaAbs(subcell,iLam)
                   else
                      kabs = grid%oneKappaAbs(iLam) * thisOctal%rho(subcell)
                   endif
                   kappaP = kappaP + dble(kabs) * &
                        dble(bnu(dble(freq(i)),dble(thisOctal%temperature(subcell))))  * dble(dnu(i))
                   norm = norm + dble(bnu(dble(freq(i)),dble(thisOctal%temperature(subcell))))  * dble(dnu(i))
                endif
             enddo
             kappaP = kappaP / norm /1.d10
             
             
             if (kappaP /= 0.) then
                newT = (pi / stefanBoltz) * aDot / (fourPi * kappaP)
                newT = newT**0.25
             else
                newT = 1.e-3
             endif
             deltaT = newT - thisOctal%temperature(subcell)
             if (deltaT > 10000.) then
                write(*,*) deltaT, thisOctal%nCrossings(subcell)
                !                write(*,*) deltaT, thisOctal%temperature(subcell), newT
                !                write(*,*) adot,kappap,thisOctal%rho(subcell)
             endif
             if (thisOctal%nCrossings(subcell) .ge. 5) then
                thisOctal%temperature(subcell) = max(1.e-3,thisOctal%temperature(subcell) + real(0.8 * deltaT))
             else
                nUnderSampled = nUndersampled + 1
             endif
             kappaP = 0.d0
             norm = 0.d0
             do i = 1, nFreq
                thisLam = (cSpeed / freq(i)) * 1.e8
                call hunt(lamArray, nLambda, real(thisLam), iLam)
                if ((iLam >=1) .and. (iLam <= nLambda)) then
                   if (.not.grid%oneKappa) then
                      kabs = thisOctal%kappaAbs(subcell,iLam)
                   else
                      kabs = grid%oneKappaAbs(iLam) * thisOctal%rho(subcell)
                   endif
                   
                   kappaP = kappaP + dble(kabs) * &
                        dble(bnu(dble(freq(i)),dble(thisOctal%temperature(subcell)))) * dble(dnu(i)) 
                   norm = norm + dble(bnu(dble(freq(i)),dble(thisOctal%temperature(subcell)))) * dble(dnu(i)) 
                endif
             enddo
             kappaP = kappaP / norm / 1.e10
             thisOctal%etaCont(subcell) = fourPi * kappaP * (stefanBoltz/pi) * (thisOctal%temperature(subcell)**4)
             totalEmission = totalEmission + thisOctal%etaCont(subcell) * V
             nDT = nDT  + 1
             meanDeltaT = meanDeltaT + deltaT
          else 
             thisOctal%etaCont(subcell) = 0.
          endif
       endif
    enddo
     end subroutine calculateTemperatureCorrections

  subroutine intersectCubeAMR(grid, posVec, direction, tval)
   use vector_mod
   use grid_mod
   implicit none
   type(GRIDTYPE) :: grid
   type(VECTOR) :: posVec, direction, norm(6), p3(6)
   type(OCTAL),pointer :: thisOctal
   type(OCTALVECTOR) :: subcen, point
   integer :: subcell
   
   real(kind=doubleKind) :: t(6),tval,denom(6)
   integer :: i,j
   logical :: ok, thisOk(6)

   point = posVec

   call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
   subcen =  subcellCentre(thisOctal,subcell)
   ok = .true.

   norm(1) = VECTOR(1., 0., 0.)
   norm(2) = VECTOR(0., 1., 0.)
   norm(3) = VECTOR(0., 0., 1.)
   norm(4) = VECTOR(-1., 0., 0.)
   norm(5) = VECTOR(0., -1., 0.)
   norm(6) = VECTOR(0., 0., -1.)

   p3(1) = VECTOR(subcen%x+thisOctal%subcellsize/2., 0., 0.)
   p3(2) = VECTOR(0., subcen%y+thisOctal%subcellsize/2. ,0.)
   p3(3) = VECTOR(0.,0.,subcen%z+thisOctal%subcellsize/2.)
   p3(4) = VECTOR(subcen%x-thisOctal%subcellsize/2., 0., 0.)
   p3(5) = VECTOR(0.,subcen%y-thisOctal%subcellsize/2.,0.)
   p3(6) = VECTOR(0.,0.,subcen%z-thisOctal%subcellsize/2.)

   thisOk = .true.

   do i = 1, 6

      denom(i) = norm(i) .dot. direction
      if (denom(i) /= 0.) then
         t(i) = (norm(i) .dot. (p3(i)-posVec))/denom(i)
      else
         thisOk(i) = .false.
         t(i) = 0.
      endif
      if (t(i) < 0.) thisOk(i) = .false.
!      if (denom > 0.) thisOK(i) = .false.
 enddo




  
  j = 0
  do i = 1, 6
    if (thisOk(i)) j=j+1
  enddo

  if (j == 0) ok = .false.
   
  if (.not.ok) then
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  tval = minval(t, mask=thisOk)
  tval = max(tval * 1.001d0,dble(thisOctal%subCellSize/1000.))


  if (tval == 0.) then
     write(*,*) posVec
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  if (tval > sqrt(3.)*thisOctal%subcellsize) then
!     write(*,*) "tval too big",tval/(sqrt(3.)*thisOctal%subcellSize)
!     write(*,*) "direction",direction
!     write(*,*) t(1:6)
!     write(*,*) denom(1:6)
  endif


  end subroutine intersectCubeAMR



 subroutine toNextEventAMR(grid, rVec, uHat,  escaped,  thisFreq, nLambda, lamArray)


   type(GRIDTYPE) :: grid
   type(VECTOR) :: rVec, uHat
   type(OCTALVECTOR) :: octVec
   type(OCTAL), pointer :: thisOctal
   type(OCTAL),pointer :: oldOctal
   integer :: subcell
   real(kind=doubleKind) :: tval, tau, r
   real :: lamArray(*)
   integer :: nLambda
   logical :: stillinGrid
   logical :: escaped
   real :: kappaScaReal, kappaAbsReal
   real(kind=doubleKind) :: thisTau
   real(kind=doubleKind) :: thisFreq
   real(kind=doubleKind) :: thisLam
   integer :: iLam

    stillinGrid = .true.
    escaped = .false.


    thisLam = (cSpeed / thisFreq) * 1.e8
    call hunt(lamArray, nLambda, real(thisLam), iLam)
    if ((ilam < 1).or.(ilam > nlambda)) then
       write(*,*) "ilam errro",ilam
    endif

    call random_number(r)
    tau = -log(r)
    call intersectCubeAMR(grid, rVec, uHat, tVal)
    octVec = rVec
    call amrGridValues(grid%octreeRoot, octVec, iLambda=iLam,  foundOctal=thisOctal, &
         foundSubcell=subcell, kappaSca=kappaScaReal, kappaAbs=kappaAbsReal, grid=grid)
    oldOctal => thisOctal
    thisTau = tVal * dble(kappaAbsReal + kappaScaReal)


    do while(stillinGrid .and. (tau > thisTau)) 
       rVec = rVec + tVal * uHat

       octVec = rVec
       if (.not.inOctal(grid%octreeRoot, octVec)) then
          stillinGrid = .false.
          escaped = .true.
       else
          thisOctal%distanceGrid(subcell) = thisOctal%distanceGrid(subcell) + tVal * dble(kappaAbsReal)
          thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1
       endif


       if (stillinGrid) then
          call random_number(r)
          tau = -log(r)
          call amrGridValues(grid%octreeRoot, octVec, iLambda=iLam,  foundOctal=thisOctal, &
               foundSubcell=subcell, kappaSca=kappaScaReal, kappaAbs=kappaAbsReal, grid=grid)
          oldOctal => thisOctal
          call intersectCubeAMR(grid, rVec, uHat, tVal)
          octVec = rVec
          thisTau = dble(tVal) * dble(kappaAbsReal + kappaScaReal)
          if (tVal == 0.) then
             escaped = .true.
             stillingrid = .false.
          endif
       endif
    enddo

    if (.not.escaped) then
       octVec = rVec
       if (.not.inOctal(grid%octreeRoot, octVec)) then
          write(*,*) "octVec",octVec
          write(*,*) "size",grid%octreeRoot%subcellsize
          stop
       endif
       if (dble(tau)/thisTau .gt. 1.d0) then
          write(*,*) "tau prob",tau,thisTau
       endif
       call intersectCubeAMR(grid, rVec, uHat, tVal)
       call amrGridValues(grid%octreeRoot, octVec, startOctal=oldOctal,iLambda=iLam, foundOctal=thisOctal, foundSubcell=subcell, &
            kappaAbs=kappaAbsReal,kappaSca=kappaScaReal, grid=grid)
          if (thisTau > 1.e-30) then
             thisOctal%distanceGrid(subcell) = thisOctal%distanceGrid(subcell) + (dble(tVal)*dble(tau)/thisTau) * dble(kappaAbsReal)
             thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1
             oldOctal => thisOctal
          endif

       if (tau > thisTau) then
          write(*,*) "tau > thistau"
          stop
       endif
       if (thisTau > 1.e-30) then
          rVec = rVec + (dble(tVal)*dble(tau)/thisTau) * uHat
       endif
    endif

 end subroutine toNextEventAMR

  recursive subroutine removeDust(thisOctal, Tthresh, nRemoved)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real :: Tthresh
  integer :: nRemoved
  integer :: subcell, i
  
  do subcell = 1, 8
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call removeDust(child, Tthresh, nRemoved)
                exit
             end if
          end do
       else
          if (thisOctal%temperature(subcell) > Tthresh) then
             write(*,*) "cell removed with T at ",thisOctal%temperature(subcell)
             thisOctal%inFlow(subcell) = .false.
             thisOctal%etaCont(subcell) = 1.e-30
             thisOctal%rho(subcell) = 1.e-20
             thisOctal%temperature(subcell) = 1.e-3
             nRemoved = nRemoved + 1
          endif
       endif
    enddo
  end subroutine removeDust

  recursive subroutine setbiasAMR(thisOctal, grid)
  type(gridtype) :: grid
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  real :: r
  
  do subcell = 1, 8
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setBiasAMR(child, grid)
                exit
             end if
          end do
       else
          r = modulus(subcellcentre(thisOctal, subcell)) / grid%rInner
          thisOctal%biasCont3D(subcell) = r
       endif
    enddo
  end subroutine setBiasAMR


end module lucy_mod

