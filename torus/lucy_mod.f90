
module lucy_mod

  use constants_mod
  use grid_mod
  use vector_mod
  use phasematrix_mod
  implicit none


contains


  subroutine lucyRadiativeEquilibrium(grid, miePhase, nMuMie, nLambda, lamArray, temperature)

    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec, uHat, uNew
    integer :: nLambda, nMuMie
    type(PHASEMATRIX):: miePhase(1:nLambda, 1:nMuMie)
    real :: lamArray(:)
    real(kind=doublekind), allocatable :: distanceGrid(:, :, :)
    real :: r
    integer, parameter :: nFreq = 1000
    integer :: i, j
    real :: freq(nFreq), dnu(nFreq), probDistPlanck(nFreq), probDistJnu(nFreq)
    real :: temperature
    integer :: nMonte = 100000, iMonte
    real :: thisFreq,  logFreqStart, logFreqEnd
    real :: albedo
    logical :: escaped
    integer :: i1, i2, i3
    real :: t1, t2, t3
    real :: thisLam
    integer :: iLam
    real :: cosTheta
    integer :: nInf
    integer :: iIter, nIter = 5
    real(kind=doubleKind) ::  kappaP
    real(kind=doubleKind) :: adot, V, epsOverDeltaT
    real :: newT, deltaT
    real :: meanDeltaT
    real :: dx, dy, tr(6), fg, bg
    real, allocatable :: tempImage(:,:)
    integer :: nDT
    real :: totalEmission

    allocate(tempImage(1:grid%nx,1:grid%ny))


    logFreqStart = 12.
    logFreqEnd = 17.
    do i = 1, nFreq
       freq(i) = logFreqStart + (real(i-1)/real(nFreq-1))*(logFreqEnd-logFreqStart)
       freq(i) = 10.**freq(i)
    enddo
    write(*,*) "Lam",(cSpeed/freq(1))*1.e8,(cSpeed/freq(nFreq))*1.e8

    do i = 2, nFreq
       dnu(i) = freq(i)-freq(i-1)
    enddo

    dnu(1) = 0.



    call setupFreqProb(temperature, freq, dnu, nFreq, ProbDistPlanck)


    write(*,'(a)') "Computing lucy radiative equilibrium..."

    allocate(distanceGrid(1:grid%nx, 1:grid%ny, 1:grid%nz))



    do iIter = 1, nIter
       distanceGrid = 0.d0
       write(*,*) "Iteration",iIter
       nInf = 0
       do iMonte = 1, nMonte
          escaped = .false.
          call random_number(r)
          call locate(probDistPlanck, nFreq, r, j)
          thisFreq = freq(j) + (freq(j+1) - freq(j))* &
               (r - probDistPlanck(j))/(probDistPlanck(j+1)-probDistPlanck(j))


          rVec  = grid%rCore * randomUnitVector()
          uHat = randomUnitVector()
         if ( (rVec .dot. uHat) < 0.) uHat = (-1.) * uHat

         do while(.not.escaped)
            call toNextEvent(grid, rVec, uHat, escaped, distanceGrid, thisFreq, nLambda, lamArray)

            if (escaped) nInf = nInf + 1

            if (.not. escaped) then

               thisLam = (cSpeed / thisFreq) * 1.e8
               call locate(lamArray, nLambda, thisLam, iLam)

               call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
               albedo = grid%kappaSca(i1,i2,i3,iLam)/(grid%kappaSca(i1,i2,i3,iLam)+grid%kappaAbs(i1,i2,i3,iLam))
               call random_number(r)
               if (r < albedo) then
                


                  uNew = newDirectionMie(uhat, thisLam, lamArray, nLambda, miePhase, nMuMie)

                  cosTheta = uHat .dot. uNew
                  uHat = uNew

               else
                  
                  probDistJnu(1) = 0.
                  do i = 2, nFreq
                     thisLam = (cSpeed / freq(i)) * 1.e8
                     call hunt(lamArray, nLambda, thisLam, iLam)
                     if (iLam < 1) iLam = 1
                     if (iLam > nLambda) iLam=nLambda
                     probDistJnu(i) = probDistJnu(i-1) + &
                          real(bnu(dble(freq(i)),dble(grid%temperature(i1,i2,i3)))) * grid%kappaAbs(i1,i2,i3,iLam) * dnu(i)
                  enddo
                  probDistJnu(1:nFreq) = probDistJnu(1:nFreq) / probDistJnu(nFreq)
                  call random_number(r)
                  call locate(probDistJnu, nFreq, r, j)
                  if (j == nFreq) j = nFreq -1
                  thisFreq = freq(j) + (freq(j+1) - freq(j))* &
                       (r - probDistJnu(j))/(probDistJnu(j+1)-probDistJnu(j))
               endif
               
            endif
         enddo
       enddo
       write(*,'(a,f7.2)') "Photons done.",real(ninf)/real(nmonte)

       epsOverDeltaT = (dble(grid%lCore)) / dble(nInf)

       V = (dble(grid%xAxis(2)-grid%xAxis(1))*dble(grid%yAxis(2)-grid%yAxis(1))*dble(grid%zAxis(2)-grid%zAxis(1)))

       meanDeltaT = 0.
       nDT = 0


       totalEmission = 0.
       do i1 = 1, grid%nx
          do i2 = 1, grid%ny
             do i3 = 1, grid%nz
                if (grid%inUse(i1,i2,i3)) then
                   adot = epsoverDeltaT * (1.d0 / v) * distancegrid(i1,i2,i3) / 1.d30
                   kappaP = 0.d0
                   do i = 1, nFreq
                      thisLam = (cSpeed / freq(i)) * 1.e8
                      call hunt(lamArray, nLambda, thisLam, iLam)
                      if (iLam < 1) ilam = 1
                      if (ilam > nLambda) ilam = nlambda
                      kappaP = kappaP + dble(grid%kappaAbs(i1,i2,i3,iLam)) * &
                           dble(bnu(dble(freq(i)),dble(grid%temperature(i1,i2,i3)))) * dble(dnu(i))
                   enddo
                   kappaP = kappaP * (pi / (stefanBoltz * grid%temperature(i1,i2,i3)**4))/1.e10

                   
                   newT = (pi / stefanBoltz) * aDot / (fourPi * kappaP)
                   newT = newT**0.25
 

                   deltaT = newT - grid%temperature(i1,i2,i3)
                   grid%temperature(i1,i2,i3) = max(1.e-3,grid%temperature(i1,i2,i3) + 0.8 * deltaT)
                   nDT = nDT  + 1
                   meanDeltaT = meanDeltaT + deltaT
                   kappaP = 0.d0
                   do i = 1, nFreq
                      thisLam = (cSpeed / freq(i)) * 1.e8
                      call hunt(lamArray, nLambda, thisLam, iLam)
                      if (iLam < 1) ilam = 1
                      if (ilam > nLambda) ilam = nlambda
                      kappaP = kappaP + dble(grid%kappaAbs(i1,i2,i3,iLam)) * &
                           dble(bnu(dble(freq(i)),dble(grid%temperature(i1,i2,i3)))) * dble(dnu(i))
                   enddo
                   kappaP = kappaP * (pi / (stefanBoltz * grid%temperature(i1,i2,i3)**4))/1.e10
                   grid%etaCont(i1,i2,i3) = fourPi * kappaP * (stefanBoltz/pi) * (grid%temperature(i1,i2,i3)**4)
                   totalEmission = totalEmission + grid%etaCont(i1,i2,i3) * V 

                endif
             enddo
          enddo
       enddo
       write(*,*) meanDeltaT / real(nDT)

       write(*,*) "Luminosity of dust / core", totalEmission * 1.e30/grid%lCore

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
  end subroutine lucyRadiativeEquilibrium


  subroutine setupFreqProb(temperature, freq, dnu, nFreq, probDist)

    ! Lucy 1999, A&A, 344, 282 Equation 3

    real :: temperature
    integer :: nFreq
    real :: freq(nFreq)
    real :: probDist(nFreq)
    real :: dnu(*)

    integer :: i


    
    probDist = 0.
    do i = 2, nFreq
       probDist(i) = probDist(i-1) + &
            real(bnu(dble(freq(i)),dble(temperature))) * dnu(i)

    enddo

    probDist(1:nFreq) = probDist(1:nFreq) / probDist(nFreq)

  end subroutine setupFreqProb


 subroutine toNextEvent(grid, rVec, uHat,  escaped, distanceGrid, thisFreq, nLambda, lamArray)


   type(GRIDTYPE) :: grid
   type(VECTOR) :: rVec, uHat
   integer :: i1, i2, i3
   real :: t1, t2, t3, tval, tau, r
   real :: lamArray(*)
   integer :: nLambda
   logical :: stillinGrid
   logical :: escaped
   real :: thisTau
   real :: thisFreq
   real(kind=doublekind) :: distanceGrid(1:grid%nx, 1:grid%ny, 1:grid%nz)
   real :: thisLam
   integer :: iLam

    stillinGrid = .true.
    escaped = .false.


    thisLam = (cSpeed / thisFreq) * 1.e8
    call hunt(lamArray, nLambda, thisLam, iLam)
    if ((ilam < 1).or.(ilam > nlambda)) then
       write(*,*) "ilam errro",ilam
    endif

    call random_number(r)
    tau = -log(r)
    call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
    call intersectCube(grid, rVec, i1, i2, i3, uHat, tVal)
    thisTau = tVal * (grid%kappaAbs(i1,i2,i3,iLam)+grid%kappaSca(i1,i2,i3,iLam))

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
          distanceGrid(i1,i2,i3) = distanceGrid(i1,i2,i3) + tVal * grid%kappaAbs(i1,i2,i3,iLam)
       endif


       if (stillinGrid) then
          call random_number(r)
          tau = -log(r)
          call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
          call intersectCube(grid, rVec, i1, i2, i3, uHat, tVal)
          thisTau = tVal * (grid%kappaAbs(i1,i2,i3,iLam)+grid%kappaSca(i1,i2,i3,iLam))
          if (tVal == 0.) then
             escaped = .true.
             stillingrid = .false.
          endif
       endif
    enddo

    if (.not.escaped) then
       call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
       call intersectCube(grid, rVec, i1, i2, i3, uHat, tVal)

       distanceGrid(i1,i2,i3) = distanceGrid(i1,i2,i3) + (tVal*tau/thisTau) * grid%kappaAbs(i1,i2,i3,iLam)
       
       if ((tVal*tau/thisTau) > 2.*(grid%xAxis(2)-grid%xAxis(1))) then
          write(*,*) "tval*tau/thistau too big", tval*tau/thisTau
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
   real :: t(6),tval,denom
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

      denom = norm(i) .dot. direction
      if (denom /= 0.) then
         t(i) = (norm(i) .dot. (p3(i)-posVec))/denom
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
  tval = max(tval * 1.000001,(grid%xAxis(2)-grid%xAxis(1))/100.)


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
     write(*,*) "tval too big",tval
  endif


  end subroutine intersectCube 


end module lucy_mod

