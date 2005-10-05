! photoionization module - started on October 4th 2005 by th

module photoion_mod

use source_mod
use grid_mod
use amr_mod
use constants_mod

implicit none

contains

  subroutine photoIonizationloop(grid, source, nSource, nLambda, lamArray)
    implicit none
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: nlambda
    real :: lamArray(:)
    integer :: nSource
    type(SOURCETYPE) :: source(:), thisSource
    integer :: iSource
    type(OCTALVECTOR) :: rVec, uHat, rHat
    real(double) :: alphaA, alpha1, v, lCore
    integer :: nMonte, iMonte
    integer :: subcell
    integer :: i
    logical :: escaped
    real(double) :: wavelength, thisFreq
    real :: thisLam
    type(OCTALVECTOR) :: outVec, octVec
    real(double) :: r
    integer :: ilam
    integer :: nInf
    real(double) :: kappaScadb, kappaAbsdb, nuHydrogen
    real(double) :: epsOverDeltaT
    real :: dummy(3)
    integer :: nIter

    lCore = 0.d0
    do i = 1, nSource
       lCore = lCore + source(i)%luminosity
    enddo
!$MPI    if(my_rank == 0) &
    write(*,'(a,1pe12.5)') "Total souce luminosity (lsol): ",lCore/lSol

    
    nMonte = 50000
    nIter = 0

    do 
       nIter = nIter + 1
       call zeroDistanceGrid(grid%octreeRoot)
       nInf=0
    call plot_AMR_values(grid, "ionization", "x-z", 0., &
         "ionization.ps/vcps", .true., .false., &
         0, dummy, dummy, dummy, real(grid%octreeRoot%subcellsize), .false.)

       call plot_AMR_values(grid, "temperature", "x-z", real(grid%octreeRoot%centre%y), &
            "lucy_temp_xz.ps/vcps", .true., .false., &
            0, dummy, dummy, dummy, real(grid%octreeRoot%subcellsize), .false.) 


       write(*,*) "Running loop with ",nmonte," photons."
    do iMonte = 1, nMonte

       call randomSource(source, nSource, iSource)
       thisSource = source(iSource)
       call getPhotonPositionDirection(thisSource, rVec, uHat,rHat)
       escaped = .false.
       call getWavelength(thisSource%spectrum, wavelength)
       thisFreq = cSpeed/(wavelength / 1.e8)

          do while(.not.escaped)

             call toNextEventPhoto(grid, rVec, uHat, escaped, thisFreq, nLambda, lamArray)

             if (.not. escaped) then

                
                thisLam = (cSpeed / thisFreq) * 1.e8
                call locate(lamArray, nLambda, real(thisLam), iLam)
                octVec = rVec 
                call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, foundsubcell=subcell,iLambda=iLam, &
                     lambda=real(thisLam), kappaSca=kappaScadb, kappaAbs=kappaAbsdb, grid=grid)

                ! recombination to continuum (equation 26 of Lucy MC trans prob II)

                alphaA = hRecombination(thisOctal%temperature(subcell))
                alpha1 = recombToGround(thisOctal%temperature(subcell))
                call random_number(r)

                nuhydrogen = 13.6 / ergtoev / hCgs
                if (r < (alpha1 / alphaA)) then                
                   call random_number(r)
                   thisFreq = nuHydrogen * (1.d0 - ((kerg*thisOctal%temperature(subcell))/(hCgs*nuHydrogen)*log(r)))
                else
                   escaped = .true.
                endif
             endif
          enddo
          nInf = nInf + 1
       end do
    epsOverDeltaT = (lCore) / dble(nInf)

    call calculateIonizationBalance(grid%octreeRoot, epsOverDeltaT)
!    call calculateThermalBalance(grid%octreeRoot, epsOverDeltaT)
    
 enddo

  end subroutine photoIonizationloop


 subroutine toNextEventPhoto(grid, rVec, uHat,  escaped,  thisFreq, nLambda, lamArray)

   type(GRIDTYPE) :: grid
   type(OCTALVECTOR) :: rVec,uHat, octVec,thisOctVec, tvec
   type(OCTAL), pointer :: thisOctal, tempOctal, sourceOctal
   type(OCTAL),pointer :: oldOctal
   type(OCTAL),pointer :: foundOctal
   integer :: subcell, isubcell, tempSubcell, sourceSubcell
   real(oct) :: tval, tau, r
   real :: lamArray(*)
   integer :: nLambda
   logical :: stillinGrid
   logical :: escaped
   logical :: twoD
   real(double) :: kappaScaDb, kappaAbsDb
   real(oct) :: thisTau
   real(oct) :: thisFreq
   real(oct) :: thisLam
   integer :: iLam
   logical ::inFlow
   integer :: imonte
   real :: diffusionZoneTemp
   logical :: photonInDiffusionZone
   logical :: leftHandBoundary
   real(double) :: prob
   real :: e, h0

    stillinGrid = .true.
    escaped = .false.
       


    thisLam = (cSpeed / thisFreq) * 1.e8
    call hunt(lamArray, nLambda, real(thisLam), iLam)
    if ((ilam < 1).or.(ilam > nlambda)) then
       write(*,*) "ilam errro",ilam
    endif

! select an initial random tau and find distance to next cell

    call random_number(r)
    tau = -log(r)
    if (grid%octreeRoot%threed) then
       call intersectCubeAMR(grid, rVec, uHat, tVal)
    else
       call intersectCubeAMR2D(grid, rVec, uHat, tVal)
    endif


    octVec = rVec
    thisOctVec = rVec

    call locate(lamArray, nLambda, real(thisLam), iLam)


    call amrGridValues(grid%octreeRoot, octVec, iLambda=iLam,  lambda=real(thisLam), foundOctal=thisOctal, &
         foundSubcell=subcell, kappaSca=kappaScadb, kappaAbs=kappaAbsdb, &
         grid=grid, inFlow=inFlow)
    oldOctal => thisOctal

    if (inFlow) then
       thisTau = dble(tVal) * (kappaAbsdb + kappaScadb)
    else
       thisTau = 1.0e-28
    end if

! if tau > thisTau then the photon has traversed a cell with no interactions

    do while(stillinGrid .and. (tau > thisTau)) 

! add on the distance to the next cell

       rVec = rVec + tVal * uHat

       octVec = rVec


! check whether the photon has escaped from the grid

       if (.not.inOctal(grid%octreeRoot, octVec)) then
          stillinGrid = .false.
          escaped = .true.
       endif


! update the distance grid

       e = thisFreq * hCgs * ergtoEV


       call phfit2(1, 1, 1 , e , h0)
       thisOctal%distanceGrid(subcell) = thisOctal%distanceGrid(subcell) &
            + tVal * dble(h0 / (hCgs * thisFreq))

          thisOctal%HIheating(subcell) = thisOctal%HIheating(subcell) &
            + dble(tVal) * dble(h0 / (hCgs * thisFreq)) &
            * (dble(hCgs * thisFreq) - dble(hCgs * nuHydrogen))



       thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1
          

! now if the photon is in the grid choose a new random tau

       if (stillinGrid) then
          call random_number(r)
          tau = -log(r)
          call amrGridValues(grid%octreeRoot, octVec, iLambda=iLam,  lambda=real(thisLam), foundOctal=thisOctal, &
               foundSubcell=subcell, kappaSca=kappaScadb, kappaAbs=kappaAbsdb, &
               grid=grid, inFlow=inFlow)
          oldOctal => thisOctal
          thisOctVec = octVec

! calculate the distance to the next cell

          if (grid%octreeRoot%threed) then
             call intersectCubeAMR(grid, rVec, uHat, tVal)
          else
             call intersectCubeAMR2D(grid, rVec, uHat, tVal)
          endif
          octVec = rVec

! calculate the optical depth to the next cell boundary

          if (inFlow) then
             thisTau = dble(tVal) * (kappaAbsdb + kappaScadb)
          else
             thisTau = 1.0e-28
          end if

          if (tVal == 0.0d0) then
             escaped = .true.
             stillingrid = .false.
          endif
       endif ! still in grid

! if photon is still in grid and  tau > tau_to_the_next_cell then loop round again
! choosing a new tau

       
    enddo
    

! the photon may have escaped the grid...

    if (.not.inOctal(grid%octreeRoot, octVec))  escaped = .true.

 ! if not the photon must interact in this cell
       

    if (.not.escaped) then
       octVec = rVec
!       if (.not.inOctal(grid%octreeRoot, octVec)) then
!          write(*,*) "Error:: Photon location is out of boundaries, but its status is not ESCAPED."
!          write(*,*) "        .... [lucy_mod::toNextEventAMR]"
!          write(*,*) "octVec-centre = ",octVec-grid%octreeRoot%centre
!          write(*,*) "cell size = ",grid%octreeRoot%subcellsize
!          stop
!       endif
       if (dble(tau)/thisTau .gt. 1.d0) then
          write(*,*) "tau prob",tau,thisTau
       endif

       call amrGridValues(grid%octreeRoot, octVec, startOctal=oldOctal,iLambda=iLam, lambda=real(thisLam),&
            foundOctal=thisOctal, foundSubcell=subcell, & 
            kappaAbs=kappaAbsdb,kappaSca=kappaScadb, grid=grid, inFlow=inFlow)
       thisOctVec = octVec


       if (thisOctal%diffusionApprox(subcell)) then
          write(*,*) "Photon in diff zone but not escaped"
       endif

       if (.not.inFlow) kappaAbsdb =0.0d0


! update the distance grid

       if (thisTau > 0.d0) then


          e = thisFreq * hCgs * ergtoEV
          call phfit2(1, 1, 1 , e , h0)
          thisOctal%distanceGrid(subcell) = thisOctal%distanceGrid(subcell) &
            + (dble(tVal)*dble(tau)/thisTau) * dble(h0 / (hCgs * thisFreq))

          thisOctal%HIheating(subcell) = thisOctal%HIheating(subcell) &
            + (dble(tVal)*dble(tau)/thisTau) * dble(h0 / (hCgs * thisFreq)) &
            * (dble(hCgs * thisFreq) - dble(hCgs * nuHydrogen))

          thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1
          

          oldOctal => thisOctal
          
       endif

       if (tau > thisTau) then
          write(*,*) "tau > thistau"
          stop
       endif

! move the requisite distance within the cell and return

       tVec = rVec
       rVec = rVec + (dble(tVal)*dble(tau)/thisTau) * uHat


! this is a workaround due to numerical problems with long pathlengths
! and small cells. needs fixing.

       call amrGridValues(grid%octreeRoot, rVec,  foundOctal=tempOctal, &
            foundSubcell=tempsubcell)

    endif

666 continue

 end subroutine toNextEventPhoto


  subroutine intersectCube(grid, posVec, i1,i2,i3,direction, tval)
   use vector_mod
   use grid_mod
   implicit none
   type(GRIDTYPE) :: grid
   type(OCTALVECTOR) :: direction
   type(OCTALVECTOR) :: posVec, norm(6), p3(6)
   real(oct) :: t(6),tval,denom(6)
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
  subroutine intersectCubeAMR(grid, posVec, direction, tval)
   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(OCTALVECTOR), intent(in) :: posVec
   type(OCTALVECTOR), intent(in) :: direction
   real(oct), intent(out) :: tval
   !
   type(OCTALVECTOR) :: norm(6), p3(6)
   type(OCTAL),pointer :: thisOctal
   type(OCTALVECTOR) :: subcen, point
   integer :: subcell
   
   real(oct) :: t(6),denom(6)
   integer :: i,j
   logical :: ok, thisOk(6)


   point = posVec

   call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
   subcen =  subcellCentre(thisOctal,subcell)
   ok = .true.

   norm(1) = OCTALVECTOR(1.0d0, 0.d0, 0.0d0)
   norm(2) = OCTALVECTOR(0.0d0, 1.0d0, 0.0d0)
   norm(3) = OCTALVECTOR(0.0d0, 0.0d0, 1.0d0)
   norm(4) = OCTALVECTOR(-1.0d0, 0.0d0, 0.0d0)
   norm(5) = OCTALVECTOR(0.0d0, -1.0d0, 0.0d0)
   norm(6) = OCTALVECTOR(0.0d0, 0.0d0, -1.0d0)

   p3(1) = OCTALVECTOR(subcen%x+thisOctal%subcellsize/2.0d0, subcen%y, subcen%z)
   p3(2) = OCTALVECTOR(subcen%x, subcen%y+thisOctal%subcellsize/2.0d0 ,subcen%z)
   p3(3) = OCTALVECTOR(subcen%x,subcen%y,subcen%z+thisOctal%subcellsize/2.0d0)
   p3(4) = OCTALVECTOR(subcen%x-thisOctal%subcellsize/2.0d0, subcen%y,  subcen%z)
   p3(5) = OCTALVECTOR(subcen%x,subcen%y-thisOctal%subcellsize/2.0d0, subcen%z)
   p3(6) = OCTALVECTOR(subcen%x,subcen%y,subcen%z-thisOctal%subcellsize/2.0d0)

   thisOk = .true.

   do i = 1, 6

      denom(i) = norm(i) .dot. direction
      if (denom(i) /= 0.0d0) then
         t(i) = (norm(i) .dot. (p3(i)-posVec))/denom(i)
      else
         thisOk(i) = .false.
         t(i) = 0.0d0
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
     write(*,*) "Error: j=0 (no intersection???) in lucy_mod::intersectCubeAMR. "
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

  subroutine intersectCubeAMR2D(grid, posVec, direction, tval)

! this is to find a cell intersection for a 2D AMR grid
! which is essentially a 2D-grid that projects into
! a 3D cylindrical coordinate system


   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(OCTALVECTOR), intent(in) :: posVec
   type(OCTALVECTOR), intent(inout) :: direction
   real(oct), intent(out) :: tval
   !
   type(OCTALVECTOR) :: norm(6), p3(6)
   type(OCTAL),pointer :: thisOctal
   type(OCTALVECTOR) :: subcen, point
   real :: dx, dz ! directions in cylindrical coords
   integer :: subcell
   real(double) :: compX, compZ,currentX,currentZ
   real(double) :: distToZBoundary, distToXboundary
   real(oct) :: r1,r2,d,cosmu,x1,x2,distTor1,distTor2, theta, mu
   integer :: i,j
   logical :: ok, thisOk(6)
   type(OCTALVECTOR) :: xHat, zHAt

   point = posVec

   call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
   subcen =  subcellCentre(thisOctal,subcell)
   ok = .true.


   r1 = subcen%x - thisOctal%subcellSize/2.d0
   r2 = subcen%x + thisOctal%subcellSize/2.d0
   d = sqrt(point%x**2+point%y**2)
   xHat = VECTOR(point%x, point%y,0.d0)
   call normalize(xHat)

   cosmu =((-1.d0)*xHat).dot.direction
   call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
   if (.not.ok) then
      write(*,*) "Quad solver failed in intersectcubeamr2d"
      direction = randomUnitVector()
      x1 = thisoctal%subcellSize/2.d0
      x2 = 0.d0
   endif
   distTor2 = max(x1,x2)

   theta = asin(max(-1.d0,min(1.d0,r1 / d)))
   cosmu = xHat.dot.direction
   mu = acos(max(-1.d0,min(1.d0,cosmu)))
   distTor1 = 1.e30
   if (mu  < theta ) then
      call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
      if (.not.ok) then
         write(*,*) "Quad solver failed in intersectcubeamr2d"
         direction = randomUnitVector()
         x1 = thisoctal%subcellSize/2.d0
         x2 = 0.d0
      endif
      distTor1 = max(x1,x2)
   endif
         
   distToXboundary = min(distTor1, distTor2)


   zHat = VECTOR(0.d0, 0.d0, 1.d0)
   compZ = zHat.dot.direction
   currentZ = point%z

   if (compZ /= 0.d0 ) then
      if (compZ > 0.d0) then
         distToZboundary = (subcen%z + thisOctal%subcellsize/2.d0 - currentZ ) / compZ
      else
         distToZboundary = abs((subcen%z - thisOctal%subcellsize/2.d0 - currentZ ) / compZ)
      endif
   else
      disttoZboundary = 1.e30
   endif

   tVal = min(distToZboundary, distToXboundary) +0.0001d0*grid%halfsmallestsubcell
   if (tVal > 1.e29) then
      write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
      write(*,*) "subcen",subcen
      write(*,*) "x,z",currentX,currentZ
   endif
   if (tval < 0.) then
      write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
      write(*,*) "subcen",subcen
      write(*,*) "x,z",currentX,currentZ
   endif

  end subroutine intersectCubeAMR2D

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
          thisOctal%incidentFlux(subcell) = 0.
          thisOCtal%nDiffusion(subcell) = 0.
          thisOctal%hIheating(subcell) = 0.d0
       endif
    enddo
  end subroutine zeroDistanceGrid

  recursive subroutine calculateIonizationBalance(thisOctal, epsOverDeltaT)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: epsOverDeltaT, v, r1, r2, ionizationCoeff, recombCoeff
  integer :: subcell, i
  type(OCTALVECTOR) :: rVec
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateIonizationBalance(child, epsOverDeltaT)
                exit
             end if
          end do
       else
          if (thisOctal%inflow(subcell)) then
             if (thisOctal%ncrossings(subcell) > 2) then
                if (thisOctal%threed) then
                   v = thisOctal%subcellSize**3
                else
                   rVec = subcellCentre(thisOctal,subcell)
                   r1 = rVec%x-thisOctal%subcellSize/2.d0
                   r2 = rVec%x+thisOctal%subcellSize/2.d0
                   v = dble(pi) * (r2**2 - r1**2) * thisOctal%subcellSize
                endif
                ionizationCoeff = (epsOverDeltaT / (v * 1.d30))*thisOctal%distanceGrid(subcell) ! equation 44 Lucy MC II
                recombCoeff = thisOctal%ne(subcell) *  hRecombination(thisOctal%temperature(subcell))
                if ((recombCoeff == 0.d0).or.((recombCoeff+ionizationCoeff) ==0.d0)) then
                   write(*,*) "!",recombCoeff,ionizationCoeff,thisOctal%temperature(subcell)
                else
                   thisOctal%nHI(subcell) = thisOctal%nh(subcell) * recombCoeff/(ionizationCoeff  + recombCoeff)
                   thisOctal%nHII(subcell) = thisOctal%nh(subcell)  - thisOctal%nhi(subcell)
                   thisOctal%ne(subcell) = thisOctal%nhii(subcell)
                endif
             else
                write(*,*) "Undersampled cell",thisOctal%ncrossings(subcell)
             endif
          endif
       endif
    enddo
  end subroutine calculateIonizationBalance

  recursive subroutine calculateThermalBalance(thisOctal, epsOverDeltaT)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: epsOverDeltaT, v, r1, r2, ionizationCoeff, recombCoeff
  integer :: subcell, i
  type(OCTALVECTOR) :: rVec
  logical :: converged
  real :: t1, t2, tm
  real(double) :: y1, y2, ym, HIheating
  real :: deltaT
  real :: underCorrection = 0.5
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateThermalBalance(child, epsOverDeltaT)
                exit
             end if
          end do
       else
          if (thisOctal%inflow(subcell)) then
             if (thisOctal%ncrossings(subcell) > 2) then
                if (thisOctal%threed) then
                   v = thisOctal%subcellSize**3
                else
                   rVec = subcellCentre(thisOctal,subcell)
                   r1 = rVec%x-thisOctal%subcellSize/2.d0
                   r2 = rVec%x+thisOctal%subcellSize/2.d0
                   v = dble(pi) * (r2**2 - r1**2) * thisOctal%subcellSize
                endif
                HIheating= thisOctal%nHi(subcell) * (epsOverDeltaT / (v * 1.d30))*thisOctal%hIheating(subcell) ! equation 21 of kenny's

                if (HIheating == 0.d0) then
                   thisOctal%temperature(subcell) = 3.
                else
                   converged = .false.
                   t1 = 1.e-5
                   t2 = 5.e6
                   y1 = (HHecooling(thisOctal%nHii(subcell), thisOctal%ne(subcell), t1) - HIheating)
                   y2 = (HHecooling(thisOctal%nHii(subcell), thisOctal%ne(subcell), t2) - HIheating)
                   if (y1*y2 > 0.d0) then
                      write(*,*) "Insufficient range"
                      write(*,*) "t1",t1,y1,HHecooling(thisOctal%nHii(subcell), thisOctal%ne(subcell), t1),HIheating
                      write(*,*) "t2",t2,y2,HHecooling(thisOctal%nHii(subcell), thisOctal%ne(subcell), t2),HIheating
                   endif
                   
                   ! find root of heating and cooling by bisection
                   
                   do while(.not.converged)
                      tm = 0.5*(t1+t2)
                      y1 = (HHecooling(thisOctal%nHii(subcell), thisOctal%ne(subcell), t1) - HIheating) * 1.d20
                      y2 = (HHecooling(thisOctal%nHii(subcell), thisOctal%ne(subcell), t2) - HIheating) * 1.d20
                      ym = (HHecooling(thisOctal%nHii(subcell), thisOctal%ne(subcell), tm) - HIheating) * 1.d20
                      !                write(*,*) t1,t2,y1,y2
                      !                write(*,*) hiheating,HHecooling(thisOctal%nHii(subcell), thisOctal%ne(subcell), tm)
                      if (y1*ym < 0.d0) then
                         t1 = t1
                         t2 = tm
                      else if (y2*ym < 0.d0) then
                         t1 = tm
                         t2 = t2
                      else
                         converged = .true.
                         tm = 0.5*(t1+t2)
                      endif
                      
                      if (abs((t1-t2)/t1) .le. 1.e-3) then
                         converged = .true.
                      endif
                      
                   enddo
                   deltaT = tm - thisOctal%temperature(subcell)
                   thisOctal%temperature(subcell) = thisOctal%temperature(subcell) + underCorrection * deltaT
                endif
             else
                write(*,*) "Undersampled cell",thisOctal%ncrossings(subcell)
             endif
          endif
       endif
    enddo
  end subroutine calculateThermalBalance


  function hRecombination(temperature) result (rate)
    real :: temperature
    real(double) :: rate, t0, t1, b

    t0 = sqrt(temperature/3.148d0)
    t1 = sqrt(temperature/7.036d5)
    b = 0.7480

    rate = 7.982d-11 / (t0*(1.d0+t0)**(1.d0-b) * (1.d0 + t1)**(1.d0+b))
  end function hRecombination

  function recombToGround(temperature) result (alpha1)
    real :: temperature
    real(double) :: alpha1

    alpha1 = 1.58d-13 * (temperature/1.d4)**(-0.53d0)  ! kenny's photo paper equation 24
  end function recombToGround

  function HHeCooling(nIon, ne, temperature) result (coolingRate)

    real(double) :: nIon, ne
    real :: temperature
    real(double) :: coolingRate
    real(double) :: gff
    real :: rootTbeta1(31) = (/ 8.287e-11, 7.821e-11, 7.356e-11, 6.982e-11, 6.430e-11, 5.971e-11, 5.515e-11, 5.062e-11, 4.614e-11, &
                               4.170e-11, 3.734e-11, 3.306e-11, 2.888e-11, 2.484e-11, 2.098e-11, 1.736e-11, 1.402e-11, 1.103e-11, &
                               8.442e-12, 6.279e-12, 4.539e-12, 3.192e-11, 2.185e-12, 1.458e-12, 9.484e-13, 6.023e-13, 3.738e-13, &
                               2.268e-13, 1.348e-13, 7.859e-14, 4.499e-14 /)

    real :: rootTbeta2(31) = (/ 1.646e-11, 1.646e-11, 1.646e-11, 1.646e-11, 1.646e-11, 1.645e-11, 1.644e-11, 1.643e-11, 1.641e-11, &
                                1.638e-11, 1.633e-11, 1.625e-11, 1.613e-11, 1.594e-11, 1.656e-11, 1.522e-11, 1.460e-11, 1.374e-11, &
                                1.260e-11, 1.119e-11, 9.571e-12, 7.844e-12, 6.146e-12, 4.601e-12, 3.295e-12, 2.262e-12, 1.494e-12, &
                                9.520e-13, 5.878e-13, 3.528e-13, 2.066e-13 /)

    real :: logT(31) = (/ 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, &
                      5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0 /)
    real(double) :: fac, thisRootTbeta1, beta, thisRootTbeta2
    real :: thisLogT
    integer :: n

    gff = 1.1d0 + 0.34d0*exp(-(5.5d0 - log10(temperature))**2 / 3)  ! Kenny equation 23

    coolingRate = 1.42d-27 * nIon * ne * gff * sqrt(temperature)  ! Kenny equation 22 (free-free)

    thisLogT = log10(temperature)
    call locate(logT, 31, thisLogT, n)
    fac = (thisLogT - logT(n))/(logT(n+1)-logT(n))
    thisRootTbeta1 = rootTbeta1(n) + (rootTbeta1(n+1)-rootTbeta1(n)) * fac

    thisRootTbeta2 = rootTbeta2(n) + (rootTbeta2(n+1)-rootTbeta2(n)) * fac

    beta = (thisrootTbeta1 + thisRootTbeta2) / sqrt(temperature)

    coolingRate = coolingRate + ne * nion * kerg * temperature * beta

  end function HHeCooling

end module

