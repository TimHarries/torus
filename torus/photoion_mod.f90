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
    real(double) :: kappaScadb, kappaAbsdb
    real(double) :: epsOverDeltaT
    real :: dummy(3)
    real :: pops(10)
    integer :: nIter
    real(double) :: crate

    lCore = 0.d0
    do i = 1, nSource
       lCore = lCore + source(i)%luminosity
    enddo
!$MPI    if(my_rank == 0) &
    write(*,'(a,1pe12.5)') "Total souce luminosity (lsol): ",lCore/lSol

    
    nMonte = 100000
    nIter = 0

!    pops = 0.
!    call solvePops(grid%ion(returnIonNumber("N I", grid%ion, grid%nIon)),pops, 100.d0, 10000.)
!    do i = 1, 5
!       write(*,*) i, pops(i)
!    enddo
!    write(*,*) "SUM: ",sum(pops(1:5))
!    call coolingRate(grid%ion(returnIonNumber("N I", grid%ion, grid%nIon)), 100.d0, 100.d0, 10000., crate)
!    write(*,*) "cooling rate:",crate
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
       call dumpLexington(grid)

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

    call calculateIonizationBalance(grid,grid%octreeRoot, epsOverDeltaT)
    call calculateThermalBalance(grid, grid%octreeRoot, epsOverDeltaT)
    
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

       call updateGrid(grid, thisOctal, subcell, thisFreq, tVal)
          

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

          call updateGrid(grid, thisOctal, subcell, thisFreq, dble(tval)*dble(tau)/thisTau)

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
          thisOctal%Hheating(subcell) = 0.d0
          thisOctal%Heheating(subcell) = 0.d0
          thisOctal%photoIonCoeff(subcell,:) = 0.d0
       endif
    enddo
  end subroutine zeroDistanceGrid

  recursive subroutine calculateIonizationBalance(grid, thisOctal, epsOverDeltaT)
    type(gridtype) :: grid
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
                call calculateIonizationBalance(grid, child, epsOverDeltaT)
                exit
             end if
          end do
       else
          if (thisOctal%inflow(subcell)) then
             if (thisOctal%ncrossings(subcell) > 10) then
                call solveIonizationBalance(grid, thisOctal, subcell, epsOverdeltaT)
             else
!                write(*,*) "Undersampled cell",thisOctal%ncrossings(subcell)
             endif
          endif
       endif
    enddo
  end subroutine calculateIonizationBalance

  recursive subroutine calculateThermalBalance(grid, thisOctal, epsOverDeltaT)
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real(double) :: epsOverDeltaT, v, r1, r2, ionizationCoeff, recombCoeff
    real(double) :: totalHeating, nHii, nHeII, nh
    integer :: subcell, i
    type(OCTALVECTOR) :: rVec
    logical :: converged
    real :: t1, t2, tm
    real(double) :: y1, y2, ym, Hheating, Heheating
    real :: deltaT
    real(double) :: junk
    real :: underCorrection = 0.8
    real :: pops(10)
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateThermalBalance(grid, child, epsOverDeltaT)
                exit
             end if
          end do
       else
          if (thisOctal%inflow(subcell)) then
             if (thisOctal%ncrossings(subcell) > 10) then
                if (thisOctal%threed) then
                   v = thisOctal%subcellSize**3
                else
                   rVec = subcellCentre(thisOctal,subcell)
                   r1 = rVec%x-thisOctal%subcellSize/2.d0
                   r2 = rVec%x+thisOctal%subcellSize/2.d0
                   v = dble(pi) * (r2**2 - r1**2) * thisOctal%subcellSize
                endif
                
                Hheating= thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,1) * grid%ion(1)%abundance &
                     * (epsOverDeltaT / (v * 1.d30))*thisOctal%Hheating(subcell) ! equation 21 of kenny's
                Heheating= thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,3) * grid%ion(3)%abundance &
                     * (epsOverDeltaT / (v * 1.d30))*thisOctal%Heheating(subcell) ! equation 21 of kenny's
                totalHeating = Hheating + HeHeating

                nHii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,2) * grid%ion(2)%abundance
                nHeii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,4) * grid%ion(4)%abundance
                nh = thisOctal%nh(subcell)

                if (totalHeating == 0.d0) then
                   thisOctal%temperature(subcell) = 3.
                else
                   converged = .false.
                   t1 = 1.e1
                   t2 = 20000.
!                   call solvePops(grid%ion(returnIonNumber("O III", grid%ion, grid%nion)), pops, thisOctal%ne(subcell), t1)
!                   do i = 1, 10
!                      write(*,*) i,pops(i)
!                   enddo
!                   stop

                   y1 = (HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t1) - totalHeating)
                   y2 = (HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t2) - totalHeating)
                   if (y1*y2 > 0.d0) then
                      write(*,*) "Insufficient range"
                      write(*,*) "t1",t1,y1,HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t1),totalHeating
                      write(*,*) "t2",t2,y2,HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t2),totalHeating
                   endif
                   
                   ! find root of heating and cooling by bisection
                   
                   do while(.not.converged)
                      tm = 0.5*(t1+t2)
                      y1 = (HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t1) - totalheating)
                      y2 = (HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t2) - totalheating)
                      ym = (HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), tm) - totalheating)

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
!                   call metalcoolingRate(grid%ion, grid%nIon, thisOctal, subcell, nh, thisOctal%ne(subcell), thisOctal%temperature(subcell), junk, .true.)
                endif
             else
!                write(*,*) "Undersampled cell",thisOctal%ncrossings(subcell)
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

  function HHeCooling(grid, thisOctal, subcell, nh, nHii, nHeii, ne, temperature) result (coolingRate)
    type(OCTAL) :: thisOctal
    integer :: subcell
    type(GRIDTYPE) :: grid
    real(double) :: nHii, nHeii, ne, nh
    real :: temperature
    real(double) :: coolingRate, crate
    real(double) :: gff
    real :: rootTbetaH(31) = (/ 8.287e-11, 7.821e-11, 7.356e-11, 6.982e-11, 6.430e-11, 5.971e-11, 5.515e-11, 5.062e-11, 4.614e-11, &
                               4.170e-11, 3.734e-11, 3.306e-11, 2.888e-11, 2.484e-11, 2.098e-11, 1.736e-11, 1.402e-11, 1.103e-11, &
                               8.442e-12, 6.279e-12, 4.539e-12, 3.192e-11, 2.185e-12, 1.458e-12, 9.484e-13, 6.023e-13, 3.738e-13, &
                               2.268e-13, 1.348e-13, 7.859e-14, 4.499e-14 /)

    real :: rootTbetaHe(18) = (/ 8.347e-11, 7.889e-11, 7.430e-11, 6.971e-11, 6.512e-11, 6.056e-11, 5.603e-11, 5.154e-11, 4.710e-11,&
                               4.274e-11, 3.847e-11, 3.431e-11, 3.031e-11, 2.650e-11, 2.291e-11, 1.960e-11, 1.660e-11, 1.394e-11 /)

    real :: logT(31) = (/ 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, &
                      5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0 /)
    real(double) :: fac, thisRootTbetaH, betaH, betaHe, thisRootTbetaHe
    real :: thisLogT
    integer :: n

    gff = 1.1d0 + 0.34d0*exp(-(5.5d0 - log10(temperature))**2 / 3)  ! Kenny equation 23

    coolingRate = 1.42d-27 * (nHii+nHeii) * ne * gff * sqrt(temperature)  ! Kenny equation 22 (free-free)

    thisLogT = log10(temperature)
    call locate(logT, 31, thisLogT, n)
    fac = (thisLogT - logT(n))/(logT(n+1)-logT(n))

    thisRootTbetaH = rootTbetaH(n) + (rootTbetaH(n+1)-rootTbetaH(n)) * fac


    if (thisLogT < logT(18)) then
       thisRootTbetaHe = rootTbetaHe(n) + (rootTbetaHe(n+1)-rootTbetaHe(n)) * fac
    else
       thisRootTbetaHe = 0.d0
    endif

    betaH = thisrootTbetaH / sqrt(temperature)
    betaHe = thisrootTbetaHe / sqrt(temperature)

    coolingRate = coolingRate + ne * nhii * kerg * temperature * betaH


    coolingRate = coolingRate + ne * nheii * kerg * temperature * betaHe

    call metalcoolingRate(grid%ion, grid%nIon, thisOctal, subcell, nh, ne, temperature, crate)
    coolingRate = coolingRate + crate
  end function HHeCooling

  subroutine updateGrid(grid, thisOctal, subcell, thisFreq, distance)
    type(GRIDTYPE) :: grid
    type(OCTAL) :: thisOctal
    integer :: subcell
    real(double) :: thisFreq, distance
    integer :: i 
    real :: e, xsec

    thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1

    e = thisFreq * hCgs * ergtoEV

    do i = 1, grid%nIon
       call phfit2(grid%ion(i)%z, grid%ion(i)%n, grid%ion(i)%outerShell , e , xsec)

       thisOctal%photoIonCoeff(subcell,i) = thisOctal%photoIonCoeff(subcell,i) &
            + distance * dble(xsec / (hCgs * thisFreq))

       ! neutral H heating

       if ((grid%Ion(i)%z == 1).and.(grid%Ion(i)%n == 1)) then
          thisOctal%Hheating(subcell) = thisOctal%Hheating(subcell) &
            + dble(distance) * dble(xsec / (thisFreq * hCgs)) &
            * (dble(hCgs * thisFreq) - dble(hCgs * grid%ion(i)%nuThresh))
       endif

       ! neutral He heating

       if ((grid%Ion(i)%z == 2).and.(grid%Ion(i)%n == 1)) then
          thisOctal%Heheating(subcell) = thisOctal%Heheating(subcell) &
            + dble(distance) * dble(xsec / (thisFreq * hCgs)) &
            * (dble(hCgs * thisFreq) - dble(hCgs * grid%ion(i)%nuThresh))
       endif

    enddo

  end subroutine updateGrid

subroutine solveIonizationBalance(grid, thisOctal, subcell, epsOverdeltaT)
  implicit none
  type(GRIDTYPE) :: grid
  type(OCTAL) :: thisOctal
  real(double) :: epsOverDeltaT, V
  integer :: subcell
  integer :: i, j, k
  real(double), allocatable :: matrixA(:,:), MatrixB(:,:)
  integer :: nIonizationStages
  type(OCTALVECTOR) :: rVec
  real(double) :: r1, r2
  integer :: iStart, iEnd, iIon
  real(double) :: chargeExchangeIonization, chargeExchangeRecombination
  logical :: ok

  if (thisOctal%threed) then
     v = thisOctal%subcellSize**3
  else
     rVec = subcellCentre(thisOctal,subcell)
     r1 = rVec%x-thisOctal%subcellSize/2.d0
     r2 = rVec%x+thisOctal%subcellSize/2.d0
     v = dble(pi) * (r2**2 - r1**2) * thisOctal%subcellSize
  endif

  k = 1
  do while(k <= grid%nIon)

     iStart = k
     iEnd = k+1
     do while(grid%ion(istart)%z == grid%ion(iEnd)%z)
        iEnd = iEnd + 1
     enddo
     iEnd = iEnd - 1
     nIonizationStages = iEnd - iStart + 1

!     write(*,*) "istart",istart,"iend",iend,"nionization",nIonizationStages

     allocate(matrixA(1:nIonizationStages,1:nIonizationStages))
     allocate(matrixB(1:nIonizationStages, 1))
     
     matrixA = 0.d0

     do i = 1, nIonizationStages - 1
        iIon = iStart+i-1
!        write(*,*) i,iion,grid%ion(iion)%species

        call getChargeExchangeRates(grid%ion(iion), thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,1),  &
             thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,2),  &
        chargeExchangeIonization, chargeExchangeRecombination)

        matrixA(i, i) = grid%ion(iIon)%abundance * thisOctal%nh(subcell)  * &
             ((epsOverDeltaT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell, iIon) + chargeExchangeIonization)! equation 44 Lucy MC II
        matrixA(i, i) = max(1.d-30, matrixA(i, i))

        matrixA(i, i+1) = -grid%ion(iIon+1)%abundance * thisOctal%nh(subcell) * &
             (recombRate(grid%ion(iIon),thisOctal%temperature(subcell)) * thisOctal%ne(subcell) + chargeExchangeRecombination)
        matrixA(i, i+1) = min(-1.d-30, matrixA(i, i+1))

        

     enddo
     do i = 1, nIonizationStages
        matrixA(nIonizationStages,i) = 1.d0
     enddo
     matrixB = 0.d0
     matrixB (nIonizationStages,1) = 1.d0

     call gaussj(matrixA, nIonizationStages, nIonizationStages, matrixB, 1, 1, ok)
     if (.not.ok) then
        write(*,*) "Ionization balance failed for element: ",trim(grid%ion(istart)%species(1:2))
        stop
     endif

     do i = 1, nIonizationStages
        iIon = istart+i-1
        thisOctal%ionFrac(subcell, iIon) = matrixB(i,1)
!        write(*,*) grid%ion(iIon)%species, thisOctal%ionFrac(subcell, iIon)
     enddo
!     write(*,*) " " 
     deallocate(matrixA, matrixB)
     k = iEnd + 1
  end do
  
  thisOctal%ne(subcell) = returnNe(thisOctal, subcell, grid%ion, grid%nion)

end subroutine solveIonizationBalance


function recombRate(thisIon, temperature) result (rate)
  type(IONTYPE) :: thisIon
  real :: temperature
  real(double) :: rate
  real :: a, b, t0, t1, c, d, f, t

! recombinations INTO this species

  select case(thisIon%z)

     case(1)
        select case(thisIon%n)
           case(1) ! H I
              a = 7.982e-11
              b = 0.7480
              t0 = 3.148e0
              t1 = 7.036e5
              rate = vernerFerland(a, temperature, t0, t1, b)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate = 0.
        end select
     case(2)
        select case(thisIon%n)
           case(2) ! He I 
              a = 3.294e-11
              b = 0.6910
              t0 = 1.544e1
              t1 = 3.676e7
              rate = vernerFerland(a, temperature, t0, t1, b)
           case(1) ! He II
              a = 1.891e-10
              b = 0.7524
              t0 = 9.370e00
              t1 = 2.774e6
              rate = vernerFerland(a, temperature, t0, t1, b)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select
     case(6)
        select case(thisIon%n)
           case(6) ! C I 
              a = 0.0108
              b = -1.075
              c = 0.2810
              d = -0.0193
              f = -0.1127
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(5) ! C II
              a = 1.8267
              b = 4.1012
              c = 4.8443
              d = 0.2261
              f = 0.5960
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(4) ! C III
              a = 2.3196
              b = 10.7328
              c = 6.8830
              d = -0.1824
              f = 0.4101
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(3) ! C IV
              a = 8.540e-11
              b = 0.5247
              t0 = 5.014e2
              t1 = 1.479e7
              rate = vernerFerland(a, temperature, t0, t1, b)
           case(2) ! C V
              a = 2.765e-10
              b = 0.6858
              t0 = 1.535e2
              t1 = 2.556e7
              rate = vernerFerland(a, temperature, t0, t1, b)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select
     case(7)
        select case(thisIon%n)
           case(7) ! N I 
              a = 0.0
              b = 0.6310
              c = 0.1990
              d = -0.0197
              f = 0.4398
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(6) ! N II
              a = 0.0320
              b = -0.6624
              c = 4.3191
              d = 0.0003
              f = 0.5946
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(5) ! N III
              a = -0.8806
              b = 11.2406
              c = 30.7066
              d = -1.1721
              f = 0.6127
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(4) ! N IV
              a = 0.4134
              b = -4.6319
              c = 25.9172
              d = -2.2290
              f = 0.2360
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(3) ! N V
              a = 1.169e-10
              b = 0.5470
              t0 = 6.793e2
              t1 = 1.650e7
              rate = vernerFerland(a, temperature, t0, t1, b)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select
     case(8)
        select case(thisIon%n)
           case(8) ! O I 
              a = 0.0
              b = 0.0238
              c = 0.0659
              d = 0.0349
              f = 0.5334
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(7) ! O II
              a = -0.0036
              b = 0.7519
              c = 1.5252
              d =-0.0838
              f = 0.2769
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(6) ! O III
              a = 0.0
              b = 21.8790
              c = 16.2730
              d = -0.8020
              f = 1.1899
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(5) ! O IV
              a = 0.00061
              b = 0.2269
              c = 32.1419
              d = 1.9939
              f = -0.0646
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(4) ! O V
              a = -2.8425
              b = 0.2283
              c = 40.4072
              d = -3.4956
              f = 1.7558
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select

     case(10)
        select case(thisIon%n) ! from nussbaumer and storey 1987
           case(10) ! Ne I 
              a = 0.0129
              b =-0.1779
              c = 0.9353
              d =-0.0682
              f = 0.4516
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(9) ! Ne II
              a = 3.6781
              b = 14.1481
              c = 17.1175
              d = -0.5017
              f = 0.2313
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(8) ! Ne III
              a =-0.0254
              b = 5.53650
              c = 17.0727
              d = -0.7725
              f = 0.1702
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(7) ! N IV
              a = 0.7469
              b =-3.2024
              c = 12.1163
              d =-1.0379
              f = 1.8482
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(6) ! N V
              a = 19.9280
              b = 235.0536
              c = 12.5096
              d = 9.1413
              f = 0.1282
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select

     case(16)
        select case(thisIon%n) ! from nussbaumer and storey 1987
           case(16) ! S I 
              rate = 3.e-13
           case(15) ! S II
              rate = 3.e-12
           case(14) ! S III
              rate = 1.5e-11
           case(13) ! S IV
              rate = 2.5e-11
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select



  end select
end function recombRate

real function vernerFerland(a, t, t0, t1, b)
  real :: a, t, t0, t1, b

! based on Verner and Ferland (1996) ApJS 103 467

  vernerFerland = a / (sqrt(t/t0) * (1.+sqrt(t/t0))**(1.-b) * (1.+sqrt(t/t1))**(1.+b))
end function vernerFerland



function nussbaumerStorey1983(t, a, b, c, d, f) result(rate)

! based on Nussbaumer & Storey 1983 AA 126 75
  real :: t, a, b, c, d, f, rate
  rate = 1.e-12 * (a/t + b + c*t + d*t**2)*(t**(-1.5))*exp(-f/t)

end function nussbaumerStorey1983


function returnNe(thisOctal, subcell, ionArray, nion) result (ne)
  real(double) :: ne, tot
  integer :: subcell
  type(OCTAL) :: thisOctal
  type(IONTYPE) :: ionArray(:)
  integer :: nion, i

  tot = 0.d0 
  do i = 1, nIon
     tot = tot + ionArray(i)%abundance * thisOctal%nh(subcell) * &
          thisOctal%ionFrac(subcell, i) * dble(ionArray(i)%z-ionArray(i)%n)
  enddo
  ne = tot
end function returnNe

subroutine dumpLexington(grid)
  type(GRIDTYPE) :: grid
  type(OCTAL), pointer :: thisOctal
  integer :: subcell
  integer :: i, j
  real(double) :: r, theta
  real :: t,hi,hei,oii,oiii,cii,ciii,civ,nii,niii,niv,neii,neiii
  type(OCTALVECTOR) :: octVec

  open(20,file="lexington.dat",form="formatted",status="unknown")

  do i = 1, 50
     r = (1.+5.d0*dble(i-1)/49.d0)*pctocm/1.e10

     t=0;hi=0; hei=0;oii=0;oiii=0;cii=0;ciii=0;civ=0;nii=0;niii=0;niv=0;neii=0;neiii=0
     do j = 1, 100
        call random_number(theta)
        theta = theta * Pi
        
        octVec = OCTALVECTOR(r*sin(theta),0.d0,r*cos(theta))
        
        call amrgridvalues(grid%octreeRoot, octVec,  foundOctal=thisOctal, foundsubcell=subcell)
        
        HI = HI + thisOctal%ionfrac(subcell,returnIonNumber("H I", grid%ion, grid%nIon))
        HeI = HeI + thisOctal%ionfrac(subcell,returnIonNumber("He I", grid%ion, grid%nIon))
        OII = OII + thisOctal%ionfrac(subcell,returnIonNumber("O II", grid%ion, grid%nIon))
        OIII = OIII + thisOctal%ionfrac(subcell,returnIonNumber("O III", grid%ion, grid%nIon))
        CII = CII + thisOctal%ionfrac(subcell,returnIonNumber("C II", grid%ion, grid%nIon))
        CIII = CIII + thisOctal%ionfrac(subcell,returnIonNumber("C III", grid%ion, grid%nIon))
        CIV = CIV + thisOctal%ionfrac(subcell,returnIonNumber("C IV", grid%ion, grid%nIon))
        NII = NII + thisOctal%ionfrac(subcell,returnIonNumber("N II", grid%ion, grid%nIon))
        NIII = NIII + thisOctal%ionfrac(subcell,returnIonNumber("N III", grid%ion, grid%nIon))
        NIV = NIV + thisOctal%ionfrac(subcell,returnIonNumber("N IV", grid%ion, grid%nIon))
        NeII = NeII + thisOctal%ionfrac(subcell,returnIonNumber("Ne II", grid%ion, grid%nIon))
        NeIII = NeIII + thisOctal%ionfrac(subcell,returnIonNumber("Ne III", grid%ion, grid%nIon))
        t  = t + thisOctal%temperature(subcell)
     enddo
     hi = hi / 100.; hei = hei/100.; oii = oii/100; oiii = oiii/100.; cii=cii/100.
     ciii = ciii/100; civ=civ/100.; nii =nii/100.; niii=niii/100.; niv=niv/100.
     neii=neii/100.; neiii=neiii/100.; t=t/100.

     write(20,'(f5.3,f12.1, 1p, 12e9.2)') r*1.e10/pctocm,t,hi,hei,oii,oiii,cii,ciii,civ,nii,niii,niv,neii,neiii
  enddo
  close(20)
end subroutine dumpLexington

subroutine getChargeExchangeRates(thisIon, nHI, nHII, IonRate, recombRate)
  type(IONTYPE) :: thisIon
  real(double) :: nHI, nHII, ionRate, recombRate


!  select case(thisIon%z)
!     case(7)
!        select case(thisIon%n)
!
! case(5)  ! N III

  ionRate = 0.d0; recombRate =0.d0
end subroutine getChargeExchangeRates


subroutine getCollisionalRates(thisIon, iTransition, temperature, excitation, deexcitation)
  real :: excitation, deexcitation
  type(IONTYPE) :: thisIon
  integer :: iTransition, i
  real :: temperature
  real :: thisGamma
  real :: t , fac

  t = max(min(20000.,temperature),5000.)
  call locate(thisIon%transition(iTransition)%t, thisIon%transition(iTransition)%ngamma, t, i)
  fac = (t - thisIon%transition(iTransition)%t(i))/(thisIon%transition(iTransition)%t(i+1) - thisIon%transition(iTransition)%t(i))
  thisGamma = thisIon%transition(iTransition)%gamma(i) + fac * (thisIon%transition(iTransition)%gamma(i+1) - thisIon%transition(iTransition)%gamma(i))

  fac = (8.63e-6 / (thisIon%level(thisIon%transition(iTransition)%i)%g  * sqrt(temperature)) ) * thisGamma
  fac = fac * exp(-thisIon%transition(iTransition)%energy / (kev*temperature))


  excitation = fac
  deexcitation = fac * thisIon%level(thisIon%transition(iTransition)%j)%g / thisIon%level(thisIon%transition(iTransition)%i)%g
end subroutine getCollisionalRates


subroutine metalcoolingRate(ionArray, nIons, thisOctal, subcell, nh, ne, temperature, total, debug)
  type(IONTYPE) :: ionArray(*)
  integer :: nIons, subcell
  type(OCTAL) :: thisOctal
  real(double) :: ne, nh
  real :: temperature
  real(double) :: rate, total
  real :: pops(10)
  integer :: i, j
  logical, optional :: debug

  total = 0.d0
  do j = 1, nIons
     if (ionArray(j)%nTransitions > 0) then
        call solvePops(ionArray(j), pops, ne, temperature)
        rate = 0.d0
        do i = 1, ionArray(j)%nTransitions
           rate = rate + pops(ionArray(j)%transition(i)%j)*ionArray(j)%transition(i)%energy*ionArray(j)%transition(i)%a/ergtoev
        enddo
        rate = rate * ionArray(j)%abundance * nh * thisOctal%ionFrac(subcell, j)
        if (present(debug)) then
           if (debug) then
              write(*,'(a,a,a,1p,e12.4,a,0p, f10.1,a,1pe12.4)') "Contribution from ",trim(ionArray(j)%species),":",rate," at T = ",temperature, &
                   " ion frac ",thisOctal%ionFrac(subcell,j)
           endif
        endif
     endif
     total = total + rate
  enddo
end subroutine metalcoolingRate
  

subroutine solvePops(thisIon, pops, ne, temperature)
  type(IONTYPE) :: thisIon
  real(double) :: ne
  real :: pops(*)
  real :: temperature
  real(double), allocatable :: matrixA(:,:), MatrixB(:,:), tempMatrix(:,:)
  integer :: n, iTrans, i, j
  real :: excitation, deexcitation, rateij, rateji
  logical :: ok

  n = thisIon%nLevels
  allocate(matrixA(1:n+1, 1:n+1), matrixB(1:n+1, 1), tempMatrix(1:n+1,1:n+1))

  matrixA = 1.d-30
  matrixB = 0.d0
  do iTrans = 1, thisIon%nTransitions
     i = thisIon%transition(iTrans)%i
     j = thisIon%transition(iTrans)%j
     call getCollisionalRates(thisIon, iTrans, temperature, excitation, deexcitation)
     rateij = max(1.e-30,excitation * ne)
     rateji = max(1.e-30,deexcitation * ne + thisIon%transition(iTrans)%a)


     matrixA(i,i) = matrixA(i,i)-rateij

     matrixA(j,j) = matrixA(j,j)-rateji

     matrixA(j,i) = matrixA(j,i) + rateij
     matrixA(i,j) = matrixA(i,j) + rateji


  enddo

  do i = 1, n
     matrixA(n+1,i) = 1.d0
  enddo
  matrixB = 0.d0
  matrixB (n+1,1) = 1.d0
  
!  do i = 1, n+1
!     write(*,'(1p,7e12.3)') matrixA(i,1:n)
!  enddo

  tempMatrix = matrixA
  call gaussj(matrixA, n+1, n+1, matrixB, 1, 1, ok)
  if (.not.ok) then
     write(*,*) "Population solver failed for: ",thisIon%species
     write(*,*) "nlevels",thisIon%nLevels,"ntrans",thisIon%nTransitions 
     do i = 1, n+1
        write(*,'(1p,9e12.3)') tempmatrix(i,1:n)
     enddo
     stop
  endif

  do i = 1, n
     pops(i) = max(1.d-30,matrixB(i,1))
  enddo

  deallocate(matrixA, matrixB, tempMatrix)

end subroutine solvePops

end module

