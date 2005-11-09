! photoionization module - started on October 4th 2005 by th

module photoion_mod

use source_mod
use grid_mod
use amr_mod
use constants_mod
implicit none

type SAHAMILNETABLE
   integer :: nFreq 
   integer :: nTemp 
   real(double), pointer :: temp(:)
   real(double), pointer :: freq(:)
   real(double), pointer :: Clyc(:, :)
end type SAHAMILNETABLE


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
    integer :: i, j
    logical :: escaped
    real(double) :: wavelength, thisFreq
    real :: thisLam
    type(OCTALVECTOR) :: outVec, octVec
    real(double) :: r
    integer :: ilam
    integer :: nInf
    real(double) :: kappaScadb, kappaAbsdb
    real(double) :: epsOverDeltaT, kappaH, kappaHe
    real :: dummy(3)
    real :: pops(10)
    integer :: nIter
    real(double) :: crate
    real :: xsec, temp, e, h0, he0
    real(double) :: luminosity1, luminosity2, luminosity3
    type(IONTYPE) :: thisIon
    real(double) :: alpha21s, alpha21p, alpha23s, photonPacketWeight
    real :: excitation, deexcitation, excitation2
    real :: fac
    type(SAHAMILNETABLE) :: hTable, heTable

    call createSahaMilneTables(hTable, heTable)

    lCore = 0.d0
    do i = 1, nSource
       lCore = lCore + source(i)%luminosity
    enddo

    write(*,'(a,1pe12.5)') "Total souce luminosity (lsol): ",lCore/lSol

    nMonte = 1000000

    nIter = 0

    do 
       nIter = nIter + 1
       nInf=0
       call plot_AMR_values(grid, "ionization", "x-z", 0., &
            "ionization.ps/vcps", .true., .false., &
            0, dummy, dummy, dummy, real(grid%octreeRoot%subcellsize), .false.)

       call plot_AMR_values(grid, "temperature", "x-z", real(grid%octreeRoot%centre%y), &
            "lucy_temp_xz.ps/vcps", .true., .false., &
            0, dummy, dummy, dummy, real(grid%octreeRoot%subcellsize), .false.) 
       epsoverdeltat = lcore/dble(nMonte)
       call dumpLexington(grid, epsoverdeltat)
       call zeroDistanceGrid(grid%octreeRoot)

       write(*,*) "Running loop with ",nmonte," photons."
       do iMonte = 1, nMonte

          call randomSource(source, nSource, iSource)
          thisSource = source(iSource)
          call getPhotonPositionDirection(thisSource, rVec, uHat,rHat)
          escaped = .false.
          call getWavelength(thisSource%spectrum, wavelength)
          photonPacketWeight = 1.d0
          thisFreq = cSpeed/(wavelength / 1.e8)

          do while(.not.escaped)

             call toNextEventPhoto(grid, rVec, uHat, escaped, thisFreq, nLambda, lamArray, photonPacketWeight)

             if (.not. escaped) then


                thisLam = (cSpeed / thisFreq) * 1.e8
                call locate(lamArray, nLambda, real(thisLam), iLam)
                octVec = rVec 
                call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, foundsubcell=subcell,iLambda=iLam, &
                     lambda=real(thisLam), kappaSca=kappaScadb, kappaAbs=kappaAbsdb, grid=grid)
                e = (hCgs * thisFreq)*ergtoev
                call phfit2(1, 1, 1 , e , h0)
                call phfit2(2, 2, 1 , e , he0)
                kappaH =  thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,1) * h0
                kappaHe = thisOctal%nh(subcell)*grid%ion(3)%abundance*thisOctal%ionFrac(subcell,3) * he0
                call random_number(r)


                if (r < (kappaH/(kappaH + kappaHe))) then

                   ! photon absorbed by H

                   alphaA = hRecombination(thisOctal%temperature(subcell))
                   alpha1 = recombToGround(thisOctal%temperature(subcell))
                   call random_number(r)
                   

                   if (r < (alpha1 / alphaA)) then                
                      ! recombination to continuum (equation 26 of Lucy MC trans prob II)
                      call getSahaMilneFreq(htable, dble(thisOctal%temperature(subcell)), thisFreq)
                      uHat =  randomUnitVector()
                   else
                      escaped = .true.
                   endif
                else
                   ! photon absorbed by He


                   call calcHeRecombs(thisOctal%temperature(subcell), alpha1, alpha21s, alpha21p, alpha23s)
                   call random_number(r)
                   alphaA = alpha1 + alpha21s + alpha21p + alpha23s
                   if ( (r <= alpha1/alphaA) ) then                   
                      ! recombination to continuum (equation 26 of Lucy MC trans prob II)
                      call getSahaMilneFreq(hetable, dble(thisOctal%temperature(subcell)), thisFreq)
                      uHat = randomUnitVector()
                   else if ((r > alpha1/alphaA).and.(r <= (alpha1+alpha21s)/alphaA)) then
                      !two photon continuum
                      call twoPhotonContinuum(thisFreq)
                      uHat =  randomUnitVector()
                            
                   else if ( (r > (alpha1+alpha21s)/alphaA).and.(r <= (alpha1+alpha21s+alpha23s)/alphaA)) then
                      ! ly alpha
                      thisFreq = (21.2 / ergtoev)/hcgs
                      uHat = randomUnitVector()
                   else 
                      thisFreq = (19.2 /ergtoev)/hcgs
                      uHat =  randomUnitVector()
                   endif
                endif
             endif
          enddo
          nInf = nInf + 1
       end do
       epsOverDeltaT = (lCore) / dble(nInf)

       write(*,*) "Calculating ionization and thermal equilibria"
       call calculateIonizationBalance(grid,grid%octreeRoot, epsOverDeltaT)
!       call calculateThermalBalance(grid, grid%octreeRoot, epsOverDeltaT)
       write(*,*) "Done."

       fac = 2.06e37

       call getForbiddenLineLuminosity(grid, "N II", 1.22e6, luminosity1)
       write(*,'(a,f12.4)') "N II (122 um):",(luminosity1+luminosity2)/fac

       call getForbiddenLineLuminosity(grid, "N II", 6584., luminosity1)
       call getForbiddenLineLuminosity(grid, "N II", 6548., luminosity2)
       write(*,'(a,f12.4)') "N II (6584+6548):",(luminosity1+luminosity2)/fac

       call getForbiddenLineLuminosity(grid, "N III", 5.73e5, luminosity1)
       write(*,'(a,f12.4)') "N III (57.3 um):",(luminosity1+luminosity2)/fac


       call getForbiddenLineLuminosity(grid, "O I", 6300., luminosity1)
       call getForbiddenLineLuminosity(grid, "O I", 6363., luminosity2)
       write(*,'(a,f12.4)') "O I (6300+6363):",(luminosity1+luminosity2)/fac
       call getForbiddenLineLuminosity(grid, "O II", 7320., luminosity1)
       call getForbiddenLineLuminosity(grid, "O II", 7330., luminosity2)
       write(*,'(a,f12.4)') "O II (7320+7330):",(luminosity1+luminosity2)/fac
       call getForbiddenLineLuminosity(grid, "O II", 3726., luminosity1)
       call getForbiddenLineLuminosity(grid, "O II", 3729., luminosity2)
       write(*,'(a,f12.4)') "O II (3726+3729):",(luminosity1+luminosity2)/fac
       call getForbiddenLineLuminosity(grid, "O III", 5007., luminosity1)
       call getForbiddenLineLuminosity(grid, "O III", 4959., luminosity2)
       write(*,'(a,f12.4)') "O III (5007+4959):",(luminosity1+luminosity2)/fac
       call getForbiddenLineLuminosity(grid, "O III", 4363., luminosity3)
       write(*,'(a,f12.4)') "O III (4363):",(luminosity3)/1.e37
       write(*,*) (luminosity1 + luminosity2)/luminosity3, 8.3*exp(3.3e4/7000.)
       call getForbiddenLineLuminosity(grid, "O III", 518145., luminosity1)
       call getForbiddenLineLuminosity(grid, "O III", 883562., luminosity2)
       write(*,'(a,f12.4)') "O III (52+88um):",(luminosity1+luminosity2)/fac


       call getForbiddenLineLuminosity(grid, "Ne II", 1.28e5, luminosity1)
       write(*,'(a,f12.4)') "Ne II (12.8um):",(luminosity1)/fac

       call getForbiddenLineLuminosity(grid, "Ne III", 1.56e5, luminosity1)
       write(*,'(a,f12.4)') "Ne III (15.5um):",(luminosity1)/fac

       call getForbiddenLineLuminosity(grid, "Ne III", 3869., luminosity1)
       call getForbiddenLineLuminosity(grid, "Ne III", 3968., luminosity2)
       write(*,'(a,f12.4)') "Ne III (3869+3968):",(luminosity1+luminosity2)/fac


       call getForbiddenLineLuminosity(grid, "S II", 6716., luminosity1)
       call getForbiddenLineLuminosity(grid, "S II", 6731., luminosity2)
       write(*,'(a,f12.4)') "S II (6716+6731):",(luminosity1+luminosity2)/fac

       call getForbiddenLineLuminosity(grid, "S II", 4068., luminosity1)
       call getForbiddenLineLuminosity(grid, "S II", 4076., luminosity2)
       write(*,'(a,f12.4)') "S II (4068+4076):",(luminosity1+luminosity2)/fac

       call getForbiddenLineLuminosity(grid, "S III", 1.87e5, luminosity1)
       write(*,'(a,f12.4)') "S III (18.7um):",(luminosity1)/fac

       call getForbiddenLineLuminosity(grid, "S III", 9532., luminosity1)
       call getForbiddenLineLuminosity(grid, "S III", 9069., luminosity2)
       write(*,'(a,f12.4)') "S III (9532+9069):",(luminosity1+luminosity2)/fac



    enddo

  end subroutine photoIonizationloop


 subroutine toNextEventPhoto(grid, rVec, uHat,  escaped,  thisFreq, nLambda, lamArray, photonPacketWeight)

   type(GRIDTYPE) :: grid
   type(OCTALVECTOR) :: rVec,uHat, octVec,thisOctVec, tvec
   type(OCTAL), pointer :: thisOctal, tempOctal, sourceOctal
   type(OCTAL),pointer :: oldOctal
   type(OCTAL),pointer :: foundOctal
   real(double) :: photonPacketWeight
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

       call updateGrid(grid, thisOctal, subcell, thisFreq, tVal, photonPacketWeight)
          

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

          call updateGrid(grid, thisOctal, subcell, thisFreq, dble(tval)*dble(tau)/thisTau, photonPacketWeight)

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
    integer :: subcell, i, j
    type(OCTALVECTOR) :: rVec
    logical :: converged, found
    real :: t1, t2, tm
    real(double) :: y1, y2, ym, ymin, Hheating, Heheating, tot
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
                totalHeating = (Hheating + HeHeating)

                nHii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,2) * grid%ion(2)%abundance
                nHeii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,4) * grid%ion(4)%abundance
                nh = thisOctal%nh(subcell)

                if (totalHeating == 0.d0) then
                   thisOctal%temperature(subcell) = 3.
                else
                   converged = .false.
!                   t1 = 1000.
!                   t2 = 100000.
!                   ymin = 1.e30
!                   found = .false.
!                   do i = 1, 1000
!                      t1 = 1000.+99000.*real(i-1)/1000.
!                      t2 = 1000.+99000.*real(i)/1000.
!                      y1 = (HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t1) - totalHeating)
!                      y2 = (HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t2) - totalHeating)
!                      if (abs(y1) < ymin) then
!                         tm = t1
!                         ymin = abs(y1)
!                      endif
!                      if (y1*y2 < 0.d0) then
!                         found = .true.
!                         exit
!                      endif
!                   enddo

                   t1 = 5000.
                   t2 = 20000.
                   found = .true.

                   if (found) then
                      y1 = (HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t1) - totalHeating)
                      y2 = (HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t2) - totalHeating)
                      if (y1*y2 > 0.d0) then

                         write(*,*) "t1", &
                              t1,y1,HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t1),totalHeating
                         write(*,*) "t2", &
                              t2,y2,HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t2),totalHeating
                         tm = 8000.
                         converged = .true.
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
                   endif
                   deltaT = tm - thisOctal%temperature(subcell)
                   thisOctal%temperature(subcell) = thisOctal%temperature(subcell) + underCorrection * deltaT
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
                               8.442e-12, 6.279e-12, 4.539e-12, 3.192e-12, 2.185e-12, 1.458e-12, 9.484e-13, 6.023e-13, 3.738e-13, &
                               2.268e-13, 1.348e-13, 7.859e-14, 4.499e-14 /)

    real :: rootTbetaHe(18) = (/ 8.347e-11, 7.889e-11, 7.430e-11, 6.971e-11, 6.512e-11, 6.056e-11, 5.603e-11, 5.154e-11, 4.710e-11,&
                               4.274e-11, 3.847e-11, 3.431e-11, 3.031e-11, 2.650e-11, 2.291e-11, 1.960e-11, 1.660e-11, 1.394e-11 /)

    real :: logT(31) = (/ 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, &
                      5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0 /)
    real(double) :: fac, thisRootTbetaH, betaH, betaHe, thisRootTbetaHe
    real :: thisLogT
    integer :: n
    
    coolingRate = 0.d0

    gff = 1.1d0 + 0.34d0*exp(-(5.5d0 - log10(temperature))**2 / 3.d0)  ! Kenny equation 23

    coolingRate = 1.42d-27 * (nHii+nHeii) * ne * gff * sqrt(temperature)  ! Kenny equation 22 (free-free)

    if (coolingRate < 0.) then
       write(*,*) "negative ff cooling",nhii,nheii,ne,gff,sqrt(temperature)
    endif

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

    if (ne * nhii * kerg * temperature * betaH < 0.) then
       write(*,*) "negative H cooling",ne,nhii,kerg,temperature,betah
    endif

    coolingRate = coolingRate + ne * nheii * kerg * temperature * betaHe

    if (ne * nheii * kerg * temperature * betaHe < 0.) then
       write(*,*) "negative He cooling",ne,nheii,kerg,temperature,betaHe
    endif

    call metalcoolingRate(grid%ion, grid%nIon, thisOctal, subcell, nh, ne, temperature, crate)
    if (crate < 0.) then
       write(*,*) "total negative metal cooling",crate
    endif
    coolingRate = coolingRate + crate


  end function HHeCooling

  subroutine updateGrid(grid, thisOctal, subcell, thisFreq, distance, photonPacketWeight)
    type(GRIDTYPE) :: grid
    type(OCTAL) :: thisOctal
    integer :: subcell
    real(double) :: thisFreq, distance
    real(double) :: photonPacketWeight
    integer :: i 
    real :: e, xsec

    thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1

    e = thisFreq * hCgs * ergtoEV

    do i = 1, grid%nIon
       call phfit2(grid%ion(i)%z, grid%ion(i)%n, grid%ion(i)%outerShell , e , xsec)
       if (xSec > 0.) then
          thisOctal%photoIonCoeff(subcell,i) = thisOctal%photoIonCoeff(subcell,i) &
               + distance * dble(xsec) / (dble(hCgs) * thisFreq) * photonPacketWeight
       endif

       ! neutral H heating

       if ((grid%Ion(i)%z == 1).and.(grid%Ion(i)%n == 1)) then
          thisOctal%Hheating(subcell) = thisOctal%Hheating(subcell) &
            + dble(distance) * dble(xsec / (thisFreq * hCgs)) &
            * (dble(hCgs * thisFreq) - dble(hCgs * grid%ion(i)%nuThresh)) * photonPacketWeight
       endif

       ! neutral He heating

       if ((grid%Ion(i)%z == 2).and.(grid%Ion(i)%n == 2)) then
          thisOctal%Heheating(subcell) = thisOctal%Heheating(subcell) &
            + dble(distance) * dble(xsec / (thisFreq * hCgs)) &
            * (dble(hCgs * thisFreq) - dble(hCgs * grid%ion(i)%nuThresh)) * photonPacketWeight
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
  real(double), allocatable :: matrixA(:,:), MatrixB(:,:), tempA(:,:)
  integer :: nIonizationStages
  type(OCTALVECTOR) :: rVec
  real(double) :: r1, r2
  integer :: iStart, iEnd, iIon
  real(double) :: chargeExchangeIonization, chargeExchangeRecombination
  logical :: ok
  real(double) :: ionRateInto, recombRateOutof

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

     allocate(matrixA(1:(nIonizationStages+1),1:(nIonizationStages+1)))
     allocate(tempA(1:(nIonizationStages+1),1:(nIonizationStages+1)))
     allocate(matrixB(1:(nIonizationStages+1), 1))
     
     matrixA = 0.d0

     do i = 1, nIonizationStages
        iIon = iStart+i-1
!        write(*,*) i,iion,grid%ion(iion)%species

        call getChargeExchangeRecomb(grid%ion(iion+1), thisOctal%temperature(subcell), &
             thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,1),  &
             thisOctal%nh(subcell)*grid%ion(2)%abundance*thisOctal%ionFrac(subcell,2),  &
             chargeExchangeRecombination)

        call getChargeExchangeIon(grid%ion(iion), thisOctal%temperature(subcell), &
             thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,1),  &
             thisOctal%nh(subcell)*grid%ion(2)%abundance*thisOctal%ionFrac(subcell,2),  &
             chargeExchangeIonization)


        if (i == 1) then
           matrixA(i, i) = -((epsOverDeltaT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell, iIon) + chargeExchangeIonization)! equation 44 Lucy MC II
           
           matrixA(i, i+1) = (recombRate(grid%ion(iIon),thisOctal%temperature(subcell)) * &
                thisOctal%ne(subcell) + chargeExchangeRecombination)
        else if (i == nIonizationStages) then
           matrixA(i,i-1) = IonRateInto
           matrixA(i,i) = -recombRateOutof
        else
!           matrixA(i,i-1) = ionRateInto
!           matrixA(i,i) = -recombRateOutof - &
!                ((epsOverDeltaT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell, iIon) + chargeExchangeIonization)
!           matrixA(i, i+1) = (recombRate(grid%ion(iIon),thisOctal%temperature(subcell)) * &
!                thisOctal%ne(subcell) + chargeExchangeRecombination)
           matrixA(i,i-1) = ionRateInto
           matrixA(i,i) = -recombRateOutof - &
                ((epsOverDeltaT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell, iIon) + chargeExchangeIonization)
           matrixA(i, i+1) = (recombRate(grid%ion(iIon),thisOctal%temperature(subcell)) * &
                thisOctal%ne(subcell) + chargeExchangeRecombination)
        endif

        
        if (i < nIonizationStages) then
           ionRateInto = ((epsOverDeltaT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell, iIon) + chargeExchangeIonization)
           recombRateOutOf =  (recombRate(grid%ion(iIon),thisOctal%temperature(subcell)) * &
                thisOctal%ne(subcell) + chargeExchangeRecombination)
        endif


     enddo

     where(matrixA < 0.d0)
        matrixA = min(matrixA,-1.d-30)
     end where

     where(matrixA >= 0.d0)
        matrixA = max(matrixA,1.d-30)
     end where

     do i = 1, nIonizationStages
        matrixA(nIonizationStages+1,i) = 1.d0
     enddo
     

     matrixB = 0.d0
     matrixB (nIonizationStages+1,1) = 1.d0


     tempA = matrixA
     if (grid%ion(istart)%species(1:1) =="x") then
        do i = 1, nIonizationStages+1
           write(*,'(1p,9e12.3)') matrixA(i,1:nIonizationStages)
        enddo
     endif



     call gaussj(matrixA, nIonizationStages+1, nIonizationStages+1, matrixB, 1, 1, ok)

     if (grid%ion(istart)%species(1:1) =="x") then
        write(*,*) " "
        write(*,'(1p,9e12.3)') matrixB(1:nIonizationStages,1)
        write(*,*) " "
     endif


     if (.not.ok) then
        write(*,*) "Ionization balance failed for element: ",trim(grid%ion(istart)%species(1:2))
        write(*,*) "temp",thisOctal%temperature(subcell)
        write(*,*) "ne",thisOctal%ne(subcell)
        do i = 1, nIonizationStages+1
           write(*,'(1p,9e12.3)') tempA(i,1:nIonizationStages)
        enddo
        write(*,*) "setting element to neutral"
        matrixB(1:nIonizationStages,1) = 0.d0
        matrixB(1,1) = 1.d0
     endif

     do i = 1, nIonizationStages
        iIon = istart+i-1
        thisOctal%ionFrac(subcell, iIon) = matrixB(i,1)
!        write(*,*) grid%ion(iIon)%species, thisOctal%ionFrac(subcell, iIon)
     enddo
!     write(*,*) " " 
     deallocate(matrixA, matrixB, tempA)
     k = iEnd + 1
  end do
  
  thisOctal%ne(subcell) = returnNe(thisOctal, subcell, grid%ion, grid%nion)

end subroutine solveIonizationBalance


function recombRate(thisIon, temperature) result (rate)
  type(IONTYPE) :: thisIon
  real :: temperature
  real(double) :: rate
  real(double) :: a, b, t0, t1, c, d, f, t

! recombinations INTO this species

  select case(thisIon%z)

     case(1)
        select case(thisIon%n)
           case(1) ! H I
              a = 7.982e-11
              b = 0.7480
              t0 = 3.148e0
              t1 = 7.036e5
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
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
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case(1) ! He II
              a = 1.891e-10
              b = 0.7524
              t0 = 9.370e00
              t1 = 2.774e6
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select
     case(6)
        select case(thisIon%n)
           case(6) ! C I 
              a = 0.0108
              b = -0.1075
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
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case(2) ! C V
              a = 2.765e-10
              b = 0.6858
              t0 = 1.535e2
              t1 = 2.556e7
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
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
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select
     case(8)
        select case(thisIon%n)
           case(8) ! O I 
              t = temperature / 10000.
              if (t < 2.) then
                 a = -0.0001
                 b = 0.0001
                 c = 0.0956
                 d = 0.0193
                 f = 0.4106
              else
                 a = 0.3715
                 b = -0.0293
                 c = -0.0597
                 d = 0.0678
                 f = 0.7993
              endif
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
              d = -0.7020
              f = 1.1899
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(5) ! O IV
              t = temperature / 10000.
              if (t < 2.) then
                 a = -0.3648
                 b =  7.2698
                 c =  17.2187
                 d =  9.8335
                 f = -0.0166
              else
                 a = -2.5053
                 b = 3.4903
                 c = 67.4128
                 d = -3.4450
                 f = 0.8501
              endif
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
              rate = 0. ! section 4 of above paper
           case(9) ! Ne II
              a = 0.0129
              b =-0.1779
              c = 0.9353
              d =-0.0682
              f = 0.4516
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(8) ! Ne III
              a = 3.6781
              b = 14.1481
              c = 17.1175
              d = -0.5017
              f = 0.2313
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(7) ! Ne IV
              a =-0.0254
              b = 5.5365
              c = 17.0727
              d = -0.7225
              f = 0.1702
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(6) ! Ne V
              a = -0.0141
              b = 33.8479
              c = 43.1608
              d =-1.6072
              f = 0.1942
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
           case(5) ! Ne VI
              a = 19.9280
              b = 235.0536
              c = 152.5096
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
              rate = 3.e-30 ! page 1344, para 2, kenny's paper
           case(14) ! S  III
              rate = 1.5e-11
           case(13) ! S IV
              rate = 2.5e-11
        end select
     case DEFAULT
        write(*,*) "No recombination rate for ",thisIon%species
        rate  = 0.
  end select
end function recombRate

real function vernerFerland(a, t, t0, t1, b)
  real(double):: a, t, t0, t1, b

! based on Verner and Ferland (1996) ApJS 103 467

  vernerFerland = a / (sqrt(t/t0) * (1.d0+sqrt(t/t0))**(1.d0-b) * (1.d0+sqrt(t/t1))**(1.d0+b))
end function vernerFerland



function nussbaumerStorey1983(t, a, b, c, d, f) result(rate)

! based on Nussbaumer & Storey 1983 AA 126 75
  real(double) :: t, a, b, c, d, f, rate
  rate = 1.d-12 * (a/t + b + c*t + d*t**2)*(t**(-1.5d0))*exp(-f/t)

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

subroutine dumpLexington(grid, epsoverdt)
  type(GRIDTYPE) :: grid
  type(OCTAL), pointer :: thisOctal
  integer :: subcell
  integer :: i, j
  real(double) :: r, theta
  real :: t,hi,hei,oii,oiii,cii,ciii,civ,nii,niii,niv,nei,neii,neiii,neiv
  real(double) :: oirate, oiirate, oiiirate, oivrate
  real(double) :: v, epsoverdt,r1,r2
  type(OCTALVECTOR) :: rvec
  type(OCTALVECTOR) :: octVec
  real :: fac

  open(20,file="lexington.dat",form="formatted",status="unknown")
  open(21,file="orates.dat",form="formatted",status="unknown")

  do i = 1, 50
     r = (1.+5.d0*dble(i-1)/49.d0)*pctocm/1.e10

     t=0;hi=0; hei=0;oii=0;oiii=0;cii=0;ciii=0;civ=0;nii=0;niii=0;niv=0;nei=0;neii=0;neiii=0;neiv=0

     oirate = 0; oiirate = 0; oiiirate = 0; oivrate = 0
     do j = 1, 100
        call random_number(theta)
        theta = theta * Pi
        
        octVec = OCTALVECTOR(r*sin(theta),0.d0,r*cos(theta))
        
        call amrgridvalues(grid%octreeRoot, octVec,  foundOctal=thisOctal, foundsubcell=subcell)
        if (thisOctal%threed) then
           v = thisOctal%subcellSize**3
        else
           rVec = subcellCentre(thisOctal,subcell)
           r1 = rVec%x-thisOctal%subcellSize/2.d0
           r2 = rVec%x+thisOctal%subcellSize/2.d0
           v = dble(pi) * (r2**2 - r1**2) * thisOctal%subcellSize
        endif
        
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
        NeI = NeI + thisOctal%ionfrac(subcell,returnIonNumber("Ne I", grid%ion, grid%nIon))
        NeII = NeII + thisOctal%ionfrac(subcell,returnIonNumber("Ne II", grid%ion, grid%nIon))
        NeIII = NeIII + thisOctal%ionfrac(subcell,returnIonNumber("Ne III", grid%ion, grid%nIon))
        NeIV = NeIV + thisOctal%ionfrac(subcell,returnIonNumber("Ne IV", grid%ion, grid%nIon))

!        fac = thisOctal%nh(subcell) * returnAbundance(8) !* thisOctal%ionfrac(subcell,returnIonNumber("O I", grid%ion, grid%nIon))
        fac = 1.
        oirate = oirate + &
             fac*((epsOverDT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell,returnIonNumber("O I", grid%ion, grid%nIon)))
        oiirate = oiirate + &
             fac*((epsOverDT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell,returnIonNumber("O II", grid%ion, grid%nIon)))
        oiiirate = oiiirate + &
             fac*((epsOverDT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell,returnIonNumber("O III", grid%ion, grid%nIon)))
!        oivrate = oivrate + &
!             fac*((epsOverDT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell,returnIonNumber("O IV", grid%ion, grid%nIon)))

        t  = t + thisOctal%temperature(subcell)
     enddo


     hi = hi / 100.; hei = hei/100.; oii = oii/100; oiii = oiii/100.; cii=cii/100.
     ciii = ciii/100; civ=civ/100.; nii =nii/100.; niii=niii/100.; niv=niv/100.
     nei=nei/100.;neii=neii/100.; neiii=neiii/100.; neiv=neiv/100.;t=t/100.

     oirate = oirate / 100.
     oiirate = oiirate / 100.
     oiiirate = oiiirate / 100.
     oivrate = oivrate / 100.

     hi = log10(max(hi, 1e-10))
     hei = log10(max(hei, 1e-10))
     oii = log10(max(oii, 1e-10))
     oiii = log10(max(oiii, 1e-10))
     cii = log10(max(cii, 1e-10))
     ciii = log10(max(ciii, 1e-10))
     civ = log10(max(civ, 1e-10))
     nii = log10(max(nii, 1e-10))
     niii = log10(max(niii, 1e-10))
     niv= log10(max(niv, 1e-10))
     nei = log10(max(nei, 1e-10))
     neii = log10(max(neii, 1e-10))
     neiii = log10(max(neiii, 1e-10))
     neiv = log10(max(neiv, 1e-10))


     write(21,'(f5.3,4e12.3)') r*1.e10/pctocm,oirate,oiirate,oiiirate,oivrate

     write(20,'(f5.3,f9.1,  14f8.3)') &
          r*1.e10/pctocm,t,hi,hei,oii,oiii,cii,ciii,civ,nii,niii,niv,nei,neii,neiii,neiv
  enddo
  close(20)
  close(21)
end subroutine dumpLexington

subroutine getChargeExchangeRecomb(parentIon, temperature, nHI, nHII, recombRate)
  type(IONTYPE) :: parentIon
  real(double) :: nHI, nHII,  recombRate
  real :: t4, a, b, c, d
  real :: temperature

  recombRate  = 0.d0

  select case(parentIon%z)
     case(7)
        select case(parentIon%n)
           case(4) ! N IV
              t4 = temperature / 1.e4
              a = 3.05e-10
              b = 0.60
              c = 2.65
              d = -0.93
              recombRate  = kingdonFerland96(t4, a, b, c, d)
        end select
     case(8)
        select case(parentIon%n)
           case(7) ! O II
              t4 = temperature / 1.e4
              a = 1.04e-9
              b = 3.15e-2
              c = -0.61
              d = -9.73
              recombRate  = kingdonFerland96(t4, a, b, c, d)
           case(6) ! OIII
              t4 = temperature / 1.e4
              a = 1.04e-9
              b = 0.27
              c = 2.02
              d = -5.92
              recombRate  = kingdonFerland96(t4, a, b, c, d)
       end select
   end select

  recombRate = recombRate * nhi

end subroutine getChargeExchangeRecomb

subroutine getChargeExchangeIon(parentIon, temperature, nHI, nHII, IonRate)
  type(IONTYPE) :: parentIon
  real(double) :: nHI, nHII, ionRate
  real :: temperature
  real :: t4, a, b, c, d

  IonRate = 0.d0

  select case(parentIon%z)
     case(7)
        select case(parentIon%n)
           case(7) ! N I charge exchange Ionization
              t4 = temperature / 1.e4
              a = 4.55e-12
              b = -0.29
              c = -0.92
              d = -8.38
              ionRate  = kingdonFerland96(t4, a, b, c, d)
        end select
     case(8)
        select case(parentIon%n)
           case(8) ! O I charge exchange ionization
              t4 = temperature / 1.e4
              a = 7.40e-11
              b = 0.47
              c = 24.37
              d = -0.74
              ionRate  = kingdonFerland96(t4, a, b, c, d)
       end select
   end select

  ionRate = ionRate * nhii
end subroutine getChargeExchangeIon


function kingdonFerland96(t4, a, b, c, d) result (alpha)
  real :: alpha
  real :: t4, a, b, c, d
  alpha = a*(t4**b)*(1.+c*exp(d*t4))
end function kingdonFerland96
  

subroutine getCollisionalRates(thisIon, iTransition, temperature, excitation, deexcitation)
  real :: excitation, deexcitation
  type(IONTYPE) :: thisIon
  integer :: iTransition, i
  real :: temperature
  real :: thisGamma
  real :: t , fac
  real :: boltzFac

  t = max(min(20000.,temperature),5000.)
  call locate(thisIon%transition(iTransition)%t, thisIon%transition(iTransition)%ngamma, t, i)
  fac = (t - thisIon%transition(iTransition)%t(i))/(thisIon%transition(iTransition)%t(i+1) - thisIon%transition(iTransition)%t(i))
  thisGamma = thisIon%transition(iTransition)%gamma(i) + &
       fac * (thisIon%transition(iTransition)%gamma(i+1) - thisIon%transition(iTransition)%gamma(i))

  boltzFac =  exp(-thisIon%transition(iTransition)%energy / (kev*temperature))

  fac = (8.63e-6 / sqrt(temperature)) * thisGamma

  deexcitation =  fac / thisIon%level(thisIon%transition(iTransition)%j)%g
  excitation =  fac / thisIon%level(thisIon%transition(iTransition)%i)%g * boltzFac

!  excitation = fac * boltzFac
!  deexcitation = fac * thisIon%level(thisIon%transition(iTransition)%j)%g &
!       / thisIon%level(thisIon%transition(iTransition)%i)%g

!  if (thisIon%species == "O III") then
!     if (iTransition ==1) then
!        write(*,*) thisGamma, boltzFac, fac, excitation, deexcitation
!     endif
!  endif


end subroutine getCollisionalRates


subroutine getForbiddenLineLuminosity(grid, species, wavelength, luminosity)
  type(GRIDTYPE) :: grid
  character(len=*) :: species
  real :: wavelength
  real(double) :: luminosity
  integer :: iIon, iTrans, i

  iTrans = 0
  iIon = returnIonNumber(species, grid%ion, grid%nIon)
  do i = 1, grid%ion(iIon)%nTransitions
     if (abs(grid%ion(iIon)%transition(i)%lambda-wavelength)/wavelength  < 0.001) then
        iTrans = i
        exit
     endif
  enddo
  if (iTrans == 0) then
     write(*,*) "No transition found at ",wavelength, "Angstroms"
     stop
  endif
  luminosity = 0.d0
  call sumLineLuminosity(grid%octreeroot, luminosity, iIon, iTrans, grid)
end subroutine getForbiddenLineLuminosity

recursive subroutine sumLineLuminosity(thisOctal, luminosity, iIon, iTrans, grid)
  type(GRIDTYPE) :: grid
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i, iIon, iTrans
  real(double) :: luminosity, v, r1, r2, rate
  real :: pops(10)
  type(OCTALVECTOR) :: rvec
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call sumLineLuminosity(child, luminosity, iIon, iTrans, grid)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal,subcell)
          if (thisOctal%threed) then
             v = thisOctal%subcellSize**3
          else
             r1 = rVec%x-thisOctal%subcellSize/2.d0
             r2 = rVec%x+thisOctal%subcellSize/2.d0
             v = dble(pi) * (r2**2 - r1**2) * thisOctal%subcellSize
          endif
          call solvePops(grid%ion(iIon), pops, thisOctal%ne(subcell), thisOctal%temperature(subcell))
          rate =  pops(grid%ion(iion)%transition(iTrans)%j) * grid%ion(iion)%transition(itrans)%energy * &
               grid%ion(iion)%transition(itrans)%a/ergtoev
          rate = rate * grid%ion(iion)%abundance * thisOctal%nh(subcell) * thisOctal%ionFrac(subcell, iion)
          luminosity = luminosity + rate * v * 1.d30
          
       endif
    enddo
  end subroutine sumLineLuminosity


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
        if (rate < 0.) then
           write(100,'(a,a,a,1p,e12.4,a,0p, f10.1,a,1pe12.4)') "Negative contribution from ", &
                trim(ionArray(j)%species),":",rate," at T = ",temperature, &
                " ion frac ",thisOctal%ionFrac(subcell,j)
        endif
        if (present(debug)) then
           if (debug) then
                 write(100,'(a,a,a,1p,e12.4,a,0p, f10.1,a,1pe12.4)') "Contribution from ", &
                      trim(ionArray(j)%species),":",rate," at T = ",temperature, &
                      " ion frac ",thisOctal%ionFrac(subcell,j)
           endif
        endif
     total = total + rate
     endif
  enddo
end subroutine metalcoolingRate
  

subroutine solvePops(thisIon, pops, ne, temperature, debug)
  type(IONTYPE) :: thisIon
  real(double) :: ne
  real :: pops(*)
  real :: temperature
  real(double), allocatable :: matrixA(:,:), MatrixB(:,:), tempMatrix(:,:)
  integer :: n, iTrans, i, j
  real :: excitation, deexcitation, rateij, rateji
  logical :: ok
  logical, optional :: debug

  n = thisIon%nLevels
  allocate(matrixA(1:n+1, 1:n+1), matrixB(1:n+1, 1), tempMatrix(1:n+1,1:n+1))

  matrixA = 1.d-30
  matrixB = 0.d0
  do iTrans = 1, thisIon%nTransitions
     i = thisIon%transition(iTrans)%i
     j = thisIon%transition(iTrans)%j
     call getCollisionalRates(thisIon, iTrans, temperature, excitation, deexcitation)
     rateij = max(1.e-40,excitation * ne)
     rateji = max(1.e-30,deexcitation * ne + thisIon%transition(iTrans)%a)

     if (PRESENT(debug)) then
        if (debug) write(*,*) i, j, rateij, rateji
     endif

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

  if (PRESENT(debug)) then
     if (debug) then
        do i = 1, n+1
           write(*,'(1p,9e12.3)') tempmatrix(i,1:n)
        enddo
     endif
  endif

  call gaussj(matrixA, n+1, n+1, matrixB, 1, 1, ok)
  if (.not.ok) then
     write(*,*) "Population solver failed for: ",thisIon%species
     write(*,*) "nlevels",thisIon%nLevels,"ntrans",thisIon%nTransitions 
     write(*,*) "temp",temperature,"ne",ne
     
     do i = 1, n+1
        write(*,'(1p,9e12.3)') tempmatrix(i,1:n)
     enddo
     matrixB = 0.d0
     matrixB(1,1) = 1.d0
     write(*,*) "Setting pops to ground state"
  endif

  do i = 1, n
     pops(i) = max(1.d-30,matrixB(i,1))
  enddo

  deallocate(matrixA, matrixB, tempMatrix)

end subroutine solvePops

subroutine calcHeRecombs(te, alpha1, alpha21s, alpha21p, alpha23s)
  real :: te, t
  real(double) :: alpha1, alpha21s, alpha21p, alpha23s

  t = te / 1.e4

  alpha1 = 1.54e-13 * t**(-0.486)

  alpha21s = 2.06e-14 * t**(-0.676)

  alpha21p = 4.17e-14 * t**(-0.861)
  
  alpha23s = 2.10e-13 * t**(-0.778)
end subroutine calcHeRecombs


subroutine createSahaMilneTables(hTable, heTable)
  type(SAHAMILNETABLE) :: hTable, heTable
  real(double) :: nu0_h, nu0_he, nufinal
  integer :: nFreq, nTemp
  integer :: i, j
  real :: e
  real(double) :: t1, t2
  real :: hxsec, hexsec
  real(double) :: dfreq, jnu

  nFreq = 100
  nTemp = 100
  hTable%nFreq = nFreq
  heTable%nFreq = nFreq
  hTable%nTemp = nTemp
  heTable%nTemp = nTemp

  allocate(hTable%freq(1:nFreq), heTable%freq(1:nFreq))
  allocate(hTable%temp(1:ntemp), heTable%temp(1:ntemp))
  allocate(hTable%Clyc(1:nTemp,1:nFreq))
  allocate(heTable%Clyc(1:nTemp,1:nFreq))

  nu0_h = 13.6d0/ergtoev/hcgs
  nu0_he = 24.59d0/ergtoev/hcgs

  nufinal = 6.d0*nu0_he 

  do i = 1, nFreq
     hTable%freq(i) = log10(nu0_h) + (log10(nuFinal)-log10(nu0_h))*dble(i-1)/dble(nFreq-1)
     heTable%freq(i) = log10(nu0_he) + (log10(nuFinal)-log10(nu0_he))*dble(i-1)/dble(nFreq-1)
  enddo
  hTable%freq(1:hTable%nFreq) = 10.d0**hTable%freq(1:hTable%nfreq)
  heTable%freq(1:hTable%nFreq) = 10.d0**heTable%freq(1:hTable%nfreq)

  t1 = 5000.d0
  t2 = 20000.d0

  do i = 1, nTemp
     hTable%temp(i) = t1 + (t2-t1)*dble(i-1)/dble(nTemp-1)
     heTable%temp(i) = t1 + (t2-t1)*dble(i-1)/dble(nTemp-1)
  enddo

  do i = 1, nTemp
     hTable%Clyc(i,1) = 0.d0
     heTable%Clyc(i,1) = 0.d0
     do j = 2, nFreq


        e = hTable%freq(j) * hcgs* ergtoev
        call phfit2(1, 1, 1 , e , hxsec)

        dFreq = hTable%freq(j)-hTable%freq(j-1)
        jnu = (hTable%freq(j)/cSpeed)**2 * hTable%temp(i)**(-1.5d0) * dble(hxsec) &
             *  exp(-hcgs*(hTable%freq(j)-nu0_h)/(kerg*hTable%temp(i)))
        hTable%Clyc(i,j) = hTable%Clyc(i,j-1) + jnu * dfreq
!        write(*,*) i,j,htable%clyc(i,j),(hTable%freq(j)/cSpeed)**2,hTable%temp(i)**1.5d0,dble(hxsec), &
!             exp(-hcgs*(hTable%freq(j)-nu0_h)/(kerg*hTable%temp(i)))
        e = heTable%freq(j) * hcgs * ergtoev
        call phfit2(2, 2, 1 , e , hexsec)

        dFreq = heTable%freq(j)-heTable%freq(j-1)
        jnu = 2.d0*heTable%freq(j)**2 * heTable%temp(i)**(-1.5d0) * dble(hexsec) &
             * exp(-hcgs*(heTable%freq(j)-nu0_h)/(kerg*heTable%temp(i)))
        heTable%Clyc(i,j) = heTable%Clyc(i,j-1) + jnu * dfreq
     enddo
     hTable%Clyc(i,1:hTable%nFreq) = hTable%Clyc(i,1:hTable%nFreq) / hTable%Clyc(i,hTable%nFreq)
     heTable%Clyc(i,1:heTable%nFreq) = heTable%Clyc(i,1:heTable%nFreq) / hTable%Clyc(i,heTable%nFreq)
  end do
end subroutine createSahaMilneTables

subroutine getSahaMilneFreq(table,temperature, thisFreq)
  type(SAHAMILNETABLE) :: table
  real(double) :: temperature, thisfreq, r, t, fac
  integer :: i, j

  t = max(5000.d0, min(20000.d0, t))
  call locate(table%temp, table%nTemp, t, i)
  call random_number(r)
  call locate(table%Clyc(i,1:table%nfreq), table%nFreq, r, j)
  fac = (r - table%Clyc(i,j))/(table%Clyc(i,j+1)-table%cLyc(i,j))
  thisFreq = table%freq(j) + fac * (table%freq(j+1)-table%freq(j))
end subroutine getSahaMilneFreq

subroutine twoPhotonContinuum(thisFreq)

! based on table ii of drake, victor, dalgarno, 1969, PhyRev Vol 180, pg 25

  real(double) :: thisFreq
  real :: y(21) = (/ 0., 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.275, 0.300, &
       0.325, 0.350, 0.375, 0.400, 0.425, 0.450, 0.475, 0.500 /)
  real :: hei(21) = (/ 0., 7.77e0, 2.52e1, 4.35e1, 5.99e1, 7.42e1, 8.64e1, 9.69e1, 1.06e2, 1.13e2, 1.20e2, 1.25e2, &
       1.30e2, 1.34e2, 1.37e2, 1.40e2, 1.42e2, 1.43e2, 1.45e2, 1.45e2, 1.45e2 /)
  real :: freq = 3.86e15, fac, r
  real :: prob(21)
  integer :: i
  prob(1) = 0.
  do i = 2, 21
     prob(i) = prob(i-1) + (y(i)-y(i-1)) * hei(i)
  enddo
  prob(1:21) = prob(1:21)/prob(21)
  call random_number(r)
  call locate(prob, 21, r, i)
  fac = y(i) + ((r - prob(i))/(prob(i+1)-prob(i)))*(y(i+1)-y(i))
  thisFreq = (1.-fac)*freq
end subroutine twoPhotonContinuum
  

end module

