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
   real(double), pointer :: emissivity(:)
end type SAHAMILNETABLE

type RECOMBTABLE
   integer :: nrho
   integer :: nTemp 
   real(double), pointer :: rho(:)
   real(double), pointer :: temp(:)
   real(double), pointer :: emissivity(:, :)
end type RECOMBTABLE


type GAMMATABLE
   integer :: ntemp
   real(double), pointer :: temp(:)
   real(double), pointer :: gamma(:)
end type GAMMATABLE

contains


  subroutine photoIonizationloop(grid, source, nSource, nLambda, lamArray)
    use input_variables, only : smoothFactor
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
    logical :: gridConverged
    real(double) :: probEsc

    type(SAHAMILNETABLE) :: hTable, heTable
    type(RECOMBTABLE) :: Hrecombtable
    type(GAMMATABLE) :: HIgammaTable

    real(double) :: hRecombemissivity, lymanContemissivity
    real(double) :: hiContEmissivity, forbiddenEmissivity
    real(double) :: probNonIonizing

    call createSahaMilneTables(hTable, heTable)

    call createRecombTable(Hrecombtable, 'e1b.d')

    call createHIGammaTable(HIgammaTable)


    lCore = 0.d0
    do i = 1, nSource
       lCore = lCore + source(i)%luminosity
    enddo

    write(*,'(a,1pe12.5)') "Total souce luminosity (lsol): ",lCore/lSol

    nMonte = 100000

    nIter = 0

    do 
       nIter = nIter + 1
       nInf=0
       call plot_AMR_values(grid, "ionization", "x-z", 0., &
            "ionization.ps/vcps", .true., .true., &
            0, dummy, dummy, dummy, real(grid%octreeRoot%subcellsize), .false.)

       call plot_AMR_values(grid, "temperature", "x-z", real(grid%octreeRoot%centre%y), &
            "lucy_temp_xz.ps/vcps", .false., .false., &
            0, dummy, dummy, dummy, real(grid%octreeRoot%subcellsize), .false.) 
       epsoverdeltat = lcore/dble(nMonte)
       call zeroDistanceGrid(grid%octreeRoot)

       write(*,*) "Running loop with ",nmonte," photons."
       mainloop: do iMonte = 1, nMonte

          call randomSource(source, nSource, iSource)
          thisSource = source(iSource)
          call getPhotonPositionDirection(thisSource, rVec, uHat,rHat)
          escaped = .false.
          call getWavelength(thisSource%spectrum, wavelength)
          photonPacketWeight = 1.d0
          thisFreq = cSpeed/(wavelength / 1.e8)
          if ((thisFreq*hcgs*ergtoev) < 13.6) then
             cycle mainloop
          endif
          do while(.not.escaped)

             call toNextEventPhoto(grid, rVec, uHat, escaped, thisFreq, nLambda, lamArray, photonPacketWeight)

             if (.not. escaped) then


                thisLam = (cSpeed / thisFreq) * 1.e8
                call locate(lamArray, nLambda, real(thisLam), iLam)
                octVec = rVec 
                call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, foundsubcell=subcell,iLambda=iLam, &
                     lambda=real(thisLam), kappaSca=kappaScadb, kappaAbs=kappaAbsdb, grid=grid)

!                do i = 1, 10
!                   thisOctal%temperature(subcell) = 8000.d0 + 4000.d0*real(i-1)/9.
                call calcEmissivities(grid, thisOctal, subcell, hTable, higammatable, hRecombTable, &
                     hRecombemissivity, lymanContemissivity, hiContEmissivity, forbiddenEmissivity)
                hRecombemissivity = hRecombemissivity / fourpi
                hiContEmissivity = hiContemissivity / fourPi

                probNonIonizing = (hRecombEmissivity + hiContEmissivity + forbiddenEmissivity)  / &
                     (hRecombEmissivity + hiContEmissivity + lymanContEmissivity + forbiddenEmissivity)

!                write(*,'(f7.1,1p,5e12.3,0p)') thisOctal%temperature(subcell), hRecombEmissivity,lymancontemissivity,&
!                     hicontemissivity,forbiddenemissivity,probNonIonizing
!                enddo
!                stop


                call random_number(r)
                if (r < probNonIonizing) then
                   escaped = .true.
                else
                   call getSahaMilneFreq(htable,dble(thisOctal%temperature(subcell)), thisFreq)
!                   thisfreq = (13.6/ergtoev)/hcgs
!                   write(*,*) hRecombemissivity,lymancontEmissivity, hicontemissivity,probNonionizing
!                   write(*,*) "forbidden",forbiddenEmissivity
                endif
             endif
          enddo
          nInf = nInf + 1
       end do mainloop
       epsOverDeltaT = (lCore) / dble(nMonte)


       do i = 1, 3
          write(*,*) "Calculating ionization and thermal equilibria",i
          call calculateIonizationBalance(grid,grid%octreeRoot, epsOverDeltaT)
!          call calculateThermalBalance(grid, grid%octreeRoot, epsOverDeltaT)
          write(*,*) "Done."
       enddo

       call dumpLexington(grid, epsoverdeltat)
       fac = 2.06e37
       
       luminosity1 = 0.d0
       call getHbetaLuminosity(grid%octreeRoot, grid, luminosity1)
       write(*,'(a,2f12.4)') "H beta :",luminosity1/1.e37,luminosity1/2.05e37
       fac = luminosity1
       

       call getForbiddenLineLuminosity(grid, "N II", 1.22e6, luminosity1)
       write(*,'(a,2f12.4)') "N II (122 um):",(luminosity1)/fac,(luminosity1)/(0.034*2.05e37)

       call getForbiddenLineLuminosity(grid, "N II", 6584., luminosity1)
       call getForbiddenLineLuminosity(grid, "N II", 6548., luminosity2)
       write(*,'(a,2f12.4)') "N II (6584+6548):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.730*2.05e37)

       call getForbiddenLineLuminosity(grid, "N III", 5.73e5, luminosity1)
       write(*,'(a,2f12.4)') "N III (57.3 um):",(luminosity1)/fac,(luminosity1+luminosity2)/(0.292*2.05e37)


       call getForbiddenLineLuminosity(grid, "O I", 6300., luminosity1)
       call getForbiddenLineLuminosity(grid, "O I", 6363., luminosity2)
       write(*,'(a,2f12.4)') "O I (6300+6363):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.0086*2.05e37)

       call getForbiddenLineLuminosity(grid, "O II", 7320., luminosity1)
       call getForbiddenLineLuminosity(grid, "O II", 7330., luminosity2)
       write(*,'(a,2f12.4)') "O II (7320+7330):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.029*2.05e37)

       call getForbiddenLineLuminosity(grid, "O II", 3726., luminosity1)
       call getForbiddenLineLuminosity(grid, "O II", 3729., luminosity2)
       write(*,'(a,2f12.4)') "O II (3726+3729):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(2.03*2.05e37)

       call getForbiddenLineLuminosity(grid, "O III", 5007., luminosity1)
       call getForbiddenLineLuminosity(grid, "O III", 4959., luminosity2)
       write(*,'(a,2f12.4)') "O III (5007+4959):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(2.18*2.05e37)

       call getForbiddenLineLuminosity(grid, "O III", 4363., luminosity3)
       write(*,'(a,2f12.4)') "O III (4363):",(luminosity3)/1.e37,luminosity3/(0.0037*2.05e37)

       call getForbiddenLineLuminosity(grid, "O III", 518145., luminosity1)
       call getForbiddenLineLuminosity(grid, "O III", 883562., luminosity2)
       write(*,'(a,2f12.4)') "O III (52+88um):,",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/((1.06+1.22)*2.05e37)


       call getForbiddenLineLuminosity(grid, "Ne II", 1.28e5, luminosity1)
       write(*,'(a,2f12.4)') "Ne II (12.8um):",(luminosity1)/fac,luminosity1/(0.195*2.05e37)

       call getForbiddenLineLuminosity(grid, "Ne III", 1.56e5, luminosity1)
       write(*,'(a,2f12.4)') "Ne III (15.5um):",(luminosity1)/fac,luminosity1/(0.322*2.05e37)

       call getForbiddenLineLuminosity(grid, "Ne III", 3869., luminosity1)
       call getForbiddenLineLuminosity(grid, "Ne III", 3968., luminosity2)
       write(*,'(a,2f12.4)') "Ne III (3869+3968):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.085*2.05e37)


       call getForbiddenLineLuminosity(grid, "S II", 6716., luminosity1)
       call getForbiddenLineLuminosity(grid, "S II", 6731., luminosity2)
       write(*,'(a,2f12.4)') "S II (6716+6731):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.147*2.05e37)

       call getForbiddenLineLuminosity(grid, "S II", 4068., luminosity1)
       call getForbiddenLineLuminosity(grid, "S II", 4076., luminosity2)
       write(*,'(a,2f12.4)') "S II (4068+4076):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.008*2.05e37)

       call getForbiddenLineLuminosity(grid, "S III", 1.87e5, luminosity1)
       write(*,'(a,2f12.4)') "S III (18.7um):",(luminosity1)/fac,luminosity1/(0.577*2.05e37)

       call getForbiddenLineLuminosity(grid, "S III", 9532., luminosity1)
       call getForbiddenLineLuminosity(grid, "S III", 9069., luminosity2)
       write(*,'(a,2f12.4)') "S III (9532+9069):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(1.22*2.05e37)



!       if ((niter > 2).and.(nIter < 8)) then
!          call locate(grid%lamArray, grid%nLambda,900.,ilam)
!          if (writeoutput) write(*,*) "Smoothing adaptive grid structure for ionization..."
!          call smoothAMRgridIonization(grid%octreeRoot,grid,gridConverged,ilam,inheritprops = .false., interpProps = .false.)
!          if (writeoutput) write(*,*) "...grid smoothing complete"
!          call locate(grid%lamArray, grid%nLambda,400.,ilam)
!          if (writeoutput) write(*,*) "Smoothing adaptive grid structure for ionization..."
!          call smoothAMRgridIonization(grid%octreeRoot,grid,gridConverged,ilam,inheritprops = .false., interpProps = .false.)
!          if (writeoutput) write(*,*) "...grid smoothing complete"
!          if (writeoutput) write(*,*) "Smoothing adaptive grid structure..."
!          gridConverged = .false.
!          do
!             call smoothAMRgrid(grid%octreeRoot,grid,smoothFactor,gridConverged,inheritprops=.false., interpProps=.false.)
!             if (gridConverged) exit
!          end do
!          if (writeoutput) write(*,*) "...grid smoothing complete"
!       endif

       nMonte = nMonte * 2

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


    call amrGridValues(grid%octreeRoot, octVec,  iLambda=iLam, lambda=real(thisLam), foundOctal=thisOctal, &
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
          call amrGridValues(grid%octreeRoot, octVec, iLambda=iLam,lambda=real(thisLam), foundOctal=thisOctal, &
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

       call amrGridValues(grid%octreeRoot, octVec, startOctal=oldOctal, iLambda=iLam,lambda=real(thisLam),&
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
          thisOCtal%nDiffusion(subcell) = 0.
          thisOctal%Hheating(subcell) = 0.d0
          thisOctal%Heheating(subcell) = 0.d0
          thisOctal%photoIonCoeff(subcell,:) = 0.d0
       endif
    enddo
  end subroutine zeroDistanceGrid

  recursive subroutine getHbetaluminosity(thisOctal, grid, luminosity)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  type(GRIDTYPE) :: grid
  real(double) :: luminosity, v, r1, r2, hbeta
  type(OCTALVECTOR) :: rVec
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call getHbetaLuminosity(child, grid, luminosity)
                exit
             end if
          end do
       else
          if (thisOctal%threed) then
             v = thisOctal%subcellSize**3
          else
             rVec = subcellCentre(thisOctal,subcell)
             r1 = rVec%x-thisOctal%subcellSize/2.d0
             r2 = rVec%x+thisOctal%subcellSize/2.d0
             v = dble(pi) * (r2**2 - r1**2) * thisOctal%subcellSize
          endif
          hbeta = (2530.d0/thisOctal%temperature(subcell)**0.833) * thisOctal%ne(subcell) * thisOctal%ionFrac(subcell, 2) * &
               thisOctal%nh(subcell)*grid%ion(1)%abundance*1.d-25
          luminosity = luminosity + hbeta * (v*1.d30)
       endif
    enddo
  end subroutine getHbetaluminosity

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
    real :: underCorrection = 1.
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


                call getHeating(grid, thisOctal, subcell, hHeating, heHeating, totalHeating, epsOverDeltaT)

                nHii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,2) * grid%ion(2)%abundance
                nHeii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,4) * grid%ion(4)%abundance
                nh = thisOctal%nh(subcell)
                

                if (totalHeating == 0.d0) then
                   thisOctal%temperature(subcell) = 1.e-3
                else
                   converged = .false.

                   t1 = 5000.
                   t2 = 20000.
                   found = .true.

                   if (found) then
                      y1 = (HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t1) - totalHeating)
                      y2 = (HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t2) - totalHeating)
                      if (y1*y2 > 0.d0) then

!                         write(*,*) "t1", &
!                              t1,y1,HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t1),totalHeating
!                         write(*,*) "t2", &
!                              t2,y2,HHecooling(grid, thisOctal, subcell, nh, nhii,nheii, thisOctal%ne(subcell), t2),totalHeating
                         tm = 5000.
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
                            write(*,*) t1, t2, y1,y2,ym
                         endif
                         
                         if (abs((t1-t2)/t1) .le. 1.e-4) then
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
    real(double) :: log10te,betarec, coolrec, betaff, coolff
    real :: ch12, ch13, ex12, ex13, th12, th13, coolcoll, te4, teused
    real :: becool
    real, parameter                :: hcRyd = &    ! constant: h*c*Ryd (Ryd at inf used) [erg]
         & 2.1799153e-11


    becool = 0.
    coolingRate = 0.d0

    gff = 1.1d0 + 0.34d0*exp(-(5.5d0 - log10(temperature))**2 / 3.d0)  ! Kenny equation 23

    coolingRate = 1.42d-27 * (nHii+nHeii) * ne * gff * sqrt(temperature)  ! Kenny equation 22 (free-free)

    thisLogT = log10(temperature)
    log10te = thisLogT
    teused = temperature


    ! cooling of gas due to FF radiation from H+ 
    ! fits to Hummer, MNRAS 268(1994) 109, Table 1. or  least square fitting to m=4

    betaFF = 1.0108464E-11 + 9.7930778E-13*log10Te - &
         & 6.6433144E-13*log10Te*log10Te + 2.4793747E-13*log10Te*log10Te*log10Te -&
         & 2.3938215E-14*log10Te*log10Te*log10Te*log10Te
    
    coolFF = Nhii*Ne*betaFF*kerg*TeUsed/sqrt(TeUsed)

    becool = becool + coolfF

    ! collisional excitation of Hydrogen
    ! Mathis, Ly alpha, beta
    ch12 = 2.47e-8
    ch13 = 1.32e-8
    ex12 = -0.228
    ex13 = -0.460
    th12 = 118338.
    th13 = 140252.
    te4 = temperature / 1.e4
    if (TeUsed > 5000.) then 
       
       coolColl = (ch12*exp(-th12/TeUsed)*Te4**ex12 + &
            & ch13*exp(-th13/TeUsed)*Te4**ex13) * &
            & hcRyd*(nh-nhii)*Ne
       
    else
       
       coolColl = 0.
       
    end if
    coolingrate = coolingrate + coolcoll

    becool = becool +coolcoll

    if (coolingRate < 0.) then
       write(*,*) "negative ff cooling",nhii,nheii,ne,gff,sqrt(temperature)
    endif

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


            ! cooling of gas due to recombination of H+
            ! fits to Hummer, MNRAS 268(1994) 109, Table 1.  
            ! least square fitting to m=4
            betaRec = 9.4255985E-11 -4.04794384E-12*log10Te &
                 & -1.0055237E-11*log10Te*log10Te +  1.99266862E-12*log10Te*log10Te*log10Te&
                 & -1.06681387E-13*log10Te*log10Te*log10Te*log10Te

            coolRec = Nhii*Ne*betaRec*kerg*TeUsed/sqrt(TeUsed)

            becool = becool+ coolrec

    coolingRate = coolingRate +  ne * nhii * kerg * temperature * betaH



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
  logical :: ok, easyWay
  real(double) :: ionRateInto, recombRateOutof
  real(double), allocatable :: xplus1overx(:)

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

     
     easyWay = .true.

     if (easyWay) then

        allocate(xplus1overx(1:nIonizationStages-1))
        do i = 1, nIonizationStages-1
           iIon = iStart+i-1
           call getChargeExchangeRecomb(grid%ion(iion+1), thisOctal%temperature(subcell), &
                thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,1),  &
                thisOctal%nh(subcell)*grid%ion(2)%abundance*thisOctal%ionFrac(subcell,2),  &
                chargeExchangeRecombination)

           call getChargeExchangeIon(grid%ion(iion), thisOctal%temperature(subcell), &
                thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,1),  &
                thisOctal%nh(subcell)*grid%ion(2)%abundance*thisOctal%ionFrac(subcell,2),  &
                chargeExchangeIonization)

        
           xplus1overx(i) = ((epsOverDeltaT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell, iIon) + chargeExchangeIonization) / &
      max(1.d-30,(recombRate(grid%ion(iIon),thisOctal%temperature(subcell)) * thisOctal%ne(subcell) + chargeExchangeRecombination))
!           if (grid%ion(iion)%species(1:1) =="O") write(*,*) i,xplus1overx(i)
        enddo
       thisOctal%ionFrac(subcell, iStart:iEnd) = 1.
        do i = 1, nIonizationStages - 1
           iIon = iStart+i-1
           thisOctal%ionFrac(subcell,iIon+1) = thisOctal%ionFrac(subcell,iIon) * xplus1overx(i)
!           if (grid%ion(iion)%species(1:1) =="O") write(*,*) i,thisOctal%ionFrac(subcell,iIon)
        enddo
        if (SUM(thisOctal%ionFrac(subcell,iStart:iEnd)) /= 0.d0) then
           thisOctal%ionFrac(subcell,iStart:iEnd) = &
                max(1.d-30,thisOctal%ionFrac(subcell,iStart:iEnd))/SUM(thisOctal%ionFrac(subcell,iStart:iEnd))
        else
           thisOctal%ionFrac(subcell,iStart:iEnd) = 1.e-30
        endif
           
!        if (grid%ion(iion)%species(1:1) =="O") then
!           write(*,*) thisOctal%ionFrac(subcell,iStart:iEnd)
!           write(*,*) " "
!        endif

        deallocate(xplus1overx)

     else

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

  endif
  k = iEnd + 1
end do
  



thisOctal%ne(subcell) = returnNe(thisOctal, subcell, grid%ion, grid%nion)

end subroutine solveIonizationBalance


function recombRate(thisIon, temperature) result (rate)
  type(IONTYPE) :: thisIon
  real :: temperature
  real(double) :: rate
  real(double) :: a, b, t0, t1, c, d, f, t, z

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
              z = 1
              a = 5.068
              b = -0.6192
              c = -0.0815
              d = 1.2910
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(5) ! C II
              a = 1.8267
              b = 4.1012
              c = 4.8443
              d = 0.2261
              f = 0.5960
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 2
              a = 5.434
              b = -0.6116
              c = 0.0694
              d = 0.7866
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(4) ! C III
              a = 2.3196
              b = 10.7328
              c = 6.8830
              d = -0.1824
              f = 0.4101
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 3
              a = 4.742
              b = -0.6167
              c = 0.2960
              d = 0.6167
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
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
              z = 1
              a = 3.874
              b = -0.6487
              c = 0.
              d = 1.
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(6) ! N II
              a = 0.0320
              b = -0.6624
              c = 4.3191
              d = 0.0003
              f = 0.5946
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 2
              a = 4.974
              b = -0.6209
              c = 0.
              d = 1.
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(5) ! N III
              a = -0.8806
              b = 11.2406
              c = 30.7066
              d = -1.1721
              f = 0.6127
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 3
              a = 4.750
              b = -0.5942
              c = 0.8452
              d = 2.8450
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(4) ! N IV
              a = 0.4134
              b = -4.6319
              c = 25.9172
              d = -2.2290
              f = 0.2360
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 4
              a = 4.626
              b = -0.9521
              c = 0.4729
              d = -0.4477
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
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
              z = 1
              a = 3.201
              b = -0.6880
              c = -0.0174
              d = 1.7070
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(7) ! O II
              a = -0.0036
              b = 0.7519
              c = 1.5252
              d =-0.0838
              f = 0.2769
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 2
              a = 4.092
              b = -0.6413
              c = 0.
              d = 1.
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(6) ! O III
              a = 0.0
              b = 21.8790
              c = 16.2730
              d = -0.7020
              f = 1.1899
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 3
              a = 4.890
              b = -0.6213
              c = 0.0184
              d = 1.5550
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
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
              z = 4
              a = 14.665
              b = -0.5140
              c = 2.7300
              d = 0.2328
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(4) ! O V
              a = -2.8425
              b = 0.2283
              c = 40.4072
              d = -3.4956
              f = 1.7558
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 5
              a = 4.626
              b = -0.9521
              c = 0.4729
              d = -0.4477
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select

     case(10)
        select case(thisIon%n) ! from nussbaumer and storey 1987
           case(10) ! Ne I
              z = 1
              a = 7.317
              b = -0.5358
              c = 2.4130
              d = 0.4176
              rate = ppb1991(z, a, b, c, d, dble(temperature))
           case(9) ! Ne II
              a = 0.0129
              b =-0.1779
              c = 0.9353
              d =-0.0682
              f = 0.4516
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 2
              a = 11.80
              b = -0.5404
              c = 3.0300
              d = 0.2050
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(8) ! Ne III
              a = 3.6781
              b = 14.1481
              c = 17.1175
              d = -0.5017
              f = 0.2313
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 3
              a = 5.841
              b = -0.5921
              c = 0.4498
              d = 0.6395
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(7) ! Ne IV
              a =-0.0254
              b = 5.5365
              c = 17.0727
              d = -0.7225
              f = 0.1702
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 4
              a = 15.550
              b = -0.4825
              c = 3.2740
              d = 0.3030
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(6) ! Ne V
              a = -0.0141
              b = 33.8479
              c = 43.1608
              d =-1.6072
              f = 0.1942
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 5
              a = 7.538
              b = -0.5540
              c = 1.2960
              d = 0.3472
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(5) ! Ne VI
              a = 19.9280
              b = 235.0536
              c = 152.5096
              d = 9.1413
              f = 0.1282
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 6
              a = 5.239
              b = -0.5966
              c = 0.7135
              d = 0.4537
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select

     case(16)
        select case(thisIon%n) ! from nussbaumer and storey 1987
           case(16) ! S I 
              rate = 3.e-13
              rate = rate + svs1982(dble(temperature), 4.10D-13, 6.30D-1)
           case(15) ! S II
              rate = 3.e-30 ! page 1344, para 2, kenny's paper
              rate = rate + svs1982(dble(temperature), 1.80D-12, 6.86D-1)
           case(14) ! S  III
              rate = 1.5e-11
              rate = rate + svs1982(dble(temperature), 2.70D-12, 7.45D-1)
           case(13) ! S IV
              rate = 2.5e-11
              rate = rate + svs1982(dble(temperature), 5.70D-12, 7.55D-1)
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
  rate = (1.d-12) * ((a/t) + b + (c*t) + d*(t**2))*(t**(-1.5d0))*exp(-f/t)

end function nussbaumerStorey1983

function ppb1991(z, a, b, c, d, temperature) result(rate)

! radiative recombination rates
! based on Pequignot, Petitjean & Boisson (1991) A&A 251 680

  real(double) :: z, a, b, c, d, temperature, t, rate

  t = 1.d-4 * temperature / z**2

  rate = 1.d-13 * z * (a*t**b)/(1.d0+c*t**d)

end function ppb1991

function svs1982(t, alpharad, xrad) result (rate)

! radiative recombination rates based on
! Shull and Van Steenberg, 1992, ApJS, 48, 95

  real(double) :: t, alpharad, xrad, rate

  rate = alpharad * (t /1.d4)**(-xrad)
end function svs1982

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
  real(double) :: hHeating, heHeating, totalHeating, heating, nh, nhii, nheii, ne
  real(double) :: cooling
  real :: netot
  open(20,file="lexington.dat",form="formatted",status="unknown")
  open(21,file="orates.dat",form="formatted",status="unknown")
  open(22,file="ne.dat",form="formatted",status="unknown")

  do i = 1, 50
     r = (1.+3.d0*dble(i-1)/49.d0)*pctocm/1.e10

     t=0;hi=0; hei=0;oii=0;oiii=0;cii=0;ciii=0;civ=0;nii=0;niii=0;niv=0;nei=0;neii=0;neiii=0;neiv=0;ne=0.

     oirate = 0; oiirate = 0; oiiirate = 0; oivrate = 0
     heating = 0.d0; cooling = 0.d0
     do j = 1, 100
        call random_number(theta)
        theta = theta * Pi
        
        octVec = OCTALVECTOR(r*sin(theta),0.d0,r*cos(theta))
        
        call amrgridvalues(grid%octreeRoot, octVec,  foundOctal=thisOctal, foundsubcell=subcell)

        nHii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,2) * grid%ion(2)%abundance
        nHeii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,4) * grid%ion(4)%abundance
        nh = thisOctal%nh(subcell)
        ne = thisOctal%ne(subcell)

        cooling = cooling + HHeCooling(grid, thisOctal, subcell, nh, nHii, nHeii, ne, thisOctal%temperature(subcell))
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
        netot = netot + thisOctal%ne(subcell)
        call getHeating(grid, thisOctal, subcell, hHeating, heHeating, totalHeating, epsOverDT)
        heating = heating + totalHeating
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
     netot = netot / 100

     oirate = oirate / 100.
     oiirate = oiirate / 100.
     oiiirate = oiiirate / 100.
     oivrate = oivrate / 100.
     heating = heating / 100.
     cooling = cooling / 100.

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
     ne = log10(max(ne,1.d-10))


     write(21,'(f5.3,1p,6e12.3,0p)') r*1.e10/pctocm,heating,cooling,oirate,oiirate,oiiirate,oivrate

     write(20,'(f5.3,f9.1,  14f8.3)') &
          r*1.e10/pctocm,t,hi,hei,oii,oiii,cii,ciii,civ,nii,niii,niv,nei,neii,neiii,neiv
     write(22,*) r*1.e10/pctocm,netot
  enddo
  close(20)
  close(21)
  close(22)
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

  t = max(min(real(thisIon%transition(iTransition)%t(thisIon%transition(iTransition)%ngamma)),temperature), &
       real(thisIon%transition(iTransition)%t(1)))
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
!     do i = 1, thisIon%transition(1)%ngamma
!        write(*,*) i,thisIon%transition(1)%t(i),thisIon%transition(1)%gamma(i)
!     enddo
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
          call solvePops(grid%ion(iIon), pops, thisOctal%ne(subcell), thisOctal%temperature(subcell), &
               thisOctal%ionFrac(subcell,iion),thisOctal%nh(subcell))
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
  do j = 3, nIons
     if (ionArray(j)%nTransitions > 0) then
        call solvePops(ionArray(j), pops, ne, temperature, thisOctal%ionFrac(subcell,j), thisOctal%nh(subcell))
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
  

subroutine solvePops(thisIon, pops, ne, temperature, ionFrac, nh, debug)
  type(IONTYPE) :: thisIon
  real(double) :: ne
  real :: pops(*)
  real :: temperature
  real(double), allocatable :: matrixA(:,:), MatrixB(:), tempMatrix(:,:), qeff(:,:),  rates(:)
  integer :: n, iTrans, i, j
  real :: excitation, deexcitation, rateij, rateji, arateji
  real(double) :: nh, ionFrac
  logical :: ok
  logical, optional :: debug

  n = thisIon%nLevels
  allocate(matrixA(1:n, 1:n), matrixB(1:n), tempMatrix(1:n,1:n), qeff(1:n,1:n), rates(1:n))

  matrixA = 1.d-30
  matrixB = 0.d0


  call getRecombs(rates, thision, dble(temperature), ne, ionFrac, nh)

  matrixA(1,:) = 1.d0
  matrixB(1) = 1.d0

  matrixB(2:n) = 0.d0 !rates(2:n) * ne * (nh * ionFrac * thisIon%abundance)
  do iTrans = 1, thisIon%nTransitions
     i = thision%transition(itrans)%i
     j = thision%Transition(itrans)%j
     call getCollisionalRates(thisIon, iTrans, temperature, excitation, deexcitation)
     qeff(i,j) = excitation
     qeff(j,i) = deexcitation
  enddo

  do i = 2, n
     do j = 1, n

        do iTrans = 1, thisIon%nTransitions
           if (((i == thision%transition(itrans)%i).and.(j == thision%Transition(itrans)%j)).or. &
               ((i == thision%transition(itrans)%j).and.(j == thision%Transition(itrans)%i))) then
              arateji =  thisIon%transition(iTrans)%a
              
              matrixA(i,j) = matrixA(i,j) + ne * qeff(j, i)
              matrixA(i,i) = matrixA(i,i) - ne * qeff(i, j)
              if (j > i) then
                 matrixA(i,j) = matrixA(i,j) + arateji
              else
                 matrixA(i,i) = matrixA(i,i) - arateji
              endif

           endif
        enddo

     enddo
  enddo

  tempMatrix = matrixA

  if (PRESENT(debug)) then
     if (debug) then
        do i = 1, n
           write(*,'(1p,9e12.3)') tempmatrix(i,1:n)
        enddo
     endif
  endif
  
  call luSlv(matrixA, matrixB, n)

  matrixB(1:n) = matrixB(1:n) / SUM(matrixB(1:n))
  
  ok = .true.

  if (.not.ok) then
     write(*,*) "Population solver failed for: ",thisIon%species
     write(*,*) matrixB(1:n)
     write(*,*) "nlevels",thisIon%nLevels,"ntrans",thisIon%nTransitions 
     write(*,*) "temp",temperature,"ne",ne
     
     do i = 1, n
        write(*,'(1p,9e12.3)') tempmatrix(i,1:n)
     enddo
     matrixB = 0.d0
     matrixB(1) = 1.d0
     write(*,*) "Setting pops to ground state"
  endif

  do i = 1, n
     pops(i) = max(1.d-30,matrixB(i))
  enddo

  deallocate(matrixA, matrixB, tempMatrix, qeff, rates)

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
  real(double) :: nu0_h, nu0_he, nufinal_h, nufinal_he
  integer :: nFreq, nTemp
  integer :: i, j
  real :: e
  real(double) :: t1, t2
  real :: hxsec, hexsec
  real(double) :: dfreq, jnu

  nFreq = 10000
  nTemp = 100
  hTable%nFreq = nFreq
  heTable%nFreq = nFreq
  hTable%nTemp = nTemp
  heTable%nTemp = nTemp

  allocate(hTable%freq(1:nFreq), heTable%freq(1:nFreq))
  allocate(hTable%emissivity(1:nTemp), heTable%emissivity(1:ntemp))
  allocate(hTable%temp(1:ntemp), heTable%temp(1:ntemp))
  allocate(hTable%Clyc(1:nTemp,1:nFreq))
  allocate(heTable%Clyc(1:nTemp,1:nFreq))

  nu0_h = 13.6d0/ergtoev/hcgs
  nu0_he = 24.59d0/ergtoev/hcgs

  nufinal_h = 2.d0*nu0_h
  nufinal_he = 2.d0*nu0_he

  do i = 1, nFreq
     hTable%freq(i) = log10(nu0_h) + (log10(nuFinal_h)-log10(nu0_h))*dble(i-1)/dble(nFreq-1)
     heTable%freq(i) = log10(nu0_he) + (log10(nuFinal_he)-log10(nu0_he))*dble(i-1)/dble(nFreq-1)
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
        jnu = ((hcgs*hTable%freq(j)**3)/(cSpeed**2)) * ((hcgs**2) /(twoPi*mElectron*Kerg*hTable%temp(i)))**(1.5d0) * &
             dble(hxsec/1.d10) *  exp(-hcgs*(hTable%freq(j)-nu0_h)/(kerg*hTable%temp(i)))
        hTable%Clyc(i,j) = hTable%Clyc(i,j-1) + jnu * dfreq


        e = heTable%freq(j) * hcgs * ergtoev
        call phfit2(2, 2, 1 , e , hexsec)

        dFreq = heTable%freq(j)-heTable%freq(j-1)
        jnu = 2.d0*((hcgs*heTable%freq(j)**3)/(cSpeed**2)) * ((hCgs**2) / (twoPi*mElectron*kerg*heTable%temp(i)))**(1.5d0) * &
             dble(hexsec/1.d10) * exp(-hcgs*(heTable%freq(j)-nu0_he)/(kerg*heTable%temp(i)))
        heTable%Clyc(i,j) = heTable%Clyc(i,j-1) + jnu * dfreq
     enddo
     hTable%emissivity(i) = hTable%Clyc(i,hTable%nFreq)
     heTable%emissivity(i) = heTable%Clyc(i,heTable%nFreq)
     hTable%Clyc(i,1:hTable%nFreq) = hTable%Clyc(i,1:hTable%nFreq) / hTable%Clyc(i,hTable%nFreq)
     heTable%Clyc(i,1:heTable%nFreq) = heTable%Clyc(i,1:heTable%nFreq) / heTable%Clyc(i,heTable%nFreq)
  end do

  do j = 1, nFreq
     write(99  ,*) htable%freq(j),htable%clyc(50,j)
  enddo


end subroutine createSahaMilneTables

subroutine getSahaMilneFreq(table,temperature, thisFreq)
  type(SAHAMILNETABLE) :: table
  real(double) :: temperature, thisfreq, r, t, fac
  integer :: i, j

  t = max(5000.d0, min(20000.d0, temperature))
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
  thisFreq = 0.
  do while((thisFreq*hcgs*ergtoev) < 13.6)
     call random_number(r)
     call locate(prob, 21, r, i)
     fac = y(i) + ((r - prob(i))/(prob(i+1)-prob(i)))*(y(i+1)-y(i))
     thisFreq = (1.-fac)*freq
  enddo
end subroutine twoPhotonContinuum

subroutine getRecombs(rates, thision, temperature, ne, ionFrac, nh)  
  real(double) :: rates(:)
  type(IONTYPE) :: thisIon
  real(double) :: temperature, ne, ionfrac, nh

  select case(thisIon%species)
     case("O II")
        rates(1) = 0.d0
        rates(2) = rateFit(7.218d0, -0.575d0, temperature)
        rates(3) = rateFit(4.812d0, -0.575d0, temperature)
        rates(4) = rateFit(3.581d0, -0.495d0, temperature)
        rates(5) = rateFit(1.790d0, -0.495d0, temperature)
     case DEFAULT
        rates = 0.d0
   end select
 end subroutine getRecombs


subroutine getHeating(grid, thisOctal, subcell, hHeating, heHeating, totalHeating, epsOverDeltaT)
  type(GRIDTYPE) :: grid
  type(OCTAL) :: thisOctal
  integer :: subcell
  real(double) :: hHeating, heHeating, totalHeating, v, r1, r2, epsOverDeltaT
  type(OCTALVECTOR) :: rVec

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
end subroutine getHeating

function rateFit(a,b,t) result (rate)
  real(double) :: a, b, t, rate

  rate = 1.d-13 * a * (t/1.d4)**b
end function rateFit

subroutine createRecombTable(table, tablefilename)
  type(RECOMBTABLE) :: table
  character(len=*) :: tablefilename
  character(len=80) :: filename, datadirectory
  integer :: ia, ib, ne, ncut, i
  real :: e(1000)

  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//tablefilename

  open(20,file=filename,status="old",form="formatted")
  read(20,*) table%nTemp, table%nrho

  allocate(table%rho(1:table%nRho))
  allocate(table%temp(1:table%ntemp))
  allocate(table%emissivity(1:table%ntemp,1:table%nrho))

  do ia=1,table%ntemp
     do ib=1,table%nrho
        read(20,'(1x,e10.3,5x,e10.3,13x,i2)') table%rho(ib),table%temp(ia),ncut
        ne=ncut*(ncut-1)/2
        read(20,'(8e10.3)') e(1:ne)
        table%emissivity(ia,ib) = SUM(e(1:ne-1))
     enddo
  enddo
  close(20)
end subroutine createRecombTable


subroutine calcEmissivities(grid, thisOctal, subcell, hTable, hiGammatable, hRecombTable, &
     hRecombemissivity, lymanContemissivity, HIcontEmissivity, forbiddenEmissivity)
  type(GRIDTYPE) :: grid
  type(OCTAL) :: thisOctal
  integer :: subcell
  type(SAHAMILNETABLE) :: hTable
  type(RECOMBTABLE) :: hRecombTable
  type(GAMMATABLE) :: hiGammaTable
  integer :: i, j
  real(double) :: hRecombEmissivity, lymanContemissivity, fac, t, u
  real(double) :: hiContemissivity, forbiddenEmissivity

  call locate(hRecombTable%temp, hRecombTable%nTemp, dble(thisOctal%temperature(subcell)), i)
  call locate(hRecombTable%rho, hRecombTable%nRho, thisOctal%ne(subcell), j)
  t = (thisOctal%temperature(subcell) - hRecombTable%temp(i)) / &
       (hRecombTable%temp(i+1) - hRecombTable%temp(i))
  u = (thisOctal%ne(subcell) - hRecombTable%rho(i)) / &
       (hRecombTable%rho(i+1) - hRecombTable%rho(i))
  hRecombEmissivity = (1.d0-t) * (1.d0-u) * log10(hRecombTable%emissivity(i,j)) + &
       (     t) * (1.d0-u) * log10(hRecombTable%emissivity(i+1,j)) + &
       (1.d0-t) * (     u) * log10(hRecombTable%emissivity(i,j+1)) + &
       (     t) * (     u) * log10(hRecombTable%emissivity(i+1,j+1)) 
  hRecombEmissivity = 10.d0**hRecombemissivity

  hRecombEmissivity = hRecombemissivity * thisOctal%ne(subcell) * &
       (thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,2) * grid%ion(1)%abundance)

  hIcontEmissivity = getGamma(HiGammaTable, dble(thisOctal%temperature(subcell))) * &
    thisOctal%ne(subcell) *(thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,2) * &
    grid%ion(1)%abundance)


  call locate(htable%temp, htable%nTemp, dble(thisOctal%temperature(subcell)), i)
  fac = (thisOctal%temperature(subcell) - hTable%temp(i)) &
       / (hTable%temp(i+1) - hTable%temp(i))

  lymanContEmissivity = log10(hTable%emissivity(i)) + fac*(log10(hTable%emissivity(i+1)/htable%emissivity(i)))
  lymanContEmissivity = 10.d0**lymanContEmissivity
  lymanContEmissivity = lymanContEmissivity &
       * thisOctal%ne(subcell) *(thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,2) * grid%ion(1)%abundance)

  call metalcoolingRate(grid%ion, grid%nIon, thisOctal, subcell, thisOctal%nh(subcell),thisOctal%ne(subcell), &
  thisOctal%temperature(subcell), forbiddenEmissivity)
  
  forbiddenEmissivity = forbiddenEmissivity/fourPi

end subroutine calcEmissivities

function getGamma(table, temp) result (gamma)
  type(GAMMATABLE) :: table
  real(double) :: temp, gamma, fac
  integer :: i
  call locate(table%temp, table%nTemp, temp, i)
  fac = (temp - table%temp(i))/(table%temp(i+1)-table%temp(i))
  gamma = log10(table%gamma(i)) + fac * log10(table%gamma(i+1)/table%gamma(i))
  gamma = 10.d0**gamma
end function getGamma
  
subroutine createHIgammaTable(table)

  type(GAMMATABLE) :: table
  real(double) :: coeff1(11,6), coeff2(11,6)
  real(double) :: temp(6) = (/ 500.d0, 1000.d0, 2500.d0, 5000.d0, 10000.d0, 20000.d0 /)
  real(double) :: sum, val, fac, freq(1000), freq1, freq2, dfreq
  integer :: i, j, k
  coeff1(1,1:6)  = (/ 2.92d-140, 3.03d-89, 5.83d-59, 4.07d-49, 2.11d-44, 3.29d-42 /)
  coeff2(1,1:6)  = (/ 2.01d-37, 7.27d-38, 1.87d-38, 6.66d-39, 2.48d-39, 1.06d-39 /)
  coeff1(2,1:6)  = (/ 6.10d-57, 7.37d-48, 1.00d-42, 3.23d-41, 1.37d-40, 2.32d-40 /)
  coeff2(2,1:6)  = (/ 6.35d-38, 2.27d-38, 5.92d-39, 2.38d-39, 1.15d-39, 6.78d-40 /)
  coeff1(3,1:6)  = (/ 6.25d-45, 4.89d-42, 1.48d-40, 3.38d-40, 4.26d-40, 4.23d-40 /)
  coeff2(3,1:6)  = (/ 2.76d-38, 9.97d-37, 3.02d-39, 1.51d-39, 9.04d-40, 6.31d-40 /)
  coeff1(4,1:6)  = (/ 1.24d-41, 1.69d-40, 5.32d-40, 6.26d-40, 5.93d-40, 5.21d-40 /)
  coeff2(4,1:6)  = (/ 1.45d-38, 5.69d-49, 21.4d-39, 1.25d-39, 8.51d-40, 6.41d-40 /)
  coeff1(5,1:6)  = (/ 1.96d-40, 5.96d-40, 8.46d-40, 7.98d-40, 6.90d-40, 5.84d-40 /)
  coeff2(5,1:6)  = (/ 9.04d-39, 4.00d-39, 1.80d-39, 1.17d-39, 1.17d-39, 8.50d-40 /)
  coeff1(6,1:6)  = (/ 6.46d-40, 1.03d-39, 1.05d-39, 9.04d-40, 7.56d-40, 6.31d-40 /)
  coeff2(6,1:6)  = (/ 6.48d-39, 3.23d-39, 1.65d-39, 1.15d-39, 8.66d-40, 6.90d-40 /)
  coeff1(7,1:6)  = (/ 1.16d-39, 1.35d-39, 1.18d-39, 9.79d-40, 8.06d-40, 6.69d-40 /)
  coeff2(7,1:6)  = (/ 5.18d-39, 2.84d-39, 1.59d-39, 1.15d-39, 8.87d-40, 7.16d-40 /)
  coeff1(8,1:6)  = (/ 1.60d-39, 1.58d-39, 1.27d-39, 1.04d-39, 8.47d-40, 7.02d-40 /)
  coeff2(8,1:6)  = (/ 4.44d-39, 2.63d-39, 1.56d-39, 1.16d-39, 9.11d-40, 7.41d-40 /)
  coeff1(9,1:6)  = (/ 1.93d-39, 1.74d-39, 1.34d-39, 1.08d-39, 8.82d-40, 7.31d-40 /)
  coeff2(9,1:6)  = (/ 4.01d-39, 2.51d-39, 1.56d-39, 1.18d-39, 9.34d-39, 7.64d-40 /)
  coeff1(10,1:6)  = (/ 2.17d-39, 1.86d-39, 1.39d-39, 1.12d-39, 9.14d-40, 7.57d-40 /)
  coeff2(10,1:6)  = (/ 3.73d-39, 2.44d-39, 1.56d-39, 1.20d-39, 9.58d-40, 7.87d-40 /)
  coeff1(11,1:6)  = (/ 2.36d-39, 1.95d-39, 1.44d-39, 1.16d-39, 9.42d-40, 7.81d-40 /)
  coeff2(11,1:6)  = (/ 3.56d-39, 2.40d-39, 1.57d-39, 1.22d-39, 9.80d-40, 8.08d-40 /)

  table%nTemp = 6
  allocate(table%gamma(6), table%temp(6))
  table%temp = temp

  do i = 1, 6
     sum = 0.d0
     do j = 1, 11
        freq1 = nuHydrogen / dble((j+1)**2)
        freq2 = nuHydrogen /dble(j**2)
        do k = 1, 1000
           freq(k) = log(freq1)+ (log(freq2)-log(freq1))*dble(k-1)/999.d0
!           freq(k) = freq1 + (freq2-freq1)*dble(k-1)/999.d0
        enddo
        freq(1:1000) = exp(freq(1:1000))
        do k = 2, 1000
           dfreq = freq(k)-freq(k-1)
           fac = (log(freq(k))-log(freq1))/(log(freq2)-log(freq1))
!           fac = dble(k-1)/999.d0
           val = log(coeff2(j,i)) + (log(coeff1(j,i))-log(coeff2(j,i)))*fac
           sum = sum + exp(val)*dFreq
        enddo
     enddo
     write(*,*) i,sum
     table%gamma(i) = sum
  enddo
end subroutine createHIgammaTable

     
end module

