module fixing_mod
  use octal_mod, ONLY: octal, subcellCentre
  use vector_mod
  use constants_mod
  use inputs_mod, ONLY: griddistancescale, cylindricalHydro, nHydroThreadsInput, smallestCellSize
  use grid_mod
  use source_mod, only: globalsourcearray

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fixOctalValues
  
contains

  subroutine fixOctalValues(grid)
    type(GRIDTYPE) :: grid
    logical, save :: flag=.true.

!    write(*,*) "fixing"    
    if (flag) call fixOctalValues_pri(grid%octreeroot, grid%currentTime)
    
  end subroutine fixOctalValues
    
  
  recursive subroutine fixOctalValues_pri(thisOctal, currentTime) 
    type(octal), pointer   :: thisOctal, child
    real(double) :: currentTime
!    real(double) :: r,theta,phi
    real(double) :: R
    type(vector) :: cen
!    real(double) :: fixedRho, fixedTemp, ethermal, uSquared
    integer :: subcell, i
 !   integer :: ri,thetai
 !   real(double), save :: rho(1000,1500), ur(1000,1500), uz(1000,1500), uphi(1000,1500)
 !   logical, save :: gotdata=.false.
 !   real(double) :: dr=2.960013665777d12 ! dr in cm and dtheta from
 !   real(double) :: dtheta=1.5723686954903871d-3 ! output array

    
 !   if (.not. gotdata) then
 !      open (unit=10, file="rho.bin", access='stream')
 !      read (10), rho
 !      close(10)
 !      open (unit=11, file="vr.bin", access='stream')
 !      read (11), ur
 !      close(11)
 !      open (unit=12, file="vz.bin", access='stream')
 !      read (12), uz
 !      close(12)
 !      open (unit=13, file="vphi.bin", access='stream')
 !      read (13), uphi
 !      close(13)
 !      gotdata=.true.
 !   endif 

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
           do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fixOctalValues_pri(child, currentTime)
                exit
             endif
           end do
       elseif (currentTime >=0) then !negative times imply we're not evolving the hydro so dont want to be fixing things

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
!          if (myrankGlobal.ne.0) cycle
          
          cen=subcellCentre(thisOctal, subcell)
!example of fixing function
             
!             if  (modulus(cen)*griddistancescale<(7.0e10*5)) then !if the distance from centre is less than 5R_sol
!                eThermal = 1.5*(1.d0/(mHydrogen))*kerg*fixedTemp
!                thisOctal%energy(subcell)=eThermal
!                thisOctal%rho(subcell)=fixedRho          
!                thisOctal%rhoe(subcell)=eThermal*fixedRho!fix the density,energy and tempature to photosphere values
!                thisOctal%temperature(subcell)=fixedTemp
!                thisOctal%rhou(subcell) = 0.0      !and the velocity to 0
!                thisOctal%rhov(subcell) = 0.0
!                thisOctal%rhow(subcell) = 0.0
!             endif
             
             !another example, temperature ceiling
!             fixedtemp = 1.0e4
!             if (thisOctal%temperature(subcell) > fixedtemp) then
!                if (thisOctal%threed) then
!                   uSquared = (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2)&
!                        /thisOctal%rho(subcell)**2
!                else if (thisOctal%twoD) then
!                   uSquared = (thisOctal%rhou(subcell)**2 +  thisOctal%rhow(subcell)**2)&
!                        /thisOctal%rho(subcell)**2
!                else
!                   uSquared=(thisOctal%rhou(subcell)/thisOctal%rho(subcell))**2
!                endif
                
!                ethermal = (1.d0/(mHydrogen))*kerg*fixedTemp
!                thisOctal%rhoe(subcell)=(ethermal+uSquared/2.0)*thisOctal%rho(subcell) !remove all energy above 10^4K e.g. for optically thin radiation 
!                thisOctal%temperature(subcell)=fixedTemp
!             endif
             
!          else !explode after 40 timesteps
!             if  (modulus(cen)*griddistancescale<(7.0e10*4.5)) then
!                ethermal = thisOctal%gamma(subcell)*(1.d0/(2.0*mHydrogen))*kerg*1.0d9
!                thisOctal%rho(subcell) =1.0d-12
!                thisOctal%rhoe(subcell)=ethermal*thisOctal%Rho(subcell)
!                thisOctal%temperature(subcell)=1.0d9
!             endif
!          endif

!applying griddata to a model

!          r     =sqrt(cen%x**2+cen%y**2+cen%z**2)*GridDistanceScale
!          theta =atan2(sqrt(cen%x**2+cen%y**2), abs(cen%z))
!          phi   =atan2(cen%y,cen%x)
!          ri    =r/dr
!          thetai=theta/dtheta
!          ri    =max(min(ri,1500),1)
!          thetai=max(min(thetai,1000),1)

!          thisOctal%rho(subcell) =rho(thetai,ri)
!          if (CylindricalHydro) then
!              thisOctal%rhou(subcell)=rho(thetai,ri)*ur(thetai,ri)
!              thisOctal%rhov(subcell)=rho(thetai,ri)*uphi(thetai,ri)*sqrt(cen%x**2+cen%y**2)
!          else
!              thisOctal%rhou(subcell)=rho(thetai,ri)*ur(thetai,ri)*cos(phi)  - uphi(thetai,ri)*sin(phi)
!              thisOctal%rhov(subcell)=rho(thetai,ri)*uphi(thetai,ri)*cos(phi)+ ur(thetai,ri)*sin(phi)
!          endif
!          thisOctal%rhow(subcell)=rho(thetai,ri)*uz(thetai,ri)

!                write (*,*) cen%x, cen%y, cen%z
!                write (*,*) thisOctal%subcellsize

          if (cylindricalHydro) then
             if (cen%x <=(thisOctal%subcellSize)) then
                thisOctal%rhou(subcell)=0!cells on the axis cant have radial momentum
                thisOctal%rhov(subcell)=0!cells on the axis cant have momentum in the phi direction 
             endif
             R=(cen%x*gridDistanceScale)
             if (abs(cen%z) <= smallestCellSize/2 ) then
                thisOctal%rhow(subcell)=0.0d0 !at finest level force disc midplane to not have z velocity
                thisOctal%rhou(subcell)=min(0.0d0,thisOctal%rhou(subcell))!        to have inwards or 0 radial velocity
                if (cen%x > 0) then
                   thisOctal%rhov(subcell)=min(thisOctal%rhov(subcell), &
                                               thisOctal%rho(subcell)*R*sqrt(globalSourceArray(1)%mass*bigG/R))!and to have (sub) keplerian rotational velocity
                endif
             endif
          endif
          
       endif
    enddo
  end subroutine fixOctalValues_pri

 end module fixing_mod
