#ifdef PDR
module pdr_mod
!Started by THaw and Bisbas 13/08/2013

use constants_mod
use messages_mod
use parallel_mod
use gridio_mod
use source_mod
use timing
use grid_mod
use vtk_mod
use amr_mod
use mpi_amr_mod
use mpi_global_mod
use unix_mod, only: unixGetenv
use definitions
!use healpix_module
use vector_mod

implicit none



! REAL(kind=dp), save :: OMEGA=0.42D0,GRAIN_RADIUS=1.0D-7,METALLICITY=1.0D0


! real(kind=dp),allocatable :: all_heating(:,:)
! real(kind=dp),allocatable :: allheating(:) 

contains

subroutine PDR_MAIN(grid)
  use unix_mod, only: unixGetenv
  use setuppdr_mod
  use nrayshealpix
  implicit none
  type(gridtype) :: grid


  call donrayshealpix()

!set up the data/allocate what needs allocated
  call writeInfo("Setting up for PDR calculation.", TRIVIAL)
  call setupPDR(grid)

!do the ray casting
  call writeInfo("Casting rays over grid.", TRIVIAL)
  call castAllRaysOverGrid(grid%octreeRoot, grid)
  call writeInfo("Done.", TRIVIAL)
  call writeInfo("Calculating dust temperature.", TRIVIAL)
  call calculate_Dust_Temperatures(grid%octreeRoot)
  call writeInfo("Done.", TRIVIAL)

  call writeVTKfile(grid, "columnDensity.vtk", valueTypeString=(/"rho       ",&
       "columnRho ", "UV        ", "uvvec     ", "dust_T    "/))

end subroutine PDR_MAIN

integer function nray_func()
  use inputs_mod, only : hlevel
  integer :: nside
  nside = 2**hlevel
  nray_func = 12*nside**2
end function nray_func
!This is the main routine that loops over the grid and does the ray casting
recursive subroutine castAllRaysOverGrid(thisOctal, grid)
  use inputs_mod, only : hlevel, maxdepthamr
  use healpix_module, only : vectors, nrays
  implicit none
  
  type(gridtype) :: grid
  integer :: subcell, ssubcell
  type(octal), pointer :: thisoctal, soctal
  type(octal), pointer :: child
  integer :: j
  type(vector) :: testPosition, rVec, startPosition, uhat, thisUVvector
  real(double) :: tVal
  integer :: nside, i, ncontributed, k



  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do j = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(j) == subcell) then
              child => thisoctal%child(j)
              call castAllRaysOverGrid(child, grid)
              exit
           end if
        end do
     else
!        if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle
!        print *, "thisOctal%ndepth", thisOctal%ndepth
 !       print *, "casting rays at", subcell!, subcellCentre(thisOctal, subcell)

!        call doSingleCellRays(thisOctal, subcell, grid)
 
        sOctal=> thisOctal
        
        startPosition = subcellCentre(thisOctal, subcell)
        testPosition = startPosition
        call findSubcellLocal(testPosition, sOctal, ssubcell)
        sOctal%columnRho(ssubcell) = 0.d0
!        thisOctal%columnRho(subcell) = 1.d-30 

        thisOctal%uv(subcell) = 0.d0
        thisOctal%av(subcell, :) = 0.d0           
        thisOctal%radsurface(subcell, :) = 0.d0           
        thisOctal%thisColRho(subcell, :, :) = 0.d0
        if(.not. thisOctal%ionfrac(subcell, 2) > 0.9d0) then
           ncontributed = 0
!        do i = 0, nrays-1
        do i = 1, nrays

           rVec = VECTOR(vectors(1, i-1), vectors(2, i-1), vectors(3, i-1))

           uhat = rVec
           call normalize(uhat)

           testPosition = startPosition
           thisOctal%columnrho(subcell)  = 0.d0
           do while (inOctal(grid%octreeRoot, testPosition))
              
              call findSubcellLocal(testPosition, sOctal, ssubcell)

              tval = 0.d0

             call distanceToCellBoundary(grid, testPosition, uHat, tVal, soctal, ssubcell)
 
              if(sOctal%ionFrac(ssubcell,2) > 0.9d0) then
                 thisUVvector = sOctal%uvVector(ssubcell)*1.d10/draine
!                 call normalize(thisUVvector)
!                 print *, "thisUVvector ", thisUVvector
                 thisOctal%radsurface(subcell, i) = - dotprod(uHat,thisUVvector)     
!                 thisOctal%radsurface(subcell, i) =  dotprod(uHat,thisUVvector)     
!                 print *, "thisOctal%radsurface(subcell, i)", thisOctal%radsurface(subcell, i)
!                 print *, "radsurface", thisOctal%radsurface(subcell, i)
                 if(thisOctal%radsurface(subcell, i) < 0.d0 ) thisOctal%radsurface(subcell, i) = 0.d0

                 exit
              end if

              
              thisOctal%columnRho(subcell) = thisOctal%columnRho(subcell) + (sOctal%rho(ssubcell)*tval)
              !loop over all species
              do k = 1, 33
                thisOctal%thisColRho(subcell, k, i)  = thisOctal%thisColRho(subcell, k, i) + &
                     thisOctal%abundance(subcell, k)*tval*1.d10*soctal%rho(subcell)/mhydrogen
                 
              enddo
!              print *, "RAY ", I, "HAS"
!              print *, thisOctal%thiscolRho(subcell, :, i)
              testPosition = testPosition + ((tVal+1.d-10*grid%halfsmallestsubcell)*uhat)
              
           end do


!          thisOctal%AV(subcell, i) = thisOctal%columnRho(subcell)*1.d10*AV_fac/mhydrogen
          thisOctal%AV(subcell, i) = thisOctal%columnRho(subcell)*1.d10*AV_fac/mhydrogen


          if(thisOctal%radsurface(subcell, i) > 0.d0) then
             thisOctal%UV(subcell) = thisOctal%UV(subcell) + &
                  thisOctal%radSurface(subcell, i) * exp(-(thisOctal%AV(subcell, i)*UV_fac))
             ncontributed = ncontributed + 1
!             print *, "UV ", thisOctal%UV(subcell)
!             print *, "exp(-(thisOctal%AV(subcell, i)*UV_fac", exp(-(thisOctal%AV(subcell, i)*UV_fac))

!             print *,thisOctal%radSurface(subcell, i) * exp(-(thisOctal%AV(subcell, i)*UV_fac))
          end if
          
       end do
!       print *, "DONE ONE "
!       print *, 'ncontributed ', ncontributed
       if(ncontributed /= 0) then
!          thisOctal%UV(subcell) = thisOctal%UV(subcell) / dble(nrays)
          thisOctal%UV(subcell) = thisOctal%UV(subcell) / dble(ncontributed)
       else
          thisOctal%UV(subcell) = 0.d0
       end if
!       thisOctal%UV(subcell) = thisOctal%UV(subcell) / dble(nrays)
       thisOctal%columnRho(subcell) = thisOctal%columnRho(subcell)/dble(nrays)
 !      print *, "AV ", thisOctal%AV(subcell,:)
  !     print *, "UV ", thisOctal%UV(subcell)
       
!       stop

       if(thisOctal%UV(subcell) < 1.d-30) thisOctal%UV(subcell) = 0.d0
    else
        thisOctal%columnRho(subcell) = 0.d0
        thisOctal%UV(subcell) = 0.d0
        thisOctal%AV(subcell,:) = 0.d0
        thisOctal%radsurface(subcell,:) = 0.d0
     end if
        
  end if
end do
  

end subroutine castAllRaysOverGrid



!A quick test of ray tracing
subroutine fireTestRays(grid)
  use inputs_mod, only : amrgridcentrex, amrgridcentrey, amrgridcentrez
  use inputs_mod, only : hlevel, maxdepthamr
  use healpix_module, only : vectors, nrays

  implicit none

  type(gridtype) :: grid
  type(vector) :: testPosition, rVec, startPosition, uhat
  integer :: subcell
  type(octal), pointer :: thisOctal
  real(double) :: x, y, z, tVal, totalLength, totalRho
  integer :: nside, i, ier


  nside = 2**hlevel
  nrays = 12*nside**2

  thisOctal=>grid%octreeRoot
  
  !choose a test location
  x = amrgridcentrex + grid%octreeroot%subcellsize*0.34d0
  y = amrgridcentrey + grid%octreeroot%subcellsize*0.34d0
  z = amrgridcentrez + grid%octreeroot%subcellsize*0.34d0
  startPosition = VECTOR(x, y, z)

  call findSubcellLocal(startPosition, thisOctal, subcell)
  
  print *, "testing rays from cell at ", startPosition
  print *, "testing with ", nrays, "rays"
  print *, "cell size is ", grid%octreeroot%subcellsize/2.**maxdepthamr/1.d8
  print *, "grid size is" , grid%octreeroot%subcellsize/1.d8

  if (thisoctal%haschild(subcell)) then
     call torus_abort("Octal has child!")
  end if

  open (1, file="testRays.dat", status="unknown", iostat=ier)
  if(ier /= 0) then
     call torus_abort("Problem opening testRays.dat")
  end if

  do i = 0, nrays-1

     !transfer healpix vectors to the test position
!     rVec = startPosition + VECTOR(vectors(1, i), vectors(2, i), vectors(3, i))
     rVec = VECTOR(vectors(1, i), vectors(2, i), vectors(3, i))
     uhat = rVec
     call normalize(uhat)
     totalLength = 0.d0
     totalRho = 0.d0
     testPosition = startPosition
     print *, "doing ray ", i!, " in direction ", uHat 
     do while (inOctal(grid%octreeRoot, testPosition))
!        print *, "walking in direction ", uHat
        call findSubcellLocal(testPosition, thisOctal, subcell)
        call distanceToCellBoundary(grid, testPosition, uHat, tVal, thisoctal)
!        print *, " ini ", testposition
        testPosition = testPosition + ((tVal+1.d-10*grid%halfsmallestsubcell)*uhat)
 !       print *, "tval", tval
  !      print *, "uhat", uhat
   !     print *, "post ", testposition
    !    stop

        write(1,'(1p,i3,3e14.2)') i, testPosition%x*1.d10/pctocm, testposition%y*1.d10/pctocm, &
             testposition%z*1.d10/pctocm
        totalLength = totalLength + (tVal)
!       print *, totallength/1.d8
        totalRho = totalRho + thisOctal%rho(subcell)*tval

     end do
     print *, "ray ", i, "travelled a distance", totalLength*1.d10/pctocm
     print *, "ray ", i, "had total density ", totalRho
     print *, "ray ", i, "had total H number density ", totalRho/mHydrogen

!     print *, "ray ", i, "had start ", totalRho
  end do

  close(1)

end subroutine fireTestRays

!subroutine


!=======================================================================
!
!  Calculate the dust temperature for each particle using the treatment
!  of Hollenbach, Takahashi & Tielens (1991, ApJ, 377, 192, eqns 5 & 6)
!  for the heating due to the incident FUV photons and the treatment of
!  Meijerink & Spaans (2005, A&A, 436, 397, eqn B.6) for heating due to
!  the incident flux of X-ray photons.
!
!  Among other things, the dust temperature can influence:
!
!     1) Cooling budget by emitting FIR photons that
!        interact with the line radiative transfer;
!     2) Gas-grain collisional heating or cooling rate;
!     3) H2 formation by changing the sticking probability;
!     4) Evaporation and condensation of molecules on grains.
!
!  The formula derived by Hollenbach, Takahashi & Tielens (1991) has
!  been modified to include the attenuation of the IR radiation. The
!  incident FUV radiation is absorbed and re-emitted in the infrared
!  by dust at the surface of the cloud (up to Av ~ 1mag). In the HTT
!  derivation, this IR radiation then serves as a second heat source
!  for dust deeper into the cloud. However, in their treatment, this
!  second re-radiated component is not attenuated with distance into
!  the cloud so it is *undiluted* with depth, leading to higher dust
!  temperatures deep within the cloud which in turn heat the gas via
!  collisions to unrealistically high temperatures. Models with high
!  gas densities and high incident FUV fluxes (e.g. n_H = 10^5 cm-3,
!  X_0 = 10^8 Draine) can produce T_gas ~ 100 K at Av ~ 50 mag!
!
!  Attenuation of the FIR radiation has therefore been introduced by
!  using an approximation for the infrared-only dust temperature from
!  Rowan-Robinson (1980, eqn 30b):
!
!  T_dust = T_0*(r/r_0)^(-0.4)
!
!  where r_0 is the cloud depth at which T_dust = T_0, corresponding
!  to an A_V of ~ 1 mag, the assumed size of the outer region of the
!  cloud that processes the incident FUV radiation and then re-emits
!  it in the FIR (see the original HTT 1991 paper for details). This
!  should prevent the dust temperature from dropping off too rapidly
!  with distance and maintain a larger warm dust region (~50-100 K).
!
!-----------------------------------------------------------------------
recursive SUBROUTINE CALCULATE_DUST_TEMPERATURES(thisOctal)

   USE HEALPIX_TYPES
!   USE MAINCODE_MODULE

   IMPLICIT NONE

   INTEGER :: J, i
   type(octal), pointer :: thisOctal, child
   integer :: subcell, nrays
   REAL(double) :: NU_0,R_0,T_0,TAU_100
   REAL(double) :: T_CMB

!  Parameters used in the HHT equations (see their paper for details)
   NU_0=2.65D15
   TAU_100=1.0D-3
   R_0=1.0D0/AV_FAC
   T_CMB=2.73D0

   nrays = nray_func()

  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do i = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(i) == subcell) then
              child => thisoctal%child(i)
              call calculate_dust_temperatures(child)
              exit
           end if
        end do
     else
!        P=IDlist_pdr(pp)
        !     Calculate the contribution to the dust temperature from the local FUV flux and the CMB background
!        PDR(P)%DUST_T=8.9D-11*NU_0*(1.71D0*PDR(P)%UVfield)+T_CMB**5
        if(.not. thisOctal%ionfrac(subcell, 2) > 0.9d0) then
           thisOctal%DUST_T(subcell)=8.9D-11*NU_0*(1.71D0*thisOctal%UV(subcell))+T_CMB**5
           
           DO J=1,NRAYS ! Loop over rays
              
              !        The minimum dust temperature is related to the incident FUV flux along each ray
              !        Convert the incident FUV flux from Draine to Habing units by multiplying by 1.7
              T_0=12.2*(1.71D0*thisOctal%RADSURFACE(subcell, j))**0.2
              
!!$!        Attenuate the FIR radiation produced in the surface layer
!!$         IF(PARTICLE(P)%TOTAL_COLUMN(J).GT.R_0) THEN
!!$            T_0=T_0*(PDR(P)%TOTAL_COLUMN(J)/R_0)**(-0.4)
!!$         END IF
              
              !        Add the contribution to the dust temperature from the FUV flux incident along this ray
              IF(T_0.GT.0) thisOctal%DUST_T(subcell)=thisOctal%DUST_T(subcell) &
                   & + (0.42-LOG(3.45D-2*TAU_100*T_0))*(3.45D-2*TAU_100*T_0)*T_0**5
              
           END DO ! End of loop over rays
           
           !     Convert from total dust emission intensity to dust temperature
           thisOctal%DUST_T(subcell)=thisOctal%DUST_T(subcell)**0.2
           
           !#ifdef XRAYS2
           !        STOP '[dust_t.F90] not coded for XRAYS=ON'
           !!       !     Calculate the contribution to the dust temperature from the local X-ray flux (assuming a fixed grain abundance of 1.6E-8)
           !       PDR(P)%DUST_T=PDR(P)%DUST_T+1.5D4*(PDR(P)%XRAY_ENERGY_DEPOSITION_RATE/1.6D-8)**0.2
           !#endif
           !     Impose a lower limit on the dust temperature, since values below 10 K can dramatically
           !     limit the rate of H2 formation on grains (the molecule cannot desorb from the surface)
           IF(thisOctal%DUST_T(subcell).LT.10) THEN
              thisOctal%DUST_T(subcell)=10.0D0
           END IF
           
           !     Check that the dust temperature is physical
           IF(thisOctal%DUST_T(subcell).GT.1000) THEN
              WRITE(6,*) 'ERROR! Calculated dust temperature exceeds 1000 K'
              STOP
           END IF
        else
           thisOctal%Dust_T(subcell) = 0.d0
        end if

     end if
  END DO ! End of loop over particles
     
END SUBROUTINE CALCULATE_DUST_TEMPERATURES
!=======================================================================




end module pdr_mod
#endif
