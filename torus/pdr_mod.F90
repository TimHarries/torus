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

implicit none

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
  call setupPDR()

!do the ray casting
  call writeInfo("Casting rays over grid.", TRIVIAL)
  call castAllRaysOverGrid(grid%octreeRoot, grid)
  call writeVTKfile(grid, "columnDensity.vtk", valueTypeString=(/"rho       ",&
       "columnRho "/))

end subroutine PDR_MAIN

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
  type(vector) :: testPosition, rVec, startPosition, uhat
  real(double) :: tVal
  integer :: nside, i 

  nside = 2**hlevel
  nrays = 12*nside**2

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
        if(.not. thisOctal%ionfrac(subcell, 2) > 0.9d0) then
        do i = 0, nrays-1
           
           !transfer healpix vectors to the test position
           !     rVec = startPosition + VECTOR(vectors(1, i), vectors(2, i), vectors(3, i))
           rVec = VECTOR(vectors(1, i), vectors(2, i), vectors(3, i))
!           mag = modulus(rVec)
           uhat = rVec
           call normalize(uhat)
!           print *, "uhat", uhat%x
!           totalLength = 0.d0
!           totalRho = 0.d0
           testPosition = startPosition
           
           do while (inOctal(grid%octreeRoot, testPosition))
              
              call findSubcellLocal(testPosition, sOctal, ssubcell)
   !           print *, " "
!              print *, "DTCB "
!              print *, "ray ", i
!              print *, testposition
!              print *, uhat              
!              print *, sOctal%temperature(ssubcell)
!              print *, " "
              tval = 0.d0
!              print *, "depth ", soctal%ndepth
             call distanceToCellBoundary(grid, testPosition, uHat, tVal, soctal, ssubcell)
 
   !              print *, "TVAL ", tval/(sqrt(2.d0)*2.d0*grid%halfsmallestsubcell)
 
              if(thisOctal%ionFrac(subcell,2) > 0.9d0) exit

              testPosition = testPosition + ((tVal+1.d-10*grid%halfsmallestsubcell)*uhat)

!              print *, "pre", sOctal%columnRho(ssubcell)
              thisOctal%columnRho(subcell) = thisOctal%columnRho(subcell) + (sOctal%rho(ssubcell)*tval)
 !             print *, "post ", sOctal%columnRho(ssubcell)
  !            print *, sOctal%Rho(ssubcell)


!              totalLength = totalLength + (tVal)

              
!              totalRho = totalRho + sOctal%rho(ssubcell)*tval              
             
           end do
           
        end do
     else
        thisOctal%columnRho(subcell) = 0.d0
     end if

!        stop
        thisOctal%columnRho(subcell) = thisOctal%columnRho(subcell)/dble(nrays)

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


end module pdr_mod
#endif
