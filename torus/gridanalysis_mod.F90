module gridanalysis_mod
  use kind_mod
  use gridtype_mod
  use amr_utils_mod
  use random_mod
  use vector_mod
  use source_mod
  implicit none

contains

  subroutine analysis(grid)
    use inputs_mod, only : imodel
    use utils_mod, only : findMultifilename
    type(GRIDTYPE) :: grid
    real(double) :: mass, mass14, mass15, mass16, mdisc
    real(double) :: flux
    character(len=80) :: thisFile

    flux = fluxThroughSphericalSurface(grid, 1500.d0*autocm/1.d10)
    flux = flux / msol * 365.25d0*24.d0*3600.d0
    mass = 0.d0
    call findMassWithBounds(grid%octreeRoot, mass, maxRadius=1500.d0*autocm/1.d10)

    mass14 = 0.d0
    mass15 = 0.d0
    mass16 = 0.d0
    mDisc = 0.d0
    call findMassWithBounds(grid%octreeRoot, mass14, minRho=1.d-14)
    call findMassWithBounds(grid%octreeRoot, mass15, minRho=1.d-15)
    call findMassWithBounds(grid%octreeRoot, mass16, minRho=1.d-16)
    call findDiscMass(grid%octreeRoot, globalSourceArray(1)%mass, mDisc, 1500.d0*autocm/1.d10)
    write(*,*) "disc mass ",mdisc/msol
    call findMultifilename("vel****.dat",iModel,thisfile)
    open(69, file=thisFile, status="unknown", form="formatted")
    call writeVelocityWithBounds(grid%octreeRoot, maxRadius=1500.d0*autocm/1.d10)
    close(69)
    call findMultifilename("kep****.dat",iModel,thisfile)
    call  writeKeplerianButterfly(thisfile,globalSourceArray(1)%mass, 26.d0*autocm, 1500.d0*autocm)
    write(*,'(i4,1p,11e12.3)') iModel, globalSourceArray(1)%mass/msol, globalSourceArray(1)%mdot/msol*356.25d0*24.d0*3600.d0, &
flux, mass/msol, mass14/msol, mass15/msol, mass16/msol, mdisc/msol
    write(48,'(i4,1p,11e12.3)') iModel, globalSourceArray(1)%mass/msol, globalSourceArray(1)%mdot/msol*356.25d0*24.d0*3600.d0, &
flux, mass/msol, mass14/msol, mass15/msol, mass16/msol, mdisc/msol


  end subroutine analysis

  subroutine calculateColumn(columnDensity, cs, inputOctal, inputSubcell, grid)
    type(VECTOR) :: rVec, direction
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, inputOctal
    real(double) :: ds,cs,totweight,columndensity,mu,pressure
    integer :: subcell, inputsubcell


    columnDensity = 0.d0
    cs = 0.d0
    totweight = 0.d0

    direction = VECTOR(0.d0, 0.d0, 1.d0)
    thisOctal => inputOctal
    subcell = inputSubcell
    rVec = subcellCentre(inputOctal, inputSubcell)
    do while(inOctal(grid%octreeRoot, rVec))
       call findSubcellLocal(rVec, thisOctal, subcell)
       call distanceToCellBoundary(grid, rVec, direction, ds, sOctal=thisOctal)
       columnDensity = columnDensity + thisOctal%subcellSize * thisOctal%rho(subcell) * 1.d10
       mu = 1.27
       pressure = thisOctal%temperature(subcell) * kerg * thisOctal%rho(subcell) / (mu * mHydrogen)
       cs = cs + sqrt(pressure / thisOctal%rho(subcell)) * thisOctal%rho(subcell)
       totweight = totweight + thisOctal%rho(subcell)
       rVec = rVec + (ds + 1.d-6)*direction
    end do
    direction = VECTOR(0.d0, 0.d0, -1.d0)
    thisOctal => inputOctal
    subcell = inputSubcell
    rVec = subcellCentre(inputOctal, inputSubcell)
    do while(inOctal(grid%octreeRoot, rVec))
       call findSubcellLocal(rVec, thisOctal, subcell)
       call distanceToCellBoundary(grid, rVec, direction, ds, sOctal=thisOctal)
       columnDensity = columnDensity + thisOctal%subcellSize * thisOctal%rho(subcell) * 1.d10
       mu = 1.27
       pressure = thisOctal%temperature(subcell) * kerg * thisOctal%rho(subcell) / (mu * mHydrogen)
       cs = cs + sqrt(pressure / thisOctal%rho(subcell)) * thisOctal%rho(subcell)
       totweight = totweight + thisOctal%rho(subcell)
       rVec = rVec + (ds + 1.d-6)*direction
    end do
    cs = cs / totweight
  end subroutine calculateColumn

  real(double) function massWithinR(thisR, grid)
    real(double) :: thisR
    type(GRIDTYPE) :: grid

    massWithinR = 0.d0
    call findMassWithBounds(grid%octreeRoot, massWithinR, maxRadius=thisR)
  end function massWithinR

  recursive subroutine calculateToomreQ(thisOctal, grid)
    use inputs_mod, only : smallestcellsize
    use ion_mod
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(VECTOR) :: rVec
    real(double) :: r, cs, omegaK, mass, sigma, toomreQ, mu, pressure
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateToomreQ(child, grid)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal, subcell)
          if (abs(abs(rVec%z) - thisOCtal%subcellSize/2.d0) < 1496.) then
             call calculateColumn(sigma, cs, thisOctal, subcell, grid)
             mass = massWithinR(modulus(rVec), grid)
             mass = mass + SUM(globalSourceArray(1:GlobalnSource)%mass)
             r = modulus(rVec)
             omegaK = sqrt(bigG * mass / (r*1d10)**3)

             mu = returnMu(thisOctal, subcell, globalIonArray, nGlobalIon)
             pressure = thisOctal%temperature(subcell) * kerg * thisOctal%rho(subcell) / (mu * mHydrogen)
             cs = sqrt(pressure / thisOctal%rho(subcell))
             toomreQ = (cs * omegaK) / (pi * bigG * sigma)
             
             if (.not.associated(thisOctal%etaCont)) allocate(thisOctal%etaCont(1:thisOctal%maxChildren))
             thisOctal%etaCont(subcell) = toomreQ
             write(*,*) "toomre Q " , toomreQ
          endif
       endif
    enddo
  end subroutine calculateToomreQ

          



  function fluxThroughSphericalSurface(grid, radius) result(totflux)
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: radius,da,flux
    integer :: i, n
    real(double) ::  totflux
    type(VECTOR) :: rVec, inVec
    n = 100000
    da = fourPi*(radius*1.d10)**2/dble(n)
    totflux = 0.d0
    thisOctal => grid%octreeRoot
    do i = 1, n

       rVec = radius*randomUnitVector()

       inVec = (-1.d0)*rVec
       call normalize(inVec)

       call findSubcellLocal(rVec, thisOctal, subcell)
       flux = invec%x * thisOctal%rhou(subcell) + &
            invec%y * thisOctal%rhov(subcell) + &
            invec%z * thisOctal%rhow(subcell)
       totFlux = totFlux + flux * da
    enddo
  end function fluxThroughSphericalSurface

  recursive subroutine findMassWithBounds(thisOctal, totalMass, minRho, maxRho, minRadius, maxRadius)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(VECTOR) :: rVec
    real(double) :: totalMass
    real(double),optional :: minRho, maxRho, minRadius, maxRadius
    real(double) :: dv,r
    integer :: subcell, i
    logical :: includeThisCell

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findMassWithBounds(child, totalMass, minRho, maxRho, minRadius, maxRadius)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
             dv = cellVolume(thisOctal, subcell)*1.d30
             rVec = subcellCentre(thisOctal, subcell) 
             r = modulus(rVec)

             includeThisCell = .true.

             if (PRESENT(minRadius)) then
                if (r < minRadius) includeThisCell = .false.
             endif

             if (PRESENT(maxRadius)) then
                if (r > maxRadius) includeThisCell = .false.
             endif

             if (PRESENT(minRho)) then
                if (thisOctal%rho(subcell) < minRho) includeThisCell = .false.
             endif

             if (PRESENT(maxRho)) then
                if (thisOctal%rho(subcell) > maxRho) includeThisCell = .false.
             endif

             
             if (includethiscell) totalMass = totalMass + thisOctal%rho(subcell) * dv
          endif
       endif
    enddo
  end subroutine findMassWithBounds

  recursive subroutine findDiscMass(thisOctal, sourceMass, totalMass,  maxRadius)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(VECTOR) :: rVec, rHat, vHat, vel, vkep
    real(double) :: totalMass
    real(double) :: maxRadius, sourceMass
    real(double) :: dv,r
    integer :: subcell, i
    logical :: includeThisCell

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findDiscMass(child, sourceMass, totalMass, maxRadius)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
             dv = cellVolume(thisOctal, subcell)*1.d30
             rVec = subcellCentre(thisOctal, subcell) 
             r = modulus(rVec)

             includeThisCell = .true.

             if (r > maxRadius) includeThisCell = .false.
             vel = VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell),thisOctal%rhov(subcell)/thisOctal%rho(subcell), 0.d0)
             rHat = rVec
             call normalize(rHat)
             vHat = rHat.cross.zhat
             call normalize(vHat)
             vKep = sqrt(bigG*sourceMass/(r*1.d10))*vHat
             if (modulus(vel-vKep)/modulus(vKep) > 0.2d0) then
                includeThisCell = .false.
             endif
             if (includethiscell) totalMass = totalMass + thisOctal%rho(subcell) * dv
          endif
       endif
    enddo
  end subroutine findDiscMass

  recursive subroutine writeVelocityWithBounds(thisOctal, minRho, maxRho, minRadius, maxRadius)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(VECTOR) :: rVec
    real(double) :: r
    real(double),optional :: minRho, maxRho, minRadius, maxRadius
    real(double) :: dv
    integer :: subcell, i
    logical :: includeThisCell

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call writeVelocityWithBounds(child, minRho, maxRho, minRadius, maxRadius)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
             dv = cellVolume(thisOctal, subcell)*1.d30
             rVec = subcellCentre(thisOctal, subcell) 

             r = modulus(rVec)

             includeThisCell = .false.

             if ((abs(rVec%z)-thisOctal%subcellSize/2.d0)  < 1.d0) then
                if ((abs(rVec%y)-thisOctal%subcellSize/2.d0)  < 1.d0) then
                   includeThisCell = .true.
                endif
             endif


             if (PRESENT(minRadius)) then
                if (r < minRadius) includeThisCell = .false.
             endif

             if (PRESENT(maxRadius)) then
                if (r > maxRadius) includeThisCell = .false.
             endif

             if (PRESENT(minRho)) then
                if (thisOctal%rho(subcell) < minRho) includeThisCell = .false.
             endif

             if (PRESENT(maxRho)) then
                if (thisOctal%rho(subcell) > maxRho) includeThisCell = .false.
             endif

             if (includethiscell) write(69,*) rVec%x/1496.d0, &
                  thisOctal%rhov(subcell)/thisOctal%rho(subcell)/1.d5
          endif
       endif
    enddo
  end subroutine writeVelocityWithBounds

  subroutine writeKeplerianButterfly(thisfile,mass, rMin, rMax)
    real(double) :: mass, rMin, rMax, r, v
    character(len=*) :: thisFile
    integer :: i
    open(66, file=thisFile, status="unknown", form="formatted")
    do i = 1, 1000
       r = rMax + (rMin-rMax)*dble(i-1)/999.d0
       v = sqrt(bigG * mass/r)
       write(66,*) -r/autocm,v/1.d5
    enddo
    do i = 1, 1000
       r = rMin + (rMax-rMin)*dble(i-1)/999.d0
       v = -sqrt(bigG * mass/r)
       write(66,*) r/autocm,v/1.d5
    enddo
    v = sqrt(bigG * mass/rMax)
    write(66,*) -rMax/autocm,v/1.d5

    close(66)
  end subroutine writeKeplerianButterfly

end module gridanalysis_mod
