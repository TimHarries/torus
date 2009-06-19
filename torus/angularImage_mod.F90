module angularImage

  use gridtype_mod, only: GRIDTYPE
  use molecular_mod, only:  moleculetype
  use kind_mod
  use vector_mod
  use datacube_mod, only: datacube
  use constants_mod
  use messages_mod
  use timing
  
  implicit none 

  public :: make_angular_image

  private 

  type(VECTOR) :: observerVelocity
  character(len=200) :: message
  type(VECTOR) :: rayposition
  real(double) :: this_gal_lat, this_gal_lon

  contains

    subroutine make_angular_image(grid)

      use datacube_mod, only: DATACUBE, initCube, addVelocityAxis, writeDataCube, freeDataCube
      use gridtype_mod, only: GRIDTYPE
      use amr_mod, only: amrGridVelocity
      use h21cm_mod, only: h21cm_lambda
      use input_variables, only: intPosX, intPosY, intPosZ, npixels, nv, minVel, maxVel

      implicit none

      TYPE(gridtype), intent(in) :: grid
      type(DATACUBE) ::  cube
      type(MOLECULETYPE) :: thisMolecule

! molecular weight is used for column density calculation
      thisMolecule%molecularWeight = mHydrogen / amu

! Set up 21cm line
      allocate( thisMolecule%transfreq(1) )
      thisMolecule%transfreq(1) = cSpeed / (h21cm_lambda)

! Set up the observer's position 
      rayposition = VECTOR(intPosX, intPosY, intPosZ)
      write(message,'(a,3(ES12.3,2x),a)') "Observer's position is ", rayposition, "(x10^10cm)" 
      call writeinfo(message, TRIVIAL)

! Get the observer's velocity from the grid
      observerVelocity = amrGridVelocity(grid%octreeRoot, rayposition, linearinterp = .false.)
      write(message,*) "Observer's velocity is ", observerVelocity * (cspeed / 1.0e5), "km/s"
      call writeinfo(message, TRIVIAL)

      call writeinfo("Initialising datacube",TRIVIAL)
      call initCube(cube, npixels, npixels, nv)
      call addvelocityAxis(cube, minVel, maxVel) 

      call writeinfo("Generating internal view", TRIVIAL)
      call createAngImage(cube, grid, thisMolecule)

      call writeinfo("Writing data cube", TRIVIAL)
      if(writeoutput) call writedatacube(cube, "theGalaxy.fits")

      call freeDataCube(cube)
     
    end subroutine make_angular_image

!-----------------------------------------------------------------------------------------------------------

    subroutine createAngImage(cube, grid, thisMolecule)

      use input_variables, only : gridDistance, npixels, nv, imageside,  &
           minVel, maxVel, nsubpixels, splitCubes
      use molecular_mod, only: calculateOctalParams
      use atom_mod, only: bnu
      use vector_mod
      use h21cm_mod, only: h21cm_lambda

#ifdef MPI
      use mpi_global_mod, only: myRankGlobal, nThreadsGlobal
#endif

      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
     type(MOLECULETYPE) :: thisMolecule
     type(GRIDTYPE) :: grid
     type(DATACUBE) :: cube
     real(double) :: deltaV
     integer :: iTrans = 1
     integer :: iv
     real(double) :: intensitysum
     real(double), save :: background
     real, allocatable :: temp(:,:,:) 
     integer :: ix1, ix2
     integer, parameter :: temp_dim=5 ! number of elements in temp array

#ifdef MPI
     ! For MPI implementations
     integer :: ierr, n           ! error flag
     integer :: itemp ! loop counter for MPI communication
     real(double), allocatable :: tempArray(:), tempArray2(:)
#endif


! Divide up the image along the x axis for MPI case, otherwise work on the whole image
#ifdef MPI
     ix1 = (myRankGlobal)   * (cube%nx / (nThreadsGlobal)) + 1
     ix2 = (myRankGlobal+1) * (cube%nx / (nThreadsGlobal))
     if (myRankGlobal == (nThreadsGlobal-1)) ix2 = cube%nx

     n = (npixels*npixels)
     allocate(tempArray(1:n), tempArray2(1:n))
#else
     ix1 = 1
     ix2 = npixels
#endif

     deltaV = minVel * 1.e5/cspeed_sgl
     
     allocate(temp(npixels,npixels,temp_dim))

     do iv = 1,nv

        deltaV = (cube%vAxis(iv)*1.e5/cSpeed_sgl) ! velocities in fraction of c
        
        if(writeoutput) then
           write(message,*) "Done ",iv," velocity"
           call tune(6, message)  ! start a stopwatch
        endif

        if(iv .eq. 1) then
           call writeinfo("Filling Octal parameters for first time",TRIVIAL)
           call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule, deltaV)
           call writeinfo("Done filling Octal parameters for first time",TRIVIAL)
        endif

        temp(:,:,:) = 0.d0
        call makeAngImageGrid(grid, cube, thisMolecule, itrans, deltaV, nSubpixels, temp, ix1, ix2)

#ifdef MPI
        do itemp = 1, temp_dim
           tempArray = reshape(temp(:,:,itemp), (/ n /))
           call MPI_ALLREDUCE(tempArray,tempArray2,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr) 
           temp(:,:,itemp) = reshape(tempArray2, (/ npixels, npixels /))
        end do
#endif        

! Intensity as brightness temperature
        cube%intensity(:,:,iv) = real(temp(:,:,1)) * (h21cm_lambda**2) / (2.0 * kErg)
        cube%tau(:,:,iv)       = real(temp(:,:,2))
        cube%nCol(:,:)         = real(temp(:,:,3)) 
        if ( splitCubes ) then 
           cube%i0_pos(:,:,iv) = real(temp(:,:,4))
           cube%i0_neg(:,:,iv) = real(temp(:,:,5))
        end if

        if(writeoutput) then
           call tune(6, message)  ! stop a stopwatch
        endif

        intensitysum = sum(temp(:,:,1)) / dble(npixels**2)

        if(iv .eq. 1) then
           background = Bnu(thisMolecule%transfreq(itrans), Tcbr)
           write(message, *) "Background Intensity: ",background
           call writeinfo(message, TRIVIAL)
        endif

        write(message,'(a,es11.4e1,tr3,a,f10.4,tr3,a,es12.4)') &
             "DELTAV(v/c):",deltaV," V (km/s):",real(cube%vAxis(iv)), "Average Intensity:",intensitysum
        call writeinfo(message,FORINFO)

     end do

#ifdef MPI
     deallocate(tempArray, tempArray2)
#endif

   end subroutine createAngImage

!-----------------------------------------------------------------------------------------------------------


    subroutine makeAngImageGrid(grid, cube, thisMolecule, itrans, deltaV, nSubpixels, imagegrid, ix1, ix2)

      use input_variables, only : npixels, imageside, centrevecx, centrevecy
      use vector_mod
      
      type(GRIDTYPE), intent(IN) :: grid
      type(datacube), intent(IN) :: cube
      type(MOLECULETYPE), intent(IN) :: thisMolecule
      integer, intent(IN) :: itrans
      real(double), intent(IN) :: deltaV
      integer, intent(IN) :: nsubpixels
      type(VECTOR) :: viewvec !, ObserverVec
      real(double) :: viewvec_x, viewvec_y, viewvec_z
      real, intent(OUT) :: imagegrid(:,:,:)
      integer, intent(in) :: ix1, ix2

      integer :: subpixels
      integer :: ipixels, jpixels

      real :: theta_min
      real :: phi_min 
      real(double) :: delta_theta, delta_phi
      real(double) :: theta_axis(npixels), phi_axis(npixels)

      if (nsubpixels .gt. 0) then ! if nsubpixels = 0 then use adaptive subpixel sampling
         subpixels = nsubpixels
      else
         subpixels = 0
      endif

      theta_min = ( 90.0 - centrevecy ) - 0.5 * imageside 
      phi_min   = ( centrevecx - 90.0 ) - 0.5 * imageside

      delta_theta = imageside / real(npixels)
      delta_phi   = imageside / real(npixels)

! Set up axis arrays
      do jpixels=1, npixels
         theta_axis(jpixels) = theta_min + ( real(npixels - jpixels + 1) * delta_theta )
      end do
      theta_axis(:) = theta_axis(:) * degToRad

! Reverse the order in which the longitude bins are populated 
      do ipixels=1, npixels
         phi_axis(ipixels) = phi_min + ( real(npixels - ipixels + 1) * delta_phi )
      end do
      phi_axis(:) = phi_axis(:) * degToRad

! Convert to galactic co-ordinates for writing out to the FITS cube.
      cube%xAxis(:) = 90.0 + ( phi_axis(:)   / degToRad )
      cube%yAxis(:) = 90.0 - ( theta_axis(:) / degToRad )

      do jpixels = 1, npixels ! raster over image
         do ipixels = ix1, ix2

            this_gal_lon = cube%xAxis(ipixels)
            this_gal_lat = cube%yAxis(jpixels)

            viewvec_x = sin( theta_axis(jpixels) ) * cos( phi_axis(ipixels) ) 
            viewvec_y = sin( theta_axis(jpixels) ) * sin( phi_axis(ipixels) ) 
            viewvec_z = cos( theta_axis(jpixels) )

            viewvec = VECTOR( viewvec_x, viewvec_y, viewvec_z )
! Take out effect of rotating galaxy away from axes
            viewvec  = rotateY( viewvec, 45.0*degToRad )
            call normalize(viewvec)

            imagegrid(ipixels,jpixels,:) = &
                 AngPixelIntensity(viewvec,grid,thisMolecule,iTrans,deltaV, subpixels)

         enddo

      enddo

    end subroutine makeAngImageGrid

 !!! Calculates the intensity for a square pixel of arbitrary size, position, orientation
 function AngPixelIntensity(viewvec,grid,thisMolecule,iTrans,deltaV,subpixels) result(out)
   
   use input_variables, only : tolerance, nsubpixels

   type(GRIDTYPE) :: grid
   type(MOLECULETYPE) :: thisMolecule
   integer :: itrans
   type(VECTOR), intent(in) :: viewVec
   real(double) :: i0, opticaldepth, i0_pos, i0_neg

   type(VECTOR) :: thisViewVec

   integer :: subpixels, minrays
   integer :: iray

   real(double) :: avgIntensityNew, avgIntensityOld
   real(double) :: varIntensityNew, varIntensityOld
   real(double) :: avgNColNew, avgNColOld
   real(double) :: nCol

   logical :: converged
   real(double) :: deltaV

   real(double) :: out(5) 

   avgIntensityOld = 0.
   varIntensityOld = 0.
   avgNColOld      = 0.0

   converged = .false. ! failed flag
     
   if(subpixels .ne. 0) then 
      converged = .true.
      minrays = nsubpixels
   else
      minrays = -1
   endif
   
   iray = 1
   
   do while((.not. converged) .or. (iray .le. minrays))  

      thisViewVec = viewVec

     call intensityalongrayRev(rayposition,thisViewVec,grid,thisMolecule,itrans,deltaV,i0,i0_pos,i0_neg, &
          tau=opticaldepth, nCol=nCol, observerVelocity=observerVelocity )

     if (isnan(i0)) then
        write(*,*) "Got nan", opticaldepth, rayposition
        i0 = 0.d0
     endif

      avgIntensityNew = ((iray - 1) * avgIntensityOld + i0) / dble(iray)
      varIntensityNew = ((iray - 1) * varIntensityOld + ((i0 - avgIntensityNew) * (i0 - avgIntensityOld))) / dble(iray)
      avgNColNew      = ((iray - 1) * avgNColOld + nCol) / dble(iray)

      if(varIntensityNew .lt. iray * (tolerance* avgIntensityNew)**2 .and. iray .gt. 1) then
         converged = .true.
         iray = iray + 1
         out(1) = avgIntensityNew
         out(2) = opticaldepth ! not averaged, just last value
         out(3) = avgNColNew
         out(4) = i0_pos
         out(5) = i0_neg
      elseif(iray .gt. 10000) then
         converged = .false.
         out(1) = avgIntensityNew
         out(2) = opticaldepth ! not averaged, just last value
         out(3) = avgNColNew
         out(4) = i0_pos
         out(5) = i0_neg
         exit
      else
         avgIntensityOld = avgIntensityNew
         varIntensityOld = varIntensityNew
         iray = iray + 1
         out(1) = avgIntensityNew
         out(2) = opticaldepth ! not averaged, just last value
         out(3) = avgNColNew
         out(4) = i0_pos
         out(5) = i0_neg
      endif

   enddo
   
 end function AngPixelIntensity

   subroutine intensityAlongRayRev(position, direction, grid, thisMolecule, iTrans, deltaV,i0,i0_pos,i0_neg,tau,tautest, &
        rhomax, i0max, nCol, observerVelocity)

     use input_variables, only : useDust, h21cm, densitysubsample, amrgridsize
     use octal_mod 
     use atom_mod, only: Bnu
     use amr_mod
     use molecular_mod, only: densite, velocity, phiprof
     use h21cm_mod, only: h21cm_lambda

     type(VECTOR) :: position, direction, dsvector, otherDirection
     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     real(double) :: disttoGrid
     integer :: itrans
     real(double) :: nMol
     real(double), intent(out) :: i0, i0_pos, i0_neg
     real(double) :: previous_i0
     real(double), optional, intent(out) :: nCol
     type(VECTOR), optional, intent(in) ::  observerVelocity
     type(OCTAL), pointer :: thisOctal
     integer :: subcell
     type(VECTOR) :: currentPosition, thisPosition, endPosition, otherSide
     type(VECTOR) :: startVel, endVel, thisVel, veldiff
     real(double) :: alphanu1, alphanu2, jnu, snu
     real(double) :: alpha
     real(double)  :: dv, deltaV, dVacrossCell
     integer :: i, icount
     real(double) :: tval
     integer :: nTau
     real(double) ::  OneOvernTauMinusOne, ds

     real(double) :: dTau, etaline, dustjnu
     real(double), intent(out), optional :: tau
     real(double),save :: BnuBckGrnd

     real(double) :: phiProfVal
     real         :: sigma_thermal
     real(double) :: deps, origdeps

     logical,save :: firsttime = .true.
     logical,save :: havebeenWarned = .false.
     character(len = 80) :: message
     logical, optional :: tautest
     logical :: dotautest

     real(double) :: dI, dIovercell, attenuateddIovercell, dtauovercell, n
     real(double), optional, intent(out) :: rhomax
     real(double), optional, intent(in) :: i0max
     real(double) :: distToObs

     if(present(tautest)) then
        dotautest = tautest
     else
        dotautest = .false.
     endif

     if(firsttime .or. dotautest) then
        BnuBckGrnd = Bnu(thisMolecule%transfreq(itrans), Tcbr)
        if(.not. dotautest) firsttime = .false.
     endif
     
     if(inOctal(grid%octreeRoot, Position)) then
        disttogrid = 0.
     else
        distToGrid = distanceToGridFromOutside(grid, position, direction) 
     endif

     if (distToGrid > 1.e29) then
        if (.not. haveBeenWarned) then
           call writewarning("ray does not intersect grid")
           write(message, *) position
           call writeinfo(message, FORINFO)
           write(message, *) direction
           call writeinfo(message, FORINFO)
           havebeenWarned = .true.

        endif
        
        i0 = 1.d-60
        tau = 1.d-60
        goto 666
     endif
     
     ! Find a position on the other side of the grid. 
     otherSide      = position + (distToGrid * direction) + (2.0 * real(amrgridsize,db) * direction)
     otherDirection%x = -1.0 * direction%x
     otherDirection%y = -1.0 * direction%y
     otherDirection%z = -1.0 * direction%z

     ! Work out distance to grid from new position 
     distToGrid = distanceToGridFromOutside(grid, otherSide, otherDirection)

     deps = 5.d-4 * grid%halfSmallestSubcell ! small value to add on to distance from grid to ensure we get onto grid
     origdeps = deps


     currentPosition = otherside + (distToGrid + deps) * otherDirection
     distToObs = (currentPosition - position) .dot. direction

     i0  = BnuBckGrnd
     tau = 0.d0
     i0_pos = BnuBckGrnd
     i0_neg = 0.d0 
     if (present(nCol)) nCol = 0.d0

     if(present(rhomax)) rhomax = 0.d0

     thisOctal => grid%octreeRoot
     icount = 0

     do while((.not. inOctal(grid%octreeRoot, currentPosition)) .and. icount .eq. 0 .and. deps .lt. 1d30)
        deps = 10.d0 * deps
        currentPosition = otherSide + (distToGrid + deps) * otherDirection
     enddo
   
     ingrid_loop: do while(inOctal(grid%octreeRoot, currentPosition) .and. distToObs > 0 )

        icount = icount + 1
           
        call findSubcelllocal(currentPosition, thisOctal, subcell)
        call distanceToCellBoundary(grid, currentPosition, otherDirection, tVal, sOctal=thisOctal)

        if(present(rhomax)) then
           rhomax = max(rhomax, thisoctal%rho(subcell))
        endif

        dtauovercell = 0.d0
        dIovercell = 0.d0
        attenuateddIovercell = 0.d0

        if(densitysubsample) then
           if ( h21cm) then
              nmol = densite(currentposition, grid) / (thisOctal%rho(subcell))
           else
              nmol = thisoctal%molabundance(subcell) * (densite(currentposition, grid) / &
                   (2.d0 * mhydrogen))
           end if
        else
           if ( h21cm ) then
              nMol = 1.0
           else
              nMol = thisOctal%molcellparam(1,subcell)
           endif
        end if

        etaline = nmol * thisOctal%molcellparam(5,subcell)

        if(usedust) then
           alphanu2 = nmol * thisOctal%molcellparam(7,subcell)
           dustjnu = nmol * thisOctal%molcellparam(8,subcell)
        else
           alphanu2 =  0.0
        endif

        thisPosition = currentPosition

        startVel = Velocity(currentPosition, grid, startoctal = thisoctal, subcell = subcell)
        endPosition = currentPosition + tval * otherDirection

        endVel = Velocity(endPosition, grid, startoctal = thisoctal, subcell = subcell)

        Veldiff = endVel - startVel

        dvAcrossCell = (veldiff.dot.direction)

        if ( h21cm ) then 
           ! Calculate line width in cm/s.
           sigma_thermal = sqrt (  (kErg * thisOctal%temperature(subcell)) / mHydrogen)
           ! Convert to Torus units (v/c)
           sigma_thermal = sigma_thermal / cspeed
           dvAcrossCell = abs(dvAcrossCell / sigma_thermal)
        else
           dvAcrossCell = abs(dvAcrossCell * thisOctal%molmicroturb(subcell))
        end if

        if(densitysubsample) then ! should replace 5 with maxdensity/mindensity * fac
           nTau = min(max(5, nint(dvAcrossCell * 5.d0)), 100) ! ensure good resolution / 5 chosen as its the magic number!
        else
           nTau = min(max(2, nint(dvAcrossCell * 5.d0)), 100) ! ensure good resolution / 5 chosen as its the magic number!
        endif
        
        OneOvernTauMinusOne = 1.d0/(nTau - 1.d0)
        
        ds = tval * OneOvernTauMinusOne
        
! Calculate column density
! Factor of 1.d10 is to convert ds to cm 
        if (present(nCol)) then 
           nCol = nCol + (thisOctal%rho(subcell) / (thisMolecule%molecularWeight * amu) ) * ds * 1.d10
        end if

        dsvector = ds * otherDirection

        previous_i0 = i0 

        ntau_loop: do i = 2, nTau 

           thisPosition = thisPosition + dsvector
           if(.not. inoctal(grid%octreeroot, thisposition)) thisPosition = thisPosition - 0.99d0 * dsvector
           thisVel = Velocity(thisPosition, grid, startoctal = thisoctal, subcell = subcell)

           if ( present(observerVelocity) ) then 
              thisVel = thisVel - observerVelocity
           end if

           dv = (thisVel .dot. direction) - deltaV

           if ( h21cm ) then 
              phiprofval = gauss (sigma_thermal, real(dv) ) / thisMolecule%transfreq(1)
           else  
              phiProfval = phiProf(dv, thisOctal%molmicroturb(subcell))
           end if

           
           if(densitysubsample .and. .not. h21cm ) then
              nmol = thisoctal%molabundance(subcell) * (Densite(thisposition, grid) / &
                     (2.d0 * mhydrogen))

              if(usedust) then
                 alphanu2 = nmol * thisOctal%molcellparam(7,subcell)
                 dustjnu =  nmol * thisOctal%molcellparam(8,subcell)
              else
                 alphanu2 = 0.0
              endif

              etaline = nmol * thisOctal%molcellparam(5,subcell)
              alphanu1 = nmol * thisOctal%molcellparam(6,subcell) * phiprofval

           else if (densitysubsample .and. h21cm) then
              nmol     = Densite(thisposition, grid) / (thisOctal%rho(subcell))
              etaline  = nmol * thisOctal%molcellparam(5,subcell)
              alphanu1 = nmol * thisOctal%molcellparam(6,subcell) * phiprofval

           else
              alphanu1 = (nmol * thisOctal%molcellparam(6,subcell)) * phiprofval
           endif

           alpha = alphanu1 + alphanu2
           dTau = alpha * ds * 1.d10
           dtauovercell = dtauovercell + dtau

           jnu = etaLine * phiProfVal

           if(useDust) jnu = jnu + dustjnu

           if (alpha .ne. 0.d0) then
              snu = jnu/alpha
           else
              snu = tiny(snu)
              i0 = i0 + tiny(i0)
           endif

           dI = (1.d0-exp(-dtau))*snu

           tau = tau + dtau

           attenuateddIovercell = attenuateddIovercell + dI

           i0 = i0 * exp(-1.0*dtau) + dI

        end do ntau_loop

        ! Change in brightness temperature over this cell
        dIovercell = (i0 - previous_i0) * (h21cm_lambda**2) / (2.0 * kErg)

        if ( dIovercell > 0 ) then 
           i0_pos = i0_pos + dIovercell
        else
           i0_neg = i0_neg + dIovercell
        end if

        if(present(i0max) .and. i0 .gt. 0.99d0 * i0max) then 
           exit
        endif

        n = thisoctal%newmolecularlevel(4,subcell)
                    
        thisoctal%newmolecularlevel(5,subcell) = (n * thisoctal%newmolecularlevel(5,subcell) + dtauovercell) / (n + 1.d0)

        ! average change in brightness temperaure over this cell
        thisoctal%newmolecularlevel(1,subcell) = (n * thisoctal%newmolecularlevel(1,subcell) + dIovercell) / (n + 1.d0) 

        ! Image co-ordinates        
        thisoctal%newmolecularlevel(2,subcell) = (n * thisoctal%newmolecularlevel(2,subcell) + this_gal_lon) / (n + 1.d0)
        thisoctal%newmolecularlevel(3,subcell) = (n * thisoctal%newmolecularlevel(3,subcell) + this_gal_lat) / (n + 1.d0)

        thisoctal%newmolecularlevel(4,subcell) = thisoctal%newmolecularlevel(4,subcell) + 1.d0

        currentPosition = currentPosition + (tval + origdeps) * otherDirection
        distToObs = (currentPosition - position) .dot. direction
        
     enddo ingrid_loop
         
666  continue

   end subroutine intensityAlongRayRev

end module angularImage
