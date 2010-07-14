module angularImage

  use kind_mod
  use vector_mod
  use constants_mod
  use messages_mod
  use gridtype_mod, only: GRIDTYPE
  use molecular_mod, only:  moleculetype
  use datacube_mod, only: datacube

  use timing, only: tune
  
  implicit none 

  public :: make_angular_image, map_dI_to_particles

  private 

  type(VECTOR) :: observerVelocity
  character(len=200) :: message
  type(VECTOR) :: rayposition
  real(double) :: this_gal_lat, this_gal_lon
! Approximate latitude range of galactic plane
  real(double), parameter :: plane_lat=1.0_db 
  integer, parameter :: vr_file_lun=81
  character(len=*), parameter :: vr_format="unformatted"
  integer :: vel_chan_num

  contains

    subroutine make_angular_image(grid)

      use datacube_mod, only: DATACUBE, initCube, addVelocityAxis, writeDataCube, freeDataCube
      use amr_mod, only: amrGridVelocity
      use h21cm_mod, only: h21cm_lambda
      use input_variables, only: intPosX, intPosY, intPosZ, npixels, nv, minVel, maxVel, intDeltaVx, intDeltaVy, intDeltaVz, &
           galaxyPositionAngle, galaxyInclination, wanttau, splitCubes, dataCubeFileName, obsVelFromGrid

      implicit none

      TYPE(gridtype), intent(in) :: grid
      type(DATACUBE) ::  cube
      type(MOLECULETYPE) :: thisMolecule
      TYPE(VECTOR) :: intVelMod  ! modify observer's velocity by this vector

! molecular weight is used for column density calculation
      thisMolecule%molecularWeight = mHydrogen / amu

! Set up 21cm line
      allocate( thisMolecule%transfreq(1) )
      thisMolecule%transfreq(1) = cSpeed / (h21cm_lambda)

! Set up the observer's position. Transform to rotated and tilted grid as required. 
      rayposition = VECTOR(intPosX, intPosY, intPosZ)
      write(message,'(a,3(ES12.3,2x),a)') "Observer's position is ", rayposition, "(x10^10cm)" 
      call writeinfo(message, TRIVIAL)
      rayposition  = rotateZ( rayposition, galaxyPositionAngle*degToRad )
      rayposition  = rotateY( rayposition, galaxyInclination*degToRad   )
      write(message,'(a,3(ES12.3,2x),a)') "Observer's position rotated to ", rayposition, "(x10^10cm)" 
      call writeinfo(message, TRIVIAL)

! Get the observer's velocity from the grid. This is already rotated and tilted.
      if ( obsVelFromGrid ) then 
         observerVelocity = amrGridVelocity(grid%octreeRoot, rayposition, linearinterp = .false.)
         write(message,*) "Observer's velocity from grid: ", observerVelocity * (cspeed / 1.0e5), "km/s"
         call writeinfo(message, TRIVIAL)
      else
         observerVelocity = VECTOR (0.0, 0.0, 0.0)
      end if

! Add velocity offset, transformed to rotated and tilted grid 
      intVelMod = VECTOR(intDeltaVx, intDeltaVy, intDeltaVz)
      intVelMod = intVelMod * (1.0e5 / cspeed)
      intVelMod = rotateZ( intVelMod, galaxyPositionAngle*degToRad )
      intVelMod = rotateY( intVelMod, galaxyInclination*degToRad   )
      observerVelocity = observerVelocity + intVelMod
      if ( obsVelFromGrid ) then
         write(message,*) "Modified observer velocity: ", observerVelocity * (cspeed / 1.0e5), "km/s"
      else
         write(message,*) "Observer's velocity set to: ", observerVelocity * (cspeed / 1.0e5), "km/s"
      end if
      call writeinfo(message, TRIVIAL)


      call writeinfo("Initialising datacube",TRIVIAL)
      call initCube(cube, npixels, npixels, nv)
      ! Reverse velocity axis 
      call addvelocityAxis(cube, maxVel, minVel) 

      call writeinfo("Generating internal view", TRIVIAL)
      call createAngImage(cube, grid, thisMolecule)

      call process_cube_for_kvis(cube)

      call writeinfo("Writing data cubes", TRIVIAL)
      if(writeoutput) then

         call writeinfo("Writing intensity to intensity_"//trim(dataCubeFileName), TRIVIAL)
         call writedatacube(cube, "intensity_"//trim(dataCubeFileName), write_Intensity=.true., &
              write_ipos=.false., write_ineg=.false., write_Tau=.false., write_nCol=.false., write_axes=.false.)

         if ( splitCubes ) then 
            call writeinfo("Writing positive intensity to intensity_pos_"//trim(dataCubeFileName), TRIVIAL)
            call writedatacube(cube, "intensity_pos_"//trim(dataCubeFileName), write_Intensity=.false., &
                 write_ipos=.true., write_ineg=.false., write_Tau=.false., write_nCol=.false., write_axes=.false.)

            call writeinfo("Writing negative intensity to intensity_neg_"//trim(dataCubeFileName), TRIVIAL)
            call writedatacube(cube, "intensity_neg_"//trim(dataCubeFileName), write_Intensity=.false., &
                 write_ipos=.false., write_ineg=.true., write_Tau=.false., write_nCol=.false., write_axes=.false.)
         end if

         if ( wanttau ) then 
            call writeinfo("Writing optical depth to tau_"//trim(dataCubeFileName), TRIVIAL)
            call writedatacube(cube, "tau_"//trim(dataCubeFileName), write_Intensity=.false., &
                 write_ipos=.false., write_ineg=.false., write_Tau=.true., write_nCol=.false., write_axes=.false.)
         end if

         call writeinfo("Writing column density to nCol_"//trim(dataCubeFileName), TRIVIAL)
         call writedatacube(cube, "nCol_"//trim(dataCubeFileName), write_Intensity=.false., &
              write_ipos=.false., write_ineg=.false., write_Tau=.false., write_nCol=.true., write_axes=.false.)

      end if

      call freeDataCube(cube)
     
    end subroutine make_angular_image

!-----------------------------------------------------------------------------------------------------------

    subroutine createAngImage(cube, grid, thisMolecule)

      use input_variables, only : npixels, nv, nsubpixels, splitCubes, wantTau
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

     allocate(temp(npixels,npixels,temp_dim))

     open (unit=vr_file_lun, file="ray_info.dat", status="replace", form=vr_format)
     do iv = 1,nv

        vel_chan_num = iv

        deltaV = (cube%vAxis(iv)*1.e5/cSpeed_sgl) ! velocities in fraction of c
        
        if(writeoutput) then
           write(message,*) "Done ",iv," velocity"
           call tune(6, message)  ! start a stopwatch
        endif

        if(iv .eq. 1) then
           call writeinfo("Filling Octal parameters for first time",TRIVIAL)
           call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule)
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
        if (wantTau ) cube%tau(:,:,iv)       = real(temp(:,:,2))
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
     close (vr_file_lun)

#ifdef MPI
     deallocate(tempArray, tempArray2)
#endif


   end subroutine createAngImage

!-----------------------------------------------------------------------------------------------------------


    subroutine makeAngImageGrid(grid, cube, thisMolecule, itrans, deltaV, nSubpixels, imagegrid, ix1, ix2)

      use input_variables, only: npixels, imageside, centrevecx, centrevecy, galaxyPositionAngle, galaxyInclination 
      use input_variables, only: nv
      use vector_mod
      
      type(GRIDTYPE), intent(IN) :: grid
      type(datacube), intent(IN) :: cube
      type(MOLECULETYPE), intent(IN) :: thisMolecule
      integer, intent(IN) :: itrans
      real(double), intent(IN) :: deltaV
      integer, intent(IN) :: nsubpixels
      type(VECTOR) :: viewvec
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

!$OMP PARALLEL default(none), private(ipixels, jpixels, this_gal_lon, this_gal_lat, viewvec),  &
!$OMP private(viewvec_x, viewvec_y, viewvec_z), &
!$OMP shared(npixels, galaxyPositionAngle, galaxyInclination, grid, thisMolecule), &
!$OMP shared(ix1, ix2, cube, theta_axis, phi_axis, iTrans, deltaV, subpixels, imagegrid, vel_chan_num, nv)
      do jpixels = 1, npixels ! raster over image
!$OMP DO
         do ipixels = ix1, ix2

            this_gal_lon = cube%xAxis(ipixels)
            this_gal_lat = cube%yAxis(jpixels)

! Write profile information in ASCII or binary as required. 
            if ( vel_chan_num == nv .and. abs(this_gal_lat) < plane_lat ) then 
               if ( vr_format == "unformatted" ) then 
                  write(vr_file_lun) real(ipixels,db), real(jpixels,db), -1.0e30_db, -1.0e30_db , -1.0e30_db , -1.0e30_db
               else
                  write(vr_file_lun,*) 
                  write(vr_file_lun,*) "#  ", ipixels, jpixels
                  write(vr_file_lun,*) 
               end if
            end if

            viewvec_x = sin( theta_axis(jpixels) ) * cos( phi_axis(ipixels) ) 
            viewvec_y = sin( theta_axis(jpixels) ) * sin( phi_axis(ipixels) ) 
            viewvec_z = cos( theta_axis(jpixels) )

            viewvec = VECTOR( viewvec_x, viewvec_y, viewvec_z )
! Take out effect of rotating galaxy away from axes in read_galaxy_sph_data
! Apply same transformation to viewvec as was carried out on particles
            viewVec  = rotateZ( viewvec, galaxyPositionAngle*degToRad )
            viewvec  = rotateY( viewvec, galaxyInclination*degToRad )
            call normalize(viewvec)

            imagegrid(ipixels,jpixels,:) = &
                 AngPixelIntensity(viewvec,grid,thisMolecule,iTrans,deltaV, subpixels)

         enddo
!$OMP END DO
      enddo
!$OMP END PARALLEL

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

   subroutine intensityAlongRayRev(position, direction, grid, thisMolecule, iTrans, deltaV,i0,i0_pos,i0_neg,tau, &
        rhomax, i0max, nCol, observerVelocity)

     use input_variables, only : useDust, h21cm, densitysubsample, nv
     use octal_mod, only: OCTAL
     use atom_mod, only: Bnu
     use amr_mod, only: inOctal, distanceToGridFromOutside, distanceToCellBoundary, findSubcelllocal
     use molecular_mod, only: interpolated_Density, velocity, phiprof
     use h21cm_mod, only: h21cm_lambda
     use utils_mod, only: gauss

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

     logical,save :: havebeenWarned = .false.
     character(len = 80) :: message

     real(double) :: dI, dIovercell, attenuateddIovercell, dtauovercell, n
     real(double), optional, intent(out) :: rhomax
     real(double), optional, intent(in) :: i0max
     real(double) :: distToObs, gridSize

     BnuBckGrnd = Bnu(thisMolecule%transfreq(itrans), Tcbr)

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
     gridsize       = grid%octreeRoot%xMax - grid%octreeRoot%xMin
     otherSide      = position + (distToGrid * direction) + (2.0 * gridsize * direction)
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
              nmol = interpolated_Density(currentposition, grid) / (thisOctal%rho(subcell))
           else
              nmol = thisoctal%molabundance(subcell) * (interpolated_Density(currentposition, grid) / &
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
           nCol = nCol + (thisOctal%rho(subcell) / (thisMolecule%molecularWeight * amu) ) * tval * 1.d10
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
              nmol = thisoctal%molabundance(subcell) * (interpolated_Density(thisposition, grid) / &
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
              nmol     = interpolated_Density(thisposition, grid) / (thisOctal%rho(subcell))
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

        ! average change in brightness temperature (per column density) over this cell
        thisoctal%newmolecularlevel(1,subcell) = (n * thisoctal%newmolecularlevel(1,subcell) + (dIovercell/nCol) ) / (n + 1.d0) 

        ! Image co-ordinates        
        thisoctal%newmolecularlevel(2,subcell) = (n * thisoctal%newmolecularlevel(2,subcell) + this_gal_lon) / (n + 1.d0)
        thisoctal%newmolecularlevel(3,subcell) = (n * thisoctal%newmolecularlevel(3,subcell) + this_gal_lat) / (n + 1.d0)

        thisoctal%newmolecularlevel(4,subcell) = thisoctal%newmolecularlevel(4,subcell) + 1.d0

        call write_los_info

        currentPosition = currentPosition + (tval + origdeps) * otherDirection
        distToObs = (currentPosition - position) .dot. direction
        
     enddo ingrid_loop
         
666  continue

     contains

       subroutine write_los_info

         implicit none

         real(double) :: this_pos, this_vel, rhoH2, H2frac

        ! Write out los profile information before update to currentPosition
         if ( vel_chan_num == nv .and. abs(this_gal_lat) < plane_lat ) then 

            this_pos = ( modulus(currentPosition - rayposition) ) * (1.0e7/pcToCm) ! in kpc
            this_vel = ( (startVel-observerVelocity) .dot. direction )  * ( cSpeed / 1.0d5) ! in km/s
            ! H2 mass density
            rhoH2    = 2.0_db * mHydrogen * thisOctal%NH2(subcell)
            ! Mass fraction of H2 (consistent with Dobbs et al 2008, MNRAS, 389, 1097)
            ! N.B. thisOctal%rho is HI mass density
            H2frac   = rhoH2 / (rhoH2 + thisOctal%rho(subcell) )

            ! Write out in ASCII or binary as required.  
            if ( vr_format == "unformatted" ) then 
               write (vr_file_lun) this_pos, this_vel, thisOctal%rho(subcell), &
                    real(thisOctal%temperature(subcell),db), thisoctal%newmolecularlevel(1,subcell), &
                    H2frac
            else
               write (vr_file_lun,'(6(ES15.6,2x))') this_pos, this_vel, thisOctal%rho(subcell), &
                    thisOctal%temperature(subcell), thisoctal%newmolecularlevel(1,subcell), &
                    H2frac
            end if

         end if

       end subroutine write_los_info

   end subroutine intensityAlongRayRev

!-----------------------------------------------------------------------------------------

   subroutine process_cube_for_kvis(cube)

     TYPE(datacube), intent(inout) :: cube 

! Set up cube to be read by kvis
!
! FITS keywords for an angular image
     cube%xUnit     = "degrees"
     cube%xAxisType = "GLON-CAR"
     cube%yAxisType = "GLAT-CAR"
     cube%vAxisType = "VELO-LSR"
     cube%intensityUnit = "K (Tb)  "

! kvis assumes that the velocity axis is in m/s 
     cube%vAxis(:) = cube%vAxis(:) * 1000.0 
     cube%vUnit    = "m/s      "

   end subroutine process_cube_for_kvis

!-----------------------------------------------------------------------------------------

   subroutine map_dI_to_particles(grid)

    use input_variables, only: sphdatafilename, galaxyPositionAngle, galaxyInclination, &
         inputFileFormat
    use sph_data_class, only: sphdata, read_galaxy_sph_data, new_read_sph_data
    use octal_mod, only: octal 
    use amr_mod, only: inOctal, findSubcellTD
    use vtk_mod, only: writeVtkFile
#ifdef MPI
    use mpi_global_mod, only: myRankGlobal
#endif

    TYPE(gridtype), intent(in) :: grid
    integer :: ipart
    TYPE(vector) :: position, positionTorus, old_position
    type(OCTAL), pointer :: thisOctal
    integer :: subcell 

    real(double) :: dI, n_sample
    real(double) :: distTotorus ! conversion factor between SPH postions and Torus positions
    real(double) :: H2_frac

#ifdef MPI
    character(len=3)    :: char_my_rank
#endif
    character(len=30)   :: outfilename
    integer, parameter  :: LUIN = 10 ! unit number of output file

    call writeVtkFile(grid, "ray_info.vtk", valueTypeString=(/"dI      ", "galLon  ", "galLat  ", "crossing"/) )

! Re-read particle data which has been deleted to save memory
    select case (inputFileFormat)
    case("binary")
       call read_galaxy_sph_data(sphdatafilename)

    case("ascii")
       call new_read_sph_data(sphdatafilename)

    case default
       call writeWarning("Unrecognised file format "//inputFileFormat)

    end select

#ifdef MPI
     write(char_my_rank, '(i3)') myRankGlobal
     outfilename="particle_dI_"//TRIM(ADJUSTL(char_my_rank))//".dat"
#else
     outfilename="particle_dI.dat"
#endif

     open (unit=LUIN, status="replace", form="formatted", file=trim(outfilename))

     distTotorus = sphdata%udist / 1.0e10_db

     particles: do ipart=1, sphdata%npart

        position = VECTOR(sphData%xn(ipart), sphData%yn(ipart), sphData%zn(ipart) )
        positionTorus = distTotorus * position

! Undo effect of rotated and tilted grid. Use this position for output
        old_position  = rotateY( position,     -1.0*galaxyInclination*degToRad   )
        old_position  = rotateZ( old_position, -1.0*galaxyPositionAngle*degToRad )

        if(inOctal(grid%octreeRoot, positionTorus)) then
           
           call findSubcellTD(positionTorus,grid%octreeRoot,thisOctal,subcell)
           
           dI       =  thisOctal%newmolecularlevel(1,subcell)
           n_sample =  thisOctal%newmolecularlevel(4,subcell)
           
           ! Calculate fraction of molecular hydrogen by number (not mass)
           ! sphdata%rhon is the mass density of HI 
           if ( associated (sphdata%rhoH2) ) then 
              H2_frac = sphdata%rhoH2(ipart) / ( (2.0_db * sphdata%rhon(ipart)) + sphdata%rhoH2(ipart) )
           else
              H2_frac = 0.0
           end if

! newmolecularlevel is a floating point number so n_sample>0.99 is a reliable way of saying one or more samples.
           if ( n_sample > 0.99 ) write(LUIN,'(9(e15.8,2x),i8)') old_position, sphdata%gasmass(ipart), sphdata%hn(ipart), &
                sphdata%rhon(ipart), dI, n_sample, H2_frac, ipart

        end if

     end do particles

     close (LUIN)

! Write splash columns file
     open(unit=LUIN, status="replace", form="formatted", file="columns")
     write(LUIN,'(a)') "x"
     write(LUIN,'(a)') "y"
     write(LUIN,'(a)') "z"
     write(LUIN,'(a)') "pmass"
     write(LUIN,'(a)') "h"
     write(LUIN,'(a)') "rho"
     write(LUIN,'(a)') "dI"
     write(LUIN,'(a)') "nsample"
     write(LUIN,'(a)') "H2frac"
     write(LUIN,'(a)') "index"
     close (LUIN)


   end subroutine map_dI_to_particles

!-----------------------------------------------------------------------------------------

end module angularImage
