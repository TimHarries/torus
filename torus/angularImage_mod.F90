#ifdef MOLECULAR
module angularImage

  use kind_mod
  use vector_mod
  use constants_mod
  use messages_mod
  use gridtype_mod, only: GRIDTYPE
  use molecular_mod, only:  moleculetype
  use datacube_mod, only: datacube
  use angularImage_utils
  use timing, only: tune
  
  implicit none 

  public :: make_angular_image
#ifdef SPH
  public :: map_dI_to_particles
#endif

  private 

  type(VECTOR) :: observerVelocity
  character(len=200) :: message
  type(VECTOR) :: rayposition
  real(double) :: this_gal_lat, this_gal_lon
! Approximate latitude range of galactic plane
  real(double), parameter :: plane_lat=1.0_db 
  integer, parameter :: vr_file_lun=81
  character(len=11),save :: vr_format="unformatted"
  integer :: vel_chan_num
  logical, save :: nColOnly 

! For storing density plot of CO vs H2 column density
  integer, save, allocatable :: CO_H2_dist(:,:)
  integer, parameter :: CO_H2_dist_size=100
! Range of data values stored
  real(double), parameter, private :: minCO = 10.0
  real(double), parameter, private :: maxCO = 18.0
  real(double), parameter, private :: minH2 = 17.0
  real(double), parameter, private :: maxH2 = 22.0

!$OMP THREADPRIVATE(this_gal_lat, this_gal_lon)
!$OMP THREADPRIVATE (vr_format)

  contains

    subroutine make_angular_image(grid)

      use datacube_mod, only: DATACUBE, initCube, addVelocityAxis, freeDataCube, convertVelocityAxis
      use amr_mod, only: amrGridVelocity
      use h21cm_mod, only: h21cm_lambda, h21cm
      use inputs_mod, only: nv, minVel, maxVel, wanttau
      use molecular_mod, only: globalMolecule
#ifdef USECFITSIO
      use datacube_mod, only : writeDataCube, dataCubeFilename
#endif

      implicit none

      TYPE(gridtype), intent(in) :: grid
      type(DATACUBE) ::  cube
      type(MOLECULETYPE) :: thisMolecule
      TYPE(VECTOR) :: intVelMod  ! modify observer's velocity by this vector

! Requesting only one velocity channel is the flag for generaing multiple column densities.
! i0_pos and i0_neg are used for storage so splitCubes needs to be true
      if (nv ==1) then
         nColOnly   = .true.
         splitCubes = .true.
         allocate( CO_H2_dist(CO_H2_dist_size,CO_H2_dist_size) )
         CO_H2_dist(:,:) = 0
      else
         nColOnly = .false.
      end if

      if ( h21cm ) then 
! molecular weight is used for column density calculation
         thisMolecule%molecularWeight = mHydrogen / amu

! Set up 21cm line
         allocate( thisMolecule%transfreq(1) )
         thisMolecule%transfreq(1) = cSpeed / (h21cm_lambda)
      end if

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
      call initCube(cube, nv, splitCubes=splitCubes, wantTau=wantTau, galacticPlaneSurvey=.true.)
      ! Reverse velocity axis 
      call addvelocityAxis(cube, maxVel, minVel) 

      call writeinfo("Generating internal view", TRIVIAL)
      if ( h21cm ) then 
         call createAngImage(cube, grid, thisMolecule)
      else
         call createAngImage(cube, grid, globalMolecule)
      end if

      call convertVelocityAxis(cube, "m/s")

#ifdef USECFITSIO
      call writeinfo("Writing data cubes", TRIVIAL)
      if(writeoutput) then

         call writeinfo("Writing intensity to intensity_"//trim(dataCubeFileName), TRIVIAL)
         call writedatacube(cube, "intensity_"//trim(dataCubeFileName), write_Intensity=.true., &
              write_ipos=.false., write_ineg=.false., write_Tau=.false., write_nCol=.false., write_axes=.false.)

         if ( nColOnly ) then
            ! We're using an intensity array to write column density so fix the units
            cube%IntensityUnit = "cm^-2     "
            call writeinfo("Writing H2 column density to nCol_H2_"//trim(dataCubeFileName), TRIVIAL)
            call writedatacube(cube, "nCol_H2_"//trim(dataCubeFileName), write_Intensity=.false., &
                 write_ipos=.true., write_ineg=.false., write_Tau=.false., write_nCol=.false., write_axes=.false.)

            call writeinfo("Writing CO column density to nCol_CO_"//trim(dataCubeFileName), TRIVIAL)
            call writedatacube(cube, "nCol_CO_"//trim(dataCubeFileName), write_Intensity=.false., &
                 write_ipos=.false., write_ineg=.true., write_Tau=.false., write_nCol=.false., write_axes=.false.)
         elseif ( splitCubes ) then 
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
#endif
      call freeDataCube(cube)

      if (nColOnly) call write_CO_vs_H2
     
    end subroutine make_angular_image

!-----------------------------------------------------------------------------------------------------------

    subroutine createAngImage(cube, grid, thisMolecule)

      use inputs_mod, only : nv, nsubpixels, wantTau, itrans
      use molecular_mod, only: calculateOctalParams
      use atom_mod, only: bnu
      use vector_mod
#ifdef MPI
      use mpi_global_mod, only: myRankGlobal, nThreadsGlobal
      use mpi
#endif

      implicit none

     type(MOLECULETYPE) :: thisMolecule
     type(GRIDTYPE) :: grid
     type(DATACUBE) :: cube
     real(double) :: deltaV
     integer :: iv
     real(double) :: intensitysum
     real(double) :: background
     real(double) :: thisLambda
     real, allocatable :: temp(:,:,:) 
     integer :: ix1, ix2
     integer, parameter :: temp_dim=5 ! number of elements in temp array

#ifdef MPI
     ! For MPI implementations
     integer :: ierr, n           ! error flag
     integer :: itemp ! loop counter for MPI communication
     real, allocatable :: tempArray(:), tempArray2(:)
#endif


! Divide up the image along the x axis for MPI case, otherwise work on the whole image
#ifdef MPI
     ix1 = (myRankGlobal)   * (cube%nx / (nThreadsGlobal)) + 1
     ix2 = (myRankGlobal+1) * (cube%nx / (nThreadsGlobal))
     if (myRankGlobal == (nThreadsGlobal-1)) ix2 = cube%nx

     n = (cube%nx*cube%ny)
     allocate(tempArray(1:n), tempArray2(1:n))
#else
     ix1 = 1
     ix2 = cube%nx
#endif

     allocate(temp(cube%nx,cube%ny,temp_dim))

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

        temp(:,:,:) = 0.0
        call makeAngImageGrid(grid, cube, thisMolecule, itrans, deltaV, nSubpixels, temp, ix1, ix2)

#ifdef MPI
        do itemp = 1, temp_dim
           tempArray = reshape(temp(:,:,itemp), (/ n /))
           call MPI_ALLREDUCE(tempArray,tempArray2,n,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr) 
           temp(:,:,itemp) = reshape(tempArray2, (/ cube%nx,cube%ny /))
        end do
#endif        

! Intensity as brightness temperature
        thisLambda = cSpeed / thisMolecule%transfreq(itrans)
        cube%intensity(:,:,iv) = real(temp(:,:,1) * (thisLambda**2) / (2.0 * kErg) )

        if (wantTau ) cube%tau(:,:,iv)       = temp(:,:,2)
        cube%nCol(:,:)         = temp(:,:,3) 
        if ( splitCubes ) then 
           cube%i0_pos(:,:,iv) = temp(:,:,4)
           cube%i0_neg(:,:,iv) = temp(:,:,5)
        end if

        if(writeoutput) then
           call tune(6, message)  ! stop a stopwatch
        endif

        intensitysum = sum(temp(:,:,1)) / real(cube%nx*cube%ny)

        if(iv .eq. 1) then
           background = Bnu(thisMolecule%transfreq(itrans), Tcbr)
           write(message, '(a,es11.4,a,es11.4,a)') &
           "Background Intensity: ",background, " at " , thisLambda, " cm"
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

      use inputs_mod, only: imageside, centrevecx, centrevecy, nv
      use vector_mod
      
      type(GRIDTYPE), intent(IN) :: grid
      type(datacube), intent(INOUT) :: cube
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
      real(double) :: theta_axis(cube%ny), phi_axis(cube%nx)

      real :: imagesideX, imagesideY

! Allow for non-square spatial axes
      imagesideX = imageside * (real(cube%nx)/real(cube%ny))
      imagesideY = imageside

      if (nsubpixels .gt. 0) then ! if nsubpixels = 0 then use adaptive subpixel sampling
         subpixels = nsubpixels
      else
         subpixels = 0
      endif

      theta_min = real(( 90.0 - centrevecy ) - 0.5 * imagesideY )
      phi_min   = real(( centrevecx - 90.0 ) - 0.5 * imagesideX)

      delta_theta = imagesideY / real(cube%ny)
      delta_phi   = imagesideX / real(cube%nx)

! Set up axis arrays
      do jpixels=1, cube%ny
         theta_axis(jpixels) = theta_min + ( real(cube%ny - jpixels + 1) * delta_theta )
      end do
      theta_axis(:) = theta_axis(:) * degToRad

! Reverse the order in which the longitude bins are populated 
      do ipixels=1, cube%nx
         phi_axis(ipixels) = phi_min + ( real(cube%nx - ipixels + 1) * delta_phi )
      end do
      phi_axis(:) = phi_axis(:) * degToRad

! Convert to galactic co-ordinates for writing out to the FITS cube.
      cube%xAxis(:) = 90.0 + ( phi_axis(:)   / degToRad )
      cube%yAxis(:) = 90.0 - ( theta_axis(:) / degToRad )

!$OMP PARALLEL default(none), private(ipixels, jpixels, viewvec),  &
!$OMP private(viewvec_x, viewvec_y, viewvec_z), &
!$OMP shared(galaxyPositionAngle, galaxyInclination, grid, thisMolecule, ncolonly), &
!$OMP shared(ix1, ix2, cube, theta_axis, phi_axis, iTrans, deltaV, subpixels, imagegrid, vel_chan_num, nv)
      do jpixels = 1, cube%ny ! raster over image
!$OMP DO
         do ipixels = ix1, ix2

            this_gal_lon = cube%xAxis(ipixels)
            this_gal_lat = cube%yAxis(jpixels)

! Write profile information in ASCII or binary as required. 
            if ( vel_chan_num == nv .and. abs(this_gal_lat) < plane_lat .and. .not. nColOnly) then 
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
                 real(AngPixelIntensity(viewvec,grid,thisMolecule,iTrans,deltaV, subpixels))

         enddo
!$OMP END DO
      enddo
!$OMP END PARALLEL

    end subroutine makeAngImageGrid

 !!! Calculates the intensity for a square pixel of arbitrary size, position, orientation
 function AngPixelIntensity(viewvec,grid,thisMolecule,iTrans,deltaV,subpixels) result(out)
   
   use inputs_mod, only : tolerance, nsubpixels

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

     if ( ncolOnly ) then 
        call intensityalongrayRev(rayposition,thisViewVec,grid,thisMolecule,itrans,deltaV,i0, &
             nCol_H2=i0_pos, nCol_CO=i0_neg, tau=opticaldepth, nCol=nCol, observerVelocity=observerVelocity )
     else
        call intensityalongrayRev(rayposition,thisViewVec,grid,thisMolecule,itrans,deltaV,i0,i0_pos,i0_neg, &
             tau=opticaldepth, nCol=nCol, observerVelocity=observerVelocity )
     end if

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
        nCol, nCol_H2, nCol_CO, observerVelocity)

     use inputs_mod, only : useDust, densitysubsample, nv, vturb
     use h21cm_mod, only: h21cm
     use octal_mod, only: OCTAL
     use atom_mod, only: Bnu
     use amr_mod, only: inOctal, distanceToGridFromOutside, distanceToCellBoundary, findSubcelllocal
     use molecular_mod, only: interpolated_Density, velocity, phiprof
     use utils_mod, only: gauss

     type(VECTOR) :: position, direction, dsvector, otherDirection
     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     real(double) :: disttoGrid
     integer :: itrans
     real(double) :: nMol
     real(double), intent(out) :: i0
     real(double), intent(out), optional :: i0_pos, i0_neg
     real(double) :: previous_i0
     real(double), optional, intent(out) :: nCol
     real(double), optional, intent(out) :: nCol_H2 ! H2 column density
     real(double), optional, intent(out) :: nCol_CO ! CO column density
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
     real(double) :: BnuBckGrnd

     real(double) :: phiProfVal
     real         :: sigma_thermal
     real(double) :: deps, origdeps

     logical,save :: havebeenWarned = .false.
     character(len = 80) :: message

     real(double) :: dI, dIovercell, attenuateddIovercell, dtauovercell, n
     real(double) :: distToObs, gridSize

     logical, parameter :: doWriteLosInfo=.false. 

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
        if ( present (tau) ) tau = 1.d-60
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
     if (present(i0_pos) ) i0_pos  = BnuBckGrnd
     if (present(i0_neg) ) i0_neg  = 0.d0 
     if (present(tau) )    tau     = 0.d0
     if (present(nCol))    nCol    = 0.d0
     if (present(nCol_H2)) nCol_H2 = 0.d0
     if (present(nCol_CO)) nCol_CO = 0.d0 

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

! molmicroturb won't be defined for 21cm case so calculate the line width here
        if ( h21cm ) then 
           ! Calculate thermal line width in cm/s.
           sigma_thermal = real(sqrt (  (kErg * thisOctal%temperature(subcell)) / (thisMolecule%molecularWeight * amu) ))
           ! Add turbulent line width in quadrature with thermal line width. vturb is in km/s so 1e5 converts to cm/s
           sigma_thermal = real( sqrt ( sigma_thermal**2 + (vturb*1.0d5)**2 ) )
           ! Convert to Torus units (v/c)
           sigma_thermal = real(sigma_thermal / cspeed)
           dvAcrossCell = abs(dvAcrossCell / sigma_thermal)
        else
           dvAcrossCell = abs(dvAcrossCell * thisOctal%molmicroturb(subcell))
        end if

        if(densitysubsample) then ! should replace 5 with maxdensity/mindensity * fac
           nTau = min(max(dssMinSample, nint(dvAcrossCell * 5.d0)), 100) ! ensure good resolution / 5 chosen as its the magic number!
        else
           nTau = min(max(2, nint(dvAcrossCell * 5.d0)), 100) ! ensure good resolution / 5 chosen as its the magic number!
        endif
        
        OneOvernTauMinusOne = 1.d0/(nTau - 1.d0)
        
        ds = tval * OneOvernTauMinusOne
        
! Calculate column density. Factor of 1.d10 is to convert ds to cm 
! For h21cm case the octal's rho is HI mass density
! For other cases use the abundance relative to H2. 
        if (present(nCol)) then 
           if (h21cm) then 
              nCol = nCol + (thisOctal%rho(subcell) / (thisMolecule%molecularWeight * amu) ) * tval * 1.d10
           else
              nCol = nCol + ( thisOctal%molabundance(subcell) * thisOctal%nH2(subcell) )* tval * 1.d10
           end if
        end if

        if (present(nCol_H2)) then 
           nCol_H2 = nCol_H2 + thisOctal%nH2(subcell) * tval * 1.d10
        end if

! molabundance stores CO abundance relative to H2
        if(present(nCol_CO)) then
           nCol_CO = nCol_CO + thisOctal%molabundance(subcell) * thisOctal%nH2(subcell) * tval * 1.d10
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

           if ( present(tau) ) tau = tau + dtau

           attenuateddIovercell = attenuateddIovercell + dI

           i0 = i0 * exp(-1.0*dtau) + dI

        end do ntau_loop

        ! Change in brightness temperature over this cell
        dIovercell = (i0 - previous_i0) * ( (cSpeed/thisMolecule%transfreq(itrans)) **2) / (2.0 * kErg)

        if ( dIovercell > 0 ) then 
           if (present(i0_pos)) i0_pos = i0_pos + dIovercell
        else
           if (present(i0_neg)) i0_neg = i0_neg + dIovercell
        end if

        n = thisoctal%newmolecularlevel(4,subcell)
                    
        thisoctal%newmolecularlevel(5,subcell) = (n * thisoctal%newmolecularlevel(5,subcell) + dtauovercell) / (n + 1.d0)

        ! average change in brightness temperature (per column density) over this cell
        thisoctal%newmolecularlevel(1,subcell) = (n * thisoctal%newmolecularlevel(1,subcell) + (dIovercell/nCol) ) / (n + 1.d0) 

        ! Image co-ordinates        
        thisoctal%newmolecularlevel(2,subcell) = (n * thisoctal%newmolecularlevel(2,subcell) + this_gal_lon) / (n + 1.d0)
        thisoctal%newmolecularlevel(3,subcell) = (n * thisoctal%newmolecularlevel(3,subcell) + this_gal_lat) / (n + 1.d0)

        thisoctal%newmolecularlevel(4,subcell) = thisoctal%newmolecularlevel(4,subcell) + 1.d0

        if (doWriteLosInfo) call write_los_info

        currentPosition = currentPosition + (tval + origdeps) * otherDirection
        distToObs = (currentPosition - position) .dot. direction
        
     enddo ingrid_loop
         
     if (nColOnly) call CO_vs_H2

666  continue

     contains

       subroutine write_los_info

         implicit none

         real(double) :: this_pos, this_vel, rhoH2, H2frac

        ! Write out los profile information before update to currentPosition
         if ( vel_chan_num == nv .and. abs(this_gal_lat) < plane_lat .and. .not. nColOnly) then 

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

       subroutine CO_vs_H2
         implicit none

         integer :: i_dist, j_dist
         real(double) :: log_nCol_H2, log_nCol_CO
         real(double) :: dx, dy 

! Don't bother with tiny column densities
         if (nCol_H2 < 1000.0 .or. nCol_CO < 1000.0 ) return

! Will be working with log(column density)
         log_nCol_H2 = log10(nCol_H2)
         log_nCol_CO = log10(nCol_CO)

! Calculate bin sizes
         dx = ( maxH2 - minH2 ) / real(CO_H2_dist_size,db)
         dy = ( maxCO - minCO ) / real(CO_H2_dist_size,db)

! Work out which bin this point is in. 
! The int function truncates towards zero so add 1
         i_dist = int( (log_nCol_H2-minH2) / dx) + 1 
         j_dist = int( (log_nCol_CO-minCO) / dy) + 1

! Increment the array of point densities if this point lies within the range considered
         if (i_dist > 0 .and. i_dist <= CO_H2_dist_size .and. &
              j_dist > 0 .and. j_dist <= CO_H2_dist_size) then
!$OMP ATOMIC
            CO_H2_dist(i_dist, j_dist) = CO_H2_dist(i_dist, j_dist) + 1
         endif

       end subroutine CO_vs_H2

   end subroutine intensityAlongRayRev

#ifdef SPH
!-----------------------------------------------------------------------------------------

   subroutine map_dI_to_particles(grid)

    use inputs_mod, only: sphdatafilename, convertRhoToHI
    use sph_data_class, only: sphdata, read_sph_data_wrapper
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

    real(double) :: dI, n_sample, l_coord, b_coord, CO_abund
    real(double) :: distTotorus ! conversion factor between SPH postions and Torus positions
    real(double) :: H2_frac, temperature

#ifdef MPI
    character(len=3)    :: char_my_rank
#endif
    character(len=30)   :: outfilename
    integer, parameter  :: LUIN = 10 ! unit number of output file

    logical :: havesphfile

! Don't do this if the user only wants column densities.
    if (nColOnly) return

! Check that the SPH file exists. If this is a run set up some other way then there won't be a valid SPH dump
    inquire (file=sphdatafilename, exist=havesphfile)
    if (.not. havesphfile) then 
       call writeInfo("sphdatafilename does not specifiy an SPH dump so will assume no mapping to particles required.", forInfo)
       return
    end if

    call writeVtkFile(grid, "ray_info.vtk", valueTypeString=(/"dI      ", "galLon  ", "galLat  ", "crossing"/) )

! Re-read particle data which has been deleted to save memory
    call read_sph_data_wrapper

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
           l_coord  =  thisOctal%newmolecularlevel(2,subcell)
           b_coord  =  thisOctal%newmolecularlevel(3,subcell)
           n_sample =  thisOctal%newmolecularlevel(4,subcell)
           
! Calculate fraction of molecular hydrogen by mass. sphdata%rhon is now the 
! total hydrogen mass density, even if convertRhoToHI is true. The conversion 
! to HI density takes place in ClusterParameter and doesn't affect the value of sphdata%rhon. 
           if ( associated (sphdata%rhoH2) ) then
              H2_frac  = sphdata%rhoH2(ipart) / sphdata%rhon(ipart)
           else
              H2_frac = 0.0
           end if

! Similarly we now calculate the CO abundance. 
           if ( associated (sphdata%rhoCO) ) then
              CO_abund = sphdata%rhoCO(ipart) / sphdata%rhon(ipart)
           else
              CO_abund = 0.0
           end if

! Convert internal energy to temperature
           temperature = sphdata%temperature(ipart) * sphdata%codeEnergytoTemperature

! newmolecularlevel is a floating point number so n_sample>0.99 is a reliable way of saying one or more samples.
           if ( n_sample > 0.99 ) write(LUIN,'(13(e15.8,2x),i8)') old_position, sphdata%gasmass(ipart), sphdata%hn(ipart), &
                sphdata%rhon(ipart), dI, n_sample, H2_frac, temperature, l_coord, b_coord, CO_abund, ipart

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
     write(LUIN,'(a)') "temperature"
     write(LUIN,'(a)') "l"
     write(LUIN,'(a)') "b"
     write(LUIN,'(a)') "CO"
     write(LUIN,'(a)') "index"
     close (LUIN)

   end subroutine map_dI_to_particles
#endif
!-----------------------------------------------------------------------------------------

subroutine write_CO_vs_H2
  implicit none

  integer :: i_dist, j_dist
  real(double) :: dx, dy, x, y

! Calculate bin sizes
  dx = ( maxH2 - minH2 ) / real(CO_H2_dist_size,db)
  dy = ( maxCO - minCO ) / real(CO_H2_dist_size,db)

  open(unit=82, status="replace", file="CO_vs_H2.dat")

  do j_dist = 1, CO_H2_dist_size
     do i_dist = 1, CO_H2_dist_size
        x = minH2 + (real(i_dist,db) - 0.5) * dx
        y = minCO + (real(j_dist,db) - 0.5) * dy
        write(82,*) x, y, CO_H2_dist(i_dist, j_dist) 
     end do
  end do

  close(82)

  deallocate (CO_H2_dist)

end subroutine write_CO_vs_H2

end module angularImage
#endif

