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

  contains

    subroutine make_angular_image(grid)

      use input_variables, only : lineImage, maxRhoCalc, densitySubSample
      use datacube_mod, only: DATACUBE, writeDataCube
      use gridtype_mod, only: GRIDTYPE
      use molecular_mod, only:  moleculetype

      implicit none

      TYPE(gridtype), intent(in) :: grid
      type(DATACUBE) ::  cube
      type(MOLECULETYPE) :: thisMolecule
      integer :: itrans, nSubpixels
      real, parameter :: thisWavelength=21.0

      lineImage        = .true.
      maxRhoCalc       = .false. 
      densitySubSample = .false.

! molecular weight is used for column density calculation
      thisMolecule%molecularWeight = mHydrogen / amu


      ! Set up 21cm line
      allocate( thisMolecule%transfreq(1) )
      thisMolecule%transfreq(1) = cSpeed / (21.0)

      call writeinfo('Generating internal view', TRIVIAL)
      call createAngImage(cube, grid, thisMolecule)
      cube%intensity(:,:,:) = cube%intensity(:,:,:) * (thisWavelength**2) / (2.0 * kErg)
      if(writeoutput) call writedatacube(cube, "theGalaxy.fits")

    end subroutine make_angular_image

!-----------------------------------------------------------------------------------------------------------

    subroutine createAngImage(cube, grid, thisMolecule)

      use input_variables, only : gridDistance, beamsize, npixels, nv, imageside, maxVel, usedust, nsubpixels
      use datacube_mod, only: telescope, initCube, addSpatialAxes, addvelocityAxis
      use molecular_mod, only: calculateOctalParams, moleculetype, intensitytoflux
      use atom_mod, only: bnu
      use vector_mod

      implicit none

     type(TELESCOPE) :: mytelescope
     type(MOLECULETYPE) :: thisMolecule
     type(GRIDTYPE) :: grid
     type(DATACUBE) :: cube
     type(VECTOR) :: viewvec, observervec, imagebasis(2)
     real(double) :: minVel
     real(double) :: deltaV
     integer :: iTrans = 1
     integer :: i
     integer :: iv
     character(len=200) :: message
     real(double) :: intensitysum, fluxsum, ddv!, dummy
     real(double), save :: background
     real(double), allocatable :: weightedfluxmap(:,:)
     real(double) :: weightedfluxsum, weightedflux
     real(double), allocatable :: fineweightedfluxmap(:,:)
     real(double) :: fineweightedfluxsum, fineweightedflux
     real(double), allocatable :: weight(:,:)
     real, allocatable :: temp(:,:,:) 

 ! SHOULD HAVE CASE(TELESCOPE) STATEMENT HERE

     mytelescope%label = 'JCMT'
     mytelescope%diameter = 15.d2 ! diameter in cm
     mytelescope%beamsize = beamsize

     minVel = (-1.d0) * maxVel

     call writeinfo("Initialising datacube",TRIVIAL)

     if(nv .eq. 0) then
        call initCube(cube, npixels, npixels, 200, mytelescope) ! Make cube
     else
        call initCube(cube, npixels, npixels, nv, mytelescope) ! Make cube
     endif

     cube%obsDistance = gridDistance * 1d10!(in cm) Additional information that will be useful
     write(message,'(a,f7.2,a)') "Observer Distance        : ",gridDistance/pctocm, " pc"
     call writeinfo(message, TRIVIAL) 
     write(message,'(a,f7.2,a)') "Finest grid resolution   : ",grid%halfsmallestsubcell*2d10/autocm, " AU"
     call writeinfo(message, TRIVIAL) 

     if(nv .ne. 0) then 
        call addvelocityAxis(cube, minVel, maxVel) ! velocities in km/s from +ve (redder, away) to -ve (bluer,towards)
     else
        cube%vAxis(1) = minVel
     endif

     deltaV = minVel * 1.e5/cspeed_sgl
     
     allocate(temp(npixels,npixels,3))

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

!        if(usedust) call adddusttoOctalParams(grid, grid%OctreeRoot, thisMolecule, deltaV)

        temp(:,:,:) = 0.d0

        call makeAngImageGrid(grid, cube, thisMolecule, itrans, deltaV, nSubpixels, temp)

        cube%intensity(:,:,iv) = real(temp(:,:,1))
        cube%tau(:,:,iv)       = real(temp(:,:,2))
        cube%nCol(:,:)         = real(temp(:,:,3)) 

        if(writeoutput) then
           call tune(6, message)  ! stop a stopwatch
        endif

        intensitysum = sum(temp(:,:,1)) / dble(npixels**2)
        fluxsum = intensitytoflux(intensitysum, dble(imageside), dble(gridDistance), thisMolecule)

        if(iv .eq. 1) then
           background = Bnu(thisMolecule%transfreq(itrans), Tcbr)
           write(message, *) "Background Intensity: ",background
           call writeinfo(message, TRIVIAL)
           background = intensitytoflux(background, dble(imageside), dble(gridDistance), thisMolecule)
           write(message, *) "Background Flux: ",background
           call writeinfo(message, TRIVIAL)
           write(message, *) ""
           call writeinfo(message, TRIVIAL)
              
        endif

        write(message,'(a,es11.4e1,tr3,a,f10.4,tr3,a,es12.4,a,es12.4,es12.4,es12.4)') &
             "DELTAV(v/c):",deltaV," V (km/s):",real(cube%vAxis(iv)), "Average Intensity:",intensitysum, &
             " FLUX: ", fluxsum, (fluxsum / thisMolecule%transfreq(itrans)) * 1e26, (fluxsum - background) &
             / thisMolecule%transfreq(itrans) * 1d26  
        open(10, file="tempfile.dat",status="unknown",form="formatted",position="append")
        write(10,'(es11.4e1,f8.4,es12.4,es12.4,es12.4,es12.4)') &
             real(cube%vAxis(iv)), deltaV, intensitysum, &
             fluxsum / thisMolecule%transfreq(itrans)* 1e26, (fluxsum - background) / thisMolecule%transfreq(itrans) * 1e26 
        close(10)
        call writeinfo(message,FORINFO)
     end do

   end subroutine createAngImage

!-----------------------------------------------------------------------------------------------------------


    subroutine makeAngImageGrid(grid, cube, thisMolecule, itrans, deltaV, nSubpixels, imagegrid)

      use input_variables, only : npixels, imageside, centrevecx, centrevecy
      use vector_mod
      
      type(GRIDTYPE), intent(IN) :: grid
      type(datacube), intent(IN) :: cube
      type(MOLECULETYPE), intent(IN) :: thisMolecule
      integer, intent(IN) :: itrans
      real(double), intent(IN) :: deltaV
      integer, intent(IN) :: nsubpixels
      type(VECTOR) :: viewvec, ObserverVec
      real(double) :: viewvec_x, viewvec_y, viewvec_z
      real, intent(OUT) :: imagegrid(:,:,:)

      real(double) :: dnpixels ! npixels as a double, save conversion
      type(VECTOR) :: pixelcorner
      integer :: subpixels
      integer :: ipixels, jpixels
      integer :: index(2)
      real(double) :: pixelside

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
         theta_axis(jpixels) = theta_min + ( real(jpixels) * delta_theta )
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
         do ipixels = 1, npixels

            viewvec_x = sin( theta_axis(jpixels) ) * cos( phi_axis(ipixels) ) 
            viewvec_y = sin( theta_axis(jpixels) ) * sin( phi_axis(ipixels) ) 
            viewvec_z = cos( theta_axis(jpixels) )

            viewvec = VECTOR( viewvec_x, viewvec_y, viewvec_z )
            call normalize(viewvec)

            imagegrid(ipixels,jpixels,:) = &
                 AngPixelIntensity(cube,viewvec,grid,thisMolecule,iTrans,deltaV, subpixels, delta_theta, delta_phi)

         enddo

      enddo

    end subroutine makeAngImageGrid

 !!! Calculates the intensity for a square pixel of arbitrary size, position, orientation
 function AngPixelIntensity(cube,viewvec,grid,thisMolecule,iTrans,deltaV,subpixels,delta_theta, delta_phi) &
      result(out)
   
   use input_variables, only : tolerance, nsubpixels
   use molecular_mod, only: intensityAlongRay, sobseq

   type(GRIDTYPE) :: grid
   type(MOLECULETYPE) :: thisMolecule
   integer :: itrans
   type(VECTOR), intent(in) :: viewVec
   real(double) :: i0, opticaldepth
   real(double), intent(in) :: delta_theta, delta_phi

   type(VECTOR) :: imagebasis(2), pixelbasis(2), pixelcorner
   type(VECTOR), parameter :: rayposition=VECTOR(0.0, 2.1e12_db, 0.0)
   type(VECTOR) :: thisViewVec

   integer :: subpixels, minrays
   integer :: i, iray

   real(double) :: avgIntensityNew, avgIntensityOld
   real(double) :: varIntensityNew, varIntensityOld
   real(double) :: avgNColNew, avgNColOld
   real(double) :: rtemp(2)
   real(double), save ::  r(10000,2)
   real(double) :: nCol

   logical :: converged
   real(double) :: deltaV
   type(DATACUBE) :: cube
   logical, save :: firsttime = .true.

   real(double) :: out(3) 

   if(firsttime) then
      
      call sobseq(rtemp,-1)
      
      do i = 1, 10000
         call sobseq(rtemp)
         r(i,:) = dble(rtemp)
      enddo
      
      firsttime = .false.
   endif

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

      ! Generate rotated viewvec to sample angular cell
      thisViewVec = viewVec
!      thisViewvec = rotateZ(thisViewVec, (r(iray,1) * delta_theta  * degtorad) )
!      thisViewvec = rotateX(thisViewVec, (r(iray,2) * delta_phi    * degtorad) )
!      call normalize(thisViewVec)

     call intensityalongray(rayposition,thisViewVec,grid,thisMolecule,itrans,deltaV,i0,tau=opticaldepth, nCol=nCol)
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
      elseif(iray .gt. 10000) then
         converged = .false.
         out(1) = avgIntensityNew
         out(2) = opticaldepth ! not averaged, just last value
         out(3) = avgNColNew
         exit
      else
         avgIntensityOld = avgIntensityNew
         varIntensityOld = varIntensityNew
         iray = iray + 1
         out(1) = avgIntensityNew
         out(2) = opticaldepth ! not averaged, just last value
         out(3) = avgNColNew
      endif

   enddo
   
 end function AngPixelIntensity

end module angularImage
