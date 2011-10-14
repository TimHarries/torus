! This module deals with blobs or clumps which can be
! used to distort the opacity grid. It includes 
! subroutines to initialize the blobs, move them, and
! to distort the grid

! written by tjh

! v1.0 on 13/08/99


module blob_mod

  use vector_mod            ! Vector maths
  use random_mod
  use constants_mod         ! Physical constants
  use gridtype_mod, only: GRIDTYPE  ! The density grid

  implicit none

  public

  ! The type definition

  type BLOBTYPE
     type(VECTOR) :: position          ! Position vector
     type(VECTOR) :: velocity          ! Velocity
     real :: contrast                  ! Density contrast
     real :: radius                    ! Blob `radius'
     logical :: inUse                  ! Is this blob in use?
  end type BLOBTYPE


contains

  ! subroutine to initialize an individual blob

  subroutine initBlob(blobs, grid, iBlob, atBase, blobContrast)
    use utils_mod, only: locate

    integer :: iBlob               ! The blob to init
    logical :: atBase              ! Start the blob at the base?
    type(BLOBTYPE) :: blobs(*)     ! array of blobs
    type(GRIDTYPE) :: grid         ! the opacity grid
    integer :: i1,i2,i3            ! indices
    real :: r1, r2, r3             ! random deviates
    real(double) :: r, mu, phi, theta      ! spherical polar coords
    real :: sinTheta              
    real :: blobContrast
 
    ! set the blob as in use

    blobs(iBlob)%inUse = .true.

    ! put the blob at the base if requested, otherwise
    ! put it between 1 and 10 stellar radii

    if (atBase) then
       r = grid%rAxis(1)*1.05
    else
       call randomNumberGenerator(getreal=r1)
       r = (1.05 + 9. * sqrt(r1)) *grid%rAxis(1)
    endif

    ! random latitude

    call randomNumberGenerator(getreal=r2)
    mu = 2.*r2 -1.
    sinTheta = real(sqrt(1.- mu*mu))

    ! random azimuth

    call randomNumberGenerator(getreal=r3)
    phi = r3 * twoPi

    ! set the position

    blobs(iBlob)%position = VECTOR(r*sinTheta*cos(phi),r*sinTheta*sin(phi),r*mu)

    ! radius random between 0.01 and 0.1 stellar radii

    call randomNumberGenerator(getreal=r1)
    blobs(iBlob)%radius = 0.1*grid%raxis(1) !(0.01 + 0.09*r1)*grid%rAxis(1)

    ! find the position of the blob within the opacity grid


    call getPolar(blobs(iBlob)%position, r, theta, phi)
    mu = cos(theta)
    call locate(grid%rAxis, grid%nr, real(r), i1)
    call locate(grid%muAxis, grid%nmu, real(mu), i2)
    call locate(grid%phiAxis, grid%nphi, real(phi), i3)

    ! set up the instantaneous blob velocity
    ! NB the velocity grid is actually unitless doppler shifts and
    ! the distances are all in 10^10 cm.

    blobs(iBlob)%velocity = (cSpeed/1.e10)*grid%velocity(i1,i2,i3)

    ! blob contrast is hardwired
    
    blobs(iBlob)%contrast = blobContrast

  end subroutine initBlob

  ! this subroutine moves the blobs about the grid between
  ! timeStart and timeEnd

  subroutine moveBlobs(maxBlobs, blobs, timeStart, timeEnd, grid)
    USE math_mod, only: interpGridVelocity
    USE grid_mod, only: getIndices

    real :: timeStart, timeEnd           ! the times
    integer :: maxBlobs                  ! the number of blobs in the array
    type(BLOBTYPE) :: blobs(*)           ! the blob array
    type(GRIDTYPE) :: grid               ! the opacity grid          
    real(double) :: t1, t2, t3                   
    integer :: i, j, i1, i2, i3          ! counters
    integer ::   nTimes = 100000    ! the number of time segments
    real :: thisdTime, timeScale, vel

    ! loop over all blobs

    do j = 1, maxBlobs

       ! is this blob active?

       if (blobs(j)%inUse) then
          
          call getIndices(grid, blobs(j)%position, i1, i2, i3, t1, t2, t3)
          vel = real(modulus(grid%velocity(i1,1,1))*cSpeed/1.e10)
          if (i1 /= grid%nr) then
             timeScale = (grid%rAxis(i1+1)-grid%rAxis(i1))/vel
          else
             timeScale = (grid%rAxis(i1)-grid%rAxis(i1-1))/vel
          endif
          
          nTimes = max(1,int((timeEnd-timeStart)/(timeScale/10.)))
          thisdTime = (timeEnd-timeStart)/real(nTimes)
         

          ! loop over all times

          do i = 1, nTimes

             ! update the blob position

             blobs(j)%position = blobs(j)%position + dble(thisdTime) * blobs(j)%velocity

             ! is the blob outside 20 core radii? If so, switch it off

             if ((modulus(blobs(j)%position)/grid%rAxis(1)) > 10.) then
                blobs(j)%inUse = .false.
                EXIT

             else

                ! if not find its new location in the grid

                call getIndices(grid, blobs(j)%position, i1, i2, i3, t1, t2, t3)

                ! update the velocity

                blobs(j)%velocity = (cSpeed/1.d10) * &
                     interpGridvelocity(grid, i1,i2,i3, t1, t2, t3)
             endif


          enddo   ! over times

       endif
    enddo   ! over all blobs

  end subroutine moveBlobs

  ! this subroutine distorts the grid by the blobs

  subroutine distortGridWithBlobs(grid, maxBlobs, blobs)

    type(GRIDTYPE) :: grid                ! the opacity grid
    integer :: maxBlobs                   ! max no of blobs
    type(BLOBTYPE) :: blobs(:)            ! blob array 
    real :: distance, sinTheta
    integer :: i, j, k, m
    real :: fac1,fac2,r
    real, allocatable :: facGrid(:,:,:)   ! enhancement factor grid
    type(VECTOR), allocatable :: posGrid(:,:,:)   ! enhancement factor grid

    ! this may take a while for a large grid/number of blobs

    write(*,*) "Distorting grid by blobs..."
 
    ! set up a temporary array 

    allocate(facGrid(1:grid%nr, 1:grid%nMu, 1:grid%nPhi))
    allocate(posGrid(1:grid%nr, 1:grid%nMu, 1:grid%nPhi))

    ! it'll be mostly ones

    facGrid = 1.

    do i = 1, grid%nr
       do j = 1, grid%nMu
          do k = 1, grid%nPhi
             
             ! determine position vector
                   
             sinTheta = sqrt(1.-grid%muAxis(j)**2)
             posGrid(i,j,k) = VECTOR(grid%rAxis(i)*cos(grid%phiAxis(k))*sinTheta, &
                  grid%rAxis(i)*sin(grid%phiAxis(k))*sinTheta, &
                  grid%rAxis(i)*grid%muAxis(j))
          enddo
       enddo
    enddo
                   

    
    ! loop over all the blobs

    do m = 1, maxBlobs


       ! only blobs in use can distort grid
       
       
       if (blobs(m)%inUse) then
!          write(*,*) "blob: ",m
          ! loop over grid
          
          do i = 1, grid%nr
             do j = 1, grid%nMu
                do k = 1, grid%nPhi
             
                   distance = real(modulus(posGrid(i,j,k) - blobs(m)%position))
                   fac1 = distance**2 / (2.*blobs(m)%radius**2)

                   if (fac1 < 10.) then
                      fac2 = blobs(m)%contrast

                      ! distort by a 3D gaussian
                   
                      
                      facGrid(i,j,k) = min(facgrid(i,j,k) + fac2 * exp(-fac1),blobs(m)%contrast)
                      r = real(modulus(posGrid(i,j,k)/dble(grid%rAxis(1))))
!                      if ((fac1 < 1.) .and. (r > 1.) .and. (r < 3.)) then
!                         write(*,*) r,facGrid(i,j,k)
!                      endif
                   endif


                enddo
             enddo
          enddo

       endif

    enddo

 ! should be able to use intrinsic arrays to do these loops,
 ! but doesn't seem to work at the moment, despite the fact
 ! the arrays are conformable

!    do i = 1, grid%nr
!       do j = 1, grid%nMu
!          do k = 1, grid%nPhi
!             if (facGrid(i,j,k) < 1.01) facGrid(i,j,k) = 1.e-20
!          enddo
!       enddo
!    enddo
             
 

    do i = 1, grid%nr
       do j = 1, grid%nMu
          do k = 1, grid%nPhi
             
             ! line opacities go as f^2
             
             grid%chiLine(i,j,k) = grid%chiLine(i,j,k)*facGrid(i,j,k)**2
             grid%etaLine(i,j,k) = grid%etaLine(i,j,k)*facGrid(i,j,k)**2

             

             ! scattering opacities go like f

             grid%kappaSca(i,j,k,1) = &
                  grid%kappaSca(i,j,k,1)*facGrid(i,j,k)

          enddo
       enddo
    enddo
    


    ! tell the users we've done

    write(*,*) "Max distortion: ",MAXVAL(facGrid)
    write(*,*) "Min distortion: ",MINVAL(facGrid)
 
    deallocate(facGrid) ! free up the array
    deallocate(posGrid) ! free up the array
 
    write(*,*) "Done."
  end subroutine distortGridWithBlobs


  ! add a new blob according to poissonion stats


  subroutine addNewBlobs(grid, maxBlobs, blobs, blobTime, dTime, nCurr, blobContrast)
    USE utils_mod, only: poidev

    integer :: maxBlobs               ! max no of blobs
    type(BLOBTYPE) :: blobs(maxBlobs) ! blob array
    integer :: nNewBlobs              ! no of new blobs
    type(GRIDTYPE) :: grid            ! opacity grid
    real :: blobTime                  ! the timescale for blob creation
    real :: dTime                     ! time interval
    logical :: freeBlob               ! found a free blob
    integer :: j,k, iBlob             ! counters
    integer :: nCurr                  ! current number of blobs
    real :: blobContrast
    

    ! the number of new blobs is found according to the random
    ! poissonion deviate.

    nNewBlobs = int(poidev(dTime/blobTime))

    nCurr = 0
    do iBlob = 1, maxBlobs
       if (blobs(iBlob)%inUse) then
          nCurr = nCurr + 1
       endif
    enddo


    ! loop over the number of new blobs

    do k = 1 , nNewBlobs

       ! find a slot for the new blob

       freeBlob = .false.
       do iBlob = 1, maxBlobs
          if (.not.blobs(iBlob)%inUse) then
             freeBlob = .true.
             j = iBlob
          endif
       enddo

       ! if a free slot is found then initialize it otherwise
       ! warn the user - this means that maxblobs must be 
       ! increased

       if (freeBlob) then
          call initBlob(blobs, grid, j, .true., blobContrast)
       else
          write(*,'(A)') "! No slot free for new blob"
          write(*,'(A)') "! increase maxBlobs"
       endif

    enddo  ! over new blobs





  end subroutine addNewBlobs


  ! subroutine to write current blob configuration to a file

  subroutine writeBlobs(filename, maxBlobs, blobs)


    character(len=*) :: filename        ! the filename
    type(BLOBTYPE) :: blobs(*)          ! the blob array
    integer :: maxBlobs                 ! max no of blobs
    integer :: i                        ! counter

    ! open file as new to prevent overwriting previous blob
    ! configurations - an easy mistake to make

    open(40, file=filename, status='new', form='unformatted', err=666)

    ! first line contains number of blobs

    write(40) maxBlobs

    ! loop over all blobs writing out details

    do i = 1, maxBlobs

       write(40) blobs(i)%position
       write(40) blobs(i)%velocity
       write(40) blobs(i)%contrast
       write(40) blobs(i)%radius
       write(40) blobs(i)% inUse

    end do ! over all blobs

    close(40) ! close output

    goto 1000 ! quit

    ! cannot open the blob file

666 continue
    write(*,'(a,a,a)') "! Blob configuration file ",trim(filename)," already exists"
    stop


1000 continue


  end subroutine writeBlobs


  ! read in the blob configuration from a file

  subroutine readBlobs(filename, maxBlobs, blobs, quiet)

    character(len=*) :: filename     ! the filename
    type(BLOBTYPE) :: blobs(:)       ! blob array
    integer :: maxBlobs              ! max number of blobs
    integer :: i                     ! counter
    integer :: nBlobs                ! actual no of blobs
    logical :: quiet

    ! info for user

    if (.not.quiet) then
       write(*,'(a,a)') "Reading blob configuration from: ",trim(filename)
    endif

    ! file must already exist

    open(40, file=filename, status='old', form='unformatted')

    read(40) nBlobs
    
    ! this might happen

    if (nBlobs > maxBlobs) then
       write(*,'(a)') "! Too many blobs in configuration file"
       write(*,'(a)') "! increase maxBlobs"
    endif

    ! read in configuration

    do i = 1, nBlobs

       read(40) blobs(i)%position
       read(40) blobs(i)%velocity
       read(40) blobs(i)%contrast
       read(40) blobs(i)%radius
       read(40) blobs(i)% inUse
    end do

    close(40) ! close file
  end subroutine readBlobs

  subroutine writeBlobPgp(viewVec, t1, t2, nphase, maxBlobs)
    integer :: i,j
    real :: t1,t2
    real :: thisTime
    character(len=80) :: blobfile, outfile
    integer :: maxBlobs
    integer :: nPhase
    type(VECTOR) :: viewVec
    type(BLOBTYPE), allocatable :: blobs(:)
    real :: oldVel
    real :: projVel
    logical :: firstTime

    allocate(blobs(1:maxBlobs))
    

    write(outfile,'(a)') "blobs.pgp"
    open(20,file=outfile,form='formatted',status='unknown')
    do j=maxBlobs,maxBlobs-100,-1
       firsttime=.true.
       oldVel = 0.
       write(*,*) j
       do i = 1, nPhase
          write(blobfile,'(a,i3.3,a)') "run",i,".blob"
          call readBlobs(blobfile, maxBlobs, blobs, .true.)
          thisTime = t1 + real(i-1)*(t2-t1)/real(nPhase-1)
          projVel = real(blobs(j)%velocity .dot. viewVec)
          if (blobs(j)%inUse) then
             if (firsttime) then
                write(20,*) "move",projVel*1.e5,thisTime
                oldVel = abs(projVel)
                firsttime = .false.
             else
                if (abs(projVel) > oldVel) then
                   if (projVel < 0.) then
                      write(20,*) "sci 4"
                   else
                      write(20,*) "sci 2"
                   endif
                   write(20,*) "draw",projVel*1.e5, thisTime
                   oldVel = abs(projVel)
                else
                   write(20,*) "move",projVel*1.e5, thisTime
                   oldVel = abs(projVel)
                endif
             endif
          endif
       enddo
    enddo
    close(20)
    deallocate(blobs)
  end subroutine writeBlobPgp

      
  
  


end module blob_mod



