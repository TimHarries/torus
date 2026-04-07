module readathena_mod
    
    use constants_mod
    use utils_mod
    use grid_mod, only: gridtype
    use inputs_mod, only : athenaFilename
    implicit none   
    contains

    subroutine readAthena(grid)
        use amr_mod, only : myScaleSmooth, findTotalMass
        type(gridtype), intent(inout) :: grid
        integer :: nblocks, nr, ntheta
        real(double) :: totalMass
        real(double), allocatable :: rarray(:,:), thetaarray(:,:), rhoarray(:,:,:,:)
        logical :: converged, gridConverged
        

        call writeInfo("Reading Athena++ output...", TRIVIAL)
       

        call read_athena(athenaFilename, nblocks, nr, ntheta, rarray, thetaarray, rhoarray)
       
        ! Here you would typically use the read arrays to initialize your grid structure
        ! For example:
        ! grid%r = rarray
        ! grid%theta = thetaarray
        ! grid%rho = rhoarray
        ! grid%ps = psarray

        converged = .false.
        do while (.not. converged)
            converged = .true.
            call splitGridAthena(grid%octreeRoot, grid, nr, ntheta, nblocks, rarray, thetaarray, rhoarray,  converged)
             do
                   gridConverged = .true.
                   call myScaleSmooth(3., grid, &
                        gridConverged,  inheritProps = .false., &
                        interpProps = .false.)
                   if (gridConverged) exit
                end do
        enddo
             ! Check for convergence criteria here (e.g., based on cell size or density gradients)
             ! If converged, set co

        call fillGridAthena(grid%octreeRoot, nr, ntheta, nblocks, rarray, thetaarray, rhoarray)
        totalMass = 0.d0
        call findTotalMass(grid%octreeRoot, totalMass)
        write(*,*) "Total mass on athena grid",totalmass/msol, " solar masses"

         ! Deallocate arrays after use
        deallocate(rarray)
        deallocate(thetaarray)
        deallocate(rhoarray)

    end subroutine readAthena

    recursive subroutine splitGridAthena(thisOctal, grid, nr, ntheta, nblock, rarray, thetaarray, rhoarray, converged)
        use inputs_mod, only: maxDepthAMR
        use amr_mod, only: addNewChild
        type(OCTAL), pointer :: thisOctal, child
        type(gridtype), intent(inout) :: grid
        real(double), intent(in) :: rarray(:,:), thetaarray(:,:), rhoarray(:,:,:,:)
        integer, intent(in) :: nblock, nr, ntheta
        logical, intent(inout) :: converged
        type(VECTOR) :: rVec
        real(double) :: r, theta, rhocell, cellsize, ang, rcylindrical
        integer :: subcell, i
        logical :: found
        do subcell = 1, thisOctal%maxChildren
            if (thisOctal%hasChild(subcell)) then
         ! find the child
             do i = 1, thisOctal%nChildren, 1
                if (thisOctal%indexChild(i) == subcell) then
               child => thisOctal%child(i)
                  call splitGridAthena(child, grid, nr, ntheta, nblock, rarray, thetaarray, rhoarray, converged)
              exit
           end if
        end do
     else
        rVec = subcellCentre(thisOctal, subcell)
        r = modulus(rVec)*1.e10
        theta = acos(rVec%z/modulus(rvec))
        rcylindrical = sqrt(rVec%x**2 + rVec%y**2)*1.e10

        ang = acos(rcylindrical/r)*180.d0/pi
        if (ang > 20.) cycle
        call return_rho_size(r, theta, nblock, nr, ntheta, rarray, thetaarray, rhoarray, rhocell, cellsize, found)
        if (found) then
            cellsize = cellsize/1.d10
        else
            cellsize = 1.e30
        endif   
        if ((thisOctal%subcellSize > cellsize).and.(thisOctal%nDepth < maxDepthAMR)) then
            call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                         inherit=.true., interp=.false.)
            converged = .false.
            exit
        endif  
    endif   
    enddo
    end subroutine splitGridAthena
        
    recursive subroutine fillGridAthena(thisOctal, nr, ntheta, nblock, rarray, thetaarray, rhoarray)
        type(OCTAL), pointer :: thisOctal, child
        real(double), intent(in) :: rarray(:,:), thetaarray(:,:), rhoarray(:,:,:,:)
        integer, intent(in) :: nblock, nr, ntheta
        real(double) :: r, theta, rhocell, cellsize
        type(VECTOR) :: rVec
        integer :: subcell,i 
        logical :: found
        do subcell = 1, thisOctal%maxChildren
            if (thisOctal%hasChild(subcell)) then
             ! find the child
                do i = 1, thisOctal%nChildren, 1
                   if (thisOctal%indexChild(i) == subcell) then
                        child => thisOctal%child(i)
                         call fillGridAthena(child, nr, ntheta, nblock, rarray, thetaarray, rhoarray)
                         exit
                    end if
                  end do
             else
                rVec = subcellCentre(thisOctal, subcell)
                r = modulus(rVec)*1.e10
                theta = acos(rVec%z/modulus(rvec))
                call return_rho_size(r, theta, nblock, nr, ntheta, rarray, thetaarray, rhoarray, rhocell, cellsize, found)
                if (found) then
                    thisOctal%rho(subcell) = max(rhocell,1.d-25)
                else
                    thisOctal%rho(subcell) = 1.0d-25  ! or some default value
                endif   
            endif
        enddo
    end subroutine fillGridAthena   

    subroutine read_athena(filename, nblocks, nr, ntheta, rarray, thetaarray, rhoarray)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: nblocks, nr, ntheta
        real(double), allocatable, intent(out) :: rarray(:,:), thetaarray(:,:), rhoarray(:,:,:,:)
        real(double), allocatable :: psarray(:,:,:,:)
        integer :: i, j, k, block, nx3, nx2, nx1
        character(len=256) :: header, block_header
        integer :: iunit = 20

        open(unit=iunit, file=filename, status='old', action='read')

        ! Read X1F_ARRAY section (r-grid)
        read(iunit, '(A)') header
        if (header /= 'X1F_ARRAY') then
            print *, 'ERROR: Expected X1F_ARRAY header'
            stop
        end if
        read(iunit, *) nblocks, nr
        allocate(rarray(nblocks, nr))
        do block = 1, nblocks
            do i = 1, nr
                read(iunit, *) rarray(block, i)
            end do
        end do

        ! Read X2F_ARRAY section (theta-grid)
        read(iunit, '(A)') header
        if (header /= 'X2F_ARRAY') then
            print *, 'ERROR: Expected X2F_ARRAY header'
            stop
        end if
        read(iunit, *) nblocks, ntheta
        allocate(thetaarray(nblocks, ntheta))
        do block = 1, nblocks
            do j = 1, ntheta
                read(iunit, *) thetaarray(block, j)
            end do
        end do

        ! Read RHO_BLOCKS section
        read(iunit, '(A)') header
        if (header /= 'RHO_BLOCKS') then
            print *, 'ERROR: Expected RHO_BLOCKS header'
            stop
        end if
        read(iunit, *) nblocks, nx3, nx2, nx1
        allocate(rhoarray(nblocks, nx3, nx2, nx1))
        do block = 1, nblocks
            read(iunit, '(A)') block_header
            do k = 1, nx3
                do j = 1, nx2
                    do i = 1, nx1
                        read(iunit, *) rhoarray(block, k, j, i)
                    end do
                end do
            end do
        end do

        ! Read PS_BLOCKS section (passive scalars)
        read(iunit, '(A)') header
        if (header /= 'PS_BLOCKS') then
            print *, 'ERROR: Expected PS_BLOCKS header'
            stop
        end if
        read(iunit, *) nblocks, nx3, nx2, nx1
        allocate(psarray(nblocks, nx3, nx2, nx1))
        do block = 1, nblocks
            read(iunit, '(A)') block_header
            do k = 1, nx3
                do j = 1, nx2
                    do i = 1, nx1
                        read(iunit, *) psarray(block, k, j, i)
                    end do
                end do
            end do
        end do

        rhoarray = psarray * rhoArray

        close(iunit)

    end subroutine read_athena

    subroutine return_rho_size(r, theta, nblocks, nr, ntheta, rarray, thetaarray, rhoarray, rho_value, cellsize, found)
        real(double), intent(in) :: r, theta
        integer, intent(in) :: nblocks, nr, ntheta
        real(double), intent(in) :: rarray(:,:), thetaarray(:,:), rhoarray(:,:,:,:)
        real(double), intent(out) :: rho_value, cellsize
        integer :: block, k, j, i
        logical :: found
        logical, save :: firstTime = .true.
        real(double) :: dr, dtheta, r_size, theta_size

        rho_value = 0.0d0
        cellsize = 0.0d0
        found = .false.
        k = 1  ! 2D simulation has only 1 phi-cell

        ! Loop over all blocks to find which one contains the point (r, theta)
        do block = 1, nblocks
            if (rarray(block, 2) < r .and. r < rarray(block, nr-2) .and. &
                thetaarray(block, 2) < theta .and. theta < thetaarray(block, ntheta-2)) then
                
                ! Find the index in theta using binary search logic
                do j = 1, ntheta-1
                    if (thetaarray(block, j) <= theta .and. theta < thetaarray(block, j+1)) then
                        exit
                    end if
                end do
                
                ! Find the index in r using binary search logic
                do i = 1, nr-1
                    if (rarray(block, i) <= r .and. r < rarray(block, i+1)) then
                        exit
                    end if
                end do
                
                rho_value = rhoarray(block, k, j, i)
                
                ! Calculate cell sizes
                dr = rarray(block, i+1) - rarray(block, i)
                dtheta = thetaarray(block, j+1) - thetaarray(block, j)
                r_size = dr
                theta_size = r * dtheta  ! Arc length in theta direction
                cellsize = max(r_size, theta_size)
                
                found = .true.
                exit
            end if
        end do

        if (.not. found) then
           if (firstTime) then
              print *, 'ERROR: Point (r, theta) = (', r, ',', theta, ') not found in any block'
              firstTime = .false.
           endif
        end if

    end subroutine return_rho_size


end module
