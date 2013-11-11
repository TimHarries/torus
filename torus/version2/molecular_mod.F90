!!!subroutine testOpticalDepth
      subroutine testOpticalDepth(grid,thisMolecule)

        use inputs_mod, only : lamstart, lamend, rinner, router, debug, useDust
        type(GRIDTYPE) :: grid
        type(MOLECULETYPE) :: thisMolecule
        type(octal), pointer   :: thisOctal
        type(VECTOR) :: currentposition(3), posvec, viewvec, unitvec, centrevec
        integer :: subcell, i, itrans
        integer :: ilamb, nlamb
        real(double) :: xmidplane, gridsize
        real :: lamb
        real(double) :: tau, dummy, kappaAbs, kappaSca, i0
        character(len=50) :: message
        itrans = 0
        centreVec= VECTOR(0.d0, 0.d0, 0.d0)
        kappaAbs = 0.d0; kappaSca = 0.d0

        call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule)
        call addDustToOctalParams(grid, grid%OctreeRoot, thisMolecule)
  
        write(message,*) "Angular dependence"
        call writeinfo(message, FORINFO)
        
        gridsize = grid%octreeroot%subcellsize
        nlamb = 500
        
        do i = 1, 90
           call intensityAlongRay(posvec, viewvec, grid, thisMolecule, iTrans, 0.d0, i0, tau)
     
           write(message,*) i, tau
           call writeinfo(message, FORINFO)
        enddo
  
        write(message,*) "Tau from above"
        call writeinfo(message, FORINFO)
        
        do i = 1, 100
     
           unitvec = VECTOR(1.d-20,1.d-20,1.d0)
           posvec = VECTOR(real(i) * grid%octreeroot%subcellsize * 0.01,0d7,2.d0 * grid%octreeroot%subcellsize)
           viewvec = (-1.d0) * unitvec
           
           call intensityAlongRay(posvec, viewvec, grid, thisMolecule, iTrans, 0.d0, i0, tau)
           
           write(message,*) i, tau
           call writeinfo(message, FORINFO)
        enddo
        
        write(message,*) "Tau from the side"
        call writeinfo(message, FORINFO)
        
        do i = 1, 100
           unitvec = VECTOR(1.d0,0.d0,0.d0)
           posvec = VECTOR(2.d0 * grid%octreeroot%subcellsize,0.d0, real(i) * grid%octreeroot%subcellsize * 0.01)
           viewvec = (-1.d0) * unitvec
           call intensityAlongRay(posvec, viewvec, grid, thisMolecule, iTrans, 0.d0, i0, tau)
           
           write(message,*) i, tau
           call writeinfo(message, FORINFO)
        enddo
        
        call intensityAlongRay(VECTOR(0.d0,0.d0,0.d0),VECTOR(1d-20,1d-20,1.d0), grid, thisMolecule, 1, 0.d0, &
                               dummy, tau, .true.)
        
        xmidplane = rinner + (router-rinner) * 0.5
  
        if(writeoutput) write(*,*) "Midplane tau!"
        if (useDust) then
           open(50, file = 'dustcheck.dat', status="unknown", form = "formatted") 
           currentposition(1) = VECTOR(xmidplane,0.d0,0.d0)
           
           call findSubcellTD(currentPosition(1), grid%octreeroot, thisOctal, subcell)
           
           do i = 0, nlamb
              
              lamb = real((10.0**(dble(i) * log10(lamend/lamstart) / 500.0)) * lamstart)
              
              call continuumIntensityAlongRay(VECTOR(1.d-10,1d-10,1d-10),VECTOR(1.0,1d-20,1.d-20), &
                                              grid, lamb, dummy, tau, .true.)
              call locate(grid%lamArray, size(grid%lamArray), real(lamb), ilamb)
              call returnKappa(grid, thisOctal, subcell, ilambda = ilamb, lambda = real(lamb),&
                   kappaAbs = kappaAbs, kappaSca = kappaSca)
              write(50, '(f8.4,tr3,f10.4,tr3,f10.4,tr3,f10.4)') lamb * 1d-4, tau, &
                   kappaAbs*1e-10 / thisOctal%rho(subcell), (kappaAbs+kappaSca)*1e-10 / thisOctal%rho(subcell)
           enddo
           
           close(50)
        endif
        
        if(debug) then
           
           currentposition(1) =  VECTOR(xmidplane,0.d0,0.d0)
           currentposition(2) =  VECTOR(xmidplane,0.d0,xmidplane)
           currentposition(3) = VECTOR(xmidplane/sqrt(2.),xmidplane/sqrt(2.),2. * xmidplane)
           
           do i = 1,3
              
              call findSubcellTD(currentPosition(i), grid%octreeroot, thisOctal, subcell)
              currentposition(i) = subcellcentre(thisOctal, subcell)
              
              write(message,*) "currentposition", currentposition(i)
              call writeinfo(message, FORINFO)
              write (message,*) "r", sqrt(currentposition(i)%x**2+currentposition(i)%y**2)
              call writeinfo(message, FORINFO)
              write(message,*) "cell velocity",modulus(thisOctal%velocity(subcell)) * cspeed / 1d5
              call writeinfo(message, FORINFO)
              write(message,*) "calc V ",modulus(keplerianVelocity(currentposition(i))) * cspeed / 1d5
              call writeinfo(message, FORINFO)
              
           enddo
   
        endif
        
        call continuumIntensityAlongRay(VECTOR(-1.d10,1.d-10,1.d-10),VECTOR(1.d0,-1d-10,-1d-10), &
                                        grid, 1e4, dummy, tau, .true.)
        write(message,*) "TAU @ 1micron", tau
        call writeinfo(message, FORINFO)
        
      end subroutine testOpticalDepth

      subroutine plotdiscValues(grid, thisMolecule)
        
        type(GRIDTYPE) :: grid
        type(MOLECULETYPE) :: thisMolecule
  real(double) :: mean(6)
  mean = 0.d0
  call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule)
  call addDustToOctalParams(grid, grid%OctreeRoot, thisMolecule)
  
  call findtempdiff(grid, grid%OctreeRoot, thisMolecule, mean, 2)
  
  call readAMRgrid("molecular_tmp.grid",.false.,grid)
  
  call findtempdiff(grid, grid%OctreeRoot, thisMolecule, mean, 0)
  
  call readAMRgrid("molecular_tmp.grid",.false.,grid)
  
  ! Set everything back to the way it was?
  call findtempdiff(grid, grid%OctreeRoot, thisMolecule, mean, 1)

end subroutine plotdiscValues

 real function nextStep(cube, ivplusone) result (step)
   
   type(DATACUBE) :: cube
   integer :: iv,i,j,ivplusone
   real(double), allocatable :: x(:), y(:), y2(:), intensitysum(:)
   real(double) :: xquad(3),yquad(3),xa(3)
   real(double) :: sumx(4),sumxy(4)
   real(double) :: a(3,3),b(3)
   real(double) :: fac
   real(double), save :: alty2
   real(double), save :: oldstep = 4.d0
   !$OMP THREADPRIVATE (alty2, oldstep)

   iv = ivplusone - 1
   
   allocate(x(iv))
   allocate(y(iv))
   allocate(y2(iv))
   allocate(intensitysum(iv))
   
   do i = 1,iv
      intensitysum(i) = SUM(cube%intensity(1:cube%nx,1:cube%ny,i))
   enddo
   
   fac = intensitysum(1)
   x = cube%vAxis(1:iv)
   y =  intensitysum / fac 
   
   call spline(x,y,iv,1.d31,1.d31,y2)
   
   xquad = x(iv-3:iv-1)
   yquad = y2(iv-3:iv-1)
   xa = xquad
   
   do i=1,4
      sumx(i) = sum(xa)
      sumxy(i) = sum(xa * yquad)
      
      xa = xa*xquad
   enddo
   
   do i = 1,3
      do j = 1,3
         if(6-i-j .ne. 0) then 
            a(i,j) = sumx(6-i-j)
         else
            a(i,j) = 3.d0
         endif
      enddo; enddo
      
      b = (/sumxy(2),sumxy(1),sum(yquad)/)
      call luSlv(a,b)

      write(*,*) b(1)*x(iv)**2+b(2)*x(iv)+b(3)
      
      alty2 = b(1)*x(iv)**2+b(2)*x(iv)+b(3)
      
      step = real((oldstep + min(max(25.d0/abs(alty2),0.5_db),4.d0))/2.d0)
      oldstep = step
      
    end function nextStep

! this routine writes a file of intensity, density, cellsize, cellvolume etc

    subroutine dumpIntensityContributions(grid, thisMolecule) 
     use mpi_global_mod, only:  myRankGlobal
     use parallel_mod
      use inputs_mod, only : itrans
      type(MOLECULETYPE) :: thisMolecule
      type(GRIDTYPE) :: grid
      type(VECTOR) :: viewVec, rayStart
      integer :: nx, ny, nv
      real(double) :: x, y, z, v, tau
      real(double) :: vmax, vmin, nuStart, nuEnd, deltaNu
      real(double) :: itot, i0, icheck
      integer :: i, j ,k, iv1, iv2
      logical, save :: firstTime=.true.
      !$OMP THREADPRIVATE(firstTime)

      viewVec = VECTOR(0.d0, 1.d0, 0.d0)
      nx = 2000
      ny = 2000
      nv = 20

      vmin = -2.d0*(1.d5)/cspeed
      vmax =  2.d0*(1.d5)/cspeed
      iv1 = 1
      iv2 = nv
      

#ifdef MPI
    iv1 = (myrankglobal) * (nv / (nThreadsGlobal)) + 1
    iv2 = (myrankglobal+1) * (nv / (nThreadsGlobal))
    if (myrankglobal == (nThreadsGlobal-1)) iv2 = nv
#endif



      nuStart = thisMolecule%transfreq(itrans) * (1.d0+vmin)
      nuEnd = thisMolecule%transfreq(itrans) * (1.d0+vmax)
      deltaNu = nuEnd - nuStart
      itot = 0.d0
      do i = 1, nx
         if (myrankglobal == 0) write(*,*) "i ",i
         do j = 1, ny 
            do k = iv1, iv2
               x = -grid%octreeRoot%subcellSize + (2.d0*grid%octreeRoot%subcellSize)*dble(i-1)/dble(nx-1)
               z = -grid%octreeRoot%subcellSize + (2.d0*grid%octreeRoot%subcellSize)*dble(j-1)/dble(ny-1)
               y = -grid%octreeRoot%subcellSize
               v = vMin + (vMax-vMin)*dble(k-1)/dble(nv-1)

               if(firstTime) then
                  call writeinfo("Filling Octal parameters for first time",TRIVIAL)
                  call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule)
                  firstTime = .false.
               endif

               rayStart = VECTOR(x, y, z)
               call intensityAlongRay(rayStart, viewvec, grid, thisMolecule, iTrans, v, i0, tau)
            enddo
         enddo
      enddo
      itot = 0.d0
      call sumDi(grid%octreeRoot, deltaNu, itot)
      if (myrankGlobal == 0) then
         open(33, file="contribs.dat", status="unknown",form="formatted")
         icheck = 0.d0
         call writeContributions(grid%octreeRoot, deltaNu, itot, icheck)
         write(*,*) "sanity check ",icheck/itot, icheck,itot
         close(33)
      endif
      call torus_mpi_barrier()
    end subroutine dumpIntensityContributions
