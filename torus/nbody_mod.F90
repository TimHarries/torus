#ifdef HYDRO
module nbody_mod

  use kind_mod
  use dimensionality_mod
  use constants_mod
  use source_mod
  use mpi_global_mod
  use mpi_amr_mod
  implicit none

contains


  subroutine calculateGasSourceInteraction(source, nSource, grid, eps)
#ifdef MPI
    use mpi
#endif
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    type(GRIDTYPE) :: grid
    real(double) :: eps
    integer :: i
#ifdef MPI
    integer :: ierr
    integer :: n
    real(double), allocatable :: temp(:), temp2(:)
#endif


    do i = 1, nSource
       source(i)%force = VECTOR(0.d0, 0.d0, 0.d0)
    enddo

!    call writeInfo("Calculating force from gas on sources...",TRIVIAL)
    if (grid%splitOverMpi) call recursiveForceFromGas(grid%octreeRoot, source, nSource, eps)

#ifdef MPI
    if (grid%splitOverMPI.and.(myrankGlobal/=0)) then
       n = nSource*3
       allocate(temp(1:n), temp2(1:n))
       temp = 0.d0
       do  i = 1, nSource
          temp((i-1)*3+1) = source(i)%force%x
          temp((i-1)*3+2) = source(i)%force%y
          temp((i-1)*3+3) = source(i)%force%z
       enddo
       call MPI_ALLREDUCE(temp, temp2, n, MPI_DOUBLE_PRECISION, MPI_SUM, amrCommunicator, ierr)
       do i = 1, nSource
          source(i)%force%x = temp2((i-1)*3+1)
          source(i)%force%y = temp2((i-1)*3+2)
          source(i)%force%z = temp2((i-1)*3+3)
       enddo
       deallocate(temp, temp2)
    endif
#endif
!    call writeInfo("Done.",TRIVIAL)

!    call writeInfo("Calculating source potential for gas", TRIVIAL)
!    call writeInfo("Done.", TRIVIAL)

!    call writeInfo("Calculating source-source forces", TRIVIAL)
    call sourceSourceForces(source, nSource, eps)
!    call writeInfo("Done.", TRIVIAL)

  end subroutine calculateGasSourceInteraction

  subroutine donBodyOnly(tEnd, dt, grid)
    use vtk_mod, only : writevtkfilenbody
    use source_mod, only : globalnSource, globalSourceArray
    character(len=80)  :: plotFile
    type(GRIDTYPE) :: grid
    real(double) :: dt, currentTime, tEnd
    real(double) :: totalEnergy, ePot, eKin
    real(double), allocatable :: tmp(:)
    integer :: it,i 

    currentTime = 0.d0
    it = 0
    if (grid%splitOverMpi.and.(myrankGlobal == 0)) goto 666

    if (writeoutput) open(44, file="energy.dat", form="formatted", status="unknown")
    if (writeoutput) write(44,'(a)') "Step      Time (s)      p.e. (J)        k.e. (J)      E_total (J)"
    if (writeoutput) open(45, file="position.dat", form="formatted", status="unknown")
    if (Writeoutput) write(45,'(a)') "Step      Time (s)     position x y z..." 
    do while (currentTime  < tEnd)
       write(plotfile,'(a,i4.4,a)') "nbody",it,".vtk"
       call writeVtkFilenBody(globalnSource, globalsourceArray, plotfile, grid)
       call sumEnergy(globalsourcearray, globalnSource, totalenergy, ePot, eKin, grid)
       if (writeoutput) write(44,'(i6,1p,20e15.5,0p)') it, currentTime, ePot, eKin, totalEnergy
       allocate(tmp(1:globalnSource*3))
       do i = 1, globalnSource
          tmp(i*3-2) = globalSourceArray(i)%position%x
          tmp(i*3-1) = globalSourceArray(i)%position%y
          tmp(i*3) = globalSourceArray(i)%position%z
       enddo
       if (writeoutput) write(45,'(i6,1p,20e15.5,0p)') it, currentTime, tmp(1:globalnsource*3)
       deallocate(tmp)
       call  nBodyStep(globalsourceArray, globalnSource, dt, grid)
       currentTime = currentTime + dt
       if (writeoutput) write(*,*) "Current time ",currentTime
       it = it + 1
    end do
    if (writeoutput) close(44)
    if (Writeoutput) close(45)
    666 continue
  end subroutine donBodyOnly
       

    

  subroutine updateSourcePositions(source, nSource, dt, grid)
    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(sourcetype) :: source(:)
    integer :: nSource, i, ia, nvar
    real(double), allocatable :: yStart(:), dydx(:)
    integer :: nok, nbad
    integer :: kmax, kount
    real(double) :: acc, eps, minDt, thisDt, thisTime
    real(double) :: dxsav, xp(2000), yp(2000,2000)
    common /path/ kmax,kount,dxsav,xp ,yp
    kmax = 0

    call calculateEps(source, nsource,eps)

    nvar = nSource * 6
    allocate(yStart(1:nvar), dydx(1:nVar))

    do i = 1,nSource
       ia = (i-1)*6 + 1
       ystart(ia+0) = returnCodeUnitLength(source(i)%position%x*1.d10)
       ystart(ia+1) = returnCodeUnitSpeed(source(i)%velocity%x)

       ystart(ia+2) = returnCodeUnitLength(source(i)%position%y*1.d10)
       ystart(ia+3) = returnCodeUnitSpeed(source(i)%velocity%y)

       ystart(ia+4) = returnCodeUnitLength(source(i)%position%z*1.d10)
       ystart(ia+5) = returnCodeUnitSpeed(source(i)%velocity%z)


    enddo

    thisTime = 0.d0
    do while (thisTime  < dt)

       call derivs(0.d0, ystart, dydx, grid, eps)

       minDt = 1.d30
       do i = 1, nSource
          ia = (i-1)*6 + 1
          acc = sqrt(dydx(ia+1)**2 + dydx(ia+3)**2 + dydx(ia+5)**2)
          minDt = min(sqrt(eps / (acc + tiny(Acc))), minDt)
       enddo
       thisDt = min(minDt, dt)
       if ((thisTime + thisDt) > dt) then
          thisDt = dt - thisTime
       endif
       if (myrankglobal == 1) write(*,*) "calling integrator with ",thisDt, thisTime, dt
       call odeint(ystart, nvar, 0.d0, thisDt, 1.d-3, thisDt, 0.d0, nok, nbad, derivs, bsstep, grid, eps)
       thisTime = thisTime + thisDt
       do i = 1, nSource
          ia = (i-1)*6 + 1
          source(i)%position%x = returnPhysicalUnitLength(ystart(ia+0))/1.d10
          source(i)%velocity%x = returnPhysicalUnitSpeed(ystart(ia+1))
          
          source(i)%position%y = returnPhysicalUnitLength(ystart(ia+2))/1.d10
          source(i)%velocity%y = returnPhysicalUnitSpeed(ystart(ia+3))
          
          source(i)%position%z = returnPhysicalUnitLength(ystart(ia+4))/1.d10
          source(i)%velocity%z = returnPhysicalUnitSpeed(ystart(ia+5))
    
       enddo
       if (writeoutput) then
          open(57, file="pos.dat",status="old",position="append")
          write(57,*)  real(source(1)%position%x*1.d10/autocm), &
               real(source(1)%position%y*1.d10/autocm), &
               real(source(2)%position%x*1.d10/autocm), &
               real(source(2)%position%y*1.d10/autocm), &
               real(source(3)%position%x*1.d10/autocm), &
               real(source(3)%position%y*1.d10/autocm)
          close(57)
       endif
    enddo
    deallocate(yStart,dydx)
!    call sumEnergy(source, nSource, energy)
!    if (Writeoutput) write(*,*) "Total energy ",energy, nok, nbad
  end subroutine updateSourcePositions

  subroutine calculateEps(source, nsource, eps)
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    real(double) :: eps, r, minsep
    integer :: i, j

    minSep = 1.d30
    do i = 1, nSource
       do j = 1, nSource
          if (i /= j) then
             r = modulus(source(i)%position-source(j)%position)*1.d10
             minSep = MIN(minSep, r)
          endif
       enddo
    enddo
    eps = 1.d-6 * minSep
    eps = max(eps, 1000.d0*autocm)

  end subroutine calculateEps

  subroutine nBodyStep(source, nSource, dt, grid)
    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(sourcetype) :: source(:)
    integer :: nSource, i, ia, nvar
    real(double), allocatable :: yStart(:), dydx(:)
    integer :: nok, nbad
    integer :: kmax, kount
    real(double) :: acc, eps, minDt, thisDt, thisTime
    real(double) :: dxsav, xp(2000), yp(2000,2000)
    common /path/ kmax,kount,dxsav,xp ,yp
    kmax = 0

    call calculateEps(source, nSource, eps)

    nvar = nSource * 6
    allocate(yStart(1:nvar), dydx(1:nVar))

    do i = 1,nSource
       ia = (i-1)*6 + 1
       ystart(ia+0) = returnCodeUnitLength(source(i)%position%x*1.d10)
       ystart(ia+1) = returnCodeUnitSpeed(source(i)%velocity%x)

       ystart(ia+2) = returnCodeUnitLength(source(i)%position%y*1.d10)
       ystart(ia+3) = returnCodeUnitSpeed(source(i)%velocity%y)

       ystart(ia+4) = returnCodeUnitLength(source(i)%position%z*1.d10)
       ystart(ia+5) = returnCodeUnitSpeed(source(i)%velocity%z)


    enddo

    thisTime = 0.d0
    do while (thisTime  < dt)

       call derivs(0.d0, ystart, dydx, grid, eps)

       minDt = 1.d30
       do i = 1, nSource
          ia = (i-1)*6 + 1
          acc = sqrt(dydx(ia+1)**2 + dydx(ia+3)**2 + dydx(ia+5)**2)
          minDt = min(sqrt(eps / (acc + tiny(Acc))), minDt)
       enddo
       thisDt = min(minDt, dt)
       if ((thisTime + thisDt) > dt) then
          thisDt = dt - thisTime
       endif
       if (writeoutput) write(*,*) "Calling ode integrator at ",100.d0*thisTime/dt, " % of total step"
       call odeint(ystart, nvar, 0.d0, thisDt, 1.d-3, thisDt, 0.d0, nok, nbad, derivs, bsstep, grid, eps)
       thisTime = thisTime + thisDt
       do i = 1, nSource
          ia = (i-1)*6 + 1
          source(i)%position%x = returnPhysicalUnitLength(ystart(ia+0))/1.d10
          source(i)%velocity%x = returnPhysicalUnitSpeed(ystart(ia+1))
          
          source(i)%position%y = returnPhysicalUnitLength(ystart(ia+2))/1.d10
          source(i)%velocity%y = returnPhysicalUnitSpeed(ystart(ia+3))
          
          source(i)%position%z = returnPhysicalUnitLength(ystart(ia+4))/1.d10
          source(i)%velocity%z = returnPhysicalUnitSpeed(ystart(ia+5))
          
       enddo
    enddo
    deallocate(yStart,dydx)
  end subroutine nBodyStep




  recursive subroutine recursiveForceFromGas(thisOctal, source, nSource, eps)
    use amr_utils_mod, only : inSubcell
    type(OCTAL), pointer :: thisOctal, child
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    integer :: iSource
    integer :: subcell, i
    real(double) :: r, eps, dv
    real(double) :: massSub
    integer :: i1, j1, k1
    integer, parameter  :: nSub = 8
    Type(VECTOR) :: rVec, cen, pos
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call recursiveForceFromGas(child, source, nSource, eps)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
          cen = subcellCentre(thisOctal, subcell)
          dv = cellVolume(thisOctal, subcell)*1.d30
          do iSource = 1, nSource
             if (.not.inSubcell(thisOctal, subcell, source(isource)%position)) then
                rVec = source(isource)%position - cen
                r = modulus(rVec)
                source(iSource)%force = source(iSource)%force - &
                     (bigG*source(isource)%mass*thisOctal%rho(subcell) * dV / (1.d20*r**2 + eps**2)) * (rVec/r)
                
             else
                massSub  = (thisOctal%rho(subcell)*dv)/dble(nSub**3)
                do i1 = 1, nSub
                   do j1 = 1, nSub
                      do k1 = 1, nSub
                         pos%x = cen%x + (dble(i1-1)+0.5d0) * thisOctal%subcellSize / dble(nSub) - thisOctal%subcellSize/2.d0
                         pos%y = cen%y + (dble(j1-1)+0.5d0) * thisOctal%subcellSize / dble(nSub) - thisOctal%subcellSize/2.d0
                         pos%z = cen%z + (dble(k1-1)+0.5d0) * thisOctal%subcellSize / dble(nSub) - thisOctal%subcellSize/2.d0
                         rVec = source(isource)%position - pos
                         r = modulus(rVec)
                         source(iSource)%force = source(iSource)%force - &
                              (bigG*source(isource)%mass*massSub / (1.d20*r**2 + eps**2)) * (rVec/r)
                      enddo
                   enddo
                enddo

             endif

          enddo
       endif
    enddo
  end subroutine recursiveForceFromGas


  recursive subroutine recursivePotentialFromGas(thisOctal, source, nSource, energy)
    use amr_utils_mod, only : inSubcell
    type(OCTAL), pointer :: thisOctal, child
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    integer :: iSource
    integer :: subcell, i
    real(double) :: r, energy, dv
    real(double) :: massSub
    integer :: i1, j1, k1
    integer, parameter  :: nSub = 8
    Type(VECTOR) :: rVec, cen, pos
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call recursivePotentialFromGas(child, source, nSource, energy)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
          cen = subcellCentre(thisOctal, subcell)
          dv = cellVolume(thisOctal, subcell)*1.d30
          do iSource = 1, nSource
             if (.not.inSubcell(thisOctal, subcell, source(isource)%position)) then
                rVec = source(isource)%position - cen
                r = modulus(rVec)
                energy = energy - &
                     (bigG*source(isource)%mass*thisOctal%rho(subcell) * dV / (1.d10*r))
                
             else
                massSub  = (thisOctal%rho(subcell)*dv)/dble(nSub**3)
                do i1 = 1, nSub
                   do j1 = 1, nSub
                      do k1 = 1, nSub
                         pos%x = cen%x + (dble(i1-1)+0.5d0) * thisOctal%subcellSize / dble(nSub) - thisOctal%subcellSize/2.d0
                         pos%y = cen%y + (dble(j1-1)+0.5d0) * thisOctal%subcellSize / dble(nSub) - thisOctal%subcellSize/2.d0
                         pos%z = cen%z + (dble(k1-1)+0.5d0) * thisOctal%subcellSize / dble(nSub) - thisOctal%subcellSize/2.d0
                         rVec = source(isource)%position - pos
                         r = modulus(rVec)
                         energy = energy - &
                              (bigG*source(isource)%mass*massSub / (1.d10*r))
                      enddo
                   enddo
                enddo

             endif

          enddo
       endif
    enddo
  end subroutine recursivePotentialFromGas


  recursive subroutine applySourcePotential(thisOctal, source, nSource, eps)
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    type(OCTAL), pointer :: thisOctal, child
    integer :: subcell, i, iSource
    type(VECTOR) :: rVec, cen
    real(double) :: r, eps
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call  applySourcePotential(child, source, nSource, eps)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
          cen = subcellCentre(thisOctal, subcell)
          do iSource = 1, nSource
             rVec = cen - source(isource)%position
             r = modulus(rVec)
             thisOctal%phi_stars(subcell) = thisOctal%phi_stars(subcell) - &
                  (bigG*source(isource)%mass/ (max(r*1.d10,eps)))
          enddo
       endif
    enddo
  end subroutine applySourcePotential

  recursive subroutine zeroSourcePotential(thisOctal)
    type(OCTAL), pointer :: thisOctal, child
    integer :: subcell, i
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call  zeroSourcePotential(child)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
          thisOctal%phi_stars(subcell) = 0.d0
       endif
    enddo
  end subroutine zeroSourcePotential


  subroutine sourceSourceForces(source, nSource, eps)
    type(SOURCETYPE) :: source(:)
    real(double) :: eps, r
    integer :: nSource
    type(VECTOR) :: rHat, rVec
    integer :: i, j

    do i = 1, nSource
       do j = 1, nSource
          if (i /= j) then
             rVec = source(j)%position - source(i)%position
             r = modulus(rVec)
             if (r /= 0.d0) then
                rHat = rVec/r
                source(i)%force = source(i)%force + (bigG * source(i)%mass*source(j)%mass / ( 1.d20*r**2 + eps**2)) * rHat
!                if ((r/eps < 10.0).and.(writeoutput)) write(*,*) "Gravity softening becoming important"
             endif
          endif
       enddo
    enddo
  end subroutine sourceSourceForces

  subroutine sumEnergy(source, nSource, totalenergy, ePot, eKin, grid)
#ifdef MPI
    use mpi
    integer :: ierr
#endif
    type(GRIDTYPE) :: grid
    type(SOURCETYPE) :: source(:)
    integer :: nSource, i, j
    real(double) :: totalenergy, epot, ekin, temp, ePotGas
    logical :: split

    split = grid%splitOverMpi
    temp = 0.d0
    totalenergy = 0.d0
    epot = 0.d0
    ekin = 0.d0
    do i = 1, nSource
       do j = 1, nSource

          if (i /= j) then
             epot = epot - 0.5d0 * bigG*source(i)%mass*source(j)%mass / &
                  (modulus(source(i)%position-source(j)%position)*1.d10)
          endif
       enddo
       ekin = ekin + 0.5d0 * source(i)%mass * modulus(source(i)%velocity)**2
    enddo

    ePotGas = 0.d0
#ifdef MPI
    if ((grid%splitOverMPI).and.(myrankGlobal/=0)) then
       call recursivePotentialFromGas(grid%octreeRoot, source, nSource, ePotGas)
       call MPI_ALLREDUCE(epotGas, temp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, amrCommunicator, ierr)
       epot = epot + temp
    endif
#endif

   totalenergy = epot + ekin
  end subroutine sumEnergy

  subroutine odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqc,grid,graveps)
    type(GRIDTYPE) :: grid
    real(double) :: eps, hdid, hnext, hmin, graveps
    integer :: nvar
    integer, parameter :: maxstp=10000, nmax=10000
    real(double), parameter :: two=2.d0,zero=0.d0,tiny=1.d-30
    real(double) :: x, x1, h, x2, h1, xsav, dxsav
    integer :: nok, nbad, kount, i, nstp, kmax
    real(double) :: xp(2000), yp(2000,2000)
    common /path/ kmax,kount,dxsav,xp ,yp
    real(double) :: ystart(nvar),yscal(nmax),y(nmax),dydx(nmax)
    external derivs, rkqc
    x=x1
    h=sign(h1,x2-x1)
    nok=0
    nbad=0
    kount=0
    do  i=1,nvar
       y(i)=ystart(i)
    enddo
    xsav=x-dxsav*two
    do nstp=1,maxstp
       call derivs(x,y,dydx,grid, graveps)
       do  i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+tiny
       enddo
       if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
             if(kount.lt.kmax-1)then
                kount=kount+1
                xp(kount)=x
                do i=1,nvar
                   yp(i,kount)=y(i)
                enddo
                xsav=x
             endif
          endif
       endif
       if((x+h-x2)*(x+h-x1).gt.zero) h=x2-x
       call rkqc(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs,grid, graveps)
       if(hdid.eq.h)then
          nok=nok+1
       else
          nbad=nbad+1
       endif
       if((x-x2)*(x2-x1).ge.zero)then
          do i=1,nvar
             ystart(i)=y(i)
          enddo
          if(kmax.ne.0)then
             kount=kount+1
             xp(kount)=x
             do i=1,nvar
                yp(i,kount)=y(i)
             enddo
          endif
          return
       endif
       if(abs(hnext).lt.hmin) then
          write(*,*) 'stepsize smaller than minimum.'
          stop
       endif

       h=hnext
    enddo
    write(*,*) 'too many steps.'
    stop
  end subroutine odeint

  subroutine bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs,grid, graveps)
    type(GRIDTYPE) :: grid
    integer :: nv
    integer, parameter :: nmax=10000,imax=11,nuse=7
    real(double) :: h, htry, x, xsav, xest, errmax, graveps
    integer :: i, j
    real(double) :: eps, hdid,  hnext
    real(double), parameter :: one=1.d0,shrink=.95d0,grow=1.2d0
    real(double) ::  y(nv),dydx(nv),yscal(nv),yerr(nmax), &
         ysav(nmax),dysav(nmax),yseq(nmax)
    integer :: nseq(imax)
    external derivs
    nseq = (/2,4,6,8,12,16,24,32,48,64,96/)
    h=htry
    xsav=x
    do i=1,nv
       ysav(i)=y(i)
       dysav(i)=dydx(i)
    enddo
1   continue
    do i=1,imax
       call mmid(ysav,dysav,nv,xsav,h,nseq(i),yseq,derivs,grid, graveps)
       xest=(h/nseq(i))**2
       call rzextr(i,xest,yseq,y,yerr,nv,nuse)
       errmax=0.
       do j=1,nv
          errmax=max(errmax,abs(yerr(j)/yscal(j)))
       enddo
       errmax=errmax/eps
       if(errmax.lt.one) then
          x=x+h
          hdid=h
          if(i.eq.nuse)then
             hnext=h*shrink
          else if(i.eq.nuse-1)then
             hnext=h*grow
          else
             hnext=(h*nseq(nuse-1))/nseq(i)
          endif
          return
       endif
    enddo
    h=0.25*h/2**((imax-nuse)/2)
    if(x+h.eq.x) then
       write(*,*) 'step size underflow.'
       stop
    endif
    goto 1
  end subroutine bsstep

  subroutine mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs,grid, graveps)
    type(GRIDTYPE)::grid
    real(double) :: xs, htot, h, h2, x, swap, graveps
    integer :: nvar
    integer :: nstep, i, n
    integer, parameter :: nmax=10000
    real(double) :: y(nvar),dydx(nvar),yout(nvar),ym(nmax),yn(nmax)
    external derivs
    h=htot/nstep
    do i=1,nvar
       ym(i)=y(i)
       yn(i)=y(i)+h*dydx(i)
    enddo
    x=xs+h
    call derivs(x,yn,yout,grid, graveps)
    h2=2.*h
    do n=2,nstep
       do i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
       enddo
       x=x+h
       call derivs(x,yn,yout,grid, graveps)
    enddo
    do i=1,nvar
       yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
    enddo
  end subroutine mmid

  subroutine derivs(x, y, dydx, grid, graveps)
    integer, parameter :: nmax = 10000
    type(GRIDTYPE) :: grid
    real(double) :: x, y(nmax), dydx(nmax), test, graveps
    integer :: i, ia
    type(SOURCETYPE), allocatable :: testsource(:)

    test = x
    ! y(1,2,3,4,5,6) is first source
    ! y(1) = x position
    ! y(2) = x velocity
    ! y(3) = y position
    ! y(4) = y velocity
    ! y(5) = z position
    ! y(6) = z velocity

    allocate(testSource(1:globalnSource))
    do i = 1, globalNSource
       testsource(i)%mass = globalSourceArray(i)%mass
       ia = 6*(i-1)+1
       testSource(i)%position = (VECTOR(returnPhysicalUnitLength(y(ia)), &
            returnPhysicalUnitLength(y(ia+2)), &
            returnPhysicalUnitLength(y(ia+4))))/1.d10
       testSource(i)%force = VECTOR(0.d0, 0.d0, 0.d0)
    enddo

    call calculateGasSourceInteraction(testsource, globalnSource, grid, graveps)

    do i = 1 , globalNSource*6-1,2
       dydx(i) = y(i+1) ! dx/dt = v_x etc
    enddo

    do i = 1 , globalNSource
       ia = (i-1)*6+1
       dydx(ia+1) = returnCodeUnitAcceleration(testsource(i)%force%x / testsource(i)%mass)
       dydx(ia+3) = returnCodeUnitAcceleration(testsource(i)%force%y / testsource(i)%mass)
       dydx(ia+5) = returnCodeUnitAcceleration(testsource(i)%force%z / testsource(i)%mass)
    enddo


!    if (writeoutput) then
!       write(*,*) "x ",x
!       do i = 1, globalnSource
!          ia = (i-1)*6+1
!          write(*,*) "source ",i
!          write(*,*) "position ",y(ia+0),y(ia+2),y(ia+4)
!          write(*,*) "velocity ",y(ia+1),y(ia+3),y(ia+5)
!       enddo
!    endif
    deallocate(testSource)
  end subroutine derivs

  SUBROUTINE RZEXTR(IEST,XEST,YEST,YZ,DY,NV,NUSE)
    integer, parameter :: IMAX=11,NMAX=10000,NCOL=7
    integer :: iest, nuse, j, m1, k , nv
    real(double) :: xest, yy, v, c, b1,b , ddy
    real(double) ::  YEST(NV),YZ(NV),DY(NV)
    real(double), save :: x(imax), D(NMAX,NCOL),FX(NCOL)
    X(IEST)=XEST
    IF(IEST.EQ.1) THEN
       DO  J=1,NV
          YZ(J)=YEST(J)
          D(J,1)=YEST(J)
          DY(J)=YEST(J)
       enddo
    ELSE
       M1=MIN(IEST,NUSE)
       DO K=1,M1-1
          FX(K+1)=X(IEST-K)/XEST
       enddo
       DO J=1,NV
          YY=YEST(J)
          V=D(J,1)
          C=YY
          D(J,1)=YY
          DO K=2,M1
             B1=FX(K)*V
             B=B1-C
             IF(B.NE.0.) THEN
                B=(C-V)/B
                DDY=C*B
                C=B1*B
             ELSE
                DDY=V
             ENDIF
             V=D(J,K)
             D(J,K)=DDY
             YY=YY+DDY
          enddo
          DY(J)=DDY
          YZ(J)=YY
       enddo
    ENDIF
  END SUBROUTINE RZEXTR



end module nbody_mod
#endif
