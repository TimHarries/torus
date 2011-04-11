module nbody_mod

  use kind_mod
  use constants_mod
  use source_mod
  use mpi_global_mod
  use mpi_amr_mod
  implicit none

contains

  subroutine calculateGasSourceInteraction(source, nSource, grid)
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    type(GRIDTYPE) :: grid
    real(double) :: eps
    integer :: i
#ifdef MPI
    include 'mpif.h'
    integer :: ierr
    integer :: n
    real(double), allocatable :: temp(:), temp2(:)
#endif

    eps = 2.d0 * grid%halfSmallestSubcell

    do i = 1, nSource
       source(i)%force = VECTOR(0.d0, 0.d0, 0.d0)
    enddo

    call writeInfo("Calculating force from gas on sources...",TRIVIAL)
    call recursiveForceFromGas(grid%octreeRoot, source, nSource, eps)

#ifdef MPI
    if (grid%splitOverMPI) then
       n = nSource*3
       allocate(temp(1:n), temp2(1:n))
       temp = 0.d0
       do  i = 1, nSource
          temp((i-1)*3+1) = source(i)%force%x
          temp((i-1)*3+2) = source(i)%force%y
          temp((i-1)*3+3) = source(i)%force%z
       enddo
    endif
    call MPI_ALLREDUCE(temp, temp2, n, MPI_DOUBLE_PRECISION, MPI_SUM, amrCommunicator, ierr)
    do i = 1, nSource
       source(i)%force%x = temp2((i-1)*3+1)
       source(i)%force%y = temp2((i-1)*3+2)
       source(i)%force%z = temp2((i-1)*3+3)
    enddo
    deallocate(temp, temp2)
#endif
    call writeInfo("Done.",TRIVIAL)

    call writeInfo("Calculating source potential for gas", TRIVIAL)
    call zeroSourcepotential(grid%octreeRoot)
    call applySourcePotential(grid%octreeRoot, source, nSource, eps)
    call writeInfo("Done.", TRIVIAL)

    call writeInfo("Calculating source-source forces", TRIVIAL)
    call sourceSourceForces(source, nSource, eps)
    call writeInfo("Done.", TRIVIAL)

  end subroutine calculateGasSourceInteraction


  subroutine updateSourcePositions(source, nSource, dt)
    real(double) :: dt
    type(VECTOR) :: acc
    type(sourcetype) :: source(:)
    integer :: nSource, i

    do i = 1, nSource
       acc = source(i)%force / source(i)%mass
       source(i)%position = source(i)%position + (source(i)%velocity*dt + (0.5d0 * dt**2) * acc)/1.d10
       source(i)%velocity = source(i)%velocity + acc * dt
       if (writeoutput) then
          write(*,*) "source ",i
          write(*,*) "velocity ",source(i)%velocity/1.d5, " km/s"
          write(*,*) "position ",source(i)%position*1.d10/pctocm, " pc"
       endif
    enddo
  end subroutine updateSourcePositions




  recursive subroutine recursiveForceFromGas(thisOctal, source, nSource, eps)
    type(OCTAL), pointer :: thisOctal, child
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    integer :: iSource
    integer :: subcell, i
    real(double) :: r, eps, dv
    Type(VECTOR) :: rVec, cen
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
             rVec = source(isource)%position - cen
             r = modulus(rVec)
             source(iSource)%force = source(iSource)%force - &
                  (bigG*source(isource)%mass*thisOctal%rho(subcell) * dV / (1.d20*(r**2 + eps**2))) * (rVec/r)
          enddo
       endif
    enddo
  end subroutine recursiveForceFromGas

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
                  (bigG*source(isource)%mass/ (1.d10*max(r,eps)))
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
             rHat = rVec/r
             source(i)%force = source(i)%force - (bigG * source(i)%mass*source(j)%mass / (( r**2 + eps**2)*1.d20)) * rHat
          endif
       enddo
    enddo
  end subroutine sourceSourceForces

  subroutine odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqc,grid)
    type(GRIDTYPE) :: grid
    real(double) :: eps, hdid, hnext, hmin
    integer :: nvar
    integer, parameter :: maxstp=10000, nmax=10
    real(double), parameter :: two=2.d0,zero=0.d0,tiny=1.d-30
    real(double) :: x, x1, h, x2, h1, xsav, dxsav
    integer :: nok, nbad, kount, i, nstp, kmax
    real(double) :: xp(200), yp(10,200)
    common /path/ kmax,kount,dxsav,xp ,yp
    real(double) :: ystart(nvar),yscal(nmax),y(nmax),dydx(nmax)
    external derivs
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
       call derivs(x,y,dydx,grid)
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
       call rkqc(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
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
       if(abs(hnext).lt.hmin) pause 'stepsize smaller than minimum.'
       h=hnext
    enddo
    pause 'too many steps.'
  end subroutine odeint

  subroutine bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs,grid)
    type(GRIDTYPE) :: grid
    integer :: nv
    integer, parameter :: nmax=10,imax=11,nuse=7
    real(double) :: h, htry, x, xsav, xest, errmax
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
       call mmid(ysav,dysav,nv,xsav,h,nseq(i),yseq,derivs,grid)
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
    if(x+h.eq.x)pause 'step size underflow.'
    goto 1
  end subroutine bsstep

  subroutine mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs,grid)
    type(GRIDTYPE)::grid
    real(double) :: xs, htot, h, h2, x, swap
    integer :: nvar
    integer :: nstep, i, n
    integer, parameter :: nmax=10
    real(double) :: y(nvar),dydx(nvar),yout(nvar),ym(nmax),yn(nmax)
    external derivs
    h=htot/nstep
    do i=1,nvar
       ym(i)=y(i)
       yn(i)=y(i)+h*dydx(i)
    enddo
    x=xs+h
    call derivs(x,yn,yout,grid)
    h2=2.*h
    do n=2,nstep
       do i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
       enddo
       x=x+h
       call derivs(x,yn,yout,grid)
    enddo
    do i=1,nvar
       yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
    enddo
  end subroutine mmid

  subroutine derivs(x, y, dydx, grid)
    type(GRIDTYPE) :: grid
    real(double) :: x, y(:), dydx(:)
    integer :: i, ia
    type(SOURCETYPE), allocatable :: testsource(:)

    ! y(1,2,3,4,5,6) is first source
    ! y(1) = x position
    ! y(2) = x velocity
    ! y(3) = y position
    ! y(4) = y velocity
    ! y(5) = z position
    ! y(6) = z velocity

    allocate(testSource(1:globalnSource))
    do i = 1, globalNSource
       testsource(i) = globalSourceArray(i)
       ia = 6*(i-1)+1
       testSource(i)%position = VECTOR(y(ia), y(ia+2), y(ia+4))
       testSource(i)%force = VECTOR(0.d0, 0.d0, 0.d0)
    enddo

    call calculateGasSourceInteraction(testsource, globalnSource, grid)

    do i = 1 , globalNSource*6-1,2
       dydx(i) = y(i+1) ! dx/dt = v_x etc
    enddo

    do i = 1 , globalNSource
       ia = (i-1)*6+1
       dydx(ia+1) = globalSourceArray(i)%force%x / testsource(i)%mass
       dydx(ia+3) = globalSourceArray(i)%force%y / testsource(i)%mass
       dydx(ia+5) = globalSourceArray(i)%force%z / testsource(i)%mass
    enddo
    deallocate(testSource)
  end subroutine derivs

  SUBROUTINE RZEXTR(IEST,XEST,YEST,YZ,DY,NV,NUSE)
    integer, parameter :: IMAX=11,NMAX=10,NCOL=7
    integer :: iest, nuse, j, m1, k , nv
    real(double) :: xest, yy, v, c, b1,b , ddy
    real(double) ::  X(IMAX),YEST(NV),YZ(NV),DY(NV),D(NMAX,NCOL),FX(NCOL)
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
