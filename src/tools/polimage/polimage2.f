      program polimage2
      implicit none
c
c     This routine is almost identical to polimage except for 
c     the ability of specfying the plot range of the map.
c
      
      INCLUDE '/star/include/sae_par'
      INCLUDE '/star/include/dat_par'
      integer iptr, qptr, uptr, vptr, velptr, aptr1, aptr2, gptr, sptr
      integer saptr
      integer nx, ny, nelm, nxg, nyg
      integer ndf1, ndf2, ndf3, ndf4, ndf5, ndfg, ndfs
      integer maxx, maxy
      integer dim(2), ndim
      logical first, dopol
      logical thisgif,oneframe
      integer status, nphase, i, ilen, ilen2
      real phase
      integer nspec
      logical circular, greyscale
      parameter( maxx = 500, maxy = 500)
      real median
      logical hardcopy, dolog

      character*80 filename, device, gdevice, imageName, specroot
      CHARACTER*(DAT__SZLOC) PLOC,DLOC,QLOC,LOC, LOCg, LOCS


      hardcopy = .false.

      WRITE(*,*) "image name?"
      READ(*,'(A)') imageName
      WRITE(*,*) "image number?"
      READ(*,*) NPHASE
      WRITE(*,*) "Circular polarization?"
      READ(*,*) CIRCULAR
      WRITE(*,*) "Logarithmic scaling?"
      READ(*,*) dolog
      WRITE(*,*) "Grayscale?"
      READ(*,*) greyscale
      if (greyscale) then
         write(*,*) "Spectrum file root?"
         read(*,'(a)') specroot
      endif

      write(*,*) "polarization images?"
      read(*,*) dopol

      write(*,*) "one frame?"
      read(*,*) oneframe

      WRITE(*,*) "device?"
      READ(*,'(a)') gdevice
      STATUS = 0
      first = .true.
      thisgif = .false.
      
      ilen = index(imageName, " ")-1
      DO I = 1, NPHASE
         phase = real(i-1)/real(nPhase-1)
      CALL NDF_BEGIN
      CALL HDS_START(STATUS)
      write(filename,'(a,i3.3)') imageName(1:ilen),i
      CALL HDS_OPEN(filename,'READ',LOC,STATUS)
      CALL NDF_IMPRT(LOC,NDF1,STATUS)
      CALL NDF_DIM(NDF1, 2, DIM, NDIM, STATUS)
      CALL NDF_MAP(NDF1,'DATA','_REAL','READ',IPTR,NELM,STATUS)
      CALL NDF_AMAP(NDF1,'Centre',1,'_REAL','READ',APTR1,NELM,STATUS)
      CALL NDF_AMAP(NDF1,'Centre',2,'_REAL','READ',APTR2,NELM,STATUS)
      CALL NDF_XLOC(NDF1,'POLARIMETRY','READ',PLOC,STATUS)
      CALL NDF_FIND(PLOC,'STOKES_Q',NDF2,STATUS)
      CALL NDF_MAP(NDF2,'DATA','_REAL','READ',QPTR,NELM,STATUS)
      CALL NDF_FIND(PLOC,'STOKES_U',NDF3,STATUS)
      CALL NDF_MAP(NDF3,'DATA','_REAL','READ',UPTR,NELM,STATUS)
      CALL NDF_FIND(PLOC,'STOKES_V',NDF4,STATUS)
      CALL NDF_MAP(NDF4,'DATA','_REAL','READ',VPTR,NELM,STATUS)
      CALL NDF_FIND(PLOC,'VEL',NDF5,STATUS)
      CALL NDF_MAP(NDF5,'DATA','_REAL','READ',VELPTR,NELM,STATUS)

      nx = dim(1)
      ny = dim(2)

      IF (GREYSCALE) THEN
         CALL HDS_OPEN("greyscale",'READ', LOCG, STATUS)
         CALL NDF_IMPRT(LOCG, NDFG, STATUS)
         CALL NDF_DIM(NDFG, 2, DIM, NDIM, STATUS)
         CALL NDF_MAP(NDFG, 'DATA', '_REAL', 'READ', 
     &                GPTR, NELM, STATUS)
         NXG = DIM(1)
         NYG = DIM(2)
         ilen2 = index(specroot," ") - 1
         write(filename,'(a,i3.3)') specroot(1:ilen2),i
         CALL HDS_OPEN(filename,'READ', LOCS, STATUS)
         CALL NDF_IMPRT(LOCS, NDFS, STATUS)
         CALL NDF_MAP(NDFS, 'DATA', '_REAL', 'READ', 
     &                SPTR, NSPEC, STATUS)
         CALL NDF_AMAP(NDFS,'Centre',1,'_REAL','READ',SAPTR,
     &                 NELM,STATUS)

      ENDIF


      if (gdevice(1:3) .eq. "/xs") then
         device = "/xs"
      endif

      if (gdevice(1:4) .eq. "/gif") then
         write(device,'(a,i3.3,a)') imagename(1:ilen),i,".gif/gif"
         thisGif = .true.
      endif
      if (gdevice(1:4) .eq. "/vps") then
         write(device,'(a,i3.3,a)') imagename(1:ilen),i,".ps/cps"
         thisGif = .true.
         hardcopy = .true.
      endif
      write(*,*) "device: ",device(1:20)
      call plotimage(%val(iptr), %val(qptr), %val(uptr), %val(vptr),
     &    %val(velptr),%val(aptr1), %val(aptr2), nx, ny,device,first,
     &    circular, thisgif, greyscale, %val(gptr), nxg, nyg, phase,
     &    %val(sptr), %val(saptr), nspec,hardcopy,specroot,oneframe,
     &    dolog, dopol)
      CALL DAT_ANNUL(PLOC,STATUS)
      CALL NDF_END(STATUS)
      CALL HDS_CLOSE(LOC,STATUS)
      CALL HDS_STOP(STATUS)
      first = .false.
      ENDDO

      end

      subroutine plotimage(iimage, qimage, uimage, vimage,velimage,
     &                     ax1, ax2, nx, ny, 
     &                     device,first,circular, thisgif, 
     &                     dogreyscale, 
     &                     greyscale, nxg, nyg, phase,
     &                     spec, aspec, nspec,hardcopy, specroot,
     &                     oneframe, dolog, dopol)
      implicit none 
      real fg1,bg1
      integer ncols
      integer nspec
      real spec(nspec), aspec(nspec)
      logical first, dogreyscale, dopol
      real fg2,bg2
      real fg3,bg3
      real phase
      integer iclo, ichi
      real wstart, wend
      integer nxg, nyg
      real greyscale(nxg, nyg)
      logical dolog
      logical circular, thisgif, oneframe
      integer nx, ny, i, j
      real fgg,bgg, range
      real iimage(nx,ny)
      real qimage(nx,ny)
      real uimage(nx,ny)
      real vimage(nx,ny)
      real velimage(nx,ny)
      real imax,imin,ir
      real ax1(nx)
      real ax2(ny)
      real tr(6)
      real tr2(6)
      real trg(6)
      real mean,sig
      real c(10)
      integer nmean
      integer maxx, maxy
      real median
      parameter( maxx = 500, maxy = 500)
      real medarray(maxx*maxy)
      real temp(maxx, maxy)
      real xArray(maxx, maxy)
      real yArray(maxx, maxy)
      real limage(maxx, maxy)
      logical hardcopy
      real junk
      real it, qt, ut
      real polmax
      real pol, pa, dt
      integer cindex(4,2)
      integer m,n,ix,iy,nx2,ny2, ci1, ci2
      character*80 device, specroot
      real cut_low, cut_high, roof
      save fg1, bg1, fg2, bg2, fg3, bg3, imax, imin


      dt = (ax1(nx) - ax1(1))/real(nx-1)
      tr(1) = ax1(1) - dt
      tr(2) = dt
      tr(3) = 0.
      dt = (ax2(ny) - ax2(1))/real(ny-1)
      tr(4) = ax2(1) - dt
      tr(5) = 0.
      tr(6) = dt

      if (first) then
         fg1 = -1.e10
         bg1 = 1.e10
      endif
      nmean = 0
      mean = 0.
      do i = 1 , nx
         do j = 1, ny
            if (first) then
               if (iimage(i,j) .gt. fg1) fg1 = iimage(i,j)
               if ((iimage(i,j) .lt. bg1).and.
     &             (iimage(i,j).gt.1e-20)) bg1 = iimage(i,j)
            endif
            if (iimage(i,j) .gt. 1.e-10) then
               mean = mean + iimage(i,j)
               nmean = nmean + 1
               medarray(nmean) = iimage(i,j)
            endif
         enddo
      enddo
c      call sort(nmean, medarray)
c      median = medarray(nmean / 2)

      polmax = -1.e10
      ix = nx/21
      iy = ny/21
      nx2 = 1
      do i = 1 , nx, ix
         ny2 = 1
         do j = 1, ny, iy
            it = 0.
            qt = 0.
            ut = 0.
            do m = i, min(i+ix-1,nx)
               do n = j, min(j+iy-1,ny)
                  it = it + iimage(m,n)
                  qt = qt + qimage(m,n)
                  ut = ut + uimage(m,n)
                enddo
            enddo
            pol = 100.*sqrt(qt*qt + ut*ut) / it
            pa = 0.5*atan2(ut,qt)
            if (it < 1.e-3) pol = 0.
            if (pa .lt. 0.) pa = pa + 3.141592654
            xArray(nx2,ny2) = pol * sin(pa)
            yArray(nx2,ny2) = -pol * cos(pa)
            if (pol .gt. polmax) polmax = pol
            ny2 = ny2 + 1
         enddo
         nx2 = nx2 + 1
      enddo
      nx2 = nx2 -1
      ny2 = ny2 - 1
      write(*,*) nx2,ny2

c      bg = bg * 1.1
c      bg = 0.1
c      fg = fg / 10.

c      fg = median * 2.
c      bg = median * 0.5
      


      write(*,*) fg1, bg1

      if (oneframe) then
         call pgbegin(0,device,1,1)
      else
         if ((.not.dogreyscale).or.(.not.dopol)) then
            call pgbegin(0,device,2,2)
         else
            call pgbegin(0, device, 3, 2)
         endif
      endif

      call pgqcir(iclo, ichi)

      ncols = ichi - iclo + 1
      ncols = ncols / 4
      cindex(1,1) = iclo
      cindex(1,2) = iclo  + ncols - 1
      cindex(2,1) = cindex(1,2) + 1
      cindex(2,2) = cindex(2,1) + nCols-1
      cindex(3,1) = cindex(2,2) + 1
      cindex(3,2) = cindex(3,1) + nCols-1
      cindex(4,1) = cindex(3,2) + 1
      cindex(4,2) = cindex(4,1) + nCols-1


      call pgqcol(ci1, ci2)

      write(*,*) "colour index range",ci1,ci2
      if (thisGif) then
         if (.not.hardcopy) then
            if (dogreyscale.and.dopol) then
               call pgpaper(5.,0.6666666)
            else
               call pgpaper(5.,1.)
            endif
         endif
      endif
      call pgpage
      call pgsch(1.3)

      if (oneframe) then

         write(*,*) "The lower cut (percentile) of mapping value?"
         read(*,*) cut_low 
         write(*,*) "The higher cut (percentile) of mapping value?"
         read(*,*) cut_high 

c     Now adjusts the plotting range

         roof = fg1
         fg1 = roof*cut_low/100.0
         bg1 = roof*cut_high/100.0

         call pgvport(0.15,0.9,0.15,0.9)
         if (dolog) then
            call pgsitf(1)
         else
            call pgsitf(0)
         endif
               
         call pgscf(2)
         call pgscir(cindex(1,1), cindex(1,2))
         call palette(3)
         call pgwnad(ax1(1), ax1(nx), ax2(1), ax2(ny))
         call pgimag(iimage, nx, ny, 1, nx, 1, ny, bg1, fg1, tr)
         call pgbox('bcnst',0,0,'bcnst',0,0)
         call pglabel("Distance (10\\u10\\d cm)", 
     &                "Distance (10\\u10\\d cm)", " ")
         call mywedge(bg1,fg1)
         goto 666
      endif


      call pgvport(0.1,0.9,0.1,0.9)
      if (dolog) then
         call pgsitf(1)
      else
         call pgsitf(0)
      endif
      call pgscir(cindex(1,1), cindex(1,2))
      call palette(2)
      call pgwnad(ax1(1), ax1(nx), ax2(1), ax2(ny))
      call pgimag(iimage, nx, ny, 1, nx, 1, ny, bg1, fg1, tr)
      call pgbox('bcnst',0,0,'bcnst',0,0)
      call pglabel(" ", " ", "Intensity image")
c      call pgwedg("R", 2., 2., fg, bg, " ")
      call mywedge(bg1,fg1)
      call pgpage


      if (dopol) then
      call pgvport(0.1,0.9,0.1,0.9)

      if (dolog) then
         call pgsitf(1)
      else
         call pgsitf(0)
      endif


      call pgscir(cindex(2,1), cindex(2,2))
      call palette(2)
      call pgwnad(ax1(1), ax1(nx), ax2(1), ax2(ny))
      call pgimag(iimage, nx, ny, 1, nx, 1, ny, bg1, fg1, tr)

      dt = (ax1(nx) - ax1(1))/(real(nx2-1))
      tr2(1) = ax1(1) - dt/2.
      tr2(2) = dt
      tr2(3) = 0.
      dt = (ax2(ny) - ax2(1))/(real(ny2-1))
      tr2(4) = ax2(1) - dt/2.
      tr2(5) = 0.
      tr2(6) = dt
      call pgsch(0.)
      call pgsci(2)
      call pgslw(2)
      write(*,*) "polmax",polmax
      call pgvect(xArray, yArray, maxx, maxy, 1, nx2, 1, ny2, 
     &            dt/80., 0, tr2, 0.)
      call pgslw(1)
      call pgsci(1)
      call pgsch(1.3)
      call pgbox('bcnst',0,0,'bcnst',0,0)
      call pglabel(" ", " ", "Intensity plus polarization vectors")

c      call pgwedg("R", 2., 2., fg, bg, " ")
      call mywedge(bg1,fg1)



      call pgpage
      endif

      if (dogreyscale) then
         wstart = aspec(1)
         wend = aspec(nspec)
         trg(2) = (wend-wstart)/real(nxg-1)
         trg(1) = wstart - trg(2)/2.
         trg(3) = 0.
         trg(6) = 1./real(nyg-1)
         trg(4) = -trg(6)/2.
         trg(5) = 0.
         call pgvport(0.1,0.9,0.1,0.9)
         call pgsitf(0)
         call pgscir(cindex(3,1),cindex(3,2))
         call palette(1)
         call pgwnad(0.,1.,0.,1.)
         fgg = -1.e30
         bgg = 1.e30
         mean = 0.
        do i = 1, nxg
            do j = 2, nyg
               fgg = max(fgg, greyscale(i,j))
               bgg = min(bgg, greyscale(i,j))
               mean = mean + greyscale(i,j)
            enddo
         enddo
         mean = mean/(real(nxg)*real(nyg-1))
         sig = 0.
        do i = 1, nxg
            do j = 2, nyg
               sig = sig + (greyscale(i,j)-mean)**2
            enddo
         enddo
         sig = sqrt(mean/(real(nxg)*real(nyg-1)))
         fgg = fgg * 0.7


         write(*,*) "GREYSCALE",fgg,bgg,nxg,nyg
         call pgwindow(wstart, wend, 0., 1.)
         write(*,*) "plotting greyscale",nxg,nyg,fgg,bgg
         call pgscir(cindex(3,1),cindex(3,2))
         call palette(1)
         call pgimag(greyscale, nxg, nyg, 1, 
     &               nxg, 1, nyg, bgg, fgg, trg)
         call pgsci(3)
         call pgmove(wstart, phase)
         call pgdraw(wend, phase)
         call pgsci(1)
         call pgbox('bcnst',0,0,'bcnst',0,0)
         call pglabel(" ", " ", "Greyscale")
         call pgpage
      endif



      if (dopol) then
      if (first) then
         fg2 =-1.e10
         bg2 = 1.e10
      endif
      do i = 1 , nx
         do j = 1, ny
            temp(i,j) = 
     &             max(1.e-10,sqrt(qimage(i,j)**2 + uimage(i,j)**2))
            if (first) then
               if (temp(i,j) .gt. fg2) fg2 = temp(i,j)
               if (temp(i,j) .lt. bg2) bg2 = temp(i,j)
            endif
         enddo
      enddo
      
      write(*,*) bg2,fg2
c      bg = bg * 1.1
c      bg = 0.1
c      fg = fg / 10.

      
      call pgvport(0.1,0.9,0.1,0.9)
      call pgwnad(ax1(1), ax1(nx), ax2(1), ax2(ny))
      call pgscir(cindex(4,1),cindex(4,2))
      call palette(4)

      if (dolog) then
         call pgsitf(1)
      else
         call pgsitf(0)
      endif


      call pgimag(temp, maxx, maxy, 1, nx, 1, ny, bg2, fg2, tr)
      call pgbox('bcnst',0,0,'bcnst',0,0)
      call pglabel(" ", " ", "Polarized intensity image")

c      call pgwedg("R", 2., 2., fg, bg, " ")
      call mywedge(bg2,fg2)

      call pgpage

      endif

      IF (CIRCULAR) THEN
         if (first) then
            fg3 =-1.e10
            bg3 = 1.e10
         endif
         do i = 1 , nx
            do j = 1, ny
               temp(i,j) =   100. * vimage(i,j) / iimage(i,j) 
               if (first) then
                  if (temp(i,j) .gt. fg3) fg3 = temp(i,j)
                  if (temp(i,j) .lt. bg3) bg3 = temp(i,j)
               endif
            enddo
         enddo

         junk = max(abs(bg3),abs(fg3))
         bg3 = -junk
         fg3 = junk
         write(*,*) "circ",bg3,fg3
         if (first) then
            if ((bg3.eq.0.).and.(fg3.eq.0.)) then
               bg3 = -100.
               fg3 = 100.
            endif
         endif
         
         call pgvport(0.1,0.9,0.1,0.9)
         call pgwnad(ax1(1), ax1(nx), ax2(1), ax2(ny))
         call pgscir(cindex(2,1),cindex(2,2))
         call palette(2)
         call pgsitf(0)
         call pgimag(temp, maxx, maxy, 1, nx, 1, ny, bg3, fg3, tr)
         call pgbox('bcnst',0,0,'bcnst',0,0)
         call pglabel(" ", " ", "Circular polarization (%)")
         
         call mywedge(fg3,bg3)
      else
         if (first) then
            fg3 =-1.e10
            bg3 = 1.e10
         endif
         mean = 0.
         do i = 1 , nx
            do j = 1, ny
               if (velimage(i,j) .ne. 0.) then
                  mean = mean + velimage(i,j)
                  n = n + 1
               endif
            enddo
         enddo
         mean = mean / real(n)
         mean = 0.
         sig = 0.
         do i = 1 , nx
            do j = 1, ny
               temp(i,j) = velimage(i,j)
               if (velimage(i,j) .ne. 0.) then
                  sig = sig +  (velimage(i,j)-mean)**2
               endif
            enddo
         enddo
         sig = sqrt (sig / real(n))


         if (first) then
            fg3 = mean + sig
            bg3 = mean - sig
         endif
         
         
         if (first) then
            if ((bg3.eq.0.).and.(fg3.eq.0.)) then
               bg3 = -100.
               fg3 = 100.
            endif
         endif

c         bg3 = -300.
c         fg3 = 300.

         write(*,*) "vel",bg3,fg3
         call pgvport(0.1,0.9,0.1,0.9)
         call pgwnad(ax1(1), ax1(nx), ax2(1), ax2(ny))
         call pgscir(cindex(2,1),cindex(2,2))
         call palette(2)
         call pgsitf(0)
         call pgimag(temp, maxx, maxy, 1, nx, 1, ny, bg3, fg3, tr)
         call pgbox('bcnst',0,0,'bcnst',0,0)
         call pglabel(" ", " ", "Velocity (km/s)")
         
         call mywedge(bg3,fg3)

      endif


      if (dogreyscale) then
         if (first) then
            imax = -1.e30
            imin =  1.e30
            do i = 1, nspec
               imax = max(imax, spec(i))
               imin = min(imin, spec(i))
            enddo
            ir = imax - imin
            imax = imax + 0.2*ir
            imin = imin - 0.2*ir
         endif
         call pgpage
         call pgvport(0.1,0.9,0.1,0.9)
         call pgwnad(0.1,0.9,0.1,0.9)
         call pgwindow(wstart, wend, imin, imax)
         call pgbox('bcnst',0,0,'bcnst',0,0)
         call pgsci(3)
         call pgbin(nspec, aspec, spec, .true.)
         call pgsci(1)
         call pglabel(" ", " ", "Spectrum")
      endif

 666  continue
      call pgend
      end

      SUBROUTINE PALETTE(TYPE)
C-----------------------------------------------------------------------
C Set a "palette" of colors in the range of color indices used by
C PGIMAG.
C-----------------------------------------------------------------------
      INTEGER TYPE
C
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
      REAL TL(4), TR(4), TG(4), TB(4)
C
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
C
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
C
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
C
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
C
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
     :         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
     :         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
     :         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
     :         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
C
      DATA TL /0.0, 0.5, 0.5, 1.0/
      DATA TR /0.2, 0.6, 0.6, 1.0/
      DATA TG /0.0, 0.0, 0.5, 1.0/
      DATA TB /1.0, 0.0, 0.0, 0.0/
C
      IF (TYPE.EQ.1) THEN
C        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, 1.0, 0.5)
      ELSE IF (TYPE.EQ.2) THEN
C        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, 1.0, 0.5)
      ELSE IF (TYPE.EQ.3) THEN
C        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, 1.0, 0.5)
      ELSE IF (TYPE.EQ.4) THEN
C        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, 1.0, 0.5)
      ELSE IF (TYPE.EQ.5) THEN
C        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, 1.0, 0.5)
      ELSE IF (TYPE.EQ.6) THEN
C        -- TJP
         CALL PGCTAB(TL, TR, TG, TB, 4, 1.0, 0.5)
      END IF
      END

      subroutine mywedge(bg,fg)
      implicit none
      real fg, bg, tr(6),dx, r
      integer nx, ny
      parameter (nx=10,ny=100)
      real wedge(nx,ny)
      integer i,j,itf
      real x1,x2,y1,y2,xh,yh

      call pgqitf(itf)

      call pgqcs(0,xh,yh)
      do i = 1, nx
         do j = 1 , ny
            r = real(j-1)/real(ny-1)
c            if (itf .eq. 0) then
               if (fg .gt. bg) then
                  wedge(i,j) = bg + (fg-bg)*r
               else
                  wedge(i,j) = bg + (fg-bg)*r
               endif
c            endif
c            if (itf .eq. 1) then
c               wedge(i,j) = 10.**(log10(bg) + (log10(fg)-log(bg))*r)
c            endif
c            if (itf .eq. 2) then
c               wedge(i,j) = (sqrt(bg) + (sqrt(fg)-sqrt(bg))*r)**2
c            endif
         enddo
      enddo
      write(*,*) "w",wedge(1,1),bg
      write(*,*) "w",wedge(1,ny),fg

      call pgqvp(0,x1,x2,y1,y2)
      
      call pgvport(x2, x2+1*xh,y1,y2)
      call pgswin(1.,10.,min(bg,fg),max(bg,fg))
      tr(1) = 0.
      tr(2) = 1.
      tr(3) = 0.
      dx = (fg - bg)/99.
      tr(4) = bg-dx
      tr(5) = 0.
      tr(6) = dx
      call pgsitf(0)
      call pgimag(wedge, nx, ny, 1, nx, 1, ny, bg, fg, tr)
      if (itf.eq.0) then
         call pgswin(1.,10.,min(bg,fg),max(bg,fg))
         call pgbox("bc", 0., 0, "bcmst", 0., 0)
      endif
      if (itf.eq.1) then
         call pgswin(1.,10.,log10(min(bg,fg)),log10(max(bg,fg)))
         call pgbox("bc", 0., 0, "bcmst", 0., 0)
      endif
      end

      SUBROUTINE SORT(N,RA)
      DIMENSION RA(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END




