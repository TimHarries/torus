      implicit none
      character*80 filename, datfilename, sdffilename
      integer i1,nSize, nPix, i, j, nphase,iphase
      real scale
      real thisImageI(-50:50,-50:50)
      real thisImageQ(-50:50,-50:50)
      real thisImageU(-50:50,-50:50)
      real thisImageV(-50:50,-50:50)
      real thisImageVel(-50:50,-50:50)

      write(*,*) "filename?"
      read(*,'(a)') filename

      write(*,*) "nphase?"
      read(*,*) nphase

      i1 = index(filename," ")
      do iphase = 1, nphase
         write(datfilename,'(a,i3.3,a)')  filename(1:i1-1),iPhase,".dat"
         write(sdffilename,'(a,i3.3)')  filename(1:i1-1),iPhase

         write(*,*) datfilename,sdffilename
         open(20, file=datfilename, status="old", form="formatted")
         
         read(20,'(2I4)') nSize
      

         if (nSize .ne. 50) then
            write(*,*) "image is wrong size."
            stop
         endif
         
         nPix = 2*nSize + 1
         
         read(20, '(1p,e13.5)') scale
         
         do i = -nSize, nSize
            do j = -nSize, nSize
               read(20, '(1p,5e13.5)') thisImageI(i,j),
     &              thisImageQ(i,j),
     &              thisImageU(i,j),
     &              thisImageV(i,j),
     &              thisImageVel(i,j)
            enddo
         enddo
         close(20)
         
         call writeImagef77(thisImageI, thisImageQ, thisImageU,
     &        thisImageV, thisImageVel, scale, nSize, nPix, 
     &        sdffilename)
      enddo

      end





      subroutine writeImagef77(thisImageI, thisImageQ, thisImageU,
     &thisImageV, thisImageVel, scale, nSize, nPix, filename)
      implicit none
      include '/star/include/dat_par'
      include '/star/include/sae_par'
      integer npix, nsize
      real thisImageI(nPix, nPix)
      real thisImageQ(nPix, nPix)
      real thisImageU(nPix, nPix)
      real thisImageV(nPix, nPix)
      real thisImageVel(nPix, nPix)
      character*(*)  filename
      character*(DAT__SZLOC) loc, ploc
      integer place, ndfq, ndfu, ndfv, ndfvel
      integer qptr, uptr, vptr, velptr
      real scale
      integer  iptr, oaxisp1, oaxisp2, ndfo
      integer status
      integer i,j
      integer  lbnd(2), ubnd(2), nelm

       status = 0

       lbnd(1) = 1
       ubnd(1) = nPix
       lbnd(2) = 1
       ubnd(2) = nPix
       nelm = npix*npix

       write(*,'(a,i3,a,i3)') "Writing image: ",npix," by ",npix

!
! Begin the ndf and hds systems
!
      call ndf_begin
      call hds_start(status)
!
! create a new tsp ndf and map it
!
      call hds_new(filename,'output','ndf',0,0,loc,status)
      call dat_new(loc,'data_array','_real',2,ubnd,status)
      call ndf_find(loc, " ", ndfo, status)
      call ndf_acre(ndfo,status)
      call ndf_acput('xaxis',ndfo,'lab',1,status)
      call ndf_acput('rstar',ndfo,'unit',1,status)
      call ndf_acput('yaxis',ndfo,'lab',2,status)
      call ndf_acput('rstar',ndfo,'unit',2,status)
      call ndf_amap(ndfo,'centre',1,'_real','update',oaxisp1,
     &ubnd(1),status)
      call ndf_amap(ndfo,'centre',2,'_real','update',oaxisp2,
     &ubnd(1),status)
      call ndf_map(ndfo,'data','_real','write',iptr,nelm,status)

      if (status .ne. sai__ok) then
         write(*,'(a)') "! Error opening output ndf image"
         stop
      endif


      call fillNdfImage(scale,%val(iptr),%val(oaxisp1), %val(oaxisp2), 
     &npix, nSize, thisImageI)

      CALL NDF_XNEW(NDFO,'POLARIMETRY','EXT',0,0,PLOC,STATUS)
      CALL NDF_PLACE(PLOC,'STOKES_Q',PLACE,STATUS)
      CALL NDF_NEW('_REAL',2,LBND,UBND,PLACE,NDFQ,STATUS)
      CALL NDF_PLACE(PLOC,'STOKES_U',PLACE,STATUS)
      CALL NDF_NEW('_REAL',2,LBND,UBND,PLACE,NDFU,STATUS)
      CALL NDF_PLACE(PLOC,'STOKES_V',PLACE,STATUS)
      CALL NDF_NEW('_REAL',2,LBND,UBND,PLACE,NDFV,STATUS)
      CALL NDF_PLACE(PLOC,'VEL',PLACE,STATUS)
      CALL NDF_NEW('_REAL',2,LBND,UBND,PLACE,NDFVEL,STATUS)

      CALL NDF_MAP(NDFQ,'DATA','_REAL','WRITE',QPTR,UBND,STATUS)
      CALL NDF_MAP(NDFU,'DATA','_REAL','WRITE',UPTR,UBND,STATUS)
      CALL NDF_MAP(NDFV,'DATA','_REAL','WRITE',VPTR,UBND,STATUS)
      CALL NDF_MAP(NDFVEL,'DATA','_REAL','WRITE',VELPTR,UBND,STATUS)

      call fillNdfImageQUV(%val(qptr), npix, nSize, thisImageQ)
      call fillNdfImageQUV(%val(uptr), npix, nSize, thisImageU)
      call fillNdfImageQUV(%val(vptr), npix, nSize, thisImageV)
      call fillNdfImageQUV(%val(velptr), npix, nSize, thisImageVel)

      if (status .ne. sai__ok) then
         write(*,'(a)') "! Error filling output ndf image"
         stop
      endif
      call ndf_end(status)
      call hds_close(loc,status)
      call hds_stop(status)


      end 
       

      subroutine fillNdfImage(scale, ndfimage, axis1, axis2, nPix, 
     &nSize, thisImage)
      implicit none
      integer  nPix, nsize
      real thisImage(nPix,nPix), scale
      real  axis1(nPix)
      real  axis2(nPix)
      real  ndfImage(nPix,nPix)
      integer  i, j

      do i = -nSize, nSize
         axis1(i+nSize+1) = real(i) * scale
         do j = -nSize, nSize
            axis2(j+nSize+1) = real(j) * scale
            ndfImage(i+nSize+1,j+nSize+1) = 
     &           thisImage(i+nSize+1,j+nSize+1)
         enddo
      enddo
      end 
 
      subroutine fillNdfImageQUV(ndfimage,  nPix, nSize, thisImage)
      implicit none
      integer  nPix, nsize
      real thisImage(nPix,nPix), scale
      real  axis1(nPix)
      real  axis2(nPix)
      real  ndfImage(nPix,nPix)
      integer  i, j

      do i = -nSize, nSize
         do j = -nSize, nSize
            ndfImage(i+nSize+1,j+nSize+1) = 
     &           thisImage(i+nSize+1,j+nSize+1)
         enddo
      enddo
      end 
          
      subroutine writePVimageF77(nx, ny, array, xaxis, yaxis,
     &                           filename)
      implicit none
      include '/star/include/dat_par'
      include '/star/include/sae_par'
      integer nx, ny
      real array(nx,ny)
      real xaxis(nx)
      real yaxis(ny)
      character*(*) filename
      integer status, lbnd(2), ubnd(2)
      integer iptr, oaxisp1, oaxisp2, ndfo, nelm
      character*(DAT__SZLOC) loc

      lbnd(1) = 1
      lbnd(2) = 1
      ubnd(1) = nx
      ubnd(2) = ny
      nelm = nx * ny

      write(*,*) "writing",nx,ny," pv image to ",filename(1:20)
!
! Begin the ndf and hds systems
!
      status = 0

      call ndf_begin
      call hds_start(status)
!
! create a new tsp ndf and map it
!
      call hds_new(filename,'output','ndf',0,0,loc,status)
      call dat_new(loc,'data_array','_real',2,ubnd,status)
      call ndf_find(loc, " ", ndfo, status)
      call ndf_acre(ndfo,status)
      call ndf_acput('v',ndfo,'lab',1,status)
      call ndf_acput('km/s',ndfo,'unit',1,status)
      call ndf_acput('p',ndfo,'lab',2,status)
      call ndf_acput('arcsec',ndfo,'unit',2,status)
      call ndf_amap(ndfo,'centre',1,'_real','update',oaxisp1,
     &ubnd(1),status)
      call ndf_amap(ndfo,'centre',2,'_real','update',oaxisp2,
     &ubnd(1),status)
      call ndf_map(ndfo,'data','_real','write',iptr,nelm,status)

      if (status .ne. sai__ok) then
         write(*,'(a)') "! Error opening output ndf image"
         stop
      endif

      call copyImage(nx, ny, array, xAxis, yAxis, %val(iptr),
     &               %val(oaxisp1), %val(oaxisp2))


      call ndf_end(status)
      call hds_close(loc,status)
      call hds_stop(status)

      end
          
      subroutine copyImage(nx, ny, inArray, inX, inY, outArray,
     &                     outX, outY)
      implicit none
      integer nx, ny
      real inArray(nx,ny), outArray(nx,ny)
      real inX(nx), outX(nx)
      real inY(ny), outY(ny)
      integer i, j
      do i = 1, nx
         do j = 1, ny
            outArray(i,j) = inArray(i,j)
         enddo
      enddo
      do i = 1, nx
         outX(i) = inX(i)
      enddo
      do j = 1, ny
         outY(j) = inY(j)
      enddo
      end

            
                         





