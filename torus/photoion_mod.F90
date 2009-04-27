! photoionization module - started on October 4th 2005 by th

module photoion_mod

use math_mod
use parallel_mod
use gridio_mod
use image_mod
use source_mod
use timing
use grid_mod
use amr_mod
use constants_mod
use messages_mod
use diffusion_mod
use photon_mod
use unix_mod, only: unixGetenv

implicit none

private
#ifdef MPI
 private :: updateGridMPIphoto
#endif
public :: photoIonizationLoop, createImage, refineLambdaArray

type SAHAMILNETABLE
   integer :: nFreq 
   integer :: nTemp 
   real(double), pointer :: temp(:)
   real(double), pointer :: freq(:)
   real(double), pointer :: Clyc(:, :)
   real(double), pointer :: emissivity(:)
end type SAHAMILNETABLE

type RECOMBTABLE
   integer :: nrho
   integer :: nTemp 
   real(double), pointer :: rho(:)
   real(double),  pointer :: temp(:)
   real(double), pointer :: emissivity(:, :)
end type RECOMBTABLE


type GAMMATABLE
   integer :: nGamma
   integer :: nFreq
   integer :: nTemp
   real(double), pointer :: freq(:)
   real(double), pointer :: temp(:)
   real(double), pointer :: gamma(:,:)
end type GAMMATABLE



  real :: heIRecombinationFit(32, 3, 3)
  real :: heIRecombinationLambda(32)
  real :: heIRecombinationNe(3)
  real :: heIIrecombinationLines(3:30, 2:16)
  type(GAMMATABLE) :: gammaTableArray(3) ! H, HeI, HeII

contains


  subroutine photoIonizationloop(grid, source, nSource, nLambda, lamArray, readlucy, writelucy, &
       lucyfileout, lucyfilein)
    use input_variables, only : nlucy, taudiff !, smoothFactor
#ifdef MPI
    use input_variables, only : blockHandout
#endif
    implicit none
#ifdef MPI
 include 'mpif.h'
#endif

    type(GRIDTYPE) :: grid
    character(len=*) :: lucyfileout, lucyfilein
    logical :: readlucy, writelucy
    type(OCTAL), pointer :: thisOctal, tempOctal
    integer :: nCellsInDiffusion
    logical :: directPhoton
    character(len=80) :: message
    integer :: tempSubcell
    integer :: nlambda
    real :: lamArray(:)
    integer :: nSource
    type(SOURCETYPE) :: source(:), thisSource
    integer :: iSource
    type(VECTOR) :: rVec, uHat, rHat
    real(double) :: lCore
    integer :: nMonte, iMonte
    integer :: subcell
    integer :: i
    logical :: escaped
    real(double) :: wavelength, thisFreq
    real :: thisLam
    type(VECTOR) :: octVec
    real(double) :: r
    integer :: ilam
    integer :: nInf
    real(double) :: kappaScadb, kappaAbsdb
    real(double) :: epsOverDeltaT
    integer :: nIter
    logical :: converged
    real :: temp
    real(double) :: luminosity1, luminosity2, luminosity3
    real(double) :: photonPacketWeight
    real(double) :: fac
!    logical :: gridConverged
    real(double) :: albedo

    logical :: ok

    type(SAHAMILNETABLE) :: hTable, heTable
    type(RECOMBTABLE) :: Hrecombtable

    real(double) :: freq(1000), dfreq(1000), spectrum(1000), nuStart, nuEnd
    real(double) :: r1, kappaAbsGas, kappaAbsDust, escat, totalFlux

    integer, parameter :: nFreq = 1000
    logical, save :: firsttime = .true.
    integer :: iMonte_beg, iMonte_end, nSCat
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: np, iOctal, iOctal_beg, iOctal_end, nOctal
    integer ::  nVoxels, nOctals

#ifdef MPI
    ! For MPI implementations
  ! For MPI implementations =====================================================
    integer ::   my_rank        ! my processor rank
    integer ::   n_proc         ! The number of processes
    integer ::   ierr           ! error flag
    integer, dimension(:), allocatable :: photonBelongsRank
    integer, parameter :: tag = 0
    logical :: rankComplete

    integer, dimension(:), allocatable :: octalsBelongRank
    real, allocatable :: tempArray(:), tArray(:)
    real(double), allocatable :: tempArrayd(:), tArrayd(:)

    ! FOR MPI IMPLEMENTATION=======================================================
    !  Get my process rank # 
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  
    ! Find the total # of precessor being used in this run
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, ierr)
    
    if (my_rank .eq. 0) then
       print *, ' '
       print *, 'Photoionization loop computed by ', n_proc, ' processors.'
       print *, ' '
    endif
    
    ! ============================================================================
#endif


    ok = .true.
    rhat = VECTOR(0.d0, 0.d0, 0.d0)

    escat = 0.d0
    nuStart = cSpeed / (1000.e4 * 1.d-8)
    nuEnd =  2.d0*maxval(grid%ion(1:grid%nIon)%nuThresh)
    heTable%nfreq = 0
    htable%nfreq = 0
    hrecombTable%nTemp = 0
    kappaAbsDb = 0.d0; kappaScaDb = 0.d0
    iSource = 0
    kappaAbsGas = 0.d0; kappaAbsDust = 0.d0
    luminosity1 = 0.d0; luminosity2 = 0.d0; luminosity3 = 0.d0
    wavelength = 0.; temp = 0.


    do i = 1, nFreq
       freq(i) = log10(nuStart) + dble(i-1)/dble(nFreq-1) * (log10(nuEnd)-log10(nuStart))
       freq(i) = 10.d0**freq(i)
    enddo
    do i = 2, nFreq-1
       dfreq(i) = (freq(i+1)-freq(i-1))/2.d0
    enddo
    dfreq(1) = 2.d0*(freq(2)-freq(1))
    dfreq(nfreq) = 2.d0*(freq(nfreq)-freq(nfreq-1))

    do i = 1, grid%nIon
       call addxSectionArray(grid%ion(i), nfreq, freq)
    enddo

    call createGammaTable(gammaTableArray(1), 'gammaHI.dat')

    call createGammaTable(gammaTableArray(2), 'gammaHeI.dat')
    call createGammaTable(gammaTableArray(3), 'gammaHeII.dat')

    call createSahaMilneTables(hTable, heTable)

    call createRecombTable(Hrecombtable, 'e1b.d')


    call readHeIRecombinationLinesFit()

    call readHeIIrecombination()

    lCore = 0.d0
    do i = 1, nSource
       lCore = lCore + source(i)%luminosity
    enddo

    write(message,'(a,1pe12.5)') "Total souce luminosity (lsol): ",lCore/lSol
    call writeInfo(message, TRIVIAL)


    nIter = 0
    
    converged = .false.

    if (readlucy) then
       call writeInfo("Reading lucy dumpfile")
       call readAmrGrid(lucyfilein,.false.,grid)
       converged = .true.
    endif


    call countVoxels(grid%octreeRoot, nOctals, nVoxels)  
    if (nLucy == 0) then
       nMonte = nVoxels * 100
    else
       nMonte = nlucy
    endif

    do while(.not.converged)
       nIter = nIter + 1
       nInf=0


       nCellsInDiffusion = 0
       call defineDiffusionOnRosseland(grid,grid%octreeRoot, taudiff, ndiff=nCellsInDiffusion)
       write(message,*) "Number of cells in diffusion zone: ", nCellsInDiffusion
       call writeInfo(message,IMPORTANT)

       epsoverdeltat = lcore/dble(nMonte)

       call zeroDistanceGrid(grid%octreeRoot)

       if (writeoutput) then
          write(*,*) "Running loop with ",nmonte," photons. Iteration: ",niter
       endif

       if (doTuning) call tune(6, "One photoionization itr")  ! start a stopwatch

       iMonte_beg = 1
       iMonte_end = nMonte

#ifdef MPI
  !====================================================================================
  ! Splitting the innerPhoton loop for multiple processors.
  if (my_rank == 0) then
     print *, ' '
     print *, 'photonLoop computed by ', n_proc-1, ' processors.'
     print *, ' '
  endif
  if (my_rank == 0) then
     ! we will use an array to store the rank of the process
     !   which will calculate each photon
     allocate(photonBelongsRank(nMonte))
    
     call mpiBlockHandout(n_proc,photonBelongsRank,blockDivFactor=10,tag=tag,&
                          setDebug=.false.)
     deallocate(photonBelongsRank) ! we don't really need this here. 
  end if
  !====================================================================================

    
    
  if (my_rank /= 0) then
    mpiBlockLoop: do  
      call mpiGetBlock(my_rank,imonte_beg, imonte_end,rankComplete,tag,setDebug=.false.)  
      if (rankComplete) exit mpiBlockLoop  
    
#endif

       mainloop: do iMonte = iMonte_beg, iMonte_end
          call randomSource(source, nSource, iSource)
          thisSource = source(iSource)
          call getPhotonPositionDirection(thisSource, rVec, uHat,rHat,grid)
          escaped = .false.

          call amrGridValues(grid%octreeRoot, rVec, foundOctal=tempOctal, &
            foundSubcell=tempsubcell)

          if (tempOctal%diffusionApprox(tempsubcell)) then
             call randomWalk(grid, tempOctal, tempSubcell, thisOctal, Subcell, temp, ok)
             if (.not.ok) cycle mainLoop
             directPhoton = .false.
             rVec = subcellCentre(thisOctal, subcell)
          endif


          call getWavelength(thisSource%spectrum, wavelength)

          photonPacketWeight = 1.d0
          thisFreq = cSpeed/(wavelength / 1.e8)


!          if ((thisFreq*hcgs*ergtoev) < 13.6) then
!             cycle mainloop
!          endif
          nScat = 0
          do while(.not.escaped)
             nScat = nScat + 1

             call toNextEventPhoto(grid, rVec, uHat, escaped, thisFreq, nLambda, lamArray, photonPacketWeight, nFreq, freq)
             if (.not. escaped) then

                thisLam = (cSpeed / thisFreq) * 1.e8
                call locate(lamArray, nLambda, real(thisLam), iLam)
                octVec = rVec 
                call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, foundsubcell=subcell,iLambda=iLam, &
                     lambda=real(thisLam), kappaSca=kappaScadb, kappaAbs=kappaAbsdb, grid=grid)

                albedo = kappaScaDb / (kappaAbsdb + kappaScadb)


                call random_number(r)
                if (r < albedo) then
                   uHat = randomUnitVector() ! isotropic scattering
                else
                   spectrum = 1.d-30

                   call returnKappa(grid, thisOctal, subcell, ilambda=ilam, kappaAbsDust=kappaAbsDust, kappaAbsGas=kappaAbsGas, &
                        kappaSca=kappaScadb, kappaAbs=kappaAbsdb, kappaScaGas=escat)


                   if ((thisFreq*hcgs*ergtoev) > 13.6) then ! ionizing photon
                      call random_number(r1)
                      if (r1 < (kappaAbsGas / (kappaAbsGas + kappaAbsDust))) then  ! absorbed by gas rather than dust
                         call addLymanContinua(nFreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
                         call addHigherContinua(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, GammaTableArray)
                         call addHydrogenRecombinationLines(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
                         call addFreeFreeContinua(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
!                        call addHeRecombinationLines(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
                         call addForbiddenLines(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
                      else
!                         call addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, nlambda, lamArray)
                      endif
                   else ! non-ionizing photon must be absorbed by dust
!                         call  addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, nlambda, lamArray)
                   endif
                   if (firsttime.and.writeoutput) then
                      firsttime = .false.
                      open(67,file="spec.dat",status="unknown",form="formatted")
                      do i = 1, nfreq
                         write(67,*) freq(i), spectrum(i)
                      enddo
                      close(67)
                   endif

                   thisFreq =  getPhotonFreq(nfreq, freq, spectrum)
!                   write(*,*) 1.d8*cspeed/thisFreq, thisFreq*hCgs*ergtoev
                   uHat = randomUnitVector() ! isotropic emission
                endif
                

             endif
             if (nScat > 100000) then
                write(*,*) "Nscat exceeded 10000, forcing escape"
                write(*,*) 1.e8*cspeed/thisFreq
                write(*,*) albedo, kappaScaDb, kappaAbsdb,escat

                thisLam = (cSpeed / thisFreq) * 1.e8
                call locate(lamArray, nLambda, real(thisLam), iLam)

                call returnKappa(grid, thisOctal, subcell, ilambda=ilam, kappaAbsDust=kappaAbsDust, kappaAbsGas=kappaAbsGas, &
                     kappaSca=kappaScadb, kappaAbs=kappaAbsdb, kappaScaGas=escat)
                write(*,*) "thislam",thislam,ilam, lamArray(ilam)
                write(*,*) lamArray(1:nLambda)
                write(*,*) "kappaAbsDust",kappaAbsDust
                write(*,*) "kappaAbsGas",kappaAbsGas
                write(*,*) "kappaSca",kappaScadb
                write(*,*) "kappaAbs",kappaAbsdb
                write(*,*) "onekappaabs",grid%oneKappaAbs(1,ilam)
                write(*,*) "onekappasca",grid%oneKappasca(1,ilam)


                escaped = .true.
             Endif
          enddo
          nInf = nInf + 1
       end do mainloop

#ifdef MPI
 if (.not.blockHandout) exit mpiblockloop        
    end do mpiBlockLoop  
  end if ! (my_rank /= 0)

       if(my_rank == 0) write(*,*) "Calling update_octal_MPI"

       call MPI_BARRIER(MPI_COMM_WORLD, ierr)


       call updateGridMPIphoto(grid)

       if(my_rank == 0) write(*,*) "Done update_octal_MPI"
#endif


       if (doTuning) call tune(6, "One photoionization itr")  ! stop a stopwatch

       epsOverDeltaT = (lCore) / dble(nMonte)

       call  identifyUndersampled(grid%octreeRoot)



    np = 1
#ifdef MPI
    my_rank = 1
#endif
    firstTime = .true.

    allocate(octalArray(grid%nOctals))
    nOctal = 0
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)
    if (nOctal /= grid%nOctals) then
       write(*,*) "Screw up in get octal array", nOctal,grid%nOctals
       stop
    endif

#ifdef MPI
    ! FOR MPI IMPLEMENTATION=======================================================
    !  Get my process rank # 
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  
    ! Find the total # of precessor being used in this run
    call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
    
    ! we will use an array to store the rank of the process
    !   which will calculate each octal's variables
    allocate(octalsBelongRank(size(octalArray)))
    
    if (my_rank == 0) then
       call mpiBlockHandout(np,octalsBelongRank,blockDivFactor=10,tag=tag,&
                            setDebug=.false.)
    
    endif
    ! ============================================================================
#endif

    
    ! default loop indices
    ioctal_beg = 1
    ioctal_end = nOctal

    if (doTuning) call tune(6, "Temperature/ion corrections")

    if (writeoutput) &
         write(*,*) "Calculating ionization and thermal equilibria"

#ifdef MPI
 if (my_rank /= 0) then
  blockLoop: do     
 call mpiGetBlock(my_rank,iOctal_beg,iOctal_end,rankComplete,tag,setDebug=.false.)
   if (rankComplete) exit blockLoop 
#endif

    do iOctal =  iOctal_beg, iOctal_end

       thisOctal => octalArray(iOctal)%content

       do i = 1 , 3
          call calculateIonizationBalance(grid,thisOctal, epsOverDeltaT)
          call calculateThermalBalance(grid, thisOctal, epsOverDeltaT)
       enddo

    enddo

#ifdef MPI
 if (.not.blockHandout) exit blockloop
 end do blockLoop        
end if ! (my_rank /= 0)


     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     ! have to send out the 'octalsBelongRank' array
     call MPI_BCAST(octalsBelongRank,SIZE(octalsBelongRank),  &
                    MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     call countVoxels(grid%octreeRoot,nOctal,nVoxels)
     allocate(tempArray(1:nVoxels))
     allocate(tArray(1:nVoxels))
     allocate(tempArrayd(1:nVoxels))
     allocate(tArrayd(1:nVoxels))
     tArray = 0.d0
     tempArray = 0.d0
     tArrayd = 0.d0
     tempArrayd = 0.d0
     call packTemperatures(octalArray, nVoxels, tArray,octalsBelongRank)
     call MPI_ALLREDUCE(tArray,tempArray,nVoxels,MPI_REAL,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
     tArray = tempArray
     call unpackTemperatures(octalArray, nVoxels, tArray)
     tArray = 0.d0
     tempArray = 0.d0
     call packne(octalArray, nVoxels, tArrayd,octalsBelongRank)
     call MPI_ALLREDUCE(tArrayd,tempArrayd,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
     tArrayd = tempArrayd
     call unpackne(octalArray, nVoxels, tArrayd)
     do i = 1, grid%nIon
       tArrayd = 0.d0
       tempArrayd = 0.d0
       call packIonFrac(octalArray, nVoxels, tArrayd, octalsBelongRank, i)
       call MPI_ALLREDUCE(tArrayd,tempArrayd,nVoxels,MPI_DOUBLE_PRECISION,&
           MPI_SUM,MPI_COMM_WORLD,ierr)
       tArrayd = tempArrayd
       call unpackIonFrac(octalArray, nVoxels, tArrayd, octalsBelongRank, i)
     enddo





     deallocate(tempArray, tArray, tempArrayd, tArrayd)
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     deallocate(octalsBelongRank)
#endif
       deallocate(octalArray)

    if (writeoutput) &
         write(*,*) "Finished calculating ionization and thermal equilibria"

       if (doTuning) call tune(6, "Temperature/ion corrections")
       call defineDiffusionOnRosseland(grid,grid%octreeRoot,taudiff)

       if (doTuning) call tune(6, "Gauss-Seidel sweeps")


       nCellsInDiffusion = 0
       call defineDiffusionOnUndersampled(grid%octreeroot, nDiff=nCellsInDiffusion)


       call solveArbitraryDiffusionZones(grid)
       call defineDiffusionOnRosseland(grid,grid%octreeRoot, taudiff, nDiff=nCellsInDiffusion)
!       call unsetOnDirect(grid%octreeRoot)
       write(message,*) "Number of cells in diffusion zone: ", nCellsInDiffusion
       call writeInfo(message,IMPORTANT)

       if (doTuning) call tune(6, "Gauss-Seidel sweeps")



       if (writeoutput) then
       call dumpLexington(grid, epsoverdeltat)

       fac = 2.06e37

       luminosity1 = 0.d0
       call getHbetaLuminosity(grid%octreeRoot, grid, luminosity1)
       write(*,'(a20,2f12.4)') "H beta :",luminosity1/1.e37,luminosity1/2.05e37
       fac = luminosity1
       

       call getForbiddenLineLuminosity(grid, "N II", 1.22d6, luminosity1)
       write(*,'(a20,2f12.4)') "N II (122 um):",(luminosity1)/fac,(luminosity1)/(0.034*2.05e37)

       call getForbiddenLineLuminosity(grid, "N II", 6584.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "N II", 6548.d0, luminosity2)
       write(*,'(a20,2f12.4)') "N II (6584+6548):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.730*2.05e37)

       call getForbiddenLineLuminosity(grid, "N III", 5.73d5, luminosity1)
       write(*,'(a20,2f12.4)') "N III (57.3 um):",(luminosity1)/fac,(luminosity1+luminosity2)/(0.292*2.05e37)


       call getForbiddenLineLuminosity(grid, "O I", 6300.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "O I", 6363.d0, luminosity2)
       write(*,'(a20,2f12.4)') "O I (6300+6363):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.0086*2.05e37)

       call getForbiddenLineLuminosity(grid, "O II", 7320.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "O II", 7330.d0, luminosity2)
       write(*,'(a20,2f12.4)') "O II (7320+7330):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.029*2.05e37)

       call getForbiddenLineLuminosity(grid, "O II", 3726.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "O II", 3729.d0, luminosity2)
       write(*,'(a20,2f12.4)') "O II (3726+3729):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(2.03*2.05e37)

       call getForbiddenLineLuminosity(grid, "O III", 5007.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "O III", 4959.d0, luminosity2)
       write(*,'(a20,2f12.4)') "O III (5007+4959):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(2.18*2.05e37)

       call getForbiddenLineLuminosity(grid, "O III", 4363.d0, luminosity3)
       write(*,'(a20,2f12.4)') "O III (4363):",(luminosity3)/1.e37,luminosity3/(0.0037*2.05e37)

       call getForbiddenLineLuminosity(grid, "O III", 518145.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "O III", 883562.d0, luminosity2)
       write(*,'(a20,2f12.4)') "O III (52+88um):,",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/((1.06+1.22)*2.05e37)


       call getForbiddenLineLuminosity(grid, "Ne II", 1.28d5, luminosity1)
       write(*,'(a20,2f12.4)') "Ne II (12.8um):",(luminosity1)/fac,luminosity1/(0.195*2.05e37)

       call getForbiddenLineLuminosity(grid, "Ne III", 1.56d5, luminosity1)
       write(*,'(a20,2f12.4)') "Ne III (15.5um):",(luminosity1)/fac,luminosity1/(0.322*2.05e37)

       call getForbiddenLineLuminosity(grid, "Ne III", 3869.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "Ne III", 3968.d0, luminosity2)
       write(*,'(a20,2f12.4)') "Ne III (3869+3968):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.085*2.05e37)


       call getForbiddenLineLuminosity(grid, "S II", 6716.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "S II", 6731.d0, luminosity2)
       write(*,'(a20,2f12.4)') "S II (6716+6731):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.147*2.05e37)

       call getForbiddenLineLuminosity(grid, "S II", 4068.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "S II", 4076.d0, luminosity2)
       write(*,'(a20,2f12.4)') "S II (4068+4076):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.008*2.05e37)

       call getForbiddenLineLuminosity(grid, "S III", 1.87d5, luminosity1)
       write(*,'(a20,2f12.4)') "S III (18.7um):",(luminosity1)/fac,luminosity1/(0.577*2.05e37)

       call getForbiddenLineLuminosity(grid, "S III", 9532.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "S III", 9069.d0, luminosity2)
       write(*,'(a20,2f12.4)') "S III (9532+9069):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(1.22*2.05e37)
    endif


!       if ((niter > 2).and.(nIter < 8)) then
!          call locate(grid%lamArray, grid%nLambda,900.,ilam)
!          if (writeoutput) write(*,*) "Smoothing adaptive grid structure for ionization..."
!          call smoothAMRgridIonization(grid%octreeRoot,grid,gridConverged,ilam,inheritprops = .false., interpProps = .false.)
!          if (writeoutput) write(*,*) "...grid smoothing complete"
!          call locate(grid%lamArray, grid%nLambda,400.,ilam)
!          if (writeoutput) write(*,*) "Smoothing adaptive grid structure for ionization..."
!          call smoothAMRgridIonization(grid%octreeRoot,grid,gridConverged,ilam,inheritprops = .false., interpProps = .false.)
!          if (writeoutput) write(*,*) "...grid smoothing complete"
!          if (writeoutput) write(*,*) "Smoothing adaptive grid structure..."
!          gridConverged = .false.
!          do
!             call smoothAMRgrid(grid%octreeRoot,grid,smoothFactor,gridConverged,inheritprops=.false., interpProps=.false.)
!             if (gridConverged) exit
!          end do
!          if (writeoutput) write(*,*) "...grid smoothing complete"
!       endif

!      if (niter < 6) nMonte = nMonte * 2


    if (writeoutput) call writeAmrGrid("photo_tmp.grid",.false.,grid)

    if (niter == 10) converged = .true. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 enddo

 thisOctal => grid%octreeRoot
 call writeInfo("Calculating continuum emissivities...",TRIVIAL)



  call  calcContinuumEmissivity(grid, thisOctal, nlambda, lamArray)
 call writeInfo("Done.",TRIVIAL)

 if (writelucy) then
    call writeAmrGrid(lucyfileout,.false.,grid)
 endif

 call writeMultiImages(grid, nSource, source, VECTOR(1.d0,0.d0,0.d0))

end subroutine photoIonizationloop


 SUBROUTINE toNextEventPhoto(grid, rVec, uHat,  escaped,  thisFreq, nLambda, lamArray, photonPacketWeight, nfreq, freq)

   type(GRIDTYPE) :: grid
   type(VECTOR) :: rVec,uHat, octVec, tvec
   type(OCTAL), pointer :: thisOctal, tempOctal
   type(OCTAL),pointer :: oldOctal
   type(OCTAL),pointer :: endOctal
   integer :: nfreq
   real(double) :: freq(:)
   integer :: endSubcell
   real(double) :: photonPacketWeight
   integer :: subcell, tempSubcell
   real(oct) :: tval, tau, r
   real :: lamArray(:)
   integer :: nLambda
   logical :: stillinGrid
   logical :: escaped
   real(double) :: kappaScaDb, kappaAbsDb
   real(oct) :: thisTau
   real(oct) :: thisFreq
   real(oct) :: thisLam
   integer :: iLam
   logical ::inFlow , ok
   real :: diffusionZoneTemp
!   real :: lambda
!   integer :: ilambda

   endSubcell = 0
    stillinGrid = .true.
    escaped = .false.
    diffusionZoneTemp = 0.
    kappaAbsDb =0.d0 ; kappaScaDb = 0.0

    ok = .true.
    thisLam = (cSpeed / thisFreq) * 1.e8
    call locate(lamArray, nLambda, real(thisLam), iLam)
    if ((ilam < 1).or.(ilam > nlambda)) then
       write(*,*) "ilam errro",ilam,thisLam
    endif

! select an initial random tau and find distance to next cell

    call random_number(r)
    tau = -log(1.0-r)

    call distanceToCellBoundary(grid, rVec, uHat, tval)
    tval = tval + 1.d-3*grid%halfSmallestSubcell


!    if (grid%octreeRoot%threed) then
!       call intersectCubeAMR(grid, rVec, uHat, tVal)
!    else
!       call intersectCubeAMR2D(grid, rVec, uHat, tVal)
!    endif


    octVec = rVec

    call locate(lamArray, nLambda, real(thisLam), iLam)


    call amrGridValues(grid%octreeRoot, octVec,  iLambda=iLam, lambda=real(thisLam), foundOctal=thisOctal, &
         foundSubcell=subcell, kappaSca=kappaScadb, kappaAbs=kappaAbsdb, &
         grid=grid, inFlow=inFlow)
    oldOctal => thisOctal

    if (inFlow) then
       thisTau = dble(tVal) * (kappaAbsdb + kappaScadb)
    else
       thisTau = 1.0e-28
    end if

! if tau > thisTau then the photon has traversed a cell with no interactions

    do while(stillinGrid .and. (tau > thisTau)) 

! add on the distance to the next cell

       rVec = rVec + tVal * uHat

       octVec = rVec


! check whether the photon has escaped from the grid

       if (.not.inOctal(grid%octreeRoot, octVec)) then
          stillinGrid = .false.
          escaped = .true.
       endif


! update the distance grid

       call updateGrid(grid, thisOctal, subcell, thisFreq, tVal, photonPacketWeight, ilam, nfreq, freq)
          


       if (stillinGrid) then
          
          call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, &
               foundSubcell=subcell)
          if (thisOctal%diffusionApprox(subcell)) then
             call randomWalk(grid, thisOctal, subcell,  endOctal, endSubcell, diffusionZoneTemp,  ok)
             if (.not.ok) goto 666
             rVec = subcellCentre(endOctal,endSubcell)
             octVec = rVec
             call amrGridValues(grid%octreeRoot, rVec, foundOctal=tempOctal, &
                  foundSubcell=tempsubcell)
             uHat = randomUnitVector()
          endif
       endif





! now if the photon is in the grid choose a new random tau

       if (stillingrid) then
          call random_number(r)
          tau = -log(1.0-r)
          call locate(lamArray, nLambda, real(thisLam), iLam)
          call amrGridValues(grid%octreeRoot, octVec, iLambda=iLam,lambda=real(thisLam), foundOctal=thisOctal, &
               foundSubcell=subcell, kappaSca=kappaScadb, kappaAbs=kappaAbsdb, &
               grid=grid, inFlow=inFlow)
          oldOctal => thisOctal

! calculate the distance to the next cell


          call distanceToCellBoundary(grid, rVec, uHat, tval, sOctal=thisOctal)

          tval = tval + 1.d-3*grid%halfSmallestSubcell

!          if (grid%octreeRoot%threed) then
!             call intersectCubeAMR(grid, rVec, uHat, tVal)
!          else
!             call intersectCubeAMR2D(grid, rVec, uHat, tVal)
!          endif
          octVec = rVec

! calculate the optical depth to the next cell boundary

          if (inFlow) then
             thisTau = dble(tVal) * (kappaAbsdb + kappaScadb)
          else
             thisTau = 1.0e-28
          end if

          if (tVal == 0.0d0) then
             escaped = .true.
             stillingrid = .false.
          endif
       endif ! still in grid

! if photon is still in grid and  tau > tau_to_the_next_cell then loop round again
! choosing a new tau

       
    enddo
    

! the photon may have escaped the grid...

    if (.not.inOctal(grid%octreeRoot, octVec))  escaped = .true.

 ! if not the photon must interact in this cell
       

    if (.not.escaped) then



       octVec = rVec
!       if (.not.inOctal(grid%octreeRoot, octVec)) then
!          write(*,*) "Error:: Photon location is out of boundaries, but its status is not ESCAPED."
!          write(*,*) "        .... [lucy_mod::toNextEventAMR]"
!          write(*,*) "octVec-centre = ",octVec-grid%octreeRoot%centre
!          write(*,*) "cell size = ",grid%octreeRoot%subcellsize
!          stop
!       endif
       if (dble(tau)/thisTau .gt. 1.d0) then
          write(*,*) "tau prob",tau,thisTau
       endif

       call amrGridValues(grid%octreeRoot, octVec, startOctal=oldOctal, iLambda=iLam,lambda=real(thisLam),&
            foundOctal=thisOctal, foundSubcell=subcell, & 
            kappaAbs=kappaAbsdb,kappaSca=kappaScadb, grid=grid, inFlow=inFlow)


       if (thisOctal%diffusionApprox(subcell)) then
          write(*,*) "Photon in diff zone but not escaped"
       endif

       if (.not.inFlow) kappaAbsdb =0.0d0


! update the distance grid

       if (thisTau > 0.d0) then

!          lambda = cSpeed*1.e8/thisFreq
!          call locate(lamArray, nLambda, lambda, ilambda)
          call updateGrid(grid, thisOctal, subcell, thisFreq, dble(tval)*dble(tau)/thisTau, photonPacketWeight, ilam, nfreq, freq)

          oldOctal => thisOctal
          
       endif

       if (tau > thisTau) then
          write(*,*) "tau > thistau"
          stop
       endif

! move the requisite distance within the cell and return

       tVec = rVec
       rVec = rVec + (dble(tVal)*dble(tau)/thisTau) * uHat

       if (.not.inOctal(grid%octreeRoot, rVec)) then  ! this is only needed due to floating point boundary issues
          escaped = .true.
          goto 666
       endif

       octVec = rVec
       call amrGridValues(grid%octreeRoot, octVec, startOctal=oldOctal, iLambda=iLam,lambda=real(thisLam),&
            foundOctal=thisOctal, foundSubcell=subcell)
       if (thisOctal%diffusionApprox(subcell)) then
          call randomWalk(grid, thisOctal, subcell,  endOctal, endSubcell, diffusionZoneTemp, ok)
          if (.not.ok) goto 666
          rVec = subcellCentre(endOctal,endSubcell)
          octVec = rVec
          call amrGridValues(grid%octreeRoot, rVec, foundOctal=tempOctal, &
            foundSubcell=tempsubcell)
          uHat = randomUnitVector()
       endif

    endif

666 continue

 end subroutine toNextEventPhoto


  subroutine intersectCube(grid, posVec, i1,i2,i3,direction, tval)
   use vector_mod
   use grid_mod
   implicit none
   type(GRIDTYPE) :: grid
   type(VECTOR) :: direction
   type(VECTOR) :: posVec, norm(6), p3(6)
   real(oct) :: t(6),tval,denom(6)
   integer :: i,j
   logical :: ok, thisOk(6)
   integer :: i1, i2, i3

   ok = .true.

   norm(1) = VECTOR(1., 0., 0.)
   norm(2) = VECTOR(0., 1., 0.)
   norm(3) = VECTOR(0., 0., 1.)
   norm(4) = VECTOR(-1., 0., 0.)
   norm(5) = VECTOR(0., -1., 0.)
   norm(6) = VECTOR(0., 0., -1.)

   p3(1) = VECTOR(grid%xAxis(i1+1), 0., 0.)
   p3(2) = VECTOR(0.,grid%yAxis(i2+1),0.)
   p3(3) = VECTOR(0.,0.,grid%zAxis(i3+1))
   p3(4) = VECTOR(grid%xAxis(i1), 0., 0.)
   p3(5) = VECTOR(0.,grid%yAxis(i2),0.)
   p3(6) = VECTOR(0.,0.,grid%zAxis(i3))

   thisOk = .true.
   
   do i = 1, 6

      denom(i) = norm(i) .dot. direction
      if (denom(i) /= 0.) then
         t(i) = (norm(i) .dot. (p3(i)-posVec))/denom(i)
      else
         thisOk(i) = .false.
         t(i) = 0.
      endif
      if (t(i) < 0.) thisOk(i) = .false.
!      if (denom > 0.) thisOK(i) = .false.
 enddo




  
  j = 0
  do i = 1, 6
    if (thisOk(i)) j=j+1
  enddo

  if (j == 0) ok = .false.
   
  if (.not.ok) then
     write(*,*) i1, i2, i3
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  tval = minval(t, mask=thisOk)
  tval = max(tval * 1.001d0,dble((grid%xAxis(2)-grid%xAxis(1))/1000.))


  if (tval == 0.) then
     write(*,*) i1, i2, i3,tval
     write(*,*) posVec
     write(*,*) grid%xAxis(i1),grid%yAxis(i2),grid%zAxis(i3)
     write(*,*) grid%xAxis(i1+1),grid%yAxis(i2+1),grid%zAxis(i3+1)
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  if (tval > 2.*(grid%xAxis(2)-grid%xAxis(1))) then
!     write(*,*) "tval too big",tval,i1,i2,i3,posvec
!     write(*,*) "direction",direction
!     write(*,*) t(1:6)
!     write(*,*) denom(1:6)
  endif


  end subroutine intersectCube 
  subroutine intersectCubeAMR(grid, posVec, direction, tval)
   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(VECTOR), intent(in) :: posVec
   type(VECTOR), intent(in) :: direction
   real(oct), intent(out) :: tval
   !
   type(VECTOR) :: norm(6), p3(6)
   type(OCTAL),pointer :: thisOctal
   type(VECTOR) :: subcen, point
   integer :: subcell
   
   real(oct) :: t(6),denom(6)
   integer :: i,j
   logical :: ok, thisOk(6)


   point = posVec

   call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
   subcen =  subcellCentre(thisOctal,subcell)
   ok = .true.

   norm(1) = VECTOR(1.0d0, 0.d0, 0.0d0)
   norm(2) = VECTOR(0.0d0, 1.0d0, 0.0d0)
   norm(3) = VECTOR(0.0d0, 0.0d0, 1.0d0)
   norm(4) = VECTOR(-1.0d0, 0.0d0, 0.0d0)
   norm(5) = VECTOR(0.0d0, -1.0d0, 0.0d0)
   norm(6) = VECTOR(0.0d0, 0.0d0, -1.0d0)

   p3(1) = VECTOR(subcen%x+thisOctal%subcellsize/2.0d0, subcen%y, subcen%z)
   p3(2) = VECTOR(subcen%x, subcen%y+thisOctal%subcellsize/2.0d0 ,subcen%z)
   p3(3) = VECTOR(subcen%x,subcen%y,subcen%z+thisOctal%subcellsize/2.0d0)
   p3(4) = VECTOR(subcen%x-thisOctal%subcellsize/2.0d0, subcen%y,  subcen%z)
   p3(5) = VECTOR(subcen%x,subcen%y-thisOctal%subcellsize/2.0d0, subcen%z)
   p3(6) = VECTOR(subcen%x,subcen%y,subcen%z-thisOctal%subcellsize/2.0d0)

   thisOk = .true.

   do i = 1, 6

      denom(i) = norm(i) .dot. direction
      if (denom(i) /= 0.0d0) then
         t(i) = (norm(i) .dot. (p3(i)-posVec))/denom(i)
      else
         thisOk(i) = .false.
         t(i) = 0.0d0
      endif
      if (t(i) < 0.) thisOk(i) = .false.
!      if (denom > 0.) thisOK(i) = .false.
 enddo

  
  j = 0
  do i = 1, 6
    if (thisOk(i)) j=j+1
  enddo

  if (j == 0) ok = .false.
   
  if (.not.ok) then
     write(*,*) "Error: j=0 (no intersection???) in lucy_mod::intersectCubeAMR. "
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  tval = minval(t, mask=thisOk)
  tval = max(tval * 1.001d0,dble(thisOctal%subCellSize/1000.))


  if (tval == 0.) then
     write(*,*) posVec
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  if (tval > sqrt(3.)*thisOctal%subcellsize) then
!     write(*,*) "tval too big",tval/(sqrt(3.)*thisOctal%subcellSize)
!     write(*,*) "direction",direction
!     write(*,*) t(1:6)
!     write(*,*) denom(1:6)
  endif


  end subroutine intersectCubeAMR

  subroutine intersectCubeAMR2D(grid, posVec, direction, tval)

! this is to find a cell intersection for a 2D AMR grid
! which is essentially a 2D-grid that projects into
! a 3D cylindrical coordinate system


   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(VECTOR), intent(in) :: posVec
   type(VECTOR), intent(inout) :: direction
   real(oct), intent(out) :: tval
   !
   type(OCTAL),pointer :: thisOctal
   type(VECTOR) :: subcen, point
   integer :: subcell
   real(double) :: compZ, currentZ
   real(double) :: distToZBoundary, distToXboundary
   real(oct) :: r1,r2,d,cosmu,x1,x2,distTor1,distTor2, theta, mu
   logical :: ok
   type(VECTOR) :: xHat, zHAt

   point = posVec

   call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
   subcen =  subcellCentre(thisOctal,subcell)
   ok = .true.


   r1 = subcen%x - thisOctal%subcellSize/2.d0
   r2 = subcen%x + thisOctal%subcellSize/2.d0
   d = sqrt(point%x**2+point%y**2)
   xHat = VECTOR(point%x, point%y,0.d0)
   call normalize(xHat)

   cosmu =((-1.d0)*xHat).dot.direction
   call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
   if (.not.ok) then
      write(*,*) "Quad solver failed in intersectcubeamr2d"
      direction = randomUnitVector()
      x1 = thisoctal%subcellSize/2.d0
      x2 = 0.d0
   endif
   distTor2 = max(x1,x2)

   theta = asin(max(-1.d0,min(1.d0,r1 / d)))
   cosmu = xHat.dot.direction
   mu = acos(max(-1.d0,min(1.d0,cosmu)))
   distTor1 = 1.e30
   if (mu  < theta ) then
      call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
      if (.not.ok) then
         write(*,*) "Quad solver failed in intersectcubeamr2d"
         direction = randomUnitVector()
         x1 = thisoctal%subcellSize/2.d0
         x2 = 0.d0
      endif
      distTor1 = max(x1,x2)
   endif
         
   distToXboundary = min(distTor1, distTor2)


   zHat = VECTOR(0.d0, 0.d0, 1.d0)
   compZ = zHat.dot.direction
   currentZ = point%z

   if (compZ /= 0.d0 ) then
      if (compZ > 0.d0) then
         distToZboundary = (subcen%z + thisOctal%subcellsize/2.d0 - currentZ ) / compZ
      else
         distToZboundary = abs((subcen%z - thisOctal%subcellsize/2.d0 - currentZ ) / compZ)
      endif
   else
      disttoZboundary = 1.e30
   endif

   tVal = min(distToZboundary, distToXboundary) +0.0001d0*grid%halfsmallestsubcell
   if (tVal > 1.e29) then
      write(*,*) tVal,compZ, distToZboundary,disttoxboundary
      write(*,*) "subcen",subcen
      write(*,*) "z",currentZ
   endif
   if (tval < 0.) then
      write(*,*) tVal,compZ, distToZboundary,disttoxboundary
      write(*,*) "subcen",subcen
      write(*,*) "z",currentZ
   endif

  end subroutine intersectCubeAMR2D

  recursive subroutine zeroDistanceGrid(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call zeroDistanceGrid(child)
                exit
             end if
          end do
       else
          thisOctal%distanceGrid(subcell) = 0.d0
          thisOctal%nCrossings(subcell) = 0
          thisOctal%undersampled(subcell) = .false.
          thisOCtal%nDiffusion(subcell) = 0.
          thisOctal%Hheating(subcell) = 0.d0
          thisOctal%Heheating(subcell) = 0.d0
          thisOctal%photoIonCoeff(subcell,:) = 0.d0
          thisOctal%distanceGrid(subcell) = 0.d0
       endif
    enddo
  end subroutine zeroDistanceGrid

  recursive subroutine getHbetaluminosity(thisOctal, grid, luminosity)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  type(GRIDTYPE) :: grid
  real(double) :: luminosity, v, hbeta
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call getHbetaLuminosity(child, grid, luminosity)
                exit
             end if
          end do
       else
          v = cellVolume(thisOctal, subcell)
          hbeta = (10.d0**(-0.870d0*log10(thisOctal%temperature(subcell))+3.57d0)) * &
               thisOctal%ne(subcell) * thisOctal%ionFrac(subcell, 2) * &
               thisOctal%nh(subcell)*grid%ion(1)%abundance*1.d-25
          luminosity = luminosity + hbeta * (v*1.d30)
       endif
    enddo
  end subroutine getHbetaluminosity

  subroutine calculateIonizationBalance(grid, thisOctal, epsOverDeltaT)
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    real(double) :: epsOverDeltaT
    integer :: subcell
    
    do subcell = 1, thisOctal%maxChildren

       if (.not.thisOctal%hasChild(subcell)) then

          if (thisOctal%inflow(subcell)) then
             if (.not.thisOctal%undersampled(subcell)) then
                call solveIonizationBalance(grid, thisOctal, subcell, thisOctal%temperature(subcell), epsOverdeltaT)
             else
                !                write(*,*) "Undersampled cell",thisOctal%ncrossings(subcell)
             endif
          endif
       endif
    enddo
  end subroutine calculateIonizationBalance

  subroutine calculateThermalBalance(grid, thisOctal, epsOverDeltaT)
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    real(double) :: epsOverDeltaT
    real(double) :: totalHeating
    integer :: subcell
    logical :: converged, found
    real :: t1, t2, tm
    real(double) :: y1, y2, ym, Hheating, Heheating, dustHeating
    real :: deltaT
    real :: underCorrection = 1.
    integer :: nIter
    heheating = 0.d0; hheating = 0.d0; dustHeating = 0.d0; totalHeating = 0.d0

    do subcell = 1, thisOctal%maxChildren

       if (.not.thisOctal%hasChild(subcell)) then
          
          if (thisOctal%inflow(subcell)) then
             if (.not.thisOctal%undersampled(subcell)) then
                
                
                call getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDeltaT)
                
                if (totalHeating == 0.d0) then
                   thisOctal%temperature(subcell) = 1.e-3
                else
                   nIter = 0
                   converged = .false.
                   
                   t1 = 1.
                   t2 = 100000.
                   found = .true.
                   
                   if (found) then
                      y1 = (HHecooling(grid, thisOctal, subcell, epsOverDeltat,t1) &
                           - totalHeating)
                      y2 = (HHecooling(grid, thisOctal, subcell, epsOverDeltat,t2) &
                           - totalHeating)
                      if (y1*y2 > 0.d0) then
                         if (HHecooling(grid, thisOctal, subcell, epsOverDeltat,t1) > totalHeating) then
                            tm = t1
                         else
                            tm  = t2
                         endif
                         converged = .true.
                      endif
                      
                      ! Find root of heating and cooling by bisection
                      
                      do while(.not.converged)
                         tm = 0.5*(t1+t2)
                         y1 = (HHecooling(grid, thisOctal, subcell, epsOverDeltat,t1) &
                              - totalheating)
                         y2 = (HHecooling(grid, thisOctal, subcell, epsOverDeltat,t2) &
                              - totalheating)
                         ym = (HHecooling(grid, thisOctal, subcell, epsOverDeltat,tm) &
                              - totalheating)
                         
                         if (y1*ym < 0.d0) then
                            t1 = t1
                            t2 = tm
                         else if (y2*ym < 0.d0) then
                            t1 = tm
                            t2 = t2
                         else
                            converged = .true.
                            tm = 0.5*(t1+t2)
                            write(*,*) t1, t2, y1,y2,ym
                         endif
                         
                         if (abs((t1-t2)/t1) .le. 1.e-2) then
                            converged = .true.
                         endif
                         niter = niter + 1
                      enddo
                   endif
                   deltaT = tm - thisOctal%temperature(subcell)
                   thisOctal%temperature(subcell) = &
                        max(thisOctal%temperature(subcell) + underCorrection * deltaT,1.e-3)
                   !             write(*,*) thisOctal%temperature(subcell), niter
                endif
             else
                !                write(*,*) "Undersampled cell",thisOctal%ncrossings(subcell)
             endif
          endif
       endif
    enddo
  end subroutine calculateThermalBalance


  function hRecombination(temperature) result (rate)
    real :: temperature
    real(double) :: rate, t0, t1, b

    t0 = sqrt(temperature/3.148d0)
    t1 = sqrt(temperature/7.036d5)
    b = 0.7480

    rate = 7.982d-11 / (t0*(1.d0+t0)**(1.d0-b) * (1.d0 + t1)**(1.d0+b))
  end function hRecombination

  function recombToGround(temperature) result (alpha1)
    real :: temperature
    real(double) :: alpha1

    alpha1 = 1.58d-13 * (temperature/1.d4)**(-0.53d0)  ! kenny's photo paper equation 24
  end function recombToGround

  
  function HHeCooling(grid, thisOctal, subcell, epsOverDeltaT, temperature, debug) result (coolingRate)
    type(OCTAL),pointer :: thisOctal
    integer :: subcell
    real(double) :: epsOverDeltaT
    type(GRIDTYPE) :: grid
    real(double) :: nHii, nHeii, ne, nh
    real :: temperature
    logical, optional :: debug
    real(double) :: coolingRate, crate, dustCooling
    real(double) :: gff
    real :: rootTbetaH(31) = (/ 8.287e-11, 7.821e-11, 7.356e-11, 6.982e-11, 6.430e-11, 5.971e-11, 5.515e-11, 5.062e-11, 4.614e-11, &
                               4.170e-11, 3.734e-11, 3.306e-11, 2.888e-11, 2.484e-11, 2.098e-11, 1.736e-11, 1.402e-11, 1.103e-11, &
                               8.442e-12, 6.279e-12, 4.539e-12, 3.192e-12, 2.185e-12, 1.458e-12, 9.484e-13, 6.023e-13, 3.738e-13, &
                               2.268e-13, 1.348e-13, 7.859e-14, 4.499e-14 /)

    real :: rootTbetaHe(18) = (/ 8.347e-11, 7.889e-11, 7.430e-11, 6.971e-11, 6.512e-11, 6.056e-11, 5.603e-11, 5.154e-11, 4.710e-11,&
                               4.274e-11, 3.847e-11, 3.431e-11, 3.031e-11, 2.650e-11, 2.291e-11, 1.960e-11, 1.660e-11, 1.394e-11 /)

    real :: logT(31) = (/ 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, &
                      5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0 /)
    real(double) :: fac, thisRootTbetaH, betaH, betaHe, thisRootTbetaHe
    real :: thisLogT
    integer :: n
    real(double) :: log10te,betarec, coolrec, betaff, coolff
    real :: ch12, ch13, ex12, ex13, th12, th13, coolcoll, te4, teused
    real(double) :: becool
    real :: kappap
    real, parameter                :: hcRyd = &    ! constant: h*c*Ryd (Ryd at inf used) [erg]
         & 2.1799153e-11


!    call solveIonizationBalance(grid, thisOctal, subcell,  temperature, epsOverdeltaT)
                
    crate = 0.d0; kappap = 0.d0
    nHii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,2) * grid%ion(2)%abundance
    nHeii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,4) * grid%ion(4)%abundance
    nh = thisOctal%nh(subcell)
    ne = thisOctal%ne(subcell)
                


    becool = 0.
    coolingRate = 0.d0

    gff = 1.1d0 + 0.34d0*exp(-(5.5d0 - log10(temperature))**2 / 3.d0)  ! Kenny equation 23

    coolingRate = 1.42d-27 * (nHii+nHeii) * ne * gff * sqrt(temperature)  ! Kenny equation 22 (free-free)

    thisLogT = log10(temperature)
    log10te = thisLogT
    teused = temperature


    ! cooling of gas due to FF radiation from H+ 
    ! fits to Hummer, MNRAS 268(1994) 109, Table 1. or  least square fitting to m=4

    betaFF = 1.0108464E-11 + 9.7930778E-13*log10Te - &
         & 6.6433144E-13*log10Te*log10Te + 2.4793747E-13*log10Te*log10Te*log10Te -&
         & 2.3938215E-14*log10Te*log10Te*log10Te*log10Te
    
    coolFF = Nhii*Ne*betaFF*kerg*TeUsed/sqrt(TeUsed)

    becool = becool + coolfF

    ! cooling of gas due to recombination of H+
    ! fits to Hummer, MNRAS 268(1994) 109, Table 1.  
    ! least square fitting to m=4
    betaRec = 9.4255985E-11 -4.04794384E-12*log10Te &
         & -1.0055237E-11*log10Te*log10Te +  1.99266862E-12*log10Te*log10Te*log10Te&
         & -1.06681387E-13*log10Te*log10Te*log10Te*log10Te
    
    coolRec = Nhii*Ne*betaRec*kerg*TeUsed/sqrt(TeUsed)
    

    becool = becool+ coolrec


    ! collisional excitation of Hydrogen
    ! Mathis, Ly alpha, beta
    ch12 = 2.47e-8
    ch13 = 1.32e-8
    ex12 = -0.228
    ex13 = -0.460
    th12 = 118338.
    th13 = 140252.
    te4 = temperature / 1.e4
    coolColl = 0.
    if (TeUsed > 5000.) then 
       
       coolColl = (ch12*exp(-th12/TeUsed)*Te4**ex12 + &
            & ch13*exp(-th13/TeUsed)*Te4**ex13) * &
            & hcRyd*(nh-nhii)*Ne
       
    else
       
       coolColl = 0.
       
    end if
    coolingrate = coolingrate + coolcoll

    becool = becool +coolcoll

    if (PRESENT(debug)) then
       if (debug) then
          if (becool /= 0.) then
             write(*,*) coolff/becool, coolrec/becool, coolcoll/becool
          endif
       endif
    endif

    if (coolingRate < 0.) then
       write(*,*) "negative ff cooling",nhii,nheii,ne,gff,sqrt(temperature)
    endif

    call locate(logT, 31, thisLogT, n)
    fac = (thisLogT - logT(n))/(logT(n+1)-logT(n))

    thisRootTbetaH = rootTbetaH(n) + (rootTbetaH(n+1)-rootTbetaH(n)) * fac


    if (thisLogT < logT(18)) then
       thisRootTbetaHe = rootTbetaHe(n) + (rootTbetaHe(n+1)-rootTbetaHe(n)) * fac
    else
       thisRootTbetaHe = 0.d0
    endif

    betaH = thisrootTbetaH / sqrt(temperature)
    betaHe = thisrootTbetaHe / sqrt(temperature)



    coolingRate = coolingRate +  ne * nhii * kerg * temperature * betaH



    if (ne * nhii * kerg * temperature * betaH < 0.) then
       write(*,*) "negative H cooling",ne,nhii,kerg,temperature,betah
    endif


    coolingRate = coolingRate + ne * nheii * kerg * temperature * betaHe

    if (ne * nheii * kerg * temperature * betaHe < 0.) then
       write(*,*) "negative He cooling",ne,nheii,kerg,temperature,betaHe
    endif


    call metalcoolingRate(grid%ion, grid%nIon, thisOctal, subcell, nh, ne, temperature, crate, debug=.true.)
    if (crate < 0.) then
       write(*,*) "total negative metal cooling",crate
    endif
!    write(*,*) coolingRate,crate,coolingRate/(coolingrate+crate)
    coolingRate = coolingRate + crate

    call returnKappa(grid, thisOctal, subcell, kappap=kappap)

    dustCooling = fourPi * kappaP * (stefanBoltz/pi) * temperature**4

    coolingRate = coolingRate + dustCooling

  end function HHeCooling

  subroutine updateGrid(grid, thisOctal, subcell, thisFreq, distance, photonPacketWeight, ilambda, nfreq, freq)
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: nFreq, iFreq
    real(double) :: freq(:)
    integer :: subcell
    real(double) :: thisFreq, distance, kappaAbs,kappaAbsDust
    integer :: ilambda
    real(double) :: photonPacketWeight
    integer :: i 
    real(double) :: fac, xSec
    kappaAbs = 0.d0; kappaAbsDust = 0.d0
    thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1


    fac = distance * photonPacketWeight / (hCgs * thisFreq)

    call locate(freq, nFreq, thisFreq, iFreq)

    do i = 1, grid%nIon

       xSec = returnxSec(grid%ion(i), thisFreq, iFreq=iFreq)
!       call phfit2(grid%ion(i)%z, grid%ion(i)%n, grid%ion(i)%outerShell , e , xsec)
       if (xSec > 0.d0) then
          thisOctal%photoIonCoeff(subcell,i) = thisOctal%photoIonCoeff(subcell,i) &
               + fac * xSec

!          thisOctal%photoIonCoeff(subcell,i) = thisOctal%photoIonCoeff(subcell,i) &
!               + distance * dble(xsec) / (dble(hCgs) * thisFreq) * photonPacketWeight
       endif

       ! neutral h heating

       if ((grid%ion(i)%z == 1).and.(grid%ion(i)%n == 1)) then
          thisoctal%hheating(subcell) = thisoctal%hheating(subcell) &
            + distance * xsec / (thisfreq * hcgs) &
            * ((hcgs * thisfreq) - (hcgs * grid%ion(i)%nuthresh)) * photonpacketweight
       endif

       ! neutral he heating

       if ((grid%ion(i)%z == 2).and.(grid%ion(i)%n == 2)) then
          thisoctal%heheating(subcell) = thisoctal%heheating(subcell) &
            + distance * xsec / (thisfreq * hcgs) &
            * ((hcgs * thisfreq) - (hcgs * grid%ion(i)%nuthresh)) * photonpacketweight
       endif

    enddo



    call returnkappa(grid, thisoctal, subcell, ilambda=ilambda, kappaabsdust=kappaabsdust, kappaabs=kappaabs)

    thisoctal%distancegrid(subcell) = thisoctal%distancegrid(subcell) &
         + dble(distance) * dble(kappaabsdust)


  end subroutine updategrid

subroutine solveIonizationBalance(grid, thisOctal, subcell, temperature, epsOverdeltaT)
  implicit none
  type(GRIDTYPE) :: grid
  type(OCTAL) :: thisOctal
  real(double) :: epsOverDeltaT, V
  real :: temperature
  integer :: subcell
  integer :: i, k
  integer :: nIonizationStages
  integer :: iStart, iEnd, iIon
  real(double) :: chargeExchangeIonization, chargeExchangeRecombination
  real(double), allocatable :: xplus1overx(:)
  chargeExchangeIonization = 0.d0; chargeExchangeRecombination = 0.d0
  v = cellVolume(thisOctal, subcell)
  k = 1
  do while(k <= grid%nIon)

     iStart = k
     iEnd = k+1
     do while(grid%ion(istart)%z == grid%ion(iEnd)%z)
        iEnd = iEnd + 1
     enddo
     iEnd = iEnd - 1
     nIonizationStages = iEnd - iStart + 1
     
     
     allocate(xplus1overx(1:nIonizationStages-1))
     do i = 1, nIonizationStages-1
        iIon = iStart+i-1
        call getChargeExchangeRecomb(grid%ion(iion+1), temperature, &
             thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,1),  &
             thisOctal%nh(subcell)*grid%ion(2)%abundance*thisOctal%ionFrac(subcell,2),  &
             chargeExchangeRecombination)
        
        call getChargeExchangeIon(grid%ion(iion), temperature, &
             thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,1),  &
             thisOctal%nh(subcell)*grid%ion(2)%abundance*thisOctal%ionFrac(subcell,2),  &
             chargeExchangeIonization)
        
        
        xplus1overx(i) = ((epsOverDeltaT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell, iIon) + chargeExchangeIonization) / &
             max(1.d-50,(recombRate(grid%ion(iIon),temperature) * thisOctal%ne(subcell) + chargeExchangeRecombination))
        !           if (grid%ion(iion)%species(1:1) =="C") write(*,*) i,xplus1overx(i)
     enddo
     thisOctal%ionFrac(subcell, iStart:iEnd) = 1.
     do i = 1, nIonizationStages - 1
        iIon = iStart+i-1
        thisOctal%ionFrac(subcell,iIon+1) = thisOctal%ionFrac(subcell,iIon) * xplus1overx(i)
        !           if (grid%ion(iion)%species(1:1) =="C") write(*,*) i,thisOctal%ionFrac(subcell,iIon)
     enddo
     if (SUM(thisOctal%ionFrac(subcell,iStart:iEnd)) /= 0.d0) then
        thisOctal%ionFrac(subcell,iStart:iEnd) = &
             max(1.d-50,thisOctal%ionFrac(subcell,iStart:iEnd))/SUM(thisOctal%ionFrac(subcell,iStart:iEnd))
     else
        thisOctal%ionFrac(subcell,iStart:iEnd) = 1.d-50
     endif
     

        deallocate(xplus1overx)

        k = iEnd + 1
     end do

     thisOctal%ne(subcell) = returnNe(thisOctal, subcell, grid%ion, grid%nion)

   end subroutine solveIonizationBalance


function recombRate(thisIon, temperature) result (rate)
  type(IONTYPE) :: thisIon
  real :: temperature
  real(double) :: rate
  real(double) :: a, b, t0, t1, c, d, f, t, z

! recombinations INTO this species

  select case(thisIon%z)

     case(1)
        select case(thisIon%n)
           case(1) ! H I
              a = 7.982e-11
              b = 0.7480
              t0 = 3.148e0
              t1 = 7.036e5
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate = 0.
        end select
     case(2)
        select case(thisIon%n)
           case(2) ! He I 
              a = 3.294e-11
              b = 0.6910
              t0 = 1.544e1
              t1 = 3.676e7
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case(1) ! He II
              a = 1.891e-10
              b = 0.7524
              t0 = 9.370e00
              t1 = 2.774e6
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select
     case(6)
        select case(thisIon%n)
           case(6) ! C I 
              a = 0.0108
              b = -0.1075
              c = 0.2810
              d = -0.0193
              f = -0.1127
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 1
              a = 5.068
              b = -0.6192
              c = -0.0815
              d = 1.2910
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(5) ! C II
              a = 1.8267
              b = 4.1012
              c = 4.8443
              d = 0.2261
              f = 0.5960
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 2
              a = 5.434
              b = -0.6116
              c = 0.0694
              d = 0.7866
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(4) ! C III
              a = 2.3196
              b = 10.7328
              c = 6.8830
              d = -0.1824
              f = 0.4101
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 3
              a = 4.742
              b = -0.6167
              c = 0.2960
              d = 0.6167
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(3) ! C IV
              a = 8.540e-11
              b = 0.5247
              t0 = 5.014e2
              t1 = 1.479e7
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case(2) ! C V
              a = 2.765e-10
              b = 0.6858
              t0 = 1.535e2
              t1 = 2.556e7
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select
     case(7)
        select case(thisIon%n)
           case(7) ! N I 
              a = 0.0
              b = 0.6310
              c = 0.1990
              d = -0.0197
              f = 0.4398
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 1
              a = 3.874
              b = -0.6487
              c = 0.
              d = 1.
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(6) ! N II
              a = 0.0320
              b = -0.6624
              c = 4.3191
              d = 0.0003
              f = 0.5946
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 2
              a = 4.974
              b = -0.6209
              c = 0.
              d = 1.
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(5) ! N III
              a = -0.8806
              b = 11.2406
              c = 30.7066
              d = -1.1721
              f = 0.6127
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 3
              a = 4.750
              b = -0.5942
              c = 0.8452
              d = 2.8450
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(4) ! N IV
              a = 0.4134
              b = -4.6319
              c = 25.9172
              d = -2.2290
              f = 0.2360
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 4
              a = 4.626
              b = -0.9521
              c = 0.4729
              d = -0.4477
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(3) ! N V
              a = 1.169e-10
              b = 0.5470
              t0 = 6.793e2
              t1 = 1.650e7
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select
     case(8)
        select case(thisIon%n)
           case(8) ! O I 
              t = temperature / 10000.
              if (t < 2.) then
                 a = -0.0001
                 b = 0.0001
                 c = 0.0956
                 d = 0.0193
                 f = 0.4106
              else
                 a = 0.3715
                 b = -0.0293
                 c = -0.0597
                 d = 0.0678
                 f = 0.7993
              endif
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 1
              a = 3.201
              b = -0.6880
              c = -0.0174
              d = 1.7070
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(7) ! O II
              a = -0.0036
              b = 0.7519
              c = 1.5252
              d =-0.0838
              f = 0.2769
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 2
              a = 4.092
              b = -0.6413
              c = 0.
              d = 1.
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(6) ! O III
              a = 0.0
              b = 21.8790
              c = 16.2730
              d = -0.7020
              f = 1.1899
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 3
              a = 4.890
              b = -0.6213
              c = 0.0184
              d = 1.5550
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(5) ! O IV
              t = temperature / 10000.
              if (t < 2.) then
                 a = -0.3648
                 b =  7.2698
                 c =  17.2187
                 d =  9.8335
                 f = -0.0166
              else
                 a = -2.5053
                 b = 3.4903
                 c = 67.4128
                 d = -3.4450
                 f = 0.8501
              endif
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 4
              a = 14.665
              b = -0.5140
              c = 2.7300
              d = 0.2328
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(4) ! O V
              a = -2.8425
              b = 0.2283
              c = 40.4072
              d = -3.4956
              f = 1.7558
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 5
              a = 4.626
              b = -0.9521
              c = 0.4729
              d = -0.4477
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select

     case(10)
        select case(thisIon%n) ! from nussbaumer and storey 1987
           case(10) ! Ne I
              z = 1
              a = 7.317
              b = -0.5358
              c = 2.4130
              d = 0.4176
              rate = ppb1991(z, a, b, c, d, dble(temperature))
           case(9) ! Ne II
              a = 0.0129
              b =-0.1779
              c = 0.9353
              d =-0.0682
              f = 0.4516
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 2
              a = 11.80
              b = -0.5404
              c = 3.0300
              d = 0.2050
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(8) ! Ne III
              a = 3.6781
              b = 14.1481
              c = 17.1175
              d = -0.5017
              f = 0.2313
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 3
              a = 5.841
              b = -0.5921
              c = 0.4498
              d = 0.6395
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(7) ! Ne IV
              a =-0.0254
              b = 5.5365
              c = 17.0727
              d = -0.7225
              f = 0.1702
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 4
              a = 15.550
              b = -0.4825
              c = 3.2740
              d = 0.3030
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(6) ! Ne V
              a = -0.0141
              b = 33.8479
              c = 43.1608
              d =-1.6072
              f = 0.1942
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 5
              a = 7.538
              b = -0.5540
              c = 1.2960
              d = 0.3472
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(5) ! Ne VI
              a = 19.9280
              b = 235.0536
              c = 152.5096
              d = 9.1413
              f = 0.1282
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 6
              a = 5.239
              b = -0.5966
              c = 0.7135
              d = 0.4537
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select

     case(16)
        select case(thisIon%n) ! from nussbaumer and storey 1987
           case(16) ! S I 
              rate = 3.e-13
              rate = rate + svs1982(dble(temperature), 4.10D-13, 6.30D-1)
           case(15) ! S II
              rate = 3.e-30 ! page 1344, para 2, kenny's paper
              rate = rate + svs1982(dble(temperature), 1.80D-12, 6.86D-1)
           case(14) ! S  III
              rate = 1.5e-11
              rate = rate + svs1982(dble(temperature), 2.70D-12, 7.45D-1)
           case(13) ! S IV
              rate = 2.5e-11
              rate = rate + svs1982(dble(temperature), 5.70D-12, 7.55D-1)
        end select
     case DEFAULT
        write(*,*) "No recombination rate for ",thisIon%species
        rate  = 0.
  end select
end function recombRate

real function vernerFerland(a, t, t0, t1, b)
  real(double):: a, t, t0, t1, b

! based on Verner and Ferland (1996) ApJS 103 467

  vernerFerland = a / (sqrt(t/t0) * (1.d0+sqrt(t/t0))**(1.d0-b) * (1.d0+sqrt(t/t1))**(1.d0+b))
end function vernerFerland



function nussbaumerStorey1983(t, a, b, c, d, f) result(rate)

! based on Nussbaumer & Storey 1983 AA 126 75
  real(double) :: t, a, b, c, d, f, rate

  if (t < 0.05d0) then
     rate = tiny(rate)
  else
     rate = (1.d-12) * ((a/t) + b + (c*t) + d*(t**2))*(t**(-1.5d0))*exp(-f/t)
  endif

end function nussbaumerStorey1983

function ppb1991(z, a, b, c, d, temperature) result(rate)

! radiative recombination rates
! based on Pequignot, Petitjean & Boisson (1991) A&A 251 680

  real(double) :: z, a, b, c, d, temperature, t, rate

  t = 1.d-4 * temperature / z**2

  rate = 1.d-13 * z * (a*t**b)/(1.d0+c*t**d)

end function ppb1991

function svs1982(t, alpharad, xrad) result (rate)

! radiative recombination rates based on
! Shull and Van Steenberg, 1992, ApJS, 48, 95

  real(double) :: t, alpharad, xrad, rate

  rate = alpharad * (t /1.d4)**(-xrad)
end function svs1982

function returnNe(thisOctal, subcell, ionArray, nion) result (ne)
  real(double) :: ne, tot
  integer :: subcell
  type(OCTAL) :: thisOctal
  type(IONTYPE) :: ionArray(:)
  integer :: nion, i

  tot = 0.d0 
  do i = 1, nIon
     tot = tot + ionArray(i)%abundance * thisOctal%nh(subcell) * &
          thisOctal%ionFrac(subcell, i) * dble(ionArray(i)%z-ionArray(i)%n)
  enddo
  ne = tot
end function returnNe

subroutine dumpLexington(grid, epsoverdt)
  type(GRIDTYPE) :: grid
  type(OCTAL), pointer :: thisOctal
  integer :: subcell
  integer :: i, j
  real(double) :: r, theta
  real :: t,hi,hei,oii,oiii,cii,ciii,civ,nii,niii,niv,nei,neii,neiii,neiv
  real(double) :: oirate, oiirate, oiiirate, oivrate
  real(double) :: v, epsoverdt
  type(VECTOR) :: octVec
  real :: fac
  real(double) :: hHeating, heHeating, totalHeating, heating, nh, nhii, nheii, ne
  real(double) :: cooling, dustHeating
  real :: netot

  dustHeating = 0.d0; heHeating = 0.d0; hHeating = 0.d0; totalHeating = 0.d0

       call writeVtkFile(grid, "lexington.vtk", &
            valueTypeString=(/"rho        ", "temperature", "HI         "/))

  open(20,file="lexington.dat",form="formatted",status="unknown")
  open(21,file="orates.dat",form="formatted",status="unknown")
  open(22,file="ne.dat",form="formatted",status="unknown")

  do i = 1, 500
     r = (1.+7.d0*dble(i-1)/499.d0)*pctocm/1.e10

     t=0;hi=0; hei=0;oii=0;oiii=0;cii=0;ciii=0;civ=0;nii=0;niii=0;niv=0;nei=0;neii=0;neiii=0;neiv=0;ne=0.

     oirate = 0; oiirate = 0; oiiirate = 0; oivrate = 0
     heating = 0.d0; cooling = 0.d0; netot = 0.d0
     do j = 1, 100
        call random_number(theta)
        theta = theta * Pi
        
        octVec = VECTOR(r*sin(theta),0.d0,r*cos(theta))
        
        call amrgridvalues(grid%octreeRoot, octVec,  foundOctal=thisOctal, foundsubcell=subcell)

        nHii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,2) * grid%ion(2)%abundance
        nHeii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,4) * grid%ion(4)%abundance
        nh = thisOctal%nh(subcell)
        ne = thisOctal%ne(subcell)

!        cooling = cooling + HHeCooling(grid, thisOctal, subcell, thisOctal%temperature(subcell))

        v = cellVolume(thisOctal, subcell)

        HI = HI + thisOctal%ionfrac(subcell,returnIonNumber("H I", grid%ion, grid%nIon))
        HeI = HeI + thisOctal%ionfrac(subcell,returnIonNumber("He I", grid%ion, grid%nIon))
        OII = OII + thisOctal%ionfrac(subcell,returnIonNumber("O II", grid%ion, grid%nIon))
        OIII = OIII + thisOctal%ionfrac(subcell,returnIonNumber("O III", grid%ion, grid%nIon))
        CII = CII + thisOctal%ionfrac(subcell,returnIonNumber("C II", grid%ion, grid%nIon))
        CIII = CIII + thisOctal%ionfrac(subcell,returnIonNumber("C III", grid%ion, grid%nIon))
        CIV = CIV + thisOctal%ionfrac(subcell,returnIonNumber("C IV", grid%ion, grid%nIon))
        NII = NII + thisOctal%ionfrac(subcell,returnIonNumber("N II", grid%ion, grid%nIon))
        NIII = NIII + thisOctal%ionfrac(subcell,returnIonNumber("N III", grid%ion, grid%nIon))
        NIV = NIV + thisOctal%ionfrac(subcell,returnIonNumber("N IV", grid%ion, grid%nIon))
        NeI = NeI + thisOctal%ionfrac(subcell,returnIonNumber("Ne I", grid%ion, grid%nIon))
        NeII = NeII + thisOctal%ionfrac(subcell,returnIonNumber("Ne II", grid%ion, grid%nIon))
        NeIII = NeIII + thisOctal%ionfrac(subcell,returnIonNumber("Ne III", grid%ion, grid%nIon))
        NeIV = NeIV + thisOctal%ionfrac(subcell,returnIonNumber("Ne IV", grid%ion, grid%nIon))
        netot = netot + thisOctal%ne(subcell)
        call getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDT)
        heating = heating + totalHeating
        fac = thisOctal%nh(subcell) * returnAbundance(8) !* thisOctal%ionfrac(subcell,returnIonNumber("O I", grid%ion, grid%nIon))
        fac = 1.
        oirate = oirate + &
             fac*((epsOverDT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell,returnIonNumber("O I", grid%ion, grid%nIon)))
        oiirate = oiirate + &
             fac*((epsOverDT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell,returnIonNumber("O II", grid%ion, grid%nIon)))
        oiiirate = oiiirate + &
             fac*((epsOverDT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell,returnIonNumber("O III", grid%ion, grid%nIon)))
        oivrate = oivrate + &
             fac*((epsOverDT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell,returnIonNumber("O IV", grid%ion, grid%nIon)))

        t  = t + thisOctal%temperature(subcell)
     enddo


     hi = hi / 100.; hei = hei/100.; oii = oii/100; oiii = oiii/100.; cii=cii/100.
     ciii = ciii/100; civ=civ/100.; nii =nii/100.; niii=niii/100.; niv=niv/100.
     nei=nei/100.;neii=neii/100.; neiii=neiii/100.; neiv=neiv/100.;t=t/100.
     netot = netot / 100

     oirate = oirate / 100.
     oiirate = oiirate / 100.
     oiiirate = oiiirate / 100.
     oivrate = oivrate / 100.
     heating = heating / 100.
     cooling = cooling / 100.

     hi = log10(max(hi, 1e-10))
     hei = log10(max(hei, 1e-10))
     oii = log10(max(oii, 1e-10))
     oiii = log10(max(oiii, 1e-10))
     cii = log10(max(cii, 1e-10))
     ciii = log10(max(ciii, 1e-10))
     civ = log10(max(civ, 1e-10))
     nii = log10(max(nii, 1e-10))
     niii = log10(max(niii, 1e-10))
     niv= log10(max(niv, 1e-10))
     nei = log10(max(nei, 1e-10))
     neii = log10(max(neii, 1e-10))
     neiii = log10(max(neiii, 1e-10))
     neiv = log10(max(neiv, 1e-10))
     ne = log10(max(ne,1.d-10))


     write(21,'(f5.3,1p,6e12.3,0p)') r*1.e10/pctocm,heating,cooling,oirate,oiirate,oiiirate,oivrate

     write(20,'(f5.3,f9.1,  14f8.3)') &
          r*1.e10/pctocm,t,hi,hei,oii,oiii,cii,ciii,civ,nii,niii,niv,nei,neii,neiii,neiv
     write(22,*) r*1.e10/pctocm,netot
  enddo
  close(20)
  close(21)
  close(22)
end subroutine dumpLexington

subroutine getChargeExchangeRecomb(parentIon, temperature, nHI, nHII, recombRate)
  type(IONTYPE) :: parentIon
  real(double) :: nHI, nHII,  recombRate
  real :: t4, a, b, c, d
  real :: temperature

  recombRate  = 0.d0

  select case(parentIon%z)
     case(7)
        select case(parentIon%n)
           case(4) ! N IV
              t4 = temperature / 1.e4
              a = 3.05e-10
              b = 0.60
              c = 2.65
              d = -0.93
              recombRate  = kingdonFerland96(t4, a, b, c, d)
        end select
     case(8)
        select case(parentIon%n)
           case(7) ! O II
              t4 = temperature / 1.e4
              a = 1.04e-9
              b = 3.15e-2
              c = -0.61
              d = -9.73
              recombRate  = kingdonFerland96(t4, a, b, c, d)
           case(6) ! OIII
              t4 = temperature / 1.e4
              a = 1.04e-9
              b = 0.27
              c = 2.02
              d = -5.92
              recombRate  = kingdonFerland96(t4, a, b, c, d)
       end select
   end select

  recombRate = recombRate * nhi

end subroutine getChargeExchangeRecomb

subroutine getChargeExchangeIon(parentIon, temperature, nHI, nHII, IonRate)
  type(IONTYPE) :: parentIon
  real(double) :: nHI, nHII, ionRate
  real :: temperature
  real :: t4, a, b, c, d

  IonRate = 0.d0

  select case(parentIon%z)
     case(7)
        select case(parentIon%n)
           case(7) ! N I charge exchange Ionization
              t4 = temperature / 1.e4
              a = 4.55e-12
              b = -0.29
              c = -0.92
              d = -8.38
              ionRate  = kingdonFerland96(t4, a, b, c, d)
        end select
     case(8)
        select case(parentIon%n)
           case(8) ! O I charge exchange ionization
              t4 = temperature / 1.e4
              a = 7.40e-11
              b = 0.47
              c = 24.37
              d = -0.74
              ionRate  = kingdonFerland96(t4, a, b, c, d)
       end select
   end select

  ionRate = ionRate * nhii
end subroutine getChargeExchangeIon


function kingdonFerland96(t4, a, b, c, d) result (alpha)
  real :: alpha
  real :: t4, a, b, c, d
  alpha = a*(t4**b)*(1.+c*exp(d*t4))
end function kingdonFerland96
  

subroutine getCollisionalRates(thisIon, iTransition, temperature, excitation, deexcitation)
  real :: excitation, deexcitation
  type(IONTYPE) :: thisIon
  integer :: iTransition, i
  real :: temperature
  real :: thisGamma
  real :: t , fac
  real :: boltzFac

  t = max(min(real(thisIon%transition(iTransition)%t(thisIon%transition(iTransition)%ngamma)),temperature), &
       real(thisIon%transition(iTransition)%t(1)))
  call locate(thisIon%transition(iTransition)%t, thisIon%transition(iTransition)%ngamma, t, i)
  fac = (t - thisIon%transition(iTransition)%t(i))/(thisIon%transition(iTransition)%t(i+1) - thisIon%transition(iTransition)%t(i))
  thisGamma = thisIon%transition(iTransition)%gamma(i) + &
       fac * (thisIon%transition(iTransition)%gamma(i+1) - thisIon%transition(iTransition)%gamma(i))

  boltzFac =  exp(-thisIon%transition(iTransition)%energy / (kev*temperature))

  fac = (8.63e-6 / sqrt(temperature)) * thisGamma

  deexcitation =  fac / thisIon%level(thisIon%transition(iTransition)%j)%g

  excitation =  fac / thisIon%level(thisIon%transition(iTransition)%i)%g * boltzFac

!  excitation = fac * boltzFac
!  deexcitation = fac * thisIon%level(thisIon%transition(iTransition)%j)%g &
!       / thisIon%level(thisIon%transition(iTransition)%i)%g

!  if (thisIon%species == "O III") then
!     do i = 1, thisIon%transition(1)%ngamma
!        write(*,*) i,thisIon%transition(1)%t(i),thisIon%transition(1)%gamma(i)
!     enddo
!     if (iTransition ==1) then
!        write(*,*) thisGamma, boltzFac, fac, excitation, deexcitation
!     endif
!  endif


end subroutine getCollisionalRates


subroutine getForbiddenLineLuminosity(grid, species, wavelength, luminosity)
  type(GRIDTYPE) :: grid
  character(len=*) :: species
  real(double) :: wavelength
  real(double) :: fac
  real(double) :: luminosity
  integer :: iIon, iTrans, i

  iTrans = 0
  iIon = returnIonNumber(species, grid%ion, grid%nIon) 
  do i = 1, grid%ion(iIon)%nTransitions
     fac = grid%ion(iIon)%transition(i)%lambda
     fac = fac - wavelength
     fac = abs(fac)
     fac = fac / wavelength
     if (fac  < 0.001d0) then
        iTrans = i
        exit
     endif
  enddo
  if (iTrans == 0) then
     write(*,*) "No transition found at ",wavelength, "Angstroms"
     stop
  endif
  luminosity = 0.d0
  call sumLineLuminosity(grid%octreeroot, luminosity, iIon, iTrans, grid)
end subroutine getForbiddenLineLuminosity

recursive subroutine sumLineLuminosity(thisOctal, luminosity, iIon, iTrans, grid)
  type(GRIDTYPE) :: grid
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i, iIon, iTrans
  real(double) :: luminosity, v, rate
  real :: pops(10)
  pops = 0.
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call sumLineLuminosity(child, luminosity, iIon, iTrans, grid)
                exit
             end if
          end do
       else
          v = cellVolume(thisOctal, subcell)

          call solvePops(grid%ion(iIon), pops, thisOctal%ne(subcell), thisOctal%temperature(subcell), &
               thisOctal%ionFrac(subcell,iion),thisOctal%nh(subcell))
          rate =  pops(grid%ion(iion)%transition(iTrans)%j) * grid%ion(iion)%transition(itrans)%energy * &
               grid%ion(iion)%transition(itrans)%a/ergtoev
          rate = rate * grid%ion(iion)%abundance * thisOctal%nh(subcell) * thisOctal%ionFrac(subcell, iion)
          luminosity = luminosity + rate * v * 1.d30
          
       endif
    enddo
  end subroutine sumLineLuminosity


subroutine metalcoolingRate(ionArray, nIons, thisOctal, subcell, nh, ne, temperature, total, debug)
  type(IONTYPE) :: ionArray(*)
  integer :: nIons, subcell
  type(OCTAL) :: thisOctal
  real(double) :: ne, nh
  real :: temperature
  real(double) :: rate, total
  real :: pops(10)
  integer :: i, j
  logical, optional :: debug
  pops = 0.
  total = 0.d0
  do j = 5, nIons
     if (ionArray(j)%nTransitions > 0) then
        call solvePops(ionArray(j), pops, ne, temperature, thisOctal%ionFrac(subcell,j), thisOctal%nh(subcell))
        rate = 0.d0
        do i = 1, ionArray(j)%nTransitions
           rate = rate + pops(ionArray(j)%transition(i)%j)*ionArray(j)%transition(i)%energy*ionArray(j)%transition(i)%a/ergtoev
        enddo
        rate = rate * ionArray(j)%abundance * nh * thisOctal%ionFrac(subcell, j)
        if (rate < 0.) then
           write(100,'(a20,a,a,1p,e12.4,a,0p, f10.1,a,1pe12.4)') "Negative contribution from ", &
                trim(ionArray(j)%species),":",rate," at T = ",temperature, &
                " ion frac ",thisOctal%ionFrac(subcell,j)
        endif
        if (present(debug)) then
           if (debug) then
                 write(100,'(a20,a,a,1p,e12.4,a,0p, f10.1,a,1pe12.4)') "Contribution from ", &
                      trim(ionArray(j)%species),":",rate," at T = ",temperature, &
                      " ion frac ",thisOctal%ionFrac(subcell,j)
           endif
        endif
     total = total + rate
     endif
  enddo
end subroutine metalcoolingRate
  

subroutine solvePops(thisIon, pops, ne, temperature, ionFrac, nh, debug)
  type(IONTYPE) :: thisIon
  real(double) :: ne
  real :: pops(*)
  real :: temperature
  real(double), allocatable :: matrixA(:,:), MatrixB(:), tempMatrix(:,:), qeff(:,:),  rates(:)
  integer :: n, iTrans, i, j
  real :: excitation, deexcitation, arateji
  real(double) :: nh, ionFrac
  logical :: ok
  logical, optional :: debug
  excitation  = 0.; deexcitation = 0.
  n = thisIon%nLevels
  allocate(matrixA(1:n, 1:n), matrixB(1:n), tempMatrix(1:n,1:n), qeff(1:n,1:n), rates(1:n))

  matrixA = 1.d-30
  matrixB = 0.d0


  call getRecombs(rates, thision, dble(temperature), ne, ionFrac, nh)

  matrixA(1,:) = 1.d0
  matrixB(1) = 1.d0

  matrixB(2:n) = 0.d0 !rates(2:n) * ne * (nh * ionFrac * thisIon%abundance)
  do iTrans = 1, thisIon%nTransitions
     i = thision%transition(itrans)%i
     j = thision%Transition(itrans)%j
     call getCollisionalRates(thisIon, iTrans, temperature, excitation, deexcitation)
     qeff(i,j) = excitation
     qeff(j,i) = deexcitation
  enddo

  do i = 2, n
     do j = 1, n

        do iTrans = 1, thisIon%nTransitions
           if (((i == thision%transition(itrans)%i).and.(j == thision%Transition(itrans)%j)).or. &
               ((i == thision%transition(itrans)%j).and.(j == thision%Transition(itrans)%i))) then
              arateji =  thisIon%transition(iTrans)%a
              
              matrixA(i,j) = matrixA(i,j) + ne * qeff(j, i)
              matrixA(i,i) = matrixA(i,i) - ne * qeff(i, j)
              if (j > i) then
                 matrixA(i,j) = matrixA(i,j) + arateji
              else
                 matrixA(i,i) = matrixA(i,i) - arateji
              endif

           endif
        enddo

     enddo
  enddo

  tempMatrix = matrixA

  if (PRESENT(debug)) then
     if (debug) then
        do i = 1, n
           write(*,'(1p,9e12.3)') tempmatrix(i,1:n)
        enddo
     endif
  endif
  
  call luSlv(matrixA, matrixB)

  matrixB(1:n) = matrixB(1:n) / SUM(matrixB(1:n))
  
  ok = .true.

  if (.not.ok) then
     write(*,*) "Population solver failed for: ",thisIon%species
     write(*,*) matrixB(1:n)
     write(*,*) "nlevels",thisIon%nLevels,"ntrans",thisIon%nTransitions 
     write(*,*) "temp",temperature,"ne",ne
     
     do i = 1, n
        write(*,'(1p,9e12.3)') tempmatrix(i,1:n)
     enddo
     matrixB = 0.d0
     matrixB(1) = 1.d0
     write(*,*) "Setting pops to ground state"
  endif

  do i = 1, n
     pops(i) = max(1.d-30,matrixB(i))
  enddo

  deallocate(matrixA, matrixB, tempMatrix, qeff, rates)

end subroutine solvePops

subroutine calcHeRecombs(te, alpha1, alpha21s, alpha21p, alpha23s)
  real :: te, t
  real(double) :: alpha1, alpha21s, alpha21p, alpha23s

  t = te / 1.e4

  alpha1 = 1.54e-13 * t**(-0.486)

  alpha21s = 2.06e-14 * t**(-0.676)

  alpha21p = 4.17e-14 * t**(-0.861)
  
  alpha23s = 2.10e-13 * t**(-0.778)
end subroutine calcHeRecombs


subroutine createSahaMilneTables(hTable, heTable)
  type(SAHAMILNETABLE) :: hTable, heTable
  real(double) :: nu0_h, nu0_he, nufinal_h, nufinal_he
  integer :: nFreq, nTemp
  integer :: i, j
  real :: e
  real(double) :: t1, t2
  real :: hxsec, hexsec
  real(double) :: dfreq, jnu

  nFreq = 10000
  nTemp = 100
  hTable%nFreq = nFreq
  heTable%nFreq = nFreq
  hTable%nTemp = nTemp
  heTable%nTemp = nTemp

  allocate(hTable%freq(1:nFreq), heTable%freq(1:nFreq))
  allocate(hTable%emissivity(1:nTemp), heTable%emissivity(1:ntemp))
  allocate(hTable%temp(1:ntemp), heTable%temp(1:ntemp))
  allocate(hTable%Clyc(1:nTemp,1:nFreq))
  allocate(heTable%Clyc(1:nTemp,1:nFreq))

  nu0_h = 13.6d0/ergtoev/hcgs
  nu0_he = 24.59d0/ergtoev/hcgs

  nufinal_h = 2.d0*nu0_h
  nufinal_he = 2.d0*nu0_he

  do i = 1, nFreq
     hTable%freq(i) = log10(nu0_h) + (log10(nuFinal_h)-log10(nu0_h))*dble(i-1)/dble(nFreq-1)
     heTable%freq(i) = log10(nu0_he) + (log10(nuFinal_he)-log10(nu0_he))*dble(i-1)/dble(nFreq-1)
  enddo
  hTable%freq(1:hTable%nFreq) = 10.d0**hTable%freq(1:hTable%nfreq)
  heTable%freq(1:hTable%nFreq) = 10.d0**heTable%freq(1:hTable%nfreq)

  t1 = 5000.d0
  t2 = 20000.d0

  do i = 1, nTemp
     hTable%temp(i) = t1 + (t2-t1)*dble(i-1)/dble(nTemp-1)
     heTable%temp(i) = t1 + (t2-t1)*dble(i-1)/dble(nTemp-1)
  enddo

  do i = 1, nTemp
     hTable%Clyc(i,1) = 0.d0
     heTable%Clyc(i,1) = 0.d0
     do j = 2, nFreq


        e = hTable%freq(j) * hcgs* ergtoev
        call phfit2(1, 1, 1 , e , hxsec)

        dFreq = hTable%freq(j)-hTable%freq(j-1)
        jnu = ((hcgs*hTable%freq(j)**3)/(cSpeed**2)) * ((hcgs**2) /(twoPi*mElectron*Kerg*hTable%temp(i)))**(1.5d0) * &
             dble(hxsec/1.d10) *  exp(-hcgs*(hTable%freq(j)-nu0_h)/(kerg*hTable%temp(i)))
        hTable%Clyc(i,j) = hTable%Clyc(i,j-1) + jnu * dfreq


        e = heTable%freq(j) * hcgs * ergtoev
        call phfit2(2, 2, 1 , e , hexsec)

        dFreq = heTable%freq(j)-heTable%freq(j-1)
        jnu = 2.d0*((hcgs*heTable%freq(j)**3)/(cSpeed**2)) * ((hCgs**2) / (twoPi*mElectron*kerg*heTable%temp(i)))**(1.5d0) * &
             dble(hexsec/1.d10) * exp(-hcgs*(heTable%freq(j)-nu0_he)/(kerg*heTable%temp(i)))
        heTable%Clyc(i,j) = heTable%Clyc(i,j-1) + jnu * dfreq
     enddo
     hTable%emissivity(i) = hTable%Clyc(i,hTable%nFreq)
     heTable%emissivity(i) = heTable%Clyc(i,heTable%nFreq)
     hTable%Clyc(i,1:hTable%nFreq) = hTable%Clyc(i,1:hTable%nFreq) / hTable%Clyc(i,hTable%nFreq)
     heTable%Clyc(i,1:heTable%nFreq) = heTable%Clyc(i,1:heTable%nFreq) / heTable%Clyc(i,heTable%nFreq)
  end do

!  do j = 1, nFreq
!     write(99  ,*) htable%freq(j),htable%clyc(50,j)
!  enddo


end subroutine createSahaMilneTables

subroutine getSahaMilneFreq(table,temperature, thisFreq)
  type(SAHAMILNETABLE) :: table
  real(double) :: temperature, thisfreq, r, t, fac
  integer :: i, j

  t = max(5000.d0, min(20000.d0, temperature))
  call locate(table%temp, table%nTemp, t, i)
  call random_number(r)
  call locate(table%Clyc(i,1:table%nfreq), table%nFreq, r, j)
  fac = (r - table%Clyc(i,j))/(table%Clyc(i,j+1)-table%cLyc(i,j))
  thisFreq = table%freq(j) + fac * (table%freq(j+1)-table%freq(j))
end subroutine getSahaMilneFreq

subroutine twoPhotonContinuum(thisFreq)

! based on table ii of drake, victor, dalgarno, 1969, PhyRev Vol 180, pg 25

  real(double) :: thisFreq
  real :: y(21) = (/ 0., 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.275, 0.300, &
       0.325, 0.350, 0.375, 0.400, 0.425, 0.450, 0.475, 0.500 /)
  real :: hei(21) = (/ 0., 7.77e0, 2.52e1, 4.35e1, 5.99e1, 7.42e1, 8.64e1, 9.69e1, 1.06e2, 1.13e2, 1.20e2, 1.25e2, &
       1.30e2, 1.34e2, 1.37e2, 1.40e2, 1.42e2, 1.43e2, 1.45e2, 1.45e2, 1.45e2 /)
  real :: freq = 3.86e15, fac, r
  real :: prob(21)
  integer :: i
  prob(1) = 0.
  do i = 2, 21
     prob(i) = prob(i-1) + (y(i)-y(i-1)) * hei(i)
  enddo
  prob(1:21) = prob(1:21)/prob(21)
  thisFreq = 0.
  do while((thisFreq*hcgs*ergtoev) < 13.6)
     call random_number(r)
     call locate(prob, 21, r, i)
     fac = y(i) + ((r - prob(i))/(prob(i+1)-prob(i)))*(y(i+1)-y(i))
     thisFreq = (1.-fac)*freq
  enddo
end subroutine twoPhotonContinuum

subroutine getRecombs(rates, thision, temperature, ne, ionFrac, nh)  
  real(double) :: rates(:)
  type(IONTYPE) :: thisIon
  real(double) :: temperature, ne, ionfrac, nh

  select case(thisIon%species)
     case("O II")
        rates(1) = 0.d0
        rates(2) = rateFit(7.218d0, -0.575d0, temperature)
        rates(3) = rateFit(4.812d0, -0.575d0, temperature)
        rates(4) = rateFit(3.581d0, -0.495d0, temperature)
        rates(5) = rateFit(1.790d0, -0.495d0, temperature)
     case DEFAULT
        rates = 0.d0
   end select
 end subroutine getRecombs


subroutine getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDeltaT)
  type(GRIDTYPE) :: grid
  type(OCTAL) :: thisOctal
  integer :: subcell
  real(double) :: hHeating, heHeating, totalHeating, v, epsOverDeltaT, dustHeating

  v = cellVolume(thisOctal, subcell)
  Hheating= thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,1) * grid%ion(1)%abundance &
       * (epsOverDeltaT / (v * 1.d30))*thisOctal%Hheating(subcell) ! equation 21 of kenny's
  Heheating= thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,3) * grid%ion(3)%abundance &
       * (epsOverDeltaT / (v * 1.d30))*thisOctal%Heheating(subcell) ! equation 21 of kenny's
  dustHeating = (epsOverDeltaT / (v * 1.d30))*thisOctal%distanceGrid(subcell) ! equation 14 of Lucy 1999
  totalHeating = (Hheating + HeHeating + dustHeating)
  
end subroutine getHeating

function rateFit(a,b,t) result (rate)
  real(double) :: a, b, t, rate

  rate = 1.d-13 * a * (t/1.d4)**b
end function rateFit

subroutine createRecombTable(table, tablefilename)
  type(RECOMBTABLE) :: table
  character(len=*) :: tablefilename
  character(len=200) :: filename, datadirectory
  integer :: ia, ib, ne, ncut, i
  real :: e(1000)

  dataDirectory = " "
  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//tablefilename

  open(20,file=filename,status="old",form="formatted")
  read(20,*) table%nTemp, table%nrho

  allocate(table%rho(1:table%nRho))
  allocate(table%temp(1:table%ntemp))
  allocate(table%emissivity(1:table%ntemp,1:table%nrho))

  do ia=1,table%ntemp
     do ib=1,table%nrho
        read(20,'(1x,e10.3,5x,e10.3,13x,i2)') table%rho(ib),table%temp(ia),ncut
        ne=ncut*(ncut-1)/2
        read(20,'(8e10.3)') e(1:ne)
        table%emissivity(ia,ib) = SUM(e(1:ne-1))
     enddo
  enddo
  close(20)
end subroutine createRecombTable



  
subroutine createGammaTable(table, thisfilename)

! Ferland 1980 PASP 92 596

  type(GAMMATABLE) :: table
  character(len=*) :: thisfilename
  character(len=200) :: dataDirectory, filename
  integer :: i
  dataDirectory = " "

  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//thisfilename

  open(40, file=filename, form="formatted", status="old")
  read(40,*) table%nTemp, table%nFreq

  allocate(table%freq(1:table%nFreq))
  allocate(table%temp(1:table%nTemp))
  allocate(table%gamma(table%nFreq, table%nTemp))

  table%temp = 0.d0
  read(40,*)  table%temp(1:table%nTemp)
  do i = 1, table%nFreq 
     read(40,*) table%freq(i), table%gamma(i,1:table%nTemp)
  enddo
  close(40)

  do i = 1, table%nFreq
     table%freq(i) = table%freq(i) * nuHydrogen
  enddo

  where (table%gamma(1:table%nFreq,1:table%nTemp) == 0.d0)
     table%gamma(1:table%nFreq,1:table%nTemp) = 1.d-30
  end where
  table%gamma(1:table%nFreq,1:table%nTemp) =log10(table%gamma(1:table%nFreq,1:table%nTemp))
  table%freq(1:table%nFreq) = log10(table%freq(1:table%nFreq))
  table%temp(1:table%nTemp) = log10(table%temp(1:table%nTemp))
end subroutine createGammaTable

real(double) function returnGamma(table, temp, freq) result(out)
  type(GAMMATABLE) :: table
  real(double) :: temp , freq
  integer :: i, j
  real(double) :: tfac, ffac, gamma1, gamma2, logT, logF

  logT = log10(temp)
  logF = log10(freq)

  call locate(table%temp, table%nTemp, logT, i)
  call locate(table%freq, table%nFreq, logF, j)

  if (logF >= table%freq(1).and.logF <= table%freq(table%nFreq)) then
     tfac  = (logT - table%temp(i))/(table%temp(i+1) - table%temp(i))
     ffac  = (logF - table%freq(j))/(table%freq(j+1) - table%freq(j))
     
     gamma1 = table%gamma(j,i) + tfac*(table%gamma(j, i+1) - table%gamma(j,i))
     gamma2 = table%gamma(j+1,i) + tfac*(table%gamma(j+1, i+1) - table%gamma(j+1,i))
     out = gamma1 + ffac*(gamma2 - gamma1)
     out = 10.d0**out
  else
     out = tiny(out)
  endif


end function returnGamma

subroutine addLymanContinua(nFreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
  type(GRIDTYPE) :: grid
  TYPE(OCTAL) :: thisOctal
  integer :: subcell
  integer :: nFreq
  integer :: i, k, iIon, n1, n2
  real(double) :: freq(:), spectrum(:), dfreq(:)
  real(double) :: jnu
  real :: e, hxsec
  real(double), parameter :: statisticalWeight(3) = (/ 2.d0, 0.5d0, 2.d0 /)


  ! do Saha-Milne continua for H, HeI and HeII

  do k = 1 , 2 !3!!!!!!!!!!!!!!!!!

     if (k == 1) iIon = 1
     if (k == 2) iIon = 3
     if (k == 3) iIon = 4

     call locate(freq, nfreq, grid%ion(iIon)%nuThresh, n1)
     n2 = nFreq
     do i = n1, n2
        
        e = freq(i) * hcgs* ergtoev


        call phfit2(grid%ion(iIon)%z, grid%ion(iIon)%n, grid%ion(iIon)%outerShell , e , hxsec)

        jnu = tiny(jnu)

        if (hxSec > 0.) then
           if (thisOctal%temperature(subcell) > 100.) then
              jnu = statisticalWeight(k) * ((hcgs*freq(i)**3)/(cSpeed**2)) * &
                   ((hcgs**2) /(twoPi*mElectron*Kerg*thisOctal%temperature(subcell)))**(1.5d0) * &
                   dble(hxsec/1.d10) *  &
                   exp(-hcgs*(freq(i)-grid%ion(iIon)%nuThresh)/(kerg*thisOctal%temperature(subcell)))
              jnu = jnu * thisOctal%ne(subcell) *(thisOctal%nh(subcell) * &
                   thisOctal%ionFrac(subcell,iIon+1) * grid%ion(iIon)%abundance)
           else
              jnu = tiny(jnu)
           endif
        endif

        spectrum(i) = spectrum(i) + jnu * dFreq(i) * fourPi   

     enddo
  enddo

end subroutine addLymanContinua


subroutine addHigherContinua(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, table)

! Ferland 1980 PASP 92 596

  type(GAMMATABLE) :: table(3)
  integer :: nFreq
  real(double) :: spectrum(:), freq(:), dfreq(:)
  type(OCTAL) :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  integer :: i, k, iEnd, iIon
  real(double) :: fac

  do k = 1 , 3

     if (k == 1) iIon = 1
     if (k == 2) iIon = 3
     if (k == 3) iIon = 4

     call locate(freq, nFreq, grid%ion(iIon)%nuThresh, iEnd)
     do i = 1, iEnd
        fac = returnGamma(table(k), dble(thisOctal%temperature(subcell)) , freq(i))
        fac = fac*1.d-40 ! units of 10^-40 erg/s/cm/cm/cm/hz
        fac = fac * thisOctal%ne(subcell) * thisOctal%nh(subcell) &
             * thisOctal%ionFrac(subcell,iIon+1) * grid%ion(iIon)%abundance
        spectrum(i) = spectrum(i) + fac*dfreq(i)
     enddo
  enddo

end subroutine addHigherContinua

subroutine addFreeFreeContinua(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)

  integer :: nFreq
  real(double) :: spectrum(:), freq(:), dfreq(:)
  type(OCTAL) :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  integer :: i
  real(double) :: fac, gauntFF, z
  real(double), allocatable :: g(:), xlf(:)
  integer :: iflag

  allocate(g(1:nFreq), xlf(1:nFreq))

  xlf = log10(hCgs*freq*ergtoev/rydbergToEv)
!hydrogen

  z = 1.d0
  call getGauntFF(z, log10(thisOctal%temperature(subcell)), xlf, g, iflag)

  do i = 1, nFreq

     fac = fourPi * 54.43 * (z**2) * g(i) * &
          exp(-freq(i)*hCgs/(kErg*thisOctal%temperature(subcell)))/sqrt(thisOctal%temperature(subcell)) * 1.d-40
     spectrum(i) = spectrum(i) + fac*dfreq(i) 
  enddo
  deallocate(g, xlf)
end subroutine addFreeFreeContinua

    ! this subroutine computes the gaunt factors for any charge
    ! it generates thermally averaged free-free non-relativistic gaunt
    ! factor for a hydrogenic ion of charge z, with a maximum relative
    ! error of 0.007, (rms fitting error = 0.001) for tmps and freqs in
    ! intervals:
    !          10^-4 <= U <= 10^1.5
    !          10^-3 <= Gams <= 10^3~
    ! where U = h*nu/k*T and gams = Z^2 * Ryd / k*T. To obtain the stated
    ! accuracy the full number of significant figures must be retained.
    !
    ! this subroutine uses a two-dimensional chebyshev epansion computed 
    ! from expressions given bu Karzas and Latter (ApJSuppl., V.6, P.167, 
    ! 1961) augmented by various limiting forms of energy-specific gaunt-
    ! factors.
    ! D.G.Hummer, Jila, May 1987. ApJ 327, 477
    ! modified with correct limits, J Ferguson, July 94
    subroutine getGauntFF(z, log10Te, xlf, g, iflag)
        implicit none

        integer, intent(out) :: iflag                ! explanation given in each area        

        real, intent(in) :: log10Te                  ! log10(Te)
        real(double), intent(in) ::  z                       ! log10 of nuclear charge
        real(double), dimension(:), intent(in) :: xlf        ! array of log10(h*nu/Ryd)
        real(double), dimension(size(xlf)), intent(out) :: g ! array of values of g

        ! local variables 
        integer ::  i, ir, j

        real, dimension(11) :: b
        real, dimension(8) :: c
        real :: con
        real(kind=8), dimension(88) :: dd
        real(kind=8), dimension(8, 11) :: d
        real :: gamma2                               ! gamma^2 (gams = Z^2 * Ryd / k*T)
        real :: slope                                ! slope
        real :: txg                                  ! hummer variable related to gamma^2
        real :: txu                                  ! hummer variable related to U
        real :: u                                    ! U = h*nu/k*T
        real :: xlrkt                                ! log(ryd/kt)

        dd = (/8.986940175e+00, -4.009515855e+00,  8.808871266e-01,& 
&          2.640245111e-02, -4.580645915e-02, -3.568055702e-03,&
&      	  2.827798067e-03,  3.365860195e-04, -8.006936989e-01,&
&      	  9.466021705e-01,  9.043402532e-02, -9.608451450e-02,&
&         -1.885629865e-02,  1.050313890e-02,  2.800889961e-03,&
&         -1.078209202e-03, -3.781305103e-01,  1.102726332e-01,&
&         -1.543619180e-02,  8.310561114e-03,  2.179620525e-02,&
&          4.259726289e-03, -4.181588794e-03, -1.770208330e-03,&
&          1.877213132e-02, -1.004885705e-01, -5.483366378e-02,& 
&         -4.520154409e-03,  8.366530426e-03,  3.700273930e-03,&
&          6.889320423e-04,  9.460313195e-05,  7.300158392e-02,&
&          3.576785497e-03, -4.545307025e-03, -1.017965604e-02,&
&         -9.530211924e-03, -3.450186162e-03,  1.040482914e-03,&
&          1.407073544e-03, -1.744671550e-03,  2.864013856e-02,&
&          1.903394837e-02,  7.091074494e-03, -9.668371391e-04,&
&         -2.999107465e-03, -1.820642230e-03, -3.874082085e-04,&
&         -1.707268366e-02, -4.694254776e-03,  1.311691517e-03,&
&          5.316703136e-03,  5.178193095e-03,  2.451228935e-03,&
&         -2.277321615e-05, -8.182359057e-04,  2.567331664e-04,&
&         -9.155339970e-03, -6.997479192e-03, -3.571518641e-03,&
&         -2.096101038e-04,  1.553822487e-03,  1.509584686e-03,&
&          6.212627837e-04,  4.098322531e-03,  1.635218463e-03,&
&         -5.918883504e-04, -2.333091048e-03, -2.484138313e-03,&
&         -1.359996060e-03, -5.371426147e-05,  5.553549563e-04,& 
&          3.837562402e-05,  2.938325230e-03,  2.393747064e-03,&
&          1.328839809e-03,  9.135013312e-05, -7.137252303e-04,&
&         -7.656848158e-04, -3.504683798e-04, -8.491991820e-04,&
&         -3.615327726e-04,  3.148015257e-04,  8.909207650e-04,&
&          9.869737522e-04,  6.134671184e-04,  1.068883394e-04,&
&         -2.046080100e-04/)
      
        d = reshape( dd, (/8, 11/) )

        ! compute temperature dependent coeffcients for U expansion

        ! xlrxt is log(ryd/kt), code note valid for Te > 10^8
        ! xlrxt is log(ryd/kt), code note valid for Te < 10^2.5
        if ( log10Te < 2.5 ) then
            xlrkt = 5.1983649 - 2.5
        else if ( log10Te > 8. ) then
            xlrkt = 5.1983649 - 8.
        else
            xlrkt = 5.1983649 - log10Te
        end if

        ! set txg
        txg = 0.66666667*(2.0*z+xlrkt)
        gamma2 = 10**(txg*1.5)

        con = 0.72727273*xlrkt+0.90909091
        do j=1,8
            ir = 9
            b(11) = d(j,11)
            b(10) = txg*b(11)+d(j,10)
            do i=1,9
                b(ir) = txg*b(ir+1)-b(ir+2)+d(j,ir)
                ir = ir-1
            end do
            c(j) = 0.25*(b(1)-b(3))
        end do
        
        ! sum U expansion
        ! loop through energy at fixed temperature
        do i = 1, size(xlf)
            txu = 0.72727273*xlf(i)+con
            u = 10**((txu - .90909091)/.72727273)
            ! criteria set by hummer limits. it is a wedge from
            ! log(hnu),log(T)) =
            ! (-5,2.5) to (-5,4) to (-1,8) to (4,8) to (-1.5,2.5).
            ! these limits correspond to the gamma^2 and U limits 
            ! given above
            if(abs(txu)<=2.0) then
                ir = 6
                b(8) = c(8)
                b(7) = txu*b(8)+c(7)
                do j=1,6
                    b(ir) = txu*b(ir+1)-b(ir+2)+c(ir)
                    ir = ir-1
                end do
                g(i) = b(1)-b(3)
                if( (log10Te>=2.5) .and. (log10Te<=8.0)) iflag = 0
                ! On the bottom side of the hummer box,u<-4 and gamma2>.33
            else if( (log10(u)<-4.0) .and. (gamma2<0.3) ) then
               g(i) = 0.551329 * log(2.24593/u)
              if( (log10Te>=2.5) .and. (log10Te<=8.0) ) iflag = 2
            ! On the bottom side of the box,u<-4 and gamma2<.33
            else if( (log10(u)<-4.0) .and. (gamma2>=0.3) ) then
                g(i) = 0.551329 * alog( 0.944931/u/sqrt(gamma2) )
                if( (log10Te>=2.5) .and.( log10Te<=8.0) ) iflag = 3
                ! Now on the bottom side of the box
                else if (txu>2.0) then
                ! Top of the box high T first
                if (log10(gamma2)<-3.) then
                    g(i) =  sqrt(0.9549297/u)
                    if(log10Te>=2.5.and.log10Te<=8.0) iflag = 4
                    ! Must interpolate between two asymptotes
                else if(log10(gamma2)<0.0) then
                    slope = ( sqrt(12*gamma2/u) - sqrt(0.9549297/u) ) / 3.0
                    g(i) =  sqrt(12*gamma2/u) + slope * log10(gamma2)
                    if((log10Te>=2.5) .and. (log10Te<=8.0)) iflag = 7
                else
                ! The top side with wedge of 1.05 where region 6 fails
                    g(i) = sqrt(12. * gamma2 / u)
                    g(i) = min(1.05,g(i))
                    if( (g(i)==1.05) .and. (log10Te>=2.5) .and. (log10Te<=8.0)) &
                         & iflag = 5
                    if( (g(i)<1.05) .and. (log10Te>=2.5) .and. (log10Te<=8.0)) &
                         & iflag = 6
                end if
            end if
            u = 0.0
        end do
    end subroutine getGauntFF



subroutine addHydrogenRecombinationLines(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)


  integer :: nFreq
  real(double) :: spectrum(:), freq(:), dfreq(:)
  type(OCTAL) :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  real :: emissivity(3:15,2:8)
  integer :: ilow, iup
  real(double) :: lymanalpha, hBeta, LineEmissivity, lineFreq, energy
  integer :: i

  real(double) :: lambdaTrans(20,20) = reshape( source=&
    (/ 000000.d-8,1215.67d-8,1025.72d-8,992.537d-8,949.743d-8,937.803d-8,930.748d-8,926.226d-8,923.150d-8,920.963d00,&
       919.352d-8,918.129d-8,917.181d-8,916.429d-8,915.824d-8,915.329d-8,914.919d-8,914.576d-8,914.286d-8,914.039d-8,&  
       0.00000d-8,0000000d-8,6562.80d-8,4861.32d-8,4340.36d-8,4101.73d-8,3970.07d-8,3889.05d-8,3835.38d-8,3797.90d-8,&
       3770.63d-8,3750.15d-8,3734.37d-8,3721.94d-8,3711.97d-8,3703.85d-8,3697.15d-8,3691.55d-8,3686.83d-8,3682.81d-8,&  
       0.00000d-8,0.00000d-8,0000000d-8,18751.0d-8,12818.1d-8,10938.1d-8,10049.4d-8,9545.98d-8,9229.02d-8,9014.91d-8,&
       8862.79d-8,8750.47d-8,8665.02d-8,8598.39d-8,8545.39d-8,8502.49d-8,8467.26d-8,8437.96d-8,8413.32d-8,8392.40d-8,&  
       0.00000d-8,0.00000d-8,0.00000d-8,0000000d-8,40512.0d-8,26252.0d-8,21655.0d-8,19445.6d-8,18174.1d-8,17362.1d-8,&
       16806.5d-8,16407.2d-8,16109.3d-8,15880.5d-8,15700.7d-8,15556.5d-8,15438.9d-8,15341.8d-8,15260.6d-8,15191.8d-8,&  
       0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0000000d-8,74578.0d-8,46525.0d-8,37395.0d-8,32961.0d-8,30384.0d-8,&
       28722.0d-8,27575.0d-8,26744.0d-8,26119.0d-8,25636.0d-8,25254.0d-8,24946.0d-8,24693.0d-8,24483.0d-8,24307.0d-8,&  
       0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,123680.d-8,75005.0d-8,59066.0d-8,51273.0d-8,&
       46712.0d-8,43753.0d-8,41697.0d-8,40198.0d-8,39065.0d-8,38184.0d-8,37484.0d-8,36916.0d-8,36449.0d-8,36060.0d-8,&  
       0.00000d-8,0000000d-8,0000000d-8,0000000d-8,0000000d-8,0.00000d-8,0000000d-8,190570.d-8,113060.d-8,87577.0d-8,&
       75061.0d-8,67701.0d-8,62902.0d-8,59552.0d-8,57099.0d-8,55237.0d-8,53783.0d-8,52622.5d-8,51679.0d-8,50899.0d-8,&  
       0.00000d-8,0000000d-8,0000000d-8,0000000d-8,0000000d-8,0.00000d-8,0000000d-8,0000000d-8,277960.d-8,162050.d-8,&
       123840.d-8,105010.d-8,93894.0d-8,86621.0d-8,81527.0d-8,77782.0d-8,74930.0d-8,72696.0d-8,70908.0d-8,69448.0d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,388590.d-8,&  
       223340.d-8,168760.d-8,141790.d-8,125840.d-8,115360.d-8,108010.d-8,102580.d-8,98443.0d-8,95191.0d-8,92579.0d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       525200.d-8,298310.d-8,223250.d-8,186100.d-8,164070.d-8,149580.d-8,139380.d-8,131840.d-8,126080.d-8,121530.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,690500.d-8,388320.d-8,288230.d-8,238620.d-8,209150.d-8,189730.d-8,176030.d-8,165900.d-8,158120.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,887300.d-8,494740.d-8,364610.d-8,300020.d-8,261610.d-8,236260.d-8,218360.d-8,205090.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,111800.d-7,619000.d-8,453290.d-8,371000.d-8,322000.d-8,289640.d-8,266740.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,138600.d-7,762300.d-8,555200.d-8,452220.d-8,390880.d-8,350300.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,169400.d-7,926100.d-8,671200.d-8,544400.d-8,468760.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,204400.d-7,111200.d-7,802300.d-8,648200.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,243800.d-7,132100.d-7,949200.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,288200.d-7,155400.d-7,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,337400.d-7,&  
       0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /), shape=(/20,20/))


  emissivity(15, 2:8) = (/ 0.156E-01, 0.541E-02, 0.266E-02, 0.154E-02, 0.974E-03, 0.657E-03, 0.465E-03 /)
  emissivity(14, 2:8) = (/ 0.192E-01, 0.666E-02, 0.328E-02, 0.190E-02, 0.120E-02, 0.811E-03, 0.573E-03 /)
  emissivity(13, 2:8) = (/ 0.240E-01, 0.832E-02, 0.411E-02, 0.237E-02, 0.151E-02, 0.102E-02, 0.719E-03 /)
  emissivity(12, 2:8) = (/ 0.305E-01, 0.106E-01, 0.524E-02, 0.303E-02, 0.193E-02, 0.130E-02, 0.917E-03 /)
  emissivity(11, 2:8) = (/ 0.397E-01, 0.138E-01, 0.683E-02, 0.396E-02, 0.252E-02, 0.170E-02, 0.119E-02 /)
  emissivity(10, 2:8) = (/ 0.530E-01, 0.184E-01, 0.914E-02, 0.531E-02, 0.338E-02, 0.228E-02, 0.158E-02 /)
  emissivity(9, 2:8) = (/ 0.731E-01, 0.254E-01, 0.127E-01, 0.737E-02, 0.469E-02, 0.312E-02, 0.204E-02 /)
  emissivity(8,2:8) = (/ 0.105E+00, 0.366E-01, 0.183E-01, 0.107E-01, 0.673E-02, 0.425E-02, 0.000E+00 /)
  emissivity(7,2:8) = (/ 0.159E+00, 0.555E-01, 0.278E-01, 0.162E-01, 0.976E-02, 0.000E+00, 0.000E+00 /)
  emissivity(6,2:8) = (/ 0.259E+00, 0.904E-01, 0.455E-01, 0.255E-01, 0.000E+00, 0.000E+00, 0.000E+00 /)
  emissivity(5,2:8) = (/ 0.468E+00, 0.163E+00, 0.802E-01, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
  emissivity(4,2:8) = (/ 0.100E+01, 0.339E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
  emissivity(3,2:8) = (/ 0.286E+01, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)


  Hbeta = 10.d0**(-0.870d0*log10(thisOctal%temperature(subcell)) + 3.57d0)
  Hbeta = Hbeta * thisOctal%ne(subcell) *  thisOctal%nh(subcell) * &
       thisOctal%ionFrac(subcell, 2) * grid%ion(1)%abundance

  ! calculate emission due to HI recombination lines [e-25 ergs/s/cm^3]
  do iup = 15, 3, -1
     do ilow = 2, min0(8, iup-1)
        lineEmissivity = emissivity(iup, ilow) * Hbeta * 1.e-25
        energy = (hydE0eV / dble(iLow**2))-(hydE0eV / dble(iUp**2))
        lineFreq = cSpeed/lambdaTrans(iup, ilow)
!        if (iup ==4 .and. ilow == 2) then
!           write(*,*) "hbeta: ",lambdatrans(iup,ilow)*1.d8, emissivity(iup, ilow)
!           stop
!        endif
        if ((lineFreq > freq(1)).and.(lineFreq < freq(nfreq))) then
           call locate(freq, nFreq, lineFreq, i)
           i = i + 1
           spectrum(i) = spectrum(i) + lineEmissivity
        endif
     end do
  end do

  LymanAlpha = 10.d0**(-0.897*log10(thisOctal%temperature(subcell)) + 5.05d0) * &
       thisOctal%ne(subcell) *(thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,2) * &
       grid%ion(1)%abundance) * 1.d-25
  lineFreq = cSpeed/1215.67D-8
  if ((lineFreq > freq(1)).and.(lineFreq < freq(nfreq))) then
     call locate(freq, nFreq, lineFreq, i)
     i = i + 1
     spectrum(i) = spectrum(i) + lymanAlpha
  endif
  

end subroutine addHydrogenRecombinationLines

real(double) function getPhotonFreq(nfreq, freq, spectrum) result(Photonfreq)
  integer :: nFreq
  real(double) :: freq(:), spectrum(:), fac, r
  real(double), allocatable :: tSpec(:)
  integer :: i
  logical, save :: firstTime = .true.

  allocate(tSpec(1:nFreq))
  
  tSpec(1:nFreq) = spectrum(1:nFreq)

  do i = 2, nFreq
     tSpec(i) = tSpec(i) + tSpec(i-1)
  enddo
  tSpec(1:nFreq) = tSpec(1:nFreq) - tSpec(1)
  if (tSpec(nFreq) > 0.d0) then
     tSpec(1:nFreq) = tSpec(1:nFreq) / tSpec(nFreq)
     call random_number(r)
     call locate(tSpec, nFreq, r, i)
     fac = (r - tSpec(i)) / (tSpec(i+1)-tSpec(i))
     photonFreq = freq(i) + fac * (freq(i+1)-freq(i))
  else
     photonFreq = cSpeed / (1.d-8 * 100.e4)
  endif
  if (firstTime) then
     firstTime = .false.
     open(32, file="pdf.dat",form="formatted",status="unknown")
     do i = 1, nFreq
        write(32,*) freq(i), tspec(i)
     enddo
     close(32)
  endif
end function getPhotonFreq


subroutine addForbiddenLines(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)

  integer :: nFreq
  real(double) :: spectrum(:), freq(:), dfreq(:)
  type(OCTAL) :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  integer :: iIon, iTrans, j
  real :: pops(10)
  real(double) :: rate, lineFreq
  pops = 0.

  do iIon = 3, grid%nIon
     do iTrans = 1, grid%ion(iIon)%nTransitions
        call solvePops(grid%ion(iIon), pops, thisOctal%ne(subcell), thisOctal%temperature(subcell), &
             thisOctal%ionFrac(subcell,iion),thisOctal%nh(subcell))
        rate =  pops(grid%ion(iion)%transition(iTrans)%j) * grid%ion(iion)%transition(itrans)%energy * &
             grid%ion(iion)%transition(itrans)%a/ergtoev
        rate = rate * grid%ion(iion)%abundance * thisOctal%nh(subcell) * thisOctal%ionFrac(subcell, iion)
        lineFreq  =  (grid%ion(iion)%transition(itrans)%energy / ergtoEV) / hCgs
        if ((lineFreq > freq(1)).and.(lineFreq < freq(nfreq))) then
           call locate(freq, nFreq, lineFreq, j)
           j = j + 1
           spectrum(j) = spectrum(j) + rate
        endif
     enddo
  enddo

end subroutine addForbiddenLines
  

subroutine addHeRecombinationLines(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)

  integer :: nFreq
  real(double) :: spectrum(:), freq(:), dfreq(:)
  type(OCTAL) :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  integer :: i, j, k
  real :: fac, t, aj ,bj, cj
  real(double) :: lineFreq !, lambda
  real :: emissivity
!  real :: heII4686
!  integer :: ilow , iup
  integer,parameter :: nHeIILyman = 4
!  real(double) :: heIILyman(4)
  real(double) :: freqheIILyman(4) = (/ 3.839530, 3.749542, 3.555121, 2.99963 /)



  ! HeI lines 

  call locate(heIrecombinationNe, 3, real(log10(thisOctal%ne(subcell))), i)
  fac = (log10(thisOctal%ne(subcell)) - heIrecombinationNe(i))/(heIrecombinationNe(i+1)-heIrecombinationNe(i))

  do j = 1, 32
     aj = heIrecombinationFit(j,i,1) + fac*(heIrecombinationfit(j,i+1,1)-heIrecombinationfit(j,i,1))
     bj = heIrecombinationFit(j,i,2) + fac*(heIrecombinationfit(j,i+1,2)-heIrecombinationfit(j,i,2))
     cj = heIrecombinationFit(j,i,3) + fac*(heIrecombinationfit(j,i+1,3)-heIrecombinationfit(j,i,3))
     t = thisOctal%temperature(subcell)/1.e4
     emissivity = aj * (t**bj) * exp(cj / t) ! Benjamin et al. 1999 ApJ 514 307
     emissivity = emissivity * thisOctal%ne(subcell) * thisOctal%nh(subcell) * &
          thisOctal%ionFrac(subcell, 3) * grid%ion(3)%abundance

     lineFreq = cspeed / (heiRecombinationLambda(j)*1.e-8)
     call locate(freq, nFreq, lineFreq, k)
     k = k + 1
     spectrum(k) = spectrum(k) + emissivity
  enddo
!

!  HeII4686 = 10.d0**(-0.997d0*log10(thisOctal%temperature(subcell))+5.16d0)
!  HeII4686 = HeII4686*thisOctal%ne(subcell)*thisOctal%nh(subcell)*thisOctal%ionFrac(subcell,5)*grid%ion(4)%abundance
!  
!  ! calculate emission due to HeII recombination lines [e-25 ergs/s/cm^3]                                                     
!  do iup = 30, 3, -1
!     do ilow = 2, min0(16, iup-1)
!        emissivity= HeIIrecombinationLines(iup, ilow)*HeII4686*1.d-25
!
!
!        lambda = 227.838 / (1./real(ilow**2)  - 1./real(iup**2))!!!!!!!!!!!!!!!!!!!
!        lineFreq = cSpeed/(lambda * 1.d-8)
!
!     call locate(freq, nFreq, lineFreq, k)
!     spectrum(k) = spectrum(k) + emissivity
!     end do
!  end do
!

  ! He II Lyman series

!  heIILyman(1:4) = (/ 0.0334, 0.0682, 0.1849, 1. /)

!  ! calculate Lyman alpha first
!  HeIILyman(4) = 10.d0**(-0.792d0*log10(thisOctal%temperature(subcell))+6.01d0)
!  HeIILyman(4) = HeIILyman(4)*thisOctal%ne(subcell)*thisOctal%nh(subcell) * &
!       grid%ion(3)%abundance * thisOctal%ionFrac(subcell, 5) * 1.d-25
!
!  do i = 1, NHeIILyman-1
!     HeIILyman(i) = HeIILyman(i)*HeIILyman(4)
!  end do
!  do i = 1, nHeIILyman
!     lineFreq = freqHeIILyman(i) * nuHydrogen
!     call locate(freq, nFreq, lineFreq, k)
!     spectrum(k) = spectrum(k) + HeIILyman(i)
!  enddo



end subroutine addHeRecombinationLines

subroutine addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, nlambda, lamArray)

  integer :: nFreq
  real(double) :: spectrum(:), freq(:), dfreq(:)
  type(OCTAL), pointer :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  integer :: nLambda
  real :: lamArray(:)
  integer :: i, iLam
  real :: thisLam
  real(double), allocatable :: kabsArray(:)


  allocate(kAbsArray(1:nlambda))

  call returnKappa(grid, thisOctal, subcell, kappaAbsArray=kAbsArray)


  do i = 1, nFreq
     thisLam = (cSpeed / freq(i)) * 1.e8
     if ((thisLam >= lamArray(1)).and.(thisLam <= lamArray(nlambda))) then
        call hunt(lamArray, nLambda, real(thisLam), iLam)
        spectrum(i) = spectrum(i) + bnu(freq(i), dble(thisOctal%temperature(subcell))) * &
             kAbsArray(iLam) *1.d-10* dFreq(i) * fourPi
     endif
  enddo

  deallocate(kAbsArray)

end subroutine addDustContinuum








subroutine readHeIRecombinationLinesFit()
  character(len=200) :: filename, datadirectory
  integer :: i
  dataDirectory = " "
  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//"bss1999.dat"

  open(40,file=filename,status="old",form="formatted")
  do i = 1, 32
     read(40, *) HeIRecombinationLambda(i), heIRecombinationFit(i,1,1:3), &
           heIRecombinationFit(i,2,1:3),  heIRecombinationFit(i,3,1:3)
  enddo
  close(40)
!  HeiRecombinationFit(1:32, 1:3, 1:3) = log10(HeiRecombinationFit(1:32, 1:3, 1:3))
  heIrecombinationNe(1) = 2.
  heIrecombinationNe(2) = 4.
  heIrecombinationNe(3) = 6.
end subroutine readHeIRecombinationLinesFit


subroutine readHeIIrecombination()
  character(len=200) :: filename, datadirectory
  integer :: iup, ilow, i
  dataDirectory = " "

  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//"r2b0100.dat"

  open(40,file=filename,status="old",form="formatted")
  
  ! read in HeII recombination lines [e-25 ergs*cm^3/s]                                                                      
  ! (Storey and Hummer MNRAS 272(1995)41)                                                                                    
  open(40, file=filename, status = "old",form="formatted")
  do iup = 30, 3, -1
     read(40,*) (HeIIRecombinationLines(iup, ilow), ilow = 2, min0(16, iup-1))
  end do
  close(40)
end subroutine readHeIIrecombination

#ifdef MPI
  subroutine updateGridMPIphoto(grid)
    implicit none
    include 'mpif.h'
    type(gridtype) :: grid
    integer :: nOctals, nVoxels, i
    real, allocatable :: nCrossings(:)
    real, allocatable :: tempRealArray(:)
    real(double), allocatable :: hHeating(:), heHeating(:)
    real(double), allocatable :: photoIonCoeff(:,:)
    real(double), allocatable :: tempDoubleArray(:)
    real(double), allocatable :: distanceGrid(:)
    integer :: np, ierr,my_rank, nIndex

    ! FOR MPI IMPLEMENTATION=======================================================
    !  Get my process rank # 
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  
    ! Find the total # of precessor being used in this run
    call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
    nOctals = 0
    nVoxels = 0
    call countVoxels(grid%octreeRoot,nOctals,nVoxels)
    allocate(nCrossings(1:nVoxels))
    allocate(hHeating(1:nVoxels))
    allocate(heHeating(1:nVoxels))
    allocate(distanceGrid(1:nVoxels))
    allocate(photoIonCoeff(1:nVoxels, 1:grid%nIon))

    nIndex = 0
    call packValues(grid%octreeRoot,nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid)


    allocate(tempDoubleArray(nVoxels))
    allocate(tempRealArray(nVoxels))

    do i = 1, grid%nIon
      tempDoubleArray = 0.d0
      call MPI_ALLREDUCE(photoIonCoeff(1:nVoxels,i),tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
          MPI_SUM,MPI_COMM_WORLD,ierr)
       photoIonCoeff(1:nVoxels, i) = tempDoubleArray 
    enddo

    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(hHeating,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    hHeating = tempDoubleArray

    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(heHeating,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    heHeating = tempDoubleArray

    tempRealArray = 0.0
    call MPI_ALLREDUCE(nCrossings,tempRealArray,nVoxels,MPI_REAL,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    nCrossings = tempRealArray 

    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(distanceGrid,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    distanceGrid = tempDoubleArray 
    
    deallocate(tempRealArray, tempDoubleArray)
     
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
    
    nIndex = 0
    call unpackValues(grid%octreeRoot, nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid)

    deallocate(nCrossings, photoIonCoeff, hHeating, heHeating, distanceGrid)

  end subroutine updateGridMPIphoto
#endif

  recursive subroutine packvalues(thisOctal,nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real :: nCrossings(:)
  real(double) :: photoIonCoeff(:,:)
  real(double) :: hHeating(:)
  real(double) :: heHeating(:)
  real(double) :: distanceGrid(:)

  integer :: nIndex
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call packvalues(child, nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid)
                exit
             end if
          end do
       else
          nIndex = nIndex + 1
          nCrossings(nIndex) = real(thisOctal%nCrossings(subcell))
          photoIonCoeff(nIndex, :) = thisOctal%photoIonCoeff(subcell, :)
          hHeating(nIndex) = thisOctal%hHeating(subcell)
          heHeating(nIndex) = thisOctal%heHeating(subcell)
          distanceGrid(nIndex) = thisOctal%distanceGrid(subcell)
       endif
    enddo
  end subroutine packvalues

  recursive subroutine unpackvalues(thisOctal,nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real :: nCrossings(:)
  real(double) :: photoIonCoeff(:,:)
  real(double) :: hHeating(:)
  real(double) :: heHeating(:)
  real(double) :: distanceGrid(:)

  integer :: nIndex
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unpackvalues(child, nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid)
                exit
             end if
          end do
       else
          nIndex = nIndex + 1
          thisOctal%nCrossings(subcell) = int(nCrossings(nIndex))
          thisOctal%photoIonCoeff(subcell, :) = photoIonCoeff(nIndex, :)
          thisOctal%hHeating(subcell) = hHeating(nIndex) 
          thisOctal%heHeating(subcell) = heHeating(nIndex) 
          thisOctal%distanceGrid(subcell) = distanceGrid(nIndex) 
       endif
    enddo
  end subroutine unpackvalues

  recursive subroutine  identifyUndersampled(thisOctal)
    use input_variables, only : minCrossings
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call identifyUndersampled(child)
                exit
             end if
          end do
       else
          if (thisOctal%nCrossings(subcell) < minCrossings) then
             thisOctal%undersampled(subcell) = .true.
          else
             thisOctal%undersampled(subcell) = .false.
          endif
       endif
    enddo
  end subroutine identifyUndersampled

  recursive subroutine  calcContinuumEmissivity(grid, thisOctal, nlambda, lamArray)
    type(GRIDTYPE) :: grid
    integer :: nFreq
    real(double), allocatable :: freq(:), spectrum(:), dfreq(:)
    integer :: nLambda
    real :: lamArray(:)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i


    nFreq = nlambda

    allocate(freq(1:nFreq), spectrum(1:nFreq), dfreq(1:nFreq))

    do i = 1, nFreq
       freq(i) = cSpeed/(lamArray(nFreq-i+1)*1.e-8)
    enddo
    do i = 2, nFreq-1
       dfreq(i) = (freq(i+1)-freq(i-1))/2.d0
    enddo
    dfreq(1) = (freq(2)-freq(1))
    dfreq(nfreq) = (freq(nfreq)-freq(nfreq-1))



  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calcContinuumEmissivity(grid, child, nlambda, lamArray)
                exit
             end if
          end do
       else
          spectrum = 1.d-30
          thisOctal%etaCont(subcell) = 1.d-40
          if (thisOctal%temperature(subcell) > 1.5d0) then
             call addLymanContinua(nFreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
             call addHigherContinua(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, GammaTableArray)
!             call addHydrogenRecombinationLines(nfreq,  freq, dfreq, spectrum, thisOctal, subcell, grid)
             !         call addHeRecombinationLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
             call addForbiddenLines(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
             call addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, nlambda, lamArray)
             
             do i = 1, nFreq
                thisOctal%etaCont(subcell) = thisOctal%etaCont(subcell) + spectrum(i)
             enddo
          endif

       endif
    enddo
  end subroutine calcContinuumEmissivity


  function getRandomWavelengthPhotoion(grid, thisOctal, subcell, lamArray, nLambda) result(thisLambda)

    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real :: thisLambda
    integer :: nFreq
    real(double), allocatable :: freq(:), spectrum(:), tspec(:),lamspec(:), dfreq(:)
    real(double) :: nuStart, nuEnd, r, fac
    integer :: nLambda
    real :: lamArray(:)
    integer :: i

    nFreq = 1000

    allocate(freq(1:nFreq), spectrum(1:nFreq), lamSpec(1:nFreq), dFreq(1:nFreq))
    nuStart = cSpeed / (1000.d4 * 1.d-8)
    nuEnd =  2.d0*maxval(grid%ion(1:grid%nIon)%nuThresh)

    do i = 1, nFreq
       freq(i) = log10(nuStart) + dble(i-1)/dble(nFreq-1) * (log10(nuEnd)-log10(nuStart))
       freq(i) = 10.d0**freq(i)
    enddo
    do i = 2, nFreq-1
       dfreq(i) = (freq(i+1)-freq(i-1))/2.d0
    enddo
    dfreq(1) = (freq(2)-freq(1))
    dfreq(nfreq) = (freq(nfreq)-freq(nfreq-1))


    spectrum = 1.d-30

    call addLymanContinua(nFreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
    call addHigherContinua(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, GammaTableArray)
    call addHydrogenRecombinationLines(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
    !                        call addHeRecombinationLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
    call addForbiddenLines(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
    call addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, nlambda, lamArray)
    

    do i = 1, nFreq
       lamSpec(i) = cspeed/freq(i)
       spectrum(i) = spectrum(i) * cspeed/(lamspec(i)**2)
    enddo
    lamSpec = lamSpec * 1.d8

    allocate(tSpec(1:nFreq))
  
    tSpec(1:nFreq) = spectrum(1:nFreq)

    do i = 2, nFreq
       tSpec(i) = tSpec(i) + tSpec(i-1)
    enddo
    tSpec(1:nFreq) = tSpec(1:nFreq) - tSpec(1)
    if (tSpec(nFreq) > 0.d0) then
       tSpec(1:nFreq) = tSpec(1:nFreq) / tSpec(nFreq)
       call random_number(r)
       call locate(tSpec, nFreq, r, i)
       fac = (r - tSpec(i)) / (tSpec(i+1)-tSpec(i))
       thisLambda = lamspec(i) + fac * (lamspec(i+1)-lamspec(i))
    else
       thisLambda = 1000.e4
  endif
    
  deallocate(freq, spectrum, lamSpec)

  end function getRandomWavelengthPhotoion

  subroutine getWavelengthBiasPhotoion(grid, thisOctal, subcell, lamArray, dlam, nLambda, ilambda, bias, useBias)

    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer :: nFreq
    real(double), allocatable :: freq(:), dfreq(:), spectrum(:), tspec(:),lamspec(:)
    real :: dlam(:), dlam2
    real(double), allocatable :: prob(:)
    real(double) :: t
    real :: thisLam
    real(double) :: thisFreq, r, fac
    integer :: nLambda
    real :: lamArray(:)
    real :: bias
    integer :: ilambda
    integer :: i
    logical :: useBias

    nFreq = nlambda

    allocate(freq(1:nFreq), spectrum(1:nFreq), lamSpec(1:nFreq), tSpec(1:nFreq), dfreq(1:nFreq))

    do i = 1, nFreq
       freq(i) = cSpeed/(lamArray(nFreq-i+1)*1.e-8)
    enddo
    do i = 2, nFreq-1
       dfreq(i) = (freq(i+1)-freq(i-1))/2.d0
    enddo
    dfreq(1) = 2.d0*(freq(2)-freq(1))
    dfreq(nfreq) = 2.d0*(freq(nfreq)-freq(nfreq-1))


    spectrum = 1.d-50

    call addLymanContinua(nFreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
!       if (iLambda == 1) then
!          open(31,file="crap.dat",status="unknown",form="formatted")
!          do i = 1, nLambda
!             write(31,*) lamArray(i), spectrum(i),dFreq(nfreq-i+1)
!          enddo
!          close(31)
!          stop
!       endif
    call addHigherContinua(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, GammaTableArray)
    call addHydrogenRecombinationLines(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
!    !                        call addHeRecombinationLines(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
    call addForbiddenLines(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)

    call addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, nlambda, lamArray)
    


    if (useBias) then
       fac = 0.d0
       do i = 1, nFreq
          fac = fac + spectrum(i)
       enddo


       do i = 1, nFreq
          lamSpec(i) = cspeed/freq(i)
       enddo
       lamSpec = lamSpec * 1.d8
       
       do i = 1, nFreq
          tSpec(i) = spectrum(nFreq-i+1)
       enddo
       spectrum(1:nFreq) = tspec
       
       do i = 1, nFreq
          tSpec(i) = lamSpec(nFreq-i+1)
       enddo
       lamSpec(1:nFreq)= tspec
       do i = 1, nFreq
          dlam2 = 1.d8 * (cSpeed / freq(nFreq-i+1)**2) * dfreq(nFreq-i+1)
          spectrum(i) = spectrum(i) / dlam2
       enddo
       bias = spectrum(iLambda)/fac
    else
       allocate(prob(1:nFreq))
       prob(1:nFreq) = spectrum(1:nFreq)
       do i = 2, nFreq
          prob(i) = prob(i-1) + prob(i)
       enddo
       prob(1:nFreq) = prob(1:nFreq)-prob(1)
       prob(1:nFreq) = prob(1:nFreq)/prob(nFreq)
       call random_number(r)
       call locate(prob, nFreq, r, i)
       t = (r - prob(i))/(prob(i+1)-prob(i))
       thisFreq = freq(i) + t * (freq(i+1)-freq(i))
       thisLam = 1.d8 * cspeed/thisFreq
       call locate(lamArray, nLambda, thisLam, iLambda)
       bias = 1.d0
    endif


    deallocate(freq, spectrum, lamSpec, tSpec)

  end subroutine getWavelengthBiasPhotoion

  subroutine refineLambdaArray(lamArray, nLambda, grid)
    type(GRIDTYPE) :: grid
    real :: lamArray(:)
    integer :: nLambda
    integer :: i, j
    integer :: iup, ilow

  real(double) :: lambdaTrans(20,20) = reshape( source=&
    (/ 000000.d-8,1215.67d-8,1025.72d-8,992.537d-8,949.743d-8,937.803d-8,930.748d-8,926.226d-8,923.150d-8,920.963d00,&
       919.352d-8,918.129d-8,917.181d-8,916.429d-8,915.824d-8,915.329d-8,914.919d-8,914.576d-8,914.286d-8,914.039d-8,&  
       0.00000d-8,0000000d-8,6562.80d-8,4861.32d-8,4340.36d-8,4101.73d-8,3970.07d-8,3889.05d-8,3835.38d-8,3797.90d-8,&
       3770.63d-8,3750.15d-8,3734.37d-8,3721.94d-8,3711.97d-8,3703.85d-8,3697.15d-8,3691.55d-8,3686.83d-8,3682.81d-8,&  
       0.00000d-8,0.00000d-8,0000000d-8,18751.0d-8,12818.1d-8,10938.1d-8,10049.4d-8,9545.98d-8,9229.02d-8,9014.91d-8,&
       8862.79d-8,8750.47d-8,8665.02d-8,8598.39d-8,8545.39d-8,8502.49d-8,8467.26d-8,8437.96d-8,8413.32d-8,8392.40d-8,&  
       0.00000d-8,0.00000d-8,0.00000d-8,0000000d-8,40512.0d-8,26252.0d-8,21655.0d-8,19445.6d-8,18174.1d-8,17362.1d-8,&
       16806.5d-8,16407.2d-8,16109.3d-8,15880.5d-8,15700.7d-8,15556.5d-8,15438.9d-8,15341.8d-8,15260.6d-8,15191.8d-8,&  
       0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0000000d-8,74578.0d-8,46525.0d-8,37395.0d-8,32961.0d-8,30384.0d-8,&
       28722.0d-8,27575.0d-8,26744.0d-8,26119.0d-8,25636.0d-8,25254.0d-8,24946.0d-8,24693.0d-8,24483.0d-8,24307.0d-8,&  
       0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,123680.d-8,75005.0d-8,59066.0d-8,51273.0d-8,&
       46712.0d-8,43753.0d-8,41697.0d-8,40198.0d-8,39065.0d-8,38184.0d-8,37484.0d-8,36916.0d-8,36449.0d-8,36060.0d-8,&  
       0.00000d-8,0000000d-8,0000000d-8,0000000d-8,0000000d-8,0.00000d-8,0000000d-8,190570.d-8,113060.d-8,87577.0d-8,&
       75061.0d-8,67701.0d-8,62902.0d-8,59552.0d-8,57099.0d-8,55237.0d-8,53783.0d-8,52622.5d-8,51679.0d-8,50899.0d-8,&  
       0.00000d-8,0000000d-8,0000000d-8,0000000d-8,0000000d-8,0.00000d-8,0000000d-8,0000000d-8,277960.d-8,162050.d-8,&
       123840.d-8,105010.d-8,93894.0d-8,86621.0d-8,81527.0d-8,77782.0d-8,74930.0d-8,72696.0d-8,70908.0d-8,69448.0d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,388590.d-8,&  
       223340.d-8,168760.d-8,141790.d-8,125840.d-8,115360.d-8,108010.d-8,102580.d-8,98443.0d-8,95191.0d-8,92579.0d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       525200.d-8,298310.d-8,223250.d-8,186100.d-8,164070.d-8,149580.d-8,139380.d-8,131840.d-8,126080.d-8,121530.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,690500.d-8,388320.d-8,288230.d-8,238620.d-8,209150.d-8,189730.d-8,176030.d-8,165900.d-8,158120.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,887300.d-8,494740.d-8,364610.d-8,300020.d-8,261610.d-8,236260.d-8,218360.d-8,205090.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,111800.d-7,619000.d-8,453290.d-8,371000.d-8,322000.d-8,289640.d-8,266740.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,138600.d-7,762300.d-8,555200.d-8,452220.d-8,390880.d-8,350300.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,169400.d-7,926100.d-8,671200.d-8,544400.d-8,468760.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,204400.d-7,111200.d-7,802300.d-8,648200.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,243800.d-7,132100.d-7,949200.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,288200.d-7,155400.d-7,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,337400.d-7,&  
       0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /), shape=(/20,20/))


  do iup = 15, 3, -1
     do ilow = 2, min0(8, iup-1)
        call insertBin(lamArray, nLambda, real(lambdaTrans(iup, ilow)*1.e8), 1.)
    enddo
 enddo
 call insertBin(lamArray, nLambda, 1215.67, 1.)


! refine for forbidden line transitions

    do i = 1 , grid%nIon
       do j = 1, grid%ion(i)%nTransitions
          call insertBin(lamArray, nLambda, &
               real(grid%ion(i)%transition(j)%lambda), 1.)
       enddo
    enddo

     
  end subroutine refineLambdaArray

#ifdef MPI

      subroutine packTemperatures(octalArray, nTemps, tArray, octalsBelongRank)
    include 'mpif.h'
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: octalsBelongRank(:)
        integer :: nTemps
        real :: tArray(:)
        integer :: iOctal, iSubcell, my_rank, ierr
        type(OCTAL), pointer :: thisOctal

       call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
       !
       ! Update the edens values of grid computed by all processors.
       !
       nTemps = 0
       do iOctal = 1, SIZE(octalArray)

          thisOctal => octalArray(iOctal)%content
          
          do iSubcell = 1, thisOctal%maxChildren
              if (.not.thisOctal%hasChild(iSubcell)) then
                 nTemps = nTemps + 1
                 if (octalsBelongRank(iOctal) == my_rank) then
                   tArray(nTemps) = thisOctal%temperature(isubcell)
                 else 
                   tArray(nTemps) = 0.d0
                 endif
              endif
          end do
       end do
     end subroutine packTemperatures

      subroutine unpackTemperatures(octalArray, nTemps, tArray)
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: nTemps
        real :: tArray(:)
        integer :: iOctal, iSubcell
        type(OCTAL), pointer :: thisOctal

       !
       ! Update the edens values of grid computed by all processors.
       !
       nTemps = 0
       do iOctal = 1, SIZE(octalArray)

          thisOctal => octalArray(iOctal)%content
          
          do iSubcell = 1, thisOctal%maxChildren
                
              if (.not.thisOctal%hasChild(iSubcell)) then
                 nTemps = nTemps + 1
                 thisOctal%temperature(iSubcell) = tArray(nTemps)
              endif
          
          end do
       end do
     end subroutine unpackTemperatures

      subroutine packNe(octalArray, nTemps, tArray, octalsBelongRank)
    include 'mpif.h'
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: octalsBelongRank(:)
        integer :: nTemps
        real(double) :: tArray(:)
        integer :: iOctal, iSubcell, my_rank, ierr
        type(OCTAL), pointer :: thisOctal

       call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
       !
       ! Update the edens values of grid computed by all processors.
       !
       nTemps = 0
       do iOctal = 1, SIZE(octalArray)

          thisOctal => octalArray(iOctal)%content
          
          do iSubcell = 1, thisOctal%maxChildren
              if (.not.thisOctal%hasChild(iSubcell)) then
                 nTemps = nTemps + 1
                 if (octalsBelongRank(iOctal) == my_rank) then
                   tArray(nTemps) = thisOctal%ne(isubcell)
                 else 
                   tArray(nTemps) = 0.d0
                 endif
              endif
          end do
       end do
     end subroutine packNe

      subroutine unpackNe(octalArray, nTemps, tArray)
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: nTemps
        real(double) :: tArray(:)
        integer :: iOctal, iSubcell
        type(OCTAL), pointer :: thisOctal

       !
       ! Update the edens values of grid computed by all processors.
       !
       nTemps = 0
       do iOctal = 1, SIZE(octalArray)

          thisOctal => octalArray(iOctal)%content
          
          do iSubcell = 1, thisOctal%maxChildren
                
              if (.not.thisOctal%hasChild(iSubcell)) then
                 nTemps = nTemps + 1
                 thisOctal%ne(iSubcell) = tArray(nTemps)
              endif
          
          end do
       end do
     end subroutine unpackne

      subroutine packIonFrac(octalArray, nTemps, tArray, octalsBelongRank, iIon)
    include 'mpif.h'
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: octalsBelongRank(:)
        integer :: nTemps
        real(double) :: tArray(:)
        integer :: iOctal, iSubcell, my_rank, ierr
        integer :: iIon
        type(OCTAL), pointer :: thisOctal

       call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
       !
       ! Update the edens values of grid computed by all processors.
       !
       nTemps = 0
       do iOctal = 1, SIZE(octalArray)

          thisOctal => octalArray(iOctal)%content
          
          do iSubcell = 1, thisOctal%maxChildren
              if (.not.thisOctal%hasChild(iSubcell)) then
                 nTemps = nTemps + 1
                 if (octalsBelongRank(iOctal) == my_rank) then
                   tArray(nTemps) = thisOctal%ionFrac(isubcell, iIon)
                 else 
                   tArray(nTemps) = 0.d0
                 endif
              endif
          end do
       end do
     end subroutine packIonFrac

      subroutine unpackIonFrac(octalArray, nTemps, tArray, octalsBelongRank, iIon)
    include 'mpif.h'
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: octalsBelongRank(:)
        integer :: nTemps
        real(double) :: tArray(:)
        integer :: iOctal, iSubcell, my_rank, ierr
        integer :: iIon
        type(OCTAL), pointer :: thisOctal

       call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
       !
       ! Update the edens values of grid computed by all processors.
       !
       nTemps = 0
       do iOctal = 1, SIZE(octalArray)

          thisOctal => octalArray(iOctal)%content
          
          do iSubcell = 1, thisOctal%maxChildren
              if (.not.thisOctal%hasChild(iSubcell)) then
                 nTemps = nTemps + 1
                 thisOctal%ionFrac(isubcell, iIon) = tArray(nTemps) 
              endif
          end do
       end do
     end subroutine unpackIonFrac
#endif

  subroutine writeMultiImages(grid, nSource, source, observerDirection)
    type(GRIDTYPE) :: grid
    integer :: nSource
    type(SOURCETYPE) :: source(:)
    type(VECTOR) :: observerDirection
    real(double) :: lambdaLine(6)
    integer :: nLine , i
    real(double) :: totalFlux

    nLine = 6
    lambdaLine(1) = 6300.d0
    lambdaLine(2) = 6363.d0
    lambdaLine(3) = 3726.d0
    lambdaLine(4) = 3729.d0
    lambdaLine(5) = 4959.d0
    lambdaLine(6) = 5007.d0
    do i = 1, nLine
       call createImage(grid, nSource, source, observerDirection, totalflux, lambdaLine(i))
    enddo
  end subroutine writeMultiImages

  subroutine createImage(grid, nSource, source, observerDirection, totalflux, lambdaLine)
    use input_variables, only : readlucy, nlambda
    type(GRIDTYPE) :: grid
    character(len=80) :: imageFilename
    integer :: nSource
    real(double) :: lambdaLine
    type(SOURCETYPE) :: source(:), thisSource
    type(PHOTON) :: thisPhoton, observerPhoton
    type(OCTAL), pointer :: thisOctal
    real(double) :: totalFlux, tempTotalFlux
    integer :: subcell
    integer :: iPhoton, nPhotons
    integer :: iLam
    integer :: iSource
    integer :: iThread
    type(VECTOR) :: rVec, uHat, rHat, observerDirection
    real(double) :: wavelength, thisFreq
    logical :: endLoop, addToiMage
    integer :: newthread
    type(IMAGETYPE) :: thisimage
    logical :: escaped, absorbed, crossedBoundary, photonsStillProcessing, stillSCattering
    real(double) :: totalEmission
    integer :: iLambdaPhoton, nInf, i
    real(double) :: lCore, probsource, r
    real(double), allocatable :: threadProbArray(:)
    integer :: np(10)
    integer :: nDone
    real(double) :: powerPerPhoton
    real(double) :: totalLineEmission, totalContEmission

    call init_random_seed()




    thisImage = initImage(50, 50, real(4.*grid%octreeRoot%subcellSize), &
         real(4.*grid%octreeRoot%subcellSize), 0., 0.)

    
    call addEmissionLine(grid, 1.d0, lambdaLine)
    call locate(grid%lamArray, grid%nlambda, real(lambdaLine), ilambdaPhoton)
    call computeProbDist(grid, totalLineEmission, totalEmission, 1.0, .false.)

    totalEmission = totalEmission * 1.d30

    if (nSource > 0) then              
       lCore = sumSourceLuminosityMonochromatic(source, nsource, dble(grid%lamArray(ilambdaPhoton)))
    else
       lcore = tiny(lcore)
    endif

    write(*,*) "Total line emission ",totalEmission
    write(*,*) "Total source emission ",lCore

    Probsource = lCore / (lCore + totalEmission)

    if (myRankGlobal == 0) then
       write(*,*) "Probability of photon from sources: ", probSource
    endif

    nPhotons = 1000000
    nInf = 0

    powerPerPhoton = (lCore + totalEmission) / dble(nPhotons)
    

    mainloop: do iPhoton = 1, nPhotons

       thisPhoton%stokes = STOKESVECTOR(1.d0, 0.d0, 0.d0, 0.d0)
       thisPhoton%iLam = iLambdaPhoton
       thisPhoton%lambda = grid%lamArray(iLambdaPhoton)
       thisPhoton%observerPhoton = .false.
       call random_number(r)


       if (r < probSource) then
          call randomSource(source, nSource, iSource)
          thisSource = source(iSource)
          call getPhotonPositionDirection(thisSource, thisPhoton%position, thisPhoton%direction, rHat,grid)         
       else
          call random_number(r)
          thisPhoton%lambda = grid%lamArray(iLambdaPhoton)
          thisPhoton%direction = randomUnitVector()
	  thisOctal => grid%octreeRoot
          call locateContProbAMR(r,thisOctal,subcell)
          thisPhoton%position = randomPositionInCell(thisOctal, subcell)
       endif

       observerPhoton = thisPhoton
       observerPhoton%observerPhoton = .true.
       observerPhoton%tau = 0.d0
       observerPhoton%direction = observerDirection
       
       call propagateObserverPhoton(grid, observerPhoton)
       call addPhotonToImageLocal(observerDirection, thisImage, observerPhoton, totalFlux)

       stillScattering = .true.
       ninf = ninf + 1
       endloop = .false.

       do while ((.not.endLoop).and.stillScattering)
          call moveToNextScattering(grid, thisPhoton, escaped, absorbed)
                
          if (escaped.or.absorbed) then
             stillScattering = .false.
             exit
          endif
          
          call scatterPhotonLocal(grid, thisPhoton)
          observerPhoton = thisPhoton
          observerPhoton%observerPhoton = .true.
          observerPhoton%tau = 0.d0
          observerPhoton%direction = observerDirection
          call propagateObserverPhoton(grid, observerPhoton)          
          call addPhotonToImageLocal(observerDirection, thisImage, observerPhoton, totalFlux)
       end do
    end do mainloop

#ifdef MPI
    call collateImages(thisImage)
     call MPI_ALLREDUCE(totalFlux, tempTotalFlux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     totalFlux = tempTotalFlux
#endif

    if (writeoutput) then
       write(imageFilename, '(a,i4.4,a)') "test_",nint(lambdaLine),".fits"
       call writeFitsImage(thisimage, imageFilename, 1.d0, "intensity")
    endif
    call freeImage(thisImage)
  end subroutine createImage


  subroutine scatterPhotonLocal(grid, thisPhoton)
    type(GRIDTYPE) :: grid
    type(PHOTON) :: thisPhoton

    thisPhoton%direction = randomUnitVector() !isotropic scattering

  end subroutine scatterPhotonLocal

  subroutine propagateObserverPhoton(grid, thisPhoton)
    type(GRIDTYPE) :: grid
    type(PHOTON) :: thisPhoton
    logical :: addToImage
    logical :: endLoop
    integer :: newThread
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: tVal
    real(double) :: kappaAbsGas, kappaScaGas, kappaExt
    integer :: ilam
    thisOctal => grid%octreeRoot
    call findSubcellLocal(thisPhoton%position, thisOctal, subcell)

    endLoop = .false.
    do while (.not.endLoop)
       call distanceToCellBoundary(grid, thisPhoton%position, thisPhoton%direction, tval, thisOctal, subcell)
       call returnKappa(grid, thisOctal, subcell, ilambda=thisPhoton%ilam, &
            kappaAbs=kappaAbsGas, kappaSca=kappaScaGas)
       kappaExt = kappaAbsGas + kappaScaGas
       thisPhoton%tau = thisPhoton%tau + tval * kappaExt
       thisPhoton%position = thisPhoton%position + (tVal + 1.d-3*grid%halfSmallestSubcell) * thisPhoton%direction
       if (.not.inOctal(grid%octreeRoot, thisPhoton%position)) then
          endLoop = .true.
       else
          call findSubcellLocal(thisPhoton%position, thisOctal, subcell)
       endif
    enddo
  end subroutine propagateObserverPhoton

  subroutine moveToNextScattering(grid, thisPhoton, escaped, absorbed)
    type(GRIDTYPE) :: grid
    type(PHOTON) :: thisphoton
    logical :: escaped, absorbed,  scattered
    integer :: newThread
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: tau, thisTau, tVal, r, albedo
    real(double) ::  kappaAbsGas, kappaScaGas,  kappaExt
    logical :: endLoop
    
    thisOctal => grid%octreeRoot
    call findSubcellLocal(thisPhoton%position, thisOctal, subcell)
    
    absorbed = .false.
    escaped = .false.
    endLoop = .false.
    scattered = .false.
    do while (.not.endLoop)
       call distanceToCellBoundary(grid, thisPhoton%position, thisPhoton%direction, tval, thisOctal, subcell)
       call returnKappa(grid, thisOctal, subcell, ilambda=thisPhoton%ilam, &
            kappaAbs=kappaAbsGas, kappaSca=kappaScaGas)
       kappaExt = kappaAbsGas + kappaScaGas
       tau = kappaExt * tVal

       call random_number(r)
       thisTau = -log(1.d0-r)

       if (thisTau > tau) then ! photon crosses to boundary
          thisPhoton%position = thisPhoton%position + (tVal + 1.d-2*grid%halfSmallestSubcell) * thisPhoton%direction
          if (.not.inOctal(grid%octreeRoot, thisPhoton%position)) then
             escaped = .true.
             endLoop = .true.
          else
             call findSubcellLocal(thisPhoton%position, thisOctal, subcell)
          endif
       else

          thisPhoton%position = thisPhoton%position + ((thisTau/tau)*tVal) * thisPhoton%direction

          endLoop = .true.
          albedo = kappaScaGas/kappaExt
          call random_number(r)
          if (r < albedo) then
             scattered = .true.
          else
             absorbed = .true.
          endif

       endif
    enddo
  end subroutine moveToNextScattering


#ifdef MPI

  subroutine collateImages(thisImage)
    include 'mpif.h'
    type(IMAGETYPE) :: thisImage
    real, allocatable :: tempRealArray(:), tempRealArray2(:)
    real(double), allocatable :: tempDoubleArray(:), tempDoubleArray2(:)
    integer :: ierr

     allocate(tempRealArray(SIZE(thisImage%pixel)))
     allocate(tempRealArray2(SIZE(thisImage%pixel)))
     allocate(tempDoubleArray(SIZE(thisImage%pixel)))
     allocate(tempDoubleArray2(SIZE(thisImage%pixel)))
     tempRealArray = 0.0
     tempRealArray2 = 0.0
     tempDoubleArray = 0.0_db
     tempDoubleArray2 = 0.0_db

     if (myrankGlobal == 1) write(*,*) "Collating images..."
     tempDoubleArray = reshape(thisImage%pixel%i,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%pixel%i = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%i))

     tempDoubleArray = reshape(thisImage%pixel%q,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%pixel%q = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%q))

     tempDoubleArray = reshape(thisImage%pixel%u,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%pixel%u = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%u))

     tempDoubleArray = reshape(thisImage%pixel%v,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%pixel%v = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%v))


     tempRealArray = reshape(thisImage%vel,(/SIZE(tempRealArray)/))
     call MPI_REDUCE(tempRealArray,tempRealArray2,SIZE(tempRealArray),MPI_REAL,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%vel = reshape(tempRealArray2,SHAPE(thisImage%vel))

     tempRealArray = reshape(thisImage%totWeight,(/SIZE(tempRealArray)/))
     call MPI_REDUCE(tempRealArray,tempRealArray2,SIZE(tempRealArray),MPI_REAL,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%totWeight = reshape(tempRealArray2,SHAPE(thisImage%totWeight))

     if (myrankGlobal == 1) write(*,*) "Done."
     deallocate(tempRealArray)
     deallocate(tempRealArray2)
     deallocate(tempDoubleArray)
     deallocate(tempDoubleArray2)

  end subroutine collateImages

#endif

   subroutine addPhotonToImageLocal(observerDirection, thisImage, thisPhoton, totalFlux)
     
     type(IMAGETYPE), intent(inout) :: thisImage
     type(PHOTON) :: thisPhoton
     type(VECTOR) :: observerDirection,  xProj, yProj, rotationAxis
     real :: xDist, yDist
     integer :: xPix, yPix
     integer :: i
     real(double) :: totalFlux

     type(VECTOR), parameter :: zAxis = VECTOR(0.d0, 0.d0, 1.d0)

     xPix = 0; yPix = 0


     xProj =  zAxis .cross. observerDirection
     call normalize(xProj)
     yProj = observerDirection .cross. xProj
     call normalize(yProj)
     xDist = (thisPhoton%position) .dot. xProj
     yDist = (thisPhoton%position) .dot. yProj
           

     call pixelLocate(thisImage, xDist, yDist, xPix, yPix)

     if ((xPix >= 1) .and. &
          (yPix >= 1) .and. &
          (xPix <= thisImage%nx) .and. &
          (yPix <= thisImage%ny)) then

              
        thisImage%pixel(xPix, yPix) = thisImage%pixel(xPix, yPix)  &
             + thisPhoton%stokes * oneOnFourPi * exp(-thisPhoton%tau)
     endif
     totalFlux = totalFlux + thisPhoton%stokes%i * oneOnFourPi * exp(-thisPhoton%tau)
           
   end subroutine addPhotonToImageLocal


   subroutine addEmissionLine(grid, dLambda, lineWavelength)
     type(GRIDTYPE) :: grid
     real(double) :: lineWavelength, dLambda
     integer :: iIon, iTransition
     logical :: ok

     call identifyTransition(grid, lineWavelength, iIon, iTransition, ok)

     if (ok) then
        call addToEmission(grid, grid%octreeRoot, dLambda, iIon, iTransition)
     endif

     
   end subroutine addEmissionLine


   recursive subroutine addToEmission(grid, thisOctal, dLambda, iIon, iTransition)
     type(GRIDTYPE) :: grid
     real(double) :: dLambda, rate
     real :: pops(10)
     integer :: iIon, iTransition
     type(octal), pointer   :: thisOctal
     type(octal), pointer  :: child 
     integer :: subcell, i


     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do i = 1, thisOctal%nChildren, 1
              if (thisOctal%indexChild(i) == subcell) then
                 child => thisOctal%child(i)
                 call addToEmission(grid, child, dLambda, iIon, iTransition)
                 exit
              end if
           end do
        else

           call solvePops(grid%ion(iIon), pops, thisOctal%ne(subcell), thisOctal%temperature(subcell), &
                thisOctal%ionFrac(subcell,iion),thisOctal%nh(subcell))
           rate =  pops(grid%ion(iion)%transition(iTransition)%j) * grid%ion(iion)%transition(itransition)%energy * &
                grid%ion(iion)%transition(itransition)%a/ergtoev
           rate = rate * grid%ion(iion)%abundance * thisOctal%nh(subcell) * thisOctal%ionFrac(subcell, iion)

           thisOctal%etaCont(subcell) =  rate/dLambda


        endif
     enddo
   end subroutine addToEmission

  subroutine identifyTransition(grid, lambda, iIon, iTransition, ok)
    type(GRIDTYPE) :: grid
    real(double) :: lambda
    integer :: iIon, iTransition, i, j
    character(len=80) :: message
    logical :: ok

    ok = .false.

    do i = 1, grid%nIon
       do j = 1, grid%ion(i)%nTransitions
          if (abs(lambda-grid%ion(i)%transition(j)%lambda) < 1.d0) then
             write(message, '(a, a, a)') "Transition of ", trim(grid%ion(i)%species), " identified."
             call writeInfo(message)
             iIon = i
             iTransition = j
             ok = .true.
          end if
       end do
    enddo
  end subroutine identifyTransition

     

end module

