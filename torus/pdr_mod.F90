module pdr_mod
#ifdef PDR
!Started by THaw and Bisbas 13/08/2013
 
use constants_mod
use messages_mod
use parallel_mod
use gridio_mod
use source_mod
use timing
use grid_mod
use vtk_mod
use amr_mod
use mpi_amr_mod
use mpi_global_mod
use unix_mod, only: unixGetenv
!use definitions
use healpix_guts
use pdr_utils_mod
!use healpix_module
use vector_mod
use setuppdr_mod
#ifdef MPI
use mpi
#endif
implicit none


!
!relative_abundance_tolerance = 1.d-8
!absolute_abundance_tolerance = 1.d-30
!nelect
!
!

! REAL(double) :: OMEGA=0.42D0,GRAIN_RADIUS=1.0D-7,METALLICITY=1.0D0


! real(kind=dp),allocatable :: all_heating(:,:)
! real(kind=dp),allocatable :: allheating(:) 

contains

subroutine PDR_MAIN(grid)
  use unix_mod, only: unixGetenv

  use nrayshealpix
  implicit none
  type(gridtype) :: grid
  real(double), allocatable :: rate(:), alpha(:), beta(:)
  real(double), allocatable :: rtmin(:), rtmax(:), gamma(:)!, duplicate(:)
  integer, allocatable :: duplicate(:)
  character(len=10), allocatable :: product(:,:), reactant(:,:)
  character(len=80) :: datfilename
  integer :: n12co, nci, ncii, noi, nelect
#ifdef MPI
  integer :: ier
#endif

  call donrayshealpix()

!set up the data/allocate what needs allocated

#ifdef MPI
  call writeInfo("Setting up for PDR calculation.", TRIVIAL)
  call setupPDR(grid, reactant, product, alpha, beta, gamma, &
       rate, duplicate, rtmin, rtmax, n12co, nci, ncii, noi, nelect)


  if(grid%octreeroot%oned) then
     write(datFilename, '(a, i4.4, a)') "setup.dat"
     call dumpValuesAlongLinePDR(grid, datFileName, VECTOR(0.d0,  0.d0, 0.d0), &
          VECTOR(5.16d0*pctocm/1.d10, 0.d0, 0.d0), 1000)
  endif

  call writeInfo("Casting rays over grid MPI.", TRIVIAL)
  call rayTraceMPI(grid, colrhoonly=.false.)

  call MPI_BARRIER(MPI_COMM_WORLD, ier)
  if(grid%octreeroot%oned) then
     write(datFilename, '(a, i4.4, a)') "postraytrace.dat"
     call dumpValuesAlongLinePDR(grid, datFileName, VECTOR(0.d0,  0.d0, 0.d0), &
          VECTOR(5.16d0*pctocm/1.d10, 0.d0, 0.d0), 1000)
  endif

  call writeVTKfile(grid, "postRayTrace.vtk", valueTypeString=(/"rho       ",&
       "columnRho " /))
!  stop
  call writeInfo("Done.", TRIVIAL)
  call writeInfo("Calculating dust temperature MPI.", TRIVIAL)

  call calculate_Dust_TemperaturesMPI(grid%octreeRoot)
  call MPI_BARRIER(MPI_COMM_WORLD, ier)
  call writeInfo("Making initial gas temperature estimates.", TRIVIAL)

  call iniTempGuessMPI(grid%octreeroot)

  call writeInfo("Calculating reaction rates in each cell.", TRIVIAL)
  call MPI_BARRIER(MPI_COMM_WORLD, ier)


  if(grid%octreeroot%oned) then
     write(datFilename, '(a, i4.4, a)') "preAbund.dat"
     call dumpValuesAlongLinePDR(grid, datFileName, VECTOR(0.d0,  0.d0, 0.d0), &
          VECTOR(5.16d0*pctocm/1.d10, 0.d0, 0.d0), 1000)
  else
     call writeVTKfile(grid, "PDR_CALC.vtk", valueTypeString=(/"rho       ",&
          "columnRho ", "UV        ", "uvvec     ", "dust_T    ", "pdrtemp   "/))
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, ier)

  call writeInfo("Doing initial abundance calculation.", TRIVIAL)


 call abundanceSweepGovernor(grid, reactant, &
       product, alpha, beta, gamma, rate, duplicate, rtmin, rtmax, n12co &
       , ncii, nci, noi, nelect, nchemiter=8, partLTE=.true.) 

!  call abundanceSweepGovernor(grid, reactant, &
!       product, alpha, beta, gamma, rate, duplicate, rtmin, rtmax, n12co &
!       , ncii, nci, noi, nelect, nchemiter=100, partLTE=.true.) 
  call MPI_BARRIER(MPI_COMM_WORLD, ier)


  call writeVTKfile(grid, "PDR_ABUND.vtk", valueTypeString=(/"rho       ",&
       "columnRho ", "UV        ", "uvvec     ", "dust_T    ", "pdrtemp   ", &
       "CO_PDR    ", "C+_PDR    ", "C_PDR     ", "CII_1     ", "CII_2     ", &
       "CII_3     ", "CII_4     "/))

  write(datFilename, '(a, i4.4, a)') "start.dat"
  call dumpValuesAlongLinePDR(grid, datFileName, VECTOR(0.d0,  0.d0, 0.d0), &
       VECTOR(5.16d0*pctocm/1.d10, 0.d0, 0.d0), 1000)

  call MPI_BARRIER(MPI_COMM_WORLD, ier)
!  stop
  call writeInfo("Initial abundance calculation done.", TRIVIAL)

  print *, " "
  call writeInfo("Doing main PDR loop.", TRIVIAL)
!  print *, "CII WEIGHTS 6", CII_WEIGHTS
  call pdr_main_loop(grid, nelect, ncii, nci, noi, n12co, reactant, &
       product, alpha, beta, gamma, rate, duplicate, rtmin, rtmax)

  call writeVTKfile(grid, "PDR_ABUND.vtk", valueTypeString=(/"rho       ",&
       "columnRho ", "UV        ", "uvvec     ", "dust_T    ", "pdrtemp   ", &
       "CO_PDR    ", "C+_PDR    ", "C_PDR     ", "cii_1to0  ", "cii_line  ",& 
       "cooling   "/))

  write(datFilename, '(a, i4.4, a)') "finish.dat"
  call dumpValuesAlongLinePDR(grid, datFileName, VECTOR(0.d0,  0.d0, 0.d0), &
       VECTOR(5.16d0*pctocm/1.d10, 0.d0, 0.d0), 1000)

#else

  call writeInfo("Casting rays over grid.", TRIVIAL)
  call castAllRaysOverGrid(grid%octreeRoot, grid)
  call writeInfo("Done.", TRIVIAL)

  call writeInfo("Calculating dust temperature.", TRIVIAL)
  call calculate_Dust_Temperatures(grid%octreeRoot)
  call writeInfo("Making initial gas temperature estimates.", TRIVIAL)
  call iniTempGuess(grid%octreeroot)
#endif

end subroutine PDR_MAIN


#ifdef MPI
subroutine pdr_main_loop(grid, nelect, ncii, nci, noi, nc12o, reactant, &
       product, alpha, beta, gamma, rate, duplicate, rtmin, rtmax)
  type(gridtype) :: grid
  real(double), allocatable :: rate(:), alpha(:), beta(:), gamma(:)
  real(double), allocatable :: rtmin(:), rtmax(:)
  integer, allocatable :: duplicate(:)
  character(len=10), allocatable :: product(:,:), reactant(:,:)
  integer :: nc12o, nci, ncii, noi, nelect!, thisThread
  logical :: overallconverged
  character(len=80) :: filename
  integer :: thisIteration, maxIter, levpop_iteration
  integer :: ier, nconverged
  integer :: k, convCounter
  logical :: buffer
  !  logical, save :: firstTime=.true.
  logical :: first_time, level_conv, anyNotConverged
    character(len=80) :: datFilename
!    logical :: inEscProb = .true.

  overallconverged = .false.
  thisIteration = 1
  levpop_iteration = 0
  maxIter = 2000
!these logicals are updated in the convergence checks
  first_time = .true.
  level_conv = .false.

  call InitLogicals(grid%octreeroot)
  buffer = .true.
  do while(.not. overallconverged)
     
     if(myrankglobal == 1) then
        print *, "! - Main PDR loop, iteration ", thisIteration
     endif
     !abundance calculation
!     if(.not. inescprob) then
     if(level_conv .and. thisIteration > 1) then
        !     if(thisIteration > 1) then
        call writeInfo("Abundance sweep...", trivial)
        call abundanceSweepGovernor(grid, reactant, &
             product, alpha, beta, gamma, rate, duplicate, rtmin, rtmax, nc12o &
             , ncii, nci, noi, nelect, nchemiter=3, partlte=.true.) 
        call writeInfo("Abundance sweep done")
!        if(myrankglobal /= 0) then
!           do k = 1, nhydrothreadsglobal
!              if(myrankglobal == k) then
!                 !calculate collisional coefficients
!                 
!                 call solvePopulations(grid%octreeroot, grid, nelect, ncii, nci, noi, nc12o, level_conv, &
!                      thisIteration, coolingOnly = .true.)
!                 call shutdownserversXRAY_PDR()
!                 !
!              else
!                 !server
!                 
!                 call raytracingserverPDR_TWO(grid)
!                 
!              endif
!              call MPI_BARRIER(amrCommunicator, ier)
!           enddo
!        endif

     else
        
        call writeInfo("Solving populations", trivial)     
        if(myrankglobal /= 0) then
           do k = 1, nhydrothreadsglobal
              if(myrankglobal == k) then
                 !calculate collisional coefficients
                 
                 call solvePopulations(grid%octreeroot, grid, nelect, ncii, nci, noi, nc12o, level_conv, &
                   thisIteration, levpop_iteration, coolingOnly = .false.)
                 call shutdownserversXRAY_PDR()
                 !
              else
                 !server
                 
                 call raytracingserverPDR_TWO(grid)
                 
              endif
              call MPI_BARRIER(amrCommunicator, ier)
           enddo
        endif
     endif
     
     call MPI_BARRIER(MPI_COMM_WORLD, ier)
     !     stop
     call writeInfo("Population solver done", trivial)
     
     call MPI_BARRIER(MPI_COMM_WORLD, ier)
     level_conv = .false.
     !     if(thisIteration > 10) then
     call writeInfo("Checking convergence", trivial)
     call MPI_BARRIER(MPI_COMM_WORLD, ier)
     !convergence checks
     !     level_conv = .true.
     !        do convCounter = 1, 2
     nconverged = 0
     anynotconverged = .false.
     convcounter = 1
!     if(1==0)then
     if(myrankglobal /= 0) then

        call checkConvergence(grid%octreeRoot, anyNotconverged, nconverged, level_conv, convcounter)
        
        call send_conv_data(nconverged, anynotconverged, overallconverged, level_conv, convcounter)
     else

        call countNConverged(nconverged, overallconverged, level_conv, convcounter)
        !              if(convCounter == 2) print *, "number converged ", nconverged
 !    endif
     endif
     call MPI_BARRIER(MPI_COMM_WORLD, ier)
     

     if(level_conv) levpop_iteration = 0
     
     !thermal balance calculation
     call writeInfo("Calculating heating rates and doing thermal balance", trivial)
     !     if(level_conv 
     if(myrankglobal /= 0) then


           call calcHeatingOverGrid(grid%octreeRoot, reactant, &
                product, alpha, beta, gamma, rate, duplicate, rtmin, rtmax, nelect, &
                nc12o, ncii, nci, noi, level_conv, first_time, grid)           
           
        endif
     call MPI_BARRIER(MPI_COMM_WORLD, ier)
     if(level_conv) first_time = .false. 
     

     nconverged = 0
     anynotconverged = .false.
     convcounter = 2
     if(myrankglobal /= 0) then
        call checkConvergence(grid%octreeRoot, anyNotconverged, nconverged, level_conv, convcounter)
        
        call send_conv_data(nconverged, anynotconverged, overallconverged, level_conv, convcounter)
     else

        call countNConverged(nconverged, overallconverged, level_conv,convcounter)
        !              if(convCounter == 2) print *, "number converged ", nconverged
     endif
     !        print *, "LEVEL CONV IS ", level_conv, myrankglobal
     !     end if
     !     endif
     call MPI_BARRIER(MPI_COMM_WORLD, ier)
     !       enddo
     call MPI_BARRIER(MPI_COMM_WORLD, ier)
     
     if(level_conv .and. .not. overallconverged) then
        call zeroLevelConv(grid%octreeroot)
        !        level_conv = .false.
     endif

!     if(modulo(thisIteration, 10) == 0) then    
        if(grid%octreeroot%oned) then
           
!           write(datFilename, '(a, i4.4, a)') "pdr",thisIteration,".dat"
!           call dumpValuesAlongLinePDR(grid, datFileName, VECTOR(0.d0,  0.d0, 0.d0), &
!                VECTOR(5.16d0*pctocm/1.d10, 0.d0, 0.d0), 1000)

           write(datFilename, '(a, i4.4, a)') "heatingcooling.dat"
           call dumpHeatingCooling(grid, datFileName, VECTOR(1.4d2,  0.d0, 0.d0), &
                VECTOR(1.d9, 0.d0, 0.d0), 1000, thisIteration)
           
        else
           write(filename,'(a, i4.4, a)') "pdr_", thisIteration,".vtk"
           call writeVTKfile(grid, filename, valueTypeString=(/"rho       ",&
                "columnRho ", "UV        ", "uvvec     ", "dust_T    ", "pdrtemp   ", &
                "CO_PDR    ", "C+_PDR    ", "C_PDR     ", "cii_1to0  ", "cii_line  ",&
                "cooling   ", "CII_1     ", "CII_2     ", "CII_3     ", "CII_4     ", &
                "heating   ", "tLast     "/))
        end if
!     endif
     
     if (thisIteration > maxIter) then
        call writeInfo("Number of PDR iterations exceeded maximum, forcing convergence", TRIVIAL)
        
        overallconverged = .true.
     endif
     thisIteration = thisIteration + 1
     print *, "RANK", myrankglobal, "THINKS ", overallconverged, level_conv
     levpop_Iteration = levpop_Iteration + 1
     call MPI_BARRIER(MPI_COMM_WORLD, ier)
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD, ier)
  call writeInfo("Writing temperatures to tLast", TRIVIAL)
  call updateTLast(grid%octreeroot)
  call MPI_BARRIER(MPI_COMM_WORLD, ier)
end subroutine pdr_main_loop


subroutine send_conv_data(nconverged, anynotconverged,overallconverged, level_conv, convcounter)
  integer :: nconverged, convcounter
  logical, intent(in) :: anynotconverged
  logical, intent(inout) :: level_conv
  integer :: status, tag1=50, ierr, tag2=100
  logical, intent(out) :: overallconverged

!  print *, "pre ", overallconverged
!
  call MPI_SEND(nconverged, 1, MPI_INTEGER, 0, tag1, MPI_COMM_WORLD,  ierr)
  call MPI_SEND(anynotconverged, 1, MPI_LOGICAL, 0, tag2, MPI_COMM_WORLD,  ierr)
  call MPI_RECV(level_conv, 1, MPI_LOGICAL, 0, tag1, MPI_COMM_WORLD,  status, ierr)
  call MPI_RECV(overallconverged, 1, MPI_LOGICAL, 0, tag2, MPI_COMM_WORLD,  status, ierr)
  
!  print *, "post ", overallconverged

!  if(.not. level_conv .and. overallconverged) then
!     level_conv = .true.
!     overallconverged =.false.
!  endif

end subroutine send_conv_data


subroutine countNConverged(nconverged, overallconverged, level_conv,convcounter)
  use inputs_mod, only : mindepthAMR, maxdepthAMR

  integer :: nconverged, i, convcounter
  integer :: totalConverged, status, ierr, tag1=50, tag2=100
  logical :: anyfailed, overallconverged, level_conv
  logical :: nfailedarray(1:nhydrothreadsglobal)
  totalConverged = 0
  nConverged = 0
  nfailedarray = .false.

  do i = 1, nhydrothreadsglobal
     call MPI_RECV(nConverged, 1, MPI_INTEGER, i, tag1, MPI_COMM_WORLD,  status, ierr)
     totalConverged = totalConverged + nconverged
     call MPI_RECV(anyfailed, 1, MPI_LOGICAL, i, tag2, MPI_COMM_WORLD,  status, ierr)
     nfailedarray(i) = anyfailed

  enddo



!getting some weird behaviour. The above scheme works but in debug mode it sometimes
!thinks its not converged. I'm putting in this additional catch for now
!  if(minDepthAMR == maxDepthAMR) then
!     if(grid%octreeroot%oned) then
!        if(int(2.d0**(maxdepthamr)) == totalConverged) nfailedArray = .false.
!     endif
!  endif

!  print *, "pre zero ", overallconverged, level_conv

  if(any(nfailedarray) == .true.) then
     if(level_conv) then
        level_conv = .true.
        overallconverged = .false.
     else
        level_conv = .false.
        overallconverged = .false.
     endif
  else
     if(level_conv) then
        level_conv = .true.
        overallconverged = .true.
     else
        level_conv = .true.
        overallconverged = .false.
     endif
  endif

  if(convcounter == 1) overallconverged = .false.

  do i = 1, nhydrothreadsglobal

     call MPI_SEND(level_conv, 1, MPI_LOGICAL, i, tag1, MPI_COMM_WORLD,  ierr)
     call MPI_SEND(overallConverged, 1, MPI_LOGICAL, i, tag2, MPI_COMM_WORLD,  ierr)
  enddo

!  if(level_conv .and. thisConverged) overallconverged = .true.

!  if(convcounter == 2) then
!     print *, "Total cells converged is ", totalconverged 
!  endif
end subroutine countNConverged

recursive subroutine checkConvergence(thisOctal, Notconverged, nconverged, level_conv, convcounter)
  type(octal), pointer :: thisOCtal, child
  integer :: j, subcell, nconverged, convcounter
  logical :: Notconverged, level_conv
!  real(double) :: relChangeCI, relChangeCII, relChangeOI, relChangeC12O
!  integer :: ilevel
!  real(double) :: cii_converged, ci_converged, oi_converged, c12o_converged

  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do j = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(j) == subcell) then
              child => thisoctal%child(j)
              call checkCOnvergence(child, Notconverged, nconverged, level_conv, convcounter)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle

        if(.not. level_conv .and. convcounter == 1) then
           if(thisOctal%level_converged(subcell)) then
              nConverged = nConverged + 1
           else
!              print *, "rank ", myrankglobal, "found a fail"
              !              anynotconverged = .true.
              !              level_conv = .false.
              Notconverged = .true.
           endif
        else
           if(thisOctal%converged(subcell)) then
              nConverged = nConverged + 1
           else
              notconverged = .true.
              !              anynotconverged = .true.
           endif
        endif
     endif
  enddo
end subroutine CheckConvergence


recursive subroutine InitLogicals(thisOctal)
  type(octal), pointer :: thisOCtal, child
  integer :: j, subcell


  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do j = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(j) == subcell) then
              child => thisoctal%child(j)
              call InitLogicals(child)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle
        thisOctal%level_converged(subcell) = .false.
        thisOctal%converged(subcell) = .false.
        thisOctal%expanded(subcell) = .false.
        thisOctal%biChop(subcell) = .false.

     endif
  enddo
end subroutine InitLogicals


recursive subroutine zeroLevelConv(thisOctal)
  type(octal), pointer :: thisOCtal, child
  integer :: j, subcell


  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do j = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(j) == subcell) then
              child => thisoctal%child(j)
              call ZeroLevelConv(child)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle
!        if(.not. thisOctal%converged(subcell)) then
           thisOctal%level_converged(subcell) = .false.
!        endif
!        thisOctal%biChop(subcell) = .false.
!        thisOctal%expanded(subcell) = .false.


!THAW -- BE CAREFUL
        thisOctal%converged(subcell) = .false.
        !        endi



     endif
  enddo
end subroutine ZeroLevelConv

recursive subroutine updateTLast(thisOctal)
  type(octal), pointer :: thisOCtal, child
  integer :: j, subcell


  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do j = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(j) == subcell) then
              child => thisoctal%child(j)
              call updateTLast(child)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle
!        if(.not. thisOctal%converged(subcell)) then
!        thisOctal%level_converged(subcell) = .false.
        thisOctal%tLast(Subcell) = max(thisOctal%tPrev(subcell), 10.0)
!           thisOctal%expanded(subcell) = .false.
!        endif
!        thisOctal%converged(subcell) = .false.

!        thisOctal%biChop(subcell) = .false.

     endif
  enddo
end subroutine UpdateTLast



recursive subroutine calcHeatingOverGrid(thisOctal, reactant, &
     product, alpha, beta, gamma, rate, duplicate, rtmin, rtmax, nelect, &
     nc12o, ncii, nci, noi, level_conv, first_time, grid)   
  use inputs_mod, only : v_turb, mindepthAMR
  type(gridtype) :: grid
  integer :: subcell
  type(octal), pointer :: thisOctal, child
  real(double) :: rate(:), alpha(:), beta(:), gamma(:)
  real(double) :: rtmin(:), rtmax(:)
  integer :: duplicate(:)
  character(len=10) :: product(:,:), reactant(:,:)
  integer :: nc12o, ncii, nci, noi, nelect
  integer :: j
  integer :: NRGR, NRH2, NRHD, NRCO, NRCI, NRSI, NSPEC, nreac, nrays
  logical, intent(in) :: level_conv
  logical, intent(inout) :: first_time
  real(double), parameter :: tmin = 10.d0
!  real(double), parameter :: tmin = 2.7d0
  real(double), parameter :: tmax = 1.d4
  real(double) :: fmean, fratio, dummytemperature, Fcrit, Tdiff
  real(double) :: totalCooling
  integer :: TRACKER, TRACKER2
  logical :: trackMe
  type(vector) :: rvec
  nreac = 329
  nspec = 33
  nrays = nray_func()
  if(thisOCtal%oneD) nrays = 1
!THAW what should fcrit and tdiff be?
  Fcrit=0.005d0
  Tdiff= 0.01d0
!  Tdiff= 0.005d0

    do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do j = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(j) == subcell) then
              child => thisoctal%child(j)
              call calcHeatingOverGrid(child, reactant, &
     product, alpha, beta, gamma, rate, duplicate, rtmin, rtmax, nelect, &
     nc12o, ncii, nci, noi, level_conv, first_time, grid)   
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle
           rVec = subcellCentre(thisOctal, subcell)
           trackMe = .false.
           if (rvec%x <= grid%halfsmallestsubcell) then

              trackMe = .true.
              if(thisOctal%converged(subcell)) then
              else

              endif
           endif

        if(thisOctal%converged(subcell)) cycle
        if(.not.thisOctal%ionfrac(subcell, 2) > 0.99d0) then
           TRACKER = 0
           TRACKER2 = 0



           call CALCULATE_REACTION_RATES(thisOctal%tLast(subcell),thisOctal%Dust_T(subcell), &
                thisOctal%radsurface(subcell, :),thisOctal%AV(subcell, :),thisOctal%thisColRho(subcell, :, :), &
                &REACTANT,PRODUCT,ALPHA,BETA,GAMMA,RATE,RTMIN,RTMAX,DUPLICATE,NSPEC,&
                &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI, nreac, nrays, nc12o, nci, debug=.true.)              

           

           totalCooling = 0.d0
           thisOctal%heatingRate(subcell,:) = 0.d0

           call calc_heating(thisOctal%rho(subcell)/mHydrogen, thisOctal%tLast(subcell), &
                thisOctal%dust_T(subcell), thisOctal%UV(subcell), v_turb, nspec, &
                thisOctal%abundance(subcell, :), nreac, rate, thisOctal%heatingRate(subcell,:), &
                NRGR, NRH2, NRHD, NRCO, NRCI, NRSI, NELECT, nrays)

           
           totalCooling = sum(thisOctal%coolingRate(subcell, :))
           fmean = thisOctal%heatingRate(subcell,12) - totalCooling
           if(totalCooling == 0.d0 .and. thisOctal%heatingRate(subcell,12) == 0.d0) then
              fRatio = 0.d0
              print *, "ALARM"
           else
              fRatio = 2.0D0*abs(Fmean)/abs(thisoCtal%heatingRate(subcell,12) + &
                   totalCooling)
           endif

             ! Store the current temperature in a dummy variable            
             ! Do not start testing the thermal balance until enough iterations have passed for the level populations to begin to converge...
           if (level_conv.and.first_time) then

              TRACKER2 = 1
              dummytemperature = thisOctal%tLast(subcell)
              ! Determine the temperature bracket to begin searching within...
              if (Fmean == 0.d0) then ! Handle the (very rare) case when the initial guess temperature is the correct value    

                 thisOctal%Tlow(subcell) = thisOctal%tLast(subcell)  ! Update the value of Tlow
                 thisOctal%Thigh(subcell) = thisOctal%tLast(subcell) ! Update the value of Thigh

              else if (Fmean > 0.d0) then !---> HEATING                  

                 !                   print *, "heating "
                 thisOctal%Tlow(subcell) = thisOctal%tLast(subcell)  ! Update the value of Tlow 
                !                   thisOctal%tLast(subcell) = 1.3D0*thisOctal%tLast(subcell) !increase 30%
                 thisOctal%tLast(subcell) = 1.3D0*thisOctal%tLow(subcell) !increase 30%
!                 thisOctal%tLast(subcell) = 1.2D0*thisOctal%tLast(subcell) !increase 30%
                 TRACKER = 1
                 !1 for heating, 0 for cooling
                 thisOctal%lastChange(subcell) = 1 !we increase                      

              else if (Fmean < 0.d0) then !---> COOLING                  

                 !                   print *, "reducing temperature"
                 thisOctal%Thigh(subcell) = thisOctal%tLast(subcell) ! Update the value of Thigh 
                 !                 !                   thisOctal%tLast(subcell) = 0.7D0*thisOctal%tLast(subcell) !decrease 30%
                 thisOctal%tLast(subcell) = 0.7D0*thisOctal%thigh(subcell) !decrease 30%
                 !                 thisOctal%tLast(subcell) = 0.8D0*thisOctal%tLast(subcell) !decrease 30%
                 thisOctal%lastChange(subcell) = 0 !we increase                      
                 TRACKER = 2
                 !                   previouschange(pp) = "C" !we decrease                     
                 !                if(myrankglobal == 1) then 

              endif

              thisOctal%tPrev(subcell) = dummytemperature
              
              
           else if (level_conv.and..not.first_time) then
              if(trackMe) print *, "STAGE 4"
              !                print *, "HEATING B", myrankglobal
              dummytemperature = thisOctal%tLast(subcell)
              ! Check for convergence in both the heating-cooling imbalance and the temperature difference between iterations
              if (Fratio < Fcrit) then
                 if(trackMe) print *, "STAGE 5"
!                 print *, "CONVERGED A", Fratio, Fcrit
                 thisOctal%converged(subcell) = .true.
                 !                   print *, "reason  A"
              endif
              if (.not.thisOctal%biChop(subcell)) then
                 if(trackMe) print *, "STAGE 6"
                 !                   print *, "ALPHA"
                 !if we *still* need to heat, increase by 30%
                 if (Fmean > 0.d0.and.thisOctal%lastChange(subcell) == 1) then
                    if(trackMe) print *, "STAGE 7"
                    thisOctal%Tlow(subcell) = thisOctal%tLast(subcell)
                    !thisOctal%tLast(subcell) = 1.2D0*thisOctal%tLast(subcell)
                    !                    thisOctal%tLast(subcell) = 1.3D0*thisOctal%tLast(subcell)
                    thisOctal%tLast(subcell) = 1.3D0*thisOctal%tLow(subcell)
                    !                    thisOctal%tLast(subcell) = 1.1D0*thisOctal%tLast(subcell)
                    thisOctal%Thigh(subcell) = thisOctal%tLast(subcell)
                    
                    thisOctal%lastChange(subcell) = 1
                    !if we *still* need to cool, decrease by 30%                                                                                                                                
                    TRACKER = 3
                 endif
                 
                 if (Fmean < 0.d0.and.thisOctal%lastChange(subcell) == 0) then
                    if(trackMe) print *, "STAGE 8"
                    !                      print *, "BETA"
                    thisOctal%Thigh(subcell) = thisOctal%tLast(subcell)
                    thisOctal%tLast(subcell) = 0.7D0*thisOctal%thigh(subcell)
                    !thisOctal%tLast(subcell) = 0.8D0*thisOctal%tLast(subcell)
                    thisOctal%Tlow(subcell) = thisOctal%tLast(subcell)
                    TRACKER = 4
                    thisOctal%lastChange(subcell) = 0
                 endif
                 
                 
                 !                endif
                 
                 !For all other cases do binary chop and flag the process as .true.
                 !Needs heating but previously it was decreased by 30%.          
                 if (Fmean > 0.d0 .and.thisOctal%lastChange(subcell) == 0) then! .and. thisOctal%tLast(subcell) /= Tmin) then
                    if(trackMe) print *, "STAGE 9"
                    !                      print *, "GAMMA"
                    !     gastemperature(pp) = 1.3D0*gastemperature(pp) !<------       
                    thisOctal%TLast(subcell) = (thisOctal%Thigh(subcell) + &
                         thisOctal%Tlow(subcell))/2.0D0
                    thisOctaL%biChop(subcell)=.true.  !from now on            
                    TRACKER2 = 10
                    !                    TRACKER = 5
                    !                      print *, "NOW DOING BINARY CHOP A"
                    !                      stop
                 endif
                 
                 !Does nothing different to the first one... THaw

                 if (Fmean < 0.d0 .and.thisOctal%lastChange(subcell) == 1) then !.and. thisOctal%tLast(subcell) /= Tmin) then
                    if(trackMe) print *, "STAGE 10"
                    !                      print *, "DELTA"
                    thisOctal%tlast(subcell) = (thisOctal%Thigh(subcell) + &
                         thisOctal%Tlow(subcell))/2.0D0
                    thisOctal%biChop(subcell) = .true.
                    TRACKER2 = 11

                 endif
                 
              else
                 
                 if (Fmean > 0.d0) then
                    if(trackMe) print *, "STAGE 11"
                    !                      print *, "EPSILON"
                    thisOctal%Tlow(subcell) = thisOctal%Tlast(subcell)
                    thisOctal%Tlast(subcell) = (thisOctal%Tlast(subcell) + &
                         thisOctal%Thigh(subcell)) / 2.0D0
!                    thisOctal%Thigh(subcell) = 1.2d0*thisOctal%Tlast(subcell)
                    !                      thisOctal%Thigh(subcell) = 1.3d0*thisOctal%Tlast(subcell)
                    TRACKER2 = 12
                    !                      TRACKER = 2
                 endif
                 
                 if (Fmean < 0.d0) then

                    !                      print *, "PHI"
                    thisOctal%Thigh(subcell) = thisOctal%Tlast(subcell)
                    thisOctal%Tlast(subcell) = (thisOctal%Tlast(subcell) + &
                         thisOctal%Tlow(subcell)) / 2.0D0
                    !                      thisOctal%Tlow(subcell) = 0.7d0*thisOctal%Tlast(subcell)
!                    thisOctal%Tlow(subcell) = 0.8d0*thisOctal%Tlast(subcell)
                    TRACKER2 = 13
                 endif
                 
              endif
              !THAW ONLY CONVERGE ON RATES
              if ((abs(thisOctal%Tlast(subcell)-thisOctal%tPrev(subcell))<=Tdiff).and.&
                   (Fratio>Fcrit)) then
                 !                  print *, "DIFF WAS ", abs(thisOctal%Tlast(subcell)-thisOctal%tPrev(subcell))
                 !                   !                   print *, "KAPPA"
                 if (thisOctal%expanded(subcell)) thisOctal%converged(subcell) = .true.

                 if (Fmean>0.d0) then
                    thisOctal%Tlow(subcell) = thisOctal%Tlast(subcell)
                    thisOctal%Tlast(subcell) = 4.0D0*thisOctal%Tlow(subcell)
                    !                        thisOctal%Tlast(subcell) = 1.5D0*thisOctal%Tlast(subcell)
                    !                         thisOctal%Thigh(subcell) = thisOctal%Tlast(subcell)
                       !                         TRACKER2 = 14
                    !                      thisOctal%Tlast(subcell) = 5.0D0*thisOctal%Tlast(subcell)
                    !           Thigh(pp) = 1.3D0*gastemperature(pp) !*************
                 endif
                 if (Fmean< 0.d0) then
                    TRACKER2 = 15
                    thisOctal%Thigh(subcell) = thisOctal%Tlast(subcell)
                    thisOctal%Tlast(subcell) = 0.25d0*thisOctal%Thigh(subcell)
                    !                        !                         thisOctal%Tlast(subcell) = (1.d0/1.5d0)*thisOctal%Tlast(subcell)
                    !                        thisOctal%Tlow(subcell) = thisOctal%Tlast(subcell)
                    !                        TRACKER = 5
                    !                        !           Tlow(pp) = 0.7D0*gastemperature(pp) !**********    
                    !                     endif
                 endif
                 thisOctal%expanded(subcell) = .true.

              endif
!           endif
              thisOctal%tPrev(subcell) = dummytemperature                
              if(trackMe) print *, "STAGE 18"
              if ((dummytemperature<=Tmin).and.&
                   (Fmean<0.d0)) then
                 if(trackMe) print *, "STAGE 19 converged"
                 !                  print *, "CONVERGED C"
                 !                  print *, "reason  C"
                 thisOctal%converged(subcell) = .true.
              endif
              if ((dummytemperature>=Tmax).and.&
                 (Fmean>0.d0))then
                 if(trackMe) print *, "STAGE 20 converged "
                 !                   print *, "CONVERGED D"
                 !                 print *, "reason  D"
                 thisOctal%converged(subcell) = .true.
                            
              endif

              if (thisOctal%converged(subcell)) then
                 if (thisOctal%tPrev(subcell)<=Tmin) &
                      thisOctal%tPrev(subcell) = Tmin
                 if (thisOctal%tPrev(subcell)>=Tmax) &
                      thisOctal%tPrev(subcell) = Tmax
              endif             
              
           endif
!           print *, "BLA"
           
          else
             if(trackMe) print *, "STAGE 21 converged"
             thisOctal%tLast(subcell) = thisOctal%temperature(subcell)
             thisOctal%converged(subcell) = .true.
          endif

       endif
       !    endif
    enddo
end subroutine calcHeatingOverGrid


recursive subroutine solvePopulations(thisOctal, grid, nelect, ncii, nci, noi, nc12o, level_conv, &
     thisIteration, levpopiteration, coolingOnly)
  use inputs_mod, only : v_turb
  type(gridtype) :: grid
  type(octal), pointer ::  thisOctal, child
  integer :: subcell
  integer :: j, i
  integer :: nrays, thisIteration, levpopiteration
  logical :: level_conv
  integer :: nelect, nproton, nh, nhe, nh2
  
  logical :: coolingOnly

!  integer :: CII_NLEV, CII_NTEMP   !CII cooling variables  
!  
  real(double) :: CII_C_COEFFS(1:5,1:5)
  real(double) :: CI_C_COEFFS(1:5,1:5)
  real(double) :: OI_C_COEFFS(1:5,1:5)
  real(double) :: C12O_C_COEFFS(1:41,1:41)
!  real(double) :: fac
  real(double) :: CI_POP(5), CII_POP(5), OI_POP(5), C12O_POP(41)!, tau_increment, tpop
  type(vector) :: startposition, testposition, rvec, uhat

  integer :: ncontributed
!  logical :: toField

  real(double) :: frac2,  tval, Hplusfrac
  real(double) :: outCii, outC12O, outCI, outOI, AC12O(41, 41), AOI(5, 5), ACI(5, 5), ACII(5, 5)
  logical :: callwrites, firstCell

  real(double) :: popStorageCII(5), popStorageCI(5), popStorageOI(5), popStorageC12O(41)

  Integer, intent(in) :: ncii, nci, noi, nc12o

  real(double) :: taucii(5,5), tauci(5,5), tauoi(5,5), tauc12o(41,41), beta!, beta(41, 41)
  real(double) :: ciisol(1:5), cisol(1:5), oisol(1:5), c12osol(1:41)
!  real(double) :: ngrain, rho_grain, emissivity, BB_ij_dust, S_ij, BB_ij, tmp2, field(41, 41)
!  real(double) :: tmp2small(5, 5), BB_small(5,5), S_small(5,5)
!  real(double) :: 
!  integer, pointer :: epray(:)        !population of evaluation points per ray  

!is gastemperature(pp) thisOCtal%temperature(subcell)?

!unspecified
  integer :: ilevel, jlevel

  nrays = nray_func()
  if(thisOctal%oneD) nrays = 1
  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do j = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(j) == subcell) then
              child => thisoctal%child(j)
              call solvePopulations(child, grid, nelect, ncii, nci, noi, nc12o, level_conv, &
                   thisIteration, levpopiteration, coolingOnly)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle
        if(.not.thisOctal%ionfrac(subcell, 2) > 0.99d0) then  
           if(thisOctal%level_converged(subcell) .or.  thisOctal%converged(subcell)) cycle

           !FIRST CALCULATE COLLISIONAL COEFFICIENTS
           taucii = 0.d0
           tauci = 0.d0
           tauoi = 0.d0
           tauc12o = 0.d0

           NH2 = 31
           NPROTON = 19
           NHE = 26
           NH = 32

           CALL FIND_CCOEFF(CII_NTEMP,CII_NLEV,thisOctal%tLast(subcell),CII_TEMPERATURES,&
                CII_H,CII_HP,CII_EL,CII_HE,CII_H2,CII_PH2,CII_OH2,&
                CII_C_COEFFS,thisOctal%abundance(subcell, NH)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NPROTON)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NELECT)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NHE)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NH2)*thisOctal%rho(subcell)/mhydrogen)

!           if(myrankglobal == 1) then
!           print *, "CII_C_COEFFS ", CII_C_COEFFS
!           print *, "CII_C_INPUTS ", CII_NTEMP,CII_NLEV,thisOctal%tLast(subcell),CII_TEMPERATURES
!           print *, " "
!           print *, CII_H,CII_HP,CII_EL,CII_HE,CII_H2,CII_PH2,CII_OH2
!           print *, " "
!           print *, CII_C_COEFFS,thisOctal%abundance(subcell, NH)*thisOctal%rho(subcell)/mhydrogen, &
!                thisOctal%abundance(subcell, NPROTON)*thisOctal%rho(subcell)/mhydrogen, &
!                thisOctal%abundance(subcell, NELECT)*thisOctal%rho(subcell)/mhydrogen, &
!                thisOctal%abundance(subcell, NHE)*thisOctal%rho(subcell)/mhydrogen, &
!!                thisOctal%abundance(subcell, NH2)*thisOctal%rho(subcell)/mhydrogen!
!
!           stop
!        endif
           !CI 
           CALL FIND_CCOEFF(CI_NTEMP,CI_NLEV,thisOctal%tLast(subcell),CI_TEMPERATURES,&
                CI_H,CI_HP,CI_EL,CI_HE,CI_H2,CI_PH2,CI_OH2,&
                CI_C_COEFFS,thisOctal%abundance(subcell, NH)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NPROTON)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NELECT)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NHE)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NH2)*thisOctal%rho(subcell)/mhydrogen)
!
!
!           !OI
           CALL FIND_CCOEFF(OI_NTEMP,OI_NLEV,thisOctal%tLast(subcell),OI_TEMPERATURES,&
                OI_H,OI_HP,OI_EL,OI_HE,OI_H2,OI_PH2,OI_OH2,&
                OI_C_COEFFS,thisOctal%abundance(subcell, NH)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NPROTON)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NELECT)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NHE)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NH2)*thisOctal%rho(subcell)/mhydrogen)
!
!           !CO
           CALL FIND_CCOEFF(C12O_NTEMP,C12O_NLEV,thisOctal%tLast(subcell),C12O_TEMPERATURES,&
                C12O_H,C12O_HP,C12O_EL,C12O_HE,C12O_H2,C12O_PH2,C12O_OH2,&
                C12O_C_COEFFS,thisOctal%abundance(subcell, NH)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NPROTON)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NELECT)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NHE)*thisOctal%rho(subcell)/mhydrogen, &
                thisOctal%abundance(subcell, NH2)*thisOctal%rho(subcell)/mhydrogen)



           

           !           !CALCULATE THE SOURCE FUNCTION
           
           !NOW CAST RAYS TO GET TAUS
           nrays = nray_func()
           if(thisOctal%oneD) nrays = 1
           
           startPosition = subcellCentre(thisOctal, subcell)
           testPosition = startPosition
           thisOctal%coolingRate(subcell,:) = 0.d0           
           ncontributed = 0         
           frac2=1.0D0/sqrt(8.0*KB*thisOctal%tLast(subcell)/PI/MP + v_turb**2)

!           print *, "frac2 ", frac2
!           print *, "vturb ", v_turb
!           print *, "thisOctal%tLast(subcell)", thisOctal%tLast(subcell)
!           print *, "x", subcellCentre(thisOctal, subcell)!

!           stop

           do i = 1, nrays              
              rVec = VECTOR(vectors(1, i-1), vectors(2, i-1), vectors(3, i-1))
              if(thisOctal%oneD) rVec = VECTOR(-1.d0, 0.d0, 0.d0) 

              firstCell = .true.
!              if(thisOctal%oneD) then
!                 if(abs(rVec%y) > 0.d0 .or. abs(rVec%z) > 0.d0 .and. rVec%x > 0.d0) then
!                    taucii = 1.d30
!                    tauci = 1.d30
!                    tauoi = 1.d30
!                    tauc12o = 1.d30
!                    goto 666
!                 endif
!              endif
           
              uhat = rVec
              call normalize(uhat)
              

              

              testPosition = startPosition                        
!              firstCell = .true.
              do while (inOctal(grid%octreeRoot, testPosition))                 

                 tval = 0.d0
                 
                 call getRayTracingValuesPDR_TWO(grid, testposition, uHat, CII_POP, CI_POP, OI_POP, C12O_POP, tVal,&
                      HplusFrac) 
                 
!                 print *, "CII_POP", cii_pop
!                 stop
                 if(Hplusfrac > 0.99d0) then                       
                    exit
                 end if
                 firstcell = .false.                 

                 !!calculate the optical depth along each ray for each species                 
                 do ilevel = 1, C12O_NLEV
                    do jlevel = 1, C12O_NLEV
 !                do ilevel = 1, 2
 !                   do jlevel = 1, 1
                       if(jlevel >= ilevel) exit
                       
                       if(ilevel <= CII_NLEV .and. jlevel <= cii_nlev) then
!                          print *, "i, j ", ilevel, jlevel, i
!                          if()
                          if(thisOctal%cii_pop(subcell, ilevel) /= 0.d0) then
                             taucii(ilevel, jlevel) = taucii(ilevel,jlevel) + updateTau(CII_A_COEFFS(ilevel,jlevel), &
                                  cii_frequencies(ilevel,jlevel), CII_POP(jlevel), CII_WEIGHTS(ilevel), CII_POP(ilevel), &
                                  CII_WEIGHTS(jlevel), frac2, tval, popStorageCII(jlevel), popStorageCII(ilevel))
                          endif
                          if(thisOctal%ci_pop(subcell, ilevel) /= 0.d0) then
                             tauci(ilevel, jlevel) = tauci(ilevel,jlevel) + updateTau(CI_A_COEFFS(ilevel,jlevel), &
                                  ci_frequencies(ilevel,jlevel), CI_POP(jlevel), CI_WEIGHTS(ilevel), CI_POP(ilevel), &
                                  CI_WEIGHTS(jlevel), frac2, tval, popStorageCI(jlevel), popStorageCI(ilevel))
                          endif
                          if(thisOctal%oi_pop(subcell, ilevel) /= 0.d0) then
                             tauoi(ilevel, jlevel) = tauoi(ilevel,jlevel) + updateTau(OI_A_COEFFS(ilevel,jlevel), &
                                  oi_frequencies(ilevel,jlevel), OI_POP(jlevel), OI_WEIGHTS(ilevel), OI_POP(ilevel), &
                                  OI_WEIGHTS(jlevel), frac2, tval, popStorageOI(jlevel), popStorageOI(ilevel))
                          endif
                       endif
                       if(thisOctal%c12o_pop(subcell, ilevel) /= 0.d0) then
                          tauc12o(ilevel, jlevel) = tauc12o(ilevel,jlevel) + updateTau(c12o_A_COEFFS(ilevel,jlevel), &
                               c12o_frequencies(ilevel,jlevel), c12o_POP(jlevel), c12o_WEIGHTS(ilevel), c12o_POP(ilevel), &
                               c12o_WEIGHTS(jlevel), frac2, tval, popStoragec12o(jlevel), popStoragec12o(ilevel))
                       endif
                       
                    end do
                 enddo
                 popStorageCII(:) = cii_pop(:)
                 popStorageCI(:) = ci_pop(:)
                 popStorageOI(:) = oi_pop(:)
                 popStorageC12o(:) = c12o_pop(:)                                  
                 
                 testPosition = testPosition + ((tVal+1.d-10*grid%halfsmallestsubcell)*uhat)              
              end do
           enddo
 
           taucii = taucii/dble(nrays)
           
           tauci = tauci/dble(nrays)
           tauoi = tauoi/dble(nrays)
           tauc12o = tauc12o/dble(nrays)
           
           !              print *, "taucii ", taucii(2,1)
!                         print *, "tauci ", tauci(2,1)
           !              print *, "tauoi ", tauoi(2,1)
!           print *, "taucii ", thisOctal%AV(subcell, 1), taucii(2,1)
!           stop

           !              taucii = taucii*1.d10
           !              tauci = tauci*1.d10
           !              tauoi = tauoi*1.d10
           !              tauc12o = tauc12o*1.d10
           
           !              print *, taucii(2, 1)
           !              print *, tauci(2, 1)
           !              print *, tauoi(2, 1)
           !              print *, tauc12o(2, 1)
           
           
666        continue
           
           !t              print *, "tauc12o", tauc12o
           thisOctal%CIItransition(subcell, :,:) = 0.d0
           thisOctal%CItransition(subcell, :,:) = 0.d0
           thisOctal%OItransition(subcell, :,:) = 0.d0
           thisOctal%C12Otransition(subcell, :,:) = 0.d0

           
           call setupEscapeParameters2(cii_FREQUENCIES, thisOctal%Dust_T(subcell), &
                thisOctal%CII_pop(subcell, :), CII_WEIGHTS, beta, CII_A_COEFFS, taucii, &
                thisOctal, subcell, CII_B_COEFFS, CII_C_COEFFS,1)
           
           call setupEscapeParameters2(ci_FREQUENCIES, thisOctal%Dust_T(subcell), &
                thisOctal%CI_pop(subcell, :), CI_WEIGHTS, beta, CI_A_COEFFS, tauci, &
                   thisOctal, subcell, CI_B_COEFFS, CI_C_COEFFS,2)
           
           call setupEscapeParameters2(oi_FREQUENCIES, thisOctal%Dust_T(subcell), &
                   thisOctal%oI_pop(subcell, :), oI_WEIGHTS, beta, oI_A_COEFFS, tauoi, &
                   thisOctal, subcell, oI_B_COEFFS, oI_C_COEFFS,3)
              
           call setupEscapeParameters2(c12o_FREQUENCIES, thisOctal%Dust_T(subcell), &
                thisOctal%C12o_pop(subcell, :), C12o_WEIGHTS, beta, C12o_A_COEFFS, tauc12o, &
                thisOctal, subcell, C12o_B_COEFFS, C12o_C_COEFFS,4)
           
!           if(.not. coolingOnly) then
              
!        enddo
!     enddo
!           print *, "BLA"

  
           ACI = 0.d0
           ACII = 0.d0
           AOI = 0.d0
           AC12O = 0.d0
!           ACI = 1.d-30
!           ACII = 1.d-30
!           AOI = 1.d-30
!           AC12O = 1.d-30
           !        ACI = 1.d-8
           !        ACII = 1.d-8
           !        AOI = 1.d-8
           !        AC12O = 1.d-8
            !do final population step
!           print *, "AOI D", AOI

           do ilevel = 1, C12O_NLEV
              outc12o = 0.d0
              outci = 0.d0
              outoi = 0.d0
              outcii = 0.d0

              do jlevel = 1, C12O_NLEV
!                 if(jlevel >= ilevel) exit              
                 
                 if(ilevel <= CII_NLEV .and. jlevel <= CII_NLEV) then
                    outcii = outcii + thisOctal%CIItransition(subcell, ilevel, jlevel)
                    ACII(ilevel, jlevel) = thisOctal%CIItransition(subcell, jlevel, ilevel)
                    
                    outci = outci + thisOctal%CItransition(subcell, ilevel, jlevel)
                    ACI(ilevel, jlevel) = thisOctal%CItransition(subcell, jlevel, ilevel)
                    
                    outoi = outoi + thisOctal%OItransition(subcell, ilevel, jlevel)
                    AOI(ilevel, jlevel) = thisOctal%OItransition(subcell, jlevel, ilevel)

                 endif
                 outc12o = outc12o + thisOctal%C12otransition(subcell, ilevel, jlevel)
                 AC12O(ilevel, jlevel) = thisOctal%C12otransition(subcell, jlevel, ilevel)


              enddo
              
              if(ilevel <= CII_NLEV) then
                 ACII(ilevel, ilevel) = -outCII
                 ACI(ilevel, ilevel) = -outCI
                 AOI(ilevel, ilevel) = -outOI
              endif

              AC12O(ilevel, ilevel) = -outC12O
!              if(ilevel == 38 .and. jlevel ==!! 38) then
!                 print *, "38 stuff !"
!                 print *, AC12O(ilev!el, ilevel)
!                 
!              endif
           enddo

!           print *, "AOI A", AOI
!
           
           DO Ilevel=1,C12O_NLEV

              if(ilevel <= CII_NLEV) then
!not sure if these should be zeroed
!THEY SHOULD BE - THaw
!                 thisOctal%cii_pop(subcell,iLevel)=0.0D0
!                 thisOctal%ci_pop(subcell,iLevel)=0.0D0
!                 thisOctal%oi_pop(subcell,iLevel)=0.0D0
                 ciisol(ilevel) = 0.d0
                 cisol(ilevel) = 0.d0
                 oisol(ilevel) = 0.d0
                 ACII(CII_NLEV,Ilevel)=1.0D-8
                 ACI(CI_NLEV,Ilevel)=1.0D-8            
                 AOI(OI_NLEV,Ilevel)=1.0D-8            
              endif!
              c12osol(ilevel) = 0.d0
!              thisOctal%c12o_pop(subcell,iLevel)=0.0D0
              AC12O(C12O_NLEV,Ilevel)=1.0D-8            
           enddo

           

           ciisol(5)=1.0D-8*thisOctal%rho(subcell)*thisOctal%abundance(subcell,ncii)/mhydrogen
           cisol(5)=1.0D-8*thisOctal%rho(subcell)*thisOctal%abundance(subcell,nci)/mhydrogen
           oisol(5)=1.0D-8*thisOctal%rho(subcell)*thisOctal%abundance(subcell,noi)/mhydrogen
           c12osol(41)=1.0D-8*thisOctal%rho(subcell)*thisOctal%abundance(subcell,nc12o)/mhydrogen

!           thisOctal%cii_pop(subcell,5)=1.0D-8*thisOctal%rho(subcell)*thisOctal%abundance(subcell,ncii)/mhydrogen
!           thisOctal%ci_pop(subcell,5)=1.0D-8*thisOctal%rho(subcell)*thisOctal%abundance(subcell,nci)/mhydrogen
!           thisOctal%oi_pop(subcell,5)=1.0D-8*thisOctal%rho(subcell)*thisOctal%abundance(subcell,noi)/mhydrogen
!           thisOctal%c12o_pop(subcell,41)=1.0D-8*thisOctal%rho(subcell)*thisOctal%abundance(subcell,nc12o)/mhydrogen

!           print *, "USING NCI = ", NCI
!           print *, "ACI ", ACI
 !          print *, "ACII ", ACII
 !          print *, "AOI ", AOI
  !         print *, "AC12O ", AC12O
!
!           do ilevel = 1, cii_nlev
!              do jlevel = 1, cii_nlev
!                 print *, "CII ", ilevel, jlevel, ACII(ilevel, jlevel)
!                 print *, "CI ", ilevel, jlevel, ACI(ilevel, jlevel)
!                 print *, "CIIt ", ilevel, jlevel, thisOctal%CIItransition(subcell, ilevel, jlevel)
!                 print *, "CIt ", ilevel, jlevel, thisOctal%CItransition(subcell, ilevel, jlevel)
!
!              enddo
!           enddo


  !         print *, "rank ", myrankglobal, "doing cii gauss-jordan"
!           CALL GAUSS_JORDAN(ACII,CII_NLEV,CII_NLEV,thisOctal%CII_POP(subcell, :),1, callwrites)
           CALL GAUSS_JORDAN(ACII,CII_NLEV,CII_NLEV, ciisol(:),1, callwrites)
   !        print *, "rank ", myrankglobal, "done"
           
!           print *, "rank ", myrankglobal, "doing ci gauss-jordan"
!           CALL GAUSS_JORDAN(ACI,CI_NLEV,CI_NLEV,thisOctal%CI_POP(subcell, :),2, callwrites)
           CALL GAUSS_JORDAN(ACI,CI_NLEV,CI_NLEV, cisol(:),2, callwrites)
!           print *, "rank ", myrankglobal, "done"
           !
!           print *, "rank ", myrankglobal, "doing oi gauss-jordan"
!           CALL GAUSS_JORDAN(AOI,OI_NLEV,OI_NLEV,thisOctal%OI_POP(subcell, :),3, callwrites)
           CALL GAUSS_JORDAN(AOI,OI_NLEV,OI_NLEV,oisol(:),3, callwrites)
!           print *, "rank ", myrankglobal, "done"
           
!           print *, "rank ", myrankglobal, "doing c12o gauss-jordan"
 !          CALL GAUSS_JORDAN(AC12O,C12O_NLEV,C12O_NLEV,thisOctal%C12O_POP(subcell, :),4, callwrites)
           CALL GAUSS_JORDAN(AC12O,C12O_NLEV,C12O_NLEV,c12osol(:),4, callwrites)
!           print *, "rank ", myrankglobal, "done"

         

           
           do ilevel = 1, C12O_NLEV
              if(Ilevel <= CII_NLEV) then

                 if (ciisol(ilevel).lt.1.D-35) then
                    ciisol(ilevel)=0.0D0!1.0D-99!then !stop 'found negative solution!'
                 endif
                 if (cisol(ilevel).lt.1.0D-35) then
                    cisol(ilevel)=0.0D0!1.0D-99!then !stop 'found negative solution!'
                 endif
                 
                 if (oisol(ilevel).lt.1.0D-35) then
                    oisol(ilevel)=0.0D0!1.0D-99!then !stop 'found negative solution!'             
                 endif
              endif

              if (c12osol(ilevel).lt.1.0D-35) then
                 C12Osol(ilevel)=0.0D0!1.0D-99!then !stop 'found negative solution!' 
              endif
              

           ENDDO

!           do ilevel = 1, C12O_NLEV
!              if(Ilevel <= CII_NLEV) then
!
!                 if (ciisol(ilevel).lt.0.0D0) then
!                    ciisol(ilevel)=0.0D0!1.0D-99!then !stop 'found negative solution!'
!                 endif
!                 if (cisol(ilevel).lt.0.0D0) then
!                    cisol(ilevel)=0.0D0!1.0D-99!then !stop 'found negative solution!'
!                 endif
!!                 
!                 if (oisol(ilevel).lt.0.0D0) then
!                    oisol(ilevel)=0.0D0!1.0D-99!then !stop 'found negative solution!'             
!                 endif
!              endif!
!
!              if (c12osol(ilevel).lt.0.0D0) then
!                 C12Osol(ilevel)=0.0D0!1.0D-99!then !stop 'found negative solution!' 
!              endif
!              
!
!           ENDDO


!THAW - USE THIS TO FORCE CONVERGENCE
           if(.not. thisOctal%oned) then
              if (levpopIteration.ge.120) then
!                 if(CIIsol(:)
                 C12Osol(:)=thisOctal%C12O_pop(subcell, :)
                 CIIsol(:)=thisOctal%CII_pop(subcell, :)
                 CIsol(:)=thisOctal%CI_pop(subcell, :)
                 OIsol(:)=thisOctal%OI_pop(subcell, :)

              elseif (levpopIteration.ge.75) then
                 CIIsol(:)=0.5*(CIIsol(:) + thisOctal%CII_pop(subcell,:))
                 CIsol(:)=0.5*(CIsol(:) + thisOctal%CI_pop(subcell,:))
                 OIsol(:)=0.5*(OIsol(:) + thisOctal%OI_pop(subcell,:))
                 C12Osol(:)=0.5*(C12Osol(:) + thisOctal%C12O_pop(subcell,:))
              endif
           endif
           thisOctal%relch(subcell, :, :) = 0.d0
           thisOctal%level_converged(subcell) = .true.
           do ilevel = 1, C12O_NLEV
              if(Ilevel <= CII_NLEV) then
                 if(ciisol(ilevel) >= thisOctal%abundance(subcell, ncii)*1.d-10) then
                    if(ciisol(ilevel) == 0.d0 .and. thisOctal%cii_pop(subcell, ilevel) == 0.d0) then
                       thisOctal%relch(subcell, 1, ilevel) = 0.d0
                    else
                       thisOctal%relch(subcell, 1, ilevel) = 2.0D0*ABS((CIIsol(ilevel)&
                            -thisOctal%CII_pop(subcell, ilevel))/(CIIsol(ilevel)+&
                            thisOCtal%CII_pop(subcell, ilevel)))
                    endif
                 endif
                 if(cisol(ilevel) >= thisOctal%abundance(subcell, nci)*1.d-10) then
                    if(cisol(ilevel) == 0.d0 .and. thisOctal%ci_pop(subcell, ilevel) == 0.d0) then
                       thisOctal%relch(subcell, 2, ilevel) = 0.d0
                    else
                       thisOctal%relch(subcell, 2, ilevel) = 2.0D0*ABS((Cisol(ilevel)&
                            -thisOctal%CI_pop(subcell, ilevel))/(CIsol(ilevel)+&
                            thisOCtal%CI_pop(subcell, ilevel)))
                    endif
                 endif
                 if(oisol(ilevel) >= thisOctal%abundance(subcell, noi)*1.d-10) then
                    if(oisol(ilevel) == 0.d0 .and. thisOctal%oi_pop(subcell, ilevel) == 0.d0) then
                       thisOctal%relch(subcell, 3, ilevel) = 0.d0
                    else
                       thisOctal%relch(subcell, 3, ilevel) = 2.0D0*ABS((oIsol(ilevel)&
                            -thisOctal%oI_pop(subcell, ilevel))/(oIsol(ilevel)+&
                            thisOCtal%oI_pop(subcell, ilevel)))
                    endif
                 endif
              else
                 thisOctal%relch(subcell, 2, ilevel) = 0.d0
                 thisOctal%relch(subcell, 3, ilevel) = 0.d0
                 thisOctal%relch(subcell, 1, ilevel) = 0.d0
              end if
              if(c12osol(ilevel) >= thisOctal%abundance(subcell, nc12o)*1.d-10) then
                 if(c12osol(ilevel) == 0.d0 .and. thisOctal%c12o_pop(subcell, ilevel) == 0.d0) then
                    thisOctal%relch(subcell, 4, ilevel) = 0.d0
                 else
                    thisOctal%relch(subcell, 4, ilevel) = 2.0D0*ABS((C12osol(ilevel)&
                         -thisOctal%C12o_pop(subcell, ilevel))/(C12osol(ilevel)+&
                         thisOCtal%C12o_pop(subcell, ilevel)))
                 endif
              endif

              
              !              level_conv = .true.
              !              thicOctal%level_converged(subcell) = .true.
           !   thisOctal%level_converged(subcell) = .true.                 
              !              print *, "checking relch", thisOctal%relch(subcell, 1, ilevel)

!              print *, "relch", thisOctal%relch(subcell, :, ilevel)
              if(ilevel <= CII_NLEV) then
                 if(thisOctal%relch(subcell, 1, ilevel) > 1.d-2 .or. &
                      thisOctal%relch(subcell, 2, ilevel) > 1.d-2 .or. &
                      thisOctal%relch(subcell, 3, ilevel) > 1.d-2 .or. &
                      thisOctal%relch(subcell, 4, ilevel) > 1.d-2) then
                    thisOctal%level_converged(subcell) = .false.                 
                    !                    print *, "thisOctal%CII_pop(subcell, ilevel)", thisOctal%CII_pop(subcell, ilevel)
!                    print *, "cii_sol ", CIIsol(ilevel)                 
!                    print *, "ilevel ",ilevel
!                    print *, thisOctal%relch(subcell, :, ilevel)
!                    print *, "GOT A FAIL "
                    !                    exit
                    !                level_conve = .false.
                 endif
              else
                 if(thisOctal%relch(subcell, 4, ilevel) > 1.d-2) then
                    thisOctal%level_converged(subcell) = .false.      
                 endif
                 
              endif
           enddo
           !        enddo
           thisOctal%cii_pop(subcell, :) = ciisol(:)
           thisOctal%ci_pop(subcell, :) = cisol(:)
           thisOctal%c12o_pop(subcell, :) = c12osol(:)
           thisOctal%oi_pop(subcell, :) = oisol(:)           

!        endif
        
        else
           thisOctal%ciiline(subcell, :, :) = 0.d0
           thisOctal%ciiTransition(subcell, :, :) = 0.d0
           thisOctal%ciline(subcell, :, :) = 0.d0
           thisOctal%ciTransition(subcell, :, :) = 0.d0
           thisOctal%oiline(subcell, :, :) = 0.d0
           thisOctal%oiTransition(subcell, :, :) = 0.d0
           thisOctal%c12oline(subcell, :, :) = 0.d0
           thisOctal%c12oTransition(subcell, :, :) = 0.d0
           thisOctal%coolingRate(subcell,:) = 0.d0
           thisOctal%CII_POP(subcell, :) = 0.d0
           thisOctal%CI_POP(subcell, :) = 0.d0
           thisOctal%OI_POP(subcell, :) = 0.d0
           thisOctal%C12O_POP(subcell, :) = 0.d0
           thisOctal%level_converged(subcell) = .true.
        end if
               
     end if
     
  enddo
end subroutine solvePopulations


real(double) function updateTau(Acoeff,freq, jPop, iWeight, iPop, &
     jWeight, frac2, tval, jpop_m1, ipop_m1)

  real(double) :: frac2, tval, frac1, frac3, ipop, iweight, jpop, jweight, freq, Acoeff
  real(double) :: jpop_m1, ipop_m1

  frac1=(Acoeff*(cspeed**3))/(8.0*pi*(freq**3))
!  print *, "frac1", frac1
!  print *, "Acoeff", Acoeff
!  print *, "freq", freq

!  stop
!  frac3=((jPop*iWeight-iPop*jWeight)+&
!       &(jPop*iweight-iPop*jWeight))&
!       /2./jWeight                       
  jpop_m1 = jpop_m1
  ipop_m1 = ipop_m1

!THAW - I AM NOT CONVINCED AS TO WHAT FRAC3 SHOULD BE...!
!  frac3 = (((jpop_m1 + jpop)*iweight)/(2.d0*jweight) - ((ipop_m1 + ipop)/2.d0))
!  frac3 = (((jpop_m1 + jpop)*iweight)/(2.d0*jweight) - ((ipop_m1 + ipop)/2.d0))
!  print *, "frac3 A", frac3
!  frac3 = (((jpop + jpop)*iweight)/(2.d0*jweight) - ((ipop + ipop)/2.d0))


  frac3 = (jpop*iweight - ipop*jweight)/jweight
!  print *, "frac3 B", frac3
!  print *, "ipop", ipop
!  print *, "jpop", jpop
!  print *, "jweight", jweight
!  print *, "iweight", iweight
  

!  print *, "frac1, ",frac1
!  print *, "frac2, ",frac2

  !  frac3 = (((jpop_m1 + jpop)*iweight)/(2.d0*jweight) - ((ipop_m1 + ipop)/2.d0))
  if(frac3 < 0.d0) frac3 = 0.d0
  !FRAC 3 is currently wrong

!THAW temporarily rming the pcToCm
! updateTau=frac1*frac2*frac3*tval*1.d10*pcToCm
! updateTau=frac1*frac2*frac3*tval*1.d10!/pctocm

 updateTau=frac1*frac2*frac3*tval*1.d10!*pctocm
! print *, "dTau, ",updateTau
! updateTau=frac1*frac2*frac3*tval!*pctocm
! updateTau=frac1*frac2*frac3*tval!*pctocm
!  updateTau=frac1*frac2*frac3*tval*pctocm

end function updateTau


!subroutine setupEscapeParameters2(freq, tDust, &
!     iPop, jPop, iweight, &
!     jweight, beta,  Acoeff, tau, thisOctal, subcell, ilevel, jlevel, &
!     BCoeffs, Ccoeffs, id)
subroutine setupEscapeParameters2(freq, tDust, &
     pop, weight, beta,  Acoeff, tau, thisOctal, subcell,  &
     BCoeffs, Ccoeffs, id)

  real(double) :: freq(:,:), tDust, pop(:), weight(:), Acoeff(:,:), beta
  real(double) :: tmp2, fac, ngrain, rho_grain, emissivity, BB_ij_dust, tpop
  real(double) :: tau(:,:)
  logical :: tofield
  logical :: fullBeta
  integer :: id
  type(octal), pointer :: thisOctal
  integer :: subcell, ilevel, jlevel
  real(double) :: BCoeffs(:,:), CCoeffs(:,:), BB, sourceFun, fac2
  real(double), allocatable :: field(:,:)
  integer :: nlev
  if(id == 1 .or. id ==2 .or. id == 3) then
     nlev = 5
  else
     nlev = 41
  endif
  allocate(field(1:nlev, 1:nlev))
  field = 0.d0

  thisOctal%coolingRate(subcell, id) = 0.d0

do ilevel = 1, nlev
   do jlevel = 1, nlev
      if(jlevel >= ilevel) exit
      tofield = .false.
      fullbeta = .true.
      TMP2=2.0D0*HP*(freq(ilevel, jlevel)**3)/(cspeed**2)
           
      !Planck function !2.7D0 is the CMBR temperature                                                        
      
      fac = HP*freq(ilevel, jlevel)/(KB*(2.7D0))
!      print *, "fac is ", fac
      if(fac < 50.d0 .and. fac /= 0.d0) then
         BB = TMP2*(1.0D0/(EXP(fac)-1.0D0))
      else
         BB = 0.d0
      endif
      !#ifdef DUST                                                                                                    
      NGRAIN=2.0D-9 !2.0D-12*densityofgas depth depented                                                     
      rho_grain=2.0D0
      EMISSIVITY=(RHO_GRAIN*NGRAIN)*(0.01*(1.3*FREQ(ilevel, jlevel)/3.0D11))
      fac2 = HP*freq(ilevel, jlevel)/KB/tDust
   !   print *, "fac2 is ", fac2
      if(fac2 < 50.d0 .and. fac2 /= 0.d0) then
         BB_ij_dust = TMP2*(1.0D0/(EXP(fac2)&
              -1.D0)*EMISSIVITY)
      else
         BB_ij_dust = 0.d0
      endif
      BB = BB + BB_ij_dust
      !#endif                                                                                                         
      if (Pop(ilevel).eq.0.d0 .or. BB == 0.d0) then
         sourceFun=0.0D0
         beta=1.0D0
         toField= .true.
         !     print *, "GOING TO FIELD", ipop
      endif
      
      if(.not. toField) then
         TPOP = 0.d0
         
         TPOP=(pop(jlevel)*weight(ilevel))/&
              (pop(ilevel)*Weight(jlevel))-1.0D0
         
         IF(abs(TPOP).lt.1.0D-50) then ! .or. thisOctal%CII_pop(subcell, ilevel) .lt. 1.d-50) then
            sourceFun=HP*FREQ(ilevel, jlevel)*&
                 Pop(ilevel)*ACOEFF(ilevel, jlevel)/4./pi
            beta=1.0D0
            !        goto 1
            fullBeta = .false.
         else
            fullBeta = .true.
            !calculation of source function (taken from UCL_PDR)                                                   
            sourceFun=TMP2/TPOP
         endif

!         print *, "sourcefun ", sourcefun
!         print *, "BB ", BB
         !     if(debug) then
         !!        print *, "sourceFun is ", sourceFun
         !        dummy = abs(sourceFun)!
         !        if(sourceFun == 0.d0) th!en
         !           print *, "zero source! function", sourcefun, tmp2, tpop
         !        endif
         !     endif
         
         if(fullBeta) then
            if(tau(ilevel, jlevel) < -5.d0) then
               beta=(1.0D0-EXP(5.0D0))/(-5.0D0)
            else if (abs(tau(ilevel, jlevel)).lt.1.0D-8) then
               beta=1.0D0
            else
               beta = (1.d0 - exp(-tau(ilevel, jlevel)))/tau(ilevel, jlevel)
            end if
         endif
      endif
      
      if(.not. toField) then
         if(sourceFun /= 0.d0) then
            if(id == 1) then
               thisOctal%ciiline(subcell,ilevel,jlevel) = Acoeff(ilevel, jlevel)*HP*freq(ilevel, jlevel) * &
                    & Pop(ilevel)*beta*(sourceFun-BB)/sourceFun
               if(thisOctal%ciiline(subcell,ilevel,jlevel) < 0.d0) then
                  print *, "line is less than zero"
                  thisOctal%ciiline(subcell,ilevel,jlevel) = 0.d0
               endif
            else if(id == 2) then
               thisOctal%ciline(subcell,ilevel,jlevel) = Acoeff(ilevel, jlevel)*HP*freq(ilevel, jlevel) * &
                    & Pop(ilevel)*beta*(sourceFun-BB)/sourceFun
               if(thisOctal%ciline(subcell,ilevel,jlevel) < 0.d0) thisOctal%ciline(subcell,ilevel,jlevel) = 0.d0
            else if(id == 3) then
               thisOctal%oiline(subcell,ilevel,jlevel) = Acoeff(ilevel, jlevel)*HP*freq(ilevel, jlevel) * &
                    & Pop(ilevel)*beta*(sourceFun-BB)/sourceFun
               if(thisOctal%oiline(subcell,ilevel,jlevel) < 0.d0) thisOctal%oiline(subcell,ilevel,jlevel) = 0.d0
            else if(id == 4) then
               thisOctal%c12oline(subcell,ilevel,jlevel) = Acoeff(ilevel, jlevel)*HP*freq(ilevel, jlevel) * &
                    & Pop(ilevel)*beta*(sourceFun-BB)/sourceFun
               if(thisOctal%c12oline(subcell,ilevel,jlevel) < 0.d0) thisOctal%c12oline(subcell,ilevel,jlevel) = 0.d0
            else
               call torus_abort("species ID not recognised in pdr_mod")
            endif
         else
            print *, "source function is zero"
            if(id == 1) then
               thisOctal%ciiline(subcell,ilevel,jlevel) = 0.d0
            else if(id == 2) then
               thisOctal%ciline(subcell,ilevel,jlevel) = 0.d0
            else if(id == 3) then
               thisOctal%oiline(subcell,ilevel,jlevel) = 0.d0
            else if(id == 4) then
               thisOctal%c12oline(subcell,ilevel,jlevel) = 0.d0
            else
               call torus_abort("species ID not recognised in pdr_mod")
            endif
            
         endif
         if(id == 1) then
            thisOctal%coolingRate(subcell,1) = thisOctal%coolingRate(subcell,1) + thisOctal%ciiline(subcell, ilevel,jlevel)        
         else if(id == 2) then
            thisOctal%coolingRate(subcell,2) = thisOctal%coolingRate(subcell,2) + thisOctal%ciline(subcell, ilevel,jlevel)        
         else if(id == 3) then
            thisOctal%coolingRate(subcell,3) = thisOctal%coolingRate(subcell,3) + thisOctal%oiline(subcell, ilevel,jlevel)        
         else if(id == 4) then
            thisOctal%coolingRate(subcell,4) = thisOctal%coolingRate(subcell,4) + thisOctal%c12oline(subcell, ilevel,jlevel)   
         endif
      endif
      
  
      if(thisoctal%coolingrate(subcell,1) == 0.d0 .and. id == 1 .and. .not. toField) then
         print *, "id", id
         print *, "cooling rate is zero ", thisOctal%ciiline(subcell, ilevel,jlevel)
         print *, "sourcefun, bb ", sourcefun, bb
         print *, "fac1, fac2 ", fac, fac2
         print *, "ipop", pop(ilevel)
 !        stop
      endif
      
      field(ilevel, jlevel) = (1.0D0-beta)*sourceFun + beta*BB
      field(jlevel, ilevel) = field(ilevel, jlevel)
      
   enddo
enddo

do ilevel = 1, nlev
   do jlevel = 1, nlev
      
      if(id == 1) then
         thisOctal%ciiTRANSITION(subcell, ilevel,jlevel)=Acoeff(ilevel, jlevel)&
              +Bcoeffs(ilevel, jlevel)*FIELD(ilevel, jlevel)&
              +Ccoeffs(ilevel, jlevel)
!         thisOctal%ciiTRANSITION(subcell, jlevel,ilevel)=Acoeff&
!              +Bcoeffs*(-FIELD) &
!              +Ccoeffs
         
         IF(ABS(thisOctal%ciiTRANSITION(subcell, ilevel,jlevel)).LT.1.0D-50) thisOctal%ciiTRANSITION(subcell,ilevel,jlevel)=&
              0.0D0
         
 !        if(thisOctal%ciiTRANSITION(subcell, ilevel,jlevel) == 0.d0) then
 !           print *, "Acoeff",Acoeff
 !           print *, "Bcoeff",Bcoeffs
 !           print *, "Ccoeff",Ccoeffs
 !           print *, "FIELD", FIELD
 !        else 
 !           print *, "TRANSITION IS ", thisOctal%ciiTRANSITION(subcell, ilevel,jlevel)
 !        end if
         
         
         
         
      elseif(id == 2) then
         thisOctal%ciTRANSITION(subcell, ilevel,jlevel)=Acoeff(ilevel, jlevel)&
              +Bcoeffs(ilevel, jlevel)*FIELD(ilevel, jlevel)&
              +Ccoeffs(ilevel, jlevel)
         
  !       thisOctal%ciTRANSITION(subcell, jlevel,ilevel)=Acoeff&
  !            +Bcoeffs*(-FIELD)&
  !            +Ccoeffs
         IF(ABS(thisOctal%ciTRANSITION(subcell, ilevel,jlevel)).LT.1.0D-50) thisOctal%ciTRANSITION(subcell,ilevel,jlevel)=&
              0.0D0
         
  !       if(thisOctal%ciiTRANSITION(subcell, ilevel,jlevel) == 0.d0) then
  !          print *, "Acoeff",Acoeff
  !          print *, "Bcoeff",Bcoeffs
  !          print *, "Ccoeff",Ccoeffs
  !          print *, "FIELD", FIELD
  !       end if
         
         
      elseif(id == 3) then
         thisOctal%oiTRANSITION(subcell, ilevel,jlevel)=Acoeff(ilevel, jlevel)&
              +Bcoeffs(ilevel, jlevel)*FIELD(ilevel, jlevel)&
              +Ccoeffs(ilevel, jlevel)
   !      thisOctal%oiTRANSITION(subcell, jlevel,ilevel)=Acoeff&
   !           +Bcoeffs*(-FIELD)&!
   !           +Ccoeffs
         IF(ABS(thisOctal%oiTRANSITION(subcell, ilevel,jlevel)).LT.1.0D-50) thisOctal%oiTRANSITION(subcell,ilevel,jlevel)=&
              0.0D0
         
         !     if(debug .and. ilevel == 3 .and. jlevel == 1) then
         !        print *,  thisOctal%oiTRANSITION(subcell, ilevel,jlevel)
         !        print *, "Acoeff ", Acoeff
         !        print *, "Bcoeffs", Bcoeffs
         !        print *, "Ccoeffs", Ccoeffs
         !        print *, "FIELD", FIELD
         !        print *, "tau", tau
         !        print *, " "
         !     endif
         !
      elseif(id == 4) then
         thisOctal%c12oTRANSITION(subcell, ilevel,jlevel)=Acoeff(ilevel, jlevel)&
              +Bcoeffs(ilevel, jlevel)*FIELD(ilevel, jlevel)&
              +Ccoeffs(ilevel, jlevel)
    !     thisOctal%c12oTRANSITION(subcell, jlevel,ilevel)=Acoeff&
    !          +Bcoeffs*(-FIELD)&
    !          +Ccoeffs
         IF(ABS(thisOctal%c12oTRANSITION(subcell, ilevel,jlevel)).LT.1.0D-50) thisOctal%c12oTRANSITION(subcell,ilevel,jlevel)=&
              0.0D0
      endif
   enddo
enddo
end subroutine setupEscapeParameters2

subroutine abundanceSweepGovernor(grid, reactant, &
     product, alpha, beta, gamma, rate, duplicate, rtmin, rtmax & 
     , n12co, ncii, nci, noi, nelect, nChemIter, partLTE)
  integer :: nChemIter, i
  type(gridtype) :: grid
  real(double) :: rate(:), alpha(:), beta(:), gamma(:)
  real(double) :: rtmin(:), rtmax(:)
  integer :: duplicate(:)
  character(len=10) :: product(:,:), reactant(:,:)
  integer :: n12co, ncii, nci, noi, nelect, ier
  !first wave
  logical :: partLTE

!  integer :: threadCounter

  do i = 1, nChemIter
     if(myrankglobal == 1) print *, "!- chemical iteration ", i
!     if(i == nChemiter) then
!        call abundanceSweepDrone(nChemIter, grid%octreeRoot, reactant, &
!             product, alpha, beta, gamma, rate, duplicate, rtmin, rtmax, nelect, &
!             n12co, ncii, nci, noi, rfile=.false.)
!     else
        call abundanceSweepDrone(nChemIter, grid%octreeRoot, reactant, &
             product, alpha, beta, gamma, rate, duplicate, rtmin, rtmax, nelect, &
             n12co, ncii, nci, noi, rfile=.false.)
!     endif
     if(myrankglobal == 1) print *, "!- abundance sweep done, ray tracing... ", i
     
     call MPI_BARRIER(MPI_COMM_WORLD, ier)
     
     !CALC 33 thisColRhos
     
     call rayTraceMPI(grid, colrhoOnly=.true.)
     if(myrankglobal == 1) print *, "!- ray tracing done ", i
     
     call MPI_BARRIER(MPI_COMM_WORLD, ier)
     !     stop
  end do
  if(partLTE) then
     if(myrankglobal == 1) print *, "!- calling partition LTE "
     call partitionLTE(grid%octreeRoot, n12co, ncii, nci, noi)
     call MPI_BARRIER(MPI_COMM_WORLD, ier)
     !     stop
     if(myrankglobal == 1) print *, "!- partition LTE done "
  endif
end subroutine abundanceSweepGovernor


recursive subroutine abundanceSweepDrone(nChemIter, thisOctal, reactant, &
     product, alpha, beta, gamma, rate, duplicate, rtmin, rtmax, nelect, &
     nc12o, ncii, nci, noi, rfile)
  type(octal), pointer :: thisOctal, child
  integer :: nchemiter, subcell, j
  real(double) :: rate(:), alpha(:), beta(:), gamma(:)
  real(double) :: rtmin(:), rtmax(:), abundanceSol(0:32)
  integer :: duplicate(:)
  character(len=10) :: product(:,:), reactant(:,:)
  integer :: NRGR, NRH2, NRHD, NRCO, NRCI, NRSI, NSPEC, nreac, nrays, nelect
  integer :: nc12o, ncii, nci, noi
  logical :: rfile
!  type(vector) :: pos
!  real(double) :: temp
!  real(double) :: abundIn(33)
  nreac = 329
  nspec = 33

  nrays = nray_func() 
  if(thisOctal%oneD) nrays = 1
!  nchemiter = nchemiter

  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do j = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(j) == subcell) then
              child => thisoctal%child(j)
              call abundanceSweepDrone(nChemIter, child, reactant, &
                   product, alpha, beta, gamma, rate, duplicate, rtmin, rtmax, nelect, &
                   nc12o, ncii, nci, noi, rfile)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle
        if(thisOctal%level_converged(subcell).or. thisOctal%converged(subcell)) cycle
!        open (101, file="ratesfile.dat", status="unknown", position="append")
        if(rfile) open(101, file="ratesfile.dat", form="formatted", status="unknown", position="append")

!     print *, "ALPHA"
!        if(.not.thisOctal%ionfrac(subcell, 2) > 0.99d0) then  !.and. .not. thisOctal%level_converged(subcell)) then
        if(.not.thisOctal%ionfrac(subcell, 2) > 0.99d0) then

           call CALCULATE_REACTION_RATES(thisOctal%tLast(subcell),thisOctal%Dust_T(subcell), &
                thisOctal%radsurface(subcell, :),thisOctal%AV(subcell, :),thisOctal%thisColRho(subcell, :, :), &
                &REACTANT,PRODUCT,ALPHA,BETA,GAMMA,RATE,RTMIN,RTMAX,DUPLICATE,NSPEC,&
                &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI, nreac, nrays, nc12o, nci, debug=.false.)   

!           pos = subcellCentre(thisOctal, subcell)
           
!           call calculate_abundances_onepart(thisOctal%abundance(subcell,:), rate, &
!                thisOctal%rho(subcell)/mhydrogen, thisOctal%tLast(subcell), nspec, nreac, nelect)
           abundancesol(0:32) = thisOctal%abundance(subcell, 1:33)
           
           if(rfile) write(101, '(1p,500e14.5)') thisOctal%AV(subcell, 1), rate

           call calculate_abundances_onepart(abundancesol(:), rate, &
                thisOctal%rho(subcell)/mhydrogen, thisOctal%tLast(subcell), nspec, nreac, nelect)

           thisOctal%abundance(subcell, 1:33) = abundancesol(0:32)
           close(101)
        end if
     end if
  end do
end subroutine abundanceSweepDrone


recursive subroutine partitionLTE(thisOctal, n12co, ncii, nci, noi)
  type(octal), pointer :: thisOctal, child
  integer :: subcell, j
  real(double) :: cii_z_function, ci_z_function, oi_z_function, c12o_z_function
!  integer :: cii_nlev, ci_nlev, oi_nlev, c12o_nlev
!!  real(double) :: cii_energies(1:5), ci_energies(1:5), oi_energies(1:5), c12o_energies(1:41)
!  real(double) :: cii_weights(1:5), ci_weights(1:5), oi_weights(1:5), c12o_weights(1:5)
  integer :: n12co, ncii, nci, noi
!  real(double) :: rate(:), alpha(:), beta(:), gamma(:)
!
!  real(double) :: rtmin(:), rtmax(:)
!  integer :: duplicate(:)
!  character(len=10) :: product(:,:), reactant(:,:)
!  integer :: NRGR, NRH2, NRHD, NRCO, NRCI, NRSI, NSPEC, nreac, nrays


  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do j = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(j) == subcell) then
              child => thisoctal%child(j)
              call partitionLTE(child, n12co, ncii, nci, noi)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle
        if(thisOctal%level_converged(subcell) .or. thisOctal%converged(subcell)) cycle     

!     cii_nlev = 5
!     ci_nlev = 5
!     oi_nlev = 5
!     c12o_nlev = 41

!        print *, "ORION"
     if(.not.thisOctal%ionfrac(subcell, 2) > 0.99d0)then 

           !       print *, "PERSEUS"
           CALL CALCULATE_PARTITION_FUNCTION(CII_Z_FUNCTION,CII_NLEV,CII_ENERGIES,CII_WEIGHTS,thisOctal%tLast(subcell))
           CALL CALCULATE_PARTITION_FUNCTION(CI_Z_FUNCTION,CI_NLEV,CI_ENERGIES,CI_WEIGHTS,thisOctal%tLast(subcell))
           CALL CALCULATE_PARTITION_FUNCTION(OI_Z_FUNCTION,OI_NLEV,OI_ENERGIES,OI_WEIGHTS,thisOctal%tLast(subcell))
           CALL CALCULATE_PARTITION_FUNCTION(C12O_Z_FUNCTION,C12O_NLEV,C12O_ENERGIES,C12O_WEIGHTS,thisOctal%tLast(subcell))
           
           !
           ! Calculate the LTE level populations
           !      print *, "cii abund ", thisOctal%abundance(subcell, Ncii)
!           print *, "pre ", thisOctal%CII_POP(subcell, :)
           CALL CALCULATE_LTE_POPULATIONS(CII_NLEV,thisOctal%CII_POP(subcell, :),CII_ENERGIES,&
             &CII_WEIGHTS,CII_Z_FUNCTION,thisOctal%abundance(subcell, Ncii)*thisOctal%rho(subcell)/mhydrogen,&
             thisOctal%tLast(subcell))
!           print *, thisOctal%av(subcell, 1), thisOctal%cii_pop(subcell, :)
           CALL CALCULATE_LTE_POPULATIONS(CI_NLEV, thisOctal%CI_POP(subcell, :), CI_ENERGIES, &
                &CI_WEIGHTS, CI_Z_FUNCTION, thisOctal%abundance(subcell, Nci)*thisOctal%rho(subcell)/mhydrogen,&
                thisOctal%tLast(subcell))
           CALL CALCULATE_LTE_POPULATIONS(OI_NLEV, thisOctal%OI_POP(subcell, :), OI_ENERGIES, &
                &OI_WEIGHTS, OI_Z_FUNCTION, thisOctal%abundance(subcell, Noi)*thisOctal%rho(subcell)/mhydrogen,&
                thisOctal%tLast(subcell))
           CALL CALCULATE_LTE_POPULATIONS(C12O_NLEV,thisOctal%C12O_POP(subcell, :),C12O_ENERGIES,&
                &C12O_WEIGHTS,C12O_Z_FUNCTION,thisOctal%abundance(subcell, N12co)*thisOctal%rho(subcell)/mhydrogen,&
                thisOctal%tLast(subcell))

 !       endif
   !     print *, "cii abund ", thisOctal%abundance(subcell, Ncii)
!     else
!        thisOctal%CII_POP(subcell, :) = 0.d0
!        thisOctal%CI_POP(subcell, :) = 0.d0
!        thisOctal%OI_POP(subcell, :) = 0.d0
!        thisOctal%C12O_POP(subcell, :) = 0.d0
!        end if
        end if
     end if
  end do
end subroutine partitionLTE






subroutine rayTraceMPI(grid, colrhoonly)

  implicit none

  type(gridtype) :: grid
!  type(sourcetype) :: thisSource
  integer :: i, ier
  logical :: colrhoonly

  if(myrankglobal /= 0) then
     do i = 1, nhydrothreadsglobal
        if(myrankglobal == i) then
!           print *, "rank ", myrankglobal, "casting rays"
           call castAllRaysOverGridMPI(grid%octreeRoot, grid, colrhoonly)
!           print *, "rank ", myrankglobal, "shutting down"
           call shutdownserversXRAY_PDR()
        else
!           print *, "rank ", myrankglobal, "serving"
           call raytracingserverPDR(grid)
        end if
!        print *, "rank ", myrankglobal, "at barrier"
        call MPI_BARRIER(amrCommunicator, ier)
     end do            
  end if
!  print *, "rank ", myrankglobal, "at barrier"
  call MPI_BARRIER(MPI_COMM_WORLD, ier)
!  call calcIonParamTemperature(grid%octreeRoot)


end subroutine rayTraceMPI


!This is the main routine that loops over the grid and does the ray casting
recursive subroutine castAllRaysOverGridMPI(thisOctal, grid, colrhoonly)
!  use mpi
  use inputs_mod, only : hlevel, maxdepthamr
!  use healpix_module, only : vectors, nrays
  implicit none
  
  type(gridtype) :: grid
  integer :: subcell
  type(octal), pointer :: thisoctal
  type(octal), pointer :: child
  integer :: j
  type(vector) :: testPosition, rVec, startPosition, uhat, thisUVvector
  real(double) :: tVal, rho, uvx, uvy, uvz, HplusFrac
  integer ::  i, ncontributed, k, nrays!, dummy
  real(double) :: abundanceArray(1:33)
  logical :: colrhoonly, didSurface

  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do j = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(j) == subcell) then
              child => thisoctal%child(j)
              call castAllRaysOverGridMPI(child, grid, colrhoonly)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle

        nrays = nray_func()
        if(thisOctal%oneD) nrays = 1
!        call findSubcellLocal(testPosition, sOctal, ssubcell)
        startPosition = subcellCentre(thisOctal, subcell)
        testPosition = startPosition

        thisOctal%thisColRho(subcell, :, :)  = 0.d0
        if(colrhoonly) then
           if(.not. thisOctal%ionfrac(subcell, 2) > 0.99d0) then
              ncontributed = 0
              
              do i = 1, nrays
                 
                 rVec = VECTOR(vectors(1, i-1), vectors(2, i-1), vectors(3, i-1))
                 if(thisOctal%oned) rvec = VECTOR(-1.d0, 0.d0, 0.d0)
                 uhat = rVec
                 call normalize(uhat)
                 thisOctal%columnRho(subcell) = 0.d0                 
                 testPosition = startPosition
!                 thisOctal%thisColRho
                 do while (inOctal(grid%octreeRoot, testPosition))
                    
                    tval = 0.d0
                    call getRayTracingValuesPDR(grid, testposition, uHat, rho, uvx, uvy, uvz, Hplusfrac, tval, &
                         abundancearray) 
!                    do dummy = 1, 33
!                       print *,"abundancearray(",dummy,")", abundancearray(dummy)
!                       !                    print *, "BETA ", sum(abundancearray(:))
!                    enddo

                    thisOctal%columnRho(subcell) = thisOctal%columnRho(subcell) + (rho*tval*1.d10/mhydrogen)          
                       
                    do k = 1, 33
                       thisOctal%thisColRho(subcell, i, k)  = thisOctal%thisColRho(subcell, i, k) + &
                            abundancearray(k)*tval*1.d10*rho/mhydrogen                 
                    enddo
                    
                    testPosition = testPosition + ((tVal+1.d-10*grid%halfsmallestsubcell)*uhat)                           
!                    if(Hplusfrac > 0.99d0 .or. .not. inoctal(grid%octreeroot, testposition)) then       
!                       testPosition = testPosition - (((tVal+1.d-10*grid%halfsmallestsubcell))*uhat)       
!                       call getRayTracingValuesPDR(grid, testPosition, uHat, rho, uvx, uvy, uvz, Hplusfrac, tval, &
!                            abundancearray) 
!                       thisOctal%columnRho(subcell) = thisOctal%columnRho(subcell) + (rho*tval*1.d10/mhydrogen)          
!                       do k = 1, 33
!!                          call getRayTracingValuesPDR(grid, testposition, uHat, rho, uvx, uvy, uvz, Hplusfrac, tval, &
!!                               abundancearray) 
!                          thisOctal%thisColRho(subcell, i, k)  = thisOctal%thisColRho(subcell, i, k) + &
!                               abundancearray(k)*tval*1.d10*rho/mhydrogen                 
!                       enddo
!                       
!!                       exit
!                       goto 555
!                    end if
!                                 

 !                   if(thisOctal%oneD .and. .not. inOctal(grid%octreeroot, testposition)) then
 !                      testPosition = testPosition - (((tVal+1.d-10*grid%halfsmallestsubcell))*uhat)       
 !                      call getRayTracingValuesPDR(grid, testposition, uHat, rho, uvx, uvy, uvz, Hplusfrac, tval, &
 !                           abundancearray) 
 !                      thisOctal%columnRho(subcell) = thisOctal%columnRho(subcell) + (rho*tval*1.d10/mhydrogen)          
 !                      do k = 1, 33
 !                         thisOctal%thisColRho(subcell, i, k)  = thisOctal%thisColRho(subcell, i, k) + &
 !                              abundancearray(k)*tval*1.d10*rho/mhydrogen                 
 !                      enddo
 !                      goto 555
 !!                      exit
 !                   endif


                    
                 end do
!555              continue
              end do
              
           else
              thisOctal%thisColRho(subcell, :, :) = 0.d0
!              thisOctal%columnRho(subcell) = 0.d0
!              thisOctal%UV(subcell) = 0.d0
!              thisOctal%AV(subcell,:) = 0.d0
!              thisOctal%radsurface(subcell,:) = 0.d0
           end if
           
        else
           thisOctal%uv(subcell) = 0.d0
           thisOctal%av(subcell, :) = 0.d0           
           thisOctal%radsurface(subcell, :) = 0.d0           
           thisOctal%thisColRho(subcell, :, :) = 0.d0
           thisOctal%columnrho(subcell)  = 0.d0
           
           if(.not. thisOctal%ionfrac(subcell, 2) > 0.99d0) then
              ncontributed = 0
              
              do i = 1, nrays
!                 print *, "doing ray ", i
                 rVec = VECTOR(vectors(1, i-1), vectors(2, i-1), vectors(3, i-1))
                 if(thisOctal%oneD) Rvec = VECTOR(-1.d0, 0.d0, 0.d0)
                 uhat = rVec
                 call normalize(uhat)
                 
                 thisOctal%columnRho(subcell) = 0.d0
                 testPosition = startPosition
 !                print *, "direction is ", uhat
 !                print *, "position is ", testPosition
                 didSurface = .false.
                 do while (inOctal(grid%octreeRoot, testPosition))

                    tval = 0.d0
                    call getRayTracingValuesPDR(grid, testPosition, uHat, rho, uvx, uvy, uvz, Hplusfrac, tval, &
                         abundancearray) 

!                    print *, "ALPHA ", sum(abundancearray(:))

                    thisOctal%columnRho(subcell) = thisOctal%columnRho(subcell) + (rho*tval*1.d10/mhydrogen)          
                    !!                    print *, "column rho is", thisOCtal%columnrho(subcell)
                    do k = 1, 33
                       thisOctal%thisColRho(subcell, i, k)  = thisOctal%thisColRho(subcell, i, k) + &
                            abundancearray(k)*tval*1.d10*rho/mhydrogen  
                    enddo




                    testPosition = testPosition + ((tVal+1.d-10*grid%halfsmallestsubcell)*uhat)         

                    if(Hplusfrac > 0.99d0 .or. .not. inOctal(grid%octreeRoot, testposition)) then
!                       testPosition = testPosition - (((tVal+1.d-10*grid%halfsmallestsubcell))*uhat)       
!                       thisOctal%columnRho(subcell) = thisOctal%columnRho(subcell) + (rho*tval*1.d10/mhydrogen)          
!                       call getRayTracingValuesPDR(grid, testPosition, uHat, rho, uvx, uvy, uvz, Hplusfrac, tval, &
!                            abundancearray) 
!                       !!                    print *, "column rho is", thisOCtal%columnrho(subcell)
!                       do k = 1, 33
!                          thisOctal%thisColRho(subcell, i, k)  = thisOctal%thisColRho(subcell, i, k) + &
!                               abundancearray(k)*tval*1.d10*rho/mhydrogen  
!                       enddo!!

                                                                   
                       if(thisOctal%oned) then
                          !THAW - not ideal that this is hardwired atm for 1D tests
                          thisUVvector = VECTOR(10.d0*draine/1.d10, 0.d0, 0.d0)
                       else                         
                          thisUVvector = VECTOR(uvx, uvy, uvz)
                       endif
                       thisUVvector = thisUVvector*1.d10/draine!

                       thisOctal%radsurface(subcell, i) = - dotprod(uHat,thisUVvector)     
                       !                      print *, "radsuface ", thisOctal%radsurface(subcell, i)
                       didSurface = .true.
                       if(thisOctal%radsurface(subcell, i) < 0.d0 ) thisOctal%radsurface(subcell, i) = 0.d0
!                       goto 666
!                       exit
                    end if

                   
                 end do
!666 continue
                 
                 !                 print *, "out of loop"
                 if(.not. didSurface) then
                    if(thisOctal%oned) then
                       !THAW - not ideal that this is hardwired atm for 1D tests
                       thisUVvector = VECTOR(10.d0*draine/1.d10, 0.d0, 0.d0)
                    else                         
                       thisUVvector = VECTOR(uvx, uvy, uvz)
                    endif
                    thisUVvector = thisUVvector*1.d10/draine
                    
                    thisOctal%radsurface(subcell, i) = - dotprod(uHat,thisUVvector)     
                    !                    print *, "radsuface ", thisOctal%radsurface(subcell, i)
                    didSurface = .true.
                    if(thisOctal%radsurface(subcell, i) < 0.d0 ) thisOctal%radsurface(subcell, i) = 0.d0
                 endif
                 
                 thisOctal%AV(subcell, i) = thisOctal%columnRho(subcell)*AV_fac
                 
                 if(thisOctal%oneD) thisOctal%AV(subcell, 1) = thisOctal%columnRho(subcell)*AV_fac
                 !                 print *, "colrho is ", thisOctal%columnRho(subcell)
                 !                 print *, "AV is ", thisOctal%AV(subcell, i)
                 
                 !                 print *, "think ive traced a distance ", thisOctal%columnRho(subcell)/(1.d3*pctocm), "pc"
                 
                 !                 print *, "rad surface is", thisOctal%radsurface(subcell, i)
                 
                 if(thisOctal%radsurface(subcell, i) > 0.d0) then
                    thisOctal%UV(subcell) = thisOctal%UV(subcell) + &
                         thisOctal%radSurface(subcell, i) * exp(-(thisOctal%AV(subcell, i)*UV_fac))
                    ncontributed = ncontributed + 1
                    !                   if(thisOctal%oned) thisOctal%UV(subcell) = thisOctal%UV(subcell) + &
                    !                     thisOctal%radSurface(subcell, i) * exp(-(thisOctal%AV(subcell, i)*UV_fac))
                 end if
                 !                 print *, "UV  is ", thisOctal%UV(subcell)                 
              end do
              if(thisOctal%threed) then
                 if(ncontributed /= 0) then
                    thisOctal%UV(subcell) = thisOctal%UV(subcell) / dble(ncontributed)
                 else
                    thisOctal%UV(subcell) = 0.d0
                 end if
                 
                 !                 thisOctal%columnRho(subcell) = thisOctal%columnRho(subcell)/dble(ncontributed)
              endif
              !              if(thisOctal%UV(subcell) < 1.d-30) thisOctal%UV(subcell) = 0.d0
           else
              print *, "stopped on ion frac"
              thisOctal%thisColRho(subcell, :, :) = 0.d0
              thisOctal%columnRho(subcell) = 0.d0
              thisOctal%UV(subcell) = 0.d0
              thisOctal%AV(subcell,:) = 0.d0
              thisOctal%radsurface(subcell,:) = 0.d0
           end if
           !  endif
        end if
     end if
  end do
  

end subroutine castAllRaysOverGridMPI
#endif


!This is the main routine that loops over the grid and does the ray casting
recursive subroutine castAllRaysOverGrid(thisOctal, grid)
  use inputs_mod, only : hlevel, maxdepthamr
!  use healpix_module, only : vectors, nrays
  implicit none
  
  type(gridtype) :: grid
  integer :: subcell, ssubcell
  type(octal), pointer :: thisoctal, soctal
  type(octal), pointer :: child
  integer :: j
  type(vector) :: testPosition, rVec, startPosition, uhat, thisUVvector
  real(double) :: tVal
  integer ::  i, ncontributed, k, nrays



  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do j = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(j) == subcell) then
              child => thisoctal%child(j)
              call castAllRaysOverGrid(child, grid)
              exit
           end if
        end do
     else
!        if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle
!        print *, "thisOctal%ndepth", thisOctal%ndepth
 !       print *, "casting rays at", subcell!, subcellCentre(thisOctal, subcell)

!        call doSingleCellRays(thisOctal, subcell, grid)

        nrays = nray_func()
        if(thisOctal%oned) nrays = 1
        sOctal=> thisOctal
        
        startPosition = subcellCentre(thisOctal, subcell)
        testPosition = startPosition
        call findSubcellLocal(testPosition, sOctal, ssubcell)
        sOctal%columnRho(ssubcell) = 0.d0
!        thisOctal%columnRho(subcell) = 1.d-30 

        thisOctal%uv(subcell) = 0.d0
        thisOctal%av(subcell, :) = 0.d0           
        thisOctal%radsurface(subcell, :) = 0.d0           
        thisOctal%thisColRho(subcell, :, :) = 0.d0
        if(.not. thisOctal%ionfrac(subcell, 2) > 0.99d0) then
           ncontributed = 0
!        do i = 0, nrays-1
        do i = 1, nrays

           rVec = VECTOR(vectors(1, i-1), vectors(2, i-1), vectors(3, i-1))
           if(thisOctal%oned) rvec = VECTOR(-1.d0, 0.d0, 0.d0)
           uhat = rVec
           call normalize(uhat)

           testPosition = startPosition
           thisOctal%columnrho(subcell)  = 0.d0
           do while (inOctal(grid%octreeRoot, testPosition))
              
              call findSubcellLocal(testPosition, sOctal, ssubcell)

              tval = 0.d0

             call distanceToCellBoundary(grid, testPosition, uHat, tVal, soctal, ssubcell)
 
              if(sOctal%ionFrac(ssubcell,2) > 0.99d0) then
                 thisUVvector = sOctal%uvVector(ssubcell)*1.d10/draine
!                 call normalize(thisUVvector)
!                 print *, "thisUVvector ", thisUVvector
                 thisOctal%radsurface(subcell, i) = - dotprod(uHat,thisUVvector)     
!                 thisOctal%radsurface(subcell, i) =  dotprod(uHat,thisUVvector)     
!                 print *, "thisOctal%radsurface(subcell, i)", thisOctal%radsurface(subcell, i)
!                 print *, "radsurface", thisOctal%radsurface(subcell, i)
                 if(thisOctal%radsurface(subcell, i) < 0.d0 ) thisOctal%radsurface(subcell, i) = 0.d0

                 exit
              end if

              
              thisOctal%columnRho(subcell) = thisOctal%columnRho(subcell) + (sOctal%rho(ssubcell)*tval)
              !loop over all species
              do k = 1, 33
                thisOctal%thisColRho(subcell, k, i)  = thisOctal%thisColRho(subcell, k, i) + &
                     sOctal%abundance(ssubcell, k)*tval*1.d10*soctal%rho(subcell)/mhydrogen
                 
              enddo
!              print *, "RAY ", I, "HAS"
!              print *, thisOctal%thiscolRho(subcell, :, i)
              testPosition = testPosition + ((tVal+1.d-10*grid%halfsmallestsubcell)*uhat)
              
           end do


!          thisOctal%AV(subcell, i) = thisOctal%columnRho(subcell)*1.d10*AV_fac/mhydrogen
          thisOctal%AV(subcell, i) = thisOctal%columnRho(subcell)*1.d10*AV_fac/mhydrogen


          if(thisOctal%radsurface(subcell, i) > 0.d0) then
             thisOctal%UV(subcell) = thisOctal%UV(subcell) + &
                  thisOctal%radSurface(subcell, i) * exp(-(thisOctal%AV(subcell, i)*UV_fac))
             ncontributed = ncontributed + 1
!             print *, "UV ", thisOctal%UV(subcell)
!             print *, "exp(-(thisOctal%AV(subcell, i)*UV_fac", exp(-(thisOctal%AV(subcell, i)*UV_fac))

!             print *,thisOctal%radSurface(subcell, i) * exp(-(thisOctal%AV(subcell, i)*UV_fac))
          end if
          
       end do
!       print *, "DONE ONE "
!       print *, 'ncontributed ', ncontributed
       if(ncontributed /= 0) then
!          thisOctal%UV(subcell) = thisOctal%UV(subcell) / dble(nrays)
          thisOctal%UV(subcell) = thisOctal%UV(subcell) / dble(ncontributed)
       else
          thisOctal%UV(subcell) = 0.d0
       end if
!       thisOctal%UV(subcell) = thisOctal%UV(subcell) / dble(nrays)
       thisOctal%columnRho(subcell) = thisOctal%columnRho(subcell)/dble(nrays)
 !      print *, "AV ", thisOctal%AV(subcell,:)
  !     print *, "UV ", thisOctal%UV(subcell)
       
!       stop

       if(thisOctal%UV(subcell) < 1.d-30) thisOctal%UV(subcell) = 0.d0
    else
        thisOctal%columnRho(subcell) = 0.d0
        thisOctal%UV(subcell) = 0.d0
        thisOctal%AV(subcell,:) = 0.d0
        thisOctal%radsurface(subcell,:) = 0.d0
     end if
        
  end if
end do
  

end subroutine castAllRaysOverGrid



!A quick test of ray tracing
subroutine fireTestRays(grid)
  use inputs_mod, only : amrgridcentrex, amrgridcentrey, amrgridcentrez
  use inputs_mod, only : hlevel, maxdepthamr
!  use healpix_module, only : vectors, nrays

  implicit none

  type(gridtype) :: grid
  type(vector) :: testPosition, rVec, startPosition, uhat
  integer :: subcell
  type(octal), pointer :: thisOctal
  real(double) :: x, y, z, tVal, totalLength, totalRho
  integer :: nside, i, ier, nrays


  nside = 2**hlevel
  nrays = 12*nside**2

  thisOctal=>grid%octreeRoot
  
  !choose a test location
  x = amrgridcentrex + grid%octreeroot%subcellsize*0.34d0
  y = amrgridcentrey + grid%octreeroot%subcellsize*0.34d0
  z = amrgridcentrez + grid%octreeroot%subcellsize*0.34d0
  startPosition = VECTOR(x, y, z)

  call findSubcellLocal(startPosition, thisOctal, subcell)
  
  print *, "testing rays from cell at ", startPosition
  print *, "testing with ", nrays, "rays"
  print *, "cell size is ", grid%octreeroot%subcellsize/2.**maxdepthamr/1.d8
  print *, "grid size is" , grid%octreeroot%subcellsize/1.d8

  if (thisoctal%haschild(subcell)) then
     call torus_abort("Octal has child!")
  end if

  open (1, file="testRays.dat", status="unknown", iostat=ier)
  if(ier /= 0) then
     call torus_abort("Problem opening testRays.dat")
  end if

  do i = 0, nrays-1

     !transfer healpix vectors to the test position
!     rVec = startPosition + VECTOR(vectors(1, i), vectors(2, i), vectors(3, i))
     rVec = VECTOR(vectors(1, i), vectors(2, i), vectors(3, i))
     uhat = rVec
     call normalize(uhat)
     totalLength = 0.d0
     totalRho = 0.d0
     testPosition = startPosition
     print *, "doing ray ", i!, " in direction ", uHat 
     do while (inOctal(grid%octreeRoot, testPosition))
!        print *, "walking in direction ", uHat
        call findSubcellLocal(testPosition, thisOctal, subcell)
        call distanceToCellBoundary(grid, testPosition, uHat, tVal, thisoctal)
!        print *, " ini ", testposition
        testPosition = testPosition + ((tVal+1.d-10*grid%halfsmallestsubcell)*uhat)
 !       print *, "tval", tval
  !      print *, "uhat", uhat
   !     print *, "post ", testposition
    !    stop

        write(1,'(1p,i3,3e14.2)') i, testPosition%x*1.d10/pctocm, testposition%y*1.d10/pctocm, &
             testposition%z*1.d10/pctocm
        totalLength = totalLength + (tVal)
!       print *, totallength/1.d8
        totalRho = totalRho + thisOctal%rho(subcell)*tval

     end do
     print *, "ray ", i, "travelled a distance", totalLength*1.d10/pctocm
     print *, "ray ", i, "had total density ", totalRho
     print *, "ray ", i, "had total H number density ", totalRho/mHydrogen

!     print *, "ray ", i, "had start ", totalRho
  end do

  close(1)

end subroutine fireTestRays

!subroutine


!=======================================================================
!
!  Calculate the dust temperature for each particle using the treatment
!  of Hollenbach, Takahashi & Tielens (1991, ApJ, 377, 192, eqns 5 & 6)
!  for the heating due to the incident FUV photons and the treatment of
!  Meijerink & Spaans (2005, A&A, 436, 397, eqn B.6) for heating due to
!  the incident flux of X-ray photons.
!
!  Among other things, the dust temperature can influence:
!
!     1) Cooling budget by emitting FIR photons that
!        interact with the line radiative transfer;
!     2) Gas-grain collisional heating or cooling rate;
!     3) H2 formation by changing the sticking probability;
!     4) Evaporation and condensation of molecules on grains.
!
!  The formula derived by Hollenbach, Takahashi & Tielens (1991) has
!  been modified to include the attenuation of the IR radiation. The
!  incident FUV radiation is absorbed and re-emitted in the infrared
!  by dust at the surface of the cloud (up to Av ~ 1mag). In the HTT
!  derivation, this IR radiation then serves as a second heat source
!  for dust deeper into the cloud. However, in their treatment, this
!  second re-radiated component is not attenuated with distance into
!  the cloud so it is *undiluted* with depth, leading to higher dust
!  temperatures deep within the cloud which in turn heat the gas via
!  collisions to unrealistically high temperatures. Models with high
!  gas densities and high incident FUV fluxes (e.g. n_H = 10^5 cm-3,
!  X_0 = 10^8 Draine) can produce T_gas ~ 100 K at Av ~ 50 mag!
!
!  Attenuation of the FIR radiation has therefore been introduced by
!  using an approximation for the infrared-only dust temperature from
!  Rowan-Robinson (1980, eqn 30b):
!
!  T_dust = T_0*(r/r_0)^(-0.4)
!
!  where r_0 is the cloud depth at which T_dust = T_0, corresponding
!  to an A_V of ~ 1 mag, the assumed size of the outer region of the
!  cloud that processes the incident FUV radiation and then re-emits
!  it in the FIR (see the original HTT 1991 paper for details). This
!  should prevent the dust temperature from dropping off too rapidly
!  with distance and maintain a larger warm dust region (~50-100 K).
!
!-----------------------------------------------------------------------
recursive SUBROUTINE CALCULATE_DUST_TEMPERATURES(thisOctal)

!   USE HEALPIX_TYPES
!   USE MAINCODE_MODULE

   IMPLICIT NONE

   INTEGER :: J, i
   type(octal), pointer :: thisOctal, child
   integer :: subcell, nrays
   REAL(double) :: NU_0,R_0,T_0,TAU_100
   REAL(double) :: T_CMB

!  Parameters used in the HHT equations (see their paper for details)
   NU_0=2.65D15
   TAU_100=1.0D-3
   R_0=1.0D0/AV_FAC
   T_CMB=2.73D0

   nrays = nray_func()
   if(thisOctal%oned) nrays = 1
  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do i = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(i) == subcell) then
              child => thisoctal%child(i)
              call calculate_dust_temperatures(child)
              exit
           end if
        end do
     else
!        P=IDlist_pdr(pp)
        !     Calculate the contribution to the dust temperature from the local FUV flux and the CMB background
!        PDR(P)%DUST_T=8.9D-11*NU_0*(1.71D0*PDR(P)%UVfield)+T_CMB**5
        if(.not. thisOctal%ionfrac(subcell, 2) > 0.99d0) then
           thisOctal%DUST_T(subcell)=8.9D-11*NU_0*(1.71D0*thisOctal%UV(subcell))+T_CMB**5
           
           DO J=1,NRAYS ! Loop over rays
              
              !        The minimum dust temperature is related to the incident FUV flux along each ray
              !        Convert the incident FUV flux from Draine to Habing units by multiplying by 1.7
              T_0=12.2*(1.71D0*thisOctal%RADSURFACE(subcell, j))**0.2
              
!!$!        Attenuate the FIR radiation produced in the surface layer
!!$         IF(PARTICLE(P)%TOTAL_COLUMN(J).GT.R_0) THEN
!!$            T_0=T_0*(PDR(P)%TOTAL_COLUMN(J)/R_0)**(-0.4)
!!$         END IF
              
              !        Add the contribution to the dust temperature from the FUV flux incident along this ray
              IF(T_0.GT.0) thisOctal%DUST_T(subcell)=thisOctal%DUST_T(subcell) &
                   & + (0.42-LOG(3.45D-2*TAU_100*T_0))*(3.45D-2*TAU_100*T_0)*T_0**5
              
           END DO ! End of loop over rays
           
           !     Convert from total dust emission intensity to dust temperature
           thisOctal%DUST_T(subcell)=thisOctal%DUST_T(subcell)**0.2
           
           !#ifdef XRAYS2
           !        STOP '[dust_t.F90] not coded for XRAYS=ON'
           !!       !     Calculate the contribution to the dust temperature from the local X-ray flux (assuming a fixed grain abundance of 1.6E-8)
           !       PDR(P)%DUST_T=PDR(P)%DUST_T+1.5D4*(PDR(P)%XRAY_ENERGY_DEPOSITION_RATE/1.6D-8)**0.2
           !#endif
           !     Impose a lower limit on the dust temperature, since values below 10 K can dramatically
           !     limit the rate of H2 formation on grains (the molecule cannot desorb from the surface)
           IF(thisOctal%DUST_T(subcell).LT.10) THEN
              thisOctal%DUST_T(subcell)=10.0D0
           END IF
           
           !     Check that the dust temperature is physical
           IF(thisOctal%DUST_T(subcell).GT.1000) THEN
              WRITE(6,*) 'ERROR! Calculated dust temperature exceeds 1000 K'
              STOP
           END IF
        else
           thisOctal%Dust_T(subcell) = 0.d0
        end if

     end if
  END DO ! End of loop over particles
     
END SUBROUTINE CALCULATE_DUST_TEMPERATURES
!=======================================================================

recursive subroutine iniTempGuess(thisOctal)
  integer :: subcell, i
  type(octal), pointer :: thisOctal, child
  real(double) :: tguess
  real(double), parameter :: tmin = 10.d0
  real(double), parameter :: tmax = 1.d4

  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do i = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(i) == subcell) then
              child => thisoctal%child(i)
              call iniTempGuess(child)
              exit
           end if
        end do
     else
        if(.not. thisOCtal%ionfrac(subcell, 2) > 0.99d0) then
           Tguess = 10.0D0*(1.0D0+(1.0D2*thisOctal%UV(subcell))**(1.0D0/3.0D0))
           thisOctal%TLast(subcell) = Tguess
           thisOctal%temperature(subcell) = real(Tguess)
           
           thisOctal%Tmin(subcell) = Tguess/2.0D0
           thisOctal%tMax(subcell) = Tguess*1.5D0

           
           
           IF (thisOctal%tmin(subcell).LT.Tmin)  thisOctal%tmin(subcell)  = Tmin
           IF (thisOctal%tmax(subcell).GT.Tmax) thisOctal%tmax(subcell) = Tmax
           
           thisOctal%Tminarray(subcell) =  thisOctal%tmin(subcell)/3.0D0 ! Bound minimum
           thisOCtal%tmaxarray(subcell) = thisOctal%tmax(subcell)*2.0D0 ! Bound maximum
           IF (thisOctal%tminarray(subcell).LT.Tmin) thisOctal%tminarray(subcell) = Tmin
           IF (thisOctal%tmaxarray(subcell).GT.Tmax) thisOctal%tmaxarray(subcell) = Tmax
        end if
     end if
  ENDDO
end subroutine iniTempGuess

#ifdef MPI
recursive subroutine iniTempGuessMPI(thisOctal)
  integer :: subcell, i
  type(octal), pointer :: thisOctal, child
  real(double) :: tguess
  real(double), parameter :: tmin = 10.d0
!  real(double), parameter :: tmin = 10.d0
  real(double), parameter :: tmax = 1.d4

  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do i = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(i) == subcell) then
              child => thisoctal%child(i)
              call iniTempGuessMPI(child)
              exit
           end if
        end do
     else

        if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle

        if(.not. thisOCtal%ionfrac(subcell, 2) > 0.99d0) then

           Tguess = 10.0D0*(1.0D0+(1.0D2*thisOctal%UV(subcell))**(1.0D0/3.0D0))


           !THAW ! ! ! ! ! ! ! ! ! ! ! !  ! ! !

!           Tguess = 70.d0

           thisOctal%TLast(subcell) = Tguess
           thisOctal%TPrev(subcell) = Tguess
           
           thisOctal%Tlow(subcell) = Tguess/2.0D0
!           thisOctal%thigh(subcell) = Tguess*1.5D0
           thisOctal%thigh(subcell) = Tguess*2.D0
           
           IF (thisOctal%tlow(subcell).LT.Tmin)  thisOctal%tlow(subcell)  = Tmin
           IF (thisOctal%thigh(subcell).GT.Tmax) thisOctal%thigh(subcell) = Tmax
           
!           thisOctal%Tminarray(subcell) =  thisOctal%tmin(subcell)/3.0D0 ! Bound minimum
!           thisOCtal%tmaxarray(subcell) = thisOctal%tmax(subcell)*2.0D0 ! Bound maximum
!           IF (thisOctal%tminarray(subcell).LT.Tmin) thisOctal%tminarray(subcell) = Tmin
!           IF (thisOctal%tmaxarray(subcell).GT.Tmax) thisOctal%tmaxarray(subcell) = Tmax
           


        end if
     end if
  ENDDO
end subroutine iniTempGuessMPI

recursive SUBROUTINE CALCULATE_DUST_TEMPERATURESMPI(thisOctal)

!   USE HEALPIX_TYPES
!   USE MAINCODE_MODULE

   IMPLICIT NONE

   INTEGER :: J, i
   type(octal), pointer :: thisOctal, child
   integer :: subcell, nrays
   REAL(double) :: NU_0,R_0,T_0,TAU_100
   REAL(double) :: T_CMB

!  Parameters used in the HHT equations (see their paper for details)
   NU_0=2.65D15
   TAU_100=1.0D-3
   R_0=1.0D0/AV_FAC
   T_CMB=2.73D0

   nrays = nray_func()
   if(thisOctal%oned) nrays=1
  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do i = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(i) == subcell) then
              child => thisoctal%child(i)
              call calculate_dust_temperaturesMPI(child)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle
!        P=IDlist_pdr(pp)
        !     Calculate the contribution to the dust temperature from the local FUV flux and the CMB background
!        PDR(P)%DUST_T=8.9D-11*NU_0*(1.71D0*PDR(P)%UVfield)+T_CMB**5
        if(.not. thisOctal%ionfrac(subcell, 2) > 0.99d0) then
           thisOctal%DUST_T(subcell)=8.9D-11*NU_0*(1.71D0*thisOctal%UV(subcell))+T_CMB**5
           
           DO J=1,NRAYS ! Loop over rays
              
              !        The minimum dust temperature is related to the incident FUV flux along each ray
              !        Convert the incident FUV flux from Draine to Habing units by multiplying by 1.7
              T_0=12.2*(1.71D0*thisOctal%RADSURFACE(subcell, j))**0.2
              
!!$!        Attenuate the FIR radiation produced in the surface layer
!!$         IF(PARTICLE(P)%TOTAL_COLUMN(J).GT.R_0) THEN
!!$            T_0=T_0*(PDR(P)%TOTAL_COLUMN(J)/R_0)**(-0.4)
!!$         END IF
              
              !        Add the contribution to the dust temperature from the FUV flux incident along this ray
              IF(T_0.GT.0) thisOctal%DUST_T(subcell)=thisOctal%DUST_T(subcell) &
                   & + (0.42-LOG(3.45D-2*TAU_100*T_0))*(3.45D-2*TAU_100*T_0)*T_0**5
              
           END DO ! End of loop over rays
           
           !     Convert from total dust emission intensity to dust temperature
           thisOctal%DUST_T(subcell)=thisOctal%DUST_T(subcell)**0.2
           
           !#ifdef XRAYS2
           !        STOP '[dust_t.F90] not coded for XRAYS=ON'
           !!       !     Calculate the contribution to the dust temperature from the local X-ray flux (assuming a fixed grain abundance of 1.6E-8)
           !       PDR(P)%DUST_T=PDR(P)%DUST_T+1.5D4*(PDR(P)%XRAY_ENERGY_DEPOSITION_RATE/1.6D-8)**0.2
           !#endif
           !     Impose a lower limit on the dust temperature, since values below 10 K can dramatically
           !     limit the rate of H2 formation on grains (the molecule cannot desorb from the surface)
           IF(thisOctal%DUST_T(subcell).LT.10) THEN
              thisOctal%DUST_T(subcell)=10.0D0
           END IF
           
           !     Check that the dust temperature is physical
           IF(thisOctal%DUST_T(subcell).GT.1000) THEN
              WRITE(6,*) 'ERROR! Calculated dust temperature exceeds 1000 K'
              STOP
           END IF
        else
           thisOctal%Dust_T(subcell) = 0.d0
        end if

     end if
  END DO ! End of loop over particles
     
END SUBROUTINE CALCULATE_DUST_TEMPERATURESMPI
!=======================================================================

#endif



#endif
end module pdr_mod

