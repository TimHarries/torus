#ifdef CHEMISTRY
module chemistry_mod

  use krome_user
  use kind_mod
  use octal_mod
  use messages_mod
  implicit none

!  real(double) :: phl(krome_nPhotoBins)
!  real(double) :: phr(krome_nPhotoBins)
!  real(double) :: phm(krome_nPhotoBins)
!  real(double) :: phJ(krome_nPhotoBins)

contains



recursive subroutine  pathlengthToIntensity(thisOctal, epsOverDt)
  use krome_main
  use krome_user
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  real(double) :: epsOverDt, v, dfreq
  integer :: i, subcell

  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call pathLengthToIntensity(child, epsOverDt)
              exit
           end if
        end do
     else

!        do i = 1, krome_nPhotoBins
!           v = cellVolume(thisOctal, subcell) * 1.d30
!           dfreq  = (phr(i)/ergtoev)/hcgs - (phl(i)/ergtoev)/hcgs

!           thisOctal%kromeintensity(subcell, i) = ((1.d0/v) * epsOverDt * &
!                thisOctal%kromeintensity(subcell, i) / fourPi)*ergtoEv / dFreq
!        enddo
     end if
  end do
end subroutine pathLengthToIntensity

recursive subroutine  zeropathlength(thisOctal)
  use krome_main
  use krome_user
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: i, subcell

  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call zeropathLength(child)
              exit
           end if
        end do
     else

        thisOctal%kromeintensity(subcell, :) = 0.d0
     end if
  end do
end subroutine zeroPathLength

recursive subroutine setMolAbundanceToKrome(thisOctal, label, isotopologueFraction)
  use krome_main
  use krome_user
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  character(len=*) :: label
  real(double) :: isotopologueFraction
  integer :: i, j, subcell

  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call setMolAbundanceToKrome(child, label, isotopologueFraction)
              exit
           end if
        end do
     else

        j = krome_get_index(label)
        thisOctal%molAbundance(subcell) = real(isotopologueFraction * &
             thisOctal%kromeSpeciesX(subcell,j) &
             / (thisOctal%rho(subcell)/(2.d0*mHydrogen)))
     end if
  end do
end subroutine setMolAbundanceToKrome


recursive subroutine  initializeChemistry(thisOctal)
  use inputs_mod, only : kromeInitialAbundances
  use krome_main
  use krome_user
  use krome_user_commons
  use constants_mod
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  real(double) :: thisX(krome_nMols), nH
  integer :: i, subcell

  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call initializeChemistry(child)
              exit
           end if
        end do
     else


        if (.not.associated(thisOctal%kromeSpeciesX)) then
           allocate(thisOctal%kromeSpeciesX(1:thisOctal%maxChildren, 1:krome_nmols))
        endif
        thisOctal%kromeSpeciesX(subcell, :) = 1.d-20

        nH = thisOctal%rho(subcell) / mHydrogen

        do i = 1, krome_nMols
           thisOctal%kromeSpeciesX(subcell,i)  = kromeInitialAbundances(i)
        enddo

       
        thisOctal%kromeSpeciesX(subcell,:) = thisOctal%kromeSpeciesX(subcell,:) * nH

        thisX = thisOctal%kromeSpeciesX(subcell,:)

        thisOctal%kromeSpeciesX(subcell,KROME_idx_e) = krome_get_electrons(thisX)

  !user commons for opacity and CR rate
  tau = 1d2 !opacity Av (#)
  zrate = 1.3d-17 !CR rate (1/s)
  gas_dust_ratio = 7.57d11 !gas/dust
  pah_size = 4d-8 !cm

   end if

  end do
end subroutine initializeChemistry

recursive subroutine doChemistryTimestepOctal(thisOctal, dt)
  use inputs_mod, only : amin, amax, qdist, dustPhysics, nDustType
  use krome_main
  use krome_user
  use krome_user_commons

  TYPE(OCTAL),pointer :: thisOctal, child
  real(double) :: dt, tgas
  real(double) :: thisX(krome_nmols)
  real(double), allocatable :: intensity(:)
  integer :: subcell,i
  
  do subcell = 1, thisOctal%maxChildren

     if (.not.thisOctal%hasChild(subcell)) then


        if ((thisOctal%temperature(subcell) < 300.d0).and.(thisOctal%rho(subcell)>1.d-29)) then
           thisX(1:krome_nmols) = thisOctal%kromeSpeciesx(subcell,1:krome_nmols)
           
           tau = thisOctal%meanAv(subcell) !opacity Av (#)
           
           
           zrate = 1.3d-17 !CR rate (1/s)
           gas_dust_ratio = 7.57d11 !gas/dust
           pah_size = 4d-8 !cm
           
           tgas = dble(thisOctal%temperature(subcell))

           
           if (dustPhysics) then
              if (ndustType > 1) then
                 call writeFatal("KROME cannot deal with more than one species of dust. Using the 1st dust type size distribution")
              endif
              !           call krome_init_dust_distribution(thisX, thisOctal%dustTypeFraction(subcell,1), alow_arg=dble(amin(1))*micronToCm, &
              !                aup_arg=dble(amax(1))*microntocm, &
              !                phi_arg=(-1.d0*dble(qdist(1))))
              !           call krome_set_Tdust(dble(thisOctal%temperature(subcell)))
           endif
           !        allocate(intensity(1:krome_nPhotoBins))
           !        intensity = tiny(intensity)
           !        intensity = intensity + thisOctal%kromeIntensity(subcell,1:krome_nPhotobins)
           
           !        call krome_set_PhotoBinJ(intensity)
           
           !        deallocate(intensity)
           if (dt > 0.d0) then
              call krome(thisX,  tgas, dt)
           endif
           thisOctal%kromeSpeciesx(subcell,1:krome_nmols) = thisX(1:krome_nmols)
        endif
     endif
  end do
end subroutine doChemistryTimestepOctal

recursive subroutine packKrome(thisOctal, index, kromeSpeciesX)
  use krome_main
  use krome_user
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: index
  real(double) :: kromeSpeciesX(:,:)
  integer :: i, subcell
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call packKrome(child, index, kromeSpeciesX)
              exit
           end if
        end do
     else
        index = index + 1
        kromeSpeciesX(index,1:krome_nmols) = thisOctal%kromeSpeciesx(subcell,1:krome_nmols)
     end if
  end do
end subroutine packKrome

recursive subroutine packAv(thisOctal, index, av)
  use krome_main
  use krome_user
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: index
  real(double) :: av(:)
  integer :: i, subcell
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call packAv(child, index, av)
              exit
           end if
        end do
     else
        index = index + 1
        av(index) = thisOctal%meanAv(subcell)
     end if
  end do
end subroutine packAv

recursive subroutine unpackKrome(thisOctal, index, kromeSpeciesX)
  use krome_main
  use krome_user
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: index
  real(double) :: kromeSpeciesX(:,:)
  integer :: i, subcell
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call unpackKrome(child, index, kromeSpeciesX)
              exit
           end if
        end do
     else
        index = index + 1
        thisOctal%kromeSpeciesx(subcell,1:krome_nmols) = kromeSpeciesX(index,1:krome_nmols) 
     end if
  end do
end subroutine unpackKrome

recursive subroutine unpackav(thisOctal, index, av)
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: index
  real(double) :: av(:)
  integer :: i, subcell
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call unpackAv(child, index, av)
              exit
           end if
        end do
     else
        index = index + 1
        thisOctal%meanAv(subcell) = av(index)
     end if
  end do
end subroutine unpackAv

#ifdef MPI
subroutine updateGridMPIkrome(grid, amrParComm)
  use amr_mod
  use krome_main
  use krome_user
  use gridtype_mod, only : gridtype
  use mpi
  implicit none
  integer :: amrParComm
  integer :: nOctals
  type(gridtype) :: grid
  integer :: nVoxels, i
  real(double), allocatable :: kromeSpeciesX(:,:)
  real(double), allocatable :: tempDoubleArray(:), thisTemp(:)
  integer :: ierr, nIndex

  ! FOR MPI IMPLEMENTATION=======================================================

  nVoxels = 0
  call countVoxels(grid%octreeRoot,nOctals,nVoxels)
  allocate(kromeSpeciesX(1:nVoxels, 1:krome_nmols))

  nIndex = 0
  call packKrome(grid%octreeRoot,nIndex, kromeSpeciesX)

  allocate(tempDoubleArray(nVoxels))
  allocate(thisTemp(nVoxels))

  do i = 1, krome_nmols
     tempDoubleArray = 0.d0
     thisTemp = kromeSpeciesX(1:nVoxels,i)
     call MPI_ALLREDUCE(thisTemp,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
          MPI_SUM, amrParComm ,ierr)
     kromeSpeciesX(1:nVoxels, i) = tempDoubleArray
  enddo


  call MPI_BARRIER(amrParComm, ierr) 

  nIndex = 0
  call unpackKrome(grid%octreeRoot, nIndex, kromeSpeciesX)
  deallocate(kromeSpeciesX, tempDoubleArray, thisTemp)

end subroutine updateGridMPIkrome

#endif

#ifdef MPI
subroutine updateGridMPIAV(grid, amrParComm)
  use amr_mod
  use gridtype_mod, only : gridtype
  use mpi
  implicit none
  integer :: amrParComm
  integer :: nOctals
  type(gridtype) :: grid
  integer :: nVoxels, i
  real(double), allocatable :: av(:)
  real(double), allocatable :: tempDoubleArray(:), thisTemp(:)
  integer :: ierr, nIndex

  ! FOR MPI IMPLEMENTATION=======================================================

  nVoxels = 0
  call countVoxels(grid%octreeRoot,nOctals,nVoxels)
  allocate(av(1:nVoxels))

  nIndex = 0
  call packAv(grid%octreeRoot,nIndex, av)

  allocate(thisTemp(nVoxels))

  call MPI_ALLREDUCE(av,thisTemp,nVoxels,MPI_DOUBLE_PRECISION,&
       MPI_SUM, amrParComm ,ierr)

  av = thisTemp
  call MPI_BARRIER(amrParComm, ierr) 

  nIndex = 0
  call unpackav(grid%octreeRoot, nIndex, av)
  deallocate(thisTemp)

end subroutine updateGridMPIAV

#endif

subroutine doChemistryOverGrid(grid, dt)
  use inputs_mod, only : splitovermpi
#ifdef MPI
  use mpi
  use mpi_global_mod
#endif
  use timing
  use gridtype_mod
  use amr_mod
  type(GRIDTYPE) :: grid
  type(OCTAL), pointer :: thisOctal
  real(double) :: dt
  integer :: nOctal, ioctal_beg, ioctal_end, ioctal
  type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals

#ifdef MPI
  integer :: np, n_rmdr, moctal
  integer :: ierr
#endif



  allocate(octalArray(grid%nOctals))
  nOctal = 0
  call getOctalArray(grid%octreeRoot,octalArray, nOctal)
  if (nOctal /= grid%nOctals) then
     write(*,*) "Screw up in get octal array", nOctal,grid%nOctals
     stop
  endif

  ! default loop indices
  ioctal_beg = 1
  ioctal_end = nOctal

  if (doTuning) call tune(6, "Chemistry step")


if (.not.splitovermpi) then
#ifdef MPI
  np = nThreadsGlobal
  n_rmdr = MOD(nOctal,np)
  mOctal = nOctal/np

  if (myRankGlobal .lt. n_rmdr ) then
     ioctal_beg = (mOctal+1)*myRankGlobal + 1
     ioctal_end = ioctal_beg + mOctal
  else
     ioctal_beg = mOctal*myRankGlobal + 1 + int(n_rmdr)
     ioctal_end = ioctal_beg + mOctal -1
  end if

#endif
endif

#ifdef MPI

  do iOctal = 1, nOctal
     thisOctal => octalArray(iOctal)%content
     if ((iOctal < iOctal_beg).or.(iOctal > iOctal_end)) then
        thisOctal%kromeSpeciesX = 0.d0
     endif
  enddo
#endif

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ioctal, thisOctal) &
  !$OMP SHARED(grid, ioctal_beg, ioctal_end, octalarray, dt, writeoutput, myrankGlobal)

  !$OMP DO SCHEDULE(DYNAMIC)
  do iOctal =  iOctal_beg, iOctal_end

     if ((ioctal_end-iOctal_beg) > 10) then
           if (mod(iOctal-iOctal_beg+1,(iOctal_end-iOctal_beg)/10)==0) then
              if (writeoutput) write(*,*) nint(100.*real(iOctal-iOctal_beg+1)/real(iOctal_end-iOctal_beg)), " % complete"
           endif
        endif

     thisOctal => octalArray(iOctal)%content

     write(*,*) myrankglobal, " doing octal ",ioctal, ioctal_beg, ioctal_end
     call doChemistryTimestepOctal(thisOctal, dt)

  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP END PARALLEL

#ifdef MPI

  if (.not.splitovermpi) then
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     call updateGridMPIkrome(grid, MPI_COMM_WORLD)
  endif

#endif
  deallocate(octalArray)

  if (doTuning) call tune(6, "Chemistry step")

end subroutine doChemistryOverGrid


subroutine initializeKrome()
  use krome_main
  use krome_user
  integer :: i

  call krome_init() !init krome (mandatory)

!  call krome_set_photoBin_J21log(5.d0, 13.6d0)
!  phl(:) = krome_get_photoBinE_left() !returns left bin points
!  phr(:) = krome_get_photoBinE_right() !returns right bin points
!  phm(:) = krome_get_photoBinE_mid() !returns left middle points
!  phJ(:) = krome_get_photoBinJ() !returns bin intensities
!  if (writeoutput) then
!     call writeBanner("Krome photo bins (left, mid, right) ev, intensity","+")
!     do i=1,krome_nPhotoBins
!        print '(I5,1p,4E10.2)',i,phl(i),phm(i),phr(i),phJ(i)
!     end do
!  endif

  call reportKromeSpecies()

end subroutine initializeKrome

subroutine storePathlength(thisOctal, subcell, frequency, pathlength)
  type(OCTAL), pointer :: thisOctal
  integer :: subcell
  real(double) :: pathlength, frequency, energy
  integer :: i

  energy = frequency * hCgs * ergtoev

!  if ((energy >= phl(1)).and.(energy <= phr(krome_nPhotoBins))) then
!     do i = 1, krome_nPhotoBins
!        if ((energy >= phl(i)).and.(energy <= phr(i))) then
!           thisOctal%kromeIntensity(subcell,i) = thisOctal%kromeintensity(subcell,i) + pathLength
!           exit
!        endif
!     enddo
!  endif

end subroutine storePathlength



subroutine doChem(grid, dt)
  use gridtype_mod
  use krome_main !use krome (mandatory)
  use krome_user !use utility (for krome_idx_* constants and others)
  use vtk_mod
  implicit none
  type(GRIDTYPE) :: grid
  integer :: i, subcell
  type(OCTAL), pointer :: thisOctal
!  character(len=16) :: species(krome_nmols)
  character(len=80) :: vtkFilename
  real(double) :: dt,totaltime

  call initializeChemistry(grid%octreeRoot)
  call calculateAv(grid)

  write(vtkFilename,"(a,i3.3,a)") "chem",0,".vtk"
  call writeVtkFile(grid, vtkFilename, &
       valueTypeString=(/"chemistry  ","rho        ","av         ","temperature"/))

  call doChemistryoverGrid(grid, dt)


  totaltime = 0.d0
!  dt = 1.d2 * yearstosecs
!  open(20,file="chem.dat",status="unknown",form="formatted")
!  do while (totalTime < 1.d8*yearstosecs)
!    call writeInfo("Doing chemistry timestep...",TRIVIAL)
!    write(*,*) "doing chemistry at ",totaltime*secstoyears
!     call doChemistryoverGrid(grid, dt)
!     thisOctal => grid%octreeRoot
!     call findSubcellLocal(vector(0.d0, 0.0d0, 0.d0),thisOctal,subcell)
!     if (totaltime > 0.) then
!        write(20,'(5f12.3)') log10(totaltime*secstoyears), log10(thisOctal%kromeSpeciesX(subcell,krome_idx_Cj)/1.d4), &
!          log10(thisOctal%kromeSpeciesX(subcell,krome_idx_HCO)/1.d4), log10(thisOctal%kromeSpeciesX(subcell,krome_idx_CO)/1.d4), &
!          log10(thisOctal%kromeSpeciesX(subcell,krome_idx_H2O)/1.d4) 
!     endif
!     totaltime =totaltime+dt
!     dt = max(1.d2*yearstosecs,totaltime/3.d0)
     call writeInfo("Done.",TRIVIAL)

!  enddo
!  close(20)
     write(vtkFilename,"(a,i3.3,a)") "chem",1,".vtk"
     call writeVtkFile(grid, vtkFilename, &
          valueTypeString=(/"chemistry  ","rho        ","av         ","temperature"/))


end subroutine doChem


subroutine reportKromeSpecies()
  use krome_user
  use inputs_mod, only : kromeInitialAbundances
  character(len=16) :: species(krome_nmols)
  integer :: i
  species = krome_get_names()
  if (writeoutput) then
     call writebanner("Krome species in reaction network","+")
     do i = 1, krome_nmols
        write(*,'(i3,1x,a16,1p,e12.3)') i, species(i), kromeinitialAbundances(i)
     enddo
  endif
end subroutine reportKromeSpecies

subroutine calculateAv(grid)
  use gridtype_mod
  use utils_mod, only : locate
  use timing
#ifdef MPI
  use mpi
#endif
  use amr_mod, only : getOctalArray
  type(GRIDTYPE) :: grid
  integer :: iLambda


  type(OCTAL), pointer :: thisOctal
  integer :: nDir
  type(VECTOR), pointer :: dir(:)
  integer :: nOctal, ioctal_beg, ioctal_end, ioctal
  type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals

#ifdef MPI
  integer :: np, n_rmdr, moctal
  integer :: ierr
#endif

  call locate(grid%lamArray, grid%nLambda, 5500., iLambda)


  call getUniformSphereDirections(1, ndir, dir)

  allocate(octalArray(grid%nOctals))
  nOctal = 0
  call getOctalArray(grid%octreeRoot,octalArray, nOctal)
  if (nOctal /= grid%nOctals) then
     write(*,*) "Screw up in get octal array", nOctal,grid%nOctals
     stop
  endif

  ! default loop indices
  ioctal_beg = 1
  ioctal_end = nOctal

  if (doTuning) call tune(6, "A_V step")


#ifdef MPI
  np = nThreadsGlobal
  n_rmdr = MOD(nOctal,np)
  mOctal = nOctal/np

  if (myRankGlobal .lt. n_rmdr ) then
     ioctal_beg = (mOctal+1)*myRankGlobal + 1
     ioctal_end = ioctal_beg + mOctal
  else
     ioctal_beg = mOctal*myRankGlobal + 1 + int(n_rmdr)
     ioctal_end = ioctal_beg + mOctal -1
  end if

#endif

#ifdef MPI

  do iOctal = 1, nOctal
     thisOctal => octalArray(iOctal)%content
     if ((iOctal < iOctal_beg).or.(iOctal > iOctal_end)) then
        if (.not.associated(thisOctal%meanAv)) allocate(thisOctal%meanAv(1:thisOctal%maxChildren))

        thisOctal%meanAv = 0.d0
     endif
  enddo
#endif

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ioctal, thisOctal) &
  !$OMP SHARED(grid, ioctal_beg, ioctal_end, octalarray, writeoutput, ilambda, dir, ndir)

  !$OMP DO SCHEDULE(DYNAMIC)
  do iOctal =  iOctal_beg, iOctal_end

     if ((ioctal_end-iOctal_beg) > 10) then
           if (mod(iOctal-iOctal_beg+1,(iOctal_end-iOctal_beg)/10)==0) then
              if (writeoutput) write(*,*) nint(100.*real(iOctal-iOctal_beg+1)/real(iOctal_end-iOctal_beg)), " % complete"
           endif
        endif

     thisOctal => octalArray(iOctal)%content

     call calculateAvOctal(grid, thisOctal, iLambda,ndir, dir)

  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP END PARALLEL

#ifdef MPI

     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     call updateGridMPIAV(grid, MPI_COMM_WORLD)

#endif
  deallocate(octalArray)
  deallocate(dir)


  if (doTuning) call tune(6, "A_V step")


end subroutine calculateAv

recursive subroutine calculateAVOctal(grid, thisOctal, ilambda,ndir,dir)
  use amr_mod, only : tauAlongPathFast
  use gridtype_mod
  type(GRIDTYPE) :: grid
  type(OCTAL), pointer :: thisOctal, child
  integer :: iDir
  integer :: index, ilambda
  integer :: ndir
  type(VECTOR) :: dir(:), cellCentre
  real(double), allocatable :: tau(:), Av(:), expMinusAv(:)

  integer :: i, subcell
  
  allocate(tau(1:nDir), av(1:nDir), expMinusAv(1:nDir))

  do subcell = 1, thisOctal%maxChildren

     if (.not.thisOctal%hasChild(subcell)) then

        cellCentre = subcellCentre(thisOctal, subcell)


        do iDir = 1, nDir
           call tauAlongPathFast(ilambda, grid, cellCentre, dir(idir), tau(idir), startOctal=thisOctal, startSubcell=subcell)
           av(idir) = 1.086d0 * tau(idir)
        enddo
        expMinusAv(1:nDir) = exp(-av(1:nDir))
        
        if (.not.associated(thisOctal%meanAv)) allocate(thisOctal%meanAv(1:thisOctal%maxChildren))
        
        thisOctal%meanAv(subcell) = -log((1.d0/dble(nDir)) * SUM(expminusav(1:nDir)))
     endif
  end do
  deallocate(tau, av, expMinusAv)
end subroutine calculateAVOctal



end module chemistry_mod
#endif

