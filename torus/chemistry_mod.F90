#ifdef CHEMISTRY
module chemistry_mod

  use kind_mod
  use octal_mod
  use messages_mod
  implicit none

contains

recursive subroutine  initializeChemistry(thisOctal)
  use krome_main
  use krome_user
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: i, subcell
  real(double) :: nh

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
        thisOctal%kromeSpeciesX(subcell, :) = 1.d-5
        thisOctal%kromeSpeciesX(subcell,krome_idx_He) = 9d-2
        thisOctal%kromeSpeciesX(subcell,krome_idx_CO) = 1d-5
        thisOctal%kromeSpeciesX(subcell,krome_idx_O) = 2.56d-4

     end if

  end do
end subroutine initializeChemistry

recursive subroutine doChemistryTimestep(thisOctal, dt)
  use krome_main
  use krome_user
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  real(double) :: dt
  real(double) :: thisX(krome_nmols)
  integer :: i, subcell
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call doChemistryTimestep(child, dt)
              exit
           end if
        end do
     else
        thisX(1:krome_nmols) = thisOctal%kromeSpeciesx(subcell,1:krome_nmols)
        call krome(thisX, thisOctal%rho(subcell), dble(thisOctal%temperature(subcell)), dt)
        thisOctal%kromeSpeciesx(subcell,1:krome_nmols) = thisX(1:krome_nmols)

     end if

  end do
end subroutine doChemistryTimestep

subroutine initializeKrome()
  use krome_main
  use krome_user

  call krome_init() !init krome (mandatory)

  call reportKromeSpecies()
end subroutine initializeKrome

subroutine doChem(grid)
  use gridtype_mod
  use krome_main !use krome (mandatory)
  use krome_user !use utility (for krome_idx_* constants and others)
  use vtk_mod
  implicit none
  type(GRIDTYPE) :: grid
  integer :: i
  character(len=16) :: species(krome_nmols)
  character(len=80) :: vtkFilename
  real(double) :: dt

  call initializeChemistry(grid%octreeRoot)

  dt = 1.d9
  do i = 0,10
     write(vtkFilename,"(a,i3.3,a)") "chem",i,".vtk"
     call writeVtkFile(grid, vtkFilename, &
          valueTypeString=(/"chemistry  ","rho        ","temperature"/))
     call writeInfo("Doing chemistry timestep...",TRIVIAL)
     call doChemistryTimestep(grid%octreeRoot,dt)
     call writeInfo("Done.",TRIVIAL)
  enddo


end subroutine doChem


subroutine reportKromeSpecies()
  use krome_user
  character(len=16) :: species(krome_nmols)
  integer :: i
  species = krome_get_names()
  call writebanner("Krome species in reaction network","+")
  do i = 1, krome_nmols
     write(*,'(i3,1x,a16)') i, species(i)
  enddo
end subroutine reportKromeSpecies


end module chemistry_mod
#endif

