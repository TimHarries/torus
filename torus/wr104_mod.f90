module wr104_mod

  use kind_mod
  use vector_mod

  type particle_list
     type(VECTOR), pointer :: position(:)
     integer :: nParticles
  end type particle_list
  

contains

  subroutine readParticles(filename, particles)
    character(len=*) :: filename
    type(particle_list) :: particles
    integer :: nParticles, i

    open(20, file=filename, status="old", form="formatted")
    nParticles = 100000
    write(*,'(a,i12,a)') "Reading ", nParticles, " particles..."
    allocate(particles%position(1:nParticles))
    do i = 1, nParticles
       read(20,*) particles%position(i)%x,particles%position(i)%y,particles%position(i)%z

!       particles%position(i) = particles%position(i) * 1000. * pctocm * 1.e3 * 1.e-10
    enddo
    close(20)
    particles%nParticles = nParticles
    write(*,'(a)') "Done."
  end subroutine readParticles


  subroutine wr104CountParticles(nParticles, particles, cellCentre, cellSize)
    type(particle_list) :: particles
    integer :: i, nParticles
    type(octalVector) :: cellCentre
    real(kind=doublekind) ::cellSize
    nParticles = 0
    do i = 1, particles%nParticles
       if ((particles%position(i)%x < cellcentre%x+cellSize/2.) .and. &
            (particles%position(i)%x > cellcentre%x-cellSize/2.) .and. &
            (particles%position(i)%y < cellcentre%y+cellSize/2.) .and. &
            (particles%position(i)%y > cellcentre%y-cellSize/2.) .and. &
            (particles%position(i)%z < cellcentre%z+cellSize/2.) .and. &
            (particles%position(i)%z > cellcentre%z-cellSize/2.)) then
          nParticles = nParticles + 1
       endif
    enddo
  end subroutine wr104CountParticles


end module wr104_mod
