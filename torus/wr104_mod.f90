module wr104_mod

  use kind_mod
  use sph_data_class
  use vector_mod

  implicit none

  type particle_list
     type(VECTOR), pointer :: position(:)
     integer :: nParticles
  end type particle_list

contains

  subroutine readWR104Particles(filename, sphData)
    type(sph_data), intent(inout) :: sphData
    character(len=*) :: filename
    integer :: npart,i


    write(*,'(a,a)') "Reading particles from: ",trim(filename)
    open(20,file=filename,status="old",form="unformatted")
    read(20) npart


    ALLOCATE(sphData%xn(npart))
    ALLOCATE(sphData%yn(npart))
    ALLOCATE(sphData%zn(npart))
    ALLOCATE(sphData%rhon(npart))
    ALLOCATE(sphData%x(1),sphdata%y(1),sphdata%z(1))
    ALLOCATE(sphData%ptmass(1))
    sphData%uDist = 200. * autocm / 300.
    sphData%uMass = 1.
    sphData%uTime = 1.e6


    sphData%npart = npart
    sphData%rhon = 1.
    read(20) sphData%xn(1:npart), sphdata%yn(1:npart), sphData%zn(1:npart)
    close(20)
    write(*,'(a)') "Done."

  end subroutine readWR104Particles





end module wr104_mod
