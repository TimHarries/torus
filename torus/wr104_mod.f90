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
    integer :: i
    double precision :: udist, umass, utime    ! Units of distance, mass, time
    integer          :: npart                  ! Number of gas particles
    double precision :: time                   ! Time of sph data dump (in units of utime)
    integer          :: nptmass                ! Number of stars/brown dwarfs
    double precision :: gaspartmass            ! Mass of each gas particles

    double precision :: dummy  ! position of gas particles



    write(*,'(a,a)') "Reading particles from: ",trim(filename)
    open(20,file=filename,status="old",form="unformatted")
    read(20) npart

    uDist = 200. * autocm / 300.
    uMass = 1.
    uTime = 1.e6
    time = 1.e6
    nptmass = 1
    gaspartmass = 1.0d-24

    ! initilaizing the sph_data object
    ! using a routine in sph_data_class.  (Allocating the memory for arrays and etc..)
    call init_sph_data(sphData, udist, umass, utime, npart, time, nptmass, gaspartmass)


    ! reading the X poistion of particles
    do i = 1, npart
       read(20) dummy
       call put_position_gas_particle(sphData, i, "x", dummy) ! saving the value
    end do

    ! reading the y poistion of particles
    do i = 1, npart
       read(20) dummy
       call put_position_gas_particle(sphData, i, "y", dummy) ! saving the value
    end do
    ! reading the z poistion of particles
    do i = 1, npart
       read(20) dummy
       call put_position_gas_particle(sphData, i, "z", dummy) ! saving the value
    end do

    close(20)

    ! setting the density to one everywhere.. 
    do i = 1, npart
       call put_rhon(sphData, i, 1.0d0)
    end do


    write(*,'(a)') "Done."

  end subroutine readWR104Particles





end module wr104_mod
