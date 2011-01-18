module wr104_mod

  use kind_mod
  use sph_data_class
  use vector_mod
  use cluster_class

  implicit none

  type particle_list
     type(VECTOR), pointer :: position(:)
     integer :: nParticles
  end type particle_list

contains

  subroutine readWR104Particles(filename, sphdata, objectDistance)
    use constants_mod, only: degtorad
    character(len=*), intent(in) :: filename
    type(sph_data), intent(inout) :: sphData
    real(double), intent(in) :: objectDistance
    !
!    integer :: i
    real(double) :: udist, umass, utime    ! Units of distance, mass, time
!    integer          :: npart                  ! Number of gas particles
    real(double) :: time                   ! Time of sph data dump (in units of utime)
    integer          :: nptmass                ! Number of stars/brown dwarfs

!    real(double) :: dummy  ! position of gas particles



    write(*,'(a,a)') "Reading particles from: ",trim(filename)
    open(20,file=filename,status="old",form="unformatted")
    read(20) npart
    write(*,*) "Reading ",npart," particles"
! assume that the units are 1 mil
    uDist = 1.e-3/3600. ! from milliarcsec to degrees
    uDist = sphData%uDist * degtorad ! to radians
    uDist = sphData%uDist * objectDistance ! cm 

    uMass = 1.
    uTime = 1.e6
    time = 1.e6
    nptmass = 1

    ! initilaizing the sph_data object
    ! using a routine in sph_data_class.  (Allocating the memory for arrays and etc..)
    call init_sph_data(udist, umass, utime, time, nptmass)


    ! reading the X poistion of particles
!    do i = 1, npart
!       read(20) dummy
!       call put_position_gas_particle(sphData, i, "x", dummy) ! saving the value
!    end do

    ! reading the y poistion of particles
!    do i = 1, npart
!       read(20) dummy
!       call put_position_gas_particle(sphData, i, "y", dummy) ! saving the value
!    end do
    ! reading the z poistion of particles
!    do i = 1, npart
!       read(20) dummy
!       call put_position_gas_particle(sphData, i, "z", dummy) ! saving the value
!    end do

    ALLOCATE(sphData%xn(npart))
    ALLOCATE(sphData%yn(npart))
    ALLOCATE(sphData%zn(npart))
    ALLOCATE(sphData%rhon(npart))
    ALLOCATE(sphData%x(1),sphdata%y(1),sphdata%z(1))
    ALLOCATE(sphData%ptmass(1))

! assume that the units are 1 mil

    sphData%uDist = 1.e-3/3600. ! from milliarcsec to degrees
    sphData%uDist = sphData%uDist * degtorad ! to radians
    sphData%uDist = 1.e10 !sphData%uDist * objectDistance ! cm 
    
    sphData%uMass = 1.
    sphData%uTime = 1.e6


    sphData%npart = npart
    sphData%rhon = 1.
    read(20) sphData%xn(1:npart), sphdata%yn(1:npart), sphData%zn(1:npart)
    close(20)
    write(*,'(a)') "Done."

  end subroutine readWR104Particles

  subroutine assign_grid_values_wr104(thisOctal, subcell)
    type(octal), pointer  :: thisOctal
    real(double) :: rhoAv
    integer :: subcell, i
    
    call find_n_particle_in_subcell(i, rhoav, thisOctal, subcell)
    thisOctal%rho(subcell) = (max(dble(i),1.d-10))/cellVolume(thisOctal, subcell)

  end subroutine assign_grid_values_wr104

end module wr104_mod
