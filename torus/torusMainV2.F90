!
!            T   O   R   U  S
!
! torus is a general purpose monte-carlo radiative transfer code used to
! produce polarized spectra of hot and cool stellar winds
!

! written by tjh
! broken by chris

! v1.0 on 16/09/99

! raman scattering stuff added 3/3/2000

! OMP parallelization calls added 1/7/2001

! TJH: 19/7/02  Jorick's spotty star stuff added
! NHS: 02/9/02  adaptive mesh code merged
! RK : 12/02/03 Parallelized major loops in lucyStateEquilibriumAMR, stateqAMR, torusMain.
! RK : 07/07/05 Added the observed flux solver in oberver's frame. 
! TJH: 10/09/08 pgplot calls removed...
program torus
  use torus_version_mod
  use input_variables         ! variables filled by inputs subroutine
  use constants_mod
  use messages_mod
  use mpi_global_mod
  use inputs_mod
  
  call inputs()

end program torus
