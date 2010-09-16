module kind_mod
  ! KIND parameters
  !
  ! Single and double are defined to be consistent with IEEE 754 definitions. 
  ! On some systems higher precisions are available (10 or 16 bytes) but 
  ! this is compiler dependent and will cause a performace hit if 
  ! implemented in software. Consequently oct is set to double by default. 

  implicit none

  integer, parameter :: single = selected_real_kind(6, 37)
  integer, parameter :: double = selected_real_kind(15, 307)
  integer, parameter :: bigInt = selected_int_kind(10)

  integer, parameter :: si = single
  integer, parameter :: db = double

! 16 byte quad precision
!  integer, parameter :: oct = selected_real_kind(33, 4931) 
  integer, parameter :: oct = double
  integer, parameter :: oc = oct

end module kind_mod

