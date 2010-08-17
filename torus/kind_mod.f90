module kind_mod
  ! KIND parameters, specific to the individual system
  ! double, octal and logic entries added by nhs

  implicit none

  integer, parameter :: single = kind(0.)
  integer, parameter :: bigInt = kind(1000000000)
  integer, parameter :: si = single
  integer, parameter :: double = kind(0.d0)
  integer, parameter :: db = double
  integer, parameter :: quad = kind(0.d0)
  integer, parameter :: qd = quad
  integer, parameter :: oct = kind(0.d0)
  integer, parameter :: oc = oct

  integer, parameter :: logic = kind(.true.)
  integer, parameter :: lg = logic

end module kind_mod
