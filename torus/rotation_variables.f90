MODULE rotation_variables

  ! in a module so that they can be shared by funcA7 (it's passed as
  !   an argument).

  IMPLICIT NONE
  PUBLIC 

    REAL :: bigR, r
    REAL :: Einitial
    REAL :: Bp
    REAL :: l
    REAL :: omega
    REAL :: omegaStar ! 
    REAL :: eta

END MODULE rotation_variables
