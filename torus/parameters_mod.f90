module parameters_mod    
   ! parameters for some of the geometry types

   use constants_mod
  
   implicit none

   public

   ! TTauri star and disk
   real, parameter :: TTauriRstar      = 2.0 * rSol 
   real, parameter :: TTauriMstar      = 0.8 * mSol
   real, parameter :: TTauriRinner     = 2.2 * TTauriRstar
   real, parameter :: TTauriRouter     = 3.0 * TTauriRstar
   real, parameter :: TTauriDiskRadius = 3.0 * TTauriRstar
   real, parameter :: TTauriMdot       = 1.e-7 * mSol / (365.25 * 24. * 3600.)
   real, parameter :: TTauriMloss      = 1.e-8 * mSol / (365.25 * 24. * 3600.)

end module parameters_mod
