module parameters_mod    
   ! parameters for some of the geometry types

   use constants_mod
  
   implicit none

   public

   ! TTauri star and disk
   ! we may need to track the minimum and maximum densities of the
   !   accretion flow 
! RK COMMENTED OUT THIS. IBM XFL compiler does not allow intrinsic function 
! in declearations.
!   real, save      :: TTauriMinRho = huge(TTauriMinRho)
!   real, save      :: TTauriMaxRho = tiny(TTauriMaxRho)
   real, save      :: TTauriMinRho = 1.0e25
   real, save      :: TTauriMaxRho = 1.0e-25

end module parameters_mod
