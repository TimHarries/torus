module torus_version_mod

  contains
    
    subroutine setVersion(v)
      use constants_mod, only : torusVersion
      include "git_version.h"
      character(len=*) :: v
      torusVersion = trim(v)//trim(gitversion)
    end subroutine setVersion

end module torus_version_mod
