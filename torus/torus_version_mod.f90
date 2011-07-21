module torus_version_mod

  contains
    
    subroutine setVersion(v)
      use constants_mod, only : torusVersion
      include "svn_version.h"
      character(len=*) :: v
      torusVersion = trim(v)//trim(svnversion)
    end subroutine setVersion

end module torus_version_mod
