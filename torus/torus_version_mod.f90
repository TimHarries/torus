module torus_version_mod

  character(len=10) :: torusVersion

  contains
    
    subroutine setVersion(v)
      character(len=*) :: v
      torusVersion = v
    end subroutine setVersion

end module torus_version_mod
