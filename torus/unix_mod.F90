module unix_mod
  
  implicit none
  
contains

  subroutine unixGetenv(envVariable, directory, length, error)

#ifdef USEUNIXENV
    use f90_unix_env
#endif
    character(len=*) :: envVariable
    character(len=*) :: directory
    integer, intent(out), optional :: length
    integer, intent(out), optional :: error

    call getenv(envVariable, directory)
    if (present(error)) error = 0
    if (present(length)) length = len(directory)

  end subroutine unixGetenv

  subroutine unixGetLogin(login)
    
!    use f90_unix_env
    
    character(len=80) :: login
    integer :: n
!    call getLogin(login)
    n = len(trim(login))
    login = " "

  end subroutine unixGetLogin

  subroutine unixGetHostname(name)

!    use f90_unix_env
    
    character(len=*) :: name
!    integer :: n
    
    name = " "
!    call getHostname(name, n)

  end subroutine unixGetHostname

    


end module unix_mod

