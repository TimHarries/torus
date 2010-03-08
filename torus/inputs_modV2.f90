module inputs_mod

  use messages_mod
  use kind_mod
  use FoX_sax

contains
  
  subroutine characters_handler(chars)
    character(len=*), intent(in) :: chars
    
    print*, chars
  end subroutine characters_handler

  subroutine  inputs()
    use FoX_sax
!    use event_handling
    type(xml_t) :: xp
    call open_xml_file(xp, 'input.xml')
    call parse(xp, characters_handler=characters_handler)
    call close_xml_t(xp)
  end subroutine inputs




end module inputs_mod

