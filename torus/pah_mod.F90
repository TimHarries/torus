module pah_mod

  implicit none
  use constants_mod

  type PAHTABLETYPE
     integer :: nu
     real, allocatable :: u(:)
  end type PAHTABLETYPE

contains

  subroutine readPAHEmissivityspectra(PAHtable, PAHtype)
    character(len=*) :: PAHtype
    type(PAHTableTYPE) :: PAHTable

    PAHtable%nU = 30
    allocate(PAHtable%U(1:PAHtable%nU))
    PAHtable%U(01) = 0.10
    PAHtable%U(02) = 0.15
    PAHtable%U(03) = 0.20
    PAHtable%U(04) = 0.30
    PAHtable%U(05) = 0.40
    PAHtable%U(06) = 0.50
    PAHtable%U(07) = 0.70
    PAHtable%U(08) = 0.80
    PAHtable%U(09) = 1.00
    PAHtable%U(10) = 1.20
    PAHtable%U(11) = 1.50
    PAHtable%U(12) = 2.00
    PAHtable%U(13) = 2.50
    PAHtable%U(14) = 3.00
    PAHtable%U(15) = 4.00
    PAHtable%U(16) = 5.00
    PAHtable%U(17) = 7.00
    PAHtable%U(18) = 8.00
    PAHtable%U(19) = 12.0
    PAHtable%U(20) = 15.0
    PAHtable%U(21) = 20.0
    PAHtable%U(22) = 25.0
    PAHtable%U(23) = 1e2
    PAHtable%U(24) = 3e2
    PAHtable%U(25) = 1e3
    PAHtable%U(26) = 3e3
    PAHtable%U(27) = 1e4
    PAHtable%U(28) = 3e4
    PAHtable%U(29) = 1e5
    PAHtable%U(30) = 3e5
   
    
    do i = 1, PAHtable%nu

       if (PAHtable%u(i) <= 10.d0) then
          write(cval,'(f4.2)') PAHtable%u(i)
       else if ((PAHtable%u(i) > 10.d0).and.(PAHtable%u(i) < 100.d0)) then
          write(cval, '(f4.1)') PAHtable%u(i)
       else
          write(cval, '(1p,e3.0)') PAHtable%u(i)
       endif
       write(filename,'(a,a,a,a,a,a,a,a,a,a,a)') &
            "U",trim(cval),"/","U",trim(cval),"_", &
            "U",trim(cval),"_",trim(PAHtype),".txt"
       write(*,*) trim(filename)
    enddo
    
  end subroutine readPAHEmissivityspectra



end module pah_mod
