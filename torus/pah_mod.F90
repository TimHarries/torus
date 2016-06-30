module pah_mod

  use constants_mod
  use unix_mod
  use messages_mod
  implicit none

  type PAHTABLETYPE
     integer :: nu
     real, allocatable :: u(:)
     real, allocatable :: lambda(:)
     real, allocatable :: jnu(:,:)
  end type PAHTABLETYPE

contains

  subroutine readPAHEmissivityTable(PAHtable, PAHtype)
    character(len=*) :: PAHtype
    type(PAHTableTYPE) :: PAHTable
    integer :: i, j 
    character(len=10) :: cval
    character(len=80) :: filename, dataDirectory, cjunk
    real :: vjunk

    call unixGetenv("TORUS_DATA", dataDirectory, i)

    PAHtable%nU = 30
    allocate(PAHtable%U(1:PAHtable%nU), PAHtable%lambda(1001), PAHtable%jnu(1:PAHtable%nu,1:1001))
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
          write(cval,'(f5.2)') PAHtable%u(i)
       else if ((PAHtable%u(i) > 10.d0).and.(PAHtable%u(i) < 100.d0)) then
          write(cval, '(f4.1)') PAHtable%u(i)
       else
          if (i == 23) cval="1e2"
          if (i == 24) cval="3e2"
          if (i == 25) cval="1e3"
          if (i == 26) cval="3e3"
          if (i == 27) cval="1e4"
          if (i == 28) cval="3e4"
          if (i == 29) cval="1e5"
          if (i == 30) cval="3e5"
       endif
       write(filename,'(a,a,a,a,a,a,a,a,a,a,a,a)') &
            trim(dataDirectory),"/PAH/", &
            "U",trim(cval),"/U",trim(cval),"_", &
            trim(cval),"_",trim(PAHtype),".txt"

       open(20,file=filename,status="old",form="formatted")
       do j = 1, 61
          read(20,'(a)') cjunk
       enddo
       do j = 1, 1001
          read(20,*) PAHtable%lambda(1002-j), vjunk, PAHtable%jnu(i,1002-j)
          if (writeoutput) write(*,*) "lambda ",PAHtable%lambda(j), PAHtable%jnu(i,1002-j)
       enddo
       close(20)
    enddo
  end subroutine readPAHEmissivityTable
  


end module pah_mod
