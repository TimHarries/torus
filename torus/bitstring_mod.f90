module bitstring_mod

  use kind_mod
  use octal_mod

  implicit none

  contains

    recursive subroutine constructBitStrings(thisOctal)
      type(octal), pointer   :: thisOctal
      type(octal), pointer  :: child 
      integer :: subcell, i
      integer :: xBit, yBit, zBit


      do subcell = 1, thisOctal%maxChildren

         if (thisOctal%nDepth > 1) then
            thisOctal%iXbitString(subcell) = thisOctal%parent%ixBitstring(thisOctal%parentSubcell)*2
            thisOctal%iYbitString(subcell) = thisOctal%parent%iyBitstring(thisOctal%parentSubcell)*2
            thisOctal%iZbitString(subcell) = thisOctal%parent%izBitstring(thisOctal%parentSubcell)*2
         endif
         ! left is zero, right is one
         ! down is zero, up is one

         if (subcell == 1) then 
            xBit = 0; yBit = 0; zBit =  0
         endif

         if (subcell == 2) then 
            xBit = 1; yBit = 0; zBit =  0
         endif

         if (subcell == 3) then 
            xBit = 0; yBit = 0; zBit =  1
         endif

         if (subcell == 4) then 
            xBit = 1; yBit = 0; zBit =  1
         endif

         if (subcell == 5) then 
            xBit = 0; yBit = 1; zBit =  0
         endif

         if (subcell == 6) then 
            xBit = 1; yBit = 1; zBit =  0
         endif

         if (subcell == 7) then 
            xBit = 0; yBit = 1; zBit =  1
         endif

         if (subcell == 8) then 
            xBit = 1; yBit = 1; zBit =  1
         endif

         if (xBit == 0) then
            thisOctal%iXbitstring(subcell) = ibclr(thisOctal%iXbitstring(subcell), 0)
         else 
            thisOctal%iXbitstring(subcell) = ibset(thisOctal%iXbitstring(subcell), 0)
         endif

         if (thisOctal%threed) then
            if (yBit == 0) then
               thisOctal%iYbitstring(subcell) = ibclr(thisOctal%iYbitstring(subcell), 0)
            else 
               thisOctal%iYbitstring(subcell) = ibset(thisOctal%iYbitstring(subcell), 0)
            endif
         endif

         if (zBit == 0) then
            thisOctal%iZbitstring(subcell) = ibclr(thisOctal%iZbitstring(subcell), 0)
         else 
            thisOctal%iZbitstring(subcell) = ibset(thisOctal%iZbitstring(subcell), 0)
         endif

         if (thisOctal%nDepth == 10) then
            write(*,*) "x ",charBitstring(thisOctal%ixBitstring(subcell), thisOctal%nDepth), &
                 " z ", charBitString(thisOctal%izBitstring(subcell), thisOctal%nDepth)
         endif
      end do




      IF ( thisOctal%nChildren > 0 ) THEN
         ! call this subroutine recursively on each of its children
         DO i = 1, thisOctal%nChildren, 1
            child => thisOctal%child(i)
            CALL constructBitStrings(child)
         END DO
      END IF

    end subroutine constructBitStrings

    character(len=32) function  charBitString(iBitString,depth) result(number)
      integer(kind=bigInt) :: iBitstring
      integer :: i, depth
      number = " "
      do i = depth-1, 0, -1
         if (btest(ibitString,i)) then
            number(32-i:32-i) = "1"
         else
            number(32-i:32-i) = "0"
         endif
      enddo
    end function charBitString

end module bitstring_mod
