module zlib_mod

  use kind_mod
  use, intrinsic :: ISO_C_BINDING
  implicit none

#ifdef USEZLIB


  interface
     function compressBound(sourceLen) bind(C, name = 'compressBound')
       use, intrinsic :: ISO_C_BINDING
       integer(c_long), value :: sourceLen
       integer(c_long) :: compressBound
     end function compressBound
  end interface

  interface
     function compress(dest, destLen, source, sourceLen) bind(C, name = 'compress')
       use, intrinsic :: ISO_C_BINDING
       type(c_ptr), value :: dest, source, destLen
       integer(c_long), value ::  sourceLen
       integer(c_int) :: compress
     end function compress
  end interface

  interface
     function uncompress(dest, destLen, source, sourceLen) bind(C, name = 'uncompress')
       use, intrinsic :: ISO_C_BINDING
       type(c_ptr), value :: dest, source, destLen
       integer(c_long), value ::  sourceLen
       integer(c_int) :: uncompress
     end function uncompress
  end interface

  public

contains

  subroutine compressBytes(inArray, nIn, outArray, nOut)
    integer(kind=1), pointer :: inArray(:)
    integer(bigint) :: nIn, nOut
    integer :: nInFourByte
    integer(kind=1), pointer :: outArray(:)
    type(c_ptr) :: c_outlen
    integer :: j, iMax
    integer(c_long), target :: outLen
    integer(kind=1), pointer :: inBuffer(:), outBuffer(:)
    type(c_ptr) :: c_inbuffer = c_null_ptr, c_outbuffer = c_null_ptr

    allocate(inBuffer(1:nIn))
    inBuffer(1:nIn) = inArray(1:nIn)

    iMax = int(compressBound(int(nIn, kind=c_long)))
    allocate(outBuffer(1:iMax))

    c_inbuffer = c_loc(inBuffer(1))
    c_outbuffer = c_loc(outBuffer(1))
    outlen = iMax
    c_outlen = c_loc(outlen)


    j = compress(c_outbuffer, c_outLen, c_inbuffer, int(nIn, kind=c_long))
       
    nOut = outLen
    allocate(outArray(1:nOut))
    outArray(1:nOut) = outbuffer(1:nOut)
    deallocate(inBuffer, outBuffer)
  end subroutine compressBytes

#endif


end module zlib_mod
