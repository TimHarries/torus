module zlib_mod

  use, intrinsic :: ISO_C_BINDING
  implicit none

#ifdef USEZLIB


  interface
     function compress(dest, destLen, source, sourceLen) bind(C, name = 'compress')
       use, intrinsic :: ISO_C_BINDING
       integer(c_int) :: compress
       integer(c_int) :: sourceLen
       integer(c_int) :: destLen
       integer(c_char), pointer :: dest, source
     end function compress
  end interface

  interface
     function compressBound(sourceLen) bind(C, name = 'compressBound')
       use, intrinsic :: ISO_C_BINDING
       integer(c_int) :: sourceLen
     end function compressBound
  end interface

  interface
     function deflateInit(stream, level) bind(C, name = 'deflateInit')
       use, intrinsic :: ISO_C_BINDING
       type, bind(C) :: streamtype
          integer(c_int) :: zalloc
          integer(c_int) :: zfree
          integer(c_int) :: opaque

          character(c_char) :: msg(*)

          integer(c_int) :: data_type
          integer(c_int) :: adler
          integer(c_int) :: reserved
       end type 
       integer(c_int) :: deflateInit
       type(STREAMTYPE) :: stream
       integer(c_int) :: level


     end function deflateInit
  end interface

  interface 
     function gzOpen(path,mode) bind(C, NAME='gzopen') 
       use, intrinsic :: ISO_C_BINDING 
       character(c_char) path(*),mode(*) 
       integer(c_int) :: gzOpen 
     end function gzOpen
  end interface
  interface 
     function gzGets(file,buf,len) bind(C, NAME='gzgets') 
       use, intrinsic :: ISO_C_BINDING 
       character(c_char) buf(*) 
       type(c_ptr) :: gzGets 
       integer(c_int),value :: file,len 
     end function gzGets
  end interface
  interface 
     function gzClose(file) bind(C, NAME='gzclose') 
       use, intrinsic :: ISO_C_BINDING 
       integer(c_int) :: gzClose 
       integer(c_int),value :: file 
     end function gzClose
  end interface



  public

contains

  subroutine gzip_test()

    type(c_ptr) :: buf 
    character(len=10) :: mode = "r"//achar(0), file ="test.gz"//achar(0) 
    integer(c_int) :: handle,sz,r 
    character(kind=C_CHAR,len=100) :: buffer 
    handle = gzOpen(file,mode) 
    print *,handle 
    sz=len(buffer) 
    buf=gzGets(handle,buffer,sz) 
    print *,buffer 
    r=gzClose(handle) 
    print *,r 

  end subroutine gzip_test

  subroutine test_compress()
    integer(bigInt) :: iBytes, iSize
    iBytes = 10000000000
    iSize = compressBound(iBytes)
    write(*,*) "max size ",isize
  end subroutine test_compress
#endif


end module zlib_mod
