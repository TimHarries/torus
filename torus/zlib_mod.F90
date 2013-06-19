module zlib_mod

  use kind_mod
  use, intrinsic :: ISO_C_BINDING
  use vector_mod
  implicit none

  integer, parameter :: maxBuffer = 2**16
  integer(kind=bigInt) :: nBuffer
  integer(kind=1) :: buffer(maxBuffer)
  logical :: uncompressedDumpFiles

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

  interface writeCompressedFile
     module procedure writeCompressedFileDoubleValue
     module procedure writeCompressedFileDoubleArray1d
     module procedure writeCompressedFileDoubleArray2d
     module procedure writeCompressedFileDoubleArray3d
     module procedure writeCompressedFileVectorValue
     module procedure writeCompressedFileVectorArray1d
     module procedure writeCompressedFileLogicalArray1d
     module procedure writeCompressedFileLogicalValue
     module procedure writeCompressedFileString
     module procedure writeCompressedFileIntegerValue
     module procedure writeCompressedFileRealValue
     module procedure writeCompressedFileRealArray1d
     module procedure writeCompressedFileRealArray2d
     module procedure writeCompressedFileIntegerArray1d
  end interface

  interface readCompressedFile
     module procedure readCompressedFileDoubleValue
     module procedure readCompressedFileDoubleArray1d
     module procedure readCompressedFileDoubleArray2d
     module procedure readCompressedFileDoubleArray3d
     module procedure readCompressedFileVectorValue
     module procedure readCompressedFileVectorArray1d
     module procedure readCompressedFileLogicalArray1d
     module procedure readCompressedFileLogicalValue
     module procedure readCompressedFileString
     module procedure readCompressedFileIntegerValue
     module procedure readCompressedFileRealValue
     module procedure readCompressedFileRealArray1d
     module procedure readCompressedFileRealArray2d
     module procedure readCompressedFileIntegerArray1d
  end interface


  public

contains

  subroutine compressBytes(inArray, nIn, outArray, nOut)
    integer(kind=1), pointer :: inArray(:)
    integer(bigint) :: nIn, nOut
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

  subroutine uncompressBytes(inArray, nIn, outArray, nOut)
    integer(kind=1), pointer :: inArray(:)
    integer(bigint) :: nIn, nOut
    integer(kind=1), pointer :: outArray(:)
    type(c_ptr) :: c_outlen
    integer :: j
    integer(c_long), target :: outLen
    integer(kind=1), pointer :: inBuffer(:), outBuffer(:)
    type(c_ptr) :: c_inbuffer = c_null_ptr, c_outbuffer = c_null_ptr

    allocate(inBuffer(1:nIn))
    inBuffer(1:nIn) = inArray(1:nIn)

    allocate(outBuffer(1:nOut))

    c_inbuffer = c_loc(inBuffer(1))
    c_outbuffer = c_loc(outBuffer(1))
    outlen = int(nOut)
    c_outlen = c_loc(outlen)


    j = uncompress(c_outbuffer, c_outLen, c_inbuffer, int(nIn, kind=c_long))
    allocate(outArray(1:nOut))
    outArray(1:nOut) = outbuffer(1:nOut)
    deallocate(inBuffer, outBuffer)
  end subroutine uncompressBytes

  subroutine openCompressedFile(lunit, thisFilename, positionStatus)
    integer :: lunit
    character(len=*), optional :: positionStatus
    character(len=*) :: thisFilename
    

    if (PRESENT(positionStatus)) then
       open(lunit, file=thisFilename, form="unformatted", status="unknown",position=positionStatus)
    else
       open(lunit, file=thisFilename, form="unformatted", status="unknown")
       buffer = 0 
       nbuffer = 0
    endif
  end subroutine openCompressedFile

  subroutine zeroZlibBuffer()
    buffer = 0
    nbuffer = 0
  end subroutine zeroZlibBuffer

  subroutine closeCompressedFile(lunit, flushBuffer)
    integer :: lunit
    integer(bigint) :: nOut
    logical, optional :: flushBuffer
    integer(kind=1), pointer :: outBuffer(:), inBuffer(:)
    integer :: i

    if (PRESENT(flushBuffer)) then
       if (flushBuffer) then
          i = int(compressBound(int(maxBuffer, kind=c_long)))
!          allocate(outBuffer(1:i))
          allocate(inBuffer(1:maxBuffer))
          inBuffer = buffer
          call compressBytes(inbuffer, int(maxBuffer, kind=bigInt), outBuffer, nOut)
          write(lunit) nout
          write(lunit) outBuffer(1:nOut)
          deallocate(outBuffer, inBuffer)
          nBuffer = 0
          buffer = 0
       endif
    else
       nBuffer = 0
       buffer = 0
    endif
    close(lunit)

  end subroutine closeCompressedFile

  subroutine compressTest()
    real(double) :: test(100000)
    character(len=4) :: testString
    integer :: i,j

    do i = 1, size(test)
       test(i) = dble(i)
    enddo
    call openCompressedFile(33, "test.zlib")
    call writeCompressedFile(33, "crap")
    call writeCompressedFile(33, 1234)
    call writeCompressedFile(33, "shit")
    call writeCompressedFile(33, 5678)
    call closeCompressedFile(33, flushBuffer=.true.)
    open(34,file="test.unc",status="unknown",form="unformatted")
    write(34) test
    close(34)
    test = 0.d0
    call openCompressedFile(33, "test.zlib")
    call readCompressedFile(33, testString)
    write(*,*) testString
    call readCompressedFile(33, j)
    write(*,*) j
    call readCompressedFile(33, testString)
    write(*,*) testString
    call readCompressedFile(33, j)
    write(*,*) j
    call closeCompressedFile(33)
  end subroutine compressTest
       

  subroutine writeCompressedFileDoubleValue(lunit, doubleValue)
    integer :: lunit
    real(double) :: doubleValue

    call writeCompressedFileGeneric(lunit, doubleValue=doubleValue)
  end subroutine writeCompressedFileDoubleValue


  subroutine writeCompressedFileRealValue(lunit, realValue)
    integer :: lunit
    real :: realValue

    call writeCompressedFileGeneric(lunit, realValue=realValue)
  end subroutine writeCompressedFileRealValue

  subroutine writeCompressedFileString(lunit, cString)
    integer :: lunit
    character(len=*) :: cString

    call writeCompressedFileGeneric(lunit, cString=cString)
  end subroutine writeCompressedFileString

  subroutine writeCompressedFileIntegerValue(lunit, intValue)
    integer :: lunit
    integer :: intValue

    call writeCompressedFileGeneric(lunit, intValue=intValue)
  end subroutine writeCompressedFileIntegerValue

  subroutine writeCompressedFileIntegerArray1d(lunit, integerArray1d)
    integer :: lunit
    integer :: integerArray1d(:)

    call writeCompressedFileGeneric(lunit, integerArray1d=integerArray1d)
  end subroutine writeCompressedFileIntegerArray1d

  subroutine writeCompressedFileDoubleArray1d(lunit, doubleArray1d)
    integer :: lunit
    real(double) :: doubleArray1d(:)

    call writeCompressedFileGeneric(lunit, doubleArray1d=doubleArray1d)
  end subroutine writeCompressedFileDoubleArray1d

  subroutine writeCompressedFileLogicalArray1d(lunit, value)
    integer :: lunit
    logical :: value(:)

    call writeCompressedFileGeneric(lunit, logicalArray1d=value)
  end subroutine writeCompressedFileLogicalArray1d

  subroutine writeCompressedFileLogicalValue(lunit, value)
    integer :: lunit
    logical :: value

    call writeCompressedFileGeneric(lunit, logicalValue=value)
  end subroutine writeCompressedFileLogicalValue

  subroutine writeCompressedFileVectorArray1d(lunit, value)
    integer :: lunit
    type(VECTOR) :: value(:)

    call writeCompressedFileGeneric(lunit, vectorArray1d=value)
  end subroutine writeCompressedFileVectorArray1d

  subroutine writeCompressedFileVectorValue(lunit, value)
    integer :: lunit
    type(VECTOR) :: value

    call writeCompressedFileGeneric(lunit, vectorValue=value)
  end subroutine writeCompressedFileVectorValue

  subroutine writeCompressedFileRealArray1d(lunit, value)
    integer :: lunit
    real :: value(:)

    call writeCompressedFileGeneric(lunit, realArray1d=value)
  end subroutine writeCompressedFileRealArray1d

  subroutine writeCompressedFileRealArray2d(lunit, value)
    integer :: lunit
    real :: value(:,:)

    call writeCompressedFileGeneric(lunit, realArray2d=value)
  end subroutine writeCompressedFileRealArray2d


  subroutine writeCompressedFileDoubleArray2d(lunit, doubleArray2d)
    integer :: lunit
    real(double) :: doubleArray2d(:,:)

    call writeCompressedFileGeneric(lunit, doubleArray2d=doubleArray2d)
  end subroutine writeCompressedFileDoubleArray2d

  subroutine writeCompressedFileDoubleArray3d(lunit, doubleArray3d)
    integer :: lunit
    real(double) :: doubleArray3d(:,:,:)

    call writeCompressedFileGeneric(lunit, doubleArray3d=doubleArray3d)
  end subroutine writeCompressedFileDoubleArray3d



  subroutine writeCompressedFileGeneric(lunit, doubleArray1d, doubleArray2d, doubleArray3d, doubleValue, &
       cString, intValue, integerArray1d, realValue, logicalArray1d, realArray1d, realArray2d, &
       vectorArray1d, logicalValue, vectorValue)
    integer :: lunit
    integer, optional :: intValue
    logical, optional :: logicalValue
    real, optional :: realValue
    type(VECTOR), optional :: vectorValue
    integer, optional :: integerArray1d(:)
    type(VECTOR), optional :: vectorArray1d(:)
    logical, optional :: logicalArray1d(:)
    real, optional :: realArray1d(:)
    character(len=*), optional :: cString
    real(double), optional :: doubleArray1d(:)
    real(double), optional :: doubleArray2d(:,:)
    real, optional :: realArray2d(:,:)
    real(double), optional :: doubleArray3d(:,:,:)
    real(double), optional :: doubleValue
    integer(kind=1), allocatable :: inputBuffer(:)
    integer(kind=1), pointer :: outBuffer(:), inBuffer(:)
    integer :: nBytes, i, iByte
    integer(kind=bigInt) :: nout
    real(double), parameter :: testDouble  = 0.d0
    real, parameter :: testReal = 0.
    integer, parameter :: testInteger = 0
    integer :: nBytesDouble, nBytesReal, nBytesInteger

    nBytesReal = 4
    nBytesDouble = 8
    nBytesInteger = 4
    
#ifdef MEMCHECK
       nbytesReal = int(sizeof(testreal))
       nbytesDouble = int(sizeof(testDouble))
       nbytesInteger = int(sizeof(testInteger))
#endif

    
    if (present(cString)) then
       nbytes = len(cString)
       allocate(inputBuffer(1:nBytes))
       inputBuffer(1:nBytes) = transfer(cString, inputBuffer(1:nBytes))
!       write(*,*) "input buffer ",inputBuffer(1:nBytes)
    endif

    if (present(intValue)) then
       nbytes = nBytesInteger
       allocate(inputBuffer(1:nBytes))
       inputBuffer(1:nBytes) = transfer(intValue, inputBuffer(1:nBytes))
    endif

    if (present(logicalArray1d)) then
#ifdef MEMCHECK
       nbytes = int(SIZE(logicalArray1d)*sizeof(logicalArray1d(1)))
#endif
       allocate(inputBuffer(1:nBytes))
       inputBuffer(1:nBytes) = transfer(logicalArray1d, inputBuffer(1:nBytes))
    endif

    if (present(logicalValue)) then
#ifdef MEMCHECK
       nbytes = int(sizeof(logicalValue))
#endif
       allocate(inputBuffer(1:nBytes))
       inputBuffer(1:nBytes) = transfer(logicalValue, inputBuffer(1:nBytes))
    endif

    if (present(realArray1d)) then
       nbytes = SIZE(realArray1d)*nBytesReal
       allocate(inputBuffer(1:nBytes))
       inputBuffer(1:nBytes) = transfer(realArray1d, inputBuffer(1:nBytes))
    endif

    if (present(vectorArray1d)) then
       nbytes = SIZE(vectorArray1d)*nBytesDouble*3
       allocate(inputBuffer(1:nBytes))
       inputBuffer(1:nBytes) = transfer(vectorArray1d, inputBuffer(1:nBytes))
    endif


    if (present(vectorValue)) then
       nbytes = nBytesDouble*3
       allocate(inputBuffer(1:nBytes))
       inputBuffer(1:nBytes) = transfer(vectorValue, inputBuffer(1:nBytes))
    endif

    if (present(realValue)) then
#ifdef MEMCHECK
       nbytes = int(sizeof(realValue))
#endif
       allocate(inputBuffer(1:nBytes))
       inputBuffer(1:nBytes) = transfer(realValue, inputBuffer(1:nBytes))
    endif


    if (present(integerArray1d)) then
       nbytes = size(integerArray1d)*nBytesInteger
       allocate(inputBuffer(1:nBytes))
       inputBuffer(1:nBytes) = transfer(integerArray1d, inputBuffer(1:nBytes))
    endif


    if (present(doubleArray1d)) then
       nbytes = size(doubleArray1d)*nBytesDouble
       allocate(inputBuffer(1:nBytes))
       inputBuffer(1:nBytes) = transfer(doubleArray1d, inputBuffer(1:nBytes))
    endif

    if (present(doubleArray2d)) then
       nbytes = size(doubleArray2d,1)*size(doubleArray2d,2)*nBytesdouble
       allocate(inputBuffer(1:nBytes))
       inputBuffer(1:nBytes) = transfer(doubleArray2d, inputBuffer(1:nBytes))
    endif

    if (present(realArray2d)) then
       nbytes = size(realArray2d,1)*size(realArray2d,2)*nBytesReal
       allocate(inputBuffer(1:nBytes))
       inputBuffer(1:nBytes) = transfer(realArray2d, inputBuffer(1:nBytes))
    endif

    if (present(doubleArray3d)) then
       nbytes = size(doubleArray3d,1)*size(doubleArray3d,2)*size(doubleArray3d,3)*nBytesDouble
       allocate(inputBuffer(1:nBytes))
       inputBuffer(1:nBytes) = transfer(doubleArray3d, inputBuffer(1:nBytes))
    endif

    if (present(doubleValue)) then
       nbytes = nBytesDouble
       allocate(inputBuffer(1:nBytes))
       inputBuffer(1:nBytes) = transfer(doubleValue, inputBuffer(1:nBytes))
    endif

    if ((nBuffer+nBytes) <= maxBuffer) then
       buffer(nBuffer+1:(nBuffer+nBytes)) = inputBuffer(1:nBytes)
       nBuffer = nBuffer + nBytes

!       write(*,*) "nbuffer ",nBuffer
!       write(*,*) "buffer ",buffer(1:nBuffer)
       if (nBuffer == maxBuffer) then
          i = int(compressBound(int(maxBuffer, kind=c_long)))
!          allocate(outBuffer(1:i))
          allocate(inBuffer(1:maxBuffer))
          inBuffer = buffer
          call compressBytes(inbuffer, int(maxBuffer, kind=bigInt), outBuffer, nOut)
          write(lunit) nout
          write(lunit) outBuffer(1:nOut)
          deallocate(outBuffer, inBuffer)
          nBuffer = 0
          buffer = 0
       endif
       goto 666
    endif
    
    do iByte = 1, nBytes
       nBuffer = nBuffer + 1
       buffer(nBuffer) = inputBuffer(iByte)
       if (nBuffer == maxBuffer) then
          i = int(compressBound(int(maxBuffer, kind=c_long)))
!          allocate(outBuffer(1:i))
          allocate(inBuffer(1:maxBuffer))
          inBuffer = buffer
          call compressBytes(inBuffer, int(maxBuffer, kind=bigInt), outBuffer, nOut)
          write(lunit) nout
          write(lunit) outBuffer(1:nOut)
          deallocate(outBuffer, inBuffer)
          nBuffer = 0
          buffer = 0
       endif
    enddo

666 continue
    if (allocated(inputBuffer)) deallocate(inputbuffer)
  end subroutine writeCompressedFileGeneric

  subroutine readCompressedFileDoubleValue(lunit, doubleValue)
    integer :: lunit
    real(double) :: doubleValue

    call readCompressedFileGeneric(lunit, doubleValue=doubleValue)
  end subroutine readCompressedFileDoubleValue

  subroutine readCompressedFileVectorValue(lunit, value)
    integer :: lunit
    type(VECTOR) :: value

    call readCompressedFileGeneric(lunit, vectorValue=value)
  end subroutine readCompressedFileVectorValue

  subroutine readCompressedFileLogicalValue(lunit, value)
    integer :: lunit
    logical :: value

    call readCompressedFileGeneric(lunit, logicalValue=value)
  end subroutine readCompressedFileLogicalValue

  subroutine readCompressedFileIntegerValue(lunit, value)
    integer :: lunit
    integer :: value

    call readCompressedFileGeneric(lunit, integerValue=value)
  end subroutine readCompressedFileIntegerValue

  subroutine readCompressedFileRealValue(lunit, value)
    integer :: lunit
    real :: value

    call readCompressedFileGeneric(lunit, realValue=value)
  end subroutine readCompressedFileRealValue

  subroutine readCompressedFileRealArray1d(lunit, value)
    integer :: lunit
    real :: value(:)

    call readCompressedFileGeneric(lunit, realArray1d=value)
  end subroutine readCompressedFileRealArray1d

  subroutine readCompressedFileIntegerArray1d(lunit, value)
    integer :: lunit
    integer :: value(:)

    call readCompressedFileGeneric(lunit, integerArray1d=value)
  end subroutine readCompressedFileIntegerArray1d

  subroutine readCompressedFileRealArray2d(lunit, value)
    integer :: lunit
    real :: value(:,:)

    call readCompressedFileGeneric(lunit, realArray2d=value)
  end subroutine readCompressedFileRealArray2d

  subroutine readCompressedFileString(lunit, value)
    integer :: lunit
    character(len=*) :: value

    call readCompressedFileGeneric(lunit, cString=value)
  end subroutine readCompressedFileString

  subroutine readCompressedFileLogicalArray1d(lunit, value)
    integer :: lunit
    logical :: value(:)

    call readCompressedFileGeneric(lunit, logicalArray1d=value)
  end subroutine readCompressedFileLogicalArray1d

  subroutine readCompressedFileVectorArray1d(lunit, value)
    integer :: lunit
    type(VECTOR) :: value(:)

    call readCompressedFileGeneric(lunit, vectorArray1d=value)
  end subroutine readCompressedFileVectorArray1d

  subroutine readCompressedFileDoubleArray1d(lunit, doubleArray1d)
    integer :: lunit
    real(double) :: doubleArray1d(:)

    call readCompressedFileGeneric(lunit, doubleArray1d=doubleArray1d)
  end subroutine readCompressedFileDoubleArray1d

  subroutine readCompressedFileDoubleArray2d(lunit, doubleArray2d)
    integer :: lunit
    real(double) :: doubleArray2d(:,:)

    call readCompressedFileGeneric(lunit, doubleArray2d=doubleArray2d)
  end subroutine readCompressedFileDoubleArray2d


  subroutine readCompressedFileDoubleArray3d(lunit, doubleArray3d)
    integer :: lunit
    real(double) :: doubleArray3d(:,:,:)

    call readCompressedFileGeneric(lunit, doubleArray3d=doubleArray3d)
  end subroutine readCompressedFileDoubleArray3d


  subroutine readCompressedFileGeneric(lunit, doubleValue, doubleArray1d, doubleArray2d, doubleArray3d, &
       vectorValue, logicalValue, integerValue, realValue, realArray1d, integerArray1d, realArray2d, cString, &
       logicalArray1d, vectorArray1d)
    integer :: lunit
    real, optional :: realValue, realArray1d(:), realArray2d(:,:)
    logical, optional :: logicalValue, logicalArray1d(:)
    type(VECTOR), optional :: vectorvalue, vectorArray1d(:)
    integer, optional :: integerValue, integerArray1d(:)
    character(len=*), optional :: cString
    real(double), optional :: doubleValue
    real(double), optional :: doubleArray1d(:)
    real(double), optional :: doubleArray2d(:,:)
    real(double), optional :: doubleArray3d(:,:,:)
    integer(kind=1), allocatable :: inputBuffer(:)
    integer(kind=1), pointer :: outBuffer(:), compressedBuffer(:)
    integer :: nBytes, iByte
    integer(kind=bigInt) :: nout, nCompressed
    
    if (present(cString)) then
       nbytes = len(cString)
       allocate(inputBuffer(1:nBytes))
    endif



    if (present(realValue)) then
       nbytes = 4
       allocate(inputBuffer(1:nBytes))
    endif

    if (present(doubleValue)) then
       nbytes = 8
       allocate(inputBuffer(1:nBytes))
    endif

    if (present(doubleArray1d)) then
       nbytes = size(doubleArray1d)*8
       allocate(inputBuffer(1:nBytes))
    endif

    if (present(realArray1d)) then
       nbytes = size(realArray1d)*4
       allocate(inputBuffer(1:nBytes))
    endif

    if (present(doubleArray2d)) then
       nbytes = size(doubleArray2d,1)*size(doubleArray2d,2)*8
       allocate(inputBuffer(1:nBytes))
    endif

    if (present(realArray2d)) then
       nbytes = size(realArray2d,1)*size(realArray2d,2)*4
       allocate(inputBuffer(1:nBytes))
    endif

    if (present(doubleArray3d)) then
       nbytes = size(doubleArray3d,1)*size(doubleArray3d,2)*size(doubleArray3d,3)*8
       allocate(inputBuffer(1:nBytes))
    endif

    if (present(vectorValue)) then
       nBytes = 3 * 8
       allocate(inputBuffer(1:nBytes))
    endif

    if (present(vectorArray1d)) then
       nBytes = size(vectorArray1d) * 3 * 8
       allocate(inputBuffer(1:nBytes))
    endif

    if (present(logicalValue)) then
#ifdef MEMCHECK
       nBytes = int(sizeof(logicalValue))
#endif
       allocate(inputBuffer(1:nBytes))
    endif

    if (present(integerValue)) then
       nBytes = 4
       allocate(inputBuffer(1:nBytes))
    endif

    if (present(integerArray1d)) then
       nBytes = size(integerArray1d)*4
       allocate(inputBuffer(1:nBytes))
    endif

    if (present(logicalArray1d)) then
#ifdef MEMCHECK
       nBytes = int(size(logicalArray1d)*sizeof(logicalArray1d(1)))
#endif
       allocate(inputBuffer(1:nBytes))
    endif
       


    do iByte = 1, nBytes

       if (nBuffer == 0) then
          read(lunit) nCompressed
          allocate(compressedBuffer(1:nCompressed))
          read(lunit) compressedBuffer
          nOut = maxBuffer
          call uncompressBytes(compressedBuffer, nCompressed, outBuffer, nOut)
          if (nOut /= maxBuffer) then
             write(*,*) "Uncompressed buffer size is ", nout
             stop
          endif
          buffer =  outbuffer
          deallocate(outBuffer, compressedBuffer)
          nBuffer = maxBuffer
       endif

       inputBuffer(iByte) = buffer(maxBuffer - nBuffer + 1)
       nBuffer = nBuffer - 1
    enddo

    if (present(doubleValue)) then
       doubleValue = transfer(inputBuffer(1:nBytes), doubleValue)
    endif

    if (present(realValue)) then
       realValue = transfer(inputBuffer(1:nBytes), realValue)
    endif

    if (present(logicalValue)) then
       logicalValue = transfer(inputBuffer(1:nBytes), logicalValue)
    endif

    if (present(integerValue)) then
       integerValue = transfer(inputBuffer(1:nBytes), integerValue)
    endif

    if (present(cString)) then
       cString = transfer(inputBuffer(1:nBytes), cString)
    endif

    if (present(vectorValue)) then
       vectorValue = transfer(inputBuffer(1:nBytes), vectorValue)
    endif


    if (present(doubleArray1d)) then
       doubleArray1d = transfer(inputBuffer(1:nBytes), doubleArray1d)
    endif

    if (present(integerArray1d)) then
       integerArray1d = transfer(inputBuffer(1:nBytes), integerArray1d)
    endif

    if (present(logicalArray1d)) then
       logicalArray1d = transfer(inputBuffer(1:nBytes), logicalArray1d)
    endif

    if (present(realArray1d)) then
       realArray1d = transfer(inputBuffer(1:nBytes), realArray1d)
    endif

    if (present(vectorArray1d)) then
       vectorArray1d = transfer(inputBuffer(1:nBytes), vectorArray1d)
    endif

    if (present(doubleArray1d)) then
       doubleArray1d = transfer(inputBuffer(1:nBytes), doubleArray1d)
    endif

    if (present(doubleArray2d)) then
       doubleArray2d = RESHAPE(transfer(inputBuffer(1:nBytes), doubleArray2d), SHAPE(doubleArray2d))
    endif

    if (present(realArray2d)) then
       realArray2d = RESHAPE(transfer(inputBuffer(1:nBytes), realArray2d), SHAPE(realArray2d))
    endif

    if (present(doubleArray3d)) then
       doubleArray3d = RESHAPE(transfer(inputBuffer(1:nBytes), doubleArray3d), SHAPE(doubleArray3d))
    endif


    deallocate(inputBuffer)
  end subroutine readCompressedFileGeneric

#endif


end module zlib_mod
