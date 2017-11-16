module citations_mod
  use unix_mod
  use messages_mod
  implicit none

  character(len=100) :: bibcodeArray(101)
  character(len=100) :: bibcodeComment(100)
  integer :: nBibcodeArray
  

  character(len=100) :: bibcodeDatabase(100)
  character(len=1000) :: bibtexDatabase(100)
  integer :: nBibcodeDatabase
  logical :: databaseUpdated

contains

  subroutine readBibtexDatabase()
    integer :: i 
    character(len=80) :: datadirectory, thisFile
    call unixGetenv("TORUS_DATA", dataDirectory, i)
    thisFile = trim(dataDirectory)//"/bibtex.dat"
    open(33, file=thisFile, status="old",form="formatted")
    read(33,*) nBibcodeDatabase 
    do i = 1, nBibcodeDatabase
       read(33, '(a100)') bibcodeDatabase(i)
       read(33, '(a1000)') bibtexDatabase(i)
    enddo
    close(33)
  end subroutine readBibtexDatabase

  subroutine writeBibtexDatabase()
    integer :: i 
    character(len=80) :: dataDirectory, thisFile
    call unixGetenv("TORUS_DATA", dataDirectory, i)
    thisFile = trim(dataDirectory)//"/bibtex.dat"
    open(33, file=thisFile, status="old",form="formatted")
    write(33,*, err=20) nBibcodedatabase
    do i = 1, nBibcodeDatabase
       write(33, '(a100)') bibcodeDatabase(i)
       write(33, '(a1000)') bibtexDatabase(i)
    enddo
20 continue
    close(33)
  end subroutine writeBibtexDatabase

  subroutine setUpBibcodesOnParametersFile()
    use inputs_mod

    if (molecularPhysics) call addBibcode("2010MNRAS.407..986R", "Molecular line transport")
    if (photoionPhysics) call addBibcode("2012MNRAS.420..562H", "Photoionisation loop")
    if (radiativeEquilibrium) then
       call addBibcode("2004MNRAS.350..565H", "Radiative equilibrium in TORUS")
       call addBibcode("1999A&A...344..282L", "Leon Lucy's radiative equilibrium algorithm")
    endif

  end subroutine setUpBibcodesOnParametersFile

  logical function inDataBase(bibcode)
    character(len=*) :: bibcode
    integer :: i
    inDatabase = .false.
    do i = 1, nBibcodedatabase
       if (trim(bibcode)==trim(bibcodedatabase(i))) then
          inDatabase = .true.
          exit
       endif
    enddo
  end function inDataBase

  function getBibtexFromDataBase(bibcode)
    character(len=*) :: bibcode
    character(len=1000) :: getBibTexFromDatabase
    integer :: i
    do i = 1, nBibcodedatabase
       if (trim(bibcode)==trim(bibcodedatabase(i))) then
          getBibtexFromDatabase = bibtexDataBase(i)
          exit
       endif
    enddo
  end function getBibtexFromDataBase

  subroutine initBibcode()
    nBibCodeArray = 0
    call readBibtexDatabase()
    databaseUpdated = .false.
  end subroutine initBibcode

  subroutine addBibcode(bibcode, comment)
    character(len=*) :: bibcode, comment

    if (uniqueCode(bibcode)) then
       nBibcodeArray = nBibCodeArray + 1
       bibCodeArray(nBibCodeArray)  = bibCode
       bibCodeComment(nBibCodeArray)  = comment
    endif
  end subroutine addBibcode

  logical function uniqueCode(bibcode)
    character(len=*) :: bibCode
    integer :: i

    uniqueCode = .true.

    if (nBibcodeArray > 0) then
       do i = 1, nBibCodeArray
          if (trim(bibcode) == trim(bibCodeArray(i))) then
             uniqueCode = .false.
             exit
          endif
       enddo
    endif
  end function uniqueCode


  subroutine writeBibtex(thisFile)
    character(len=*) :: thisFile
    character(len=1000) :: bibtexString
    integer :: i
    logical :: ok
    open(33, file=thisfile, status="unknown", form="formatted")
    do i = 1, nBibCodeArray
       call getBibtex(bibCodeArray(i), bibtexString, ok)
       write(33,'(a)') "%"
       write(33,'(a,a)') "% ", trim(bibCodeComment(i))
       write(33,'(a)') "%"

       if (ok) then
          call addnewLines(bibtexString)
          write(33,'(a)') trim(bibtexString)
       else
          write(33,'(a,a)') "Cannot search ADS: You need to cite ",trim(bibCodeArray(i))
       endif
    enddo
    close(33)
    call writeBanner("Please remember to cite the papers in the torus.bib file","!")
    if (databaseUpdated) call writeBibtexDatabase()
  end subroutine writeBibtex

  

  subroutine getBibtex(bibcode, returnedString, ok) 
    character(len=*) :: bibcode, returnedString
    character(len=1000) :: curlString, tempString
    integer :: i
    logical :: ok

    ok = .true.

    if (inDatabase(bibcode)) then
       returnedString = getBibtexFromDatabase(bibcode)
    else

       write(curlstring,'(a,a,a,a)') 'curl -sH "Authorization: Bearer CICjAcmSaIcj7jMssIBx41Eq63UuZtCqQOMtWPf0" -H', &
            ' "Content-Type: text/plain" -X POST -d ''{"bibcode":["',trim(bibcode), &
            '"]}'' https://api.adsabs.harvard.edu/v1/export/bibtex > adsreturn'
!       call execute_command_line(curlString, wait=.true., exitstat=i)
       i = 99
       if (i /= 0) then
          call writeWarning("Cannot execute ADS search")
          ok = .false.
          goto 666
       endif
       open(34, file="adsreturn", form="formatted",status="old")
       read(34,'(a1000)') returnedString
       close(34)
!       call execute_command_line("rm adsreturn", wait=.true.)
       call parseStringToBibtex(returnedString, tempString)
       returnedString = tempString
       nBibcodeDatabase = nBibcodeDatabase + 1
       bibcodeDatabase(nBibcodeDatabase) = bibcode
       bibtexDatabase(nBibcodeDatabase) = returnedString
       databaseUpdated = .true.
    endif
666 continue
  end subroutine getBibtex

  subroutine parseStringToBibtex(longstring, bibtexString)
    character(len=*) :: longString, bibtexString
    integer :: iStart, iEnd
    iStart = index(longString, "@ARTICLE")
    iEnd = index(longString,"\n\n")-1
    bibtexString = longString(istart:iend)
  end subroutine parseStringToBibtex


  subroutine addNewlines(bibtexString)
    character(len=*) :: bibtexString
    character(len=1) :: newline
    integer :: i
    newline = new_line(newline)
    do i = 1, len(bibtexString)-1
       if (bibtexString(i:(i+1))=="\n") then
          bibtexString(i:(i+1)) = " "//newline
       endif
    enddo
  end subroutine addNewlines


end module citations_mod
