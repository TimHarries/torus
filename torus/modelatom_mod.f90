module modelatom_mod

  ! module for model atoms created by tjh 23/10/06

  use kind_mod
  use constants_mod
  use unix_mod
  use messages_mod
  use hyd_col_coeff
  use utils_mod
  use stateq_mod
  use utils_mod

  type MODELATOM
     character(len=10) :: name
     integer :: charge
     real(double) :: mass ! amu
     integer :: nLevels
     character(len=10), pointer :: level(:)
     character(len=10) :: continuumLevel
     real(double), pointer :: energy(:)  ! erg
     real(double), pointer :: g(:)
     integer, pointer :: ionStage(:)
     integer, pointer :: nQuantum(:)
     integer :: nTrans
     integer :: nRBBtrans
     integer, pointer :: indexRBB(:)
     character(len=3), pointer :: transType(:) ! RBB, RBF, CBB, CBF etc
     real(double), pointer :: transFreq(:)
     integer, pointer :: iLower(:)
     integer, pointer :: iUpper(:)
     integer, pointer :: equation(:)
     integer, pointer :: nParams(:)
     real(double), pointer :: params(:,:)
  end type MODELATOM


contains

  subroutine readAtom(thisAtom, atomfilename)
    character(len=*) :: atomfilename
    character(len=200) :: dataDirectory, thisfilename
    type(MODELATOM) :: thisAtom
    integer :: i, nChunks, j
    character(len=120) :: junk
    character(len=20) :: chunk(20)
    character(len=80) :: message
    real(double) :: a, blu, bul
    call unixGetenv("TORUS_DATA", dataDirectory, i)
    thisfilename = trim(dataDirectory)//"/"//atomfilename

    open(30, file=thisfilename, status="old", form="formatted")

    read(30,*) thisAtom%name
    write(message,'(a,a)') "Reading atom model from: ",trim(atomfilename)
    call writeInfo(message,TRIVIAL)
    read(30,*) thisAtom%charge
    read(30,*) thisAtom%mass
    read(30,*) thisAtom%nLevels


    allocate(thisAtom%energy(1:thisAtom%nLevels))
    allocate(thisAtom%g(1:thisAtom%nLevels))
    allocate(thisAtom%level(1:thisAtom%nLevels))

    do i = 1, thisAtom%nLevels
       read(30,'(a120)') junk
       call splitIntoChunks(junk, nChunks, chunk)
       thisAtom%level(i) = chunk(1)
       thisAtom%continuumLevel = chunk(2)
       read(chunk(3), *) thisAtom%energy(i)
       read(chunk(4), *) thisAtom%g(i)
    enddo

    thisAtom%energy = hydE0eVdb - (thisAtom%energy*hCgs*ergtoev) ! eV

    read(30,*) thisAtom%nTrans


    allocate(thisAtom%transType(1:thisAtom%nTrans))
    allocate(thisAtom%transFreq(1:thisAtom%nTrans))
    allocate(thisAtom%iLower(1:thisAtom%nTrans))
    allocate(thisAtom%iUpper(1:thisAtom%nTrans))
    allocate(thisAtom%equation(1:thisAtom%nTrans))
    allocate(thisAtom%nparams(1:thisAtom%nTrans))
    allocate(thisAtom%params(1:thisAtom%nTrans,10))


    do i = 1, thisAtom%nTrans
       read(30,'(a120)') junk
       call splitIntoChunks(junk, nChunks, chunk)
       thisAtom%transType(i) = chunk(1)
       thisAtom%iLower(i) = getLevel(thisAtom, chunk(2))
       thisAtom%iUpper(i) = getLevel(thisAtom, chunk(3))
       read(chunk(4),*) thisAtom%equation(i)
       read(chunk(5),*) thisAtom%nParams(i)
       do j = 1, thisAtom%nParams(i)
          read(chunk(5+j),*) thisAtom%params(i, j)
       enddo
       if (thisAtom%iUpper(i) <= thisAtom%nLevels) then
          thisAtom%transFreq(i) = ((thisAtom%energy(thisAtom%iUpper(i)) - thisAtom%energy(thisAtom%iLower(i)))*evToErg)/Hcgs
       endif
    enddo
    close(30)
    call writeInfo("Done.",TRIVIAL)

    do i = 1, thisAtom%nTrans
       if (thisAtom%transType(i) == "RBB") then
          call returnEinsteinCoeffs(thisAtom, i, a, Bul, Blu)
          write(*,*) i,1.d8*cSpeed/thisAtom%transFreq(i),a, blu, bul
       endif
    enddo

    thisAtom%nRBBtrans = 0
    do i = 1, thisAtom%nTrans
       if (thisAtom%transType(i) == "RBB") thisAtom%nRBBtrans = thisAtom%nRBBtrans + 1
    enddo
    write(*,*) "Number of radiative BB transitions",thisAtom%nRBBtrans

    allocate(thisAtom%indexRBB(1:thisAtom%nRBBtrans))

    j = 0
    do i = 1, thisAtom%nTrans
       if (thisAtom%transType(i) == "RBB") then
           j = j + 1
           thisAtom%indexRBB(j) = i
        endif
    enddo
       


  end subroutine readAtom







  subroutine splitIntoChunks(cLine, nChunks, chunk)
    character(len=*) :: cLine
    integer :: nChunks
    character(len=20) :: chunk(:)
    integer :: i
    character(len=200) :: tLine

    tLine = cLine
    i = 99
    nChunks = 0
    call stripLeadingSpaces(tLine)

    do while (i /= 0)

       i = index(tLine, " ")
       if (i /= 0) then
          nChunks = nChunks + 1
          chunk(nChunks) = tLine(1:i-1)
          tLine(1:) = tLine(i+1:) 
          call stripLeadingSpaces(tLine)
          if (len(trim(tLine)) == 0) i = 0
       endif
    enddo
  end subroutine splitIntoChunks

  subroutine stripLeadingSpaces(cString)
    character(len=*) :: cString
    logical :: finished

    finished = .false.
    do while (.not.finished)
       if (cString(1:1) == " ") then
          cString(1:) = cString(2:)
       else
          finished = .true.
       endif
       if (len(trim(cString)) == 0) finished = .true.
    enddo
  end subroutine stripLeadingSpaces

  function getLevel(thisAtom, clevel) result(level)
    type(modelAtom) :: thisAtom
    character(len=*) :: cLevel
    integer :: level
    integer :: i
    level = 0
    do i = 1, thisAtom%nLevels
       if (trim(cLevel) == trim(thisAtom%level(i))) then
          level = i
          exit
       endif
    enddo
    if (trim(cLevel) == trim(thisAtom%continuumLevel)) then
       level = thisAtom%nLevels + 1
    endif

    if (level == 0) then
       write(*,*) "Level not found: ",trim(cLevel)
       stop
    endif
  end function getLevel



  subroutine returnEinsteinCoeffs(thisAtom, iTrans, aEinstein, BulEinstein, BluEinstein)
    type(MODELATOM) :: thisAtom
    integer :: iTrans
    real(double) :: aEinstein, BulEinstein, BluEinstein, f

    if (iTrans > thisAtom%nTrans) then
       call writeFatal("returnEinsteinCoeffs: iTrans greater than number of transitions")
       stop
    endif
    if (thisAtom%transType(iTrans) /= "RBB") then
       call writeFatal("returnEinsteinCoeffs: iTrans is not a radiative bound-bound transition")
       stop
    endif
    select case(thisAtom%equation(iTrans))
    case(4,1,3)
       f = thisAtom%params(iTrans,1)
       iUpper = thisAtom%iUpper(iTrans)
       iLower = thisAtom%iLower(iTrans)
       aEinstein = (8.d0*thisAtom%transFreq(iTrans)**2 * pi**2 * eCharge**2) / &
            (mElectron*cSpeed**3) * f * thisAtom%g(iLower) / thisAtom%g(iUpper)
       BluEinstein = (fourPi * pi * eCharge**2)/(mElectron * hCgs * thisAtom%transFreq(iTrans) * cspeed) * f
       BulEinstein = BluEinstein * thisAtom%g(iLower) / thisAtom%g(iUpper)
    case DEFAULT
       call writeFatal("returnEinsteinCoeffs: unrecognised equation for RBB transition")
       stop
    end select
  end subroutine returnEinsteinCoeffs


  function photoCrossSection(thisAtom, iLevel, nu)
    type(MODELATOM) :: thisAtom
    integer :: iLevel
    real(double) :: nu
    
    if (iLevel > thisAtom%nLevels) then
       call writeFatal("photocrosssection: Level greater than nlevels")
       stop
    endif
    select case(thisAtom%name)
       case("H")
          photoCrossSection = annu_hyd(iLevel, nu)
       case DEFAULT
          call writeFatal("photocrosssection: atom not recognised")
          stop
     end select
   end function photoCrossSection


  function collisionRate(thisAtom, iTrans, temperature) result(rate)
    type(MODELATOM) :: thisAtom
    integer :: iTrans
    real(double) :: temperature, ne, rate, u0, u1, u2
    real(double) :: logGamma
    real(double) :: sigma0

    if (iTrans > thisAtom%nTrans) then
       call writeFatal("returnEinsteinCoeffs: iTrans greater than number of transitions")
       stop
    endif
    if (thisAtom%transType(iTrans)(1:1) /= "C") then
       call writeFatal("collisionRate: iTrans is not a collisional transition")
       stop
    endif


    select case(thisAtom%transType(iTrans))
    case("CBB")
       select case(thisAtom%name)
       case("H") 
          rate = cijt_hyd_hillier(thisAtom%iLower(iTrans), thisAtom%iUpper(iTrans), temperature, thisAtom)

       case("He") 
          select case(thisAtom%equation(itrans))
             case(2)
                if (thisAtom%nParams(iTrans) == 1) then
                   u0 = Hcgs*thisAtom%transFreq(iTrans)/(kErg*temperature)
                   rate = 8.631d-6/(thisAtom%g(thisAtom%iLower(itrans))*sqrt(temperature)) * exp(-u0)
                else
                   call writeFatal("equn not implemented")
                   stop
                endif
             case(5)
                EH = hydE0eVdb * evtoErg
                u0 = Hcgs*thisAtom%transFreq(iTrans)/(kErg*temperature)
                rate = 5.465d-11 * sqrt(temperature)*thisAtom%params(itrans,1)
                rate = rate * (EH/(hCgs*thisAtom%transFreq(iTrans)))**2 * u0 * expint(1, u0)
             case(6)
                EH = hydE0eVdb * evtoErg
                u0 = Hcgs*thisAtom%transFreq(iTrans)/(kErg*temperature)
                u1 = u0 + 0.2d0
                rate = 5.465d-11 * sqrt(temperature)*thisAtom%params(itrans,1)
                rate = rate * (EH/(hCgs*thisAtom%transFreq(iTrans)))**2 * u0 * (expint(1, u0) - &
                     (u0/u1)*exp(-0.2)*expint(1,u1))
            case(8)
               logGamma = thisAtom%params(itrans,1) + thisAtom%params(iTrans,2) * &
                    log10(temperature)+thisAtom%params(itrans,3)/(log10(temperature)**2)
               rate = 5.465d-11 * sqrt(temperature) * exp(-(hCgs*thisAtom%transFreq(iTrans)/(kerg*temperature))) * logGamma
             case(9)
                u0 = Hcgs*thisAtom%transFreq(iTrans)/(kErg*temperature)
                rate = 5.465d-11 * sqrt(temperature) * exp(-u0)*(1.d0+u0)
             case DEFAULT
                call writeFatal("CBB rate equation not implemented for He")
                stop
           end select
               

       case DEFAULT
          call writeFatal("collisionRate: bound-bound collision type not implemented")
          stop
       end select
    case("CBF")
       select case(thisAtom%name)
       case("H") 
          rate = cikt_hyd_hillier(thisAtom%iLower(iTrans), temperature, thisAtom)
       case("He")
          if (thisAtom%iLower(iTrans) <= 15) then
             sigma0 = 1.64d0
             u0 = Hcgs*thisAtom%transFreq(iTrans)/(kErg*temperature)
             u1 = u0 + 0.27d0
             u2 = u0 + 1.43d0
             rate =  5.465d-11 * sqrt(temperature) *sigma0 * (u0*expint(1,u0) - (0.728d0*u0**2/u1)*expint(1,u1) - 0.189d0 * u0**2 &
                  * exp(-u0)*((2.d0+u2)/u2*3))
          else
             write(*,*) "????"
             sigma0 = 666.d0!!!!
             u0 = Hcgs*thisAtom%transFreq(iTrans)/(kErg*temperature)
             rate = 1.55d13*thisAtom%params(itrans,2)*sigma0*exp(-u0)/u0 /sqrt(temperature)
          endif
       case DEFAULT
          call writeFatal("collisionRate: bound-bound collision type not implemented2")
          stop
       end select
    case DEFAULT
       call writeFatal("collisionRate: collision type not implemented")
       stop
    end select
  end function collisionRate

  function cijt_hyd_hillier(i,j,t,thisAtom) result (cijt)
    implicit none
    !
    ! this function calculates the excitation rate using values from Hillier !
    type(MODELATOM) :: thisAtom
    integer, intent(in)              :: i,j  ! lower/upper level
    real(double),intent(in) :: t    ! temperature
    real(double) :: cijt
    real(double) :: t1              ! temperature
    real(double) :: chi             ! excitation energy
    real(double) :: bsum(9)         ! summation
    integer               :: i1,j1           ! lower/upper level
    integer :: lower, upper
    real(double) :: factor
    real(double),parameter :: avals(7) =                                    &    
         (/ 2.579997e-10_db, -1.629166e-10_db, 7.713069e-11_db, -2.668768e-11_db, &
         6.642513e-12_db, -9.422885e-13_db, 0.e0_db                            /)
    integer :: level
    real(double), parameter :: eTrans(23) =                      &
         (/ (hydE0eVdb*(1.0d0 - 1.0d0/level**2),level=1,SIZE(eTrans)) /) 
    integer :: r                             ! counter


    t1 = min(t,1.e5_db)
    i1 = min(i,j)
    j1 = max(i,j)
    chi=abs(eTrans(j)-eTrans(i))
    chi=chi/(kev*t1)

    call locate(tempTable,size(tempTable),t,lower)
    if (lower == 0 .or. lower == size(tempTable)) then
       print *, 'In cijt, temperature is out of range! (',t1,')'
       stop
    end if
    upper = lower + 1
    factor = (t1 - tempTable(lower)) / (tempTable(upper) - tempTable(lower))
    cijt = (3.23_db*8.63e-08_db) * (omegaij(lower,j1,i1)*(1.0_db-factor) + omegaij(upper,j1,i1)*(factor)) * &
         exp(-chi) / thisAtom%g(i1) / sqrt(t1*10e-4_db)


    if (j < i) then
       cijt = cijt / (exp(-chi) * 2.)
    endif

  end function cijt_hyd_hillier

  function cikt_hyd_hillier(i,t,thisAtom) result(cikt)
    !
    ! this function calculates the collisional ionization rate
    ! (using Hillier coefficients for most levels).
    ! 
    type(MODELATOM) :: thisAtom
    integer, intent(in)              :: i        ! the level
    real(double),intent(in) :: t        ! the temperature
    real(double) cikt
    real(double) :: t1                  
    real(double) :: gammait             ! see k&c
    real(double) :: lgt
    real(double) :: chi                 ! potential from i to k
    integer :: lower, upper
    real(double) :: factor
    real(double), parameter :: eTrans(23) =                      &
         (/ (hydE0eVdb*(1.0d0 - 1.0d0/level**2),level=1,SIZE(eTrans)) /) 

    ! making cint a PARAMETER may cause problems with XL Fortran
    !  real(double) :: cint(5,10)
    !  cint = reshape(source=                                                                   &
    real(double), parameter :: cint(5,10) = reshape(source=                         &
         (/-0.4350000e0_db, 0.3000000e00_db,0.0000000e0_db, 0.0000000e0_db,0.00000000e0_db,  &
         1.9987261e1_db,-5.8906298e-5_db,0.0000000e0_db,-2.8185937e4_db,5.44441600e7_db,  &
         1.3935312e3_db,-1.6805859e02_db,0.0000000e0_db,-2.5390000e3_db,0.00000000e0_db,  &
         2.0684609e3_db,-3.3415820e02_db,0.0000000e0_db, 0.0000000e0_db,-7.6440625e3_db,  &
         3.2174844e3_db,-5.5882422e02_db,0.0000000e0_db, 0.0000000e0_db,-6.8632500e3_db,  &
         5.7591250e3_db,-1.5163125e03_db,8.1750000e1_db, 0.0000000e0_db,0.00000000e0_db,  &
         1.4614750e4_db,-4.8283750e03_db,3.9335938e2_db, 0.0000000e0_db,0.00000000e0_db,  &
         2.8279250e4_db,-1.0172750e04_db,9.1967968e2_db, 0.0000000e0_db,0.00000000e0_db,  &
         4.6799250e4_db,-1.7627500e04_db,1.6742031e3_db, 0.0000000e0_db,0.00000000e0_db,  &
         -7.4073000e4_db, 6.8599375e03_db,0.0000000e0_db, 2.0169800e5_db,0.00000000e0_db/),&
         shape=(/5,10/))

    t1 = min(t,1.5e5_db)
    lgt=log10(t1)
    chi=hydE0eV-eTrans(i)
    !    if (i .ne. 2) then 
    call locate(tempTable,size(tempTable),t,lower)
    if (lower == 0 .or. lower == size(tempTable)) then
       print *, 'In cikt, temperature is out of range! (',t1,')'
       stop
    end if
    upper = lower + 1
    factor = (t1 - tempTable(lower)) / (tempTable(upper) - tempTable(lower))
    cikt = 3.23_db*8.63e-08_db * (omegaik(lower,i)*(1.0_db-factor) + omegaik(upper,i)*factor) * &
         exp(-chi/(kev*t1)) / thisAtom%g(i) / sqrt(t1*10e-4_db)

    !    else
    !       gammait=cint(1,i)+cint(2,i)*t1+(cint(4,i)/t1)+(cint(5,i)/t1**2)
    !       cikt=(5.465e-11)*sqrt(t1)*exp(-chi/(kev*t1))*gammait
    !    endif

  end function cikt_hyd_hillier

  real(double) pure function annu_hyd(n,nu)
    !
    ! this subroutine calculates the photoionization x-section
    ! for hydrogen from the n-th level for a given freq photon.
    !
    integer, intent(in)              :: n     ! the level
    real(double), intent(in):: nu    ! the photon frequency
    real(double)            :: lam,e ! the photon wavelength
    real(double), parameter :: eTrans(23) =                      &
         (/ (hydE0eVdb*(1.0d0 - 1.0d0/level**2),level=1,SIZE(eTrans)) /) 

    lam=cSpeed/nu
    lam=lam * 1.e8_db

    e = hCgs * nu * ergToEv

    if (e > (hydE0eV - eTrans(n))) then
       annu_hyd=1.044e-26_db*gii(n,1.e0_db,lam)*(lam**3)/dble(n)**5
    else
       annu_hyd = 1.e-30_db
    endif

  end function annu_hyd

end module modelatom_mod
