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

  implicit none

  type MODELATOM
     character(len=10) :: name
     integer :: charge
     integer :: nz
     real(double) :: mass ! amu
     real(double) :: abundance
     integer :: nLevels
     character(len=10), pointer :: level(:)
     real(double), pointer :: energy(:)  ! erg
     real(double), pointer :: g(:)
     real(double) :: iPot
     integer, pointer :: ionStage(:)
     integer, pointer :: nQuantum(:)
     integer :: nTrans
     integer :: nRBBTrans
     integer :: indexRBBTrans(200)
     character(len=3), pointer :: transType(:) ! RBB, RBF, CBB, CBF etc
     real(double), pointer :: transFreq(:)
     integer, pointer :: iLower(:)
     integer, pointer :: iUpper(:)
     integer, pointer :: equation(:)
     integer, pointer :: nParams(:)
     real(double), pointer :: params(:,:)
     real(double), pointer :: fMatrix(:,:)
  end type MODELATOM

  type TOPBASETYPE
     character(len=10) :: atom
     integer, pointer :: i(:)
     integer :: nz
     integer :: ne
     integer,pointer :: islp(:)
     integer, pointer  :: ilv(:)
     real(double), pointer :: e(:)
     integer :: nLevels
     integer, pointer :: nFreq(:)
     real(double), pointer :: freq(:,:)
     real(double), pointer :: xSection(:,:)
  end type TOPBASETYPE


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
    read(30,*) thisAtom%nz
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
       read(chunk(3), *) thisAtom%energy(i)
       read(chunk(4), *) thisAtom%g(i)
    enddo

    thisAtom%iPot = abs(thisAtom%energy(1))*hCgs*ergtoev
    thisAtom%energy = (thisAtom%energy(1) - thisAtom%energy)*hCgs*ergtoev ! eV

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

    allocate(thisAtom%fMatrix(1:thisAtom%nLevels,1:thisAtom%nLevels))
    thisAtom%fMatrix =0.d0
    do i = 1, thisAtom%nTrans
       if (thisAtom%transType(i) == "RBB") then
           thisAtom%fMatrix(thisAtom%iLower(i), thisAtom%iUpper(i)) = thisAtom%params(i,1)
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

    if (level == 0) then
       write(*,*) "Level not found: ",trim(cLevel)
       stop
    endif
  end function getLevel



  subroutine returnEinsteinCoeffs(thisAtom, iTrans, aEinstein, BulEinstein, BluEinstein)
    type(MODELATOM) :: thisAtom
    integer :: iTrans
    real(double) :: aEinstein, BulEinstein, BluEinstein, f
    integer :: iUpper, iLower

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


  real(double) function photoCrossSection(thisAtom, iLevel, nu)
    type(MODELATOM) :: thisAtom
    integer :: iLevel
    real(double) :: nu
    character(len=2) :: shell
    integer :: is
    real :: x
    
    if (iLevel > thisAtom%nLevels) then
       call writeFatal("photocrosssection: Level greater than nlevels")
       stop
    endif
    select case(thisAtom%name)
       case("HI")
          photoCrossSection = annu_hyd(iLevel, nu)
!          call phfit2(1,1,1,real(nu*hCgs*ergtoev), x)
!          photoCrossSection  = x * 1.d-10
       case("HeI")
          call phfit2(2,2,1,real(nu*hCgs*ergtoev), x)
          photoCrossSection  = x * 1.d-10
       case("HeII")
          call phfit2(2,1,1,real(nu*hCgs*ergtoev), x)
          photoCrossSection  = x * 1.d-10
       case DEFAULT
          call writeFatal("photocrosssection: atom not recognised")
          stop
     end select
   end function photoCrossSection


  function collisionRate(thisAtom, iTrans, temperature,label) result(rate)
    type(MODELATOM) :: thisAtom
    integer :: iTrans
    real(double) :: temperature, ne, rate, u0, u1, u2
    character(len=*), optional :: label
    real(double) :: logGamma
    real(double) :: sigma0
    real :: x
    real(double) :: fij, eh, gamma, logt, g
    integer :: i, j

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
       case("HI") 
          rate = cijt_hyd_hillier(thisAtom%iLower(iTrans), thisAtom%iUpper(iTrans), temperature, thisAtom)

       case("HeI") 
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
                i = thisAtom%iLower(iTrans)
                j = thisAtom%iUpper(iTrans)
                EH = hydE0eVdb * evtoErg
                u0 = Hcgs*thisAtom%transFreq(iTrans)/(kErg*temperature)
                u1 = u0 + 0.2d0
                rate = 5.465d-11 * sqrt(temperature) * 4.d0 * thisAtom%params(itrans,1) 
                rate = rate * (EH/(hCgs*thisAtom%transFreq(iTrans)))**2 * u0 * (expint(1, u0) - &
                     (u0/u1)*exp(-0.2d0)*expint(1,u1))
            case(8)
               logGamma = thisAtom%params(itrans,1) + thisAtom%params(iTrans,2) * &
                    log10(temperature)+thisAtom%params(itrans,3)/(log10(temperature)**2)
               rate = 5.465d-11 * sqrt(temperature) * exp(-(hCgs*thisAtom%transFreq(iTrans)/(kerg*temperature))) * (10.d0**logGamma)
             case(9)
                u0 = Hcgs*thisAtom%transFreq(iTrans)/(kErg*temperature)
                rate = 5.465d-11 * sqrt(temperature) * exp(-u0)*(1.d0+u0)
                rate = rate * 1.d-2!!!!!!!!!!

             case DEFAULT
                call writeFatal("CBB rate equation not implemented for He")
                stop
           end select
       case("HeII") 
          EH = hydE0eVdb * evtoErg
          u0 = Hcgs*thisAtom%transFreq(iTrans)/(kErg*temperature)
          fij = thisAtom%fMatrix(thisAtom%iLower(iTrans), thisAtom%iUpper(iTrans))
          rate =  5.465d-11 * sqrt(temperature) * (EH/(Hcgs*thisAtom%transFreq(iTrans)))**2 &
               * u0 * fij * exp(1.d0)* (u0*exp(-u0)*log(2.d0) + u0*expint(1,u0))
               

       case DEFAULT
          call writeFatal("collisionRate: bound-bound collision type not implemented")
          stop
       end select
    case("CBF")
       select case(thisAtom%name)
       case("HI") 
          rate = cikt_hyd_hillier(thisAtom%iLower(iTrans), temperature, thisAtom)
       case("HeI")
          if (thisAtom%iLower(iTrans) <= 15) then
             sigma0 = 1.64d0
             u0 = (thisAtom%iPot-thisAtom%energy(thisAtom%iLower(iTrans)))/(kEv*temperature)
             u1 = u0 + 0.27d0
             u2 = u0 + 1.43d0
             rate =  5.465d-11 * sqrt(temperature) *sigma0 * (u0*expint(1,u0) - (0.728d0*u0**2/u1)*expint(1,u1) - 0.189d0 * u0**2 &
                  * exp(-u0)*((2.d0+u2)/u2**3))
          else
             call phfit2(2,2,1,real(thisAtom%ipot*1.01), x)
             g = thisAtom%params(iTrans,2)
             sigma0  = x * 1.d-10
             u0 = (thisAtom%iPot - thisAtom%energy(thisAtom%iLower(itrans)))/(kEv*temperature)
             rate = 1.55d13*thisAtom%params(itrans,2)*g*sigma0*exp(-u0)/u0 /sqrt(temperature)
          endif
       case("HeII")
          if (thisAtom%iLower(iTrans) <= 10) then
             if (thisAtom%iLower(iTrans) <= 3) then
                gamma = thisAtom%params(iTrans,1)  + thisAtom%params(iTrans,2)*temperature + &
                     thisAtom%params(iTrans,4)/temperature + &
                     thisAtom%params(iTrans,5)/temperature**2
             else
                logt = log10(temperature)
                gamma = thisAtom%params(iTrans,1)  + thisAtom%params(iTrans,2)*logt + thisAtom%params(iTrans,3)*logt**2 +&
                     thisAtom%params(iTrans,5)/logt**2
             endif
             rate =  5.465d-11 * sqrt(temperature) * &
                  exp(-(thisAtom%iPot-thisAtom%energy(thisAtom%iLower(iTrans)))/(Kev*temperature)) * Gamma
             rate = 0.d0
          else
             call phfit2(2,1,1,real(thisAtom%ipot*1.01), x)
             sigma0  = x * 1.d-10
             u0 = (thisAtom%iPot - thisAtom%energy(thisAtom%iLower(itrans)))/(kEv*temperature)
             rate = 1.55d13*thisAtom%params(itrans,2)*sigma0*exp(-u0)/u0 /sqrt(temperature)
             rate = 0.d0
          endif
       case DEFAULT
          call writeFatal("collisionRate: bound-free collision type not implemented")
          stop
       end select
    case DEFAULT
       call writeFatal("collisionRate: collision type not implemented")
       stop
    end select
    if (rate <0.d0) then
       write(*,*) "negative CBF rate for atom ",thisAtom%name,thisAtom%iLower(iTrans),thisAtom%iUpper(iTrans)
       rate = 1.d-30
    endif
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
    integer :: level
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

  real(double) function annu_hyd(n,nu)
    !
    ! this subroutine calculates the photoionization x-section
    ! for hydrogen from the n-th level for a given freq photon.
    !
    integer, intent(in)              :: n     ! the level
    real(double), intent(in):: nu    ! the photon frequency
    real(double)            :: lam,e ! the photon wavelength
    integer :: level
    real(double), parameter :: eTrans(23) =                      &
         (/ (hydE0eVdb*(1.0d0 - 1.0d0/level**2),level=1,SIZE(eTrans)) /) 

    lam=cSpeed/nu
    lam=lam * 1.e8_db

    e = hCgs * nu * ergToEv

    if (e > (hydE0eV - eTrans(n))) then
       annu_hyd=1.044e-26_db*gii(n,1.e0_db,lam)*(lam**3)/dble(n)**5
    else
       annu_hyd = 1.e-50_db
    endif

  end function annu_hyd

  subroutine readTopbase(base, basefilename)
    character(len=*) :: basefilename
    character(len=200) :: dataDirectory, thisFilename
    type(TOPBASETYPE) :: base
    integer :: i, j
    character(len=80) :: junk
    integer, parameter :: maxFreq = 4000
    
    call writeInfo("Read TOPBASE data from: "//trim(baseFilename),TRIVIAL)

    call unixGetenv("TORUS_DATA", dataDirectory, i)
    thisfilename = trim(dataDirectory)//"/"//basefilename

    open(30, file=thisfilename, status="old", form="formatted")


    read(30,*) base%nLevels
    read(30,'(a80)') junk
    read(30,'(a80)') junk
    read(30,'(a80)') junk

    allocate(base%i(1:base%nLevels))
    allocate(base%islp(1:base%nLevels))
    allocate(base%ilv(1:base%nLevels))
    allocate(base%e(1:base%nLevels))
    allocate(base%nFreq(1:base%nLevels))
    allocate(base%freq(1:base%nLevels, 1:maxFreq))
    allocate(base%xSection(1:base%nLevels, 1:maxFreq))

    do i = 1, base%nLevels
       read(30,*) base%i(i), base%nz, base%ne, base%islp(i), base%ilv(i), base%e(i), base%nFreq(i)
       if (base%nFreq(i) > maxFreq) then
          write(*,*) "Not enough freq space ",base%nfreq(i),maxFreq
          stop
       endif
       do j = 1, base%nFreq(i)
          read(30,*) base%freq(i,j), base%xSection(i,j)
       enddo
    enddo
    close(30)

    base%e = abs(base%e)*rydbergtoev
    base%e = base%e(1) - base%e
    call writeInfo("Done.",TRIVIAL)



  end subroutine readTopbase

  function BoltzSahaGeneral(thisAtom, nion, level, Ne, t) result(ratio)
  
    type(MODELATOM) :: thisAtom
    integer              :: nIon, level
    real(double) :: Ne, t, ratio, nk
    real(double) ::  Ucoeff(5,5)
    real(double) :: N2, N1, N0, u0, u1, u2, N1overN0, N2overN1, pe, tot
    uCoeff(1,1:5) = (/0.30103d0, -0.00001d0, 0.d0, 0.d0, 0.d0 /)
    uCoeff(2,1:5) = (/0.00000d0,  0.00000d0, 0.d0, 0.d0, 0.d0 /)
    uCoeff(3,1:5) = (/0.30103d0,  0.00000d0, 0.d0, 0.d0, 0.d0 /)

    pe = ne * kerg * t
    select case(thisAtom%name)
       case("HI")
          u0 = getUT(t, uCoeff(1,1:5))
          u1 = 1.d0
          N1overN0 = ((-5040.d0/t)*thisAtom%iPot + 2.5d0*log10(t) + log10(u1/u0)-0.1762d0)
          N1overN0 = (10.d0**N1overN0)/pe

!          n1overn0 = (2.d0*u1)/(ne*u0) * ((twoPi * melectron * kerg * t)/(hCgs**2))**1.5d0 * exp(-thisAtom%iPot/(kev*t))

          ratio = (thisAtom%g(level)*exp(-thisAtom%energy(level)/(kev * t))) / u0 / N1overn0

       case("HeI")
          u0 = 1.d0
          u1 = getUT(t, uCoeff(3,1:5))
          N1overN0 = ((-5040.d0/t)*24.59d0 + 2.5d0*log10(t) + log10(u1/u0)-0.1762d0)
          N1overN0 = (10.d0**N1overN0)/pe
          ratio = (thisAtom%g(level)*exp(-thisAtom%energy(level)/(kev * t))) / u0 / N1overn0
       case("HeII")
          u0 = 1.d0
          u1 = getUT(t, uCoeff(3,1:5))
          u2 = 1.d0
          N1overN0 = ((-5040.d0/t)*24.59d0 + 2.5d0*log10(t) + log10(u1/u0)-0.1762d0)
          N1overN0 = (10.d0**N1overN0)/pe

          N2overN1 = ((-5040.d0/t)*54.42d0 + 2.5d0*log10(t) + log10(u2/u1)-0.1762d0)
          N2overN1 = (10.d0**N2overN1)/pe
          ratio = (thisAtom%g(level)*exp(-thisAtom%energy(level)/(kev * t))) / u1 / N2overn1
       case DEFAULT
          write(*,*) "atom not recognised in boltzsahageneral ",thisAtom%name
          stop
     end select
   end function BoltzSahaGeneral
          

   function getUT(t, coeff) result (ut)
     real(double) :: t, coeff(:),  ut
     real(double) :: logphi

     logphi = 5040.d0/t

     ut = coeff(1) + coeff(2)*logphi + coeff(3) * logphi**2 + &
          coeff(4) * logphi**3 + coeff(5) * logphi**4
     ut = 10.d0**ut
   end function getUT

   subroutine createRRBarrays(nAtom, thisAtom, nRBBtrans, indexAtom, indexRBBTrans)
     integer :: nAtom
     type(MODELATOM) :: thisAtom(:)
     integer :: nRBBTrans
     integer :: indexAtom(:)
     integer :: indexRBBTrans(:)
     integer :: iATom, iTrans

     nRBBTrans = 0

     do iAtom = 1, nAtom
        thisAtom(iAtom)%nRBBTrans = 0
        do iTrans = 1, thisAtom(iAtom)%nTrans
           if (thisAtom(iAtom)%transType(iTrans) == "RBB") then
              nRBBTrans = nRBBTrans + 1
              indexAtom(nRBBTrans) = iAtom
              indexRBBTrans(nRBBTrans) = iTrans
              thisAtom(iAtom)%nRBBTrans = thisAtom(iAtom)%nRBBTrans + 1
              thisAtom(iAtom)%indexRBBtrans(iTrans) = thisAtom(iAtom)%nRBBTrans
           endif
        enddo
     enddo
   end subroutine createRRBarrays


  real(double) function oldcikt_ma(i,t, thisAtom) 
    !
    ! this function calculates the collisional ionization rate (see k&c and
    ! mihalas 1967 (apj 149 169) ).
    ! 
    integer, intent(in)              :: i        ! the level
    real(double),intent(in) :: t        ! the temperature
    type(MODELATOM) :: thisAtom
    real(double) :: crap
    real(double) :: t1                  
    real(double) :: gammait             ! see k&c
    real(double) :: lgt
    real(double) :: chi                 ! potential from i to k
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
    if (i .ne. 2) then 
       gammait=cint(1,i)+cint(2,i)*lgt+cint(3,i)*(lgt**2)+ &
            (cint(4,i)/lgt)+(cint(5,i)/(lgt**2))
    else
       gammait=cint(1,i)+cint(2,i)*t1+(cint(4,i)/t1)+(cint(5,i)/t1**2)
    endif
    chi=hydE0eV-thisAtom%energy(i)
    oldcikt_ma=(5.465e-11)*sqrt(t1)*exp(-chi/(kev*t1))*gammait

  end function oldcikt_ma

end module modelatom_mod
