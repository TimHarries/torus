#ifdef CMFATOM
module modelatom_mod

  ! module for model atoms created by tjh 23/10/06

  use kind_mod
  use constants_mod
  use messages_mod
  use unix_mod, only: unixGetenv
  use hyd_col_coeff, only: tempTable, omegaij, omegaik
  use stateq_mod, only:  expint
  use utils_mod, only: locate, sort
  use source_mod, only: SOURCETYPE
  use phfit_mod, only : phfit2

  implicit none

  type MODELATOM
     character(len=10) :: name
     integer :: charge
     integer :: nz
     real(double) :: mass ! amu
     integer :: nLevels
     character(len=10), pointer :: level(:)
     real(double), pointer :: energy(:)  ! erg
     real(double), pointer :: g(:)
     real(double) :: iPot
     integer :: nPhotoFreq
     real(double), pointer :: photoFreq(:) => null()
     real(double), pointer :: photoSec(:,:) => null()
     integer, pointer :: ionStage(:)
     integer, pointer :: nQuantum(:)
     integer :: nTrans
     integer :: nRBBTrans
     integer :: indexRBBTrans(200)
     integer :: nRBFTrans
     integer :: indexRBFTrans(200)
     character(len=3), pointer :: transType(:) ! RBB, RBF, CBB, CBF etc
     real(double), pointer :: transFreq(:)
     integer, pointer :: iLower(:)
     integer, pointer :: iUpper(:)
     logical, pointer :: inDetailedBalance(:)
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

  type(MODELATOM), pointer :: globalAtomArray(:)

contains


  subroutine identifyTransitionCmf(lamLine, atomArray, iAtom, iTrans)
    real(double) :: lamLine, lamTrans
    type(MODELATOM) :: atomArray(:)
    integer :: iAtom, iTrans
    logical :: found
    integer :: i, j
    character(len=80) :: message
    found = .false.
    do i = 1, SIZE(atomArray)
       do j = 1, atomArray(i)%nTrans
          if (atomArray(i)%transType(j)=="RBB") then
             lamTrans = ((cspeed/atomArray(i)%transfreq(j))*1d8)/nAir !!tjgw201 (04/02/19) this will not match wavelengths in H.atm becuse of correction for air (rather than vacuum)
             if ( (abs(lamTrans-lamLine)) < 1.d0) then
                found = .true.
                iAtom = i
                iTrans = j
                exit
             endif
          endif
       enddo
       if (found) exit
    enddo
    if (found) then
       write(message, '(a, a, a, i2)') "Transition found: name=",trim(AtomArray(iAtom)%name), ", iTrans=",iTrans
       call writeInfo(message,TRIVIAL)
    else
       write(message,'(a,f10.1)') "Transition not found in identifyTransition: ",lamline
       call writeFatal(message)
       stop
    endif
  end subroutine identifyTransitionCmf

  subroutine readAtom(thisAtom, atomfilename)
    character(len=*) :: atomfilename
    character(len=200) :: dataDirectory, thisfilename
    type(MODELATOM) :: thisAtom
    integer :: i, nChunks, j
    character(len=120) :: junk
    character(len=20) :: chunk(20)
    character(len=80) :: message
    real(double) :: lamTrans
    dataDirectory = " "; chunk = " "; nChunks = 0

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
       thisAtom%level(i) = chunk(1)(1:10)
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
    allocate(thisAtom%inDetailedBalance(1:thisAtom%nTrans))

    thisAtom%inDetailedBalance = .false.

    do i = 1, thisAtom%nTrans
       read(30,'(a120)') junk
       call splitIntoChunks(junk, nChunks, chunk)
       thisAtom%transType(i) = chunk(1)(1:3)
       thisAtom%iLower(i) = getLevel(thisAtom, chunk(2))
       thisAtom%iUpper(i) = getLevel(thisAtom, chunk(3))
       read(chunk(4),*) thisAtom%equation(i)
       read(chunk(5),*) thisAtom%nParams(i)
       do j = 1, thisAtom%nParams(i)
          read(chunk(5+j),*) thisAtom%params(i, j)
       enddo
       if (thisAtom%iUpper(i) <= thisAtom%nLevels) then
          thisAtom%transFreq(i) = ((thisAtom%energy(thisAtom%iUpper(i)) - thisAtom%energy(thisAtom%iLower(i)))*evToErg)/Hcgs
          lamTrans = ((cspeed/thisAtom%transfreq(i))*1d8)/nAir
          !print*,"freq",thisAtom%transFreq(i)," iupper",thisAtom%iUpper(i)," ilower",thisatom%iLower(i)," Wavelength",lamTrans


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

  function returnNe(pops, thisAtom) result(ne)
    real(double) :: pops(:,:), ne
    type(MODELATOM) :: thisAtom(:)
    integer :: nAtom, iAtom

    ne = 0.d0
    nAtom = size(thisAtom)
    do iAtom = 1,nAtom
       ne = ne + pops(iAtom,thisatom(iatom)%nlevels) * (dble(thisAtom(iatom)%charge)+1.d0)
    enddo
  end function returnNe


  subroutine returnEinsteinCoeffs(thisAtom, iTrans, aEinstein, BulEinstein, BluEinstein)
    type(MODELATOM) :: thisAtom
    integer :: iTrans
    real(double) :: aEinstein, BulEinstein, BluEinstein, f
    integer :: iUpper, iLower
    real(double), parameter :: fac = (8.d0 * pi**2 * eCharge**2)/(mElectron*cSpeed**3)
    real(double), parameter :: fac2 = (fourPi * pi * eCharge**2)/(mElectron * hCgs *  cspeed)
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
!       aEinstein = (8.d0*thisAtom%transFreq(iTrans)**2 * pi**2 * eCharge**2) / &
!            (mElectron*cSpeed**3) * f * thisAtom%g(iLower) / thisAtom%g(iUpper)
!       BluEinstein = (fourPi * pi * eCharge**2)/(mElectron * hCgs * thisAtom%transFreq(iTrans) * cspeed) * f

       aEinstein = fac * thisAtom%transFreq(iTrans)**2 * f * thisAtom%g(iLower) / thisAtom%g(iUpper)


       BluEinstein = fac2 /  thisAtom%transFreq(iTrans) * f
       BulEinstein = BluEinstein * thisAtom%g(iLower) / thisAtom%g(iUpper)
    case DEFAULT
       call writeFatal("returnEinsteinCoeffs: unrecognised equation for RBB transition")
       stop
    end select
  end subroutine returnEinsteinCoeffs


  real(double) function photoCrossSection(thisAtom, iLevel, nu)
    type(MODELATOM) :: thisAtom
    integer :: iLevel
    real(double) :: nu, e
    real :: x
    photoCrossSection = tiny(photoCrossSection)


    if (iLevel > thisAtom%nLevels) then
       call writeFatal("photocrosssection: Level greater than nlevels")
       stop
    endif
    select case(thisAtom%nz)
       case(1)
          photoCrossSection = annu_hyd(iLevel, nu)
!          call phfit2(1,1,1,real(nu*hCgs*ergtoev), x)
!          photoCrossSection  = x * 1.d-10
       case(2)
          select case(thisAtom%charge)
             case(0)
                photoCrossSection = annu_He_I(ilevel,nu)
             case(1)
                photoCrossSection =  annu_He_II(ilevel,nu)
          end select
       case(20)
          e = nu * hcgs * ergtoev
          select case(thisAtom%charge)
             case(0)
                call phfit2(20,20,1,real(e), x)
                photoCrossSection = x * 1.d-10
             case(1)
             call phfit2(20,19,1,real(e), x)
                photoCrossSection = x * 1.d-10
          end select
      case DEFAULT
          call writeFatal("photocrosssection: atom not recognised")
          stop
     end select
   end function photoCrossSection

   real(double) function gauntTest(nu, gIIx, gIIy, gIIz)
     real(double) :: nu, gIIx, gIIy, gIIz
     gauntTest = gIIx + gIIy / nu + GIIz / nu**2
   end function gauntTest

  real(double) function quickPhotoCrossSection(thisAtom, iRBF, iFreq)
    type(MODELATOM) :: thisAtom
    integer :: iRBF, iFreq
    quickPhotoCrossSection  = thisAtom%photoSec(iRBF, iFreq)
  end function quickPhotoCrossSection




   subroutine addCrossSectionstoAtom(thisAtom, nFreq, freq)
     type(MODELATOM) :: thisAtom
     integer :: nFreq
     real(double) :: freq(:)
     integer :: iFreq, j, iTrans, iLevel

     if (associated(thisAtom%photoFreq)) deallocate(thisAtom%photoFreq)
     if (associated(thisAtom%photoSec)) deallocate(thisAtom%photoSec)

     thisAtom%nPhotoFreq = nFreq
     allocate(thisAtom%photoFreq(1:nfreq))
     allocate(thisAtom%photoSec(1:thisAtom%nRBFtrans, 1:nFreq))

     do j = 1, thisAtom%nRBFTrans
        iTrans = thisAtom%indexRBFtrans(j)
        iLevel = thisAtom%iLower(iTrans)
        do iFreq = 1, nFreq
           thisAtom%photoSec(j, iFreq) = photoCrossSection(thisAtom, iLevel, freq(ifreq))
        enddo
     enddo



   end subroutine addCrossSectionstoAtom




  function collisionRate(thisAtom, iTrans, temperature) result(rate)
    type(MODELATOM) :: thisAtom
    integer :: iTrans
    real(double) :: temperature, rate, u0, u1, u2
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

             case DEFAULT
                call writeFatal("CBB rate equation not implemented for He")
                stop
           end select
       case("HeII")
          EH = hydE0eVdb * evtoErg
          u0 = Hcgs*thisAtom%transFreq(iTrans)/(kErg*temperature)
          fij = thisAtom%fMatrix(thisAtom%iLower(iTrans), thisAtom%iUpper(iTrans))
          rate =  5.465d-11 * sqrt(temperature) * (EH/(Hcgs*thisAtom%transFreq(iTrans)))**2 &
               * fij * exp(1.d0)* (u0*exp(-u0)*log(2.d0) + u0*expint(1,u0))


          rate = 5.465d-11 * sqrt(temperature) * exp(1.d0) * fij * (1.d0/16.d0) * &
               (1.d0/dble(thisAtom%ilower(itrans))**2 - 1.d0/dble(thisAtom%iUpper(itrans))**2) * &
               (u0 * exp(-u0) * log(2.d0) + u0 * expint(1,u0))


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
          else
             call phfit2(2,1,1,real(thisAtom%ipot*1.01), x)
             sigma0  = x * 1.d-10
             u0 = (thisAtom%iPot - thisAtom%energy(thisAtom%iLower(itrans)))/(kEv*temperature)
             rate = 1.55d13*thisAtom%params(itrans,2)*sigma0*exp(-u0)/u0 /sqrt(temperature)
          endif
!          rate =  HeIICollisionalIonRates(thisAtom%ilower(itrans), &
!               thisAtom%iPot - thisAtom%energy(thisAtom%iLower(itrans)) , temperature)
!         write(*,*) "rate ",rate, HeIICollisionalIonRates(thisAtom%ilower(itrans), &
!               thisAtom%iPot - thisAtom%energy(thisAtom%iLower(itrans)) , temperature)
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
    integer               :: i1,j1           ! lower/upper level
    integer :: lower, upper
    real(double) :: factor
    integer :: level
    real(double), parameter :: eTrans(23) =                      &
         (/ (hydE0eVdb*(1.0d0 - 1.0d0/level**2),level=1,SIZE(eTrans)) /)

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
!    real(double) :: gammait             ! see k&c
    real(double) :: chi                 ! potential from i to k
    integer :: lower, upper
    real(double) :: factor
    integer :: level
    real(double), parameter :: eTrans(23) =                      &
         (/ (hydE0eVdb*(1.0d0 - 1.0d0/level**2),level=1,SIZE(eTrans)) /)

    ! making cint a PARAMETER may cause problems with XL Fortran
    !  real(double) :: cint(5,10)
    !  cint = reshape(source=                                                                   &
!!$    real(double), parameter :: cint(5,10) = reshape(source=                         &
!!$         (/-0.4350000e0_db, 0.3000000e00_db,0.0000000e0_db, 0.0000000e0_db,0.00000000e0_db,  &
!!$         1.9987261e1_db,-5.8906298e-5_db,0.0000000e0_db,-2.8185937e4_db,5.44441600e7_db,  &
!!$         1.3935312e3_db,-1.6805859e02_db,0.0000000e0_db,-2.5390000e3_db,0.00000000e0_db,  &
!!$         2.0684609e3_db,-3.3415820e02_db,0.0000000e0_db, 0.0000000e0_db,-7.6440625e3_db,  &
!!$         3.2174844e3_db,-5.5882422e02_db,0.0000000e0_db, 0.0000000e0_db,-6.8632500e3_db,  &
!!$         5.7591250e3_db,-1.5163125e03_db,8.1750000e1_db, 0.0000000e0_db,0.00000000e0_db,  &
!!$         1.4614750e4_db,-4.8283750e03_db,3.9335938e2_db, 0.0000000e0_db,0.00000000e0_db,  &
!!$         2.8279250e4_db,-1.0172750e04_db,9.1967968e2_db, 0.0000000e0_db,0.00000000e0_db,  &
!!$         4.6799250e4_db,-1.7627500e04_db,1.6742031e3_db, 0.0000000e0_db,0.00000000e0_db,  &
!!$         -7.4073000e4_db, 6.8599375e03_db,0.0000000e0_db, 2.0169800e5_db,0.00000000e0_db/),&
!!$         shape=(/5,10/))

    t1 = min(t,1.5e5_db)
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
         (/ (hydE0eVdb*(1.0d0 - 1.0d0/dble(level)**2),level=1,SIZE(eTrans)) /)

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
    dataDirectory = " "
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

  function BoltzSahaGeneral(thisAtom, level, Ne, t) result(ratio)

    type(MODELATOM) :: thisAtom
    integer              ::  level
    real(double) :: Ne, t, ratio
    real :: tReal
    real, parameter ::  Ucoeff(5,5) =  reshape(source=  &
        (/0.30103e0, -0.00001e0, 0.00000e0,  0.00000e0, 0.00000e0, &
          0.00000e0,  0.00000e0, 0.00000e0,  0.00000e0, 0.00000e0, &
          0.30103e0,  0.00000e0, 0.00000e0,  0.00000e0, 0.00000e0, &
          0.07460e0, -0.75759e0, 2.58494e0, -3.53170e0, 1.65240e0, &
          0.34383e0, -0.41472e0, 1.01550e0,  0.31930e0, 0.00000e0  &
           /), shape=(/5,5/), order=(/2,1/))
    real(double) :: u0, u1, u2, N1overN0, N2overN1, pe, ci
    real(double) :: n0overn1,n1overn2
       pe = ne * kerg * t
       tReal = real(t)
       select case(thisAtom%nz)
       case(1)
          u0 = getUT(treal, uCoeff(1,:))
          u1 = 1.d0
          N1overN0 = ((-5040.d0/t)*thisAtom%iPot + 2.5d0*log10(t) + log10(u1/u0)-0.1762d0)
          if ((t > 3.).and.(pe /= 0.d0)) then
             N1overN0 = (10.d0**N1overN0)/pe
          else
             N1overN0 = 0.d0
          endif

          !          n1overn0 = (2.d0*u1)/(ne*u0) * ((twoPi * melectron * kerg * t)/(hCgs**2))**1.5d0 * exp(-thisAtom%iPot/(kev*t))

          if (N1overN0 /= 0.d0) then
             ratio = (thisAtom%g(level)*exp(-thisAtom%energy(level)/(kev * t))) / u0 / N1overn0
          else
             if (level == 1) then
                ratio = 1.d0
             else
                ratio = 0.d0
             endif
          endif
       case(2)

          ci = 2.07d-16
          u0 = getUT(treal, uCoeff(2,:))
          u1 = getUT(treal, uCoeff(3,:))
          u2 = 1.d0


          if (t > 3.) then
!             write(*,*) "t ",t, " ne ",ne
             N0overN1 = ne * (u0/u1) * ci * exp (24.59d0/(kev*t)) / t**1.5
             n1overn2 = ne * (u1/u2) * ci * exp (54.42d0/(kev*t)) / t**1.5
          else
             n0overn1 = 1.d10
             n1overn2 = 1.d10
          endif

          select case(thisAtom%charge)
          case(0)
             ratio = ((thisAtom%g(level)*exp(-thisAtom%energy(level)/(kev * t))) / u0) * N0overn1
          case(1)
             ratio = ((thisAtom%g(level)*exp(-thisAtom%energy(level)/(kev * t))) / u1) * N1overn2
          case DEFAULT
             write(*,*) "wrong charge for He"
             stop
          end select
       case(20)
          select case(thisAtom%charge)
          case(0)
             u0 = 1.d0
             u1 = getUT(treal, uCoeff(4,:))
             N1overN0 = ((-5040.d0/t)*6.11d0 + 2.5d0*log10(t) + log10(u1/u0)-0.1762d0)
             N1overN0 = (10.d0**N1overN0)/pe
             ratio = (thisAtom%g(level)*exp(-thisAtom%energy(level)/(kev * t))) / u0 / N1overn0
          case(1)
             u1 = getUT(treal, uCoeff(4,:))
             u2 = getUT(treal, uCoeff(5,:))
             N2overN1 = ((-5040.d0/t)*11.87d0 + 2.5d0*log10(t) + log10(u2/u1)-0.1762d0)
             N2overN1 = (10.d0**N2overN1)/pe

             ratio = (thisAtom%g(level)*exp(-thisAtom%energy(level)/(kev * t))) / u1 / N2overn1
          case DEFAULT
             write(*,*) "wrong charge for Ca"
             stop
          end select
       case DEFAULT
          write(*,*) "atom not recognised in boltzsahageneral ",thisAtom%name
          stop
       end select
!       write(*,*) ratio,ne,t
!    if (isnan(ratio)) then
!       write(*,*) "bug in boltz",ne,t
!       stop
!    endif
   end function BoltzSahaGeneral

  function BoltzSahaEquationAbsolute(thisAtom, nTotal, level, Ne, t) result(pop)

    type(MODELATOM) :: thisAtom
    real(double) :: Ne, t, pop, nTotal
    real :: tReal
    integer :: level
    real(double) :: c
    real, parameter ::  Ucoeff(5,5) =  reshape(source=  &
        (/0.30103e0, -0.00001e0, 0.00000e0,  0.00000e0, 0.00000e0, &
          0.00000e0,  0.00000e0, 0.00000e0,  0.00000e0, 0.00000e0, &
          0.30103e0,  0.00000e0, 0.00000e0,  0.00000e0, 0.00000e0, &
          0.07460e0, -0.75759e0, 2.58494e0, -3.53170e0, 1.65240e0, &
          0.34383e0, -0.41472e0, 1.01550e0,  0.31930e0, 0.00000e0  &
           /), shape=(/5,5/),order=(/2,1/))
    real(double) :: u0, u1, u2, N1overN0,  pe, ci
    real(double) :: n0overn1,n1overn2
       pe = ne * kerg * t
       tReal = real(t)
       select case(thisAtom%nz)
       case(1)
          u0 = getUT(treal, uCoeff(1,:))
          u1 = 1.d0
          N1overN0 = ((-5040.d0/t)*thisAtom%iPot + 2.5d0*log10(t) + log10(u1/u0)-0.1762d0)
          if ((t > 3.).and.(pe > 1.d-30)) then
             N1overN0 = (10.d0**N1overN0)/pe
          else
             N1overN0 = 0.d0
          endif

          pop = nTotal / (1.d0+n1OverN0)
          pop = pop * thisAtom%g(level)*exp(-thisAtom%energy(level)/(kev * t))/u0
       case(2)

          ci = 2.07d-16
          u0 = getUT(treal, uCoeff(2,:))
          u1 = getUT(treal, uCoeff(3,:))
          u2 = 1.d0


          if (t > 3.) then
             N0overN1 = min(1.d20,ne * (u0/u1) * ci * exp (24.59d0/(kev*t)) / t**1.5)
             n1overn2 = min(1.d20,ne * (u1/u2) * ci * exp (54.42d0/(kev*t)) / t**1.5)

             c = 1.d0 / ( 1.d0 + n1overn2 + n0overn1 * n1overn2)
          else
             n0overn1 = 1.d10
             n1overn2 = 1.d10
             c = 1.d-20
          endif

          select case(thisAtom%charge)
          case(0)
             pop = (nTotal * c) * n1overn2 * n0overn1
             pop = pop * thisAtom%g(level)*exp(-thisAtom%energy(level)/(kev * t))/u0
          case(1)
             pop = (nTotal * c) * n1overn2
             pop = pop * thisAtom%g(level)*exp(-thisAtom%energy(level)/(kev * t))/u1
          case DEFAULT
             write(*,*) "wrong charge for He"
             stop
          end select
       case(20)
             write(*,*) "wrong charge for Ca"
             stop
       case DEFAULT
          write(*,*) "atom not recognised in sahaEquationAbsolute ",thisAtom%name
          stop
       end select
   end function BoltzSahaEquationAbsolute

  function SahaEquationNextIon(thisAtom, nTotal, Ne, t) result(pop)

    type(MODELATOM) :: thisAtom
    real(double) :: Ne, t, pop, nTotal
    real :: tReal
    real(double) :: c
    real, parameter ::  Ucoeff(5,5) =  reshape(source=  &
        (/0.30103e0, -0.00001e0, 0.00000e0,  0.00000e0, 0.00000e0, &
          0.00000e0,  0.00000e0, 0.00000e0,  0.00000e0, 0.00000e0, &
          0.30103e0,  0.00000e0, 0.00000e0,  0.00000e0, 0.00000e0, &
          0.07460e0, -0.75759e0, 2.58494e0, -3.53170e0, 1.65240e0, &
          0.34383e0, -0.41472e0, 1.01550e0,  0.31930e0, 0.00000e0  &
           /), shape=(/5,5/), order=(/2,1/))
    real(double) :: u0, u1, u2, N1overN0,  pe, ci
    real(double) :: n0overn1,n1overn2
       pe = ne * kerg * t
       tReal = real(t)
       select case(thisAtom%nz)
       case(1)
          u0 = getUT(treal, uCoeff(1,:))
          u1 = 1.d0
          N1overN0 = ((-5040.d0/t)*thisAtom%iPot + 2.5d0*log10(t) + log10(u1/u0)-0.1762d0)
          if ((t > 3.).and.(pe >  1.d-30)) then
             N1overN0 = (10.d0**N1overN0)/pe
          else
             N1overN0 = 1.d-30
          endif

          pop = nTotal / (1.d0 + 1.d0/ max(1.d-20,n1OverN0))
       case(2)

          ci = 2.07d-16
          u0 = getUT(treal, uCoeff(2,:))
          u1 = getUT(treal, uCoeff(3,:))
          u2 = 1.d0


          if (t > 3.) then
             N0overN1 = min(1.d20,ne * (u0/u1) * ci * exp (24.59d0/(kev*t)) / t**1.5)
             n1overn2 = min(1.d20,ne * (u1/u2) * ci * exp (54.42d0/(kev*t)) / t**1.5)
             c = 1.d0 / ( 1.d0 + n1overn2 + n0overn1 * n1overn2)
          else
             n0overn1 = 1.d10
             n1overn2 = 1.d10
             c = 1.d-20
          endif

          select case(thisAtom%charge)
          case(0)
             pop = (nTotal * c) * n1overn2
          case(1)
             pop = (nTotal * c)
          case DEFAULT
             write(*,*) "wrong charge for He"
             stop
          end select
       case(20)
             write(*,*) "wrong charge for Ca"
             stop
       case DEFAULT
          write(*,*) "atom not recognised in sahaEquationAbsolute ",thisAtom%name
          stop
       end select
   end function SahaEquationNextIon


   function getUT(t, coeff) result (ut)
     real :: t, coeff(:),  ut
     real :: logphi

     logphi = log10(5040.0/t)

     ut = coeff(1) + coeff(2)*logphi
     ut = 10.**ut
   end function getUT

   subroutine createRBBarrays(nAtom, thisAtom, nRBBtrans, indexAtom, indexRBBTrans)
     integer :: nAtom
     type(MODELATOM) :: thisAtom(:)
     integer :: nRBBTrans
     integer :: indexAtom(:)
     integer :: indexRBBTrans(:)
     integer :: iATom, iTrans

     nRBBTrans = 0

     do iAtom = 1, nAtom
        thisAtom(iAtom)%nRBBTrans = 0
        thisAtom(iAtom)%nRBFTrans = 0
        do iTrans = 1, thisAtom(iAtom)%nTrans
           if (thisAtom(iAtom)%transType(iTrans) == "RBB") then
              nRBBTrans = nRBBTrans + 1
              indexAtom(nRBBTrans) = iAtom
              indexRBBTrans(nRBBTrans) = iTrans
              thisAtom(iAtom)%nRBBTrans = thisAtom(iAtom)%nRBBTrans + 1
              thisAtom(iAtom)%indexRBBtrans(itrans) = nRBBTrans
           endif
           if (thisAtom(iAtom)%transType(iTrans) == "RBF") then
              thisAtom(iAtom)%nRBFTrans = thisAtom(iAtom)%nRBFTrans + 1
              thisAtom(iAtom)%indexRBFtrans(thisAtom(iAtom)%nRBFTrans) = iTrans
           endif
        enddo
!        write(*,*) "nRBBtrans",thisAtom(iAtom)%nRBBTrans, &
!             thisAtom(iatom)%indexRBBTrans(1:thisAtom(iatom)%nrbbtrans)
!        write(*,*) "nRBFtrans",thisAtom(iAtom)%nRBFTrans, &
!             thisAtom(iatom)%indexRBFTrans(1:thisAtom(iatom)%nrbftrans)
     enddo
   end subroutine createRBBarrays


  real(double) function oldcikt_ma(i,t, thisAtom)
    !
    ! this function calculates the collisional ionization rate (see k&c and
    ! mihalas 1967 (apj 149 169) ).
    !
    integer, intent(in)              :: i        ! the level
    real(double),intent(in) :: t        ! the temperature
    type(MODELATOM) :: thisAtom
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


  function etaContHydrogen(freq, thisAtom, pops, nLevels, ne, temperature) result(eta)

! bound-free and free-free continuum emissivity

    type(MODELATOM) :: thisAtom
    real(double) :: freq
    real(double) :: pops(:)
    integer :: nLevels
    real(double) :: ne
    real(double) :: temperature
    real(double) :: nStar
    real(double) :: eta
    integer :: j
    real(double) :: thresh, photonEnergy

    eta=0.d0
!    bound-free

    do j=1, nLevels-1

       thresh=(thisAtom%iPot - thisAtom%energy(j))
       photonEnergy = freq * hCgs * ergtoEv
       if (photonEnergy.ge.thresh) then
          nStar = BoltzSahaGeneral(thisAtom, j, Ne, temperature) * pops(nLevels)
          eta = eta + nStar * annu_hyd(j,freq) * exp(-(hcgs*freq)/(kerg*temperature))
       endif
    enddo

! free-free

    eta = eta + (ne**2) * alpkk_hyd(freq,temperature) * exp(-(hcgs*freq)/(kerg*temperature))/fourPi

    eta=eta*real((2.0*dble(hcgs)*dble(freq)**3)/(dble(cspeed)**2))

  end function etaContHydrogen

  Real(double) pure function alpkk_hyd(freq,t)
     !
     ! this function returns the free-free absorption coefficient for hydrogen
     !
     real(double), intent(in) :: freq,t
     real(double)             :: wav,gauntf

     wav=1.e8_db*cSpeed/freq
     gauntf=giii_hyd(1.e0_db,t,wav)
     alpkk_hyd=gauntf*real(3.6d8/((dble(freq)**3)*sqrt(dble(t))))

   end function alpkk_hyd


  real(double) pure function giii_hyd (z, t, wl)
     !
     !   ferland's fabulous functional fits
     !

     real(double), intent(in) :: wl, t, z
     real(double) :: c, u, ulog, gam2
     integer               :: i,j,k, m
     real(double) :: b2
     real(double) :: frac, sum1, sum2, d
     ! making coeff and a2 PARAMETERs may cause problems with XL Fortran
     !  real(double) :: coeff(28)
     !  real(double) :: a2(7)
     !  coeff =                                                           &
     !     (/1.102d0       ,-0.1085d0     ,0.09775d0     ,-0.01125d0     ,&
     !       1.2d0         ,-0.24016667d0 ,0.07675d0     ,-0.01658333d0  ,&
     !       1.26d0        ,-0.313166667d0,0.15075d0     ,0.00241667d0   ,&
     !       1.29d0        ,-0.4518333d0  ,0.12925d0     ,0.00258333d0   ,&
     !       1.27d0        ,-0.579d0      ,0.092d0       ,-0.003d0       ,&
     !       1.16d0        ,-0.707333d0   ,0.112d0       ,0.0053333333d0 ,&
     !       0.883d0       ,-0.76885d0    ,0.190175d0    ,0.022675d0     /)
     !  a2 = (/100.d0, 10.d0, 3.d0, 1.d0, 0.3d0, 0.1d0, 0.001d0/)

     real(double), parameter :: coeff(28) =                   &
        (/1.102d0       ,-0.1085d0     ,0.09775d0     ,-0.01125d0     ,&
          1.2d0         ,-0.24016667d0 ,0.07675d0     ,-0.01658333d0  ,&
          1.26d0        ,-0.313166667d0,0.15075d0     ,0.00241667d0   ,&
          1.29d0        ,-0.4518333d0  ,0.12925d0     ,0.00258333d0   ,&
          1.27d0        ,-0.579d0      ,0.092d0       ,-0.003d0       ,&
          1.16d0        ,-0.707333d0   ,0.112d0       ,0.0053333333d0 ,&
          0.883d0       ,-0.76885d0    ,0.190175d0    ,0.022675d0     /)
     real(double), parameter :: a2(7) =                       &
          (/100.d0, 10.d0, 3.d0, 1.d0, 0.3d0, 0.1d0, 0.001d0/)

       u = 1.44e+8 / (wl*t)
       ulog = log10(u)
       gam2 = 1.58e+5 * z*z/t
       if (gam2.gt.a2(7)) go to 10
         i = 7
         j = 7
         k = 7
         frac = 0.5
         go to 60
   10  continue
       if (gam2.lt.a2(1)) go to 20
         i = 1
         j = 1
         k = 1
         frac = 0.5
         go to 60
   20  continue
       do 30 i = 2, 7
           if (gam2.gt.a2(i)) go to 40
   30  continue
   40  continue
       k = i - 1

       b2 = log10(a2(k))
       c = log10(a2(i))
       gam2 = log10(gam2)
       frac = abs ((gam2-b2) / (b2-c))
   60  continue
       k = (k-1)*4
       sum1 = coeff(k+1)
       d = 1.0
       do 70 m = 2, 4
           d = d*ulog
           sum1 = sum1 + coeff(k+m)*d
   70  continue
       sum1 = sum1 * (1.0 - frac)
       i = (i-1)*4
       sum2 = coeff(i+1)
       d = 1.0
       do 80 m = 2, 4
           d = d*ulog
           sum2 = sum2 + coeff(i+m)*d
   80  continue
       sum2 = sum2 * frac

       giii_hyd = sum1 + sum2
     end function giii_hyd


     function bfOpacity(freq, nAtom, thisAtom, pops, nstar, ne, temperature, ifreq) result(kappa)
       real(double) :: freq
       integer, optional :: ifreq
       type(MODELATOM) :: thisAtom(:)
       integer :: nAtom
       real(double) :: pops(:,:), nstar(:,:),test
       integer :: iAtom
       integer :: iTrans
       real(double) :: kappa, fac, ne, temperature
       integer :: ilower, j
       logical, save :: firstTime = .true.
       test = nstar(1, 1)
       kappa = 0.d0
       do iAtom = 1, nAtom
          do j = 1, thisAtom(iAtom)%nRBFtrans
             iTrans = thisAtom(iAtom)%indexRBFtrans(j)
          if (thisAtom(iatom)%transType(itrans) /= "RBF") then
             write(*,*) "transtype bug in bfemissivity"
             stop
          endif
             ilower = thisAtom(iAtom)%iLower(iTrans)
             if (ilower < 10) then
                fac = exp(-hCgs*freq / (kerg * temperature))
                if (present(iFreq)) then
                   kappa = kappa + quickPhotoCrossSection(thisAtom(iAtom), j, iFreq) * &
                        (pops(iAtom, ilower))! -nstar(iatom,ilower)*fac)

                else
                   kappa = kappa + photoCrossSection(thisAtom(iAtom),  ilower, freq) * &
                        (pops(iAtom, ilower))! - nstar(iatom,ilower)*fac)

                endif
             endif
          enddo
       enddo

       kappa = kappa +  ne * pops(1,thisAtom(1)%nlevels) * &
            alpkk_hyd(freq,dble(temperature)) * (1.d0-exp(-(hcgs*freq)/(kerg*temperature)))
       kappa = kappa + ne * sigmae


       if (kappa < 0.d0) then
          if (firstTime) then
             write(*,*) "negative bfopacity!!!!"
             firstTime = .false.
          endif
          kappa = 0.d0
       endif

     end function bfOpacity

  function bfEmissivity(freq, nAtom, thisAtom, pops, nStar, temperature, ne,  ifreq) result(eta)


! bound-free and free-free continuum emissivity

    type(MODELATOM) :: thisAtom(:)
    integer :: nAtom
    real(double) :: nStar(:,:)
    real(double) :: freq
    integer, optional :: iFreq
    integer :: iTrans
    real(double) :: temperature
    real(double) :: eta, fac
    integer :: iAtom, i, iLower, iUpper
    real(double) :: ne
    real(double) :: thresh, photonEnergy, expFac
    logical, save :: firsttime = .true.
    real(double) :: pops(:,:)

    eta=0.d0

!    bound-free

    fac = (2.d0 * hCgs * freq**3)/(cSpeed**2)
    expFac = exp(-(hcgs*freq)/(kerg*temperature))

    do iAtom = 1, nAtom

       do  i = 1, thisAtom(iAtom)%nRBFtrans
          iTrans = thisAtom(iAtom)%indexRBFtrans(i)
          iLower = thisAtom(iAtom)%iLower(iTrans)
          iUpper = thisAtom(iAtom)%iUpper(iTrans)
          if (thisAtom(iatom)%transType(itrans) /= "RBF") then
             write(*,*) "transtype bug in bfemissivity"
             stop
          endif
          if (iLower < 10) then
             thresh=(thisAtom(iAtom)%iPot - thisAtom(iAtom)%energy(iLower))
             photonEnergy = freq * hCgs * ergtoEv
             if (photonEnergy.ge.thresh) then
                if (present(ifreq)) then
                   eta = eta + fac * nStar(iAtom,iLower) * quickPhotoCrossSection(thisAtom(iAtom), i, iFreq) * expFac
                else
                   eta = eta + fac * nStar(iAtom,iLower) * photoCrossSection(thisAtom(iAtom), iLower, freq) * expFac
                endif
             endif
          endif
       enddo
    enddo


    eta = eta + fac *  ne *  pops(1,thisAtom(1)%nlevels) * alpkk_hyd(freq,dble(temperature)) * expFac


    if (eta < 0.d0) then
       if (firstTime)  then
          write(*,*) "negative bf emissivity"
          firstTime = .false.
       endif
       eta = 0.d0
    endif


  end function bfEmissivity


  subroutine createContFreqArray(nFreq, freq, nAtom, thisAtom, nsource, source)
    integer :: nAtom
    type(MODELATOM) :: thisAtom(:)
    integer :: nsource
    type(SOURCETYPE) :: source(:)
    integer :: nfreq
    real(double) :: freq(:)
    integer :: iAtom, iLevel
    integer :: i
    real(double) :: nuThresh, lamMin, lamMax
    real(double) :: nuStart, nuEnd
    integer :: nEven
    integer :: iTrans

    nFreq = 0

    nEven = 10

    nuStart = 25.d0 * evtoerg/hcgs
    nuEnd = 10.d0 * 1000.d0 * evtoerg/hcgs

    do i = 1, nEven
       nFreq = nFreq + 1
       freq(nFreq) = log10(nuStart) + (dble(i-1)/dble(nEven-1))*(log10(nuEnd)-log10(nuStart))
    enddo

    nuStart = cSpeed / (8.d5 * 1.d-8)
    nuEnd = cSpeed / (5.d0 * 1.d-8)

    nEven = 100


    do i = 1, nEven
       nFreq = nFreq + 1
       freq(nFreq) = log10(nuStart) + (dble(i-1)/dble(nEven-1))*(log10(nuEnd)-log10(nuStart))
    enddo

    freq(1:nFreq) = 10.d0**freq(1:nFreq)


! bound-free edges

    do iAtom = 1, nAtom

       do iLevel = 1, min(5,thisAtom(iAtom)%nLevels-1)
          nuThresh = (thisAtom(iAtom)%iPot - thisAtom(iAtom)%energy(iLevel))*evtoerg/hCgs
          nfreq = nfreq + 1
          freq(nfreq) = nuthresh * 0.99d0
          nfreq = nfreq + 1
          freq(nfreq) = nuthresh * 1.01d0
       enddo
    enddo
    lamMin = 1.d30
    lamMax = -1.d30
    do i = 1, nSource
       if (source(i)%spectrum%lambda(1) < lamMin) lamMin = source(i)%spectrum%lambda(1)
       if (source(i)%spectrum%lambda(source(i)%spectrum%nlambda) > lamMax) &
            lamMax = source(i)%spectrum%lambda(source(i)%spectrum%nlambda)
    enddo


    do iAtom = 1, nAtom
       do iTrans = 1, thisAtom(iAtom)%nTrans
          if (thisAtom(iAtom)%transType(itrans) == "RBB") then
             nfreq  = nFreq + 1
             freq(nFreq) = thisAtom(iAtom)%transFreq(iTrans)
          endif
       enddo
    enddo

    nfreq = nfreq + 1
    freq(nfreq) = cspeed/ (lamMin * 1.d-8)
    nfreq = nfreq + 1
    freq(nfreq) = cspeed/ (lamMax * 1.d-8)
    call sort(nFreq, Freq)


    if (writeoutput) Write(*,*) "Number of frequency points in continuum: ",nfreq
  end subroutine createContFreqArray

  subroutine stripAtomLevels(thisAtom, maxBoundLevels)
    type(MODELATOM) :: thisAtom
    integer :: maxBoundLevels
    integer :: iOldContinuum, iNewContinuum, iTrans

    if (maxBoundLevels > thisAtom%nLevels-1) then
       call writeFatal("maxBoundLevels exceeds bound levels in atom")
    endif

    iOldContinuum = thisAtom%nLevels
    iNewContinuum = maxBoundLevels + 1
    thisAtom%nLevels = maxBoundLevels+1

    thisAtom%level(maxBoundLevels+1) = thisAtom%level(iOldContinuum)
    thisAtom%g(maxBoundLevels+1) = thisAtom%g(iOldContinuum)
    thisAtom%energy(maxBoundLevels+1) = thisAtom%energy(iOldContinuum)

    do iTrans = 1, thisAtom%nTrans


       if (thisAtom%transType(iTrans)(2:3) == "BB") then
          if ((thisAtom%iLower(iTrans) > maxBoundLevels).or. &
               (thisAtom%iUpper(iTrans) > maxBoundLevels)) then
             thisAtom%transType(iTrans) = "OFF"
          endif
       endif


       if (thisAtom%transType(iTrans)(2:3) == "BF") then
          if (thisAtom%ilower(iTrans) <= maxBoundLevels) then
             thisAtom%iUpper(iTrans) = iNewContinuum
          else
             thisAtom%transType(iTrans) = "OFF"
          endif
       endif
!       write(*,*) itrans,thisAtom%transType(iTrans),thisAtom%iLower(iTrans), thisAtom%iupper(iTrans)
    enddo
  end subroutine stripAtomLevels

  function HeIICollisionalIonRates(i, iPot, temp) result (xSec)
! page 922 of Klein and Castor (1978) ApJ 220 902
    integer :: i
    real(double) :: iPot, temp
    real(double) :: xSec
    real(double) :: ai(12) = (/-4.59d-6, 2.20d-4, 1.18d-4, 1.76d-4, 2.70d-4, 3.95d-4, 5.51d-4, &
         7.38d-4, 9.55d-4, 1.20d-3, 1.48d-3, 1.79d-3/)
    real(double) :: bi(12) = (/ -1002190.d0, 2611480.d0, 262714.d0, 118656.d0, 72355.8d0, 50622.2d0, 38299.9d0, &
         30461.2d0, 25074.7d0, 21164.8d0, 18209.2d0, 15904.0d0/)
    xSec = (ai(i) / (temp + bi(i))) *sqrt(temp)*exp(-iPot/(kev*temp))
  end function HeIICollisionalIonRates

  !
  ! Bound-free (photoionization X-section from n-th energy level) for a
  ! photon with its frequency nu.
  !
  ! See sub_phot_gen.f in CMFGEN source files distrubution
  !
  ! The data tabulation is from Hillier.. Don't know where the original
  ! data is from.... Need to ask him
  ! Units in cm^2.
  function annu_He_I(n,nu) result(out)
    implicit none
    real(double) :: out
    integer, intent(in)      :: n   ! [-] energy level
    real(double), intent(in) :: nu  ! [Hz]  frequency of photon
    !
    integer, parameter :: nmax= 19  ! max level in the data file.. skipping the SuperLevels
    integer, parameter :: npar=8    ! # of parameter in the fitting function
    !
    real(double),save :: a(nmax,npar)     ! fitting parameters
    logical, save :: first_time = .true.
    !
    real(double) :: alpha
    character(LEN=80) :: msg
    !
    real(double), parameter :: E_1k=198310.7600d0 !  [Ry] ionization energy from ground
    ! Conversion factor for unit change from [Ry] to [Hz]
    real(double), parameter :: conv_fac = (13.6056923/109737.316 * 1.60217646d-12/6.626205d-27)
    ! Conversion factor for unit change from [Ry] to [eV]
!    real(double), parameter :: Ry2eV = 13.6056923/109737.316  !
    !

    ! Energy levels with respect the the ground state... in Hillier's data
    ! Units are in [Ry]
    real(double), parameter :: E_i(nmax) =  &
         (/ 0.00000000, 159850.32000000, 166271.70000000, 169081.26000000, 171129.15000000, &
         183231.08000000, 184859.06000000, 185559.01000000, 186095.90000000, 186099.22000000, &
         186203.62000000, 190292.46000000, 190934.50000000, 191211.42000000, 191438.83000000, &
         191440.71000000, 191446.61000000, 191447.24000000, 191486.95000000 /)

    real(double) :: Eion_from_n, U, X, T1, EDGE

!    ! JUST FOR TESTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    out = annu(n,nu)   !hydrogen!!
!    return
    !



    if (first_time) then
       first_time =.false.
       !
       ! Assiging the data (just copying from heiphot_a4.dat of Hillier's CMFGEN data file)
       !
       ! For level 1 ==> 1s1s1Se
       a(1,1:npar) = (/ 8.734D-01, -1.545D+00, -1.093D+00, 5.918D-01, 0.3262d0, 0.695319d0, -1.29300d0, 0.384251d0 /)
       ! For level 2 ==> 1s2s3Se
       a(2,1:npar) = (/ 7.377D-01, -9.327D-01, -1.466D+00, 6.891D-01, 0.7228d0, 0.901802d0, -1.85905d0, 0.854216d0 /)
       ! For level 3 ==> 1s2s1Se
       a(3,1:npar) = (/ 9.771D-01, -1.567D+00, -4.739D-01, -1.302D-01,0.6135d0, 1.13101d0,  -2.15771d0, 0.989081d0 /)
       ! For level 4 ==> 1s2p3Po
       a(4,1:npar) = (/ 1.204D+00, -2.809D+00, -3.094D-01, 1.100D-01, 0.9302d0, 1.16294d0,  -2.95758d0, 1.06043d0  /)
       ! For level 5 ==> 1s2p1Po
       a(5,1:npar) = (/ 1.129D+00, -3.149D+00, -1.910D-01, -5.244D-01,0.7360d0, 1.79027d0,  -4.47151d0, 0.986915d0 /)
       ! For level 6 ==> 1s3s3Se
       a(6,1:npar) = (/ 9.031D-01, -1.157D+00, -7.151D-01, 1.832D-01, 1.076d0,  1.25389d0,  -2.04057d0, 1.25914d0  /)
       ! For level 7 ==> 1s3s1Se
       a(7,1:npar) = (/ 1.174D+00, -1.638D+00, -2.831D-01, -3.281D-02, 0.9233d0, 1.36313d0, -2.13263d0, 1.34758d0  /)
       ! For level 8 ==> 1s3p3Po
       a(8,1:npar) = (/ 1.455D+00, -2.254D+00, -4.795D-01,  6.872D-02, 1.1440d0, 1.86467d0, -3.07110d0, 1.39913d0  /)
       ! For level 9 ==> 1s3d3De
       a(9,1:npar) = (/ 1.267D+00, -3.417D+00, -5.038D-01, -1.797D-02, 0.9585d0, 1.74181d0, -4.41218d0, 1.36716d0  /)
       ! For level 10 ==> 1s3d1De
       a(10,1:npar)= (/ 1.258D+00, -3.442D+00, -4.731D-01, -9.522D-02, 0.9586d0, 1.83599d0, -4.85608d0, 1.38737d0  /)
       ! For level 11 ==> 1s3p1Po
       a(11,1:npar)= (/ 1.431D+00, -2.511D+00, -3.710D-01, -1.933D-01, 1.0410d0, 2.23543d0, -4.3993d0,  1.35996d0  /)
       ! For level 12 ==> 1s4s3Se
       a(12,1:npar)= (/ 1.031D+00, -1.313D+00, -4.517D-01,  9.207D-02, 1.206d0,  1.39033d0, -2.02189d0, 1.52406d0  /)
       ! For level 13 ==> 1s4s1Se
       a(13,1:npar)= (/ 1.324D+00, -1.692D+00, -2.916D-01,  9.027D-02, 0.8438d0, 1.51684d0, -2.10272d0, 1.60093d0  /)
       ! For level 14 ==> 1s4p3Po
       a(14,1:npar)= (/ 1.619D+00, -2.109D+00, -3.357D-01, -2.532D-02, 1.0280d0, 2.02110d0, -2.87157d0, 1.64584d0  /)
       ! For level 15 ==>1s4d3De
       a(15,1:npar)= (/ 1.565D+00, -2.781D+00, -6.497D-01, -5.979D-03, 1.0410d0, 2.25756d0, -4.12940d0, 1.60889d0  /)
       ! For level 16 ==>1s4d1De
       a(16,1:npar)= (/ 1.553D+00, -2.781D+00, -6.841D-01, -4.083D-03, 1.1870d0, 2.50403d0, -4.40022d0, 1.63108d0  /)
       ! For level 17 ==>1s4f3Fo  ===> Modefied Seaton Fit for this level
       a(17,1:npar)= (/ 4.000D+00,  3.000D+00,  3.000D+00,  0.000D+00, 0.0000d0, 0.0000d0,   0.0000d0,  0.0000d0   /)
       ! For level 18 ==>1s4f1Fo  ===> Modefied Seaton Fit for this level
       a(18,1:npar)= (/ 4.000D+00,  3.000D+00,  3.000D+00,  0.000D+00, 0.0000d0, 0.0000d0,   0.0000d0,  0.0000d0   /)
       ! For level 19 ==>1s4p1Po
       a(19,1:npar)= (/ 1.620D+00, -2.303D+00, -3.045D-01, -1.391D-01, 1.2720d0, 2.63942d0, -3.71668d0, 1.60525d0  /)
    end if

    ! Checking

    if (nu <= 0.0d0) then
       write(msg,*) "Error:: nu <= 0.0d0 in annu_He_I!!!"
       call writeInfo(msg)
       write(msg,*) "   nu = ", nu
       call writefatal(msg)
       alpha = sqrt(nu)
       stop
    end if


    !
    ! Computes the cross section...

    if (n == 17 .or. n==18) then
       ! This is TYPE =3 in Hillier
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       out = 1.0d-20  ! just returning for small number for now
    elseif (n > nmax) then
       ! No data exists
       out = 1.0d-20  ! just returning for small number for now
    elseif (n<1) then
       ! out of range
       write(msg,*) "Error:: n is negative in photoion_xsect::annu_He_I."
       call writeinfo(msg)
       write(msg,*) "n, nmax = ", n, nmax
       call writefatal(msg)
       stop
    else
       !  This is TYPE =6 in Hillier
       !
       ! This type is for the Hummer fits to the Opacity cross sections of HeI.
       ! See HEI_PHOT_OPAC.
       !
       ! We issue the SPACING command to ensure that rounding error does not
       ! cause X to be < 0 at the bound-free edge.
       !
       ! U+SPACING(U) is the smallest number different from U (and lager).
       !

       ! Ionization energy from level n
       Eion_from_n = E_1k - E_i(n)    ! [Ry]

       EDGE=  Eion_from_n*conv_fac ! edge freqenecy  [Hz]


       if (EDGE <= 0.0d0) then
          write(msg,*) "Error:: EDGE <= 0.0d0 in annu_He_I!!!"
          call writeInfo(msg)
          write(msg,*) "   EDGE, n, nu = ",  EDGE, n, nu
          call writefatal(msg)
          stop
       end if

       U= nu/EDGE  ! dimensionless

       !
       X=LOG10(U+3.0D0*SPACING(U)) ! slightly larger than U
!       if (n == 1) then
!          print *, "n            = ", n
!          print *, "EDGE [in eV] = ", Eion_from_n*Ry2eV
!          print *, "EDGE [in Hz] = ", Eion_from_n*conv_fac
!          print *, "U [-]        = " , U
!          print *, "X [-]        = " , X
!          print *, "a(n,5)       = " , a(n,5)
!       end if

       IF(X .GE. 0)THEN
          IF(X .LT. a(n,5))THEN
             T1=((a(n,4)*X+a(n,3))*X + a(n,2))*X + a(n,1)
          ELSE
             T1=a(n,6)+a(n,7)*X    ! This is as in Hillier (sub_photo_gen.f
          END IF
          alpha = 10.0D0**(T1)  ! this should be in megabarn

          out = alpha*megabarn  ! should be in cm^2 here now.

       ELSE ! this is below the edge
          out = 1.0d-30
       END IF

    end if

   end function annu_He_I



  !
  !
  ! Bound-Free abs. coeff. for Helium (He II)
  !
  real(double) function annu_He_II(n,nu)
    !
    ! this subroutine calculates the photoionization x-section
    ! for He II from the n-th level for a given freq photon.
    !
    !  See P.132 of the book entitled:
    !  "The observation and analysis of stellar photosphere" by D. F. Gray
    !
    ! We use the fact that He II is hydrogenic ion and use propertionality
    ! from the proton number (Z).  (see Mihalas's book)
    !
    integer, intent(in)              :: n     ! the level
    real(double), intent(in):: nu    ! the photon frequency
    !
    !
    real(double), parameter :: Z = 2.0d0 ! proton number for He

    ! Using b-f x-section of hyrdogen and rescale..


    annu_He_II = annu_hgi(n,z,nu)   ! using a fucntion in this module


  end function annu_He_II



  !
  !
  ! Bound-Free abs. coeff. hydrogenic ion with nuclear charge z
  !
  real(double) function annu_hgi(n,z,nu)
    !
    ! this subroutine calculates the photoionization x-section
    ! for hydrogenenic ion from the n-th level for a given freq photon.
    !
    !  See P.132 of the book entitled:
    !  "The observation and analysis of stellar photosphere" by D. F. Gray
    !
    ! Note: Mihals's Stellar Atomosphere book on P105 states
    !       1. En (energy level) scales as Z^2
    !       2. B-F cross section scales as Z^4
    !       3. F-F cross section scales as Z^2
    !
    integer, intent(in)     :: n      ! the level
    real(double), intent(in):: z      ! Nuclear charge
    real(double), intent(in):: nu     ! the photon frequency
    real(double)            :: lam, E ! the photon wavelength
    real(double)            :: E_ion  ! ionization energy
    real(double)            :: En     ! Energy level
    real(double) :: z2, annu_hyd
    !
    ! Enegry levels of H I
    real(double), parameter :: En_H_I(30) =   &
       & (/ 0.0d0, 10.19882475d0, 12.087496d0, 12.74853094d0, 13.05449568d0,  &
       &    13.22069875d0, 13.32091396d0, 13.38595748d0, 13.43055111d0, 13.46244867d0,  &
       &    13.48604926d0, 13.50399944d0, 13.5179689d0, 13.52905324d0,  13.53799552d0,  &
       &    13.54531412d0, 13.5513796d0,  13.55646253d0, 13.56076421d0, 13.56443692d0,  &
       &    13.56759755d0, 13.57033706d0, 13.57272708d0, 13.57482461d0, 13.57667551d0,  &
       &    13.57831697d0, 13.57977946d0, 13.58108806d0, 13.58226364d0, 13.58332363d0   /)

    lam=cSpeed/nu
    lam=lam * 1.d8

    z2 = z*z
    E_ion = hydE0eV*z2
    En    = En_H_I(n)*z2

    E = hCgs * nu * ergToEv  ! energy of photon

    ! Check if the energy is the edge energy/frequency
    if (E > (E_ion - En)) then

!       annu_hyd = 1.044d-26*gii(n,1.0d0,lam)*(lam**3)/dble(n)**5

       annu_hyd = 1.044d-26*gii(n,z,lam)*(lam**3)/dble(n)**5

       annu_hgi = annu_hyd*z2*z2  ! Z^4 factor
    else
       annu_hgi = 1.d-30
    endif


  end function annu_hgi



  !
  !
  !
  real(double) pure function gii(n,z,wl)
    !
    ! returns bound-free gaunt factors (nicked from idh)
    !
    ! the original code had real/integer divisions, I'm asuming these are
    !   intentional and leaving them in. nhs.
    integer, intent(in)   :: n
    real(double), intent(in) ::  z,wl
    real(double) ::  efree,sum,ag,alam
    integer               :: i

    real(double) :: nDouble
    ! making coeff  a PARAMETER may cause problems with XL Fortran
    real(double),dimension(6), parameter :: coeff = &
          (/-0.338276d0, 0.808398d0, -0.59941d0, 0.104292d0, &
              -6.61998d-3, 1.15609d-4 /)
    nDouble = real(n,kind=double)

    efree=(911.76e0_db/wl - 1.e0_db/(nDouble*nDouble)) / (z*z)

    if (efree.le.(1.e0_db+2.e0_db/n)) then
       gii = giia(nDouble,z,wl)
       return
    elseif (efree.gt.10.e0_db) then
       sum=coeff(1)
       ag = 1.e0_db
       efree = log10(efree)
       do i = 2, 6
          ag = ag * efree
          sum = sum + coeff(i) * ag
       enddo
       gii = 10.e0_db**sum
       return
    else
       alam = 911.76e0_db / (z*z*(1.e0_db+2.e0_db/n)+1.e0_db/(n*n))
       alam = giia(nDouble,z,alam)
       sum = log10(1.e0 + 2.e0/n)
       efree = log10(efree)
       sum = (efree-sum) / (1.e0_db-sum)
       gii = (1.e0-sum) * alam + 0.93e0_db*sum+0.2e0_db*sum*(1.e0_db-sum)
       return
    endif

  end function gii



  !
  ! Used in gii function above.
  !
  real pure function giia(nDouble,z,wl)
    real(double), intent(in) :: nDouble
    real(double), intent(in) :: z, wl

    real(double)             :: u, gii, term
    real(double), parameter  :: twoThirds = 2.e0_db/3.e0_db
    real(double), parameter  :: fourThirds = 2.0_db * twoThirds

    u = nDouble*nDouble*911.76e0_db / (wl*z*z)-1.e0_db
    gii = 1.e0_db + 0.1728e0_db * (u-1.e0_db) / ((nDouble*(u+1.e0_db))**twoThirds)
    term = 0.0496e0_db * (1.e0_db+u*(u+1.333e0_db))/(nDouble*(u+1.e0_db)**(fourThirds))
    if ((term/gii).le.0.25e0_db) then
       giia=real(gii-term)
       return
    else
       giia = 1.0_db
       return
    endif

  end function giia


end module modelatom_mod
#endif
