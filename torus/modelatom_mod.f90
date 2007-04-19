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
     integer :: nLevels
     character(len=10), pointer :: level(:)
     real(double), pointer :: energy(:)  ! erg
     real(double), pointer :: g(:)
     real(double) :: iPot
     integer :: nPhotoFreq
     real(double), pointer :: photoFreq(:)
     real(double), pointer :: photoSec(:,:)
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
    allocate(thisAtom%inDetailedBalance(1:thisAtom%nTrans))

    thisAtom%inDetailedBalance = .false.

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


  real(double) function photoCrossSection(thisAtom, iTrans, iLevel, nu)
    type(MODELATOM) :: thisAtom
    integer :: iTrans
    integer :: iLevel
    real(double) :: nu, e
    character(len=2) :: shell
    integer :: is
    real(double) :: photocrosssection2
    real :: x
    real(double) :: photonEnergy, lamMicrons, logKappa
    real(double) :: sigma0, nuThresh, alpha, s, gIIx, gIIy, gIIz
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
                photonEnergy = nu * hCgs * ergtoEv
                if (photonEnergy > (thisAtom%iPot - thisAtom%energy(iLevel))) then
                   select case(iLevel)  ! from Gingerich 1964
                     case(1)
                        logKappa = 14.47d0 - 2.d0 * log10(nu)
                        photoCrossSection = (10.d0**logKappa)
                        call phfit2(2,2,1,real(nu*hCgs*ergtoev), x)
!                        write(*,*) "xsec", x * 1.d-10,photocrosssection
                        photocrosssection = x * 1.d-10
                     case(2)
                        logKappa = -17.387d0 + 0.679*log10(nu*1.d-14) - 0.727d0*log10(nu*1.d-14)
                        photoCrossSection = (10.d0**logKappa)
                     case(3)
                        logKappa = 11.65d0 - 1.91d0*log10(nu)
                        photoCrossSection = (10.d0**logKappa)
                     case(4)
                        logKappa = 26.57d0 - 2.9d0*log10(nu)
                        photoCrossSection = (10.d0**logKappa)
                     case(5)
                        logKappa = 35.31d0 - 3.5d0*log10(nu)
                        photoCrossSection = (10.d0**logKappa)
                     case(6)
                        logKappa = 31.06d0 - 3.3d0*log10(nu)
                        photoCrossSection = (10.d0**logKappa)
                     case(7)
                        logKappa = 35.48d0 - 3.6d0*log10(nu)
                        photoCrossSection = (10.d0**logKappa)
                     case DEFAULT
                        lamMicrons = (cspeed/nu)/microntocm
                        photoCrossSection = (1.0499d-14 * lamMicrons**3) / dble(iLevel**5) 
                   end select
                endif
             case(1)
                select case(thisAtom%equation(iTrans))
                   case(1)
                      nuThresh = (thisAtom%iPot-thisAtom%energy(thisAtom%iLower(iTrans)))*evtoerg/hCgs
                      if (nu > nuThresh) then
                         sigma0 = thisAtom%params(iTrans,1)
                         alpha = thisAtom%params(iTrans,2)
                         s = thisAtom%params(iTrans,3)
                         photoCrossSection = sigma0 * (nuThresh/nu)**s * (alpha + (1.d0-alpha)*(nuThresh/nu))
                      endif
                   case(2)
                      nuThresh = (thisAtom%iPot-thisAtom%energy(thisAtom%iLower(iTrans)))*evtoerg/hCgs
                      if (nu > nuThresh) then
                         sigma0 = thisAtom%params(iTrans,1)
                         alpha = thisAtom%params(iTrans,2)
                         s = thisAtom%params(iTrans,3)
                         gIIx = thisAtom%params(iTrans,4)
                         gIIy = thisAtom%params(iTrans,5)
                         gIIz = thisAtom%params(iTrans,6)
                         photoCrossSection = sigma0 * (nuThresh/nu)**s * (alpha + (1.d0-alpha)*(nuThresh/nu))!*gauntII(gIIx, GIIy, gIIz)
                      endif
                  end select
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


  real(double) function quickPhotoCrossSection(thisAtom, iRBF, iFreq)
    type(MODELATOM) :: thisAtom
    integer :: iRBF, iFreq
    quickPhotoCrossSection  = thisAtom%photoSec(iRBF, iFreq)
  end function quickPhotoCrossSection

    


   subroutine addCrossSectionstoAtom(thisAtom, nFreq, freq)
     type(MODELATOM) :: thisAtom
     integer :: nFreq
     real(double) :: freq(:)
     integer :: iFreq, i, j, iTrans, iLevel, iLineTrans

     if (associated(thisAtom%photoFreq)) deallocate(thisAtom%photoFreq)
     if (associated(thisAtom%photoSec)) deallocate(thisAtom%photoSec)

     thisAtom%nPhotoFreq = nFreq
     allocate(thisAtom%photoFreq(1:nfreq))
     allocate(thisAtom%photoSec(1:thisAtom%nRBFtrans, 1:nFreq))

     do j = 1, thisAtom%nRBFTrans
        iTrans = thisAtom%indexRBFtrans(j)
        iLevel = thisAtom%iLower(iTrans)
        do iFreq = 1, nFreq
           thisAtom%photoSec(j, iFreq) = photoCrossSection(thisAtom, iTrans, iLevel, freq(ifreq))
        enddo
     enddo
        


   end subroutine addCrossSectionstoAtom




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
    real :: tReal
    real, parameter ::  Ucoeff(5,5) =  reshape(source=  &
        (/0.30103e0, -0.00001e0, 0.00000e0,  0.00000e0, 0.00000e0, &
          0.00000e0,  0.00000e0, 0.00000e0,  0.00000e0, 0.00000e0, &
          0.30103e0,  0.00000e0, 0.00000e0,  0.00000e0, 0.00000e0, &
          0.07460e0, -0.75759e0, 2.58494e0, -3.53170e0, 1.65240e0, &
          0.34383e0, -0.41472e0, 1.01550e0,  0.31930e0, 0.00000e0  &
           /), shape=(/5,5/))
    real(double) :: N2, N1, N0, u0, u1, u2, N1overN0, N2overN1, pe, tot
       pe = ne * kerg * t
       tReal = t
       select case(thisAtom%nz)
       case(1)
          u0 = getUT(treal, uCoeff(1,:))
          u1 = 1.d0
          N1overN0 = ((-5040.d0/t)*thisAtom%iPot + 2.5d0*log10(t) + log10(u1/u0)-0.1762d0)
          N1overN0 = (10.d0**N1overN0)/pe
          
          !          n1overn0 = (2.d0*u1)/(ne*u0) * ((twoPi * melectron * kerg * t)/(hCgs**2))**1.5d0 * exp(-thisAtom%iPot/(kev*t))
          
          if (N1overN0 /= 0.d0) then
             ratio = (thisAtom%g(level)*exp(-thisAtom%energy(level)/(kev * t))) / u0 / N1overn0
          else
             ratio = 1.d10
          endif
       case(2)
          select case(thisAtom%charge)
          case(0)
             u0 = getUT(treal, uCoeff(2,:))
             u1 = getUT(treal, uCoeff(3,:))
             N1overN0 = ((-5040.d0/t)*24.59d0 + 2.5d0*log10(t) + log10(u1/u0)-0.1762d0)
             N1overN0 = (10.d0**N1overN0)/pe
             ratio = (thisAtom%g(level)*exp(-thisAtom%energy(level)/(kev * t))) / u0 / N1overn0
          case(1)
             u1 = getUT(treal, uCoeff(3,:))
             u2 = 1.d0
             N2overN1 = ((-5040.d0/t)*54.42d0 + 2.5d0*log10(t) + log10(u2/u1)-0.1762d0)
             N2overN1 = (10.d0**N2overN1)/pe
             ratio = (thisAtom%g(level)*exp(-thisAtom%energy(level)/(kev * t))) / u1 / N2overn1
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
          

   function getUT(t, coeff) result (ut)
     real :: t, coeff(:),  ut
     real :: logphi

     logphi = 5040.0/t

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
              thisAtom(iAtom)%indexRBBtrans(itrans) = thisAtom(iAtom)%nRBBTrans
           endif
           if (thisAtom(iAtom)%transType(iTrans) == "RBF") then
              thisAtom(iAtom)%nRBFTrans = thisAtom(iAtom)%nRBFTrans + 1
              thisAtom(iAtom)%indexRBFtrans(thisAtom(iAtom)%nRBFTrans) = iTrans
           endif
        enddo
        write(*,*) "nRBBtrans",thisAtom(iAtom)%nRBBTrans, &
             thisAtom(iatom)%indexRBBTrans(1:thisAtom(iatom)%nrbbtrans)
        write(*,*) "nRBFtrans",thisAtom(iAtom)%nRBFTrans, &
             thisAtom(iatom)%indexRBFTrans(1:thisAtom(iatom)%nrbftrans)
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
          nStar = BoltzSahaGeneral(thisAtom, 1, j, Ne, temperature) * pops(nLevels)
          eta = eta + nStar * annu_hyd(j,freq) * exp(-(hcgs*freq)/(kerg*temperature))
       endif
    enddo
    
! free-free

    eta = eta + (ne**2) * alpkk_hyd(freq,temperature) * exp(-(hcgs*freq)/(kerg*temperature))
  
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
   50  continue
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


     function bfOpacity(freq, nAtom, thisAtom, pops, ne, nStar, temperature, ifreq) result(kappa)
       real(double) :: freq
       integer, optional :: ifreq
       integer :: nAtom
       type(MODELATOM) :: thisAtom(:)
       real(double) :: pops(:,:)
       integer :: iAtom
       integer :: iTrans
       real(double) :: kappa, fac, ne, temperature
       real(double) :: nStar(:,:)
       integer :: i, j

       kappa = 0.d0
       do iAtom = 1, nAtom
          do j = 1, thisAtom(iAtom)%nRBFtrans
             iTrans = thisAtom(iAtom)%indexRBFtrans(j)
             i = thisAtom(iAtom)%iLower(iTrans)
             if (i < 7) then
!                fac = exp(-hCgs*thisAtom(iAtom)%transFreq(iTrans) / (kerg * temperature))
                fac = exp(-hCgs*freq / (kerg * temperature))
                if (present(iFreq)) then
                   kappa = kappa + quickPhotoCrossSection(thisAtom(iAtom), j, iFreq) * &
                        (pops(iAtom, i) - nStar(iatom,i)*fac)
                else
                   kappa = kappa + photoCrossSection(thisAtom(iAtom), iTrans, i, freq) * &
                        (pops(iAtom, i) - nStar(iatom,i) *fac) 
                endif
             endif
          enddo
       enddo

       kappa = kappa + ne * sigmae

       if (kappa < 0.d0) kappa = 0.d0 

       
!write(*,*) "kappa",kappa

     end function bfOpacity

  function bfEmissivity(freq, nAtom, thisAtom, nStar, temperature, ne, Jnu, ifreq) result(eta)

! bound-free and free-free continuum emissivity

    type(MODELATOM) :: thisAtom(:)
    integer :: nAtom
    real(double) :: nStar(:,:)
    real(double) :: freq
    integer, optional :: iFreq
    integer :: iTrans
    integer :: nLevels
    real(double) :: temperature
    real(double) :: eta
    integer :: i, j, iAtom
    real(double) :: Jnu ! mean intensity
    real(double) :: ne
    real(double) :: thresh, photonEnergy, expFac

    eta=0.d0

!    bound-free

    expFac = exp(-(hcgs*freq)/(kerg*temperature)) * (2.d0 * hCgs * freq**3)/(cSpeed**2)

    do iAtom = 1, nAtom


       do  i = 1, thisAtom(iAtom)%nRBFtrans
          iTrans = thisAtom(iAtom)%indexRBFtrans(i)
          j = thisAtom(iAtom)%iLower(iTrans)
          if ( j < 7) then
             thresh=(thisAtom(iAtom)%iPot - thisAtom(iAtom)%energy(j))
             photonEnergy = freq * hCgs * ergtoEv
             if (photonEnergy.ge.thresh) then
                if (present(ifreq)) then
                   eta = eta + nStar(iAtom,j) * quickPhotoCrossSection(thisAtom(iAtom), i, iFreq) * expFac
                else
                   eta = eta + nStar(iAtom,j) * photoCrossSection(thisAtom(iAtom), iTrans, j, freq) * expFac
                endif
             endif
          endif
       enddo
    enddo
    eta = eta +  ne * sigmaE * Jnu ! coherent electron scattering

  end function bfEmissivity


  subroutine createContFreqArray(nFreq, freq, nAtom, thisAtom, nsource, source, maxFreq)
    integer :: nAtom
    type(MODELATOM) :: thisAtom(:)
    integer :: nsource
    type(SOURCETYPE) :: source(:)
    integer :: nfreq
    real(double) :: freq(:)
    integer :: iAtom, iLevel
    integer :: maxFreq, i
    real(double) :: nuThresh, lamMin, lamMax
    real(double) :: nuStart, nuEnd
    integer :: nEven
    integer :: iRBB, iTrans

    nuStart = cSpeed / (8.d5 * 1.d-8)
    nuEnd = cSpeed / (10.d0 * 1.d-8)

    nEven = 50


    do i = 1, nEven
       freq(i) = log10(nuStart) + (dble(i-1)/dble(nEven-1))*(log10(nuEnd)-log10(nuStart))
    enddo
    freq(1:nEven) = 10.d0**freq(1:nEven)
    nfreq = nEven


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


    Write(*,*) "Number of frequency points in continuum: ",nfreq
  end subroutine createContFreqArray

  subroutine stripAtomLevels(thisAtom, maxBoundLevels)
    type(MODELATOM) :: thisAtom
    integer :: maxBoundLevels
    integer :: i, iOldContinuum, iNewContinuum, iTrans
  
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

end module modelatom_mod
