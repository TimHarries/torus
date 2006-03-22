module density_mod

  !
  ! Generic function to return a density for a given point
  ! and geometry type.
  !
  ! If you add a new geometry, please don't forget to add the name
  ! in the list in 
  !
  
  use utils_mod
  use vector_mod
  use gridtype_mod
  use jets_mod
  use wr104_mod
  use ostar_mod
  use luc_cir3d_class
  use cmfgen_class

  implicit none
  
  public :: density,spiralWindDensity
  ! the specific definition of density functions should really be private, 
  ! but for now they are public... Better yet, they should be in their own module.
  ! See jets_mod for example. 
  
  private :: &
       print_geometry_list, &
       TTauriDensity
  
contains

  !
  !  For a given point (as a vector) with each component in
  !  10^10 cm and gridtype object, it will return the density in g/cm^3
  function density(r_vec, grid) RESULT(out)
    implicit none
    real(double) :: out
    type(octalVector), intent(in) :: r_vec
    type(gridtype), intent(in) :: grid


    select case (grid%geometry)
       
    case ("jets")
       ! using a routine in jets_mod.f90
       out = JetsDensity(r_vec, grid)/1000.d0 
       !                              ^^^^^^^
       !                          Converting kg/m^3 to g/cm^3

!       ! test function....
!       out = density_inverse_sq(r_vec)

    case ("ttauri")
       !  using a routine this module
       out = TTauriDensity(r_vec, grid) ! [g/cm^3]
!       out = TTauriDensityFuzzy(r_vec, grid) ! [g/cm^3]
       ! Double check the units

    case("testamr")
       out = testDensity(r_vec,grid)

    case("proto")
       out = protoDensity(r_vec,grid)

    case("wrshell")
       out = wrshellDensity(r_vec,grid)

    case("spiralwind")
       out = spiralWindDensity(r_vec, grid)

    case("benchmark")
       out = benchmarkDensity(r_vec, grid)

    case("melvin")
       out = melvinDensity(r_vec, grid)

    case("shakara","aksco")
       out = shakaraSunyaevDisc(r_vec, grid)

    case("warpeddisc")
       out = warpedDisc(r_vec, grid)

    case("ppdisk")
       out = ppdiskDensity(r_vec, grid)

    case("clumpydisc")
       out = clumpydisc(r_vec, grid)
    case("luc_cir3d")
       out = luc_cir3d_density(r_vec)  ! [g/cm^3]

    case("cmfgen")
       out = cmfgen_density(r_vec)  ! [g/cm^3]
       
    case default
       print *, 'Error:: Geometry option passed to [density_mod::density] '
       print *, '       (through grid) has not been implemented yet. Sorry.'
       print *, 'Your geometry = ', grid%geometry
       print *, 'Avilable gemetries are: '
       call print_geometry_list()
       print *, 'Terminating the program ...'
       stop 
    end select

  end function density
    


  !
  ! List the available geometry
  ! 
  subroutine print_geometry_list()    
    implicit none
    ! # of option available currently.
    integer, parameter :: n = 7
    character(LEN=30) :: name(n)
    integer :: i

    ! List here
    name(1) =  'jets'
    name(2) =  'wr104'
    name(3) =  'ttauri'
    name(4) =  'testamr'
    name(5) =  'spiralwind'
    name(6) =  'luc_cir3d'
    name(7) =  'cmfgen'
    
    do i = 1, n
       write(*, *)  '   ', i, '. ', name(i)
    end do
       
  end subroutine print_geometry_list




  LOGICAL PURE FUNCTION TTauriInFlow(point,grid,ignoreDisk)
    ! tests if a point lies within the T Tauri accretion flow  
    
    use input_variables, only: TTauriRinner, TTauriRouter, TTauriRstar, &
                               TTauriDiskHeight, TTauriMstar

    IMPLICIT NONE

    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid
    LOGICAL, OPTIONAL, INTENT(IN) :: ignoreDisk
    
    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec
    REAL              :: phi
    REAL              :: r, rM, theta, y
    

    starPosn = grid%starPos1
    pointVec = (point - starPosn) * 1.e10_oc
    r = modulus( pointVec ) 
    
    theta = ACOS(MIN(ABS(pointVec%z/r),0.998_oc))
    rM  = r / SIN(theta)**2
    y = SIN(theta)**2 

    ! test if the point lies within the star
    IF ( r < TTauriRstar ) THEN
      TTauriInFlow = .FALSE.
      
    ! test if the point lies too close to the disk
    ELSE IF ( ABS(pointVec%z) < TTauriDiskHeight) THEN
      IF (PRESENT(ignoreDisk)) THEN
        IF (ignoreDisk) THEN 
          TTauriInFlow = .TRUE.
        END IF 
      ELSE 
         TTauriInFlow = .FALSE.
      END IF 
   
    ! test if the point lies outside the accretion stream
    ELSE IF ((rM > TTauriRinner*0.99d0) .AND. (rM < TTauriRouter*1.01d0 )) THEN
      TTauriInFlow = .TRUE.

    ELSE
      TTauriInFlow = .FALSE.

    END IF

  END FUNCTION TTauriInFlow


  real function TTauriDensity(point,grid,ignoreDisk) result(rho) 
    ! given a position in the grid, we calculate the elapsed free-fall duration
    ! from the disc surface, and then use the mass accretion rate at the time
    ! when the material left the disc. 

    use input_variables, only: TTauriRinner, TTauriRouter, TTauriRstar, &
                               TTauriMstar
    use flowSpeedVariables
    
    type(GRIDTYPE), intent(in)    :: grid
    type(octalVector), intent(in) :: point
    logical, optional, intent(in) :: ignoreDisk
    real :: y

    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec

    real(oct) :: r, theta
    real :: TTauriMdotLocal

    starPosn = grid%starPos1
    pointVec = (point - starPosn) * 1.e10_oc
    r = modulus( pointVec ) 
    theta = acos( pointVec%z  / r )
    y = SIN(theta)**2 

    IF (TTauriInFlow(point,grid,ignoreDisk)) then 
    
      ! we call the accretion rate function to determine the appropriate value.
      TTauriMdotLocal = TTauriVariableMdot(point,grid,ignoreDisk)

      rho = (TTauriMdotLocal * TTauriRstar) / (4.0 * pi * &
            (TTauriRStar/TTauriRinner - TTauriRstar/TTauriRouter)) &
            * (r**(-5.0/2.0) / SQRT( 2.0 * bigG * TTauriMstar )) &
            * (SQRT( 4.0 - 3.0*y) / SQRT( 1.0 - y)) 
 
      rho = max(rho,1.e-25)
    ELSE
      rho = 1.e-25
      RETURN
    END IF
    
  end function TTauriDensity



  pure function TTauriFlowSpeedFunc(bigR)
    ! returns the component of the T Tauri flow speed that is in the disc plane

    use flowSpeedVariables
    use input_variables, only: TTauriMstar, TTauriRstar
    
    real, intent(in), dimension(:) :: bigR 
    real, dimension(size(bigR)) :: radius
    real, dimension(size(bigR)) :: TTauriFlowSpeedFunc
    real, dimension(size(bigR)) :: y

    radius = (flowPointRm * bigR**2)**(1./3.) 
    y = radius / flowPointRm

    TTauriFlowSpeedFunc = SQRT((2.0 * bigG * TTauriMstar / TTauriRstar) * &
                          (TTauriRstar / radius(:) - TTauriRstar/flowPointRm)) * &
                          (3.0 * SQRT(y) * SQRT(1.0-y) / SQRT(4.0 - (3.0*y))) 
    TTauriFlowSpeedFunc = 1. / TTauriFlowSpeedFunc 
    
  end function TTauriFlowSpeedFunc

  real function TTauriTimeInFlow(rM,rStart)
    ! lookup table
      
    use flowSpeedVariables
    use input_variables, only: TTauriRinner, TTauriRouter, TTauriRstar, mDotType
    use utils_mod, only: qsimp, locate

    real, intent(in) :: rM, rStart
    real, dimension(:,:), allocatable, save :: timeTable
    real, dimension(:), allocatable, save   :: rMindex
    real, dimension(:), allocatable, save   :: bigRindex
    integer, parameter :: numRm = 80
    integer, parameter :: numBigR = 700
    integer :: iRm, iBigR
    logical, save :: warned_once_already = .false.

    if (mdottype == "constantcurtains" .or. mdottype == "constant") then
       ! no need to create the timeTable !
       TTauriTimeInFlow = 1.0
    else

       if (.not. allocated(timeTable)) then
          print *, 'Creating lookup table in ''TTauriTimeInFlow'''
          
          allocate(timeTable(numBigR,numRm))
          allocate(rMindex(numRm))
          allocate(bigRindex(numbigR))
          
          do iRm = 1, numRm, 1
             rMindex(iRm) = (TTauriRinner*0.99) + (((TTauriRouter*1.01)-(TTauriRinner*0.99))/REAL(numRm-1) * REAL(iRm-1))
          end do
          
          do iBigR = 1, numBigR
             bigRindex(iBigR) = TTauriRstar*0.5 + ((TTauriRouter-TTauriRstar*0.5)/REAL(numBigR-1) * REAL(iBigR-1))
          end do
          
          do iRm = 1, numRm, 1
             flowPointRm = rmIndex(iRm)
             do iBigR = 1, numBigR, 1
                
                if (bigRindex(iBigR) < rmIndex(iRm)*0.9995) then 
                   timeTable(iBigR,iRm) = qsimp(TTauriFlowSpeedFunc,bigRindex(iBigR),rmIndex(iRm)*0.9995)
                else
                   timeTable(iBigR,iRm) = 0.0
                end if
                
             end do
          end do
          print *, 'Lookup table complete'

       end if

       call locate(rMindex,SIZE(rMindex),rM,iRm)
       call locate(bigRindex,SIZE(bigRindex),rStart,ibigR)
       if ((iRm == 0) .or. (iRm == numRm) .or. (iBigR == 0) .or. (ibigR == numBigR)) then
          if (.not. warned_once_already) then
             print *, '''Locate'' returned out of range value in TTauriTimeInFlow'
             print *, rMindex(1), rMindex(SIZE(rMindex)), rM, iRm
             print *, "=========WARNING *** WARNING *** WARNING ================"
             print *, "Program has set TTauriTimeFlow to 1.0, and continuing..."
             print *, "=========WARNING *** WARNING *** WARNING ================"
             warned_once_already = .true.
          end if
          TTauriTimeInFlow = 1.0
       else
          TTauriTimeInFlow = timeTable(iBigR,iRm)
       end if

    end if  ! (mdottype == "constantcurtains" .or. mdottype == "constant")
    
  end function TTauriTimeInFlow

  real function TTauriVariableMdot(point,grid,ignoreDisk)
    ! returns the mass accretion rate at a given point at the start of an
    !   accretion flow, at a certain time.
    ! the value is in units of grams per second, *IF* this value was constant
    !   across the entire flow base. see hartmann, hewett and calvet 1994.
    ! each different type of accretion uses values from the input file.
   
    use flowSpeedVariables
    use utils_mod
    use clump_mod
    use input_variables, only: MdotParameter1, MdotParameter2, &
                               MdotParameter3, MdotParameter4, &
                               MdotParameter5, MdotParameter6, &
                               MdotType, curtainsPhi1s, curtainsPhi1e, &
                               curtainsPhi2s, curtainsPhi2e, TTauriRinner,&
                               curtain_number, curtain_width, &
                               TTauriRouter, TTauriDiskHeight, TTauriRstar, &
                               phaseTime, nStartPhase, nPhase

    type(GRIDTYPE), intent(in)    :: grid
    type(octalVector), intent(in) :: point
    logical, optional, intent(in) :: ignoreDisk
    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec
    real :: r, rM, theta, y
    real :: Rstart
    !real :: Rend, thetaDisk
    real :: timeFlowStart, timeInFlow
    real :: azimuth    ! azimuth angle (radians)

    ! clump variables
    real    :: timeBeforeSimulation, totalTime
    integer :: arraySizeNeeded, iClump
    real    :: uniformRandom, timeToNextClump
    real    :: runningTime, meanTime
    ! the 'clumps' variable is now in clump_mod 
    real(oct) :: r_oct
    
    starPosn = grid%starPos1
    pointVec = (point - starPosn) * 1.e10_oc
    r_oct = modulus( pointVec ) 
    r = r_oct
    
    if (TTauriInFlow(point,grid,ignoreDisk)) then
      
      theta = acos( pointVec%z  / r )
      if (abs(modulo(theta,real(pi))) > 1.e-10 ) then 
        rM  = r / sin(theta)**2
      else
        rM = huge(rM)
      end if
      
      y = SIN(theta)**2 

      !radius = (rM-TTauriRinner) / (TTauriRouter-TTauriRinner)
      !radius = min(0.999,radius) ! to correct numerical errors
      !radius = max(0.001,radius)

      ! calculate the start and end points for the integration

      !rEnd = rM * 0.9999 
      rStart = sqrt(pointVec%x**2 + pointVec%y**2)

      flowPointRm = rM ! this gets passed to TTauriFlowSpeedFunc through a module
      !timeInFlow = qsimp(TTauriFlowSpeedFunc,rStart,rEnd)

      timeInFlow = TTauriTimeInFlow(rM,rStart)

      timeFlowStart = grid%timeNow - timeInFlow

      azimuth = ATAN2(pointVec%y,pointVec%x)

      select case (MdotType)
 
        case ('constant')
          TTauriVariableMdot = mdotparameter1
 
        case ('pulse1')
          ! after half an hour, step-change the accretion rate for a given time
          
          !  mDotParameter1 = normal accretion rate (Msol/year)
          !  mDotParameter2 = exceptional accretion rate (Msol/year)
          !  mDotParameter3 = time from start before exceptional accretion rate (s)
          !  mDotParameter4 = duration of exceptional accretion rate event (s)
          
          if (timeFlowStart < mDotParameter3) then
            TTauriVariableMdot = mdotparameter1
          else if (timeFlowStart < (mDotParameter3 + mDotParameter4)) then 
            TTauriVariableMdot = mdotparameter2
          else 
            TTauriVariableMdot = mdotparameter1
          end if
 
        case ('clumpy1')
          ! a steady stream with some some patches of increased accretion
          
          ! mdotparameter1 = steady 'background' level (Msol/year)
          ! mdotparameter2 = effective accretion level of each clump
          !   this is quantified as the total accretion rate *if* the 
          !   accretion stream was entirely at the 'clumpy' rate
          !   (Msol/year)
          ! mdotparameter3 = mean time between clumps forming (s)
          ! mdotparameter4 = duration of each clump (s)
          ! mdotparameter5 = angular size of each clump (deg)
         
          if (.not. allocated(clumps)) then 
          
            ! we first set up an array storing the times when *all* of the clumps
            !   (present and future) will be launched. we make this big enough 
            !   to last for all the phases of the simulation, and also we start
            !   the clumps well in the past. 

            ! overestimate how far back in time we have to go         
            timeBeforeSimulation = 200000. !(s)

            ! calculate the total time range over which we need to track clumps
            totalTime = timeBeforeSimulation + phaseTime * REAL(nStartPhase+nPhase) * 1.5

            ! estimate how many clumps we will have during this time
            meanTime = mDotParameter3
            arraySizeNeeded = INT(totalTime / meanTime + 1.)
            arraySizeNeeded = MAX(arraySizeNeeded,1)
            allocate(clumps(arraySizeNeeded))

            runningTime = -1.0 * timeBeforeSimulation

            do iClump = 1, size(clumps), 1 

              clumps(iClump)%startTime = runningTime
              clumps(iClump)%duration  = mDotParameter4
              clumps(iClump)%mDot  = mDotParameter2
              call random_number(uniformRandom)
              clumps(iClump)%azimuth = (uniformRandom * twoPi) - pi
              clumps(iClump)%angularSize = mDotParameter5 * degToRad 

              call random_number(uniformRandom)
              timeToNextClump = -meanTime * log(1-uniformRandom)

              runningTime = runningTime + timeToNextClump
            end do 

            print *, 'Using mDot model clumpy1'
            print *, 'Mean time between clumps (s): ',meanTime
            print *, 'Duration of each clump (s): ', mDotParameter4
            print *, 'Total number of clumps (an overestimate): ',SIZE(clumps)
            print *, 'Angular size (in azimuth) of a clump (deg): ', mDotParameter5

            do iClump = 1, size(clumps), 1
              print *, iClump, clumps(iClump)
            end do
                     
          end if

          ! now we actually use the clumps...
          
          where (((azimuth > clumps%azimuth-(clumps%angularSize/2.0)) .or.        & 
                  (azimuth-twoPi > clumps%azimuth-(clumps%angularSize/2.0))).and. &
                 ((azimuth < clumps%azimuth+(clumps%angularSize/2.0)) .or.        & 
                  (azimuth+twoPi < clumps%azimuth+(clumps%angularSize/2.0))).and. &
                  (timeFlowStart > clumps%startTime)                  .and.       &
                  (timeFlowStart < clumps%startTime + clumps%duration)    )
            clumps%activeFlag = .true.
          elsewhere 
            clumps%activeFlag = .false.
          end where

          if (any(clumps%activeFlag)) then
            TTauriVariableMdot = MAXVAL(clumps%mDot,MASK=clumps%activeFlag)
          else
            TTauriVariableMdot = mDotParameter1
          end if
          
        case ('linearincrease')
          ! mdotparameter1 gives value at t=0
          ! mdotparameter2 specifies rate of change (in Msol per year per second)
          TTauriVariableMdot = mdotparameter1 + mdotparameter2*timeFlowStart

        case ('constantcurtains')  
          if (azimuth < 0.) azimuth = azimuth + twoPi
          if (in_curtain(azimuth)) then ! using a function in this module
            TTauriVariableMdot = mdotparameter1
          else 
            TTauriVariableMdot = mDotParameter2 
          end if
          
        case ('semicurtains')  
          ! one curtain is only above the disc, and one is only below.
          
          if (azimuth < 0.) azimuth = azimuth + twoPi
          if (((azimuth > curtainsPhi1s).and.(azimuth < curtainsPhi1e)) .and.  &
              pointVec%z > 0.0) then 
            TTauriVariableMdot = mdotparameter1
          else if (((azimuth > curtainsPhi2s).and.(azimuth < curtainsPhi2e)) .and. &
              pointVec%z < 0.0) then 
            TTauriVariableMdot = mdotparameter1
          else 
            TTauriVariableMdot = mDotParameter2 
          end if
          
        case ('test')
          if (azimuth < 0.) azimuth = azimuth + twoPi
          azimuth = azimuth / twoPi
          azimuth = azimuth * 20 
          if (MOD(azimuth,2.0) > 1.0) then 
            IF (MOD(timeFlowStart,3600.) > 1700.) then
              TTauriVariableMdot = 0.
            ELSE 
              TTauriVariableMdot = mdotparameter1
            END IF
          else 
            IF (MOD(timeFlowStart,3600.) > 1700.) then
              TTauriVariableMdot = mdotparameter1
            ELSE 
              TTauriVariableMdot = 0.
            END IF
          end if

        case DEFAULT
          print *, 'Accretion type not recognized in TTauriVariableMdot'
          stop
      end select
      
      TTauriVariableMdot = TTauriVariableMdot * secsToYears * mSol
!print *, 'TTauriVariableMdot = ',TTauriVariableMdot
      
    else
      TTauriVariableMdot = 1.e-25 
    end if

  contains

    real function discGetYfunction(y)
    
      USE functionVariables
      
      real, intent(in) :: y
      
      discGetYfunction = 3.*sqrt(y)*sqrt(1.-y)/sqrt(4.-3.*y) - zValue
      
    end function discGetYfunction
    
  end function TTauriVariableMdot


  ! 
  ! 
  ! Function to check if a given phi (azimuth angle) is within one of curtain in 
  ! constant curtain. 
  function in_curtain(phi) RESULT(out)
    use input_variables, only: curtain_number, curtain_width
    implicit none
    logical :: out 
    !
    real, intent(in)  :: phi  ! should be in radians
    logical, save :: first_time = .true.
    real :: pi = 3.141592654
    ! this offset is need to ensure that there is  density 
    ! at phi=45 degrees since this plane is use to 2D calulations
    !  and then to map to 3D density structure. 
    real, save :: offset
    real,save :: gap  ! between curtains (radians)
    real, allocatable, save :: beg_c(:), end_c(:)  ! of curtains (in radians)
    integer :: i
    real :: phi_local

    if (first_time) then
       ! check the validity of input parameters
       ! (curtain_number, curtain_width).
       if (real(curtain_number)*curtain_width > 2.0*Pi) then
          write(*,*) "Error:: Too many curtains or the curtain size is too wide! [density_mod::in_curtain]."
          write(*,*) "curtain_number=", curtain_number
          write(*,*) "curtain_width=", curtain_width*180.0/pi, "  (degrees)"
          write(*,*) " ---> Adjust these parameters in your input file!"
          stop
       end if

       ALLOCATE(beg_c(curtain_number), end_c(curtain_number))
       offset = pi/4.0
       gap = (2.0*Pi - real(curtain_number)*curtain_width ) &            
            / real(curtain_number)

       ! set up the curtain ranges
       ! these arrays are saved for subseqence uses.
       beg_c(1) = offset - curtain_width/2.0
       end_c(1) = beg_c(1) + curtain_width
       do i = 2, curtain_number
          beg_c(i) = end_c(i-1) + gap
          end_c(i) = beg_c(i) + curtain_width
       end do

       ! add 2Pi if any of them are less than 0
       do i = 1, curtain_number
          if (beg_c(i) < 0.0) beg_c(i) = beg_c(i)+2*pi
          if (end_c(i) < 0.0) end_c(i) = end_c(i)+2*pi
       end do

       first_time=.false.
    end if

    ! Now check if phi is within range.
    phi_local = phi
    if (phi < 0.0) phi_local=phi+2.0*Pi ! for safty
    out = .false.
    do i =1, curtain_number
       if ( beg_c(i) < phi_local .and. phi_local < (beg_c(i) + curtain_width) ) then
          out = .true.
          goto 100
       end if
       if ( end_c(i) > phi_local .and. phi_local > (end_c(i) - curtain_width) ) then
          out = .true.
          goto 100
       end if
    end do
100 continue

    
  end function in_curtain

  !
  !
  !
!  REAL PURE FUNCTION TTauriDensityOld(point,grid,ignoreDisk) RESULT(rho)
!    ! calculates the density at a given point for a model of a T Tauri 
!    !   star with magnetospheric accretion
!    !   see Hartman, Hewett & Calvet 1994ApJ...426..669H 
!
!    use parameters_mod
!    use input_variables, only: TTauriMstar, TTauriRstar, TTauriRinner, TTauriRouter
!
!    IMPLICIT NONE
!
!    TYPE(octalVector), INTENT(IN) :: point
!    TYPE(gridtype), INTENT(IN)    :: grid
!    LOGICAL, OPTIONAL, INTENT(IN) :: ignoreDisk
!
!    TYPE(octalVector) :: starPosn
!    TYPE(octalVector) :: pointVec
!
!    REAL :: r, rM, theta, y, ang
!
!    starPosn = grid%starPos1
!    pointVec = (point - starPosn) * 1.e10_oc
!    r = modulus( pointVec ) 
!
!    IF (TTauriInFlow(point,grid,ignoreDisk)) THEN 
!    
!      theta = ACOS( pointVec%z  / r )
!      IF (ABS(MODULO(theta,pi)) > 1.e-10 ) THEN 
!        rM  = r / SIN(theta)**2
!      ELSE
!        rM = HUGE(rM)
!      END IF
!
!      ang = ATAN2(pointVec%y,pointVec%x)
!      IF (ang < 0.) ang = ang + twoPi
!
!      ! test if the point lies outside the accretion stream
!      IF (((rM > TTauriRinner) .AND. (rM < TTauriRouter )) .AND. &
!         (.NOT. curtains .OR. (curtains .AND.                   & 
!         (((ang > curtainsPhi1s).and.(ang < curtainsPhi1e)).or.  &
!         ((ang > curtainsPhi2s).and.(ang < curtainsPhi2e)))))) THEN 
!
!         y = SIN(theta)**2 
!
!         rho = (TTauriMdot * TTauriRstar) / (4.0 * pi * &
!         (TTauriRStar/TTauriRinner - TTauriRstar/TTauriRouter))
!         rho = rho * r**(-5.0/2.0) / SQRT( 2.0 * bigG * TTauriMstar ) 
!         rho = rho * SQRT( 4.0 - 3.0*y) / SQRT( 1.0 - y) 
!
!       ELSE
!         rho = 1.e-25
!         RETURN
!       END IF
!  
!    ELSE
!      rho = 1.e-25
!      RETURN
!    END IF
!
!  END FUNCTION TTauriDensityOld



  function testDensity(point, grid)
    use input_variables
    real :: testDensity
    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid
    real :: r
    r = modulus(point)
    testDensity = 1.e-30
    if ((r > grid%rInner).and.(r < grid%rOuter)) then
       testDensity = rho * (grid%rInner / r)**2
    endif
  end function testDensity

  function protoDensity(point, grid) result(testdensity)
    use input_variables
    real :: testDensity
    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid
    real :: r
    r = modulus(point)
    testDensity = 1.e-30
    if ((r > grid%rInner).and.(r < grid%rOuter)) then
       testDensity = grid%densityScaleFac*rho * (grid%rInner / r)**rPower
    endif
  end function protoDensity

  function wrshellDensity(point, grid) result(testdensity)
    use constants_mod
    use input_variables
    real :: testDensity
    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid
    real :: r
    r = modulus(point)
    testDensity = tiny(testDensity)
    if ((r > grid%rInner).and.(r < grid%rOuter)) then
       testDensity = mdot / (fourPi * r**2 * vTerm *1.e20)
    endif
  end function wrshellDensity

  real function benchmarkDensity(point, grid)

    use input_variables
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(octalVector), INTENT(IN) :: point
    real :: r, hr, rd
    
    rd = rOuter / 2.
    benchmarkDensity = 1.e-30
    r = sqrt(point%x**2 + point%y**2)
    if ((r > rInner).and.(r < rOuter)) then
       hr = height * (r/rd)**1.125
       benchmarkDensity = rho * ((r / rd)**(-1.))*exp(-pi/4.*(point%z/hr)**2)
    endif
!    if ((r < rInner).and.(abs(point%z) < rInner)) then
!       r = rInner
!       hr = height * (r/rd)**1.125
!       benchmarkDensity = rho * ((r / rd)**(-1.))*exp(-pi/4.*(point%z/hr)**2)*100.
!    endif
  end function benchmarkDensity

  function shakaraSunyaevDisc(point, grid) result (rhoOut)
    use input_variables
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(octalVector), INTENT(IN) :: point
    real(double) :: r, h, rhoOut, warpHeight
    real(double) :: kspiral
    real(double) :: xpoint,ypoint,rscale,r1,fac
    integer :: nspiral1
    real(double) :: phase(10)
    integer :: i
    real(double) :: phi

    nSpiral1 = 3
    do i = 1, nspiral1
       phase(i)=twoPi*real(i-1)/real(nSpiral1)
    enddo

    kspiral = grid%rOuter/Pi


    rhoOut = tiny(rhoOut)
    r = sqrt(point%x**2 + point%y**2)
    phi = atan2(point%y,point%x)
    warpHeight = 0. !cos(phi) * rInner * sin(30.*degtorad) * sqrt(rinner / r)
    if ((r > rinner).and.(r < rOuter)) then
       h = height * (r / (100.d0*autocm/1.d10))**betaDisc
       rhoOut = dble(rho0) * (dble(rInner)/r)**dble(alphaDisc) * exp(-0.5d0 * (dble(point%z-warpheight)/h)**2)

       fac =  1.d0-min(dble(r - rInner)/(0.01d0*rinner),1.d0)
       fac = exp(-fac*10.d0)
       rhoOut = rhoOut * fac
       rhoOut = max(rhoOut, tiny(rhoOut))
    endif


  end function shakaraSunyaevDisc

  function warpedDisc(point, grid) result (rhoOut)
    use input_variables
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(octalVector), INTENT(IN) :: point
    real(double) :: r, h, rhoOut, warpHeight
    real(double) :: kspiral
    real(double) :: xpoint,ypoint,rscale,r1,fac
    integer :: nspiral1
    real(double) :: phase(10)
    integer :: i
    real(double) :: phi

    rhoOut = tiny(rhoOut)
    r = sqrt(point%x**2 + point%y**2)
    phi = atan2(point%y,point%x)
    warpheight = 0.3 * rOuter * (r / rOuter)**2 * cos(phi)

    if ((r > rinner).and.(r < rOuter)) then
       h = height * (r / (100.d0*autocm/1.d10))**betaDisc
       rhoOut = dble(rho0) * (dble(rInner)/r)**dble(alphaDisc) * exp(-0.5d0 * (dble(point%z-warpheight)/h)**2)
       fac =  1.d0-min(dble(r - rInner)/(0.02d0*rinner),1.d0)
       fac = exp(-fac*10.d0)
       rhoOut = rhoOut * fac
       rhoOut = max(rhoOut, tiny(rhoOut))
    endif


  end function warpedDisc

  real(double) function melvinDensity(point, grid) 

    use input_variables
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(octalVector), INTENT(IN) :: point
    real :: r,  mStar, rCav, r_c, r_d
    real :: openingAngle, mu, mu_0, bigR
    real :: z, H_R, H_0, alpha
    real :: rhoEnv, rhoDisc
    real :: cavZ
    real :: rstar

! Envelope with cavity plus alpha disc model, as presented
! by Alvarez, Hoare and Lucas (2004).


    
    rStar = 10. * rSol
    H_0 = 10. * autocm
    mStar = 10. * mSol
    mDot = 1.11e-4 * mSol / (365.25*24.*3600.)
    rCav = 50.*auTocm
    r_c = 50. * auTocm
    r_d = 1000. * auTocm
    openingAngle = 20. * degtoRad
    alpha = 1.875
    beta  = 1.125
    rho0 = 2.e-5

    melvinDensity = 1.e-30
    r = modulus(point) * 1.e10
    mu = point%z * 1.e10 / r
    z = point%z * 1.e10

    bigR = sqrt(point%x**2 + point%y**2)*1.e10
    H_R = H_0 * (bigR / (100.*auToCm))**beta

    rhoEnv = 1.d-30
    if (r > r_c) then
       mu_0 = rtnewt(-0.2 , 1.5 , 1.e-4, r/r_c, abs(mu))
       rhoEnv = (Mdot/(8. * pi * r_c * sqrt(bigG * Mstar)))  ! Equation 1
       rhoEnv = rhoEnv * 1./sqrt(1.+mu_0) * 1./sqrt(r)
    endif
    rhoDisc = 1.d-30
    if ((bigR < r_d).and.(bigR > 10.*rStar)) then
       rhoDisc = rho0 * ((bigR/Rstar)**(-alpha))*exp(-(z**2 / (2.*H_R**2))) ! Eq 3
    endif
    rhoDisc = 1.e-30


    melvinDensity = max(rhoDisc, rhoEnv)

    if (bigR < Rcav) melvinDensity = 1.e-30

    if (bigR > Rcav) then
       cavZ = tan(pi/2.-0.5*openingAngle)*(bigR - Rcav) ! Eq 4
       if (abs(z) > cavZ) then
          melvinDensity = 1.e-30
       endif
    endif


  end function melvinDensity

  real function rtnewt(x1,x2,xacc, p1, p2) result(junk)

    real :: x1, x2, xacc, p1, p2
    integer :: jmax, j
    real ::  dx, f, df
    parameter (jmax=20)
    junk = 0.5 * (x1+x2)
    do j=1,jmax
       call equation2(junk,f,df,p1,p2)
       dx=f/df
       junk=junk-dx
       if((x1-junk)*(junk-x2).lt.0.) then
          write(*,*) 'RTNEWT: jumped out of brackets',p1,p2,junk
          stop
       endif
       if(abs(dx).lt.xacc) return
    enddo
    write(*,*) 'rtnewt exceeding maximum iterations'
  end function rtnewt
  

!  real function Equation2(mu0, eq2, deq2, r, mu)
  subroutine Equation2(mu0, eq2, deq2, r, mu)
    real :: r, mu, mu0
    real :: eq2, deq2

    eq2 = mu0**3 + (r-1.)*mu0 -r*mu
    deq2 = 3.*mu0**2 + r - 1.

  end subroutine Equation2
!  end function Equation2

  function clumpyDisc(rVec, grid) result (rhoOut)
    type(OCTALVECTOR) :: rVec
    type(GRIDTYPE) :: grid
    real :: rhoOut
    
    call findFactor(rhoOut, rVec, grid%gArray, grid%ng)
    rhoOut = rhoOut *   grid%densityScalefac
    if (modulus(rVec) < grid%rInner) rhoOut = 1.e-33

  end function clumpyDisc

  ! chris (26/05/04)
  ! This returns the *dust* density (not the gas density) and assumes a fixed
  ! gas:dust ratio.
  ! chris (19/04/05)
  ! Returns real density now (and has for some time).
  function ppdiskDensity(point, grid)
    use input_variables
    use constants_mod
    implicit none
    REAL(double) :: ppdiskDensity
    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid

    REAL(double) :: R, rInnerAU, rOuterAU
    REAL(double) :: scaleHeight, smoothScaleLength

    R = sqrt((point%x)**2 + (point%y)**2) * (1.d10/auToCm)
    rInnerAU = grid%rInner * (1.d10/auToCm)
    rOuterAU = grid%rOuter * (1.d10/auToCm)

    if ((R .lt. rInnerAU) .or. (R .gt. rOuterAU)) then
       ppdiskDensity = tiny(ppdiskDensity)
    else
!!       h = 0.15
!!       R0 = 1. (AU) 1.5e3 (10^10 cm)
!!       flaringPower = 1.5
!!       rho0 = 1.d-4 (Msol/AU**3) 5.941284951d-11 (g/cm**3)
!!       sigmaPower = 0.5
!!

!!       scaleHeight = h * (R/R0)**flaringPower
!       scaleHeight = 0.15 * (R/1.d0)**1.5
       scaleHeight = height * rHeight * (R/rHeight)**flaringPower

!!       asdf = 1.d-4 or gap
!!       ppdiskDensity = asdf * rho0/sqrt(pi)/scaleHeight * exp(-(R*z/scaleHeight)**2)/R**(sigmaPower)
! rho0 is what it's supposed to be now. May include some test reduction
! in density.
! rho0 is actually the surface density at R0 (or rHeight)
       ppdiskDensity = rho0/sqrt(2.*pi)/scaleHeight &
                       * exp(-0.5*(point%z*(1.d10/auToCm)/scaleHeight)**2)*(rHeight/R)**sigmaPower

! This test reduction is absorbed into rho0.
!       ! Reduce disc density by some amount (for testing)
!       ! The parameter 'rho' is a scaling factor here
!       ppdiskDensity = rho0 * ppdiskDensity
!!       ppdiskDensity = 1.d-2 * ppdiskDensity

       ! Add a gap into the disc
       ppdiskDensity = ppdiskDensity * fractgap(R)

       ! Smooth the inner edge of the disc (at 0.4AU)
!       ppdiskDensity = ppdiskDensity / (1 + exp(-100.*(R-0.4)))
       smoothScaleLength = height * rHeight * (rSmooth/rHeight)**flaringPower
!       ppdiskDensity = ppdiskDensity / (1. + exp(-(4.38/smoothScaleLength)*dble(R-rSmooth)))
       ppdiskDensity = ppdiskDensity / (1. + 81.d0**(dble(rSmooth-R)/smoothScaleLength))

       ! Convert the density to torus units (Msol/AU^3 --> g/cm^3)
       ppdiskDensity = ppdiskDensity * 5.941284951d-7

! We work in real mass now, rather than dust mass
!       ! We assume the gas:dust ratio is 100:1
!       ppdiskDensity = 1.d-2 * ppdiskDensity
    end if
  end function ppdiskDensity

  function fractgap(R)
      use constants_mod
      use input_variables, only : rGap, height, mPlanet, gapViscAlpha
      implicit none
      REAL(double) :: x_nu, xmu, visc, const, xx, arg, fractgap, gapfloor, x, gapalph
      REAL(double), INTENT(IN) :: R

!====================== Matt's gap ========================
!    x_nu = 2 *(MAX(3e-6,xmu)/3.0)**(1.0/3.0)
!    xmu = planetmass
!    visc = h_over_r**2 * gapalph
!    gapalph = 0.004
!
!    function  fractgap(x,x_nu,xmu,visc,gapfloor)
!      real fractgap
!      if (gapfloor .lt. 0.) then
!        fractgap=1.0
!        return
!      else
!        const=-xmu/(3d0*3.1415926535d0*visc)
!        xx=abs(x/x_nu)
!        arg=const*(1d0+xx*(3d0+4.5d0*xx+4.5d0*xx*xx))
!        fractgap=exp(arg*exp(-3d0*xx)) + gapfloor
!      endif
!      return
!    end
!==========================================================
      gapfloor = tiny(gapfloor)
!      gapfloor = 1.d-4
!      xmu = 1.d-3
      xmu = mPlanet
!      gapalph = 0.0001
      gapalph = gapViscAlpha
      x_nu = 2 *(MAX(3d-6,xmu)/3.0)**(1.0/3.0)
!      visc = 0.15**2 * gapalph
      visc = height**2 * gapalph

!      x = R - 1.d0
      x = R - rGap

      if (gapfloor .lt. 0.) then
        fractgap=1.0
      else
        const=-xmu/(3.d0*pi*visc)
        xx=abs(x/x_nu)
        arg=const*(1d0+xx*(3d0+4.5d0*xx+4.5d0*xx*xx))
        fractgap=exp(arg*exp(-3d0*xx)) + gapfloor
      endif
  end function fractgap



end module density_mod 
