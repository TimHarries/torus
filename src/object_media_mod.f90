module object_media_mod
    use kind_mod
    use unix_mod
    implicit none

    type mediumProps
        real(double) :: muAbsorb
        real(double) :: muScatter
        real(double) :: gFactor
        real(double) :: nRefract
    end type mediumProps

    type medium
        character(len=80) :: label
        !integer :: nlambda
        real(double), pointer :: lambdaList(:) => null()
        type(mediumProps), pointer :: propList(:) => null()
    end type medium

!    integer :: nTissue
!    type(TISSUETYPE), pointer :: tissueArray(:)
!    integer :: nh2o, nHb
!    real(double), pointer :: lamH2O(:), muH2O(:)
!    real(double), pointer :: lamHb(:), muHb(:), muHbo2(:)


contains


!  real(double)  function interpArray(x, nx, xArray, yArray) result (out)
!    real(double) :: x, xArray(:), yArray(:)
!    integer :: nx, i
!
!    call locate(xArray, nx, x, i)
!    out = yArray(i) + (yArray(i+1)-yArray(i)) * (x - xArray(i))/(xArray(i+1)-xArray(i))
!  end function interpArray

    function createNewMedium(nameString, lambdaList) result(med)
        character(len=*), intent(in) :: nameString
        !real(double), allocatable, intent(in) :: lambdaList(:)
        real(double), intent(in) :: lambdaList(:)
        type(medium) :: med

        integer :: nLambda
        nLambda = size(lambdaList)

        med%label = trim(nameString)
        allocate(med%lambdaList(nLambda), &
                 med%propList(nLambda))

        med%lambdaList = lambdaList
        med%propList = mediumProps(0,0,0,1)

        return
    end function createNewMedium


    subroutine writeProps(filename, thismedium)
      character(len=*) :: filename
      type(medium) :: thisMedium
      type(mediumProps) :: props
      real(double) :: lam
      integer :: i

      open(33,file=filename, form="formatted", status="unknown")
      do i = 1, 100
         lam = 400.d0 + 500.d0 * dble(i-1)/99.d0
         props = getMediumProps(thisMedium, lam)
         write(33, '(1p,4e12.3)') lam, props%muAbsorb, props%muScatter, props%gfactor
      enddo
      close(33)
    end subroutine writeProps



    function getMediumProps(med, lambda) result(props)
        ! For given medium, uses linear interpolation to find to find the material properties
        ! at the given wavelength lambda.
        ! Returns a type(mediumProps)
        use messages_mod, only: writeWarning
        use utils_mod, only: locate

        implicit none
        type(medium), intent(in) :: med
        real(double), intent(in) :: lambda
        type(mediumProps) :: props

        integer :: i
        real(double) :: w
        type(mediumProps) :: p0, p1

        if(.not. associated(med%propList)) then
            call writeWarning('Attempted access an undefined medium')
            stop
        end if

        props = med%propList(1)

        if(.not. associated(med%lambdaList)) return

        ! Linear interpolation
        call locate(med%lambdaList, size(med%lambdaList), lambda, i)
        w = (lambda - med%lambdaList(i))/(med%lambdaList(i+1) - med%lambdaList(i))
        p0 = med%propList(i)
        p1 = med%propList(i+1)

        props%muAbsorb = (1.d0-w)*p0%muAbsorb + w*p1%muAbsorb
        props%muScatter = (1.d0-w)*p0%muScatter + w*p1%muScatter
        props%gFactor = (1.d0-w)*p0%gFactor + w*p1%gFactor
        props%nRefract = (1.d0-w)*p0%nRefract + w*p1%nRefract

    end function getMediumProps

    function defineAir() result(air)

      type(medium) :: air

      allocate(air%propList(1:1))
      air%propList(1)%gfactor = 0.d0
      air%propList(1)%muScatter = 0.d0
      air%propList(1)%muAbsorb = 0.d0
      air%propList(1)%nRefract = 1.000277
    endfunction defineAir

    function defineRubber() result(rubber)

      type(medium) :: rubber

      allocate(rubber%propList(1:1))
      rubber%propList(1)%gfactor = 0.
      rubber%propList(1)%muScatter = 0.
      rubber%propList(1)%muAbsorb = 1e10
      rubber%propList(1)%nRefract = 1.000277
      rubber%label = "rubber"
    endfunction defineRubber

    function defineFatSolution(solutionType) result(fatSolution)

! fits from Michels, Foschm, Kienle, 2008, Optics Express, 16, 5907

      type(medium) :: fatSolution
      type(medium) :: water
      type(mediumprops) :: waterProps
      character(len=*) :: solutionType
      integer, parameter :: nLambda = 100
      real(double), parameter :: lamStart = 400.d0, lamEnd = 900.d0 ! wavelenghs in nm
      real(double) :: muAbsWater, muAbsSoy
      real(double) :: a_g, a_s, a_soy, b_s, b_soy, sigma_water, sigma_soy, x0_soy, y0_g, sfac
      integer :: i

      fatSolution%label = trim(solutionType)
      water  = defineWater(20.d0)

      allocate(fatSolution%propList(1:nLambda))
      allocate(fatSolution%lambdaList(1:nLambda))
      do i = 1, nLambda
         fatSolution%lambdaList(i) = lamStart + (lamEnd-lamStart)*dble(i-1)/dble(nLambda-1)
      enddo

      a_soy = 1.171d5
      b_soy = -3.659d1
      x0_soy = -3.210d2
      

      select case(solutionType)
      case("intralipid10")
         y0_g = 1.018d0
         a_g = -8.82d-4
         a_s = 4.857d8
         b_s = -2.644d0
         sigma_soy = 0.12
         sigma_water = 0.88
         sfac = 1.d0
      case("intralipid01")
         y0_g = 1.018d0
         a_g = -8.82d-4
         a_s = 4.857d8
         b_s = -2.644d0
         sigma_soy = 0.012
         sigma_water = 0.988
         sfac = 0.1d0
      case("intralipid001")
         y0_g = 1.018d0
         a_g = -8.82d-4
         a_s = 4.857d8
         b_s = -2.644d0
         sigma_soy = 0.0012
         sigma_water = 0.9988
         sfac = 0.01d0
      case DEFAULT
         write(*,*) "Optical properties not defined for this concentration"
         stop
      end select
      do i = 1, nLambda
         fatSolution%propList(i)%gfactor = y0_g + a_g * fatSolution%lambdaList(i)
         fatSolution%propList(i)%muScatter = sfac * a_s *  fatSolution%lambdaList(i)**b_s * 10.d0 ! per mm to per cm

         muAbsSoy = (a_soy / (1.d0+exp(-(fatSolution%lambdaList(i)-x0_soy)/b_soy))) * 10.d0

         waterProps = getMediumProps(water, fatSolution%lambdalist(i))
         muAbsWater = waterProps%muAbsorb 
         fatSolution%propList(i)%muAbsorb = sigma_soy * muAbsSoy + sigma_water * muAbsWater
         fatSolution%propList(i)%nRefract = waterProps%nRefract * &
              0.9d0 + 0.1d0*(waterProps%nRefract+0.14d0) 
      enddo
    end function defineFatSolution

          
    function defineWater(tempIn) result(water)
        ! Defines the medium properties of water at 25 degrees C, or given (optional) temperature
        ! Absorbtion rate is obtained by reading from file
        ! Refractive index is computed using an approximation
        implicit none
        real(double), intent(in), optional :: tempIn
        type(medium) :: water

        real(double) :: temp
        integer :: i, nLambda
        character(len=120) :: dataDirectory, charLine, filename

        ! constants taken from
        ! The International Association for the Properties of Water and Steam
        ! Release on the Refractive Index of Ordinary Water Substance
        !     as a Function of Wavelength, Temperature and Pressure
        ! http://www.iapws.org/relguide/rindex.pdf
        real(double) :: TBar, rBar, lBar, rhs
        real(double) :: a0, a1, a2, a3, a4, a5, a6, a7, lUV, lIR
        parameter( a0 = 0.244258d0, &
                   a1 = 9.746345d-3, &
                   a2 = -3.732350d-3, &
                   a3 = 2.686785d-4, &
                   a4 = 1.589206d-3, &
                   a5 = 2.459343d-3, &
                   a6 = 0.900705d0, &
                   a7 = -1.666262d-2, &
                   lUV = 0.229202d0, &
                   lIR = 5.432937d0)


        temp = 25 ! degrees Celsius
        if(present(tempIn)) temp = tempIn

        write (water%label, '(a,f5.1,a)') 'water (T=',temp,' C)'

        ! read absorbtion from file
        call unixGetenv("TORUS_DATA", dataDirectory, i)
        filename = trim(dataDirectory)//"/"//"buiteveld1994.dat"

        ! find size first
        open(21, file=filename,status="old", form="formatted")
        ! FIXME: issues error/warning if file is not found
        nLambda = 0
10      continue
        read(21,*,end=20) charLine
        if (charLine(1:1)=="#") goto 10
        nLambda = nLambda + 1
        goto 10
20      continue
        close(21)

        allocate(water%lambdaList(nLambda), &
                 water%propList(nLambda))

        water%propList = mediumProps(0,0,0,1)

        ! now read data
        open(21, file=filename,status="old", form="formatted")
        nLambda = 0
30      continue
        read(21,'(a)',end=40) charLine
        if (charLine(1:1)=="#") goto 30
        nLambda = nLambda + 1
        read(charLine,*) water%lambdaList(nLambda), water%propList(nLambda)%muAbsorb
        goto 30
40      continue
        close(21)

        ! compute refractive index
        TBar = (temp+273.15d0)/273.15d0
        rBar = 1.d0

        do i=1,nLambda
            lBar = water%LambdaList(i)/589.d0

            rhs = a0 + a1*rBar + a2*TBar + a3*lBar**2*TBar + a4/lBar + &
                  a5/(lBar**2-lUV**2) + a6/(lBar**2-lIR**2) + a7*rBar**2

            water%propList(i)%nRefract = sqrt( (1.d0 + 2*rhs)/(1.d0 - rhs) )
        end do

        return

    end function defineWater


    function defineGlass(refIndexfilename) result(glass)
      use constants_mod
      use utils_mod
      character(len=*) :: refIndexfilename
      character(len=80) :: filename, junk, dataDirectory
      type(medium) :: glass
      real(double) :: lambda_n(1000), lambda_k(1000), n(1000), k(1000), alpha(1000)
      integer :: i, j, in, ik
      call unixGetenv("TORUS_DATA", dataDirectory, i)
      filename = trim(dataDirectory)//"/"//trim(refIndexFilename)

      open(20, file=filename, status="old", form="formatted")

      read(20,*) junk
      read(20,*) junk
      in = 1
    10 read(20,*, err=20) lambda_n(in), n(in)
      in = in+1
      goto 10
 20  continue
      in = in - 1
      ik = 1
      read(20,*) junk
 30    read(20,*, end=40) lambda_k(ik), k(ik)
      ik = ik + 1
      goto 30
 40   continue
      ik = ik - 1
      close(20)

      alpha(1:ik) = 0.!fourPi*k(1:ik)/lambda_k(1:ik)
      lambda_n(1:in) = lambda_n(1:in) * 1000.d0 ! microns to nm
      lambda_k(1:ik) = lambda_k(1:ik) * 1000.d0 ! microns to nm

      glass%label = refIndexfilename
      allocate(glass%lambdaList(1:in))
      allocate(glass%propList(1:in))
      glass%lambdaList(1:in) = lambda_n(1:in)

      do i = 1, in
         glass%propList(i)%nRefract = n(i)
         glass%propList(i)%muScatter = 0.d0
         glass%propList(i)%gFactor = 0.d0
         call locate(lambda_k, ik, lambda_n(i), j)
         glass%propList(i)%muAbsorb = alpha(j) + (alpha(j+1)-alpha(j)) * &
              (lambda_n(i) - lambda_k(j))/(lambda_k(j+1)-lambda_k(j))
      enddo

    end function defineGlass

    function defineHemoglobin(bound) result(Hb)        ! Defines the medium properties of hemoglobin, either bound or unbound
        ! Absorbtion rates are obtained by reading from file
        implicit none
        logical, intent(in) :: bound    ! Flag to read either oxygen (un)bound hemoglobin
        type(medium) :: Hb

        real(double) :: mua_hbo2, mua_hb
        integer :: i, nLambda
        character(len=120) :: dataDirectory, charLine, filename

        hb%label = "blood"

        call unixGetenv("TORUS_DATA", dataDirectory, i)
        filename = trim(dataDirectory)//"/"//"prahl1999.dat"


        ! find size
        open(21, file=filename,status="old", form="formatted")
        ! FIXME: issues error/warning if file is not found
        nLambda = 0
10      continue
        read(21,*,end=20) charLine
        if (charLine(1:1)=="#") goto 10
        nLambda = nLambda + 1
        goto 10
20      continue
        close(21)

        ! allocate and initialize
        allocate(Hb%lambdaList(nLambda), &
                 Hb%propList(nLambda))

        Hb%propList = mediumProps(0,0,0,1)

        open(21, file=filename,status="old", form="formatted")
        nLambda = 0
30      continue
        read(21,'(a)',end=40) charLine
        if (charLine(1:1)=="#") goto 30
        nLambda = nLambda + 1
        read(charLine,*) Hb%lambdaList(nLambda), mua_hbo2, mua_hb
        if (bound) then
            Hb%propList(nLambda)%muAbsorb = mua_hbo2
        else
            Hb%propList(nLambda)%muAbsorb = mua_hb
        end if
        goto 30
40      continue
        close(21)

        !mua_hbo2 = mua_hbo2 * 2.303d0 * 150.d0 / 64500.d0
        !mua_hb   = mua_hbo2 * 2.303d0 * 150.d0 / 64500.d0 ! FIXME: typo?!

        Hb%propList%muAbsorb = Hb%propList%muAbsorb * 2.303d0 * 150.d0 / 64500.d0
        Hb%propList%nRefract = 1.35 ! Zhernovaya et al. 2011 Phys Med Biol 56 4013

    end function defineHemoglobin


end module object_media_mod



