
! written by tjh


! v1.0 on 13/01/05

module gaussian_mod


  use kind_mod
  use constants_mod
  use vector_mod

  implicit none

  public

  ! The definition of the vector type

  type GAUSSIAN
     real(kind=doubleKind) :: sigma
     type(OCTALVECTOR) :: centre
     real(kind=doubleKind) :: amplitude
  end type GAUSSIAN

contains

  subroutine setGaussian(thisGaussian, centre, sigma, amplitude)
    type(GAUSSIAN) :: thisGaussian
    type(VECTOR) :: centre
    real :: sigma, amplitude, thisamplitude

!    thisamplitude = amplitude  / sqrt(twoPi * sigma**2)
    thisamplitude = amplitude  /  sigma**3
    thisGaussian%centre = centre
    thisGaussian%sigma = sigma
    thisGaussian%amplitude = thisamplitude
  end subroutine setGaussian

  subroutine findFactor(fac, position, gArray, ng)
    real :: fac
    integer :: ng, i
    type(OCTALVECTOR) :: position
    type(GAUSSIAN) :: gArray(ng)

    fac = 0.
    do i = 1, ng
       fac = fac + eval(gArray(i),position)
    end do
  end subroutine findFactor

  real function eval(thisgaussian, position)
    type(GAUSSIAN) :: thisGaussian
    type(OCTALVECTOR) :: position
    real(kind=octalKind) :: x

    x = modulus(thisGaussian%centre - position)
    eval = thisGaussian%amplitude * exp(-x**2 / (2.*thisGaussian%sigma**2))
  end function eval

  subroutine createDiscGaussians(ng, gArray)
    use input_variables
    integer :: ng, n
    real :: z, ang
    type(VECTOR) :: position
    type(GAUSSIAN) :: gArray(ng)
    real :: r1, r, h, alphaDisc, betaDisc, sigma
    real :: rAxis(1000),rProb(1000), massWeight(1000), clumpMass
    integer :: i
    alphaDisc = 2.25
    betaDisc = 1.25
    
    do i = 1, 1000
       rAxis(i) = log10(rInner) + real(i-1)/999.*(log10(rOuter)-log10(rInner))
       rAxis(i) = 10.**rAxis(i)
    enddo

    massWeight(1:1000) = rAxis(1)/rAxis(1:1000)

    rProb(1) = 0.
    do i = 2, 1000
       rProb(i) = massWeight(i) * (rAxis(i)/rAxis(1))**(betaDisc-alphaDisc)*(twoPi*rAxis(i))*(rAxis(i)-rAxis(i-1))
    enddo
    do i = 2, 1000
       rProb(i) = rProb(i) + rProb(i-1)
    enddo
    rProb(1:1000) = rProb(1:1000)/rProb(1000)
    open(99,file="prob.dat",status="unknown",form="formatted")
    do i = 1, 1000
       write(99,*) rAxis(i),rProb(i)
    enddo
    close(99)

    do i = 1, ng
       call random_number(r1)
       call locate(rProb, 1000, r1, n)
       r = rAxis(n)+(r1-rProb(n))/(rProb(n+1)-rProb(n))*(rAxis(n+1)-rAxis(n))
       clumpMass = massWeight(n) + (r1-rProb(n))/(rProb(n+1)-rProb(n))*(massWeight(n+1)-massWeight(n))
       clumpMass = 1./clumpMass
       h = height * (r / (100.*autocm/1.e10))**betaDisc
       z = gasdev() * h
       call random_number(ang)
       ang = ang * twoPi
       position%x = r * cos(ang)
       position%y = r * sin(ang)
       position%z = z * 1.e-6
       sigma = h
       call setGaussian(gArray(i), position, sigma, clumpMass)
    enddo
      

  end subroutine createDiscGaussians
    


end module gaussian_mod
