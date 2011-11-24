  subroutine createDiscGaussians(ng, gArray)
    use inputs_mod
    use utils_mod, only: locate, gasdev
    integer :: ng, n
    real :: z, ang
    type(VECTOR) :: position
    type(GAUSSIAN) :: gArray(ng)
    real :: r1, r, h, sigma
    real :: rAxis(1000),rProb(1000), massWeight(1000), clumpMass
    integer :: i
    
    do i = 1, 1000
       rAxis(i) = log10(rInner) + real(i-1)/999.*(log10(rOuter)-log10(rInner))
       rAxis(i) = 10.**rAxis(i)
    enddo

    massWeight(1:1000) = rAxis(1)/rAxis(1:1000)

    rProb(1) = 0.
    do i = 2, 1000
       rProb(i) = real(massWeight(i) * (rAxis(i)/rAxis(1))**(betaDisc-alphaDisc)*(twoPi*rAxis(i))*(rAxis(i)-rAxis(i-1)))
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
       call randomNumberGenerator(getReal=r1)
       call locate(rProb, 1000, r1, n)
       r = rAxis(n)+(r1-rProb(n))/(rProb(n+1)-rProb(n))*(rAxis(n+1)-rAxis(n))
       clumpMass = massWeight(n) + (r1-rProb(n))/(rProb(n+1)-rProb(n))*(massWeight(n+1)-massWeight(n))
       clumpMass = 1./clumpMass
       h = real(height * (r / (100.*autocm/1.e10))**betaDisc)
       z = gasdev() * h
       call randomNumberGenerator(getReal=ang)
       ang = ang * real(twoPi)
       position%x = r * cos(ang)
       position%y = r * sin(ang)
       position%z = z * 1.e-6
       sigma = h
       call setGaussian(gArray(i), position, sigma, clumpMass)
    enddo
      

  end subroutine createDiscGaussians
    

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

