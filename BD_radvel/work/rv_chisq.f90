program rv_chisq

  use chi_square

  implicit none
  ! Given a set of RV measurements (with their uncertainties)
  ! from two epochs, this program computes the chisq prob of 
  ! the RV being constant. 

  integer, parameter :: nobj=17  ! number of objects in the list
  integer, parameter :: ntime=2  ! number of observations for each object
  real*8 :: RV(ntime, nobj)      ! radial velocities in km/s
  real*8 :: err_RV(ntime, nobj)  ! uncertainties in RV in km/s
  character(LEN=10) :: object_id(nobj)   ! name of the objects in the list.
  character(LEN=10) :: dum_a     ! 
  !
  ! File unti numbers
  integer, parameter :: LUIN =23
  integer, parameter :: LUOUT =24
  integer, parameter :: LUHIST =25
  integer, parameter :: LUEXP =26
  integer, parameter :: LUdiv =27
  !
  integer :: i, j
  ! used for chisq
  real*8  :: RV_data(ntime)  ! RV from data
  real*8  :: RV_fit(ntime)   ! RV of fitting function
  real*8  :: error(ntime)    ! uncertainties in RV 
  real*8  :: RV_ave(nobj)  ! average RV from two epochs.
  integer :: nu ! number of degree of freedom
  real*8  :: p(nobj)   ! chi_squre probability
  real*8  :: xs(nobj)  ! chi_squre
  ! for histogram
   integer, parameter :: nbin = 7
!  integer, parameter :: nbin = 14
  real*8 :: N(nbin)  ! number of object
  real*8 :: binsize 
  real*8 :: mlp, norm

  ! expected distribution chi-sq distibution
  integer,parameter :: nxs = 100
  real*8 :: p_exp(nxs)
  real*8 :: xs_exp(nxs)
  real*8 :: mlp_exp(nxs)
  real*8 :: freq(nxs)
  real*8 :: dlog10_xs, area, xs_min, xs_max, d_xs
  real*8 :: err_extra

  ! reading the data file.
  open(unit=LUIN, file="rv.dat", status="old")
  read(LUIN, "(a5)") dum_a  
  do j = 1, nobj
     do i = 1, ntime
        read(LUIN, *) object_id(j), RV(i,j), err_RV(i,j)
     end do
  end do
  close(LUIN)


  !
  ! Adding extra error values for testting
!  err_extra = 0.15d0

  write(*,*) "Enter the value additional systematic error in RV (km/s): "
  read(*,*) err_extra
  !err_extra = 0.1d0
  err_RV(:,:) = err_RV(:,:)+ err_extra

  ! Now computes the chi-squiare probability
  ! Note that the number of fitting parameter is npara = 1 (constant) here. 
  ! # of deg. of freedom (nu)  is = ntime - npara.
  ! Since ntime = 2, npara=1, we have nu=1
  nu = ntime - 1
  RV_ave(:) = 0.0d0
  do j = 1, nobj
     do i = 1, ntime
!        RV_ave(j) = RV_ave(j) + RV(i,j)  ! not weighted
        RV_ave(j) = RV_ave(j) + RV(i,j)* (1.0d0/err_RV(i,j))**2  ! weighted mean
        RV_data(i) = RV(i,j)
        error(i) = err_RV(i,j)
     end do
!     RV_ave(j) = RV_ave(j)/dble(ntime) / (1.0d0/err_RV(1,j) + 1.0d0/err_RV(2,j))
     RV_ave(j) = RV_ave(j)/ (1.0d0/err_RV(1,j)**2 + 1.0d0/err_RV(2,j)**2)
!     RV_ave(j) = RV_ave(j)/dble(ntime) 
     RV_fit(:) = RV_ave(j)  ! assuming contant at the average of two RV measuremnets

     ! now computes the chi-sq probability
     ! -- using a function in chi_square module
     xs(j) = chi_sq(RV_data, error, RV_fit, ntime)
     xs(j) = min(xs(j), 50.0d0)

     p(j) = chi_sq_prob(xs(j), nu)
     ! just for safty
     p(j) = max(p(j),1.0d-99)
  end do


  !
  ! Write the results in a file.
  open(unit=LUOUT, file="chisq_prob_const.dat", status="replace")

  write(LUOUT, "(a, a10, a14, 2x, a14, 2x, a14)") "#", "object ID", "RV mean",  "chisq (-)", "chisq prob (-)"
  do j = 1, nobj     
     write(LUOUT, "(a10, 1P, E14.3, 2x, E14.3, 2x, E14.3)") object_id(j), RV_ave(j), xs(j), p(j)
  end do

  close(LUOUT)


  !
  ! Now makes histogram for N vs -log(p) as in Maxted and Jeffries (2005)
  !
  N(:) = 0.0d0
!  binsize  = 0.5d0
  binsize  = 3.0d0/dble(nbin-1)
  do j = 1, nobj
     mlp = -1.0d0*log10(p(j))
     i = mlp/binsize
     i = i + 1
!     print *, "i = ", i 
     if (i>nbin) i = nbin
     N(i) = N(i) + 1
  end do

  ! writing result in a file
  open(unit=LUHIST, file="chisq_prob_histgram.dat", status="replace")
  write(LUHIST, "(a, a14, 2x, a14)") "#", "-log(p)", "N"
  do j = 1, nbin 
     mlp = binsize*dble(j-1)
     write(LUHIST, "(1P, E14.3, 2x, E14.3)") mlp,  N(j) 
  end do
  close(LUHIST)


  ! 
  ! For expected chi-square probabilities
  
  ! chi equally spaced between xs_min and xs_max
  xs_max = 18.0d0
  xs_min = 1.0d-5
  d_xs = (xs_max - xs_min)/dble(nxs-1)
  
  do i = 1, nxs
     xs_exp(i) = xs_min + d_xs*dble(i-1)
  end do

  ! now computes the probability
  area = 0.0d0
  nu = ntime - 0 ! parameter free (see Maxted & Jeffries 2005)
  do i = 1, nxs
     p_exp(i) = chi_sq_prob(xs_exp(i), nu)
     mlp_exp(i) = -1.0d0*log10(p_exp(i))
  end do

  ! finding the normailization
  area = 0.0d0
  do i = 1, nxs-1
     area = area + 0.5d0*(p_exp(i)+p_exp(i+1))*(mlp_exp(i+1)-mlp_exp(i))  ! not very accurate
  end do
  area = area + p_exp(nxs)*(mlp_exp(nxs) - mlp_exp(nxs-1))
  norm = dble(nobj)*binsize

!  do i = 1, nxs
!     print*, "p_exp(",i,") =", p_exp(i) 
!  end do
!  do i = 1, nxs
!     print*, "mlp_exp(",i,") =", mlp_exp(i) 
!  end do

  print *, "area = ", area
  print *, "norm = ", norm

  ! writing result in a file
  open(unit=LUEXP, file="chisq_prob_expected.dat", status="replace")
  write(LUEXP, "(a, a14, 2x, a14)") "#", "-log(p)", "N_expected"
  do i = 1, nxs
     write(LUEXP, "(1P, E14.3, 2x, E14.3)") mlp_exp(i),  p_exp(i)*norm/area
!     write(LUEXP, "(1P, E14.3, 2x, E14.3)") mlp_exp(i),  p_exp(i)
     ! The last term above is scaleing the ditribution to the total area of the histrogram (data) 
     ! just depends on the total number of object and the bin size.
     
  end do
  close(LUEXP)



  !
  ! writing the deviation from the mean of each data points
  !
  open(unit=LUDIV, file="rv_wrt_mean.dat", status="replace")
  write(LUDIV, "(a, a14, 2x, a14)") "obj id", "RV-<RV>"
  do j = 1, nobj
     do i = 1, ntime
        write(LUDIV,"(a14, 2x, 1P, E14.3)") object_id(j), RV(i,j)-RV_ave(j)
     end do
     write(LUDIV,"(a)") " "
  end do
  close(LUDIV)
  

end program rv_chisq
