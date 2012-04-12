!Another naff program storing a bunch of menial calculations
program virial_calc

implicit none

!parameters
double precision, parameter :: pi = 3.1415926535897932
double precision, parameter :: mH = 1.6733d-24
double precision, parameter :: pc = 3.0856776d18
double precision, parameter :: msol = 1.9891d33

double precision, parameter :: Tdust_low = 15.80d0
double precision, parameter :: Mcloud_low = 0.11757d2*msol
double precision, parameter :: Rcloud_low = 1.336d0*pc
double precision, parameter :: Tdust_med = 25.74d0
double precision, parameter :: Mcloud_med = 0.69212d1*msol
double precision, parameter :: Rcloud_med = 0.517d0*pc
double precision, parameter :: Tdust_hi = 20.67d0
double precision, parameter :: Mcloud_hi = 0.43506d1*msol
double precision, parameter :: Rcloud_hi = 1.222d0*pc

double precision, parameter :: T_radio = 1.d4
double precision, parameter :: T_rat_lo = 8497.d0
double precision, parameter :: T_rat_med = 8498.d0
double precision, parameter :: T_rat_hi = 8379.d0

double precision, parameter :: ne_lr = 15.94d0
double precision, parameter :: ne_lB = 16.74d0
double precision, parameter :: ne_lC = 15.48d0
double precision, parameter :: ne_lD = 14.77d0
double precision, parameter :: ne_mr = 102.76d0
double precision, parameter :: ne_mB = 102.59d0
double precision, parameter :: ne_mC = 101.39d0
double precision, parameter :: ne_mD = 92.06d0
double precision, parameter :: ne_hr = 37.88d0
double precision, parameter :: ne_hB = 40.50d0
double precision, parameter :: ne_hC = 38.19d0
double precision, parameter :: ne_hD = 39.45d0
double precision, parameter :: ne_rat_lo = 19.03d0
double precision, parameter :: ne_rat_med = 39.63d0
double precision, parameter :: ne_rat_lhi= 55.08d0

!IBL properties
double precision :: c_s_radio
double precision :: c_s_rat_lo
double precision :: c_s_rat_med
double precision :: c_s_rat_hi
double precision :: P_lr
double precision :: P_lB
double precision :: P_lC
double precision :: P_lD
double precision :: P_mr
double precision :: P_mB
double precision :: P_mC
double precision :: P_mD
double precision :: P_hr
double precision :: P_hB
double precision :: P_hC
double precision :: P_hD
double precision :: P_rat_lo
double precision :: P_rat_med
double precision :: P_rat_hi


!Neutral cloud properties

double precision :: c_cloud_low
double precision :: c_cloud_med
double precision :: c_cloud_hi
double precision :: P_cloud_low
double precision :: P_cloud_med
double precision :: P_cloud_hi


!neutral clouds

!low flux
c_cloud_low = cloudSoundSpeed(Tdust_low)
c_cloud_med = cloudSoundSpeed(Tdust_med)
c_cloud_hi = cloudSoundSpeed(Tdust_hi)
P_cloud_low = supportPressure(Mcloud_low, c_cloud_low, Rcloud_low)
P_cloud_med = supportPressure(Mcloud_med, c_cloud_med, Rcloud_med)
P_cloud_hi = supportPressure(Mcloud_hi, c_cloud_hi, Rcloud_hi)


print *, "**************************"
print *, " NEUTRAL CLOUD PROPERTIES "
print *, "**************************"
print *, "P low: ", P_cloud_low
print *, "P med: ", P_cloud_med
print *, "P hi: ", P_cloud_hi
print *, " "
print *, " "

c_s_radio = IBLSoundSpeed(T_radio)
c_s_rat_lo = IBLSoundSpeed(T_rat_lo)
c_s_rat_med = IBLSoundSpeed(T_rat_med)
c_s_rat_hi = IBLSoundSpeed(T_rat_hi)

P_lr = IBLpressure(ne_lr, c_s_radio)
P_lB = IBLpressure(ne_lB, c_s_radio)
P_lC = IBLpressure(ne_lC, c_s_radio)
P_lD = IBLpressure(ne_lD, c_s_radio)
P_mr = IBLpressure(ne_mr, c_s_radio)
P_mB = IBLpressure(ne_mB, c_s_radio)
P_mC = IBLpressure(ne_mC, c_s_radio)
P_mD = IBLpressure(ne_mD, c_s_radio)
P_hr = IBLpressure(ne_hr, c_s_radio)
P_hB = IBLpressure(ne_hB, c_s_radio)
P_hC = IBLpressure(ne_hC, c_s_radio)
P_hD = IBLpressure(ne_hD, c_s_radio)
P_rat_lo = IBLpressure(ne_rat_lo, c_s_rat_lo)
P_rat_med = IBLpressure(ne_rat_med, c_s_rat_med)
P_rat_hi = IBLpressure(ne_rat_hi, c_s_rat_hi)

print *, "**************************"
print *, "      IBL PRESSURES       "
print *, "**************************"
print *, "P_lr ", P_lr
print *, "P_lB ", P_lB
print *, "P_lC ", P_lC
print *, "P_lD ", P_lD
print *, "P_rat_lo ", P_rat_lo
print *, "P_mr ", P_mr
print *, "P_mB ", P_mB
print *, "P_mC ", P_mC
print *, "P_mD ", P_mD
print *, "P_rat_med ", P_rat_med
print *, "P_hr ", P_hr
print *, "P_hB ", P_hB
print *, "P_hC ", P_hC
print *, "P_hD ", P_hD
print *, "P_rat_hi ", P_rat_hi
print *, " "
print *, " "

contains

double precision function supportPressure(M, c, R)
  double precision, parameter :: pi = 3.1415926535897932
  double precision, parameter :: G = 6.67259d-8 
  double precision :: M
  double precision :: c
  double precision :: R

  supportPressure = ((3.d0*M*c**2)/(4.d0*pi*R**3)) - ((3.d0*G*M)/(20.d0*pi*R**4))

end function supportPressure

double precision function IBLpressure(ne, c_s)
  double precision, parameter :: mH = 1.6733d-24
  double precision :: c_s
  double precision :: ne

  IBLpressure = 2.d0*ne*mH*c_s**2
  
end function IBLpressure


double precision function cloudSoundSpeed(T)
  double precision, parameter :: mH = 1.6733d-24
  double precision, parameter :: kerg = 1.380626d-16
  double precision :: T

  isoSoundSpeed = sqtr(kerg*T/(1.36*mH))
  
end function cloudSoundSpeed

double precision function IBLSoundSpeed(T)
  double precision, parameter :: mH = 1.6733d-24
  double precision, parameter :: kerg = 1.380626d-16
  double precision :: T

  isoSoundSpeed = sqtr(kerg*T/(0.6*mH))

end function IBLSoundSpeed

end program virial_calc
