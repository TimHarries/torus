!A brutally ugly program to do the radio slog calculations for me
program radio_analyzer

implicit none

!constants
double precision, parameter :: pi = 3.1415926535897932
double precision, parameter :: pcToCm = 3.0856776d18
double precision, parameter :: auToCm = 1.495979d13
double precision, parameter :: degToRad = pi/180.d0
double precision, parameter :: radToDeg = 180.d0/pi
double precision, parameter :: radiansToArcSec = radToDeg * 60.d0 * 60.d0
double precision, parameter :: ArcSecsToRadians = pi/6.48d5 ! 6.48d5 = 180*60*60
double precision, parameter :: arcsec = 1.d0/3600.d0 * degtorad


!measured average fluxes (basically inputs)
double precision, parameter :: lr = 0.009384
double precision, parameter :: lB = 0.01035
double precision, parameter :: lC = 0.008854
double precision, parameter :: lD = 0.008057
double precision, parameter :: mr = 0.1509
double precision, parameter :: mB = 0.1504
double precision, parameter :: mC = 0.1469
double precision, parameter :: mD = 0.1211
double precision, parameter :: hr = 0.04847
double precision, parameter :: hB = 0.0554
double precision, parameter :: hC = 0.04926
double precision, parameter :: hD = 0.05257

!grid properties
double precision, parameter :: theta_pix_D = 2.3067336 !pixel angular diameter
double precision, parameter :: R_lo_pix = 119.45       !size of masked areas
double precision, parameter :: R_med_pix = 46.21
double precision, parameter :: R_hi_pix = 109.24
double precision, parameter :: gridL = 9.25e5          !grid size (in au)
integer, parameter :: npix=401                         !num pixels on image

!calculation parameters
double precision, parameter :: A = 1.24d10
double precision, parameter :: B = 122.41d0
double precision, parameter :: C= 4.4d-3
double precision, parameter :: nu_GHz = 1.5 !GHz
double precision, parameter :: Te = 1.d4
double precision, parameter :: eta = 0.2d0
double precision, parameter :: R_l = 1.336 !cloud radii in pc
double precision, parameter :: R_m = 0.517
double precision, parameter :: R_h = 1.222

!other
double precision :: gridSize, dx, distance
double precision :: theta_l_R_rad, theta_m_R_rad, theta_h_R_Rad    !angular radii of masked areas radians
double precision :: theta_l_R, theta_m_R, theta_h_R    !angular radii of masked areas as
double precision :: omega_l, omega_m, omega_h          !solid angle subtended by masked areas

!Integrated fluxes
double precision :: F_lr, F_lB, F_lC, F_lD
double precision :: F_mr, F_mB, F_mC, F_mD
double precision :: F_hr, F_hB, F_hC, F_hD

!Incident ionizing fluxes
double precision :: Phi_lr, Phi_lB, Phi_lC, Phi_lD
double precision :: Phi_mr, Phi_mB, Phi_mC, Phi_mD
double precision :: Phi_hr, Phi_hB, Phi_hC, Phi_hD

!n
double precision :: ne_lr, ne_lB, ne_lC, ne_lD
double precision :: ne_mr, ne_mB, ne_mC, ne_mD
double precision :: ne_hr, ne_hB, ne_hC, ne_hD

!Mass loss rates
double precision :: M_lr, M_lB, M_lC, M_lD
double precision :: M_mr, M_mB, M_mC, M_mD
double precision :: M_hr, M_hB, M_hC, M_hD

!photo-evaporative strength
double precision :: S_lr, S_lB, S_lC, S_lD
double precision :: S_mr, S_mB, S_mC, S_mD
double precision :: S_hr, S_hB, S_hC, S_hD

gridSize = gridL*auToCm
dx = gridSize/dble(npix)
distance = 1000.*pctocm

theta_l_R = R_lo_pix*theta_pix_D
theta_m_R = R_med_pix*theta_pix_D
theta_h_R = R_hi_pix*theta_pix_D
theta_l_R_rad = R_lo_pix*theta_pix_D*ArcSecsToRadians
theta_m_R_rad = R_med_pix*theta_pix_D*ArcSecsToRadians
theta_h_R_rad = R_hi_pix*theta_pix_D*ArcSecsToRadians

omega_l = 2.d0*pi*(1.d0-cos(theta_l_R_rad))
omega_m = 2.d0*pi*(1.d0-cos(theta_m_R_rad))
omega_h = 2.d0*pi*(1.d0-cos(theta_h_R_rad))

print *, "Unmasked solid angles"
print *, "Low flux: ", omega_l, "strad"
print *, "Medium flux: ", omega_m, "strad"
print *, "High flux: ", omega_h, "strad"

print *, "theta_l_R ", theta_l_R
print *, "theta_m_R ", theta_m_R
print *, "theta_h_R ", theta_h_R

!Calculate the integrated fluxes
F_lr = lr*omega_l*1.d9
F_lB = lB*omega_l*1.d9
F_lC = lC*omega_l*1.d9
F_lD = lD*omega_l*1.d9
F_mr = mr*omega_m*1.d9
F_mB = mB*omega_m*1.d9
F_mC = mC*omega_m*1.d9
F_mD = mD*omega_m*1.d9
F_hr = hr*omega_h*1.d9
F_hB = hB*omega_h*1.d9
F_hC = hC*omega_h*1.d9
F_hD = hD*omega_h*1.d9

print *, " "
print *, "Integrated fluxes (mJy)"
print *, "-----------------------"
print *, "Low flux: "
print *, "r ", F_lr
print *, "B ", F_lB
print *, "C ", F_lC
print *, "D ", F_lD
print *, "Medium flux: "
print *, "r ", F_mr
print *, "B ", F_mB
print *, "C ", F_mC
print *, "D ", F_mD
print *, "High flux: "
print *, "r ", F_hr
print *, "B ", F_hB
print *, "C ", F_hC
print *, "D ", F_hD

!calculate the ionizing fluxes

Phi_lr = calcPhi(F_lr, 2.d0*theta_l_R)
Phi_lB = calcPhi(F_lB, 2.d0*theta_l_R)
Phi_lC = calcPhi(F_lC, 2.d0*theta_l_R)
Phi_lD = calcPhi(F_lB, 2.d0*theta_l_R)
Phi_mr = calcPhi(F_mr, 2.d0*theta_m_R)
Phi_mB = calcPhi(F_mB, 2.d0*theta_m_R)
Phi_mC = calcPhi(F_mC, 2.d0*theta_m_R)
Phi_mD = calcPhi(F_mB, 2.d0*theta_m_R)
Phi_hr = calcPhi(F_hr, 2.d0*theta_h_R)
Phi_hB = calcPhi(F_hB, 2.d0*theta_h_R)
Phi_hC = calcPhi(F_hC, 2.d0*theta_h_R)
Phi_hD = calcPhi(F_hB, 2.d0*theta_h_R)

print *, " "
print *, "Incident ionizing fluxes (x10^9)"
print *, "-----------------------"
print *, "Low flux: "
print *, "r ", Phi_lr/1.d9
print *, "B ", Phi_lB/1.d9
print *, "C ", Phi_lC/1.d9
print *, "D ", Phi_lD/1.d9
print *, "Medium flux: "
print *, "r ", Phi_mr/1.d9
print *, "B ", Phi_mB/1.d9
print *, "C ", Phi_mC/1.d9
print *, "D ", Phi_mD/1.d9
print *, "High flux: "
print *, "r ", Phi_hr/1.d9
print *, "B ", Phi_hB/1.d9
print *, "C ", Phi_hC/1.d9
print *, "D ", Phi_hD/1.d9



!calculate the electron densitites

ne_lr = calcne(F_lr, 2.d0*theta_l_R,R_l)
ne_lB = calcne(F_lB, 2.d0*theta_l_R,R_l)
ne_lC = calcne(F_lC, 2.d0*theta_l_R,R_l)
ne_lD = calcne(F_lB, 2.d0*theta_l_R,R_l)
ne_mr = calcne(F_mr, 2.d0*theta_m_R,R_m)
ne_mB = calcne(F_mB, 2.d0*theta_m_R,R_m)
ne_mC = calcne(F_mC, 2.d0*theta_m_R,R_m)
ne_mD = calcne(F_mB, 2.d0*theta_m_R,R_m)
ne_hr = calcne(F_hr, 2.d0*theta_h_R,R_h)
ne_hB = calcne(F_hB, 2.d0*theta_h_R,R_h)
ne_hC = calcne(F_hC, 2.d0*theta_h_R,R_h)
ne_hD = calcne(F_hB, 2.d0*theta_h_R,R_h)

print *, " "
print *, "Electron densities"
print *, "-----------------------"
print *, "Low flux: "
print *, "r ", ne_lr
print *, "B ", ne_lB
print *, "C ", ne_lC
print *, "D ", ne_lD
print *, "Medium flux: "
print *, "r ", ne_mr
print *, "B ", ne_mB
print *, "C ", ne_mC
print *, "D ", ne_mD
print *, "High flux: "
print *, "r ", ne_hr
print *, "B ", ne_hB
print *, "C ", ne_hC
print *, "D ", ne_hD



!calculate the mass loss rates

M_lr = C*(Phi_lr**0.5)*(R_l**1.5)
M_lB = C*(Phi_lB**0.5)*(R_l**1.5)
M_lC = C*(Phi_lC**0.5)*(R_l**1.5)
M_lD = C*(Phi_lD**0.5)*(R_l**1.5)
M_mr = C*(Phi_mr**0.5)*(R_m**1.5)
M_mB = C*(Phi_mB**0.5)*(R_m**1.5)
M_mC = C*(Phi_mC**0.5)*(R_m**1.5)
M_mD = C*(Phi_mD**0.5)*(R_m**1.5)
M_hr = C*(Phi_hr**0.5)*(R_h**1.5)
M_hB = C*(Phi_hB**0.5)*(R_h**1.5)
M_hC = C*(Phi_hC**0.5)*(R_h**1.5)
M_hD = C*(Phi_hD**0.5)*(R_h**1.5)

print *, " "
print *, "Mass loss rates (Mokyr^-1)"
print *, "-----------------------"
print *, "Low flux: "
print *, "r ", M_lr/1.d3
print *, "B ", M_lB/1.d3
print *, "C ", M_lC/1.d3
print *, "D ", M_lD/1.d3
print *, "Medium flux: "
print *, "r ", M_mr/1.d3
print *, "B ", M_mB/1.d3
print *, "C ", M_mC/1.d3
print *, "D ", M_mD/1.d3
print *, "High flux: "
print *, "r ", M_hr/1.d3
print *, "B ", M_hB/1.d3
print *, "C ", M_hC/1.d3
print *, "D ", M_hD/1.d3


!calculate the photo-evaporative strenghts

S_lr = M_lr/(4.d0*pi*R_l**2)
S_lB = M_lB/(4.d0*pi*R_l**2)
S_lC = M_lC/(4.d0*pi*R_l**2)
S_lD = M_lD/(4.d0*pi*R_l**2)
S_mr = M_mr/(4.d0*pi*R_m**2)
S_mB = M_mB/(4.d0*pi*R_m**2)
S_mC = M_mC/(4.d0*pi*R_m**2)
S_mD = M_mD/(4.d0*pi*R_m**2)
S_hr = M_hr/(4.d0*pi*R_h**2)
S_hB = M_hB/(4.d0*pi*R_h**2)
S_hC = M_hC/(4.d0*pi*R_h**2)
S_hD = M_hD/(4.d0*pi*R_h**2)


print *, " "
print *, "Photo-evaporative strengths (Mokyr^-1pc^-2)"
print *, "-----------------------"
print *, "Low flux: "
print *, "r ", S_lr/1.d3
print *, "B ", S_lB/1.d3
print *, "C ", S_lC/1.d3
print *, "D ", S_lD/1.d3
print *, "Medium flux: "
print *, "r ", S_mr/1.d3
print *, "B ", S_mB/1.d3
print *, "C ", S_mC/1.d3
print *, "D ", S_mD/1.d3
print *, "High flux: "
print *, "r ", S_hr/1.d3
print *, "B ", S_hB/1.d3
print *, "C ", S_hC/1.d3
print *, "D ", S_hD/1.d3

contains

double precision function calcPhi(flux, theta)

  double precision, parameter :: A = 1.24d10
  double precision, parameter :: nu_GHz = 1.5 !GHz
  double precision, parameter :: Te = 1.d4
  double precision :: flux
  double precision :: theta

  calcPhi = A*flux*(Te**0.35)*(nu_GHz**0.1)/(theta**2)

end function calcPhi


double precision function calcNe(flux, theta, R)
  double precision, parameter :: B = 122.41
  double precision, parameter :: nu_GHz = 1.5 !GHz
  double precision, parameter :: Te = 1.d4
  double precision, parameter :: eta=0.2d0
  double precision :: flux
  double precision :: theta
  double precision :: R

  
  calcNe = B*((flux*(Te**0.35)*(nu_GHz**0.1)*(theta**(-2)))/(eta*R))**0.5

end function calcNe
end program radio_analyzer
