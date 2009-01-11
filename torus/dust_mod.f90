module dust_mod

  use gridtype_mod
  use phasematrix_mod
  use grid_mod
  use constants_mod
  use octal_mod
  use amr_mod
  use messages_mod

  implicit none
  public

contains

  subroutine getRefractiveIndex(lambda, nLambda, graintype, mReal, mImg)
    use unix_mod, only: unixGetenv

    ! gr1_dl - graphite orth  .25
    ! gr2_dl - graphite para  .125

    integer, parameter :: npLnk = 98
    real(single) :: n_sil_ow(npLnk),k_sil_ow(npLnk),n_sil_oc(npLnk), &
         k_sil_oc(npLnk), n_sil_dl(npLnk), k_sil_dl(npLnk),             &
         n_amc_hn(npLnk), k_amc_hn(npLnk), n_sic_pg(npLnk),             &
         k_sic_pg(npLnk), n_gr1_dl(npLnk), k_gr1_dl(npLnk),             &
         n_gr2_dl(npLnk), k_gr2_dl(npLnk), lam_nk(npLnk)                
    DATA lam_nk/1.00E-02,2.00E-02,3.00E-02,4.00E-02,5.00E-02,             &
         6.00E-02,8.00E-02,1.00E-01,1.20E-01,1.50E-01,2.00E-01,2.50E-01,    &
         3.00E-01,3.60E-01,4.40E-01,5.50E-01,7.00E-01,8.50E-01,1.00E+00,    &
         1.15E+00,1.30E+00,1.70E+00,2.00E+00,2.20E+00,2.70E+00,3.00E+00,    &
         3.50E+00,4.00E+00,4.50E+00,5.00E+00,5.50E+00,6.00E+00,6.50E+00,    &
         7.00E+00,7.50E+00,8.00E+00,8.50E+00,9.00E+00,9.40E+00,9.55E+00,    &
         9.70E+00,9.85E+00,1.00E+01,1.05E+01,1.10E+01,1.13E+01,1.16E+01,    &
         1.20E+01,1.25E+01,1.30E+01,1.35E+01,1.40E+01,1.45E+01,1.50E+01,    &
         1.60E+01,1.70E+01,1.80E+01,1.90E+01,2.00E+01,2.20E+01,2.30E+01,    &
         2.40E+01,2.50E+01,2.60E+01,2.70E+01,2.80E+01,3.00E+01,3.50E+01,    &
         4.00E+01,4.50E+01,5.00E+01,5.50E+01,6.00E+01,6.50E+01,7.00E+01,    &
         7.50E+01,8.00E+01,8.50E+01,9.00E+01,9.50E+01,1.00E+02,1.05E+02,    &
         1.10E+02,1.20E+02,1.30E+02,1.40E+02,1.50E+02,2.00E+02,2.50E+02,    &
         3.00E+02,4.00E+02,5.00E+02,7.00E+02,1.30E+03,4.00E+03,1.30E+04,    &
         2.00E+04,3.60E+04/                                                 
    DATA n_sil_ow/1.81E+00,1.81E+00,1.81E+00,1.81E+00,1.81E+00,           &
         1.81E+00,1.81E+00,1.81E+00,1.81E+00,1.81E+00,1.81E+00,1.81E+00,    &
         1.81E+00,1.81E+00,1.82E+00,1.84E+00,1.85E+00,1.85E+00,1.85E+00,    &
         1.85E+00,1.85E+00,1.87E+00,1.88E+00,1.88E+00,1.89E+00,1.89E+00,    &
         1.88E+00,1.87E+00,1.87E+00,1.85E+00,1.83E+00,1.81E+00,1.79E+00,    &
         1.75E+00,1.71E+00,1.64E+00,1.55E+00,1.45E+00,1.40E+00,1.43E+00,    &
         1.44E+00,1.46E+00,1.51E+00,1.73E+00,1.91E+00,2.00E+00,2.11E+00,    &
         2.21E+00,2.21E+00,2.16E+00,2.09E+00,2.06E+00,2.03E+00,1.98E+00,    &
         1.91E+00,1.90E+00,1.96E+00,2.04E+00,2.12E+00,2.27E+00,2.34E+00,    &
         2.36E+00,2.38E+00,2.40E+00,2.42E+00,2.43E+00,2.45E+00,2.60E+00,    &
         2.66E+00,2.70E+00,2.74E+00,2.76E+00,2.77E+00,2.79E+00,2.80E+00,    &
         2.81E+00,2.82E+00,2.83E+00,2.85E+00,2.86E+00,2.87E+00,2.87E+00,    &
         2.87E+00,2.87E+00,2.88E+00,2.88E+00,2.88E+00,2.89E+00,2.90E+00,    &
         2.90E+00,2.90E+00,2.90E+00,2.91E+00,2.91E+00,2.91E+00,2.91E+00,    &
         2.91E+00,2.91E+00/                                                 
    DATA k_sil_ow/9.99E-02,9.99E-02,9.99E-02,9.99E-02,9.99E-02,           &
         9.99E-02,9.99E-02,9.99E-02,9.99E-02,9.99E-02,9.99E-02,9.99E-02,    &
         9.99E-02,9.99E-02,9.08E-02,7.17E-02,6.03E-02,5.64E-02,5.46E-02,    &
         6.12E-02,6.75E-02,7.23E-02,7.33E-02,6.96E-02,6.19E-02,5.88E-02,    &
         5.48E-02,5.24E-02,5.07E-02,5.14E-02,5.24E-02,5.50E-02,5.79E-02,    &
         6.36E-02,7.07E-02,9.03E-02,1.31E-01,2.63E-01,4.31E-01,4.96E-01,    &
         5.60E-01,6.30E-01,7.12E-01,8.10E-01,7.94E-01,7.86E-01,7.68E-01,    &
         6.13E-01,4.93E-01,3.98E-01,3.74E-01,3.69E-01,3.66E-01,3.72E-01,    &
         4.39E-01,5.63E-01,6.68E-01,7.23E-01,7.51E-01,7.55E-01,7.34E-01,    &
         7.10E-01,6.85E-01,6.61E-01,6.37E-01,6.44E-01,6.58E-01,6.28E-01,    &
         5.78E-01,5.23E-01,4.69E-01,4.34E-01,4.13E-01,3.92E-01,3.70E-01,    &
         3.49E-01,3.28E-01,3.06E-01,2.85E-01,2.64E-01,2.42E-01,2.37E-01,    &
         2.32E-01,2.22E-01,2.12E-01,2.02E-01,1.92E-01,1.42E-01,1.00E-01,    &
         9.09E-02,7.24E-02,5.39E-02,3.79E-02,2.15E-02,7.26E-03,2.39E-03,    &
         2.39E-03,2.39E-03/                                                 
    DATA n_sil_oc/1.77E+00,1.77E+00,1.77E+00,1.77E+00,1.77E+00,           &
         1.77E+00,1.77E+00,1.77E+00,1.77E+00,1.77E+00,1.77E+00,1.77E+00,    &
         1.77E+00,1.77E+00,1.79E+00,1.82E+00,1.83E+00,1.82E+00,1.81E+00,    &
         1.81E+00,1.81E+00,1.84E+00,1.85E+00,1.86E+00,1.87E+00,1.88E+00,    &
         1.88E+00,1.87E+00,1.86E+00,1.85E+00,1.83E+00,1.80E+00,1.78E+00,    &
         1.74E+00,1.69E+00,1.62E+00,1.51E+00,1.39E+00,1.34E+00,1.36E+00,    &
         1.37E+00,1.39E+00,1.44E+00,1.67E+00,1.86E+00,1.96E+00,2.08E+00,    &
         2.22E+00,2.23E+00,2.17E+00,2.09E+00,2.04E+00,2.01E+00,1.92E+00,    &
         1.78E+00,1.73E+00,1.80E+00,1.98E+00,2.17E+00,2.40E+00,2.48E+00,    &
         2.53E+00,2.58E+00,2.63E+00,2.69E+00,2.70E+00,2.74E+00,2.90E+00,    &
         2.95E+00,2.98E+00,3.01E+00,3.02E+00,3.03E+00,3.03E+00,3.03E+00,    &
         3.04E+00,3.04E+00,3.05E+00,3.05E+00,3.05E+00,3.06E+00,3.06E+00,    &
         3.06E+00,3.06E+00,3.06E+00,3.06E+00,3.06E+00,3.07E+00,3.07E+00,    &
         3.07E+00,3.07E+00,3.08E+00,3.08E+00,3.08E+00,3.08E+00,3.08E+00,    &
         3.08E+00,3.08E+00/                                                 
    DATA k_sil_oc/8.95E-02,8.95E-02,8.95E-02,8.95E-02,8.95E-02,           &
         8.95E-02,8.95E-02,8.95E-02,8.95E-02,8.95E-02,8.95E-02,8.95E-02,    &
         8.95E-02,8.95E-02,8.23E-02,6.64E-02,5.09E-02,4.46E-02,4.69E-02,    &
         6.37E-02,7.89E-02,9.62E-02,1.03E-01,1.00E-01,9.27E-02,8.97E-02,    &
         8.50E-02,8.14E-02,7.85E-02,7.80E-02,7.77E-02,8.06E-02,8.40E-02,    &
         8.91E-02,9.51E-02,1.10E-01,1.48E-01,2.95E-01,4.78E-01,5.49E-01,    &
         6.19E-01,6.96E-01,7.87E-01,9.11E-01,9.10E-01,9.14E-01,9.07E-01,    &
         7.40E-01,5.94E-01,4.74E-01,4.39E-01,4.28E-01,4.20E-01,4.35E-01,    &
         5.35E-01,7.36E-01,9.52E-01,1.10E+00,1.13E+00,1.07E+00,1.04E+00,    &
         1.01E+00,0.97E+00,0.94E+00,8.98E-01,8.77E-01,8.42E-01,7.52E-01,    &
         6.60E-01,5.68E-01,4.75E-01,4.25E-01,4.04E-01,3.82E-01,3.60E-01,    &
         3.39E-01,3.17E-01,2.96E-01,2.74E-01,2.52E-01,2.31E-01,2.26E-01,    &
         2.21E-01,2.12E-01,2.02E-01,1.92E-01,1.83E-01,1.35E-01,9.44E-02,    &
         8.56E-02,6.82E-02,5.07E-02,3.57E-02,2.02E-02,6.82E-03,2.29E-03,    &
         2.29E-03,2.29E-03/                                                 
    DATA n_sil_dl/8.66E-01,8.66E-01,8.66E-01,8.25E-01,7.81E-01,           &
         8.90E-01,1.29E+00,1.60E+00,1.87E+00,2.26E+00,1.93E+00,1.80E+00,    &
         1.76E+00,1.74E+00,1.73E+00,1.72E+00,1.72E+00,1.71E+00,1.71E+00,    &
         1.71E+00,1.71E+00,1.71E+00,1.71E+00,1.71E+00,1.70E+00,1.70E+00,    &
         1.69E+00,1.68E+00,1.66E+00,1.64E+00,1.60E+00,1.57E+00,1.52E+00,    &
         1.47E+00,1.34E+00,1.21E+00,1.16E+00,1.11E+00,1.20E+00,1.24E+00,    &
         1.30E+00,1.35E+00,1.39E+00,1.57E+00,1.75E+00,1.82E+00,1.90E+00,    &
         2.00E+00,2.04E+00,2.09E+00,2.04E+00,2.00E+00,1.91E+00,1.82E+00,    &
         1.73E+00,1.69E+00,1.71E+00,1.78E+00,1.89E+00,2.04E+00,2.11E+00,    &
         2.18E+00,2.26E+00,2.30E+00,2.35E+00,2.40E+00,2.50E+00,2.63E+00,    &
         2.76E+00,2.89E+00,3.02E+00,3.07E+00,3.13E+00,3.18E+00,3.23E+00,    &
         3.25E+00,3.27E+00,3.29E+00,3.31E+00,3.33E+00,3.35E+00,3.35E+00,    &
         3.36E+00,3.37E+00,3.38E+00,3.39E+00,3.40E+00,3.41E+00,3.42E+00,    &
         3.42E+00,3.43E+00,3.43E+00,3.43E+00,3.43E+00,3.43E+00,3.43E+00,    &
         3.43E+00,3.43E+00/                                                 
    DATA k_sil_dl/1.39E-01,1.39E-01,1.39E-01,2.55E-01,4.29E-01,           &
         6.73E-01,8.78E-01,9.26E-01,7.22E-01,5.30E-01,5.32E-02,2.77E-02,    &
         2.84E-02,2.88E-02,2.91E-02,2.94E-02,2.97E-02,3.00E-02,3.03E-02,    &
         3.06E-02,3.09E-02,3.21E-02,3.31E-02,3.39E-02,3.60E-02,3.72E-02,    &
         3.94E-02,4.11E-02,4.25E-02,4.40E-02,4.72E-02,5.05E-02,5.36E-02,    &
         5.66E-02,1.14E-01,1.71E-01,3.68E-01,5.66E-01,7.65E-01,8.29E-01,    &
         8.71E-01,8.97E-01,9.24E-01,9.63E-01,1.00E+00,9.65E-01,9.28E-01,    &
         8.78E-01,7.71E-01,6.63E-01,5.70E-01,4.77E-01,4.83E-01,4.88E-01,    &
         5.83E-01,7.14E-01,8.55E-01,9.76E-01,1.05E+00,1.08E+00,1.10E+00,    &
         1.11E+00,1.13E+00,1.12E+00,1.12E+00,1.12E+00,1.11E+00,1.06E+00,    &
         1.01E+00,9.66E-01,9.19E-01,8.66E-01,8.13E-01,7.59E-01,7.06E-01,    &
         6.72E-01,6.38E-01,6.03E-01,5.69E-01,5.35E-01,5.00E-01,4.83E-01,    &
         4.67E-01,4.33E-01,3.99E-01,3.65E-01,3.32E-01,2.48E-01,2.06E-01,    &
         1.65E-01,1.32E-01,9.87E-02,7.04E-02,3.83E-02,2.46E-02,2.46E-02,    &
         2.46E-02,2.46E-02/                                                 
    DATA n_amc_hn/8.40E-01,8.40E-01,8.40E-01,8.40E-01,7.40E-01,           &
         6.90E-01,9.30E-01,1.53E+00,1.74E+00,1.55E+00,1.22E+00,1.40E+00,    &
         1.60E+00,1.71E+00,1.78E+00,1.85E+00,1.94E+00,2.01E+00,2.11E+00,    &
         2.21E+00,2.31E+00,2.50E+00,2.63E+00,2.70E+00,2.82E+00,2.86E+00,    &
         2.95E+00,3.03E+00,3.04E+00,3.04E+00,3.09E+00,3.15E+00,3.16E+00,    &
         3.16E+00,3.25E+00,3.35E+00,3.38E+00,3.42E+00,3.47E+00,3.49E+00,    &
         3.51E+00,3.53E+00,3.55E+00,3.59E+00,3.63E+00,3.65E+00,3.67E+00,    &
         3.70E+00,3.73E+00,3.75E+00,3.79E+00,3.84E+00,3.86E+00,3.88E+00,    &
         3.96E+00,3.99E+00,4.10E+00,4.14E+00,4.18E+00,4.26E+00,4.31E+00,    &
         4.34E+00,4.36E+00,4.42E+00,4.47E+00,4.50E+00,4.57E+00,4.77E+00,    &
         4.94E+00,5.11E+00,5.25E+00,5.32E+00,5.44E+00,5.49E+00,5.67E+00,    &
         5.71E+00,5.85E+00,5.90E+00,5.99E+00,5.94E+00,6.32E+00,6.70E+00,    &
         6.79E+00,6.99E+00,7.17E+00,7.37E+00,7.55E+00,8.50E+00,9.32E+00,    &
         1.01E+01,1.15E+01,1.25E+01,1.41E+01,1.65E+01,1.65E+01,1.65E+01,    &
         1.65E+01,1.65E+01/                                                 
    DATA k_amc_hn/1.08E-01,1.08E-01,1.08E-01,1.08E-01,1.77E-01,           &
         3.80E-01,9.00E-01,8.40E-01,5.60E-01,1.77E-01,3.21E-01,7.40E-01,    &
         7.20E-01,6.86E-01,6.70E-01,6.95E-01,7.70E-01,8.25E-01,9.00E-01,    &
         9.38E-01,9.60E-01,9.95E-01,1.02E+00,1.01E+00,9.96E-01,9.90E-01,    &
         1.01E+00,1.04E+00,1.04E+00,1.03E+00,1.09E+00,1.15E+00,1.20E+00,    &
         1.25E+00,1.34E+00,1.42E+00,1.44E+00,1.47E+00,1.50E+00,1.51E+00,    &
         1.52E+00,1.53E+00,1.54E+00,1.57E+00,1.60E+00,1.61E+00,1.62E+00,    &
         1.63E+00,1.66E+00,1.69E+00,1.71E+00,1.74E+00,1.76E+00,1.78E+00,    &
         1.83E+00,1.89E+00,1.95E+00,1.95E+00,1.98E+00,2.05E+00,2.08E+00,    &
         2.13E+00,2.19E+00,2.23E+00,2.27E+00,2.32E+00,2.40E+00,2.60E+00,    &
         2.77E+00,2.92E+00,3.00E+00,3.18E+00,3.31E+00,3.50E+00,3.70E+00,    &
         3.80E+00,4.00E+00,4.15E+00,4.30E+00,4.48E+00,4.59E+00,4.70E+00,    &
         4.79E+00,4.97E+00,5.15E+00,5.33E+00,5.51E+00,6.41E+00,7.17E+00,    &
         7.92E+00,9.43E+00,1.04E+01,1.25E+01,1.41E+01,1.41E+01,1.41E+01,    &
         1.41E+01,1.41E+01/                                                 
    DATA n_sic_pg/7.06E-01,7.06E-01,7.06E-01,7.06E-01,7.06E-01,           &
         7.06E-01,7.06E-01,7.06E-01,9.69E-01,2.07E+00,4.29E+00,5.16E+00,    &
         4.75E+00,3.42E+00,2.59E+00,2.55E+00,2.51E+00,2.50E+00,2.47E+00,    &
         2.49E+00,2.50E+00,2.51E+00,2.55E+00,2.58E+00,2.65E+00,2.69E+00,    &
         2.73E+00,2.76E+00,2.76E+00,2.77E+00,2.74E+00,2.71E+00,2.67E+00,    &
         2.63E+00,2.57E+00,2.51E+00,2.42E+00,2.33E+00,2.17E+00,2.11E+00,    &
         2.04E+00,1.96E+00,1.88E+00,1.56E+00,1.25E+00,1.47E+00,1.69E+00,    &
         1.99E+00,3.07E+00,4.14E+00,4.21E+00,4.28E+00,4.11E+00,3.93E+00,    &
         3.71E+00,3.59E+00,3.52E+00,3.49E+00,3.46E+00,3.43E+00,3.42E+00,    &
         3.40E+00,3.39E+00,3.39E+00,3.39E+00,3.38E+00,3.38E+00,3.38E+00,    &
         3.39E+00,3.39E+00,3.39E+00,3.40E+00,3.41E+00,3.41E+00,3.42E+00,    &
         3.42E+00,3.43E+00,3.43E+00,3.44E+00,3.44E+00,3.45E+00,3.45E+00,    &
         3.46E+00,3.46E+00,3.47E+00,3.48E+00,3.49E+00,3.51E+00,3.52E+00,    &
         3.52E+00,3.52E+00,3.52E+00,3.52E+00,3.52E+00,3.52E+00,3.52E+00,    &
         3.52E+00,3.52E+00/                                                 
    DATA k_sic_pg/1.53E+00,1.53E+00,1.53E+00,1.53E+00,1.53E+00,           &
         1.53E+00,1.53E+00,1.53E+00,1.43E+00,1.28E+00,1.06E+00,3.50E-01,    &
         2.50E-01,1.62E-01,1.04E-01,1.13E-01,1.21E-01,1.30E-01,1.44E-01,    &
         1.62E-01,1.80E-01,2.17E-01,2.63E-01,2.76E-01,2.77E-01,2.78E-01,    &
         2.51E-01,2.23E-01,1.88E-01,1.53E-01,1.22E-01,9.08E-02,8.07E-02,    &
         7.06E-02,6.69E-02,6.31E-02,6.59E-02,6.87E-02,9.61E-02,1.05E-01,    &
         1.09E-01,1.24E-01,1.38E-01,5.43E-01,9.48E-01,1.40E+00,1.85E+00,    &
         2.45E+00,2.45E+00,2.45E+00,1.69E+00,9.19E-01,7.00E-01,4.80E-01,    &
         3.62E-01,2.90E-01,2.80E-01,2.82E-01,2.70E-01,2.57E-01,2.51E-01,    &
         2.45E-01,2.38E-01,2.35E-01,2.32E-01,2.29E-01,2.23E-01,2.13E-01,    &
         2.04E-01,1.94E-01,1.84E-01,1.80E-01,1.75E-01,1.71E-01,1.66E-01,    &
         1.63E-01,1.60E-01,1.57E-01,1.54E-01,1.51E-01,1.48E-01,1.46E-01,    &
         1.43E-01,1.39E-01,1.34E-01,1.29E-01,1.25E-01,1.01E-01,9.22E-02,    &
         8.30E-02,8.30E-02,8.30E-02,8.30E-02,8.30E-02,8.30E-02,8.30E-02,    &
         8.30E-02,8.30E-02/                                                 
    DATA n_gr1_dl/9.81E-01,9.25E-01,8.48E-01,7.72E-01,7.23E-01,           &
         7.65E-01,9.91E-01,1.39E+00,2.57E+00,1.94E+00,1.61E+00,1.56E+00,    &
         1.93E+00,2.17E+00,2.35E+00,2.34E+00,2.29E+00,2.25E+00,2.23E+00,    &
         2.21E+00,2.20E+00,2.19E+00,2.18E+00,2.17E+00,2.17E+00,2.16E+00,    &
         2.16E+00,2.15E+00,2.15E+00,2.14E+00,2.13E+00,2.13E+00,2.12E+00,    &
         2.11E+00,2.10E+00,2.09E+00,2.08E+00,2.07E+00,2.06E+00,2.06E+00,    &
         2.05E+00,2.05E+00,2.04E+00,2.03E+00,2.01E+00,1.98E+00,2.05E+00,    &
         2.01E+00,1.99E+00,1.97E+00,1.96E+00,1.95E+00,1.93E+00,1.92E+00,    &
         1.89E+00,1.87E+00,1.84E+00,1.82E+00,1.80E+00,1.76E+00,1.74E+00,    &
         1.73E+00,1.72E+00,1.72E+00,1.71E+00,1.71E+00,1.71E+00,1.76E+00,    &
         1.83E+00,1.92E+00,2.02E+00,2.11E+00,2.21E+00,2.30E+00,2.40E+00,    &
         2.49E+00,2.57E+00,2.66E+00,2.74E+00,2.82E+00,2.90E+00,2.97E+00,    &
         3.05E+00,3.19E+00,3.33E+00,3.46E+00,3.58E+00,4.15E+00,4.65E+00,    &
         5.10E+00,5.90E+00,6.59E+00,7.80E+00,9.33E+00,9.33E+00,9.33E+00,    &
         9.33E+00,9.33E+00/                                                 
    DATA k_gr1_dl/3.08E-03,2.78E-02,9.28E-02,2.09E-01,3.70E-01,           &
         5.39E-01,9.54E-01,1.65E+00,6.55E-01,1.94E-01,2.75E-01,6.38E-01,    &
         8.16E-01,6.57E-01,4.79E-01,2.55E-01,1.41E-01,9.41E-02,6.96E-02,    &
         5.51E-02,4.60E-02,3.33E-02,2.87E-02,2.68E-02,2.42E-02,2.35E-02,    &
         2.32E-02,2.36E-02,2.44E-02,2.57E-02,2.72E-02,2.91E-02,3.13E-02,    &
         3.41E-02,3.74E-02,4.14E-02,4.61E-02,5.16E-02,5.65E-02,5.85E-02,    &
         6.05E-02,6.26E-02,6.47E-02,7.23E-02,8.11E-02,8.83E-02,9.49E-02,    &
         9.92E-02,1.10E-01,1.22E-01,1.34E-01,1.47E-01,1.61E-01,1.75E-01,    &
         2.07E-01,2.42E-01,2.79E-01,3.20E-01,3.63E-01,4.55E-01,5.04E-01,    &
         5.55E-01,6.06E-01,6.59E-01,7.11E-01,7.64E-01,8.68E-01,1.11E+00,    &
         1.33E+00,1.52E+00,1.69E+00,1.85E+00,1.99E+00,2.11E+00,2.23E+00,    &
         2.34E+00,2.44E+00,2.54E+00,2.63E+00,2.72E+00,2.81E+00,2.89E+00,    &
         2.97E+00,3.12E+00,3.27E+00,3.40E+00,3.53E+00,4.12E+00,4.62E+00,    &
         5.08E+00,5.88E+00,6.58E+00,7.79E+00,9.32E+00,9.32E+00,9.32E+00,    &
         9.32E+00,9.32E+00/                                                 
    DATA n_gr2_dl/9.80E-01,9.20E-01,8.39E-01,7.64E-01,5.55E-01,           &
         5.43E-01,1.01E+00,2.08E+00,1.87E+00,1.27E+00,6.88E-01,1.25E+00,    &
         2.39E+00,2.58E+00,2.66E+00,2.74E+00,2.87E+00,3.03E+00,3.19E+00,    &
         3.34E+00,3.47E+00,3.76E+00,3.94E+00,4.05E+00,4.29E+00,4.43E+00,    &
         4.64E+00,4.82E+00,4.99E+00,5.14E+00,5.29E+00,5.43E+00,5.57E+00,    &
         5.69E+00,5.80E+00,5.90E+00,6.01E+00,6.11E+00,6.20E+00,6.22E+00,    &
         6.25E+00,6.28E+00,6.31E+00,6.39E+00,6.47E+00,6.52E+00,6.56E+00,    &
         6.62E+00,6.67E+00,6.74E+00,6.81E+00,6.89E+00,6.96E+00,7.03E+00,    &
         7.15E+00,7.26E+00,7.36E+00,7.46E+00,7.53E+00,7.83E+00,7.99E+00,    &
         8.13E+00,8.25E+00,8.47E+00,8.73E+00,9.01E+00,9.60E+00,1.14E+01,    &
         1.35E+01,1.59E+01,1.88E+01,2.13E+01,2.33E+01,2.48E+01,2.59E+01,    &
         2.66E+01,2.70E+01,2.70E+01,2.66E+01,2.62E+01,2.56E+01,2.46E+01,    &
         2.34E+01,2.05E+01,1.75E+01,1.49E+01,1.33E+01,1.33E+01,1.68E+01,    &
         2.11E+01,3.02E+01,3.89E+01,5.47E+01,7.40E+01,7.40E+01,7.40E+01,    &
         7.40E+01,7.40E+01/                                                 
    DATA k_gr2_dl/3.09E-03,2.79E-02,9.38E-02,1.36E-01,3.62E-01,           &
         7.20E-01,1.48E+00,1.24E+00,3.38E-01,2.77E-01,1.09E+00,2.27E+00,    &
         2.08E+00,1.67E+00,1.56E+00,1.56E+00,1.68E+00,1.82E+00,1.94E+00,    &
         2.03E+00,2.12E+00,2.33E+00,2.48E+00,2.58E+00,2.82E+00,2.96E+00,    &
         3.18E+00,3.39E+00,3.60E+00,3.79E+00,4.00E+00,4.19E+00,4.37E+00,    &
         4.55E+00,4.72E+00,4.91E+00,5.09E+00,5.27E+00,5.41E+00,5.46E+00,    &
         5.50E+00,5.55E+00,5.60E+00,5.77E+00,5.94E+00,6.04E+00,6.14E+00,    &
         6.27E+00,6.44E+00,6.63E+00,6.82E+00,7.00E+00,7.17E+00,7.35E+00,    &
         7.69E+00,8.05E+00,8.42E+00,8.80E+00,9.24E+00,1.01E+01,1.05E+01,    &
         1.10E+01,1.15E+01,1.20E+01,1.24E+01,1.29E+01,1.38E+01,1.59E+01,    &
         1.76E+01,1.90E+01,1.95E+01,1.92E+01,1.83E+01,1.73E+01,1.62E+01,    &
         1.50E+01,1.38E+01,1.26E+01,1.16E+01,1.08E+01,9.98E+00,9.24E+00,    &
         8.73E+00,8.70E+00,9.94E+00,1.24E+01,1.54E+01,2.94E+01,3.95E+01,    &
         4.78E+01,6.08E+01,7.09E+01,8.61E+01,1.03E+02,1.03E+02,1.03E+02,    &
         1.03E+02,1.03E+02/                                                     
    real :: lambda(*)
    integer :: nLambda, nRef
    real, allocatable :: tempIm(:), tempReal(:), lamRef(:)
    real(double) :: dydx
    integer :: i,j 
    real :: t
    character(len=*), intent(in) :: graintype
    character(len=200) :: filename, dataDirectory
    real :: mReal(*), mImg(*)
    dataDirectory = " "
    select case(graintype)
    case("sil_ow")
       allocate(tempIm(1:npLnk))
       allocate(tempReal(1:npLnk))
       allocate(lamRef(1:npLnk))
       lamRef = lam_nk
       nRef = npLnk
       tempReal = n_sil_ow
       tempIm = k_sil_ow

    case("sil_oc")
       allocate(tempIm(1:npLnk))
       allocate(tempReal(1:npLnk))
       allocate(lamRef(1:npLnk))
       lamRef = lam_nk
       nRef = npLnk
       tempReal = n_sil_oc
       tempIm = k_sil_oc

    case("sil_dl")
       allocate(tempIm(1:npLnk))
       allocate(tempReal(1:npLnk))
       allocate(lamRef(1:npLnk))
       lamRef = lam_nk
       nRef = npLnk
       tempReal = n_sil_dl
       tempIm = k_sil_dl

    case("amc_hn")
       allocate(tempIm(1:npLnk))
       allocate(tempReal(1:npLnk))
       allocate(lamRef(1:npLnk))
       lamRef = lam_nk
       nRef = npLnk
       tempReal = n_amc_hn
       tempIm = k_amc_hn

    case("sic_pg")
       allocate(tempIm(1:npLnk))
       allocate(tempReal(1:npLnk))
       allocate(lamRef(1:npLnk))
       lamRef = lam_nk
       nRef = npLnk
       tempReal = n_sic_pg
       tempIm = k_sic_pg

    case("gr1_dl")
       allocate(tempIm(1:npLnk))
       allocate(tempReal(1:npLnk))
       allocate(lamRef(1:npLnk))
       lamRef = lam_nk
       nRef = npLnk
       tempReal = n_gr1_dl
       tempIm = k_gr1_dl

    case("gr2_dl")
       allocate(tempIm(1:npLnk))
       allocate(tempReal(1:npLnk))
       allocate(lamRef(1:npLnk))
       lamRef = lam_nk
       nRef = npLnk
       tempReal = n_gr2_dl
       tempIm = k_gr2_dl

    case("amc_zb")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       filename = trim(dataDirectory)//"/"//"amC-zb2.nk"
       if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
       open(20,file=filename,status="old",form="formatted")
       nRef = 1245
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       do i = 1, nRef
          read(20,*) lamRef(i), tempReal(i), tempIm(i)
       enddo
       close(20)

    case("draine_sil")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       filename = trim(dataDirectory)//"/"//"Draine_Si_sUV.dat"
       if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
       open(20,file=filename,status="old",form="formatted")
       nRef = 1201
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       do i = nRef, 1, -1
          read(20,*) lamRef(i), tempReal(i), tempIm(i)
       enddo
       close(20)

    case("pinteISM")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       filename = trim(dataDirectory)//"/"//"pinteISM.dust"
       if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
       open(20,file=filename,status="old",form="formatted")
       nRef = 38
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       do i = nRef, 1, -1
          read(20,*) lamRef(i), tempReal(i), tempIm(i)
       enddo
       close(20)

    case DEFAULT
       if (writeoutput) write(*,'(a,a,a)') "! Grain type ", trim(graintype)," not recognised"
       stop
    end select

    do i = 1, nLambda 
       if (lambda(i)*real(angsToMicrons) <= lamRef(nRef)) then
          call locate(lamRef, nRef, lambda(i)*real(angsToMicrons), j)
          t = (lambda(i)*angsToMicrons - lamRef(j))/(lamRef(j+1) - lamRef(j))
          mReal(i) = tempReal(j) + t * (tempReal(j+1) - tempReal(j))
          mImg(i) = tempIm(j) + t * (tempIm(j+1) - tempIm(j))         
       else
          call writeWarning("Extrapolating grain properties")
          dydx = (log10(tempReal(nref)) - log10(tempReal(nRef-1))) / &
               (log10(lamRef(nref))-log10(lamRef(nRef-1)))
          mReal(i) = log10(tempReal(nref)) + dydx * &
               (log10(lambda(i)*angsToMicrons) - log10(lamRef(nRef)))
          mReal(i) = 10.d0**mreal(i)
          dydx = (log10(tempIm(nref)) - log10(tempIm(nRef-1))) / &
               (log10(lamRef(nref))-log10(lamRef(nRef-1)))
          mImg(i) = log10(tempIm(nref)) + dydx * &
               (log10(lambda(i)*angsToMicrons) - log10(lamRef(nRef)))
          mImg(i) = 10.d0**mImg(i)
       endif

    enddo

  end subroutine getRefractiveIndex

  subroutine fillGridMie(grid, scale, aMin, aMax, a0, qDist, pDist, &
       ngrain, abundance, grainname, thisDust)

    implicit none
    type(GRIDTYPE) :: grid
    integer :: thisDust
    real :: aMin, aMax,a0, qDist, pDist
    real, allocatable :: sigmaAbs(:), sigmaSca(:), sigmaExt(:)
    real :: scale
    real, allocatable :: mReal(:), mImg(:)          ! size = nlamda
    real, allocatable :: mReal2D(:,:), mImg2D(:,:)  ! size = ngrain x nlambda
    real :: meanParticleMass
    real :: getMeanMass2
    real :: rayleigh, gsca
    external getMeanMass2
    integer, intent(in) :: ngrain  ! number of grain types
    real, intent(in) :: abundance(:)   ! relative abundance of grains
    character(len=*) :: grainname(:)   ! names of grains available
    real :: sig_ext, sig_scat, sig_abs
    real :: total_abundance
    character(len=80) :: albedoFilename

    integer :: i, j, k


    scale = 1.

    allocate(sigmaAbs(1:grid%nLambda))
    allocate(sigmaSca(1:grid%nLambda))
    allocate(sigmaExt(1:grid%nLambda))


    call writeInfo("NEW: Filling grid with mie cross-sections...", TRIVIAL)
    call writeInfo("Dust law: n(a) = const * a^-q * Exp( -(a/a0)^p )", TRIVIAL)
    call writeInfo("          where  amin < a < amax", TRIVIAL)
    call writeFormatted("(a,f10.3)","    amin  = ",  aMin, TRIVIAL)
    call writeFormatted("(a,f10.3)","    amax  = ",  aMax, TRIVIAL)
    call writeFormatted("(a,e12.3)","    a0    = ",  a0, TRIVIAL)
    call writeFormatted("(a,e12.3)","    qDist = ",  qDist, TRIVIAL)
    call writeFormatted("(a,e12.3)","    pDist = ",  pDist, TRIVIAL)


    allocate(mReal(1:grid%nLambda))
    allocate(mImg(1:grid%nLambda))

    ! Synthetic grains

    ! quick test for zero total abundance.
    total_abundance = SUM(abundance)
    if ( total_abundance <= 0.0 ) then
       if (writeoutput) write(*,*) "Error:: total_abundance <= 0.0 in  grain_mod::fillGridMie."
       if (writeoutput) write(*,*) "  ==> You probably forgot to assign abundance in your parameter file!"
       if (writeoutput) write(*,*) "  ==> Exiting the prograim ... "
       stop 
    end if

    ! allocate mem for temp arrays
    allocate(mReal2D(1:ngrain, 1:grid%nLambda))
    allocate(mImg2D(1:ngrain, 1:grid%nLambda))
    ! initializing the values
    mReal2D(:,:) = 0.0; mImg2D(:,:) = 0.0

    ! Find the index of refractions for all types of grains available
    do j = 1, ngrain
       call getRefractiveIndex(grid%lamArray, grid%nLambda, grainname(j), mReal, mImg)
       mReal2D(j,:) = mReal(:)  ! copying the values to a 2D maxtrix
       mImg2D(j,:)  = mImg(:)   ! copying the values to a 2D maxtrix            
    end do

    ! finding the cross sections
    sigmaExt(:) = 0.0; sigmaAbs(:)=0.0; sigmaSca(:)=0.0 ! initializing the values

!    if (writeoutput) open(20,file="albedo.dat",form="formatted",status="unknown")
!    if (writeoutput) open(21,file="gfactor.dat",form="formatted",status="unknown")
    do i = 1, grid%nLambda
       do j = 1, ngrain
          call mieDistCrossSection(aMin, aMax, a0, qDist, pDist, grid%lamArray(i), &
               mReal2D(j,i), mImg2D(j,i), sig_ext, sig_scat, sig_abs, gsca)

          ! Weighting the cross section according to their abundance...            
          sigmaExt(i) = sig_ext*abundance(j)+ sigmaExt(i)
          sigmaAbs(i) = sig_abs*abundance(j)+ sigmaAbs(i)
          sigmaSca(i) = sig_scat*abundance(j)+ sigmaSca(i)
       end do
       sigmaExt(i) =    sigmaExt(i)/total_abundance 
       sigmaAbs(i) =    sigmaAbs(i)/total_abundance 
       sigmaSca(i) =    sigmaSca(i)/total_abundance 

!       if (writeoutput) write(21,*) grid%lamArray(i), gsca
!       if (writeoutput) write(20,*) grid%lamArray(i),sigmaExt(i),sigmaAbs(i),sigmaSca(i),sigmaSca(i)/sigmaExt(i)
    end do

!    if (writeoutput) close(20)
!    if (writeoutput) close(21)

    if (.not.grid%oneKappa) then
       if (grid%adaptive) then
          if (writeoutput) write(*,'(a,i3)') "Filling AMR grid with mie cross sections...",grid%nLambda
          call fillAMRgridMie(grid%OctreeRoot, sigmaSca, sigmaAbs, grid%nLambda)
       endif

       if (grid%cartesian) then

          do i = 1, grid%nx
             do j = 1, grid%ny
                do k = 1, grid%nz


                   if (grid%inUse(i,j,k)) then
                      grid%kappaAbs(i,j,k,1:grid%nLambda) = sigmaAbs  * grid%rho(i,j,k)
                      grid%kappaSca(i,j,k,1:grid%nLambda) = sigmaSca  * grid%rho(i,j,k)
                   endif

                   !      write(*,*) grid%kappaAbs(i,j,k,1:grid%nLambda),grid%kappaSca(i,j,k,1:grid%nLambda), grid%rho(i,j,k),scale
                enddo
             enddo
          enddo
       endif

       if (grid%polar) then
          do i = 1, grid%nr
             do j = 1, grid%nmu
                do k = 1, grid%nPhi

                   grid%kappaAbs(i,j,k,1:grid%nLambda) = sigmaAbs * grid%rho(i,j,k)
                   grid%kappaSca(i,j,k,1:grid%nLambda) = sigmaSca * grid%rho(i,j,k)

                enddo
             enddo
          enddo

       endif

       where(grid%kappaAbs < 1.e-25) grid%kappaAbs = 1.e-25
       where(grid%kappaSca < 1.e-25) grid%kappaSca = 1.e-25


       grid%kappaAbs = grid%kappaAbs * 1.e10
       grid%kappaSca = grid%kappaSca * 1.e10
    else
       call writeFormatted("(a,i4)", "Filling the oneKappa arrays: ",grid%nLambda, TRIVIAL)

       meanParticleMass = 0.
       do i = 1, ngrain
          meanParticleMass = meanParticleMass + getMeanMass2(aMin, aMax, a0, qDist, pDist, grainname(i))*abundance(i)
       enddo
       grid%oneKappaAbs(thisDust,1:grid%nLambda) = (sigmaAbs(1:grid%nLambda) * 1.e10)/meanParticleMass
       grid%oneKappaSca(thisDust,1:grid%nLambda) = (sigmaSca(1:grid%nLambda) * 1.e10)/meanParticleMass

       write(albedoFilename,'(a,i2.2,a)') "albedo",thisDust,".dat"
       if (writeoutput) open(20,file=albedoFilename,form="formatted",status="unknown")
       if (writeoutput) write(20,'(a120)') "# Columns are: wavelength (microns), kappa ext (cm^2 g^-1), &
           &    kappa abs (cm^2 g^-1), kappa sca (cm^2 g^-1), albedo"
       if (writeoutput) write(20,*) "# Note that the opacities are per gram of dust"

       do i = 1, grid%nLambda
          rayleigh = (8.*pi**2)/(grid%lamArray(i)*angstromtocm)* &
               aimag((cmplx(mreal(i),mimg(i))**2-cmplx(1.,0.))/(cmplx(mreal(i),mimg(i))**2+cmplx(2.,0.)))*(amin*microntocm)**3
          rayleigh = rayleigh / meanParticleMass
          if (writeoutput) &
               write(20,*) grid%lamArray(i)*angstomicrons,&
               (grid%oneKappaAbs(thisdust,i)+grid%oneKappaSca(thisdust,i))/1.e10, &
               grid%oneKappaAbs(thisdust,i)/1.e10,grid%oneKappaSca(thisdust,i)/1.e10, &
               grid%oneKappaSca(thisdust,i)/(grid%oneKappaAbs(thisdust,i)+grid%oneKappaSca(thisdust,i))
       enddo
       if (writeoutput) close(20)




    endif
    deallocate(sigmaAbs, sigmaSca)
    call writeInfo("mie cross-sections done. Note 10^10 factor", TRIVIAL)
  end subroutine fillGridMie

  subroutine setKappaTest(grid, scale, aMin, aMax, a0, qDist, pDist, grainType, &
       ngrain, abundance, grainname, lambdaTau)

    implicit none
    type(GRIDTYPE) :: grid
    real :: aMin, aMax, a0, qDist, pDist
    real :: sigmaAbs, sigmaSca, sigmaExt
    real :: scale
    real, allocatable :: mReal(:), mImg(:)
    real, allocatable :: mReal2D(:,:), mImg2D(:,:)  ! size = ngrain x nlambda
    character(len=*) :: grainType
    integer, intent(in) :: ngrain  ! number of grain types
    real, intent(in) :: abundance(ngrain)   ! relative abundance of grains
    character(len=*) :: grainname(ngrain)   ! names of grains available
    real :: sig_ext, sig_scat, sig_abs
    real :: total_abundance, gsca
    real :: meanParticleMass
    real :: lambdaTau
    real :: getMeanMass2

    integer :: i, j 

    scale = 1.
    allocate(mReal(1:grid%nLambda))
    allocate(mImg(1:grid%nLambda))
    call locate(grid%lamArray, grid%nLambda, lambdaTau, i)


    if (writeoutput) write(*,*) "kappa test set for: ",grid%lamarray(i)


    scale = 1.
    if (writeoutput) write(*,'(a)') "NEW: Filling grid with mie cross-sections..."

    if (graintype(1:5) == "mixed") then
       ! Synthetic grains

       ! quick test for zero total abundance.
       total_abundance = SUM(abundance)
       if ( total_abundance <= 0.0 ) then
          write(*,*) "Error:: total_abundance <= 0.0 in  grain_mod::setKappaTest."
          write(*,*) "  ==> You probably forgot to assign abundance in your parameter file!"
          write(*,*) "  ==> Exiting the prograim ... "
          stop 
       end if

       ! allocate mem for temp arrays
       allocate(mReal2D(1:ngrain, 1:grid%nLambda))
       allocate(mImg2D(1:ngrain, 1:grid%nLambda))
       ! initializing the values
       mReal2D(:,:) = 0.0; mImg2D(:,:) = 0.0

       ! Find the index of refractions for all types of grains available
       do j = 1, ngrain
          call getRefractiveIndex(grid%lamArray, grid%nLambda, grainname(j), mReal, mImg)
          mReal2D(j,:) = mReal(:)  ! copying the values to a 2D maxtrix
          mImg2D(j,:)  = mImg(:)   ! copying the values to a 2D maxtrix            
       end do

       ! finding the cross sections
       sigmaExt = 0.0; sigmaAbs=0.0; sigmaSca=0.0 ! initializing the values

       do j = 1, ngrain
          call mieDistCrossSection(aMin, aMax, a0, qDist, pDist, grid%lamArray(i), &
               mReal2D(j,i), mImg2D(j,i), sig_ext, sig_scat, sig_abs ,gsca)

          ! Weighting the cross section according to their abundance...            
          sigmaExt = sig_ext*abundance(j)+ sigmaExt
          sigmaAbs = sig_abs*abundance(j)+ sigmaAbs
          sigmaSca = sig_scat*abundance(j)+ sigmaSca
       end do
       sigmaExt =    sigmaExt/total_abundance 
       sigmaAbs =    sigmaAbs/total_abundance 
       sigmaSca =    sigmaSca/total_abundance 

    else 
       ! Do a single grain calculations...       
       call getRefractiveIndex(grid%lamArray, grid%nLambda, graintype, mReal, mImg)         
       call mieDistCrossSection(aMin, aMax, a0, qDist, pDist,grid%lamArray(i),  &
            mReal(i), mImg(i), sigmaExt, sigmaSca, sigmaAbs, gsca)
    end if
    meanParticleMass = getMeanMass2(aMin, aMax, a0, qDist, pDist, graintype)

    grid%kappaTest = sigmaExt * 1.e10 / meanParticleMass

    if (writeoutput) write(*,*) "kappa test is: ",grid%kappatest
  end subroutine setKappaTest

  subroutine setKappa(kappaAbs, kappaSca, lambda, nLambda, aMin, aMax, a0, qDist, pDist, grainType)

    implicit none
    real :: aMin, aMax, qDist, a0, pDist
    real :: sigmaAbs, sigmaSca, sigmaExt
    real, allocatable :: mReal(:), mImg(:)
    character(len=*) :: grainType
    real :: lambda(:), kappaAbs(:), kappaSca(:), gSca
    integer :: nLambda
    integer :: i

    if (writeoutput) write(*,'(a)') "Setting kappas for: ",trim(grainType)

    allocate(mReal(1:nLambda))
    allocate(mImg(1:nLambda))
    call getRefractiveIndex(lambda, nLambda, graintype, mReal, mImg)
    do i = 1, nLambda
       call mieDistCrossSection(aMin, aMax, a0, qDist, pDist, lambda(i),  mReal(i), mImg(i), sigmaExt, &
            sigmaSca, sigmaAbs, gSca)
       kappaAbs(i) = sigmaAbs * 1.e10
       kappaSca(i) = sigmaSca * 1.e10
    enddo
  end subroutine setKappa

  !
  ! Computes Kappa at a single fgrequency/wavelength.
  !
  subroutine MieCrossSection(sigmaExt, sigmaAbs, sigmaSca, &
       aMin, aMax, a0, qDist, pDist, grainType, &
       ngrain, abundance, grainname, lambda)

    implicit none
    real, intent(out) :: sigmaExt, sigmaAbs, sigmaSca ! total, absorption and scattering  x-sections
    real, intent(in) :: aMin, aMax, a0, qDist, pDist
    character(len=*), intent(in) :: grainType
    integer, intent(in) :: ngrain  ! number of grain types
    real, intent(in) :: abundance(ngrain)   ! relative abundance of grains
    character(len=*) :: grainname(ngrain)   ! names of grains available
    real, intent(in) :: lambda  ! wavelength at which sigma are cmoputed
    !
    real :: mReal(ngrain), mImg(ngrain)  ! automatic arrays
    !
    real :: sig_ext, sig_scat, sig_abs
    real :: total_abundance, gsca
!    real :: meanParticleMass
!    real :: getMeanMass2
    integer, parameter :: nlambda = 1
    integer :: j 
    ! dummy array for interface with subroutines
    real :: lamArray(nlambda)   
    real :: mRealArray(nlambda), mImgArray(nlambda)  ! automatic arrays


    lamArray(:) = lambda  ! same for all elemets
    if (graintype(1:5) == "mixed") then
       ! Synthetic grains

       ! quick test for zero total abundance.
       total_abundance = SUM(abundance)
       if ( total_abundance <= 0.0 ) then
          write(*,*) "Error:: total_abundance <= 0.0 in  dust_mod::getKappa."
          write(*,*) "  ==> You probably forgot to assign abundance in your parameter file!"
          write(*,*) "  ==> Exiting the prograim ... "
          stop 
       end if

       ! initializing the values
       mReal(:) = 0.0; mImg(:) = 0.0

       ! Find the index of refractions for all types of grains available
       do j = 1, ngrain
          call getRefractiveIndex(lamArray, nLambda, grainname(j), mRealArray, mImgArray)
          mReal(j) = mRealArray(1) 
          mImg(j)  = mImgArray(1)  
       end do

       ! finding the cross sections
       sigmaExt = 0.0; sigmaAbs=0.0; sigmaSca=0.0 ! initializing the values

       do j = 1, ngrain
          call mieDistCrossSection(aMin, aMax, a0, qDist, pDist, lamArray(1), &
               mReal(j), mImg(j), sig_ext, sig_scat, sig_abs ,gsca)

          ! Weighting the cross section according to their abundance...            
          sigmaExt = sig_ext*abundance(j)+ sigmaExt
          sigmaAbs = sig_abs*abundance(j)+ sigmaAbs
          sigmaSca = sig_scat*abundance(j)+ sigmaSca
       end do
       sigmaExt =    sigmaExt/total_abundance 
       sigmaAbs =    sigmaAbs/total_abundance 
       sigmaSca =    sigmaSca/total_abundance 

    else 
       ! Do a single grain calculations... 
       call getRefractiveIndex(lamArray, nLambda, graintype, mRealArray, mImgArray)         
       call mieDistCrossSection(aMin, aMax, a0, qDist, pDist,lamArray(1),  &
            mRealArray(1), mImgArray(1), sigmaExt, sigmaSca, sigmaAbs, gsca)
    end if

    !      meanParticleMass = getMeanMass2(aMin, aMax, a0, qDist, pDist, graintype)
  end subroutine MieCrossSection

  recursive subroutine fillAMRgridMie(thisOctal, sigmaSca, sigmaAbs, nLambda)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: nLambda
    real :: sigmaSca(*), sigmaAbs(*)
    integer :: subcell, i
!    write(*,*) subcell,sigmasca(1),sigmaabs(1),nlambda

    do subcell = 1, 8
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fillAMRgridMie(child, sigmaSca, sigmaAbs, nLambda)
                exit
             end if
          end do
       else
          if (thisOctal%inFlow(subcell)) then
             thisOctal%kappaAbs(subcell,1:nLambda)   &
                  = sigmaAbs(1:nLambda) * thisOctal%rho(subcell) * 1.e10
             thisOctal%kappaSca(subcell,1:nLambda)   &
                  = sigmaSca(1:nLambda) * thisOctal%rho(subcell) * 1.e10
          else
             thisOctal%kappaAbs(subcell,1:nLambda)   &
                  = 0.0
             thisOctal%kappaSca(subcell,1:nLambda)   &
                  = 0.0
          end if

       endif
    enddo
  end subroutine fillAMRgridMie

  subroutine dustPropertiesfromFile(filename, nlambda, lambda, kappaAbs, kappaSca)
    implicit none
    character(len=*) :: filename
    integer :: nlambda
    real :: lambda(*)
    real :: kappaAbs(*), kappaSca(*)
    real :: sigmaExt(1000),sigmaSca(1000), kappa(1000), albedo(1000), tlam(1000)
    real :: tSca(1000), tAbs(1000), sigmaAbs(1000)
    character(len=40) :: filetype
    character(len=80) :: message
    integer :: npts, i, j

    write(message,'(a,a)') "Reading dust properties from: ",trim(filename)
    call writeInfo(message, TRIVIAL)
    open(20, file=filename, status="old", form="formatted")
    read(20,'(a)') filetype

    select case (filetype)

    case("kenny")
       npts = 1
10     read(20,*,end=20) tlam(npts), sigmaExt(npts),sigmaSca(npts),kappa(npts)
       npts = npts + 1
       goto 10
20     continue
       npts = npts - 1
       close(20)
       albedo(1:npts) = sigmaSca(1:npts) / sigmaExt(1:npts)

       tAbs(1:npts) = (1.-albedo(1:npts))*kappa(1:npts)
       tSca(1:npts) = albedo(1:npts)*kappa(1:npts)

       tlam(1:npts) = tlam(1:npts) * 1.e4 ! microns to angstrom

    case("jenny")
       npts = 1
30     read(20,*,end=40) tlam(npts), sigmaAbs(npts),sigmaSca(npts)
       npts = npts + 1
       goto 30
40     continue
       npts = npts - 1
       close(20)
       tAbs(1:npts) = sigmaAbs(1:npts)
       tSca(1:npts) = sigmaSca(1:npts)
       tlam(1:npts) = tlam(1:npts) * 1.e4 ! microns to angstrom

    case DEFAULT
       write(*,'(a)') "! Dust properties file has unknown type",trim(filetype)
       stop

    end select

    do i = 1, nLambda
       call locate(tlam,npts,lambda(i),j)
       kappaAbs(i) = logint(lambda(i), tlam(j), tlam(j+1), tAbs(j), tAbs(j+1))*1.e10
       kappaSca(i) = logint(lambda(i), tlam(j), tlam(j+1), tSca(j), tSca(j+1))*1.e10
    enddo

    !    write(*,*) "Correcting xsections by dust-to-gas ratio!!!!!!!!"
    !    kappaAbs(1:nLambda) = kappaAbs(1:nLambda) * 0.01
    !    kappaSca(1:nLambda) = kappaSca(1:nLambda) * 0.01

  end subroutine dustPropertiesfromFile

  recursive subroutine fillDustShakara(grid, thisOctal)

    use input_variables, only : rInner, rOuter
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    type(VECTOR) :: rVec
    real :: x, z
    real :: height
    real(double) :: fac
    integer :: subcell, i
    height = 0.d0

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fillDustShakara(grid, child)
                exit
             end if
          end do
       else

          thisOctal%DustTypeFraction(subcell,1:2) = 0.d0
          thisOctal%DustTypeFraction(subcell,1) = 1.d0
          rVec = subcellCentre(thisOctal, subcell)
          x = rVec%x
          z = rVec%z
          if ( (x > rInner).and.(x < rOuter)) then
             call returnScaleHeight(grid, x, height)
             fac = exp(-abs(z/height))
             thisOctal%dustTypeFraction(subcell,1) = 1.d0 - fac
             thisOctal%dustTypeFraction(subcell,2) = fac
          endif

       end if
    end do

  end subroutine fillDustShakara





  recursive subroutine fillDustUniform(grid, thisOctal)

    use input_variables, only : maxDustTypes, nDustType, grainFrac
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fillDustUniform(grid, child)
                exit
             end if
          end do
       else
          thisOctal%dustTypeFraction(subcell,1:nDustType) = grainFrac(1:nDustType)
       end if
    end do

  end subroutine fillDustUniform


  recursive subroutine sublimateDustWR104(thisOctal)
    use input_variables, only : tThresh
    type(OCTAL), pointer :: thisOctal, child
    integer :: subcell, i
    real(double) :: temperature
    real(double), parameter :: subRange = 1.d-3
    real(double) :: frac
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call sublimateDustWR104(child)
                exit
             end if
          end do
       else
          temperature = thisOctal%temperature(subcell)

          frac = 1.d0
          if (temperature > tThresh) then
             frac = exp(-dble((temperature-tThresh)/subRange))
          endif
          thisOctal%dustTypeFraction(subcell,1) = frac
       endif
       thisOctal%dustTypeFraction(subcell,1) = max(1.d-20, frac)
    end do
  end subroutine sublimateDustWR104

  recursive subroutine sublimateDust(grid, thisOctal, totFrac, nFrac, tauMax)

    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real :: totFrac
    real :: tauMax
    integer :: nFrac
    real(double) :: frac, newFrac, deltaFrac, normFac, thistau
    real ::  temperature, sublimationTemp, subrange
    real :: underCorrect = 1.
    integer :: ilambda
    real(double) :: kappaSca, kappaAbs
    integer :: subcell, i

    kappaSca = 0.d0; kappaAbs = 0.d0
    subrange = 1.

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call sublimateDust(grid, child, totFrac,nFrac, tauMax)
                exit
             end if
          end do
       else

          temperature = thisOctal%temperature(subcell)
          sublimationTemp = max(700.d0,2000.d0 * thisOctal%rho(subcell)**(1.95d-2))
          if (temperature < sublimationTemp) newFrac = 1.

          if (temperature >= sublimationTemp) then
             newfrac = exp(-dble((temperature-sublimationtemp)/subRange))
          endif
          newfrac = max(newfrac,1.d-20)

          deltaFrac = newFrac - thisOctal%oldFrac(subcell)

          frac = thisOctal%oldFrac(subcell) + underCorrect * deltaFrac
          frac = max(frac, 1.d-20)

          normfac = SUM(thisOctal%dustTypeFraction(subcell,:))


          if (normFac /= 0.d0) then
             thisOctal%dustTypeFraction(subcell,:) = thisOctal%dustTypeFraction(subcell,:) * frac /normFac
          else
             thisOctal%dustTypeFraction(subcell,:) = 0.d0
             thisOctal%dustTypeFraction(subcell,1) = frac
          end if


          call locate(grid%lamArray, grid%nLambda, 5500., iLambda)
          call returnKappa(grid, thisOctal, subcell, ilambda, kappaSca=kappaSca, kappaAbs=kappaAbs)

          thisTau = (kappaAbs+kappaSca)*thisOctal%subcellSize
          if (thisTau > tauMax) then
             frac = frac *tauMax / thisTau 
             normfac = SUM(thisOctal%dustTypeFraction(subcell,:))
             thisOctal%dustTypeFraction(subcell,:) = thisOctal%dustTypeFraction(subcell,:) * frac /normFac
          endif


          if (deltaFrac /= 0.) then
             nfrac = nfrac + 1
             totFrac = totFrac + abs(deltaFrac)
          endif

          thisOctal%oldFrac(subcell) = frac

       end if
    end do

  end subroutine sublimateDust

  subroutine returnScaleHeight(grid, x,  height)
    type(gridtype) :: grid
    integer :: nz
    real :: x, z(10000), height
    real :: subcellsize(10000), temperature(10000)
    real(double) :: rho(10000)
    real(double) :: rho_over_e
    integer :: i, j
    nz  = 0
    temperature = 0.; subcellSize = 0.; rho = 0.0; z=0.d0
    call getTemperatureDensityRundust(grid, z, subcellsize, rho, temperature, x, 0., nz, 1.)


    if (rho(2) >= rho(1)) then
       height = 1.e20
       goto 666
    endif

    rho_over_e = rho(1) / exp(1.d0)

    if (rho_over_e < rho(nz)) then
       height = z(nz) /1.e10
    else
       j = 1
       do i = 1, nz
          if ((rho_over_e < rho(i)).and.(rho_over_e > rho(i+1))) then
             j = i
             exit
          endif
       enddo
       height = z(j) + (z(j+1)-z(j))*(rho_over_e - rho(j))/(rho(j+1)-rho(i))
       height = height / 1.e10
    endif
666 continue
  end subroutine returnScaleHeight

  subroutine getTemperatureDensityRundust(grid, zAxis, subcellsize, rho, temperature, xPos, yPos, nz, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    integer :: nz
    real(double) :: rho(:)
    real ::temperature(:), zAxis(:), subcellsize(:)
    real :: xPos, yPos
    integer :: subcell
    real(double) :: rhotemp
    real :: temptemp
    real :: direction
    type(VECTOR) :: currentPos, temp
    real :: halfSmallestSubcell

    nz = 0
    halfSmallestSubcell = grid%halfSmallestSubcell

    currentPos = VECTOR(xPos, yPos, direction*halfSmallestSubcell)

    do while(abs(currentPos%z) < grid%ocTreeRoot%subcellsize)
       call amrGridValues(grid%octreeRoot, currentPos, foundOctal=thisOctal, &
            foundSubcell=subcell, rho=rhotemp, temperature=temptemp)
       thisOctal%chiLine(subcell) = 1.e-30
       !       if (thisOctal%inFlow(subcell)) then
       nz = nz + 1
       temperature(nz) = temptemp
       rho(nz) = rhotemp
       temp = subCellCentre(thisOctal, subcell)
       zAxis(nz) = temp%z
       subcellsize(nz) = thisOctal%subcellsize
       !       endif
       currentPos = VECTOR(xPos, yPos, zAxis(nz)+0.5*direction*thisOctal%subcellsize+direction*halfSmallestSubcell)
       !       else
       !          currentPos = VECTOR(xPos, yPos, grid%octreeRoot%subcellsize+halfSmallestSubcell)
       !       endif
    end do
    zAxis(1:nz) = abs(zAxis(1:nz)) * 1.e10  ! convert to cm
  end subroutine getTemperatureDensityRundust

  recursive subroutine stripDustAway(thisOctal, fac, radius)

    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real(double) :: fac
    real(double) :: radius, r
    integer :: subcell, i


    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call stripDustAway(child, fac, radius)
                exit
             end if
          end do
       else
          r = modulus(subcellCentre(thisOctal, subcell))
          if (r < radius) then
             thisOctal%dustTypeFraction(subcell,:) = thisOctal%dustTypeFraction(subcell,:) * fac
             thisOctal%etaCont(subcell) = tiny(thisOctal%etaCont(subcell))
             thisOctal%oldFrac(subcell) = fac
          endif
       end if
    end do

  end subroutine stripDustAway

  subroutine createRossArray(grid)
    use input_variables, only : nDustType
    type(GRIDTYPE) :: grid
    integer :: i, j, k
    real :: bNuTot, rosselandKappa, temperature
    real(double) :: dFreq, Freq
    real :: maxTemp

    maxTemp = 3000.

    if (grid%geometry == "wr104") maxTemp = 30000.

    if (grid%nTempRossArray == 0) then
       call writeFatal("nTempRossArray is zero on call to createRossArray")
    endif

    do k = 1, grid%nTempRossArray
       temperature = 3. + (maxTemp-3.)*real(k-1)/real(grid%nTempRossArray-1)
       do j = 1, nDustType
          rosselandKappa = 0.
          Bnutot = 0.
          do i =  grid%nLambda,2,-1
             freq = cSpeed / (grid%lamArray(i)*1.e-8)
             dfreq = cSpeed / (grid%lamArray(i)*1.e-8) - cSpeed / (grid%lamArray(i-1)*1.e-8)
!             rosselandKappa = rosselandKappa + bnu(freq, dble(temperature)) * dFreq / &
!                  (grid%oneKappaabs(j,i)+grid%oneKappaSca(j,i))
!             bnutot = bnutot + bnu(freq, dble(temperature)) * dfreq

             rosselandKappa = rosselandKappa + dbNubydT(freq, dble(temperature)) * dFreq / &
                  (grid%oneKappaabs(j,i)+grid%oneKappaSca(j,i))
             bnutot = bnutot + dbNubydT(freq, dble(temperature))*dfreq

          enddo
          if (rosselandkappa /= 0.) then
             rosselandKappa = (bnutot / rosselandKappa)/1.d10
          endif
          grid%kappaRossArray(j,k) = rosselandKappa
       enddo
       grid%tempRossArray(k) = temperature
    enddo
  end subroutine createRossArray



  recursive subroutine fillDustWhitney(grid, thisOctal)
    use input_variables
    TYPE(VECTOR) :: point
    TYPE(gridtype), INTENT(IN)    :: grid
    real :: r, mu, mu_0, muCavity, rhoEnv, r_c
    real :: h, rhoDisc, alpha
    real(double) :: fac
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fillDustWhitney(grid, child)
                exit
             end if
          end do
       else

          point = subcellCentre(thisOctal, subcell)
          muCavity = cos(cavAngle/2.)
          r = modulus(point)*1.e10

          mu = (point%z*1.e10) /r

          r_c = erInner
          alpha = 2.25
          beta = 1.25

          rhoEnv = cavdens * mHydrogen

          ! by default use envelope grains

          thisOctal%dustTypeFraction(subcell,1:4) = 0.d0
          thisOctal%dustTypeFraction(subcell,3) = 1.d0


          if ((r > erInner).and.(r < erOuter)) then
             mu_0 = rtnewtdust(-0.2 , 1.5 , 1.e-4, r/r_c, abs(mu))
             ! equation 1 for Whitney 2003 ApJ 591 1049 has a mistake in it
             ! this is from Momose et al. 1998 ApJ 504 314

             rhoEnv = (mdotenv / fourPi) * (bigG * mCore)**(-0.5) * r**(-1.5) * &
                  (1. + abs(mu)/mu_0)**(-0.5) * &
                  (abs(mu)/mu_0 + (2.*mu_0**2 * r_c/r))**(-1.)

             fac =  1.d0-min(dble(r - erInner)/(0.02d0*erinner),1.d0)
             fac = exp(-fac*10.d0)
             rhoEnv = rhoEnv * fac
             rhoEnv = max(rhoEnv, tiny(rhoEnv))
          endif

          if (mu_0 > muCavity) then
             !       fac =  1.d0-min(dble(r - drInner)/(0.02d0*erinner),1.d0)
             !       fac = exp(-fac*10.d0)
             rhoEnv = cavdens * 2.*mHydrogen

             ! outflow cavity grains

             thisOctal%dustTypeFraction(subcell,1:4) = 0.d0
             thisOctal%dustTypeFraction(subcell,4) = 1.d0
          endif

          if (r < erInner) then
             thisOctal%dustTypeFraction(subcell,1:4) = 0.d0
             thisOctal%dustTypeFraction(subcell,4) = 1.d0
          endif

          rho0  = mDisc *(beta-alpha+2.) / ( twoPi**1.5 * 0.01*rStellar * rStellar**(alpha-beta) * ( &
               (drouter**(beta-alpha+2.)-drInner**(beta-alpha+2.))) )

          r = sqrt(point%x**2 + point%y**2)*1.e10
          h = 0.01 * rStellar * (r/rStellar)**beta
          rhoDisc = 1.e-30
          if ((r > drInner).and.(r < drOuter)) then
             rhoDisc = rho0 * (rStellar/r)**alpha  * exp(-0.5*((point%z*1.e10)/h)**2)
             fac =  1.d0-min(dble(r - drInner)/(0.02d0*drinner),1.d0)
             fac = exp(-fac*10.d0)
             rhoDisc = rhoDisc * fac
             rhoDisc = max(rhoDisc, tiny(rhoDisc))

             if (rhoDisc > 1.e-9) then ! large midplane dust grains
                thisOctal%dusttypeFraction(subcell,1:4) = 0.d0
                thisOctal%dusttypeFraction(subcell,1) = 1.d0
             else if (rhoDisc > rhoEnv) then  ! normal disc dust
                thisOctal%dusttypeFraction(subcell,1:4) = 0.d0
                thisOctal%dusttypeFraction(subcell,2) = 1.d0
             endif
          endif


       endif
    enddo

  end subroutine fillDustWhitney

  real function rtnewtdust(x1,x2,xacc, p1, p2) result(junk)

    real :: x1, x2, xacc, p1, p2
    integer :: jmax, j
    real ::  dx, f, df
    parameter (jmax=20)
    df =0.; f = 0.; dx = 0.
    junk = 0.5 * (x1+x2)
    do j=1,jmax
       call equation2dust(junk,f,df,p1,p2)
       dx=f/df
       junk=junk-dx
       if((x1-junk)*(junk-x2).lt.0.) then
          write(*,*) 'RTNEWT: jumped out of brackets',p1,p2,junk
          stop
       endif
       if(abs(dx).lt.xacc) return
    enddo
    write(*,*) 'rtnewt exceeding maximum iterations'
  end function rtnewtdust


  subroutine Equation2dust(mu0, eq2, deq2, r, mu)
    real :: r, mu, mu0
    real :: eq2, deq2

    eq2 = mu0**3 + (r-1.)*mu0 -r*mu
    deq2 = 3.*mu0**2 + r - 1.

  end subroutine Equation2dust


  subroutine parseGrainType(grainString, nTypes, name, abundance)

    character(len=*) :: grainString
    integer :: nTypes
    character(len=*) :: name(:)
    character(len=80) :: tempString
    real :: abundance(:)
    integer :: i, j

    nTypes = 0

    i = index(grainString,":")

    ! not a mixed grain type

    if (i == 0) then
       nTypes = 1
       name(1) = trim(grainString)
       abundance(1) = 1.
       goto 999
    endif
    tempString = grainString
    do while (index(tempString,":") /=0)
       i = index(tempString,":")
       nTypes = nTypes + 1
       name(nTypes) = tempString(1:i-1)
       tempString(1:) = tempString(i+1:)
       if (index(tempString,":") /=0 ) then
          j = index(tempString,":")-1
          read(tempString(1:j), *, err=666) abundance(nTypes)
       else
          j = len(trim(tempString))   
          read(tempString(1:j), *,err=666) abundance(nTypes)
          goto 999
       endif
       tempString(1:) = tempString(j+2:)
    end do

666 continue
    write(*,*) "Error parsing graintype string"
    write(*,*) "Grain string: ",trim(grainString)
    write(*,*) "temp string: ",trim(tempString)
    stop
999 continue
  end subroutine parseGrainType

  subroutine createDustCrossSectionPhaseMatrix(grid, xArray, nLambda, miePhase, nMuMie)

    use mieDistPhaseMatrix_mod
    use input_variables, only : mie, useDust, dustFile, nDustType, graintype, ngrain, &
         grainname, x_grain, amin, amax, a0, qdist, pdist, dustToGas, scale, &
         dustfilename, isotropicScattering, readmiephase, writemiephase, useOldMiePhaseCalc, &
         ttau_disc_on, grainFrac
    real, allocatable :: mReal(:,:), mImg(:,:), tmReal(:), tmImg(:)
    real, allocatable :: mReal2D(:,:), mImg2D(:,:)
    type(PHASEMATRIX),pointer :: miePhase(:,:,:)
    integer :: nMuMie
    real :: mu, total_dust_abundance
    real :: kAbs, kSca
    integer :: i, j, k
    character(len=80) :: message


    type(GRIDTYPE) :: grid
    real :: xArray(:)
    integer :: nLambda

    ! Note: the first index should be either lambda or mu
    !       in order to speedup the array operations!!!  (RK) 
    if (associated(miePhase)) then
       deallocate(miePhase)
       nullify(miePhase)
    endif

    allocate(miePhase(1:nDustType,1:nLambda,1:nMumie)) 



    if (mie .or. useDust) then
       if (associated(grid%onekappaAbs)) deallocate(grid%onekappaAbs)
       if (associated(grid%onekappaSca)) deallocate(grid%onekappaSca)
       allocate(grid%oneKappaAbs(1:nDustType, 1:nLambda), grid%oneKappaSca(1:nDustType, 1:nLambda))

       if (.not.dustfile) then
          write(message,'(a,f5.2)') "Multiplying the opacities by the dust-to-gas ratio of: ",dusttogas
          call writeInfo(message, FORINFO)
          do i = 1, nDustType
             call parseGrainType(graintype(i), ngrain, grainname, x_grain)
             call fillGridMie(grid, scale, aMin(i), aMax(i), a0(i), qDist(i), pDist(i), &
                  ngrain, X_grain, grainname, i)
             grid%oneKappaAbs(i,1:grid%nLambda) =  grid%oneKappaAbs(i,1:grid%nLambda) * dustToGas
             grid%oneKappaSca(i,1:grid%nLambda) =  grid%oneKappaSca(i,1:grid%nLambda) * dustToGas
          enddo
       else
          do i = 1, nDustType
             call dustPropertiesfromFile(dustfilename(i), grid%nlambda, xArray, &
                  grid%onekappaAbs(i,1:grid%nlambda), grid%onekappaSca(i,1:grid%nLambda))
          enddo
       endif
       if (writeoutput) then
          open(20, file="albedo.dat", status="unknown", form="formatted")
          write(20,'(a120)') "# Columns are: wavelength (microns), kappa ext (cm^2 g^-1), &
             &  kappa abs (cm^2 g^-1), kappa sca (cm^2 g^-1), albedo"
          write(20,*) "# Note that the opacities are per gram of dust"
          do i = 1, nLambda
             kAbs = SUM(grid%oneKappaAbs(1:nDustType,i)*grainFrac(1:nDustType))/1.e10/dusttogas
             kSca = SUM(grid%oneKappaSca(1:nDustType,i)*grainFrac(1:nDustType))/1.e10/dusttogas
             write(20,*) xArray(i)*angstomicrons, kAbs+kSca, kAbs, kSca, kSca/(kAbs+kSca)
          enddo
          close(20)
       endif



       call writeInfo("Creating Rosseland opacity lu table",TRIVIAL)
       call createRossArray(grid)
       call writeInfo("Done.",TRIVIAL)
    endif


    if (mie .or. (grid%geometry == "ttauri" .and. ttau_disc_on)) then
       ! construct the mie phase matrix
       call writeInfo("Computing Mie phase grid...",TRIVIAL)

       if (isotropicScattering) then
          miePhase = fillIsotropic()
          call writeInfo("Completed.",TRIVIAL)
          return
       endif

       if (readMiePhase) then
          open(144, file='miephasefile', status="old", form="unformatted")
          read(unit=144) miePhase
          close(144)

       else

          allocate(mReal(1:nDusttype,1:nLambda))
          allocate(mImg(1:nDustType,1:nLambda))
          allocate(tmReal(1:nLambda))
          allocate(tmImg(1:nLambda))

          ! Set up refractive indices
          do i = 1, nDustType
             call parseGrainType(graintype(i), ngrain, grainname, x_grain)
             ! quick test for zero total dust abundance.
             total_dust_abundance = SUM(x_grain)
             if (total_dust_abundance <= 0.0) then
                write(*,*) "Error:: total_dust_abundance <= 0.0 in torusMain."
                write(*,*) "  ==> You probably forgot to assign dust abundance in your parameter file!"
                write(*,*) "  ==> Exiting the program ... "
                stop 
             end if

             ! allocate mem for temp arrays
             allocate(mReal2D(1:ngrain, 1:nLambda))
             allocate(mImg2D(1:ngrain, 1:nLambda))
             ! initializing the values
             mReal2D(:,:) = 0.0; mImg2D(:,:) = 0.0

             ! Find the index of refractions for all types of grains available
             do k = 1, ngrain
                call getRefractiveIndex(xArray, nLambda, grainname(k), tmReal, tmImg)
                mReal2D(k,:) = tmReal(:)  ! copying the values to a 2D maxtrix
                mImg2D(k,:)  = tmImg(:)   ! copying the values to a 2D maxtrix            
             end do

             ! Finding the weighted average of the refractive index.
             mReal(i,:) = 0.0; mImg(i,:) = 0.0
             do j = 1, nLambda
                do k = 1, ngrain
                   mReal(i,j) = mReal(i,j) + mReal2D(k,j)*X_grain(k)
                   mImg(i,j)  = mImg(i,j) + mImg2D(k,j) *X_grain(k)
                end do
                mReal(i,j) = mReal(i,j) / total_dust_abundance
                mImg(i,j)  = mImg(i,j)  / total_dust_abundance
             end do

             deallocate(mReal2D)
             deallocate(mImg2D)
          end do

          deallocate(tmReal)
          deallocate(tmImg)

          ! You should use the new wrapper as it is much faster.
          ! The old and new methods give exactly the same result when optimisations
          ! are turned off.
          ! When optimisations are turned on, the methods give different result,
          ! and *neither* matches the result obtained when optimisations are off.
          if (useOldMiePhaseCalc) then
             do i = 1, nDustType
                do j = 1, nLambda
                   do k = 1, nMumie
                      mu = 2.*real(k-1)/real(nMumie-1)-1.
                      call mieDistPhaseMatrixOld(aMin(i), aMax(i), a0(i), qDist(i), pDist(i), &
                           xArray(j), mu, miePhase(i,j,k), mReal(i,j), mImg(i,j))
                   enddo
                   call normalizeMiePhase(miePhase(i,j,1:nMuMie), nMuMie)
                end do
             end do
          else
             call mieDistPhaseMatrixWrapper(nDustType, nLambda, nMuMie, xArray, mReal, mImg, miePhase)
          end if

          deallocate(mReal)
          deallocate(mImg)
       endif

       if (writeMiePhase) then
          open(144, file='miephasefile', status="replace", form="unformatted")
          write(unit=144) miePhase
          close(144)
       end if


       call fixMiePhase(miePhase, nDustType, nLambda, nMuMie)

       call writeInfo("Completed.",TRIVIAL)
    endif
  end subroutine createDustCrossSectionPhaseMatrix


  recursive subroutine allocateMemoryForDust(thisOctal)
    use input_variables, only : nDustType
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call allocateMemoryForDust(child)
                exit
             end if
          end do
       else
          call allocateAttribute(thisOctal%oldFrac, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%dustType, thisOctal%maxChildren)
          ALLOCATE(thisOctal%dusttypefraction(thisOctal%maxchildren,  nDustType))
          thisOctal%dustTypeFraction(thisOctal%maxchildren,1:nDustType) = 0.d0
          thisOctal%dustTypeFraction(thisOctal%maxchildren,1) = 1.d0
       endif
    enddo
  end subroutine allocateMemoryForDust
  

end module dust_mod
