!
! a module to do 3d statistical equilibrium according to Klein & Castor (1978)
! - lots of it based on my "stateq" PhD code
!
module stateq_mod

  use gridtype_mod
  use grid_mod
  use math_mod
  use vector_mod
  use constants_mod
  use path_integral
  use jets_mod
  use opacity_lte_mod
  use hyd_col_coeff
				   
  implicit none
!  public

!  public initGridStateq, expint, gii, giia, beta_mn, cijt
 
  integer, private :: level

  real(double) :: bEinstein(20, 20) 
  
  real(double), parameter :: eTrans(23) =                      &
       (/ (hydE0eVdb - hydE0eVdb/level**2,level=1,SIZE(eTrans)) /)
  
  real(double), parameter :: gDegen(23) = &
       (/ (2.0_db*level**2,level=1,SIZE(gDegen)) /)
       
  ! data from: 
  !   W. L. Wiese, M. W. Smith, and B. M. Glennon
  !   Atomic transition probabilities : a critical data compilation
  !   National standard reference data series NSRDS-NBS / 
  !     United States. National Bureau of Standards ; 4, 22
  
  real(double), parameter :: aEinstein(20,20) = reshape( source=&                 
    (/ 0.000d0, 4.699d8, 5.575d7, 1.278d7, 4.125d6, 1.644d6, 7.568d5, 3.869d5, 2.143d5, 1.263d5, & 
       7.834d4, 5.066d4, 3.393d4, 2.341d4, 1.657d4, 1.200d4, 8858.d0, 6654.d0, 5077.d0, 3928.d0, & 
       0.000d0, 0.000d0, 4.410d7, 8.419d6, 2.530d6, 9.732d5, 4.389d5, 2.215d5, 1.216d5, 7.122d4, &
       4.397d4, 2.834d4, 1.893d4, 1.303d4, 9210.d0, 6658.d0, 4910.d0, 3685.d0, 2809.d0, 2172.d0, &
       0.000d0, 0.000d0, 0.000d0, 8.986d6, 2.201d6, 7.783d5, 3.358d5, 1.651d5, 8.905d4, 5.156d4, &
       3.156d4, 2.021d4, 1.343d4, 9211.d0, 6490.d0, 4680.d0, 3444.d0, 2580.d0, 1964.d0, 1517.d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 2.699d6, 7.711d5, 3.041d5, 1.424d5, 7.459d4, 4.235d4, &
       2.556d4, 1.620d4, 1.069d4, 7288.d0, 5110.d0, 3671.d0, 2693.d0, 2013.d0, 1529.d0, 1178.d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 1.025d6, 3.253d5, 1.388d5, 6.908d4, 3.800d4, &
       2.246d4, 1.402d4, 9148.d0, 6185.d0, 4308.d0, 3079.d0, 2249.d0, 1675.d0, 1268.d0, 975.1d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 4.561d5, 1.561d5, 7.065d4, 3.688d4, &
       2.110d4, 1.288d4, 8271.d0, 5526.d0, 3815.d0, 2707.d0, 1966.d0, 1457.d0, 1099.d0, 842.4d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 2.272d5, 8.237d4, 3.905d4, &
       2.117d4, 1.250d4, 7845.d0, 5156.d0, 3516.d0, 2471.d0, 1781.d0, 1312.d0, 984.9d0, 751.7d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 1.233d5, 4.676d4, &
       2.301d4, 1.287d4, 7804.d0, 5010.d0, 3359.d0, 2331.d0, 1664.d0, 1216.d0, 906.9d0, 688.6d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 7.151d4, &
       2.813d4, 1.427d4, 8192.d0, 5080.d0, 3325.d0, 2268.d0, 1598.d0, 1156.d0, 855.5d0, 645.2d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, &
       4.377d4, 1.774d4, 9231.d4, 5417.d0, 3424.d0, 2280.d0, 1578.d0, 1127.d0, 825.2d0, 617.3d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, &
       0.000d0, 2.799d4, 1.163d0, 6186.d0, 3699.d0, 2377.d0, 1606.d0, 1127.d0, 814.1d0, 602.6d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, &
       0.000d0, 0.000d0, 1.857d4, 7884.d0, 4271.d0, 2596.d0, 1693.d0, 1159.d0, 822.3d0, 600.5d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, &
       0.000d0, 0.000d0, 0.000d0, 1.271d4, 5496.d0, 3026.d0, 1866.d0, 1232.d0, 853.2d0, 611.9d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 8933.d0, 3926.d0, 2192.d0, 1369.d0, 914.4d0, 639.7d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 6429.d0, 2864.d0, 1620.d0, 1023.d0, 690.3d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 4720.d0, 2130.d0, 1217.d0, 776.7d0, & 
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 3530.d0, 1610.d0, 929.6d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 2680.d0, 1235.d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, &
       0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 2067.d0, & 
       0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /), shape=(/20,20/))
  
  real(double) :: fStrength(20,20) = reshape( source=& 
    (/ 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,                    &
       4.162d-1, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       7.910d-2, 6.407d-1, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       2.899d-2, 1.193d-1, 8.421d-1, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       1.394d-2, 4.467d-2, 1.506d-1, 1.038d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       7.799d-3, 2.209d-2, 5.584d-2, 1.793d-1, 1.231d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       4.814d-3, 1.270d-2, 2.768d-2, 6.549d-2, 2.069d-1, 1.424d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, & 
       0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       3.183d-3, 8.036d-3, 1.604d-2, 3.230d-2, 7.448d-2, 2.340d-1, 1.616d00, 0.000d00, 0.000d00, 0.000d00, &
       0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       2.216d-3, 5.429d-3, 1.023d-2, 1.870d-2, 3.645d-2, 8.315d-2, 2.609d-1, 1.807d00, 0.000d00, 0.000d00, &
       0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       1.605d-3, 3.851d-3, 6.980d-3, 1.196d-2, 2.104d-2, 4.038d-2, 9.163d-2, 2.876d-1, 1.999d00, 0.000d00, &
       0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       1.201d-3, 2.835d-3, 4.996d-3, 8.187d-3, 1.344d-2, 2.320d-2, 4.416d-2, 1.000d-1, 3.143d-1, 2.190d00, &
       0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       9.214d-4, 2.151d-3, 3.711d-3, 5.886d-3, 9.209d-3, 1.479d-2, 2.525d-2, 4.787d-2, 1.083d-1, 3.408d-1, &
       2.381d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       7.227d-4, 1.672d-3, 2.839d-3, 4.393d-3, 6.631d-3, 1.012d-2, 1.605d-2, 2.724d-2, 5.152d-2, 1.166d-1, &
       3.673d-1, 2.572d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       5.774d-4, 1.326d-3, 2.224d-3, 3.375d-3, 4.959d-3, 7.289d-3, 1.097d-2, 1.726d-2, 2.918d-2, 5.513d-2, &
       1.248d-1, 3.938d-1, 2.763d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       4.686d-4, 1.070d-3, 1.776d-3, 2.656d-3, 3.821d-3, 5.455d-3, 7.891d-3, 1.177d-2, 1.843d-2, 3.109d-2, &
       5.872d-2, 1.330d-1, 4.202d-1, 2.954d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       3.856d-4, 8.764d-4, 1.443d-3, 2.131d-3, 3.014d-3, 4.207d-3, 5.905d-3, 8.456d-3, 1.254d-2, 1.958d-2, &
       3.298d-2, 6.228d-2, 1.412d-1, 4.467d-1, 3.145d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, &
       3.211d-4, 7.270d-4, 1.188d-3, 1.739d-3, 2.425d-3, 3.324d-3, 4.556d-3, 6.323d-3, 8.995d-3, 1.328d-2, &
       2.070d-2, 3.486d-2, 6.584d-2, 1.494d-1, 4.731d-1, 3.336d00, 0.000d00, 0.000d00, 0.000d00, 0.000d00, & 
       2.702d-4, 6.099d-4, 9.916d-4, 1.439d-3, 1.984d-3, 2.679d-3, 3.602d-3, 4.877d-3, 6.719d-3, 9.515d-3, &
       1.402d-2, 2.182d-2, 3.672d-2, 6.938d-2, 1.575d-1, 4.995d-1, 3.527d-1, 0.000d00, 0.000d00, 0.000d00, &
       2.296d-4, 5.167d-4, 8.361d-4, 1.204d-3, 1.646d-3, 2.196d-3, 2.905d-3, 3.856d-3, 5.180d-3, 7.099d-3, &
       1.002d-2, 1.474d-2, 2.292d-2, 3.858d-2, 7.292d-2, 1.657d-1, 5.259d-1, 3.718d00, 0.000d00, 0.000d00, &
       1.967d-4, 4.416d-4, 7.118d-4, 1.019d-3, 1.382d-3, 1.825d-3, 2.383d-3, 3.112d-3, 4.094d-3, 5.468d-3, &
       7.468d-3, 1.052d-2, 1.545d-2, 2.402d-2, 4.043d-2, 7.644d-2, 1.738d-1, 5.523d-1, 3.909d00, 0.000d00 /), shape=(/20,20/))
 
  real(double) :: lambdaTrans(20,20) = reshape( source=&
    (/ 000000.d-8,1215.67d-8,1025.72d-8,992.537d-8,949.743d-8,937.803d-8,930.748d-8,926.226d-8,923.150d-8,920.963d00,&
       919.352d-8,918.129d-8,917.181d-8,916.429d-8,915.824d-8,915.329d-8,914.919d-8,914.576d-8,914.286d-8,914.039d-8,&  
       0.00000d-8,0000000d-8,6562.80d-8,4861.32d-8,4340.36d-8,4101.73d-8,3970.07d-8,3889.05d-8,3835.38d-8,3797.90d-8,&
       3770.63d-8,3750.15d-8,3734.37d-8,3721.94d-8,3711.97d-8,3703.85d-8,3697.15d-8,3691.55d-8,3686.83d-8,3682.81d-8,&  
       0.00000d-8,0.00000d-8,0000000d-8,18751.0d-8,12818.1d-8,10938.1d-8,10049.4d-8,9545.98d-8,9229.02d-8,9014.91d-8,&
       8862.79d-8,8750.47d-8,8665.02d-8,8598.39d-8,8545.39d-8,8502.49d-8,8467.26d-8,8437.96d-8,8413.32d-8,8392.40d-8,&  
       0.00000d-8,0.00000d-8,0.00000d-8,0000000d-8,40512.0d-8,26252.0d-8,21655.0d-8,19445.6d-8,18174.1d-8,17362.1d-8,&
       16806.5d-8,16407.2d-8,16109.3d-8,15880.5d-8,15700.7d-8,15556.5d-8,15438.9d-8,15341.8d-8,15260.6d-8,15191.8d-8,&  
       0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0000000d-8,74578.0d-8,46525.0d-8,37395.0d-8,32961.0d-8,30384.0d-8,&
       28722.0d-8,27575.0d-8,26744.0d-8,26119.0d-8,25636.0d-8,25254.0d-8,24946.0d-8,24693.0d-8,24483.0d-8,24307.0d-8,&  
       0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,123680.d-8,75005.0d-8,59066.0d-8,51273.0d-8,&
       46712.0d-8,43753.0d-8,41697.0d-8,40198.0d-8,39065.0d-8,38184.0d-8,37484.0d-8,36916.0d-8,36449.0d-8,36060.0d-8,&  
       0.00000d-8,0000000d-8,0000000d-8,0000000d-8,0000000d-8,0.00000d-8,0000000d-8,190570.d-8,113060.d-8,87577.0d-8,&
       75061.0d-8,67701.0d-8,62902.0d-8,59552.0d-8,57099.0d-8,55237.0d-8,53783.0d-8,52622.5d-8,51679.0d-8,50899.0d-8,&  
       0.00000d-8,0000000d-8,0000000d-8,0000000d-8,0000000d-8,0.00000d-8,0000000d-8,0000000d-8,277960.d-8,162050.d-8,&
       123840.d-8,105010.d-8,93894.0d-8,86621.0d-8,81527.0d-8,77782.0d-8,74930.0d-8,72696.0d-8,70908.0d-8,69448.0d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,388590.d-8,&  
       223340.d-8,168760.d-8,141790.d-8,125840.d-8,115360.d-8,108010.d-8,102580.d-8,98443.0d-8,95191.0d-8,92579.0d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       525200.d-8,298310.d-8,223250.d-8,186100.d-8,164070.d-8,149580.d-8,139380.d-8,131840.d-8,126080.d-8,121530.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,690500.d-8,388320.d-8,288230.d-8,238620.d-8,209150.d-8,189730.d-8,176030.d-8,165900.d-8,158120.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,887300.d-8,494740.d-8,364610.d-8,300020.d-8,261610.d-8,236260.d-8,218360.d-8,205090.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,111800.d-7,619000.d-8,453290.d-8,371000.d-8,322000.d-8,289640.d-8,266740.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,138600.d-7,762300.d-8,555200.d-8,452220.d-8,390880.d-8,350300.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,169400.d-7,926100.d-8,671200.d-8,544400.d-8,468760.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,204400.d-7,111200.d-7,802300.d-8,648200.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,243800.d-7,132100.d-7,949200.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,288200.d-7,155400.d-7,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,337400.d-7,&  
       0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /), shape=(/20,20/))

contains

  real function beta_mn(m, n, rVec, i1, i2, i3, grid, thisOctal, thisSubcell)

    type(GRIDTYPE), intent(in)   :: grid
    type(VECTOR), intent(in)     :: rVec
    integer, intent(in)          :: i1, i2, i3
    integer, intent(in)          :: m,n
    type(OCTAL), pointer,optional:: thisOctal
    integer, intent(in),optional :: thisSubcell 
    
    type(VECTOR)      :: direction
    integer           :: i, j
    real              :: theta, phi
    integer,parameter :: nTheta = 10
    integer           :: nPhi
    real              :: dTheta, dPhi, dOmega
    real              :: totomega
    real(double)              :: escprob,  tau_mn
    type(octalVector) :: rVecOctal
    type(OCTAL),pointer:: octalCopy
    
    dtheta = pi / real(ntheta-1)
    beta_mn = 0.
    totomega = 0.
    
    do i = 1, ntheta
       theta = pi*real(i-1)/real(ntheta-1)
       nphi = max(2,nint(real(ntheta)*sin(theta)))
       dphi = twopi / real(nphi-1)
       do j = 1, nphi-1
          phi = twopi * real(j-1)/real(nphi-1)
          direction = vector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta))
          domega = sin(theta)*dtheta*dphi
          totomega = totomega + domega

!          call integratePath(real(lambdaTrans(m,n))*1.e8, real(lambdaTrans(m,n))*1.e8, &
!                   grid%velocity(i1,i2,i3), &
!                   rVec, direction, grid, &
!                   lambdaArray, tauExt, tauAbs, tauSca, maxTau, nTau, .false., &
!                   escProb, .false., real(lambdaTrans(m,n))*1.e8, real(lambdaTrans(m,n))*1.e8, &
!                   1, contTau, hitCore, .false.,.false., 0.,&
!                   VECTOR(0.,0.,0.), 1., .true., m, n, real(fstrength(m,n)), real(gDegen(m)), real(gDegen(n)), &
!                    tau_mn)

          tau_mn = (pi*echarge**2)/(melectron*cspeed)
          tau_mn = tau_mn * gDegen(m) * fStrength(m,n)
          if (grid%adaptive) then
             if (grid%geometry(1:4) == "jets")  then
                tau_mn = tau_mn * abs((thisOctal%N(thisSubcell,m)/gDegen(m)) - &
                     (thisOctal%N(thisSubcell,n)/gDegen(n)))
                tau_mn = tau_mn/dV_dn_jets(rVec, direction)
             else
                tau_mn = tau_mn * abs((thisOctal%N(thisSubcell,m)/gDegen(m)) - &
                                      (thisOctal%N(thisSubcell,n)/gDegen(n))) ! eq 5.
                rVecOctal = rVec
                octalCopy => thisOctal
                tau_mn = tau_mn / (amrGridDirectionalDeriv(grid,rVecOctal,direction, &
                                                        startOctal=octalCopy) / 1.e10)
             end if
          else
             tau_mn = tau_mn * abs((grid%n(i1,i2,i3,m)/gDegen(m)) - &
                                   (grid%n(i1,i2,i3,n)/gDegen(n))) ! eq 5.
             tau_mn = tau_mn / (directionalderiv(grid,s2o(rvec),i1,i2,i3,direction)/1.e10)
          end if
          
          tau_mn = tau_mn *  lambdaTrans(m,n) / cSpeed_dbl

          if (tau_mn < 1.e-20) then
             tau_mn = 1.e-20 
          endif

          if (tau_mn < 0.1) then
             escProb = 1.0-tau_mn*0.5*(1.0 - tau_mn/3.0*(1. - tau_mn*0.25*(1.0 - 0.20*tau_mn)))
          else if (tau_mn < 15.) then
             escProb = (1.0-exp(-tau_mn))/tau_mn
          else
             escProb = 1./tau_mn
          end if
          beta_mn = beta_mn + escprob * domega

       enddo
    enddo
    beta_mn = beta_mn / totOmega 

!   if ((m == 1).and.(n == 2)) beta_mn = beta_mn * 100.
!    write(*,*) totOmega/fourPi
!   write(*,*) m,n,beta_mn

  end function beta_mn


  real function beta_cmn(m,n,rVec,i1,i2,i3,grid,nstar,thisOctal,thisSubcell)

    integer, intent(in)         :: m,n
    type(gridtype), intent(in)  :: grid
    type(vector), intent(in)    :: rVec
    integer, intent(in)         :: i1, i2, i3
    integer, parameter          :: ntheta=10
    integer, parameter          :: nphi=10
    integer, intent(in)         :: nstar
    type(OCTAL),pointer,optional:: thisOctal
    integer,intent(in),optional :: thisSubcell 
    
    real :: r, thetatostar, phiTostar
    real :: h, totomega, disttooccult
    integer :: i, j
    type(vector) :: tostar, occultposition, starposition
    type(vector) :: tooccult, direction
    real :: starradius, occultradius
    real :: disttostar
    real :: cosang, sinang, ang
    logical :: occulted, full
    real :: tau_cmn
    real :: dtheta, dphi
    real :: theta, phi, domega
    real :: escprob, dotprod
    integer, parameter :: maxTau = 4000
    integer :: nTau
    real, allocatable :: contTau(:,:)
    real :: tauExt(maxTau)
    real :: tauAbs(maxTau)
    real :: tauSca(maxTau)
    real :: theta0, phi0
    real :: lambdaArray(maxTau)
    !logical :: hitCore
    real :: tauLocal
    type(octalVector) :: rVecOctal
    type(OCTAL),pointer :: octalCopy

    full = .false.

    if (nstar == 1) then
       starradius = grid%rstar1
       occultradius = grid%rstar2
       starposition = grid%starpos1
       occultposition = grid%starpos2
    else
       starradius = grid%rstar2
       occultradius = grid%rstar1
       starposition = grid%starpos2
       occultposition = grid%starpos1
    endif

    tostar = starposition - rVec
    disttostar = modulus(tostar)
    tostar = tostar / disttostar

    tooccult = occultposition - rVec
    disttooccult = modulus(tooccult)
    tooccult = tooccult / disttooccult

    sinang = starradius / disttostar
    ang = asin(min(1.,max(-1.,sinAng)))
    call getPolar(toStar, r, thetaTostar, phiToStar)

    dtheta = pi / real(ntheta-1)
    dphi = twopi / real(nphi-1)
    beta_cmn = 0.
    totomega = 0.
    do i = 1, ntheta
       theta0 = (2.*real(i-1)/real(nTheta-1)-1.)*ang
       theta = thetaToStar + theta0
       dTheta = 2.*ang/real(nTheta-1)
       do j = 1, nphi
          phi0 = (2.*real(j-1)/real(nPhi-1)-1.)*ang
          phi = phiToStar + phi0
!          if (phi < 0.) phi = phi + twoPi
!          if (phi > twoPi) phi = phi - twoPi
!       if (theta > pi) theta = theta - pi
!       if (theta < 0.) theta = theta + pi

          dphi = 2.*ang/real(nPhi-1)
          occulted = .false.
          direction = vector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta))
          if (theta > pi) theta = theta - pi
          if (theta < 0.) theta = theta + pi
          dOmega = dtheta*dphi

          dotprod = direction .dot. tostar
          escProb = 0.
          if ((dotprod > 0.)  .and. (acos(min(1.,max(-1.,dotprod))) < 1.01*ang)) then

             dotprod = direction .dot. tooccult
             if ((dotprod > 0.).and.(distToOccult < distToStar)) then
                sinang = sqrt(max(0.,(1.-dotprod*dotprod)))
                if ((sinang*disttooccult) < occultradius) then 
                   occulted = .true.
                endif
             endif
!    write(*,*) "calc",pi*starradius**2/disttostar**2

             if (.not.occulted) then

                if (full) then
!                   call integratePath(real(lambdaTrans(m,n))*1.e8, real(lambdaTrans(m,n))*1.e8, &
!                        grid%velocity(i1,i2,i3), &
!                        rVec, direction, grid, &
!                        lambdaArray, tauExt, tauAbs, tauSca, maxTau, nTau, .true., &
!                        escProb, .false., real(lambdaTrans(m,n))*1.e8, real(lambdaTrans(m,n))*1.e8, &
!                        1, contTau, hitCore, .false.,.false., 0.,&
!                        VECTOR(0.,0.,0.), 1., .true., m, n, real(fstrength(m,n)), real(gDegen(m)), real(gDegen(n)), &
!                        tauLocal)
!                   tau_cmn = tauLocal + tauExt(nTau-1)
                else


                   tau_cmn = (pi*echarge**2)/(melectron*cspeed)
                   tau_cmn = tau_cmn * gDegen(m) * fStrength(m,n)
                   if (grid%adaptive) then
                      tau_cmn = tau_cmn * abs((thisOctal%N(thisSubcell,m)/gDegen(m)) - &
                                              (thisOctal%N(thisSubcell,n)/gDegen(n))) ! eq 5.
                   else 
                      tau_cmn = tau_cmn * abs((grid%n(i1,i2,i3,m)/gDegen(m)) - &
                                              (grid%n(i1,i2,i3,n)/gDegen(n))) ! eq 5.
                   endif
                   if (tau_cmn < 0.) then
                      tau_cmn = 0. !abs(tau_cmn)
                   endif
                   tau_cmn = tau_cmn  / (cSpeed / lambdaTrans(m,n))
                   if (grid%adaptive) then
                      if (grid%geometry(1:4) == "jets")  then
                         tau_cmn = tau_cmn/dV_dn_jets(rVec, direction)
                      else
                         rVecOctal = rVec
                         
                         octalCopy => thisOctal
                         tau_cmn = tau_cmn / (amrGridDirectionalDeriv(grid,rVecOctal,direction, &
                                                                   startOctal=octalCopy) / 1.e10)
                      end if
                   else
                      tau_cmn = tau_cmn / (directionalderiv(grid,s2o(rVec),i1,i2,i3,direction)/1.e10)
                   end if
                endif

                if (tau_cmn < 0.1) then
                   escProb = 1.d0-tau_cmn*0.5d0*(  1.0d0 - tau_cmn/3.0d0*( 1.0d0 - &
                        tau_cmn*0.25d0*(1.0d0 -0.20d0*tau_cmn) )  )
                else
                   escprob = (1.d0 - exp(-tau_cmn))/tau_cmn
                endif
                beta_cmn = beta_cmn + escprob * domega
                totOmega = totOmega + dOmega
             endif
          endif
       enddo
    enddo
    beta_cmn =   beta_cmn / fourpi
!    write(*,*) "calc",pi*starradius**2/disttostar**2
!    write(*,*) "found",totOmega
  end function beta_cmn

  subroutine beta_cmn_sub(m,n,i1,i2,i3,grid,nstar)

    integer :: m,n
    type(gridtype) :: grid
    type(vector) :: rvec
    integer :: i1, i2, i3
    integer :: i, j, nphi=5, ntheta=5
    real :: r, thetatostar, phiTostar,betacmn
    real :: h, totomega, disttooccult
    integer nstar
    type(vector) :: tostar, occultposition, starposition
    type(vector) :: tooccult, direction
    real :: starradius, occultradius
    real :: disttostar
    real :: cosang, sinang, ang
    logical :: occulted, full
    real :: tau_cmn
    real :: dtheta, dphi
    real :: theta, phi, domega
    real :: escprob, dotprod
    integer, parameter :: maxTau = 4000
    integer :: nTau
    real, allocatable :: contTau(:,:)
    real :: tauExt(maxTau)
    real :: tauAbs(maxTau)
    real :: tauSca(maxTau)
    real :: lambdaArray(maxTau)
    logical :: hitCore
    real :: tauLocal

    full = .false.

    rVec = VECTOR(grid%xAxis(i1),grid%yAxis(i2),grid%zAxis(i3))

    if (nstar == 1) then
       starradius = grid%rstar1
       occultradius = grid%rstar2
       starposition = grid%starpos1
       occultposition = grid%starpos2
    else
       starradius = grid%rstar2
       occultradius = grid%rstar1
       starposition = grid%starpos2
       occultposition = grid%starpos1
    endif

    tostar = starposition - rvec
    disttostar = modulus(tostar)
    tostar = tostar / disttostar

    tooccult = occultposition - rvec
    disttooccult = modulus(tooccult)
    tooccult = tooccult / disttooccult

    h  = sqrt(max(0.,(disttostar**2 - starradius**2)))
    cosang = h / disttostar
    ang = acos(min(1.,max(-1.,cosAng)))

    call getPolar(toStar, r, thetaTostar, phiToStar)

    dtheta = pi / real(ntheta-1)
    dphi = twopi / real(nphi-1)
    betacmn = 0.
    totomega = 0.
    do i = 1, ntheta
       theta = thetaToStar + (2.*real(i-1)/real(nTheta-1)-1.)*ang
       if (theta > pi) theta = theta - pi
       if (theta < 0.) theta = theta + pi
       dTheta = 2.*ang/real(nTheta-1)
       do j = 1, nphi
          phi = phiToStar + (2.*real(j-1)/real(nPhi-1)-1.)*ang
          if (phi < 0.) phi = phi + twoPi
          if (phi > twoPi) phi = phi + twoPi
          dphi = 2.*ang/real(nPhi-1)
          occulted = .false.
          direction = vector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta))
          domega = sin(theta)*dtheta*dphi
          dotprod = direction .dot. tostar
          escProb = 0.
          if ((dotprod > 0.) .and. (acos(min(1.,max(-1.,dotprod))) < ang)) then

             dotprod = direction .dot. tooccult
             if ((dotprod > 0.).and.(distToOccult < distToStar)) then
                sinang = sqrt(max(0.,(1.-dotprod*dotprod)))
                if ((sinang*disttooccult) < occultradius) then 
                   occulted = .true.
                endif
             endif


             if (.not.occulted) then

                if (full) then
!                   call integratePath(real(lambdaTrans(m,n))*1.e8, real(lambdaTrans(m,n))*1.e8, &
!                        grid%velocity(i1,i2,i3), &
!                        rVec, direction, grid, &
!                        lambdaArray, tauExt, tauAbs, tauSca, maxTau, nTau, .true., &
!                        escProb, .false., real(lambdaTrans(m,n))*1.e8, real(lambdaTrans(m,n))*1.e8, &
!                        1, contTau, hitCore, .false.,.false., 0.,&
!                        VECTOR(0.,0.,0.), 1., .true., m, n, real(fstrength(m,n)), real(gDegen(m)), real(gDegen(n)), &
!                        tauLocal)
!                   tau_cmn = tauLocal + tauExt(nTau-1)
                else


                   tau_cmn = (pi*echarge**2)/(melectron*cspeed)
                   tau_cmn = tau_cmn * gDegen(m) * fStrength(m,n)
                   tau_cmn = tau_cmn * abs((grid%n(i1,i2,i3,m)/gDegen(m)) - (grid%n(i1,i2,i3,n)/gDegen(n))) ! eq 5.
                   if (tau_cmn < 0.) then
                      tau_cmn = 0. !abs(tau_cmn)
                   endif
                   tau_cmn = tau_cmn  / (cSpeed / lambdaTrans(m,n))
                   tau_cmn = tau_cmn / (directionalderiv(grid,s2o(rvec),i1,i2,i3,direction)/1.e10)
                endif

                if (tau_cmn < 0.1) then
                   escProb = 1.-tau_cmn*0.5*(  1.0 - tau_cmn/3.0*( 1.0 - &
                        tau_cmn*0.25*(1.0 -0.20*tau_cmn) )  )
                else
                   escprob = (1. - exp(-tau_cmn))/tau_cmn
                endif
                betacmn = betacmn + escprob * domega
                totOmega = totOmega + dOmega
             endif
          endif
       enddo
    enddo
    betacmn =   betacmn / fourpi
!    write(*,*) "calc",fourPi*(pi*starradius**2)/(fourPi * disttostar**2)
!    write(*,*) "found",totOmega
  end subroutine beta_cmn_sub



  ! this subroutine fills the grid with nlte level populations
  ! the grid should have temperatures/densities at all points

  subroutine initgridstateq(grid, contfile1, contfile2, popFileName, &
       readPops, writePops, lte, nLower, nUpper)
    implicit none
    type(gridtype) :: grid
    integer, parameter :: maxLevels = statEqMaxLevels
    integer :: nLower, nUpper
    logical :: lte, debugInfo
    real :: sinTheta
    logical :: readPops, writePops
    character(len=*) :: popFileName
    real :: visFrac1, visFrac2
    logical :: isBinary
    integer :: i, j, k, m
    real(double) :: freq
    real :: hnu1(2000), hnu2(2000)
    real :: nuarray1(2000), nuarray2(2000)
    real :: percentDone
    integer :: nnu1, nnu2
    character(len=*) :: contfile1, contfile2
    real(double), parameter :: tolx = 1.d-3
    real(double), parameter :: tolf = 1.d-3
    real(double) :: x(maxLevels+1), nTot
    real :: t1, t2, t3
    integer :: i1, i2, i3
    type(vector) :: rvec, rHat, thisVec
    real :: departCoeff(maxLevels), ang, r
    real, allocatable :: departCoeffall(:)
    logical :: oneD, twoD, threeD

    integer :: nIter, iIter
    real :: oldLevels(maxLevels+1)
    real(double), allocatable :: xall(:)
    real (kind=double):: ne1, ne2
    real :: temp, crate
    real(double) :: phiT
    logical :: ok
    logical :: aboutX, aboutZ
    integer :: ierr


    debugInfo = .true.
    oneD = .false.
    twoD = .false.
    threeD = .true.
    aboutX = .false.
    aboutZ = .true.


    isBinary = .false.
    if (grid%geometry == "binary") then
       isBinary = .true.
       threeD = .false.
       twoD = .true.
       aboutX = .true.
       aboutZ = .false.
    endif

    if (grid%geometry == "wind") then
       threeD = .false.
       oneD = .true.
       aboutZ = .false.
    endif
    
    if (grid%geometry == "donati") then
       twoD = .true.
       threeD = .false.
       aboutZ = .true.
       aboutX = .false.
    endif

    if (grid%geometry == "ttauri") then
       twoD = .true.
       threeD = .false.
       aboutZ = .true.
       aboutX = .false.
    endif

    if (grid%geometry == "luc_cir3d") then
       twoD = .true.
       threeD = .false.
       aboutZ = .true.
       aboutX = .false.
    endif

    if (grid%geometry == "puls") then
       twoD = .true.
       threeD = .false.
       aboutZ = .true.
       aboutX = .false.
    endif

    if (lte) then
       nIter = 1
    else
       nIter = 1
    endif

    do k = 2, maxlevels
       do i=1, k-1
          lambdaTrans(i, k) = lambdaTrans(k, i)
          ! aEinstein(k, i) = aEinstein(k, i) commented out by nhs
          freq = cspeed / lambdaTrans(i, k)
          bEinstein(k, i) = ((dble(aEinstein(k,i))*dble(cspeed)**2) / (2.d0*dble(hcgs)*dble(freq)**3))
          bEinstein(i, k) = bEinstein(k,i) * gDegen(k) / gDegen(i)
          fStrength(k, i) = -gDegen(i) * fStrength(i,k) / gDegen(k)
       enddo
    enddo





    allocate(grid%n(1:grid%na1, 1:grid%na2, 1:grid%na3,1:maxlevels))
    allocate(grid%ne(1:grid%na1, 1:grid%na2, 1:grid%na3))
    allocate(grid%nTot(1:grid%na1, 1:grid%na2, 1:grid%na3))

    grid%n = 1.e-20
    grid%ne = 1.

    if (readPops) then
       call readGridPopulations(popFilename, grid, maxLevels)
    else
       do i = 1, grid%na1
          do j = 1, grid%na2 
             do k = 1, grid%na3
                if ((.not.grid%instar(i,j,k)).and.grid%inUse(i,j,k)) then
                   nTot = dble(grid%rho(i,j,k) / mhydrogen)


                   if (grid%temperature(i,j,k) < 1.e5) then
                      phiT = 1.5d0*log10(2.d0*pi*mElectron*kErg*grid%temperature(i,j,k)) - &
                           3.d0*log10(hCgs) + log10(exp(-1.0*hydE0eVdb/(kEV*grid%temperature(i,j,k))))
                     
                      phiT = 10.d0**phiT
                      call solveQuadDble(1.d0, phiT, &
                           -1.d0*phiT*dble(nTot), ne1, ne2, ok)
                      grid%ne(i,j,k) = min(max(ne1,ne2),nTot)   
                      grid%ne(i,j,k) = max(grid%ne(i,j,k),1.d0)
                   else
                      grid%ne(i,j,k) = dble(ntot)
                   endif

                   grid%nTot(i,j,k) = nTot 
                   do m = 1, maxlevels
                      grid%n(i,j,k,m) = boltzsaha(m, grid%ne(i,j,k), &
                           dble(grid%temperature(i,j,k)))
                   enddo
                endif
             enddo
          enddo
       enddo




       open(20,file=contfile1,status="old",form="formatted")
       nnu1 = 1
10     continue
       read(20,*,end=20) nuarray1(nnu1), hnu1(nnu1)
       nnu1 = nnu1 + 1
       goto 10
20     continue
       nnu1 = nnu1  - 1
9      close(20)
       hnu1(1:nnu1) = hnu1(1:nnu1) / fourPi   ! Converts from flux to flux momemnt

       if (grid%geometry == "binary") then
          open(20,file=contfile2,status="old",form="formatted")
          nnu2 = 1
30        continue
          read(20,*,end=40) nuarray2(nnu2), hnu2(nnu2)
          nnu2 = nnu2 + 1
          goto 30
40        continue
          nnu2 = nnu2  - 1
          close(20)
          hnu2(1:nnu2) = hnu2(1:nnu2) / fourPi   ! Converts from flux to flux momemnt
       endif

       if (threeD) then

          percentDone = 0.
          departCoeff = 1.

          !$OMP PARALLEL DO &
          !$OMP DEFAULT(NONE) &
          !$OMP PRIVATE(i1, i2, i3) &
          !$OMP SHARED(grid) &
          !$OMP SHARED(lte, hnu1, hnu2, nuArray1, nnu1, nuarray2, nnu2) &
          !$OMP SHARED(isBinary, debugInfo) &
          !$OMP PRIVATE(i, visFrac1, visFrac2, rVec, nTot, percentDone) &
          !$OMP PRIVATE(departCoeffAll, xAll) 
          do i1 = 1, grid%na1
             allocate(departCoeffAll(1:maxlevels))
             departCoeffAll = 1.d0
             do i2 = 1, grid%na2
                do i3 = 1, grid%na3
                   allocate(xall(1:maxLevels+1))
                   if (grid%inUse(i1,i2,i3)) then
                      rvec = vector(grid%xaxis(i1), grid%yaxis(i2), grid%zaxis(i3))
                      if (.not.grid%inStar(i1,i2,i3)) then
                         do i = 1, maxlevels
                            xall(i) = grid%n(i1,i2,i3,i) * dble(departCoeffall(i))
                         enddo
                         xall(maxlevels+1) = grid%ne(i1,i2,i3)
                         if (.not.lte) then

                            if (grid%geometry == "binary") then
                               call occultTest(grid, i1, i2, i3, grid%starPos1, grid%rStar1, &
                                    grid%starPos2, grid%rStar2, visFrac1)
                               call occultTest(grid, i1, i2, i3, grid%starPos2, grid%rStar2, &
                                    grid%starPos1, grid%rStar1, visFrac2)
                            else
                               visFrac1 = 1.
                               visFrac2 = 0.
                            endif


                            call mnewt_stateq(grid, maxLevels+1,xall,maxlevels+1,tolx,tolf, hnu1, &
                                 nuarray1, nnu1, hnu2, nuarray2, nnu2, rvec, i1, i2, i3, visFrac1, visFrac2, isBinary)
                            if ((.not.lte).and.debugInfo) write(*,'(a,3i3,f8.1)') "grid cell",i1,i2,i3,grid%temperature(i1,i2,i3)
                            write(*,*) visFrac1, visFrac2
                            grid%ne(i1,i2,i3) = xall(maxlevels+1)
                            do i = 1 , maxLevels
                               departCoeffall(i) = real(xall(i))/boltzSaha(i, grid%Ne(i1,i2,i3), dble(grid%temperature(i1,i2,i3)))
                               grid%n(i1,i2,i3,i) = xall(i)
                               if ((.not.lte).and.debugInfo) then
                                  write(*,'(i3,1p,e12.3,e12.3)') &
                                       i,departCoeffall(i), xall(i)
                               endif

                            enddo
                            if ((.not.lte).and.debugInfo) write(*,'(a,f6.2)') "log Ne: ",log10(grid%ne(i1,i2,i3))
                            nTot = dble(grid%rho(i1,i2,i3) / mhydrogen)
                            if ((.not.lte).and.debugInfo) write(*,'(a,1p,e12.5)') "Ne / Ntot: ",grid%ne(i1,i2,i3)/nTot
                         endif
                      endif
                   endif
                   deallocate(xall)
                enddo
             enddo
             percentDone = 100.*(real(i1)/real(grid%nx))
             if (.not.lte) write(*,'(a,f5.1)') "Percentage complete: ",percentDone
             deallocate(departCoeffAll)

          enddo
          !$OMP END PARALLEL DO

       endif

       ! if we are only computing in 2D then compute the populations in the z=0 plane
       ! and remap the populations
       ! onto the full 3D grid via rotation about the x-axis

       if (twoD) then


          if (aboutX) then
             do iIter = 1, nIter
                do i1 = 1, grid%nx
                   departCoeff = 1.
                   do i2 = 1, grid%ny
                      rvec = vector(grid%xaxis(i1), grid%yaxis(i2), 0.)
                      call hunt(grid%zAxis, grid%nz, 0., i3)

                      if (.not.grid%inStar(i1,i2,i3)) then
                         if (modulus(rVec-grid%starpos1) /= 0.) then
                            if (.not.lte) then
                               departCoeff(1) = 1./(0.5*(1.-sqrt(max(0.,(1.-grid%rStar1**2/modulus(rVec-grid%starpos1)**2)))))
                            endif
                         endif
                         if (iIter == 1) then
                            do i = 1, maxlevels
                               x(i) = grid%n(i1,i2,i3,i) * dble(departCoeff(i))
                            enddo
                            x(maxlevels+1) = grid%ne(i1,i2,i3)
                            oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                            oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
                         else
                            oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                            oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
                            do i = 1, maxlevels
                               x(i) = grid%n(i1,i2,i3,i)
                            enddo
                            x(maxlevels+1) = grid%ne(i1,i2,i3)
                         endif


                         if (.not.lte) then
                            if (grid%geometry == "binary") then
                               call occultTest(grid, i1, i2, i3, grid%starPos1, grid%rStar1, &
                                    grid%starPos2, grid%rStar2, visFrac1)
                               call occultTest(grid, i1, i2, i3, grid%starPos2, grid%rStar2, &
                                    grid%starPos1, grid%rStar1, visFrac2)
                            else
                               visFrac1 = 1.
                               visFrac2 = 0.
                            endif

                            call mnewt_stateq(grid, maxLevels+1,x,maxlevels+1,tolx,tolf, hnu1, &
                                 nuarray1, nnu1, hnu2, nuarray2, nnu2, rvec, i1, i2, i3, visFrac1, visFrac2, isBinary)
                         endif
                         if (.not.lte) write(*,*) "Grid: ",i1,i2,i3
                         do i = 1 , maxLevels
                            departCoeff(i) = real(x(i))/boltzSaha(i, grid%Ne(i1,i2,i3), dble(grid%temperature(i1,i2,i3)))

                            if (.not.lte) write(*,'(i3,1p,e12.3,e12.3)') i,departCoeff(i), x(i)/oldLevels(i)
                            grid%n(i1,i2,i3,i) = x(i)
                         enddo
                         grid%ne(i1,i2,i3) = x(maxLevels+1)
                         write(*,'(a,1p,e12.3,e12.3)') "Ne ",grid%ne(i1,i2,i3),grid%ne(i1,i2,i3)/oldLevels(maxLevels+1)
                      endif
                   enddo
                   percentDone = 100.*(real(i1)/real(grid%nx))
                   if (.not.lte) write(*,'(a,f5.1)') "Percentage complete: ",percentDone
                enddo
             enddo

             ! now perform the remapping

             do i = 1, grid%nx
                do j = 1, grid%ny
                   do k = 1, grid%nz
                      if (k /= i3) then
                         rVec = vector(grid%xaxis(i), grid%yaxis(j), grid%zaxis(k))
                         ang = atan2(rVec%z, rVec%y)
                         thisVec = rotateX(rVec, ang)
                         call getIndices(grid, thisVec, i1, i2, i3, t1, t2, t3)
                         if (.not.grid%inStar(i1,i2,i3)) then
                            grid%ne(i,j,k) = grid%ne(i1,i2,i3)
                            do m = 1, maxLevels
                               grid%n(i,j,k,m) = grid%n(i1,i2,i3,m)
                            enddo
                         endif
                      endif
                   enddo
                enddo
             enddo

          endif

          if (aboutZ) then
             if (.not.grid%cartesian) then


                WRITE(*,*) "about Z in polar"
                do iIter = 1, nIter
                   !$OMP PARALLEL DO &
                   !$OMP DEFAULT(NONE) &
                   !$OMP PRIVATE(i1, i2, i3) &
                   !$OMP SHARED(grid) &
                   !$OMP SHARED(lte, hnu1, hnu2, nuArray1, nnu1, nuarray2, nnu2) &
                   !$OMP SHARED(isBinary, debugInfo, iiter) &
                   !$OMP PRIVATE(i, visFrac1, visFrac2, rVec, nTot, percentDone) &
                   !$OMP PRIVATE(departCoeffAll, xAll, oldLevels, sinTheta, ang) 
                   do i1 = 1, grid%nr
                      allocate(departCoeffAll(1:maxlevels))
                      allocate(xall(1:maxLevels+1))
                      do i2 = 1, grid%nmu
                         departCoeffAll = 1.d0
                         ang = 0.
                         sinTheta = sqrt(1.d0 - grid%muAxis(i2)**2)
                         rvec = vector(grid%raxis(i1)*sinTheta*cos(ang), &
                              grid%raxis(i1)*sinTheta*sin(ang), &
                              grid%rAxis(i1)*grid%muAxis(i2))
                         call hunt(grid%phiAxis, grid%nphi, ang, i3)



                         if (.not.grid%inStar(i1,i2,i3).and.grid%inUse(i1,i2,i3)) then
                            if (iIter == 1) then
                               do i = 1, maxlevels
                                  xall(i) = grid%n(i1,i2,i3,i) * dble(departCoeffAll(i))
                               enddo
                               xall(maxlevels+1) = grid%ne(i1,i2,i3)
                               oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                               oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
                            else
                               oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                               oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
                               do i = 1, maxlevels
                                  xall(i) = grid%n(i1,i2,i3,i)
                               enddo
                               xall(maxlevels+1) = grid%ne(i1,i2,i3)
                            endif


                            if (.not.lte) then
                               if (grid%geometry == "binary") then
                                  call occultTest(grid, i1, i2, i3, grid%starPos1, grid%rStar1, &
                                       grid%starPos2, grid%rStar2, visFrac1)
                                  call occultTest(grid, i1, i2, i3, grid%starPos2, grid%rStar2, &
                                       grid%starPos1, grid%rStar1, visFrac2)
                               else
                                  visFrac1 = 1.
                                  visFrac2 = 0.
                               endif


                               call mnewt_stateq(grid, maxLevels+1,xall,maxlevels+1,tolx,tolf, hnu1, &
                                    nuarray1, nnu1, hnu2, nuarray2, nnu2, rvec, i1, i2, i3, visFrac1, visFrac2, isBinary)
                            endif
                            if ((.not.lte).and.debuginfo) write(*,*) "Grid: ",i1,i2,i3,grid%temperature(i1,i2,i3)
                            do i = 1 , maxLevels
                               departCoeffall(i) = real(xall(i))/boltzSaha(i, xall(maxlevels+1), dble(grid%temperature(i1,i2,i3)))

                               if ((.not.lte).and.debugInfo) then
                                  write(*,'(i3,1p,e12.3,e12.3, e12.3)') &
                                       i,departCoeffall(i), xall(i), abs(xall(i)-oldLevels(i))/xall(i)
                               endif
                               grid%n(i1,i2,i3,i) = xall(i)
                            enddo
                            grid%ne(i1,i2,i3) = xall(maxLevels+1)
                            if ((.not.lte).and.debugInfo) then
                               write(*,'(a,1p,e12.3,e12.3)') "Ne ",grid%ne(i1,i2,i3),grid%ne(i1,i2,i3)/oldLevels(maxLevels+1)
                            endif
                         if (minval(xall(1:maxlevels)) < 0.) then
                            write(*,'(a)') "Negative population..."
                            stop
                         endif
                      endif
                   enddo
                   percentDone = 100.*(real(i1)/real(grid%nr))
                   if ((.not.lte).and.debugInfo) write(*,'(a,f5.1)') "Percentage complete: ",percentDone
                   deallocate(departCoeffAll)
                   deallocate(xall)

                enddo
!$OMP END PARALLEL DO

             enddo

          else

             do iIter = 1, nIter
                do i1 = 1, grid%nx
                   departCoeff = 1.
                   do i3 = 1, grid%nz
                      rvec = vector(grid%xaxis(i1), 0., grid%zaxis(i3))
                      call hunt(grid%yAxis, grid%ny, 0., i2)

                      if (.not.grid%inStar(i1,i2,i3).and.grid%inUse(i1,i2,i3)) then
                         departCoeff = 1. !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                         if (iIter == 1) then
                            do i = 1, maxlevels
                               x(i) = grid%n(i1,i2,i3,i) * dble(departCoeff(i))
                            enddo
                            x(maxlevels+1) = grid%ne(i1,i2,i3)
                            oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                            oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
                         else
                            oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                            oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
                            do i = 1, maxlevels
                               x(i) = grid%n(i1,i2,i3,i)
                            enddo
                            x(maxlevels+1) = grid%ne(i1,i2,i3)
                         endif


                         if (.not.lte) then
                            if (grid%geometry == "binary") then
                               call occultTest(grid, i1, i2, i3, grid%starPos1, grid%rStar1, &
                                    grid%starPos2, grid%rStar2, visFrac1)
                               call occultTest(grid, i1, i2, i3, grid%starPos2, grid%rStar2, &
                                    grid%starPos1, grid%rStar1, visFrac2)
                            else
                               visFrac1 = 1.
                               visFrac2 = 0.
                            endif

                            call mnewt_stateq(grid, maxLevels+1,x,maxlevels+1,tolx,tolf, hnu1, &
                                 nuarray1, nnu1, hnu2, nuarray2, nnu2, rvec, i1, i2, i3, visFrac1, visFrac2, isBinary)
                         endif
                         if (.not.lte) write(*,*) "Grid: ",i1,i2,i3,grid%temperature(i1,i2,i3)
                         do i = 1 , maxLevels
                            departCoeff(i) = real(x(i))/boltzSaha(i, x(maxlevels+1), dble(grid%temperature(i1,i2,i3)))

                            if (.not.lte) then
                               write(*,'(i3,1p,e12.3,e12.3)') &
                                    i,departCoeff(i), x(i)
                            endif
                            grid%n(i1,i2,i3,i) = x(i)
                         enddo
                         grid%ne(i1,i2,i3) = x(maxLevels+1)
                         if (.not.lte) then
                            write(*,'(a,1p,e12.3,e12.3)') "Ne ",grid%ne(i1,i2,i3),grid%ne(i1,i2,i3)/oldLevels(maxLevels+1)
                         endif
                      endif
                   enddo
                   percentDone = 100.*(real(i1)/real(grid%nx))
                   if (.not.lte) write(*,'(a,f5.1)') "Percentage complete: ",percentDone
                enddo
             enddo
          endif
             
             ! now perform the remapping

          if (.not.lte) then
          if (grid%cartesian) then
             do i = 1, grid%nx
                do j = 1, grid%ny
                   do k = 1, grid%nz
                      if (j /= i2) then
                         rVec = vector(grid%xaxis(i), grid%yaxis(j), grid%zaxis(k))
                         ang = atan2(rVec%y, rVec%x)
                         thisVec = rotateZ(rVec, ang)
                         call getIndices(grid, thisVec, i1, i2, i3, t1, t2, t3)
                         if (.not.grid%inStar(i1,i2,i3)) then
                            grid%ne(i,j,k) = grid%ne(i1,i2,i3)
                            do m = 1, maxLevels
                               grid%n(i,j,k,m) = grid%n(i1,i2,i3,m)
                            enddo
                         endif
                      endif
                   enddo
                enddo
             enddo
             else

                write(*,'(a)') "Remapping about Z axis..."
                do k = 2, grid%nphi
                   grid%ne(1:grid%nr,1:grid%nMu,k) = &
                        grid%ne(1:grid%nr,1:grid%nMu,1)
                   grid%n(1:grid%nr,1:grid%nmu,k,1:maxLevels) = &
                        grid%n(1:grid%nr,1:grid%nMu,1,1:maxLevels)
                enddo
             endif
       endif
    endif
    
 endif



 if (oneD) then
    
    nIter = 1.
    
    departCoeff = 1.
    do iIter = 1, nIter
       do i1 = 1,grid%na1
          i2 = 1
          i3 = 1
          sinTheta = sqrt(1.d0 - grid%muAxis(i2)**2)
          rvec = vector(grid%raxis(i1)*sinTheta*cos(ang), &
               grid%raxis(i1)*sinTheta*sin(ang), &
               grid%rAxis(i1)*grid%muAxis(i2))
          
          where (departCoeff < 0.) 
             departCoeff = 1.
          end where
          departCoeff(1) = 1./(0.5*(1.-sqrt(max(0.,(1.-grid%rStar1**2/modulus(rVec-grid%starpos1)**2)))))
          
          
          if (.not.grid%inStar(i1,i2,i3).and.grid%inUse(i1,i2,i3)) then
             
             if (iIter == 1) then
                do i = 1, maxlevels
                   x(i) = grid%n(i1,i2,i3,i) * dble(departCoeff(i))
                enddo
                x(maxlevels+1) = grid%ne(i1,i2,i3)
                oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
             else
                oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
                do i = 1, maxlevels
                   x(i) = grid%n(i1,i2,i3,i)
                enddo
                x(maxlevels+1) = grid%ne(i1,i2,i3)
             endif
             
             
             if (.not.lte) then
                if (grid%geometry == "binary") then
                   call occultTest(grid, i1, i2, i3, grid%starPos1, grid%rStar1, &
                        grid%starPos2, grid%rStar2, visFrac1)
                   call occultTest(grid, i1, i2, i3, grid%starPos2, grid%rStar2, &
                        grid%starPos1, grid%rStar1, visFrac2)
                else
                   visFrac1 = 1.
                   visFrac2 = 0.
                endif
                
                call mnewt_stateq(grid, maxLevels+1,x,maxlevels+1,tolx,tolf, hnu1, &
                     nuarray1, nnu1, hnu2, nuarray2, nnu2, rvec, i1, i2, i3, visFrac1, visFrac2, isBinary)
             endif
             if (.not.lte) write(*,*) "Grid: ",i1,i2,i3,grid%temperature(i1,i2,i3)
             do i = 1 , maxLevels
                departCoeff(i) = real(x(i))/boltzSaha(i, x(maxlevels+1), dble(grid%temperature(i1,i2,i3)))
                
                if (.not.lte) then
                   write(*,'(i3,1p,e12.3,e12.3,e12.3,e12.3)') &
                        i,departCoeff(i), x(i),log10(departCoeff(i)),log10(x(2)/x(1))
                endif
                grid%n(i1,i2,i3,i) = x(i)
             enddo
             grid%ne(i1,i2,i3) = x(maxLevels+1)
             if (.not.lte) then
                write(*,'(a,1p,e12.3,e12.3)') "Ne ",grid%ne(i1,i2,i3),grid%ne(i1,i2,i3)/grid%nTot(i1,i2,i3)
             endif

             endif
             enddo
          enddo
             ! now perform the remapping
             rHat = VECTOR(1., 1., 1.)
             call normalize(rHat)

             do i = 1, grid%na1
                do j = 1, grid%na2
                   do k = 1, grid%na3
                      i1 = i
                      i2 = 1
                      i3 = 1
                      if (.not.grid%inStar(i1,i2,i3)) then
                         grid%ne(i,j,k) = grid%ne(i1,i2,i3)
                         do m = 1, maxLevels
                            grid%n(i,j,k,m) = grid%n(i1,i2,i3,m)
                         enddo
                      endif
                   enddo
                enddo
             enddo

          endif







          if (writePops) then
             call writeGridPopulations(popFilename, grid, maxLevels)
          endif
    endif

    write(*,'(a,f8.1)') "Generating opacities for ",lambdaTrans(nLower, nUpper)*1.e8


!    write(*,*) "!!!!!!!!!!!!!!! depart coeff of level 3 is 2"
!    grid%n(:,:,:,3) = grid%n(:,:,:,3) * 100.
!    grid%n(:,:,:,5) = grid%n(:,:,:,5) * 30.
    call generateOpacities(grid, nLower, nUpper)

    deallocate(grid%n)
    deallocate(grid%ne)
    deallocate(grid%nTot)



  end subroutine initGridStateq



  real(double) elemental function BoltzSaha(m, Ne, t)
  
    integer, intent(in)              :: m
    real(double), intent(in):: Ne, t
    
    real(double), parameter :: ci = 2.07d-16
    real(double), parameter :: iPot = hydE0eV

    BoltzSaha = Ne**2 * gDegen(m) * ci * exp( (iPot-eTrans(m))/(kev * t) ) / t**1.5
    
  end function BoltzSaha


  real(double) pure function boltzmann(m, N0, t)
  
    integer, intent(in)              :: m
    real(double), intent(in):: N0, t
    real(double)            :: z0

    z0 = SUM(gDegen*exp(-eTrans/(kev*t)))
    boltzmann = (gDegen(m)/z0)*n0*exp(-eTrans(m)/(kev*t))

  end function boltzmann


!  real(double) pure function oldcikt(i,t)
!    !
!    ! this function calculates the collisional ionization rate (see k&c and
!    ! mihalas 1967 (apj 149 169) ).
!    ! 
!    integer, intent(in)              :: i        ! the level
!    real(double),intent(in) :: t        ! the temperature
!    
!    real(double) :: t1                  
!    real(double) :: gammait             ! see k&c
!    real(double) :: lgt
!    real(double) :: chi                 ! potential from i to k
!    ! making cint a PARAMETER may cause problems with XL Fortran
!    !  real(double) :: cint(5,10)
!    !  cint = reshape(source=                                                                   &
!    real(double), parameter :: cint(5,10) = reshape(source=                         &
!         (/-0.4350000e0_db, 0.3000000e00_db,0.0000000e0_db, 0.0000000e0_db,0.00000000e0_db,  &
!            1.9987261e1_db,-5.8906298e-5_db,0.0000000e0_db,-2.8185937e4_db,5.44441600e7_db,  &
!            1.3935312e3_db,-1.6805859e02_db,0.0000000e0_db,-2.5390000e3_db,0.00000000e0_db,  &
!            2.0684609e3_db,-3.3415820e02_db,0.0000000e0_db, 0.0000000e0_db,-7.6440625e3_db,  &
!            3.2174844e3_db,-5.5882422e02_db,0.0000000e0_db, 0.0000000e0_db,-6.8632500e3_db,  &
!            5.7591250e3_db,-1.5163125e03_db,8.1750000e1_db, 0.0000000e0_db,0.00000000e0_db,  &
!            1.4614750e4_db,-4.8283750e03_db,3.9335938e2_db, 0.0000000e0_db,0.00000000e0_db,  &
!            2.8279250e4_db,-1.0172750e04_db,9.1967968e2_db, 0.0000000e0_db,0.00000000e0_db,  &
!            4.6799250e4_db,-1.7627500e04_db,1.6742031e3_db, 0.0000000e0_db,0.00000000e0_db,  &
!           -7.4073000e4_db, 6.8599375e03_db,0.0000000e0_db, 2.0169800e5_db,0.00000000e0_db/),&
!                                                                             shape=(/5,10/))
!
!    t1 = min(t,1.5e5_db)
!    lgt=log10(t1)
!    if (i .ne. 2) then 
!       gammait=cint(1,i)+cint(2,i)*lgt+cint(3,i)*(lgt**2)+ &
!            (cint(4,i)/lgt)+(cint(5,i)/(lgt**2))
!    else
!       gammait=cint(1,i)+cint(2,i)*t1+(cint(4,i)/t1)+(cint(5,i)/t1**2)
!    endif
!    chi=hydE0eV-eTrans(i)
!    oldcikt=(5.465e-11)*sqrt(t1)*exp(-chi/(kev*t1))*gammait
!
!  end function oldcikt
  
  real(double) function cikt(i,t)
    !
    ! this function calculates the collisional ionization rate
    ! (using Hillier coefficients for most levels).
    ! 
    integer, intent(in)              :: i        ! the level
    real(double),intent(in) :: t        ! the temperature
    
    real(double) :: t1                  
    real(double) :: gammait             ! see k&c
    real(double) :: lgt
    real(double) :: chi                 ! potential from i to k
    integer :: lower, upper
    real(double) :: factor
    
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
              exp(-chi/(kev*t1)) / gDegen(i) / sqrt(t1*10e-4_db)
        
!    else
!       gammait=cint(1,i)+cint(2,i)*t1+(cint(4,i)/t1)+(cint(5,i)/t1**2)
!       cikt=(5.465e-11)*sqrt(t1)*exp(-chi/(kev*t1))*gammait
!    endif

  end function cikt


!  real(double) pure function oldcijt(i,j,t)
!    !
!    ! this function calculates the excitation rate according to crandall et al
!    ! (apj 191 789).
!    !
!    integer, intent(in)              :: i,j  ! lower/upper level
!    real(double),intent(in) :: t    ! temperature
!    
!    integer :: r                             ! counter
!    real(double) :: t1              ! temperature
!    real(double) :: bsum(9)         ! summation
!    real(double) :: chi             ! excitation energy
!    real(double) :: x
!    integer               :: i1,j1           ! lower/upper level
!    ! making avals a PARAMETER may cause problems with XL Fortran
!    !   real(double) :: avals(7)     
!    !   avals = (/ 2.579997e-10_db, -1.629166e-10_db, 7.713069e-11_db, -2.668768e-11_db, &
!    !              6.642513e-12_db, -9.422885e-13_db, 0.e0_db                            /)
!    real(double),parameter :: avals(7) =                                    &    
!            (/ 2.579997e-10_db, -1.629166e-10_db, 7.713069e-11_db, -2.668768e-11_db, &
!               6.642513e-12_db, -9.422885e-13_db, 0.e0_db                            /)
!
!    t1 = min(t,1.e5_db)
!    i1 = min(i,j)
!    j1 = max(i,j)
!
!    x=log10(t1)-4.e0
!    chi=abs(eTrans(j)-eTrans(i))
!    chi=chi/(kev*t1)
!    if ((i1.eq.1).and.(j1.eq.2)) then
!       bsum(9)=0.0
!       bsum(8)=0.0
!       do r=7,1,-1
!          bsum(r)=2.0*x*bsum(r+1)-bsum(r+2)+avals(r)
!       enddo
!       oldcijt=4.e0*sqrt(t1)*exp(-chi)*0.5e0*(bsum(1)-bsum(3))
!    else
!       oldcijt=abs(5.465e-11_db*sqrt(t1)*4.e0_db*fStrength(i1,j1)/ &
!            (( (1.e0/real(i1)**2) - (1.0/real(j1)**2))**2))*chi* &
!            (expint(1,chi)+0.148e0_db*chi*expint(5,chi))
!    endif
!
!    if (j < i) then
!       oldcijt = oldcijt / (exp(-chi) * 2.)
!    endif
!
!
!  end function oldcijt
  
  real(double) function cijt(i,j,t)
    implicit none
    !
    ! this function calculates the excitation rate using values from Hillier !
    integer, intent(in)              :: i,j  ! lower/upper level
    real(double),intent(in) :: t    ! temperature
    
    real(double) :: t1              ! temperature
    real(double) :: chi             ! excitation energy
    real(double) :: bsum(9)         ! summation
    integer               :: i1,j1           ! lower/upper level
    integer :: lower, upper
    real(double) :: factor
    real(double),parameter :: avals(7) =                                    &    
            (/ 2.579997e-10_db, -1.629166e-10_db, 7.713069e-11_db, -2.668768e-11_db, &
               6.642513e-12_db, -9.422885e-13_db, 0.e0_db                            /)
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
    exp(-chi) / gDegen(i1) / sqrt(t1*10e-4_db)

      
    if (j < i) then
       cijt = cijt / (exp(-chi) * 2.)
    endif


  end function cijt


  real(double) pure function annu(n,nu)
    !
    ! this subroutine calculates the photoionization x-section
    ! for hydrogen from the n-th level for a given freq photon.
    !
    integer, intent(in)              :: n     ! the level
    real(double), intent(in):: nu    ! the photon frequency
    real(double)            :: lam,e ! the photon wavelength
    
    lam=cSpeed/nu
    lam=lam * 1.e8_db

    e = hCgs * nu * ergToEv

    if (e > (hydE0eV - eTrans(n))) then
       annu=1.044e-26_db*gii(n,1.e0_db,lam)*(lam**3)/dble(n)**5
    else
       annu = 1.e-30_db
    endif

  end function annu


  real pure function expint(n,x1)
    ! 
    ! exponential integrals. numerical approximations given by gray.
    !
    integer, intent(in)              :: n
    real(double),intent(in) :: x1
    real(double)            :: x,ep(5)
    integer                          :: i
    ! making a1 and b1 PARAMETERs may cause problems with XL Fortran
    real(double),parameter  :: a1(4) =                                      &
       (/ 0.2677737343e0_db,8.6347608925e0_db, 18.059016973e0_db, 8.5733287401e0_db /)
    real(double),parameter  :: b1(4) =                                      &
       (/ 3.9584969228e0_db,21.0996530827e0_db,25.6329561486e0_db,9.5733223454e0_db /)

    x=abs(x1)
    if (x.le.1.0) then
       ep(1)=-log(x)-0.57721566e0_db+0.99999193e0_db*x-0.24991055e0_db*x**2+ &
            0.05519968e0_db*x**3-0.00976004_db*x**4+0.00107857_db*x**5
    else
       ep(1)=(x**4+a1(4)*x**3+a1(3)*x**2+a1(2)*x+a1(1))
       ep(1)=ep(1)/(x**4+b1(4)*x**3+b1(3)*x**2+b1(2)*x+b1(1))
       ep(1)=ep(1)/(x*exp(x))
    endif
    if (n.gt.1) then
       do i=1,n-1
          ep(i+1)=(exp(-x)-x*ep(i))/real(i)
       enddo
    endif
    expint=ep(n)
  end function expint



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
  

  real pure function giia(nDouble,z,wl)
    real(double), intent(in) :: nDouble
    real(double), intent(in) :: z, wl

    real(double)             :: ag, u, gii, term
    real(double), parameter  :: twoThirds = 2.e0_db/3.e0_db
    real(double), parameter  :: fourThirds = 2.0_db * twoThirds 

    u = nDouble*nDouble*911.76e0_db / (wl*z*z)-1.e0_db
    gii = 1.e0_db + 0.1728e0_db * (u-1.e0_db) / ((nDouble*(u+1.e0_db))**twoThirds)
    term = 0.0496e0_db * (1.e0_db+u*(u+1.333e0_db))/(nDouble*(u+1.e0_db)**(fourThirds))
    if ((term/gii).le.0.25e0_db) then
       giia=gii-term
       return
    else
       giia = 1.0_db
       return
    endif
    
  end function giia


  real(double) function equation8(n, nPop, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, &
       rVec, i1, i2, i3, grid, visFrac1, visFrac2, binary, thisOctal, thisSubcell)

    use input_variables, only : LyContThick
    type(GRIDTYPE), intent(in)  :: grid
    integer, intent(in)         :: nPop, nNu1, nNu2
    logical, intent(in)         :: binary
    real, intent(in)            :: visFrac1, visFrac2
    real, intent(in), dimension(:) :: Hnu1, nuArray1
    real, intent(in), dimension(:) :: Hnu2, nuArray2
    integer, intent(in)         :: i1, i2, i3
    type(VECTOR), intent(in)    :: rVec
    type(OCTAL),pointer,optional:: thisOctal
    integer,intent(in),optional :: thisSubcell
    integer :: m,n
    real(double) :: nstar
    real(double) :: freq
    integer, parameter :: debug = 0
    integer, save :: i = 0
    real(double) :: Inu
    real(double) :: tot
    real(double) :: fac1, fac2

    if (n == debug) write(*,*) " "

    ! we have one IF statement to decide whether we use the block containin the
    !   AMR version of the code. 
   if (grid%adaptive) then
   
       tot = 0.e0_db

       do m = 1, n - 1
          freq = cSpeed/lambdaTrans(m,n)
          if ((freq >= nuArray1(1)).and.(freq <= nuArray1(nnu1))) then
             call hunt(nuArray1, nNu1, real(freq), i)
             Inu = 4.*logint(real(freq), nuArray1(i),nuArray1(i+1),hnu1(i), hnu1(i+1))
          else
             inu = 0.
          endif
   
          tot = tot  + (thisOctal%N(thisSubcell,m)*bEinstein(m,n) &
               - thisOctal%N(thisSubcell,n)*bEinstein(n,m)) * &
               beta_cmn(m, n, rVec, i1, i2, i3, grid, 1, thisOctal, thisSubcell) * Inu
   
          if (binary) then
             if ((freq >= nuArray2(1)).and.(freq <= nuArray2(nnu2))) then
                call hunt(nuArray2, nNu2, real(freq), i)
                Inu = 4.*logint(real(freq), nuArray2(i),nuArray2(i+1),hnu2(i), hnu2(i+1))
             else
                inu = 0.
             endif
             tot = tot  + (thisOctal%N(thisSubcell,m)*bEinstein(m,n) - &
                  thisOctal%N(thisSubcell,n)*bEinstein(n,m)) * &
                  beta_cmn(m, n, rVec, i1, i2, i3, grid, 2, thisOctal, thisSubcell) * Inu
          endif
   
          tot = tot - thisOctal%N(thisSubcell,n)*aEinstein(n,m)*&
                   beta_mn(m, n, rVec, i1, i2, i3, grid, thisOctal, thisSubcell)
          tot = tot + thisOctal%Ne(thisSubcell) * &
               (thisOctal%N(thisSubcell,m)*cijt(m,n,real(thisOctal%temperature(thisSubcell),kind=db)) &
               - thisOctal%N(thisSubcell,n) * cijt(n,m,real(thisOctal%temperature(thisSubcell),kind=db)))
       enddo
   
       if (n == debug) then
          write(*,*) "m < n",tot
       endif
   
   
       do m = n+1, nPop
          tot = tot + thisOctal%N(thisSubcell,m)*aEinstein(m,n)*&
                   beta_mn(n,m,rVec, i1, i2, i3,grid, thisOctal, thisSubcell)
          if (n == debug) then
             write(*,*) "spont",m,tot
          endif
   
          freq = cSpeed/lambdaTrans(m,n)
          if ((freq >= nuArray1(1)).and.(freq <= nuArray1(nnu1))) then
             call hunt(nuArray1, nNu1, real(freq), i)
             Inu = 4.*logint(real(freq), nuArray1(i),nuArray1(i+1),hnu1(i), hnu1(i+1))
          else
             inu = 0.
          endif

          tot = tot  - (thisOctal%N(thisSubcell,n)*bEinstein(n,m) - &
               thisOctal%N(thisSubcell,m)*bEinstein(m,n)) * &
               beta_cmn(n, m, rVec, i1, i2, i3, grid, 1, thisOctal, thisSubcell) * Inu
   
          if (n == debug) then
             write(*,*) m,"star 1",tot
          endif
   
   
          if (binary) then
             if ((freq >= nuArray2(1)).and.(freq <= nuArray2(nnu2))) then
                call hunt(nuArray2, nNu2, real(freq), i)
                Inu = 4.*logint(real(freq), nuArray2(i),nuArray2(i+1),hnu2(i), hnu2(i+1))
             else
                inu = 0.
             endif
             tot = tot  -(thisOctal%N(thisSubcell,n)*bEinstein(n,m) - &
                  thisOctal%N(thisSubcell,m)*bEinstein(m,n)) * &
                  beta_cmn(n, m, rVec, i1, i2, i3, grid, 2, thisOctal, thisSubcell) * Inu
          endif
   
          if (n == debug) then
             write(*,*) m,"star 2",tot
          endif
   
   
          tot = tot + thisOctal%Ne(thisSubcell) * (thisOctal%N(thisSubcell,m) &
           *cijt(m,n,real(thisOctal%temperature(thisSubcell),kind=db)) - thisOctal%N(thisSubcell,n)&
           *cijt(n,m,real(thisOctal%temperature(thisSubcell),kind=db)))
          if (n == debug) then
             write(*,*) m,"collisional",tot
          endif
   
   
       enddo
   
   
       if ( n == debug) then
          write(*,*) "m > n",tot
       endif

       NStar = boltzSaha(n, thisOctal%Ne(thisSubcell),dble(thisOctal%temperature(thisSubcell)))

       fac1 = integral1(n,hnu1, nuArray1, nNu1, rVec, i1, i2, i3, grid, 1, thisOctal, thisSubcell)*visFrac1 
       if (binary) then
          fac2 = integral1(n,hnu2, nuArray2, nNu2, rVec, i1, i2, i3, grid, 2, thisOctal, thisSubcell)*visFrac2
       else
          fac2 = 0.
       endif
       if (.not. (n == 1 .and. LyContThick)) then
          tot = tot + Nstar * (fac1 + fac2 + thisOctal%Ne(thisSubcell)* &
          cikt(n,real(thisOctal%temperature(thisSubcell),kind=db)))
       end if
   
       if (n == debug) then
          write(*,*) "Recombination: ",tot
          write(*,*) nstar, fac1, fac2, thisOctal%Ne(thisSubcell)* &
             cikt(n,real(thisOctal%temperature(thisSubcell),kind=db))
       endif
       
       
       fac1 = integral2(n,hnu1, nuArray1, nNu1, rVec, grid, 1)*visFrac1 
 
       if (binary) then
          fac2 = integral2(n,hnu2, nuArray2, nNu2, rVec, grid, 2)*visFrac2
       else
          fac2 = 0.
       endif
       
       if (.not. (n == 1 .and. LyContThick)) then
         tot = tot - thisOctal%N(thisSubcell,n)*(fac1 + fac2 + thisOctal%Ne(thisSubcell) * &
          cikt(n, real(thisOctal%temperature(thisSubcell),kind=db)))
       end if
   
       if (n == debug) then
          write(*,*) "Ionization: ", &
               thisOctal%N(thisSubcell,n)*(fac1 + fac2 + thisOctal%Ne(thisSubcell) * &
                  cikt(n, real(thisOctal%temperature(thisSubcell),kind=db)))
          write(*,*) thisOctal%N(thisSubcell,n),fac1,fac2,thisOctal%Ne(thisSubcell) * &
             cikt(n, real(thisOctal%temperature(thisSubcell),kind=db))
       endif 
       
       
    else ! grid is not adaptive
       
       
       tot = 0.e0_db
       do m = 1, n - 1
          freq = cSpeed/lambdaTrans(m,n)
          if ((freq >= nuArray1(1)).and.(freq <= nuArray1(nnu1))) then
             call hunt(nuArray1, nNu1, real(freq), i)
             Inu = 4.*logint(real(freq), nuArray1(i),nuArray1(i+1),hnu1(i), hnu1(i+1))
          else
             inu = 0.
          endif
   
   
          tot = tot  + (grid%N(i1,i2,i3,m)*bEinstein(m,n) &
               - grid%N(i1,i2,i3,n)*bEinstein(n,m)) * &
               beta_cmn(m, n, rVec, i1, i2, i3, grid, 1) * Inu
   
          if (binary) then
             if ((freq >= nuArray2(1)).and.(freq <= nuArray2(nnu2))) then
                call hunt(nuArray2, nNu2, real(freq), i)
                Inu = 4.*logint(real(freq), nuArray2(i),nuArray2(i+1),hnu2(i), hnu2(i+1))
             else
                inu = 0.
             endif
             tot = tot  + (grid%N(i1,i2,i3,m)*bEinstein(m,n) - &
                  grid%N(i1,i2,i3,n)*bEinstein(n,m)) * &
                  beta_cmn(m, n, rVec, i1, i2, i3, grid, 2) * Inu
          endif
   
          tot = tot - grid%N(i1,i2,i3,n)*aEinstein(n,m)*beta_mn(m, n, rVec, i1, i2, i3, grid)
          tot = tot + grid%Ne(i1,i2,i3) * &
               (grid%N(i1,i2,i3,m)*cijt(m,n,dble(grid%temperature(i1,i2,i3))) - grid%N(i1,i2,i3,n) * &
                cijt(n,m,dble(grid%temperature(i1,i2,i3))))
       enddo
   
       if (n == debug) then
          write(*,*) "m < n",tot
       endif
   
   
       do m = n+1, nPop
          tot = tot + grid%N(i1,i2,i3,m)*aEinstein(m,n)*beta_mn(n,m,rVec, i1, i2, i3,grid)
          if (n == debug) then
             write(*,*) "spont",m,tot
          endif
   
          freq = cSpeed/lambdaTrans(m,n)
          if ((freq >= nuArray1(1)).and.(freq <= nuArray1(nnu1))) then
             call hunt(nuArray1, nNu1, real(freq), i)
             Inu = 4.*logint(real(freq), nuArray1(i),nuArray1(i+1),hnu1(i), hnu1(i+1))
          else
             inu = 0.
          endif
   
   
          tot = tot  - (grid%N(i1,i2,i3,n)*bEinstein(n,m) - &
               grid%N(i1,i2,i3,m)*bEinstein(m,n)) * &
               beta_cmn(n, m, rVec, i1, i2, i3, grid, 1) * Inu
   
          if (n == debug) then
             write(*,*) m,"star 1",tot
          endif
   
   
          if (binary) then
             if ((freq >= nuArray2(1)).and.(freq <= nuArray2(nnu2))) then
                call hunt(nuArray2, nNu2, real(freq), i)
                Inu = 4.*logint(real(freq), nuArray2(i),nuArray2(i+1),hnu2(i), hnu2(i+1))
             else
                inu = 0.
             endif
             tot = tot  -(grid%N(i1,i2,i3,n)*bEinstein(n,m) - &
                  grid%N(i1,i2,i3,m)*bEinstein(m,n)) * &
                  beta_cmn(n, m, rVec, i1, i2, i3, grid, 2) * Inu
          endif
   
          if (n == debug) then
             write(*,*) m,"star 2",tot
          endif
   
   
          tot = tot + grid%Ne(i1,i2,i3) * (grid%N(i1,i2,i3,m) &
           *cijt(m,n,dble(grid%temperature(i1,i2,i3))) - grid%N(i1,i2,i3,n)*cijt(n,m,dble(grid%temperature(i1,i2,i3))))
          if (n == debug) then
             write(*,*) m,"collisional",tot
          endif
   
   
       enddo
   
   
       if ( n == debug) then
          write(*,*) "m > n",tot
       endif
   
       NStar = boltzSaha(n, grid%Ne(i1,i2,i3),dble(grid%temperature(i1,i2,i3)))
   
       fac1 = integral1(n,hnu1, nuArray1, nNu1, rVec, i1, i2, i3, grid, 1, thisOctal, thisSubcell)*visFrac1
       if (binary) then
          fac2 = integral1(n,hnu2, nuArray2, nNu2, rVec, i1, i2, i3, grid, 2, thisOctal, thisSubcell)*visFrac2
       else
          fac2 = 0.
       endif
       tot = tot + Nstar * (fac1 + fac2 + grid%Ne(i1,i2,i3)*cikt(n,dble(grid%temperature(i1,i2,i3))))
   
       if (n == debug) then
          write(*,*) "Recombination: ",tot
          write(*,*) nstar, fac1, fac2, grid%Ne(i1,i2,i3)*cikt(n,dble(grid%temperature(i1,i2,i3)))
       endif
       
       
   
       fac1 = integral2(n,hnu1, nuArray1, nNu1, rVec, grid, 1)*visFrac1
       if (binary) then
          fac2 = integral2(n,hnu2, nuArray2, nNu2, rVec, grid, 2)*visFrac2
       else
          fac2 = 0.
       endif
       tot = tot - grid%N(i1,i2,i3,n)*(fac1 + fac2 + grid%Ne(i1,i2,i3) * cikt(n, dble(grid%temperature(i1,i2,i3))))
   
       if (n == debug) then
          write(*,*) "Ionization: ", &
               grid%N(i1,i2,i3,n)*(fac1 + fac2 + grid%Ne(i1,i2,i3) * cikt(n, dble(grid%temperature(i1,i2,i3))))
          write(*,*) grid%N(i1,i2,i3,n),fac1,fac2,grid%Ne(i1,i2,i3) * cikt(n, dble(grid%temperature(i1,i2,i3)))
       endif 
       
    end if ! (grid%adaptive)

    if (n == debug) write(*,*) "Final total:",tot
    equation8 = tot


  end function equation8
 

  real (kind=double) pure function equation14(nPop, grid, i1, i2, i3, thisOctal, thisSubcell)

    type(GRIDTYPE), intent(in)   :: grid
    integer, intent(in)          :: i1, i2, i3, nPop
    type(OCTAL),pointer,optional :: thisOctal
    integer,intent(in),optional  :: thisSubcell
    
    real(double) :: tot
    integer :: i
    integer, parameter :: maxLevels = statEqMaxLevels

    if (grid%adaptive) then

      tot = 0.e0_db
       do i = 1, nPop
          tot = tot + dble(thisOctal%N(thisSubcell,i))
       enddo
       
       do i = nPop+1, maxLevels+3
          tot = tot + boltzSaha(i, dble(thisOctal%Ne(thisSubcell)), dble(thisOctal%temperature(thisSubcell)))
       enddo
   !    write(*,*) "Neutral",tot
       tot = tot + thisOctal%Ne(thisSubcell)
   !    write(*,*) "Ne",grid%ne(i1,i2,i3)
       tot = tot - dble(thisOctal%rho(thisSubcell) / mHydrogen)
   !    write(*,*) "Ntot",grid%rho(i1,i2,i3) / mHydrogen
       equation14 = tot
   !    write(*,*) "tot",tot

    else ! grid not adaptive

      tot = 0.e0_db
       do i = 1, nPop
          tot = tot + dble(grid%n(i1,i2,i3,i))
       enddo
       
       do i = nPop+1, maxLevels+3
          tot = tot + boltzSaha(i, grid%ne(i1,i2,i3), dble(grid%temperature(i1,i2,i3)))
       enddo
   !    write(*,*) "Neutral",tot
       tot = tot + grid%ne(i1,i2,i3)
   !    write(*,*) "Ne",grid%ne(i1,i2,i3)
       tot = tot - dble(grid%rho(i1,i2,i3) / mHydrogen)
   !    write(*,*) "Ntot",grid%rho(i1,i2,i3) / mHydrogen
       equation14 = tot
   !    write(*,*) "tot",tot

    end if ! (grid%adaptive)

  end function equation14



  real function integral1(n,hnu, nuArray, nNu, rVec, i1, i2, i3, grid, nStar, thisOctal, thisSubcell)
    real, intent(in), dimension(:) :: hnu, nuArray
    integer, intent(in)      :: nNu, nStar, n
    type(VECTOR), intent(in) :: rVec
    type(GRIDTYPE),intent(in):: grid
    integer,intent(in)       :: i1, i2, i3
    type(OCTAL),pointer,optional:: thisOctal
    integer,intent(in),optional :: thisSubcell
    
    real    :: r
    integer :: i, imin
    real    :: fac1, fac2
    real    :: w, x1, freq, tot, jnu

    if (nStar == 1) then
       r = modulus(rVec-grid%starPos1)
       x1 = sqrt(max(0.,(1. - grid%rStar1**2 / r**2)))
       w = 0.5*(1. - x1)
    else
       r = modulus(rVec-grid%starPos2)
       x1 = sqrt(max(0.,1. - grid%rStar2**2 / r**2))
       w = 0.5*(1. - x1)
    endif

    freq = ((hydE0eV-eTrans(n))*1.602192e-12)/hcgs
    if ((freq < nuArray(1)) .or.(freq > nuArray(nNu))) then
       write(*,*) " "
       write(*,*) "Error in [steteq_mod::integral1]. Freq out of range!"
       write(*,*) " =====>  nNu = ", nNu
       write(*,*) " =====> freq = ", freq
       write(*,*) " =====> nuArray(1) =  ", nuArray(1)
       write(*,*) " =====> nuArray(nNu) =  ", nuArray(nNu)
       write(*,*) " =====> You should try to fix this problem !"
       write(*,*) " =====> HINT: Add a point with low flux value at the end of "
       write(*,*) " =====>       your input flux data file."
       write(*,*) " "
       jnu = 1.e-28
       iMin = 1
       stop
    else
       iMin = -1
       call hunt(nuArray, nNu, freq, iMin)
       jnu = 4.*w*logint(freq,nuArray(iMin), nuArray(iMin+1), hnu(imin), hnu(imin+1))
    endif
    tot = 0.

    if (grid%adaptive) then 
            
       fac1  = (fourPi/(hCgs*freq))*annu(n,dble(freq))* &
            ((2.*real(dble(hCgs)*dble(freq)**3) / (cSpeed*cSpeed)) + jnu) &
            *exp(-hCgs*freq/(kErg*thisOctal%temperature(thisSubcell)))
       jnu = 4.*w*hnu(imin+1)
       fac2  = (fourPi/(hcgs*nuArray(imin+1)))*annu(n,dble(nuArray(imin+1)))* &
            ((2.*real(dble(hcgs)*dble(nuArray(imin+1))**3) / (cSpeed*cSpeed)) + jnu) &
            *exp(-hcgs*nuArray(imin+1)/(kErg*thisOctal%temperature(thisSubcell)))
       tot = tot + 0.5*(fac1 + fac2)*(nuArray(imin+1)-freq)
       do i = iMin+1, nNu-1
          Jnu = 4.*w*hnu(i)
          fac1  = (fourPi/(hCgs*nuArray(i)))*annu(n,dble(nuArray(i)))* &
               ((2.*real(dble(hCgs)*dble(nuArray(i))**3) / (cSpeed*cSpeed)) + jnu) &
               *exp(-hCgs*nuArray(i)/(kErg*thisOctal%temperature(thisSubcell)))
          Jnu = 4.*w*hnu(i+1)
          fac2  = (fourPi/(hcgs*nuArray(i+1)))*annu(n,dble(nuArray(i+1)))* &
               ((2.*real(dble(hcgs)*dble(nuArray(i+1))**3) / (cSpeed*cSpeed)) + jnu) &
               *exp(-hcgs*nuArray(i+1)/(kErg*thisOctal%temperature(thisSubcell)))
          tot = tot + 0.5*(fac1 + fac2)*(nuArray(i+1)-nuArray(i))
       enddo

    else ! not adaptive
            
       fac1  = (fourPi/(hCgs*freq))*annu(n,dble(freq))* &
            ((2.*real(dble(hCgs)*dble(freq)**3) / (cSpeed*cSpeed)) + jnu) &
            *exp(-hCgs*freq/(kErg*grid%temperature(i1,i2,i3)))
       jnu = 4.*w*hnu(imin+1)
       fac2  = (fourPi/(hcgs*nuArray(imin+1)))*annu(n,dble(nuArray(imin+1)))* &
            ((2.*real(dble(hcgs)*dble(nuArray(imin+1))**3) / (cSpeed*cSpeed)) + jnu) &
            *exp(-hcgs*nuArray(imin+1)/(kErg*grid%temperature(i1,i2,i3)))
       tot = tot + 0.5*(fac1 + fac2)*(nuArray(imin+1)-freq)
       do i = iMin+1, nNu-1
          Jnu = 4.*w*hnu(i)
          fac1  = (fourPi/(hCgs*nuArray(i)))*annu(n,dble(nuArray(i)))* &
               ((2.*real(dble(hCgs)*dble(nuArray(i))**3) / (cSpeed*cSpeed)) + jnu) &
               *exp(-hCgs*nuArray(i)/(kErg*grid%temperature(i1,i2,i3)))
          Jnu = 4.*w*hnu(i+1)
          fac2  = (fourPi/(hcgs*nuArray(i+1)))*annu(n,dble(nuArray(i+1)))* &
               ((2.*real(dble(hcgs)*dble(nuArray(i+1))**3) / (cSpeed*cSpeed)) + jnu) &
               *exp(-hcgs*nuArray(i+1)/(kErg*grid%temperature(i1,i2,i3)))
          tot = tot + 0.5*(fac1 + fac2)*(nuArray(i+1)-nuArray(i))
       enddo

    end if ! (grid%adaptive)
    
    integral1 = tot
    
  end function integral1

  
  real function integral2(n,hnu, nuArray, nNu, rVec, grid, nStar)
    real, dimension(:) :: hnu, nuArray
    integer :: nNu, nStar, n
    type(VECTOR) :: rVec
    type(GRIDTYPE) :: grid
    real :: r
    integer :: i, iMin
    real :: fac1, fac2
    real :: w, x1, freq, tot, jnu

    if (nStar == 1) then
       r = modulus(rVec-grid%starPos1)
       x1 = sqrt(max(0.,(1. - grid%rStar1**2 / r**2)))
       w = 0.5*(1. - x1)
    else
       r = modulus(rVec-grid%starPos2)
       x1 = sqrt(max(0.,1. - grid%rStar2**2 / r**2))
       w = 0.5*(1. - x1)
    endif

    tot = 0.
    freq = ((hydE0eV-eTrans(n))*1.602192e-12)/hcgs
    if ((freq < nuArray(1)).or.(freq > nuArray(nNu))) then
       write(*,*) " "
       write(*,*) "Error in [steteq_mod::integral2]. Freq out of range!", &
            nNu,freq,nuArray(1),nuArray(nNu)
       write(*,*) " =====>  nNu = ", nNu
       write(*,*) " =====> freq = ", freq
       write(*,*) " =====> nuArray(1) =  ", nuArray(1)
       write(*,*) " =====> nuArray(nNu) =  ", nuArray(nNu)
       write(*,*) " =====> You should try to fix this problem !"
       write(*,*) " =====> HINT: Add a point with low flux value at the end of "
       write(*,*) " =====>       your input flux data file."
       write(*,*) " "
       iMin = 1
       jnu = 1.e-28
       stop
    else
       iMin = -1
       call hunt(nuArray, nNu, freq, iMin)
       jnu = 4.*w*logint(freq,nuArray(imin), nuArray(iMin+1), hnu(imin), hnu(imin+1))
    endif
    fac1  = (fourPi/(hcgs*freq))*annu(n,dble(freq))*jnu
    jnu = 4.*w*hnu(imin+1)
    fac2  = (fourPi/(hcgs*nuArray(imin+1)))*annu(n,dble(nuArray(imin+1)))*jnu
    tot = tot + 0.5*(fac1 + fac2)*(nuArray(imin+1)-freq)
    do i = iMin+1, nNu-1
       Jnu = 4.*w*hnu(i)
       fac1  = (fourPi/(hcgs*nuArray(i)))*annu(n,dble(nuArray(i)))*jnu
       Jnu = 4.*w*hnu(i+1)
       fac2  = (fourPi/(hcgs*nuArray(i+1)))*annu(n,dble(nuArray(i+1)))*jnu
       tot = tot + 0.5*(fac1 + fac2)*(nuArray(i+1)-nuArray(i))
    enddo
    integral2 = tot
  end function integral2

  real function integral1_new(n,hnu, nuArray, nNu, nuArray2, hnu2, nNu2, rVec, i1, i2, i3, grid, nStar, &
       thisOctal, thisSubcell, photoOmega, hotOMega)

    real, intent(in), dimension(:) :: hnu, nuArray
    real, intent(in), dimension(:) :: hnu2, nuArray2
    integer, intent(in)      :: nNu, nNu2, nStar, n
    type(VECTOR), intent(in) :: rVec
    type(GRIDTYPE),intent(in):: grid
    integer,intent(in)       :: i1, i2, i3
    type(OCTAL),pointer,optional:: thisOctal
    integer,intent(in),optional :: thisSubcell
    real    :: r
    integer :: i, imin, imin2
    real    :: fac1, fac2
    real    :: w, x1, freq, tot, jnu
    real :: hotOmega, photoOmega
    real :: jnuPhoto, jnuHot
    type(SURFACETYPE) :: starSurface


    if (nStar == 1) then
       r = modulus(rVec-grid%starPos1)
       x1 = sqrt(max(0.,(1. - grid%rStar1**2 / r**2)))
       w = 0.5*(1. - x1)
    else
       r = modulus(rVec-grid%starPos2)
       x1 = sqrt(max(0.,1. - grid%rStar2**2 / r**2))
       w = 0.5*(1. - x1)
    endif

    freq = ((hydE0eV-eTrans(n))*1.602192e-12)/hcgs
    if ((freq < nuArray(1)) .or.(freq > nuArray(nNu))) then
       write(*,*) " "
       write(*,*) "Error in [steteq_mod::integral1_new]. Freq out of range!", &
            nNu,freq,nuArray(1),nuArray(nNu)
       write(*,*) " =====>  nNu = ", nNu
       write(*,*) " =====> freq = ", freq
       write(*,*) " =====> nuArray(1) =  ", nuArray(1)
       write(*,*) " =====> nuArray(nNu) =  ", nuArray(nNu)
       write(*,*) " =====> You should try to fix this problem !"
       write(*,*) " =====> HINT: Add a point with low flux value at the end of "
       write(*,*) " =====>       your input flux data file."
       write(*,*) " "
       jnuPhoto = 1.e-28
       iMin = 1
       stop
    else
       iMin = -1
       call hunt(nuArray, nNu, freq, iMin)
       jnuPhoto = 4.*w*logint(freq,nuArray(iMin), nuArray(iMin+1), hnu(imin), hnu(imin+1))
    endif

    if ((freq < nuArray2(1)) .or.(freq > nuArray2(nNu2))) then
       write(*,*) " "
       write(*,*) "Error in [steteq_mod::integral1_new]. Freq out of range!"
       write(*,*) " =====>  nNu = ", nNu
       write(*,*) " =====> freq = ", freq
       write(*,*) " =====> nuArray(1) =  ", nuArray2(1)
       write(*,*) " =====> nuArray(nNu) =  ", nuArray2(nNu)
       write(*,*) " =====> You should try to fix this problem !"
       write(*,*) " =====> HINT: Add a point with low flux value at the end of "
       write(*,*) " =====>       your input flux data file."
       write(*,*) " "
       jnuHot = 1.e-28
       iMin2 = 1
       stop
    else
       iMin2 = -1
       call hunt(nuArray2, nNu2, freq, iMin2)
       jnuHot = 4.*w*logint(freq,nuArray2(iMin2), nuArray2(iMin2+1), hnu2(imin2), hnu2(imin2+1))
    endif

    jnu = (photoOmega * jnuPhoto + hotOmega * jnuHot) / (photoOmega + hotOmega)

    tot = 0.

    if (grid%adaptive) then 
            
       fac1  = (fourPi/(hCgs*freq))*annu(n,dble(freq))* &
            ((2.*real(dble(hCgs)*dble(freq)**3) / (cSpeed*cSpeed)) + jnu)&
            *exp(-hCgs*freq/(kErg*thisOctal%temperature(thisSubcell)))
       jnu = 4.*w* (photoOmega *hnu(imin+1) + hotOmega * hnu2(imin+1)) / (photoOmega + hotOmega)
       fac2  = (fourPi/(hcgs*nuArray(imin+1)))*annu(n,dble(nuArray(imin+1)))* &
            ((2.*real(dble(hcgs)*dble(nuArray(imin+1))**3) / (cSpeed*cSpeed)) + jnu) &
            *exp(-hcgs*nuArray(imin+1)/(kErg*thisOctal%temperature(thisSubcell)))
       tot = tot + 0.5*(fac1 + fac2)*(nuArray(imin+1)-freq)
       do i = iMin+1, nNu-1
          Jnu = 4.*w* (photoOmega * hnu(i) + hotOmega * hnu2(i)) / (photoOmega + hotOmega)
          fac1  = (fourPi/(hCgs*nuArray(i)))*annu(n,dble(nuArray(i)))* &
               ((2.*real(dble(hCgs)*dble(nuArray(i))**3) / (cSpeed*cSpeed)) + jnu) &
               *exp(-hCgs*nuArray(i)/(kErg*thisOctal%temperature(thisSubcell)))
          Jnu = 4.*w*(photoOmega * hnu(i+1) + hotOmega * hnu2(i+1)) / (photoOmega + hotOmega)
          fac2  = (fourPi/(hcgs*nuArray(i+1)))*annu(n,dble(nuArray(i+1)))* &
               ((2.*real(dble(hcgs)*dble(nuArray(i+1))**3) / (cSpeed*cSpeed)) + jnu) &
               *exp(-hcgs*nuArray(i+1)/(kErg*thisOctal%temperature(thisSubcell)))
          tot = tot + 0.5*(fac1 + fac2)*(nuArray(i+1)-nuArray(i))
       enddo

    else ! not adaptive
            
       fac1  = (fourPi/(hCgs*freq))*annu(n,dble(freq))* &
            ((2.*real(dble(hCgs)*dble(freq)**3) / (cSpeed*cSpeed)) + jnu) &
            *exp(-hCgs*freq/(kErg*grid%temperature(i1,i2,i3)))
       jnu = 4.*w*hnu(imin+1)
       fac2  = (fourPi/(hcgs*nuArray(imin+1)))*annu(n,dble(nuArray(imin+1)))* &
            ((2.*real(dble(hcgs)*dble(nuArray(imin+1))**3) / (cSpeed*cSpeed)) + jnu) &
            *exp(-hcgs*nuArray(imin+1)/(kErg*grid%temperature(i1,i2,i3)))
       tot = tot + 0.5*(fac1 + fac2)*(nuArray(imin+1)-freq)
       do i = iMin+1, nNu-1
          Jnu = 4.*w*hnu(i)
          fac1  = (fourPi/(hCgs*nuArray(i)))*annu(n,dble(nuArray(i)))* &
               ((2.*real(dble(hCgs)*dble(nuArray(i))**3) / (cSpeed*cSpeed)) + jnu) &
               *exp(-hCgs*nuArray(i)/(kErg*grid%temperature(i1,i2,i3)))
          Jnu = 4.*w*hnu(i+1)
          fac2  = (fourPi/(hcgs*nuArray(i+1)))*annu(n,dble(nuArray(i+1)))* &
               ((2.*real(dble(hcgs)*dble(nuArray(i+1))**3) / (cSpeed*cSpeed)) + jnu) &
               *exp(-hcgs*nuArray(i+1)/(kErg*grid%temperature(i1,i2,i3)))
          tot = tot + 0.5*(fac1 + fac2)*(nuArray(i+1)-nuArray(i))
       enddo

    end if ! (grid%adaptive)
    
    integral1_new = tot
    
  end function integral1_new

  
  real function integral2_new(n,hnu, nuArray, nNu, hnu2, nuArray2, nNu2, rVec, grid, nStar, photoOmega, hotOmega)
    real, dimension(:) :: hnu, nuArray
    real, dimension(:) :: hnu2, nuArray2
    integer :: nNu, nNu2, nStar, n
    type(VECTOR) :: rVec
    type(GRIDTYPE) :: grid
    real :: r
    integer :: i, iMin, imin2
    real :: fac1, fac2
    real :: w, x1, freq, tot, jnu
    real :: jnuPhoto, jnuHot
    real :: photoOmega, hotOmega

    tot = 0.

    if (nStar == 1) then
       r = modulus(rVec-grid%starPos1)
       x1 = sqrt(max(0.,(1. - grid%rStar1**2 / r**2)))
       w = 0.5*(1. - x1)
    else
       r = modulus(rVec-grid%starPos2)
       x1 = sqrt(max(0.,1. - grid%rStar2**2 / r**2))
       w = 0.5*(1. - x1)
    endif

    freq = ((hydE0eV-eTrans(n))*1.602192e-12)/hcgs
    if ((freq < nuArray(1)) .or.(freq > nuArray(nNu))) then
       write(*,*) " "
       write(*,*) "Error in [steteq_mod::integral2_new]. Freq out of range!"
       write(*,*) " =====>  nNu = ", nNu
       write(*,*) " =====> freq = ", freq
       write(*,*) " =====> nuArray(1) =  ", nuArray(1)
       write(*,*) " =====> nuArray(nNu) =  ", nuArray(nNu)
       write(*,*) " =====> You should try to fix this problem !"
       write(*,*) " =====> HINT: Add a point with low flux value at the end of "
       write(*,*) " =====>       your input flux data file."
       write(*,*) " "
       jnuPhoto = 1.e-28
       iMin = 1
       stop
    else
       iMin = -1
       call hunt(nuArray, nNu, freq, iMin)
       jnuPhoto = 4.*w*logint(freq,nuArray(iMin), nuArray(iMin+1), hnu(imin), hnu(imin+1))
    endif
    if ((freq < nuArray2(1)) .or.(freq > nuArray2(nNu2))) then
       write(*,*) " "
       write(*,*) "Error in [steteq_mod::integral2_new]. Freq out of range!"
       write(*,*) " =====>  nNu = ", nNu
       write(*,*) " =====> freq = ", freq
       write(*,*) " =====> nuArray(1) =  ", nuArray2(1)
       write(*,*) " =====> nuArray(nNu) =  ", nuArray2(nNu)
       write(*,*) " =====> You should try to fix this problem !"
       write(*,*) " =====> HINT: Add a point with low flux value at the end of "
       write(*,*) " =====>       your input flux data file."	  
       write(*,*) " "
       jnuHot = 1.e-28
       iMin = 1
       stop
    else
       iMin = -1
       call hunt(nuArray2, nNu2, freq, iMin)
       jnuHot = 4.*w*logint(freq,nuArray2(iMin), nuArray2(iMin+1), hnu2(imin), hnu2(imin+1))
    endif

    jnu = (photoOmega * jnuPhoto + hotOmega * jnuHot) / (photoOmega + hotOmega)
    fac1  = (fourPi/(hcgs*freq))*annu(n,dble(freq))*jnu


    jnu = 4.*w* (photoOmega * hnu(imin+1) + hotOmega * hnu2(imin+1)) / (photoOmega + hotOmega)
    fac2  = (fourPi/(hcgs*nuArray(imin+1)))*annu(n,dble(nuArray(imin+1)))*jnu

    tot = tot + 0.5*(fac1 + fac2)*(nuArray(imin+1)-freq)
    do i = iMin+1, nNu-1
       Jnu = 4.*w*(photoOmega *hnu(i) + hotOmega * hnu2(i))/(photoOmega + hotOmega)
       fac1  = (fourPi/(hcgs*nuArray(i)))*annu(n,dble(nuArray(i)))*jnu
       Jnu = 4.*w*(photoOmega *hnu(i+1) + hotOmega * hnu2(i+1))/(photoOmega + hotOmega)
       fac2  = (fourPi/(hcgs*nuArray(i+1)))*annu(n,dble(nuArray(i+1)))*jnu
       tot = tot + 0.5*(fac1 + fac2)*(nuArray(i+1)-nuArray(i))
    enddo
    integral2_new = tot
  end function integral2_new


  subroutine setupMatrices(x, alpha, beta, np, rVec, i1, i2, i3, grid, &
       hnu1, nuArray1, nnu1, hnu2, nuarray2, nnu2, visFrac1, visFrac2, &
       isBinary, thisOctal, thisSubcell,negativeErrors)
    integer, intent(in)                :: np
    type(GRIDTYPE), intent(inout)      :: grid
    logical, intent(in)                :: isBinary
    real, intent(in)                   :: visFrac1, visFrac2
    real(double),intent(in),dimension(:) :: x
    real(double), intent(out) :: alpha(:,:)
    real(double), intent(out) :: beta(:)
    real, intent(in), dimension(:)     :: hnu1, hnu2
    real, intent(in), dimension(:)     :: nuarray1, nuarray2
    integer, intent(in)                :: nNu1, nNu2
    type(VECTOR), intent(in)           :: rVec
    integer, intent(in)                :: i1, i2, i3
    type(octal), pointer, optional     :: thisOctal 
    integer, intent(in),optional       :: thisSubcell 
    integer, intent(inout), optional   :: negativeErrors
    
    real(double)              :: tmp, incr
    real(double), parameter   :: fac = 1.e-3_db
    integer                            :: i, j
    integer, parameter                 :: maxLevels = statEqMaxLevels
    
    if (grid%adaptive) then 

       thisOctal%N(thisSubcell,1:maxLevels) = x(1:maxLevels)
       thisOctal%Ne(thisSubcell) = x(maxLevels+1)
   
       do i = 1, maxLevels
          beta(i) = -equation8(i, maxLevels, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, &
               rVec, i1, i2, i3, grid, visFrac1, visFrac2, isBinary, thisOctal, thisSubcell)
       enddo
       beta(maxLevels+1) = -equation14(maxLevels, grid, i1, i2, i3, thisOctal, thisSubcell)
   
   !    write(*,*) "beta",real(beta(1:7))
       do i = 1, maxLevels
          do j = 1, maxLevels
             tmp = thisOctal%N(thisSubcell,j)
             thisOctal%N(thisSubcell,j) = thisOctal%N(thisSubcell,j) * (1.e0_db + fac)
             incr = equation8(i, maxLevels, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, rVec, &
                  i1, i2, i3, grid, visFrac1, visFrac2, isBinary, thisOctal, thisSubcell)
             alpha(i,j) = (beta(i) + incr) / (tmp*fac)
             if (alpha(i,j) == 0.) then
   !             write(*,*) i,j,beta(i),incr,tmp,grid%n(i1,i2,i3,j),grid%temperature(i1,i2,i3)
             endif
             thisOctal%N(thisSubcell,j) = tmp
          enddo
       enddo
   
   !    do j = 1, maxLevels
   !       tmp = grid%N(i1,i2,i3,j)
   !       grid%N(i1,i2,i3,j) = grid%N(i1,i2,i3,j) * (1.d0 + fac)
   !       incr = equation14(maxLevels, grid, i1, i2, i3)
   !       alpha(maxLevels+1,j) = 1.d0 !(beta(maxLevels+1) + incr) / (tmp*fac) ! 1.d0
   !       grid%N(i1,i2,i3,j) = tmp
   !    enddo
   
       alpha(maxLevels+1,1:maxLevels) = 1.e0_db
   
       tmp = thisOctal%Ne(thisSubcell)
       thisOctal%Ne(thisSubcell) = thisOctal%Ne(thisSubcell) * (1.e0_db+fac)
       do i = 1, maxLevels
          incr =  equation8(i, maxLevels, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, rVec, &
               i1, i2, i3, grid, visFrac1, visFrac2, isBinary, thisOctal, thisSubcell)
   !       alpha(i, maxLevels+1) = (beta(maxLevels+1)+incr)/(tmp*fac)
          alpha(i, maxLevels+1) = (beta(i)+incr)/(tmp*fac)
       enddo
   
       incr = equation14(maxLevels, grid, i1, i2, i3, thisOctal, thisSubcell)
       incr = incr + beta(maxLevels+1)
       alpha(maxLevels+1,maxLevels+1) = incr/(tmp*fac)
       thisOctal%Ne(thisSubcell) = tmp
   
   !    do i = 1, maxLevels+1
   !       do j = 1, maxLevels+1
   !          write(*,*) i,j,real(alpha(i,j))
   !       enddo
   !    enddo
   
    else ! grid not adaptive
           
       do i = 1, maxLevels
          grid%N(i1,i2,i3,i) = x(i)
       enddo
       grid%Ne(i1,i2,i3) = x(maxLevels+1)
   
       do i = 1, maxLevels
          beta(i) = -equation8(i, maxLevels, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, &
               rVec, i1, i2, i3, grid, visFrac1, visFrac2, isBinary)
       enddo
       beta(maxLevels+1) = -equation14(maxLevels, grid, i1, i2, i3)
   
   !    write(*,*) "beta",real(beta(1:7))
       do i = 1, maxLevels
          do j = 1, maxLevels
             tmp = grid%N(i1,i2,i3,j)
             grid%N(i1,i2,i3,j) = grid%N(i1,i2,i3,j) * (1.e0_db + fac)
             incr = equation8(i, maxLevels, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, rVec, &
                  i1, i2, i3, grid, visFrac1, visFrac2, isBinary)
             alpha(i,j) = (beta(i) + incr) / (tmp*fac)
             if (alpha(i,j) == 0.) then
   !             write(*,*) i,j,beta(i),incr,tmp,grid%n(i1,i2,i3,j),grid%temperature(i1,i2,i3)
             endif
             grid%N(i1,i2,i3,j) = tmp
          enddo
       enddo
   
   !    do j = 1, maxLevels
   !       tmp = grid%N(i1,i2,i3,j)
   !       grid%N(i1,i2,i3,j) = grid%N(i1,i2,i3,j) * (1.d0 + fac)
   !       incr = equation14(maxLevels, grid, i1, i2, i3)
   !       alpha(maxLevels+1,j) = 1.d0 !(beta(maxLevels+1) + incr) / (tmp*fac) ! 1.d0
   !       grid%N(i1,i2,i3,j) = tmp
   !    enddo
   
       alpha(maxLevels+1,1:maxLevels) = 1.e0_db
   
       tmp = grid%Ne(i1,i2,i3)
       grid%Ne(i1,i2,i3) = grid%Ne(i1,i2,i3) * (1.e0_db+fac)
       do i = 1, maxLevels
          incr =  equation8(i, maxLevels, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, rVec, &
               i1, i2, i3, grid, visFrac1, visFrac2, isBinary)
   !       alpha(i, maxLevels+1) = (beta(maxLevels+1)+incr)/(tmp*fac)
          alpha(i, maxLevels+1) = (beta(i)+incr)/(tmp*fac)
       enddo
   
       incr = equation14(maxLevels, grid, i1, i2, i3)
       incr = incr + beta(maxLevels+1)
       alpha(maxLevels+1,maxLevels+1) = incr/(tmp*fac)
       grid%Ne(i1,i2,i3) = tmp
   
   !    do i = 1, maxLevels+1
   !       do j = 1, maxLevels+1
   !          write(*,*) i,j,real(alpha(i,j))
   !       enddo
   !    enddo

    end if ! (grid%adaptive)
  
    if (any(x < 0.0_db)) then 
      print *, x
      print *, 'Warning: in setupMatrices, negative value detected'
      if (present(negativeErrors)) then
        negativeErrors = negativeErrors + 1
      else
        print *, 'Negative value handling doesn''t seem to be implemented in'
        print *, '  this version of the code?'
      end if
    end if

  end subroutine setupMatrices

  subroutine mnewt_stateq(grid, ntrial,x,n,tolx,tolf, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, &
                   rVec, i1, i2, i3, visFrac1, visFrac2, isBinary, thisOctal, thisSubcell, &
                   needRestart, maxNegatives)
    implicit none

    type(GRIDTYPE), intent(inout)  :: grid
    integer, parameter             :: maxLevels = statEqMaxLevels
    integer, intent(in)            :: i1, i2, i3
    logical, intent(in)            :: isBinary
    real, intent(in)               :: visFrac1, visFrac2
    type(VECTOR), intent(in)       :: rVec
    real, intent(in), dimension(:) :: nuarray1, nuarray2
    real, intent(in), dimension(:) :: hnu1, hnu2
    integer, intent(in)            :: nNu1, nNu2
    real(double),intent(inout),dimension(:) :: x
    integer, intent(in)            :: nTrial
    integer, intent(in)            :: n
    integer                        :: i, k
    integer                        :: indx(n)
    real(double),intent(in):: tolx, tolf
    real(double)          :: errf, d, errx
    real(double)          :: alpha(n,n),beta(n)
    type(octal), pointer, optional :: thisOctal 
    integer, intent(in),optional   :: thisSubcell 
    logical, parameter             :: debug = .true.
    integer                        :: negativeErrors
    logical, intent(out), optional :: needRestart
    integer, intent(in), optional  :: maxNegatives
    
    if (size(x) /= n) then 
      print *, 'In subroutine mnewt_stateq, we are assuming argumnent ''n'' = SIZE(x),',&
               ' but in this case it is not.'
      stop
    end if
    
    negativeErrors = 0
    needRestart = .false.

      do k = 1, ntrial
        if (debug) print *, 'Trial ',k      
        
        call setupMatrices(x,alpha,beta,n,rVec, i1, i2, i3, grid,Hnu1, nuArray1, nNu1, Hnu2, &
             nuArray2, nNu2, visFrac1, visFrac2, isBinary, thisOctal=thisOctal, &
             thisSubcell=thisSubcell, negativeErrors=negativeErrors)

        if (present(needRestart) .and. present(maxNegatives)) then
          if (negativeErrors >= maxNegatives) then 
            needRestart = .true.
            return
          end if
        end if
             
        errf=0.
        do i = 1, n
          errf=errf+abs(beta(i))
        end do 
        !if (debug) print *, "errf,tolf",errf,tolf
        if ((errf.le.tolf) .and. all(x > 0.0_db)) then 
           if (debug) print *, k, 'trials.'
           return
        end if

        call ludcmp_f77(alpha,n,n,indx,d) 

        call lubksb_f77(alpha,n,n,indx,beta)
        errx=0.
        do i = 1, n
          errx=errx+abs(beta(i))
          x(i)=x(i)+beta(i)
        end do 
        !if (debug) print *, "errx,tolx",errx,tolx
        if ((errx.le.tolx) .and. all(x > 0.0_db)) then
           if (debug) print *, k, 'trials.'
           return
        end if
      end do

      if (any(x < 0.0_db)) then 
        needRestart = .true.
      end if

      print *, ntrial, 'trials' 

  end subroutine mnewt_stateq

  real(double) pure function alpkk(freq,t)
     !
     ! this function returns the free-free absorption coefficient for hydrogen
     !
     real(double), intent(in) :: freq,t
     real(double)             :: wav,gauntf
      
     wav=1.e8_db*cSpeed/freq
     gauntf=giii(1.e0_db,t,wav)
     alpkk=gauntf*real(3.6d8/((dble(freq)**3)*sqrt(dble(t))))
     
  end function alpkk


  real(double) pure function giii (z, t, wl)
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

       giii = sum1 + sum2
  end function giii


    subroutine occultTest(grid, i1, i2, i3, starPos, starRadius, &
         occultPos, occultRadius, visFrac)

      type(GRIDTYPE) :: grid
      integer :: i1, i2, i3
      real :: visFrac
      type(VECTOR) :: toStar, toOccult, direction, rVec, starPos, occultPos
      real :: distToStar, distToOccult, h
      real :: starRadius, occultRadius
      real :: thetaToStar, phiToStar
      logical :: occulted
      real :: dotProd, sinAng
      real :: r, occultedFrac
      integer :: nTheta = 10, nPhi = 10
      real :: cosAng, ang, theta, phi, dtheta, dphi, domega
      integer :: i, j

      rVec = VECTOR(grid%xAxis(i1), grid%yAxis(i2), grid%zAxis(i3))
      toStar = starPos - rVec
      toOccult = occultPos - rVec
      distToStar = modulus(toStar)
      distToOccult = modulus(toOccult)
      toStar = toStar / distToStar
      toOccult = toOccult / distToOccult
      
      if (distToOccult > distToStar) then
         visFrac = 1.
         goto 666
      endif

      
      visFrac = 0.
      occultedFrac = 0.

      h  = sqrt(max(0.,distToStar**2 - starRadius**2))
      cosang = h / distToStar
      ang = acos(min(1.,max(-1.,cosAng)))
      call getPolar(toStar, r, thetaTostar, phiToStar)
      dtheta = pi / real(ntheta-1)
      dphi = twopi / real(nphi-1)
      do i = 1, ntheta
       theta = thetaToStar + (2.*real(i-1)/real(nTheta-1)-1.)*ang
       if (theta > pi) theta = theta - pi
       if (theta < 0.) theta = theta + pi
       dTheta = 2.*ang/real(nTheta-1)
       do j = 1, nphi
          phi = phiToStar + (2.*real(j-1)/real(nPhi-1)-1.)*ang
          if (phi < 0.) phi = phi + twoPi
          if (phi > twoPi) phi = phi + twoPi
          dphi = 2.*ang/real(nPhi-1)

          occulted = .false.
          direction = vector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta))
          domega = sin(theta)*dtheta*dphi
          dotprod = direction .dot. tostar
          if ((dotprod > 0.) .and. (acos(min(1.,max(-1.,dotprod))) < ang)) then

             visFrac = visFrac + dOmega
             dotprod = direction .dot. tooccult
             if (dotprod > 0.) then
                sinang = sqrt(max(0.,(1.-dotprod*dotprod)))
                if ((sinang*disttooccult) < occultradius) then 
                   occultedFrac = occultedFrac + dOmega
                endif
             endif
          endif
       enddo
    enddo
    visFrac = (visFrac-occultedFrac)/visFrac
666 continue
  end subroutine occultTest

  
  subroutine generateOpacities(grid, m, n)

    type(GRIDTYPE) :: grid
    integer, parameter :: maxLevels = statEqMaxLevels
    integer :: m, n
    integer :: i1, i2, i3
    real :: chil, fac
    real :: transe, thresh
    integer :: i,j,k
    real :: chi, eta
    real(double) :: freq

    write(*,'(a)') "Generating opacities..."

    do k = 2, maxlevels
       do i=1, k-1
          lambdaTrans(i, k) = lambdaTrans(k, i)
          ! aEinstein(k, i) = aEinstein(k, i) commented out by nhs
          freq = cspeed / lambdaTrans(i, k)
          bEinstein(k, i) = ((dble(aEinstein(k,i))*dble(cspeed)**2) / (2.d0*dble(hcgs)*dble(freq)**3))
          bEinstein(i, k) = bEinstein(k,i) * gDegen(k) / gDegen(i)
          fStrength(k, i) = -gDegen(i) * fStrength(i,k) / gDegen(k)
       enddo
    enddo

    transe = abs(eTrans(n)-eTrans(m))
    freq = cspeed/lambdaTrans(m,n)

    grid%etaLine = 1.e-20
    grid%etaCont = 1.e-20
    grid%kappaAbs = 1.e-20
    grid%kappaSca = 1.e-20
    grid%chiLine = 1.e-20

    do i1 = 1, grid%na1
       do i2 = 1, grid%na2
          do i3 = 1, grid%na3
             if (.not.grid%inStar(i1,i2,i3).and.grid%inUse(i1,i2,i3)) then
                grid%kappasca(i1,i2,i3,1) = grid%ne(i1,i2,i3) * sigmaE * 1.e10
                !
                ! calculate the line opacity and emissivity
                !
                chil=((pi*eCharge**2)/(mElectron*cSpeed))*fStrength(m,n)
                chil = chil* (grid%n(i1,i2,i3,m)-( ( gDegen(m) / gDegen(n)) * grid%n(i1,i2,i3,n)))
                grid%chiLine(i1,i2,i3) = 1.e10*chil
                
                if (grid%n(i1,i2,i3,n) == 0.d0) then
                   write(*,*) i1,i2,i3,n
                   write(*,*) grid%n(i1,i2,i3,1:maxLevels)
                   write(*,*) grid%Ne(i1,i2,i3)
                   write(*,*) grid%temperature(i1,i2,i3)
                   stop
                endif
                fac=( ((grid%n(i1,i2,i3,m)*gDegen(n))/(grid%n(i1,i2,i3,n)*gDegen(m)))-1.d0)
                grid%etaLine(i1,i2,i3)=1.e10*chil*real((2.d0*dble(hcgs)*dble(freq)**3)/(dble((cSpeed*cSpeed))))/fac

!                write(*,*) i1,i2,i3,grid%etaline(i1,i2,i3)/grid%chiline(i1,i2,i3)
                !
                ! continuous opacity.. bound-free and free-free processes (+es)
                !
                chi=0.d0
                do j=1,maxLevels
                   thresh=(hydE0eVdb-eTrans(j))
                   if (transe.ge.thresh) then
                      chi=chi+(grid%n(i1,i2,i3,j)- &
                           boltzsaha(j, grid%ne(i1,i2,i3), dble(grid%temperature(i1,i2,i3)))* &
                           exp((-hcgs*freq)/(kerg*grid%temperature(i1,i2,i3))))* annu(j,dble(freq))
                   endif
                enddo
                
                chi=chi+real(grid%ne(i1,i2,i3))**2*alpkk(freq,dble(grid%temperature(i1,i2,i3)))*&
                     (1.d0-exp((-hcgs*freq)/(kerg*grid%temperature(i1,i2,i3))))
!                chi=chi+grid%ne(i1,i2,i3)*sigmaE
                
                grid%kappaabs(i1,i2,i3,1) = chi * 1.e10
                !
                ! continuous emissivity...bf and ff
                ! 
                eta=0.d0
                do j=1,20
                   thresh=(hydE0eVdb-eTrans(j))
                   if (transe.ge.thresh) then
                      eta=eta+boltzsaha(j, grid%ne(i1,i2,i3),dble(grid%temperature(i1,i2,i3))) &
                           *annu(j,freq)*exp(-(hcgs*freq)/(kerg*grid%temperature(i1,i2,i3)))
                   endif
                enddo
                
                eta=eta + (grid%ne(i1,i2,i3)**2) * alpkk(freq,dble(grid%temperature(i1,i2,i3)))* &
                     exp(-(hcgs*freq)/(kerg*grid%temperature(i1,i2,i3)))
                
                eta=eta*real((2.0*dble(hcgs)*dble(freq)**3)/(dble(cspeed)**2))
                
                grid%etacont(i1,i2,i3) = eta*1.e10

                if (grid%chiLine(i1,i2,i3) < 0.) then
                   write(*,*) i1, i2, i3, chil
                   grid%chiLine(i1,i2,i3) = 1.e-20
                   grid%etaLine(i1,i2,i3) = 1.e-20
                   grid%etaCont(i1,i2,i3) = 1.e-20
                   grid%kappaAbs(i1,i2,i3,1) = 1.e-20
                   grid%kappaSca(i1,i2,i3,1) = 1.e-20
                endif

             endif
          enddo
       enddo
    enddo
    write(*,'(a)') "Done."


    where (grid%kappaabs < 0.) grid%kappaabs = 1.e-20
  end subroutine generateOpacities

  subroutine amrStateq(grid, contfile, lte, nLower, nUpper, starSurface,&
                       recalcPrevious, ion_name, ion_frac) 
    ! calculate the statistical equilibrium for the subcells in an
    !   adaptive octal grid.

    USE input_variables, ONLY: LyContThick, statEq1stOctant
  
    use parallel_mod
    implicit none
    include 'mpif.h'
    type(GRIDTYPE),intent(inout):: grid      
    character(len=*),intent(in) :: contfile      ! filename for continuum flux
    logical,intent(in)          :: lte           ! true if lte conditions
    integer,intent(in)          :: nLower, nUpper! level populations
    logical, intent(in)         :: recalcPrevious ! whether to improve some previous results
    ! Name of the ion (see opacity_lte_mod.f90 for the list of a valid name.)
    character(LEN=*),intent(in), optional :: ion_name
    real, intent(in), optional            :: ion_frac      ! n_ion/n_specie
    type(SURFACETYPE), intent(in) :: starSurface
    integer                     :: i, m       ! loop counters
    integer, parameter          :: maxLevels = statEqMaxLevels ! number of levels to compute
    integer, dimension(maxLevels), parameter :: allLevels = (/ (i,i=1,maxLevels) /)
    integer                     :: iOctal        ! loop counter
    integer                     :: iSubcell      ! loop counter
    integer                     :: nOctal        ! number of octals in grid
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    real(double),parameter :: tolx = 1.e6_db   ! tolerance
    real(double),parameter :: tolf = 1.e6_db  ! tolerance
    integer                     :: returnVal     ! status value
    real(double),dimension(1:maxLevels+1) :: xall
    real(double)       :: nTot
    real(double)       :: phiT
    real                        :: visFrac1      ! visible fraction of star 1
    real                        :: visFrac2      ! visible fraction of star 2
    logical                     :: isBinary      ! true for binary system
    integer                     :: nNu1, nNu2
    real                        :: hNu1(2000)
    real                        :: hNu2(2000)
    real                        :: nuarray1(2000)
    real                        :: nuarray2(2000)
    type(vector)                :: rVec
    real(double)       :: NeLTE
    type(octal), pointer        :: thisOctal => null()
    real, dimension(maxLevels)  :: departCoeff
    real, dimension(maxLevels)  :: previousXall
    real                        :: previousNeRatio      
    logical, parameter          :: debugInfo = .true.
    type(octalVector)           :: starCentre 
    real(double), parameter :: CI = 2.07d-16   ! in cgs units
    real(double)            :: T,ne, x    
    real(double), parameter :: NeFactor = 20.0_db ! NeFactor
    real :: lam
    integer                     :: numLTEsubcells 

    ! 2-d case variables
    type(octalListElement),pointer  :: listHead => NULL() ! linked list of octals
    type(octalWrapper), allocatable :: planeOctals(:) ! array of pointers to octals in 2d plane
    type(octalWrapper), allocatable :: subsetOctals(:) ! array of pointers to octals

    real(double) :: dummy
    integer       ::   ioctal_beg, ioctal_end  
    logical :: mapped_already
    ! For MPI implementations
    integer       ::   my_rank        ! my processor rank
    integer       ::   np             ! The number of processes
    integer       ::   ierr           ! error flag

    numLTEsubcells = 0

    if (LyContThick) then
      print *, '**************************************************'
      print *, 'Lyman Continuum Thick!'
      print *, '**************************************************'
    end if
    ! Initialize the data arrays (lambdaTrans, bEinstein, fStrength) defined at the top of this module.
    call map_transition_arrays(maxLevels)
    
    starCentre = grid%starPos1

    if (grid%geometry == "ttauri" .or. grid%geometry == "luc_cir3d") then 
      nNu1 = starSurface%nNuHotFlux
      nNu2 = nNu1
      hnu2(1:nNu2) = 0.0
       !do i = 1, nNu2
       !   lam = (cSpeed / nuArray2(i)) * 1.e8
       !   hnu2(i) = pi*blackbody(tAccretion, lam)/fourPi
       !enddo

    else   
      open(20,file=contfile,status="old",form="formatted")
      nnu1 = 1
      do
         read(20,*,iostat=returnVal) nuarray1(nNu1), hnu1(nNu1)
         if (returnVal /= 0) exit
         nNu1 = nNu1 + 1
      end do
      nNu1 = nNu1  - 1
      close(20)
      hNu1(1:nNu1) = hnu1(1:nNu1) / fourPi   ! Converts from flux to flux momemnt
      nNu2 = 0
    endif


    ! 2-D case for rotationally symmetric geometries
        print *, "grid%statEq2d=", grid%statEq2d
        print *, "grid%amr2dOnly=", grid%amr2dOnly
        print *, "statEq1stOctant=", statEq1stOctant
        print *, "lte=", lte
        print *, "recalcPrevious=", recalcPrevious
    mapped_already =.false.
    if (grid%statEq2d .and. (.not. lte) ) then
       ! we first get a list of the octals that lie in a plane across the
       !   space.
       nOCtal = 0
       if (.not. recalcPrevious) then
          call getIntersectedOctals(grid%octreeRoot,listHead,grid,nOctal)
          print *, nOctal, ' octals intersect the plane for 2-d statEq'
       else
          call getIntersectedOctals(grid%octreeRoot,listHead,grid,nOctal,onlyChanged=.true.)
          print *, nOctal, ' octals in the 2-d plane have been updated'
       end if
       
       ! allocate the array of pointers to octals that intersect the plane 
       allocate(planeOctals(nOctal))
       call moveOctalListToArray(listHead,planeOctals)
       
       ! now calculate the non-LTE populations in that plane.
       ! we use 'setupCoeffs' to store the departure coefficients
       ! we use 'propogateCoeffs' to use the previous (radius-sorted) cell's
       !   coefficients as a starting point.
       ! 'firstTime' must be on
       call calcAMRstatEq(planeOctals, setupCoeffs=.true., propogateCoeffs=.true., &
            firstTime=.true.,NeFactor=NeFactor)
       deallocate(planeOctals) ! no longer need the pointers to that plane.
       
       print *, '   Mapping 2-D cells to 3-D...'
       call map2DstatEq(grid%octreeRoot,grid)
       print *, '   ...3-D mapping done'


    !
    !
    elseif (statEq1stOctant) then
       ! we first get a list of the octals that lie in a plane across the
       !   space.
       nOCtal = 0
       if (.not. recalcPrevious) then
          call getOctalsInFirstOctant(grid%octreeRoot,listHead,grid,nOctal)
          print *, nOctal, ' octals (in the first octant) used in amrstatEq'
       else
          call getOctalsInFirstOctant(grid%octreeRoot,listHead,grid,nOctal,onlyChanged=.true.)
          print *, nOctal, ' octals (in the first octant) have been updated/changed'
       end if
       
       ! allocate the array of pointers to octals in first octant
       allocate(subsetOctals(nOctal))
       call moveOctalListToArray(listHead,subsetOctals)
       
       ! now calculate the non-LTE populations 
       ! we use 'setupCoeffs' to store the departure coefficients
       ! we use 'propogateCoeffs' to use the previous (radius-sorted) cell's
       !   coefficients as a starting point.
       ! 'firstTime' must be on
       call calcAMRstatEq(subsetOctals, setupCoeffs=.true., propogateCoeffs=.true., &
            firstTime=.true.,NeFactor=NeFactor)
       deallocate(subsetOctals) ! no longer need the pointers to that plane.
       
       print *, '   Mapping the first octant octal cells to other octants...'
       call mapOctantStatEq(grid%octreeRoot,grid)
       print *, '   ...3-D mapping done'


    ! if we're going to calculate solutions for all the cells (not just those in 
    !   a plane): 
    elseif ((.not. grid%amr2dOnly) .or. lte .or. recalcPrevious) then    


       write(*,*) "Doing all cells..."
      ! get an array of octals comprising the entire tree
      allocate(octalArray(grid%nOctals))
      call getOctalArray(grid%octreeRoot,octalArray, nOctal)

      if ((grid%statEq2d .and. (.not. lte)) .or. recalcPrevious) then

!        if ((grid%statEq2d .and. (.not. lte)) .and. recalcPrevious) then
          ! map the departure coefficients from the 2-d plane to the new cells
          print *, 'Mapping 2-d plane coefficients to modified cells...'
          do iOctal = 1, size(octalArray), 1
            if (any(octalArray(iOctal)%inUse)) then
              call map2dStatEq(octalArray(iOctal)%content,grid,               &
                 subcellMask=(octalArray(iOctal)%content%changed .and. octalArray(iOctal)%inUse))
            end if
          end do
!        end if

          ! we do not want to 'setupCoeffs'
          ! we don't use 'propogateCoeffs' because we have already mapped good
          !   starting values to each cell.
          ! 'firstTime' must be off because we have already computed LTE values.
        call calcAMRstatEq(octalArray, setupCoeffs=.false., propogateCoeffs=.false.,&
                           firstTime=.false.,NeFactor=NeFactor)
        !print *, 'and again, just for fun...'                   
        !call calcAMRstatEq(octalArray, setupCoeffs=.false., propogateCoeffs=.false.,&
        !                   firstTime=.false.,NeFactor=NeFactor)
      else
                
        ! we solve from scratch for all the cells
        call calcAMRstatEq(octalArray, setupCoeffs=.false., propogateCoeffs=.true., &
                           firstTime=.true.,NeFactor=NeFactor)
      end if 
      deallocate(octalArray)
    end if
    ! if we are going to use them later, we must store all the departure
    ! coefficients.
    if (.not. lte) then
      
      if (.not.(allocated(octalArray))) then
        allocate(octalArray(grid%nOctals))
        call getOctalArray(grid%octreeRoot,octalArray, nOctal)
      end if
      
      call saveAllDepartCoeffs(octalArray)
    end if
      
    !! we no longer need to store the departure coefficients
    !if (grid%statEq2d .and. (.not. lte) .and. (.not. recalcPrevious)) &
    !  call removeDepartCoeffs(grid%octreeRoot)

    call generateOpacitiesAMR(grid, nLower, nUpper)
    
  contains
  
  subroutine calcAMRstatEq(octalArray,setupCoeffs,propogateCoeffs,firstTime,NeFactor)
    ! takes an array of pointers to octals, and calculates the statistical 
    !   equilibrium. 

    include 'mpif.h'
    type(octalWrapper), intent(inout), dimension(:) :: octalArray
    logical, intent(in) :: setupCoeffs  ! whether to store the departure 
                                        !   coefficients in the octal structure.
    logical, intent(in) :: propogateCoeffs ! whether to use previous cell's
                                           !   departure coefficients as
                                           !   starting values.
                                           
    logical, intent(in) :: firstTime ! if this is the first run, we have to 
                                     !   calculate LTE values first.
    real(double), intent(in) :: NeFactor ! fudge factor for Ne
    real :: photoOmega, hotOmega
    logical :: needRestart ! non-physical solution flag
    logical, dimension(8) :: canUse ! flag for valid alternative subcells
    integer :: previousSubcell
    real(double),dimension(maxLevels) :: tempLevels
    real(double) :: tempNe
    real,dimension(SIZE(starSurface%nuArray)) :: photoFlux1
     logical :: dcAllocated
     integer, dimension(:), allocatable :: octalsBelongRank
     logical :: rankComplete
     integer :: iRank
     integer :: tag = 0
     integer :: tempInt

    if ( firstTime .and.  (  &
         grid%geometry(1:8) == 'windtest'  .or.   &
         grid%geometry(1:6) == 'ttauri'    .or.   &
         grid%geometry(1:6) == 'luc_cir3d')  )  then
       ! sort the octals by radius
       print *, '   Sorting octalArray by radius...'
       call sortOctalArray(octalArray,grid)
       print *, '   ...sorting complete'
    end if
   
    do i = 1, size(octalArray), 1 
      if (recalcPrevious) then
        ! if we're just updating a few cells, set the mask
        octalArray(i)%inUse = octalArray(i)%content%changed
      end if
      octalArray(i)%content%changed = .false.
    end do
    
    if (firstTime .or. (recalcPrevious .and. .not. grid%statEq2d)) then 
       print *, "   calculating LTE values..." ,size(octalArray)
       do iOctal = 1, SIZE(octalArray), 1
          thisOctal => octalArray(iOctal)%content
          do iSubcell = 1, thisOctal%maxChildren
             if (octalArray(iOctal)%inUse(iSubcell).and.thisOctal%inflow(isubcell)) then 
                call LTElevels(thisOctal%temperature(iSubcell),  &
                               thisOctal%rho(iSubcell),thisOctal%Ne(iSubcell),&
                               thisOctal%nTot(iSubcell),thisOctal%N(iSubcell,:))
                if (.not. associated(thisOctal%departCoeff)) then
                  allocate(thisOctal%departCoeff(8,maxLevels+1))
                  thisOctal%departCoeff = 1.0
                end if

                if (thisOctal%Ne(iSubcell)<0) then
                   write(*,*) 'Error:: ne (electron density) < 0  in stateq_mod::amrStateq!!!'
                   write(*,*) "thisOctal%Ne(",iSubcell,") = ", thisOctal%Ne(iSubcell)
                   write(*,*) "thisOctal%centre = ", thisOctal%centre
                   write(*,*) "iOctal = ", iOctal
                   stop
                end if

             endif
          enddo
       enddo
       print *, "   ...LTE calculations complete." 
    end if

    if (.not.lte) then
       print *, "   calculating non-LTE values..." 
            
       visFrac1 = 1.
       visFrac2 = 0.
       isBinary = .false.
       
       ! default loop indecies
       ioctal_beg = 1
       ioctal_end = SIZE(octalArray)       


!$OMP    ! FOR OpenMP IMPLEMENTATION=======================================================
!$OMP    ! NOT WORKING YET WITH INTEL COMPILER!!!!!!
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(iOctal, thisOctal, iSubcell, NeLTE, xAll, previousXall, rVec) &
!$OMP PRIVATE(photoFlux1, hNu1, nuArray1, previousSubcell, previousNeRatio) &
!$OMP PRIVATE(tempNe, tempLevels, nTot, needRestart, departCoeff, i, starCentre) &
!$OMP SHARED(iOctal_beg, iOctal_end, octalArray, setupCoeffs, starSurface) & 
!$OMP SHARED(recalcPrevious, firstTime, propogateCoeffs, grid) & 
!$OMP SHARED(visFrac1, visFrac2, isBinary) & 
!$OMP SHARED(nNu1, nuArray2, hNu2, nNu2) & 
!$OMP SHARED(numLTEsubcells,  NeFactor) 


      
    ! FOR MPI IMPLEMENTATION=======================================================
    !  Get my process rank # 
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  
    ! Find the total # of precessor being used in this run
    call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
    
    ! we will use an array to store the rank of the process
    !   which will calculate each octal's variables
    allocate(octalsBelongRank(size(octalArray)))
    
    if (my_rank == 0) then
       print *, ' '
       print *, 'calcAMRstatEq routine  computed by ', np-1, ' processors.'
       print *, ' '
       call mpiBlockHandout(np,octalsBelongRank,blockDivFactor=1,tag=tag,&
                            maxBlockSize=10,setDebug=.false.)
    
    endif
    ! ============================================================================

       
       ! initialize some variables
       departCoeff = 1.0
       previousXall = -9.9
       previousNeRatio = 1.0

 if (my_rank /= 0) then
  blockLoop: do     
 call mpiGetBlock(my_rank,iOctal_beg,iOctal_end,rankComplete,tag,setDebug=.false.)
   if (rankComplete) exit blockLoop 

!$OMP DO SCHEDULE(DYNAMIC)
       do iOctal =  iOctal_beg, iOctal_end, 1

          if (debugInfo) print *, 'Octal #',iOctal
          thisOctal => octalArray(iOctal)%content
          if (setupCoeffs) allocate (thisOctal%departCoeff(8,maxLevels+1))
            ! departCoeff(maxLevels+1) is the LTE/nonLTE Ne ratio
          
          do iSubcell = 1,  thisOctal%maxChildren, 1
             if (octalArray(iOctal)%inUse(iSubcell) .and. ((thisOctal%inFlow(iSubcell) .and. &
                thisOctal%nTot(iSubcell) > 1.0) )) then
                if (recalcPrevious .and. .not. firstTime) then
                  ! calculate LTE values, and then scale by old coefficients
                  
                  call LTElevels(thisOctal%temperature(iSubcell),  &
                                 thisOctal%rho(iSubcell),thisOctal%Ne(iSubcell),&
                                 thisOctal%nTot(iSubcell),thisOctal%N(iSubcell,:))

                  thisOctal%Ne(iSubcell) = &
                    thisOctal%Ne(iSubcell) * thisOctal%departCoeff(iSubcell,maxLevels+1)      
                  thisOctal%N(iSubcell,:) = &
                    thisOctal%N(iSubcell,:) * thisOctal%departCoeff(iSubcell,1:maxLevels)      
                  
                end if
                
                ! store the electron density so that we can use it for
                !   comparison later.
                NeLTE = thisOctal%Ne(iSubcell)
                
                ! if we are using starting values from the previous cell, apply
                !   them now
                if (propogateCoeffs) then
                   thisOctal%Ne(iSubcell) = thisOctal%Ne(iSubcell) * previousNeRatio     
                   thisOctal%N(iSubcell,1:maxLevels) = departCoeff(1:maxLevels) * &
                                      boltzSaha(allLevels, thisOctal%Ne(iSubcell),&
                                      real(thisOctal%temperature(iSubcell),kind=db))
                end if
                
                xAll(1:maxLevels) = thisOctal%N(iSubcell,1:maxLevels)
                xAll(maxLevels+1) = thisOctal%Ne(iSubcell)
                
                if (.not. propogateCoeffs .and. .not. firstTime) &
                       previousXall = xAll(1:maxLevels)
                       
                rVec = subcellCentre(thisOctal,iSubcell)

                !call photoSolidAngle(rVec, starSurface, hotOmega, photoOmega)
                photoFlux1 = photoFluxIntegral(rVec, starSurface, SIZE(starSurface%nuArray))
                photoFlux1 = photoFlux1 / fourPi   ! Converts from flux to flux moment
                hNu1(1:nNu1) = photoFlux1(1:nNu1)
                nuArray1(1:nNu1) = starSurface%nuArray(:)

                if (isBinary) then
                  print *, 'Binary geometry not implemented for amrStateq surface models'
                  stop
                end if
                  
                call mnewt_stateq(grid, maxLevels+1, xAll, maxlevels+1, tolx, tolf, hNu1(1:nNu1), &
                           nuArray1, nNu1, &
                           hNu2, nuArray2, nNu2, rVec, 1, 1, 1, visFrac1, visFrac2,&
                           isBinary, thisOctal, iSubcell,   &
                           needRestart=needRestart,maxNegatives=10)

                if (needRestart) then

                  ! first try departure coefficients from a neighbouring cell

                  do previousSubcell = 1, iSubcell-1, 1

                    if (thisOctal%inFlow(previousSubcell)           .and.  &
                        (thisOctal%nTot(previousSubcell) > 1.0)     .and.  &
                        octalArray(iOctal)%inUse(previousSubcell)) then
                      print *, 'Using starting coefficients from previous',&
                               ' subcell: ',previousSubcell
                      ! calculate the departure coefficients for the previous
                      !   subcell 
                      call LTElevels(thisOctal%temperature(previousSubcell),  &
                                     thisOctal%rho(previousSubcell),tempNe,nTot,tempLevels)
                      tempNe = thisOctal%Ne(previousSubcell) / tempNe
                      tempLevels = thisOctal%N(previousSubcell,1:maxLevels) / tempLevels
                     
                      ! recalculate LTE for the current subcell 
                      call LTElevels(thisOctal%temperature(iSubcell),  &
                                     thisOctal%rho(iSubcell),xAll(maxLevels+1),&
                                     thisOctal%nTot(iSubcell),xAll(1:maxLevels))
                      ! apply coefficients
                      xAll(maxLevels+1) = xAll(maxLevels+1) * tempNe
                      xAll(1:maxLevels) = xAll(1:maxLevels) * tempLevels
                      
                      call mnewt_stateq(grid, maxLevels+1, xAll, maxlevels+1, tolx, tolf, &
                           hNu1(1:nNu1), nuArray1, nNu1, &
                           hNu2, nuArray2, nNu2, rVec, 1, 1, 1, visFrac1, visFrac2,&
                           isBinary, thisOctal, iSubcell, needRestart,maxNegatives=8)

                      if (.not. needRestart) exit 
                    end if
                  
                  end do
                
                  if (needRestart) then
                    print *, 'None of the neighbouring subcells provided',&
                             ' suitable starting points.'

                    ! try starting from LTE
                    print *, 'Trying LTE starting point.'
                    call LTElevels(thisOctal%temperature(iSubcell),  &
                                   thisOctal%rho(iSubcell),xAll(maxLevels+1),&
                                   thisOctal%nTot(iSubcell),xAll(1:maxLevels))

                    call mnewt_stateq(grid, maxLevels+1, xAll, maxlevels+1, tolx, &
                         tolf, hNu1(1:nNu1), nuArray1, nNu1, &
                         hNu2, nuArray2, nNu2, rVec, 1, 1, 1, visFrac1, visFrac2,&
                         isBinary, thisOctal, iSubcell, needRestart,maxNegatives=6)

                    if (needRestart) then
                      print *, 'Starting from LTE in stateq did not stop numerical problem...'
                      print *, 'Trying LTE with Ne multiplied by ',NeFactor

                      call LTElevels(thisOctal%temperature(iSubcell),  &
                                     thisOctal%rho(iSubcell),xAll(maxLevels+1),&
                                     thisOctal%nTot(iSubcell),xAll(1:maxLevels),NeFactor=NeFactor)
                      call mnewt_stateq(grid, maxLevels+1, xAll, maxlevels+1, tolx, tolf,  &
                         hNu1(1:nNu1), nuArray1, nNu1, &
                         hNu2, nuArray2, nNu2, rVec, 1, 1, 1, visFrac1, visFrac2,&
                         isBinary, thisOctal, iSubcell, needRestart,maxNegatives=12)
                      
                      if (needRestart) then
                        print *, 'Multiplying Ne by a fudge factor did not stop numerical problem...'
                        print *, 'Fixing subcell (at depth ',thisOctal%nDepth,') at LTE values'
                        call LTElevels(thisOctal%temperature(iSubcell),  &
                                       thisOctal%rho(iSubcell),xAll(maxLevels+1),&
                                       thisOctal%nTot(iSubcell),xAll(1:maxLevels))
                        numLTEsubcells = numLTEsubcells + 1
                      end if
                    end if

                  end if
                end if
                           
                ! store the results
                thisOctal%N(iSubcell,1:maxLevels) = xall(1:maxLevels)
                thisOctal%Ne(iSubcell) = xall(maxLevels+1)
                
                if (propogateCoeffs .or. setupCoeffs .or. associated(thisOctal%departCoeff)) then

                   ! we work out the Ne ratio of non-LTE to LTE...
                   previousNeRatio = thisOctal%Ne(iSubcell) / NeLTE
                  
                   ! and the departure coefficients from LTE at the current (non-LTE) Ne
                   departCoeff(1:maxLevels) = real(xall(1:maxLevels))/                         &
                                                 boltzSaha(allLevels, thisOctal%Ne(iSubcell),&
                                                 real(thisOctal%temperature(iSubcell),kind=db))
                   ! store these in the octal
                   thisOctal%departCoeff(iSubcell,1:maxLevels) = departCoeff(:)
                   thisOctal%departCoeff(iSubcell,maxLevels+1) = previousNeRatio
                end if
                           
                if (debugInfo) then
                   write (*,'(a12,i1,a15,f7.0,a16,e10.1)') '   subcell #', iSubcell, '  temperature: ', & 
                                      thisOctal%temperature(iSubcell),'  num. density: ',thisOctal%rho(iSubcell)/mHydrogen 
                   write (*,'(a,f8.4,a)') '    radius = ',&
                      modulus(subcellCentre(thisOctal,iSubcell) - starCentre)/grid%rStar1,' rStar'
                   write(*,*) '  level   dep.coeff      N          N LTE     start N  '
                   do i = 1 , maxLevels
                      departCoeff(i) = real(xall(i))/boltzSaha(i, thisOctal%Ne(iSubcell),          &
                                                real(thisOctal%temperature(iSubcell),kind=db))
                      write(*,'(a5,i3,1p,e12.3,e12.3,e12.3,e12.3)') '     ',i,departCoeff(i),xall(i),&
                       boltzSaha(i, thisOctal%Ne(iSubcell),real(thisOctal%temperature(iSubcell),kind=db)),&
                       previousXall(i)
                   enddo
                   if (.not. propogateCoeffs .and. .not. firstTime) then 
                     write(*,*) 'Log Ne = ',REAL(log10(thisOctal%Ne(iSubcell))),&
                                '  Ne final / Ne starting value = ',REAL(thisOctal%Ne(iSubcell)/NeLTE)
                   else
                     write(*,*) 'Log Ne = ',REAL(log10(thisOctal%Ne(iSubcell))), &
                                '  Ne non-LTE/LTE = ',REAL(thisOctal%Ne(iSubcell)/NeLTE)
                   end if
                   
                   !write(*,'(a,f6.2)') 'Hot spot viewing fraction(%): ',100.*hotOmega/(photoOmega+hotOmega)
                end if   
                
             else if (setupCoeffs) then
                
                if (.not. associated(thisOctal%departCoeff)) allocate (thisOctal%departCoeff(8,maxLevels+1))
                thisOctal%departCoeff(iSubcell,:) = 1.0
                
             endif
          enddo
       enddo
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL

 end do blockLoop        
 end if ! (my_rank /= 0)

     print *,'Process ',my_rank,' waiting to update values in calcAMRstatEq...' 
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     ! have to send out the 'octalsBelongRank' array
     call MPI_BCAST(octalsBelongRank,SIZE(octalsBelongRank),  &
                    MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

       !
       ! Update the values (N and Ne) of grid computed by all processors.
       !
       do iOctal = 1, SIZE(octalArray)
          !print *,'Process ',my_rank,' starting octal ',iOctal 

          !if (my_rank==0)   print *,'Root reports rank ',octalsBelongRank(ioctal), 'for octal',ioctal
          thisOctal => octalArray(iOctal)%content
          
          ! we may need to allocate departure coefficients
          if (my_rank == octalsBelongRank(iOctal)) then
            dcAllocated = associated(thisOctal%departCoeff)
            !print *, 'rank ',my_rank,'broadcasting dcAllocated:',dcAllocated
          end if
          
          call MPI_BCAST(dcAllocated, 1, &
               MPI_LOGICAL, octalsBelongRank(iOctal), MPI_COMM_WORLD, ierr)
          
          if (dcAllocated .and. .not. associated(thisOctal%departCoeff)) then 
            allocate(thisOctal%departCoeff(8,maxLevels+1))
            !print *, 'process ',my_rank,'allocating dc'
          end if
          
          do iSubcell = 1, thisOctal%maxChildren
             if (octalArray(iOctal)%inUse(iSubcell)) then
                
                call MPI_BCAST(thisOctal%nTot(iSubcell), 1, MPI_DOUBLE_PRECISION,&
                     octalsBelongRank(iOctal), MPI_COMM_WORLD, ierr)
                call MPI_BCAST(thisOctal%Ne(iSubcell), 1, MPI_DOUBLE_PRECISION, &
                     octalsBelongRank(iOctal), MPI_COMM_WORLD, ierr)
                call MPI_BCAST(thisOctal%N(iSubcell, 1:maxlevels), maxlevels, &
                     MPI_DOUBLE_PRECISION, octalsBelongRank(iOctal), MPI_COMM_WORLD, ierr)
             end if
          
             if (dcAllocated) then
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
                call MPI_BCAST(thisOctal%departCoeff(iSubcell, 1:maxlevels+1), maxlevels+1, &
                     MPI_REAL, octalsBelongRank(iOctal), MPI_COMM_WORLD, ierr)
             end if   
          end do
       end do
          
     tempInt = 0
     call MPI_REDUCE(numLTEsubcells,tempInt,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     numLTEsubcells = tempInt
          
     print *,'Process ',my_rank,' finished updating values in calcAMRstatEq...' 
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
      
      if (my_rank==0) then
      print *, "   non-LTE calculations complete..." 
      print *, numLTEsubcells,' subcells were fixed at LTE values'
      end if
       
    endif

    end subroutine calcAMRstatEq

  end subroutine amrStateq

  subroutine saveAllDepartCoeffs(octalArray)
    ! populate the 'departcoeff' variable for all of the octals.
    ! we need these when the 'recalcPrevious' flag is on.
    
    type(octalWrapper), intent(inout), dimension(:) :: octalArray
    
    integer :: iOctal, iSubcell
    real(double) :: NeLTE, dummy
    real(double), dimension(stateqMaxLevels) :: nLte
    type(octal), pointer :: thisOctal
    integer :: i
    integer, dimension(stateqMaxLevels), parameter :: allLevels = (/ (i,i=1,stateqMaxLevels) /)
    
    do iOctal = 1, size(octalArray), 1
      thisOctal => octalArray(iOctal)%content
    
      if (.not. associated(thisoctal%departCoeff)) &
        allocate(thisOctal%departCoeff(8,(stateqMaxLevels+1)))
      
      thisOctal%departCoeff = 1.0  
      
      do iSubcell = 1, thisOctal%maxChildren, 1
      
        if (thisOctal%inFlow(iSubcell) .and. (thisOctal%nTot(iSubcell) > 1.0) ) then
          call LTElevels(thisOctal%temperature(iSubcell),            &
                         thisOctal%rho(iSubcell), NeLTE, dummy, nLTE)
          thisOctal%departCoeff(iSubcell,stateqMaxLevels+1) = &
             thisOctal%Ne(iSubcell) / NeLTE
          thisOctal%departCoeff(iSubcell,1:stateqMaxLevels) = real(thisOctal%N(iSubcell,1:stateqMaxLevels)/  &
             boltzSaha(allLevels, thisOctal%Ne(iSubcell),real(thisOctal%temperature(iSubcell),kind=db)))
        else 
          thisOctal%departCoeff = 1.0
        end if
      end do
    end do

  end subroutine saveAllDepartCoeffs

                
  subroutine generateOpacitiesAMR(grid, nLower, nUpper)
    ! for all the cells in the grid.
  
    type(GRIDTYPE),intent(inout):: grid      
    integer,intent(in)          :: nLower, nUpper! level populations
    real(double)                :: chil
    real(double)                :: chi
    real(double)                :: fac
    real(double)                :: eta
    real(double)                :: thresh
    integer                     :: j             ! loop counter
    integer                     :: iOctal        ! loop counter
    integer                     :: iSubcell      ! loop counter
    integer                     :: nOctal
    real(double)                :: transe
    integer, parameter          :: maxLevels = statEqMaxLevels ! number of levels to compute
    real(double)       :: freq
    logical, parameter          :: debugInfo = .true.
    real(double), parameter :: CI = 2.07d-16   ! in cgs units
    type(octalWrapper), dimension(:), allocatable :: octalArray ! array containing pointers to octals
type(octal), pointer :: testOctal
    call map_transition_arrays(maxLevels)
    freq = cspeed/lambdaTrans(nLower,nUpper)

    nOctal = 0 
    allocate(octalArray(grid%nOctals))
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)
  
    transe = abs(eTrans(nUpper)-eTrans(nLower))
    write(*,'(a,f8.1)') "   generating opacities for ",lambdaTrans(nLower, nUpper)*1.e8

    do iOctal = 1, SIZE(octalArray), 1
       do iSubcell = 1, octalArray(iOctal)%content%maxChildren
          if (octalArray(iOctal)%inUse(iSubcell).and.octalArray(iOctal)%content%inFlow(isubcell)) then
             testOctal => octalArray(iOctal)%content
             octalArray(iOctal)%content%kappaSca(iSubcell,1) = &
                octalArray(iOctal)%content%Ne(iSubcell) * sigmae * 1.e10
             
             if (octalArray(iOctal)%content%kappaSca(iSubcell,1) < 0.0) then              
                octalArray(iOctal)%content%kappaSca(iSubcell,1) = 1.e-20
                print *, ' in amrStatEq, negative kappaAbs value fixed!'
             end if

             ! calculate the line opacity and emissivity
             
             chil=( (pi*eCharge**2) / (mElectron*cSpeed) ) * fStrength(nLower,nUpper)
             chil = chil * (octalArray(iOctal)%content%N(iSubcell,nLower) - &
                       ((gDegen(nLower) / gDegen(nUpper)) * octalArray(iOctal)%content%N(iSubcell,nUpper)))
             octalArray(iOctal)%content%chiLine(iSubcell) = 1.e10 * chil
             
             if (octalArray(iOctal)%content%n(iSubcell,nUpper) == 0.e0_db) then
                write(*,*) 'In amrStatEq, octalArray(iOctal)%content%n(iSubcell,nUpper) == 0.d0'
                write(*,*) nUpper
                write(*,*) octalArray(iOctal)%content%N(iSubcell,1:maxLevels)
                write(*,*) octalArray(iOctal)%content%Ne(iSubcell)
                write(*,*) octalArray(iOctal)%content%temperature(iSubcell)
                stop
             endif
             
             fac=(((octalArray(iOctal)%content%N(iSubcell,nLower)* gDegen(nUpper)) / &
                 (octalArray(iOctal)%content%N(iSubcell,nUpper)*gDegen(nLower)))-1.e0_db)
             octalArray(iOctal)%content%etaLine(iSubcell) = &
                  1.e10*chil*real((2.e0_db*dble(hcgs)*dble(freq)**3)/(dble((cSpeed*cSpeed))))/fac
             
             if (octalArray(iOctal)%content%etaLine(iSubcell) < 0.0) then              
                octalArray(iOctal)%content%etaLine(iSubcell) = 1.e-20
                print *, ' in amrStatEq, negative kappaAbs value fixed!'
             end if
             
             ! continuous opacity.. bound-free and free-free processes (+es)
             
             chi=0.e0_db
             do j=1,maxLevels
                thresh=(hydE0eVdb-eTrans(j))
                if (transe.ge.thresh) then
                   chi=chi+(octalArray(iOctal)%content%N(iSubcell,j)- &
                        boltzsaha(j, real(octalArray(iOctal)%content%Ne(iSubcell),kind=db), &
                                  real(octalArray(iOctal)%content%temperature(iSubcell),kind=db))* &
                        exp((-hcgs*freq)/(kerg*octalArray(iOctal)%content%temperature(iSubcell)))) * &
                        annu(j,real(freq,kind=db))
                endif
             enddo
             
             chi=chi+real(octalArray(iOctal)%content%Ne(iSubcell))**2 * &
                    alpkk(freq,real(octalArray(iOctal)%content%temperature(iSubcell),kind=db))*&
                    (1.e0_db-exp((-hcgs*freq)/(kerg*octalArray(iOctal)%content%temperature(iSubcell))))
             !                chi=chi+grid%ne(i1,i2,i3)*sigmaE
             if (chi < 0.0) then              
                chi = 1.e-30
                print *, ' in amrStatEq, negative kappaAbs value fixed!'
             end if
             octalArray(iOctal)%content%kappaAbs(iSubcell,1) = chi * 1.e10

             
             ! continuous emissivity...bf and ff
              
             eta=0.e0_db
             do j=1,20
                thresh=(hydE0eVdb-eTrans(j))
                if (transe.ge.thresh) then
                   eta=eta+boltzsaha(j, real(octalArray(iOctal)%content%Ne(iSubcell),kind=db),&
                                     real(octalArray(iOctal)%content%temperature(iSubcell),kind=db)) * &
                        annu(j,freq)*exp(-(hcgs*freq)/(kerg*octalArray(iOctal)%content%temperature(iSubcell)))
                endif
             enddo
             
             eta=eta + (octalArray(iOctal)%content%Ne(iSubcell)**2) * &
                        alpkk(freq,real(octalArray(iOctal)%content%temperature(iSubcell),kind=db))* &
                        exp(-(hcgs*freq)/(kerg*octalArray(iOctal)%content%temperature(iSubcell)))
             
             eta=eta*real((2.0*dble(hcgs)*dble(freq)**3)/(dble(cspeed)**2))
             
             octalArray(iOctal)%content%etaCont(iSubcell) = eta*1.e10
             
             if (octalArray(iOctal)%content%etaCont(iSubcell) < 0.0) then              
                octalArray(iOctal)%content%etaCont(iSubcell) = 1.e-20
                print *, ' in amrStatEq, negative kappaAbs value fixed!'
             end if
             
             if (octalArray(iOctal)%content%chiLine(iSubcell) < 0.) then
                write(*,*) iOCtal, chil
                octalArray(iOctal)%content%chiLine(iSubcell)    = 1.e-20
                octalArray(iOctal)%content%etaLine(iSubcell)    = 1.e-20
                octalArray(iOctal)%content%etaCont(iSubcell)    = 1.e-20
                octalArray(iOctal)%content%kappaAbs(iSubcell,1) = 1.e-20
                octalArray(iOctal)%content%kappaSca(iSubcell,1) = 1.e-20
             endif

          endif
       enddo
    enddo

    deallocate(octalArray)

  end subroutine generateOpacitiesAMR
  

  recursive subroutine map2DstatEq(thisOctal,grid,subcellMask)
    ! applies the departure coefficients from the 2-d plane to the rest of the
    !   grid (for rotationally symmetric geometries).

    type(octal), pointer        :: thisOctal
    type(GRIDTYPE),intent(inout):: grid      
    logical(kind=logic), dimension(:), intent(in), optional :: subcellMask

    real(double) :: nTot
    real(double) :: phiT
    integer               :: iSubcell
    type(octalVector)     :: inPlanePoint
    !real(double) :: Ne1, Ne2
    !logical               :: ok
    real(oct)  :: distance
    type(octalVector)     :: starCentre
    type(octalVector)     :: rHat
    type(octalVector)     :: cellCentre
    type(octal),pointer   :: inPlaneOctal
    integer               :: inPlaneSubcell
    real, dimension(1:grid%maxLevels+1) :: departCoeff 
      ! the last value is the electron density coefficent

    integer               :: i
    integer, parameter    :: maxLevels = statEqMaxLevels 
    integer, dimension(maxLevels), parameter :: allLevels = (/ (i,i=1,maxLevels) /)
    type(octal), pointer  :: child
    integer               :: m, iChild
    
    starCentre = grid%starPos1
    rHat = octalVector(1.0_oc/SQRT(2.0_oc),1.0_oc/SQRT(2.0_oc),0.0_oc)
      
    DO iSubcell = 1, thisOctal%maxChildren, 1
      if (present(subcellMask)) then
        if (.not. subcellMask(iSubcell)) cycle
      end if
    
      ! find the corresponding subcell in the statEq-populated plane
      
      cellCentre = subcellCentre(thisOctal,iSubcell)

      ! calculate the radial distance from the star's centre in the x-y plane 
      distance = modulus(cellCentre -                                         &
                         octalVector(starCentre%x,starCentre%y,cellCentre%z))
      
      inPlanePoint = starCentre + (distance * rHat)
      inPlanePoint%z = cellCentre%z 

      ! get the parameters from the in-plane octal
      call amrGridValues(grid%octreeRoot,inPlanePoint,                       &
                         foundOctal=inPlaneOctal,foundSubcell=inPlaneSubcell,&
                         departCoeff=departCoeff)

      ! if the current subcell is in the plane, we can ignore it
      if (associated(thisOctal,inPlaneOctal) .and. (iSubcell==inPlaneSubcell)) &
        cycle 
        
      ! calculate LTE values for the current subcell (to get Ne)
      call LTElevels(thisOctal%temperature(iSubcell),  &
                     thisOctal%rho(iSubcell),thisOctal%Ne(iSubcell),&
                     thisOctal%nTot(iSubcell),thisOctal%N(iSubcell,:))
      thisOctal%Ne(iSubcell) =  thisOctal%Ne(iSubcell) *                        &
                                   real(departCoeff(maxLevels+1),KIND=double)
                                   
      ! calculate level pops. for the current subcell (with modified Ne)
      thisOctal%N(iSubcell,1:maxLevels) = departCoeff(1:maxLevels) * &
                         boltzSaha(allLevels, thisOctal%Ne(iSubcell),&
                         real(thisOctal%temperature(iSubcell),kind=db))
             
      if (.not. associated(thisOctal%departCoeff)) &
        allocate(thisOctal%departCoeff(8,maxLevels+1))

      thisOctal%departCoeff(iSubcell,:) = departCoeff
      thisOctal%changed(iSubcell) = .false.
    
    end do 
    
    ! call recursively on any children
    if (thisOctal%nChildren > 0) then 
      do iChild = 1, thisOctal%nChildren, 1
        child => thisOctal%child(iChild)
        call map2DstatEq(child,grid)
      end do
    end if
  
  end subroutine map2dStatEq
  

  ! This is base on map2dStatEq. 
  ! Maps the departure coefficients to the first octant (x>0, y>0, z>0) 
  ! to other 7 octants using symmetries ( symmetry about z=0 plane, 
  ! phi = pi/2 rotational symmetry around z axis)
  ! 
 recursive subroutine mapOctantStatEq(thisOctal,grid,subcellMask)
    ! applies the departure coefficients from the 2-d plane to the rest of the
    !   grid (for rotationally symmetric geometries).

    type(octal), pointer        :: thisOctal
    type(GRIDTYPE),intent(inout):: grid      
    logical(kind=logic), dimension(:), intent(in), optional :: subcellMask

    real(double) :: nTot
    real(double) :: phiT
    integer               :: iSubcell
    !real(double) :: Ne1, Ne2
    !logical               :: ok
    real(oct)  :: distance
    type(octalVector)     :: starCentre
    type(octalVector)     :: cellCentre, r, r2
    type(octal),pointer   :: foundOctal
    integer               :: foundSubcell
    real(oct) :: phi
    real, dimension(1:grid%maxLevels+1) :: departCoeff 
      ! the last value is the electron density coefficent

    integer               :: i
    integer, parameter    :: maxLevels = statEqMaxLevels 
    integer, dimension(maxLevels), parameter :: allLevels = (/ (i,i=1,maxLevels) /)
    type(octal), pointer  :: child
    integer               :: m, iChild
    logical   :: skip_rotation
    
    starCentre = grid%starPos1
      
    MAINLOOP: DO iSubcell = 1, thisOctal%maxChildren, 1
      if (present(subcellMask)) then
        if (.not. subcellMask(iSubcell)) cycle
      end if
    
      ! find the corresponding subcell in the statEq-populated plane      
      cellCentre = subcellCentre(thisOctal,iSubcell) - starCentre
      ! 
      if (cellCentre%x > 0.0 .and. cellCentre%y > 0.0 .and. cellCentre%z > 0.0) then
         ! the point is in the first octant; hence, no need to map the value.
         CYCLE  MAINLOOP
      end if

      ! Now transfer this vector to a corresponding vector in the first
      ! octant.
      r = cellCentre
      if (r%z < 0.0) r%z = -(r%z)  ! flips sign

      ! then rotate around the z axis.
      skip_rotation=.false.
      if (r%x<0.0) then
         if (r%y<0.0) then
            ! in third quadrant
            phi = pi
         else
            ! in the second quadrant
            phi = 0.5d0*pi
         end if         
      else
         if (r%y<0.0) then
            ! in fourth quadrant
            phi = 1.5d0*pi
         else
            ! in the first quadrant
            phi = 0.0d0
            skip_rotation=.true.
         end if                  
      end if
      if (skip_rotation)   then
         r2 = r
      else
         r2 = rotateZ(r, phi)   
      end if

      ! get the parameters from the octals in the first octant
      call amrGridValues(grid%octreeRoot, r2,                       &
                         foundOctal=foundOctal,foundSubcell=foundSubcell,&
                         departCoeff=departCoeff)
        
      ! calculate LTE values for the current subcell (to get Ne)
      call LTElevels(thisOctal%temperature(iSubcell),  &
                     thisOctal%rho(iSubcell),thisOctal%Ne(iSubcell),&
                     thisOctal%nTot(iSubcell),thisOctal%N(iSubcell,:))
      thisOctal%Ne(iSubcell) =  thisOctal%Ne(iSubcell) *                        &
                                   real(departCoeff(maxLevels+1),KIND=double)
                                   
      ! calculate level pops. for the current subcell (with modified Ne)
      thisOctal%N(iSubcell,1:maxLevels) = departCoeff(1:maxLevels) * &
                         boltzSaha(allLevels, thisOctal%Ne(iSubcell),&
                         real(thisOctal%temperature(iSubcell),kind=db))
             
      if (.not. associated(thisOctal%departCoeff)) &
        allocate(thisOctal%departCoeff(8,maxLevels+1))

      thisOctal%departCoeff(iSubcell,:) = departCoeff
      thisOctal%changed(iSubcell) = .false.
    
   end do MAINLOOP
    
    ! call recursively on any children
    if (thisOctal%nChildren > 0) then 
      do iChild = 1, thisOctal%nChildren, 1
        child => thisOctal%child(iChild)
        call mapOctantStatEq(child,grid)
      end do
    end if
  
  end subroutine mapOctantStatEq
  



  ! Initialize the data arrays (lambdaTrans, bEinstein, fStrength) defined at the top of this module.
  subroutine map_transition_arrays(maxLevels)
    implicit none
    integer, intent(in)   :: maxLevels      
    integer               :: i, k
    real(double) :: freq    
    logical, save         :: alreadyDone = .false.
    
    if (alreadyDone) return
    
    if (maxLevels > 20) then
      print *, 'Number of levels is > 20. Stopping in stateq_mod'
      stop
    end if
    do k = 2, maxlevels
       do i=1, k-1
          lambdaTrans(i, k) = lambdaTrans(k, i)
          ! aEinstein(k, i) = aEinstein(k, i) commented out by nhs
          freq = cspeed / lambdaTrans(i, k)
          bEinstein(k, i) = ((aEinstein(k,i)*(dble(cspeed))**2)  &
                             / (2.e0_db*dble(hcgs)*(dble(freq))**3))
          bEinstein(i, k) = bEinstein(k,i) * gDegen(k) / gDegen(i)
          fStrength(k, i) = -gDegen(i) * fStrength(i,k) / gDegen(k)
       end do
    end do

    alreadyDone = .true.

  end subroutine map_transition_arrays

  
  !
  ! Computes the partition function of hydrogen (HI)
  ! 
  
  pure function Z_HI(nmax, T) RESULT(out)
    implicit none
    real(double)  :: out 
    integer, intent(in) :: nmax       ! maximum level used for the hydrogen
    real(double), intent(in) :: T ! Temperature in Kelvin
    !    
    integer :: i 
    real(double) :: sum
    real(double), parameter :: k = 8.617342d-5      ! [eV/K]      ! Boltzmann constant
    real(double)  :: x
    real(double)  :: exp_minus_x

    sum = 0.0d0
    do i = 1, nmax
       x = hydE0eV*(1.0d0-1.0d0/dble(i*i))/(k*T)
       if (x < 0.001) then ! use approximation
          exp_minus_x = 1.0d0 - x + 0.5d0*x*x - x*x*x/3.0d0
       else
          exp_minus_x = EXP(-x)
       end if
       sum = sum + 2.0d0*dble(i*i) * exp_minus_x
    end do
    
    out = sum

  end function Z_HI

  pure subroutine LTElevels(Tsingle,rho,Ne,nTot,levels,NeFactor)

    real, intent(in)                  :: Tsingle 
    real(double)             :: T
    real, intent(in)                  :: rho
    real(double), intent(out):: Ne 
    real(double), intent(out):: nTot
    real(double), dimension(:), intent(out) :: levels
    real(double), intent(in), optional :: NeFactor 
    real(double), parameter  :: CI = 2.07d-16   ! in cgs units
    integer, parameter                :: maxLevels = statEqMaxLevels 
    integer :: i
    integer, dimension(maxLevels), parameter :: allLevels = (/ (i,i=1,maxLevels) /)
    real(double)             :: phiT
    
    !==========================================================================
    !
    ! Note:           1         h^2
    !         CI =   --- ------------------  = 2.07x10^-16 (in cgs)
    !                 2     2 Pi m k T
    !
    ! where m is the mass of electron.
    ! See page 49-50 of "Fundation of Radiation Hydrodynamic" by Mihalas & Mihalas.

    T = Tsingle
    
    nTot = rho/mHydrogen
    phiT = 0.5d0*CI*Z_HI(maxLevels,T)*EXP(real(hydE0eV,kind=double)/(kev*T))

    ! Solving for phi(T)*ne^2 + 2ne -n =0 for ne, and choosing the physical
    ! solution ... 
    Ne = (sqrt(nTot*phiT+1.0_db) -1.0_db)/phiT
    Ne = min(Ne, nTot)  ! to avoid unphysical solution.
    if (Ne<=0) Ne =nTot ! to avoid unphysical solution.
   
    if (present(NeFactor)) Ne = Ne * NeFactor
    
    levels = boltzsaha(allLevels,Ne,T) 
  end subroutine LTElevels

end module stateq_mod

!!! vim:set filetype=fortran :                                !!!
!!! otherwise vim won't recognize a file with the suffix .raw !!!

