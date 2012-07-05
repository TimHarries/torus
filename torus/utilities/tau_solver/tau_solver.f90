program tau_solver

implicit none


double precision, parameter :: T_13 = 1.2169526815414429
double precision, parameter :: T_18 = 0.8931993842124939
double precision, parameter :: X = 8.d0
double precision, parameter :: dtau = 1.d-3
double precision :: ratio
double precision :: thisTau, startTau, endTau
double precision :: testVal, iniTest, endTest, lowestDiff, lowestTau
integer :: i, nsteps

ratio = T_13/T_18

startTau = 0.d0
endTau = 100.d0

iniTest = testRat(startTau)
endTest = testRat(endTau)

print *, "----------------------"
print *, "TAU SOLVER"
print *, "----------------------"
print *, "- Thomas Haworth 2012"
print *, " " 
print *, " " 
print *, "searching for ratio = ", ratio
print *, "between ", startTau, "and ", endTau
print *, " "

!print *, "diff 1 = ", ratio - iniTest
!print *, "diff 2 = ", ratio - endTest

nsteps = int((endTau-startTau)/dTau)

   thisTau = startTau
   lowestDiff = 1.d10
do i = 1, nSteps
   testVal = abs(ratio -  testRat(thisTau))
   if(testVal < lowestDiff) then
      lowestDiff = testVal
      lowestTau = thisTau
   end if
   thisTau = thisTau + dTau
end do

print *, "lowest tau is: ", lowestTau
print *, "lowest diff is: ", lowestDiff

contains

double precision function testRat(thisTau)
  double precision, parameter :: X = 16.d0
  double precision :: thisTau

  testRat = (1.d0 - exp(-thisTau))/(1.d0 - exp(-(thisTau/X)))
!  print *, "top ", dble((1.d0 - exp(-(thisTau/X))))
!  print *, "bottom ", dble((1.d0 - exp(-thisTau)))
!  print *, "exp(-thisTau)", exp(-thisTau)
!  print *, "exp(-thisTau/X)", exp(-thisTau/X)
end function testRat




end program tau_solver
