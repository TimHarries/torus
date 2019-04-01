echo "Compiling check.f90"
gfortran -o check_nbody check.f90
echo "Running check"
./check_nbody

