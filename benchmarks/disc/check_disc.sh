#!/bin/bash

# Check SEDs
echo "Compiling comparespec code"
gfortran -o comparespec comparespec.f90
# Check first SED
echo Comparing the 12.5 degree model...
cp test_inc013.dat speca.dat
cp sed100_125.dat specb.dat
./comparespec
# Check second SED
echo Comparing the 77.5 degree model...
cp test_inc077.dat speca.dat
cp sed100_775.dat specb.dat
./comparespec

echo "Compiling check_disc_image"
gfortran -o check_disc_image check_disc_image.f90 -lcfitsio
if [[ -x check_disc_image ]]; then
    echo "Checking image"
else
    echo "Could not build check_disc_image. Skipping this test ..."
fi

