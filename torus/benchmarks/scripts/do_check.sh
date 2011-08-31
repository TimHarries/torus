#!/bin/ksh

if [[ -x comparespec ]]; then
    echo "Found comparespec"
else
    echo "Building comparespec"
    gfortran -o comparespec comparespec.f90
fi

echo Comparing the 12.5 degree model...
cp test_inc013.dat speca.dat
cp sed100_125.dat specb.dat
comparespec

echo Comparing the 77.5 degree model...
cp test_inc077.dat speca.dat
cp sed100_775.dat specb.dat
comparespec

