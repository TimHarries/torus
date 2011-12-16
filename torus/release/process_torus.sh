#!/usr/bin/ksh

echo "Processing Torus source code"
rm -f *90
ln -s ../*.F90 . 
for srcFile in *.F90; do
    outfile=`echo $srcFile | sed s/F90/f90/`
    fpp -P $srcFile $outfile
# Add OpenMP key if you want to build with OpenMP. 
#    fpp -P -D_OPENMP $srcFile $outfile
done

rm *.F90
cp ../*.f90 .

# Get rid of files we don't need
rm hydrodynamics_mod.f90 photoionAMR_mod.f90 photoion_mod.f90
rm photoion_utils_mod.f90 angularImage_mod.f90 molecular_mod.f90 
rm torusMod.f90 ion_mod.f90 nbody_mod.f90

ln -s ../makedepend

echo -n "Number of lines "
wc -l *90 | tail -1 

echo "All done"
