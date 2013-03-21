#!/usr/bin/ksh

echo "Processing Torus source code"
rm -f *90
ln -s ../*.F90 . 
for srcFile in *.F90; do
    outfile=`echo $srcFile | sed s/F90/f90/`
#    fpp -P $srcFile $outfile
# Add OpenMP key if you want to build with OpenMP. 
    fpp -P -D_OPENMP -DUSEZLIB -DMEMCHECK $srcFile $outfile
done

rm *.F90
cp ../*.f90 .

# Get rid of files we don't need
rm hydrodynamics_mod.f90 photoionAMR_mod.f90 photoion_mod.f90
rm photoion_utils_mod.f90 angularImage_mod.f90 molecular_mod.f90 
rm torusMod.f90 ion_mod.f90 nbody_mod.f90 qShep*90 timedep_mod.f90
rm cluster_class.f90 sph_data_class.f90 
rm phfit2.f90 cmf_mod.f90 modelatom_mod.f90 h21cm_mod.f90
rm isochrone_class.f90 viscosity_mod.f90 stateq_mod.f90 mpi_amr_mod.f90

# Only used by stateq_mod
rm math_mod2.f90

rm -f makedepend
ln -s ../makedepend

echo -n "Number of lines in reduced code: "
num_red=`wc -l *90 | tail -1 | awk '{print $1}'`
echo $num_red

echo -n "Number of lines in full code: "
num_full=`wc -l ../*90 | tail -1 | awk '{print $1}`
echo $num_full

frac=`echo $num_red $num_full | awk '{print $1/$2}`
echo "Fraction in reduced code: ${frac}"


echo "All done"
