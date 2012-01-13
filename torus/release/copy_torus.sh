#!/bin/ksh

cp ../*.[fF]90 . 

# Get rid of files we don't need
rm hydrodynamics_mod.?90 photoionAMR_mod.?90 photoion_mod.?90
rm photoion_utils_mod.?90 angularImage_mod.?90 molecular_mod.?90 
rm torusMod.?90 ion_mod.?90 nbody_mod.?90 qShep*90 timedep_mod.?90
rm cluster_class.?90 sph_data_class.?90 wr104_mod.?90
rm phfit2.?90 cmf_mod.?90 modelatom_mod.?90 h21cm_mod.?90
rm isochrone_class.?90 stateq_mod.?90 math_mod2.?90 hyd_col_coeff.?90

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

