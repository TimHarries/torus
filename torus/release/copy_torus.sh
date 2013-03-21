#!/bin/ksh

cd ..
svnstring=`svnversion`
export TORUS_TAR_FILE=torus2.0-${svnstring}_public_alpha.tar
echo "Making ${TORUS_TAR_FILE}"
cd release

echo "Copying Torus source code"
cp ../*.[fF]90 . 

echo "Removing files not for release"
# Get rid of files we don't need
rm hydrodynamics_mod.?90 photoionAMR_mod.?90 photoion_mod.?90
rm photoion_utils_mod.?90 angularImage_mod.?90 molecular_mod.?90 
rm torusMod.?90 ion_mod.?90 nbody_mod.?90 qShep*90 timedep_mod.?90
rm cluster_class.?90 sph_data_class.?90 
rm phfit2.?90 cmf_mod.?90 modelatom_mod.?90 h21cm_mod.?90
rm isochrone_class.?90 viscosity_mod.?90 stateq_mod.?90

# Only used by stateq_mod
rm math_mod2.?90

cp ../makedepend . 

echo 
echo -n "Number of lines in reduced code: "
num_red=`wc -l *90 | tail -1 | awk '{print $1}'`
echo $num_red

echo -n "Number of lines in full code: "
num_full=`wc -l ../*90 | tail -1 | awk '{print $1}`
echo $num_full

frac=`echo $num_red $num_full | awk '{print $1/$2}`
echo "Fraction in reduced code: ${frac}"
echo 

echo "Making tar file"
make depends > /dev/null 2>&1
mkdir torus
mv *90 make.depends makedepend torus
cp svn_version.h Makefile torus
tar cf ${TORUS_TAR_FILE} torus

echo "Testing build"
rm -rf build
mkdir build 
cd build
ln -s ../${TORUS_TAR_FILE}
tar xf ${TORUS_TAR_FILE}
mv torus/* .
make SYSTEM=ifort > compile_log 2>&1

if [[ -x torus.ifort ]]; then 
    echo "Found torus.ifort: build OK"
    echo "Cleaning up"
    cd ..
    rm -r torus 
    rm -r build
    gzip ${TORUS_TAR_FILE}
    exit 0
else
    echo "Found torus.ifort not found: build failed"
    exit 1
fi
