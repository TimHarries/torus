#!/bin/ksh

echo "Extracting release version of Torus"

which ifort > /dev/null 2>&1
if [[ $? -eq 0 ]]; then
    echo "Will use ifort for build test"
    export TESTSYS=ifort
else
    echo "Will use gfortran for build test"
    export TESTSYS=gfortran
fi 

# Get the subversion revision string from the top level of the repository
cd ..
svnstring=`svnversion`
export TORUS_TAR_FILE=torus2.0-${svnstring}_public.tar
echo "Making ${TORUS_TAR_FILE}"
cd release

echo "Copying Torus source code"
cp ../*.[fF]90 . 

# Get rid of files we don't want to release
echo "Removing files not for release"
rm hydrodynamics_mod.?90 photoionAMR_mod.?90 photoion_mod.?90
rm photoion_utils_mod.?90 angularImage_mod.?90 molecular_mod.?90 
rm torusMod.?90 ion_mod.?90 nbody_mod.?90 qShep*90 timedep_mod.?90
rm datacube_mod.?90 phfit2.?90 cmf_mod.?90 modelatom_mod.?90
rm viscosity_mod.?90 stateq_mod.?90
rm math_mod2.?90

# Fix permissions
chmod -x *90

# Report some stats
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
# Use the main repository makefile to generate make.depends
ln -s ../Makefile
ln -s ../makedepend 
make depends > /dev/null 2>&1
rm Makefile makedepend

# Now generate the release version Makefile
cat Makefile_release make.depends > Makefile

mkdir torus
mv *90 Makefile torus
cp svn_version.h  torus
tar cf ${TORUS_TAR_FILE} torus

echo "Testing build"
rm -rf build
mkdir build 
cd build
ln -s ../${TORUS_TAR_FILE}
tar xf ${TORUS_TAR_FILE}
mv torus/* .
make SYSTEM=${TESTSYS} > compile_log 2>&1

if [[ -x torus.${TESTSYS} ]]; then 
    echo "Found torus.${TESTSYS}: build OK"
    echo "Cleaning up"
    cd ..
    rm -r torus 
    rm -r build
    rm make.depends
    gzip ${TORUS_TAR_FILE}
    exit 0
else
    echo "torus.${TESTSYS} not found: build failed"
    exit 1
fi
