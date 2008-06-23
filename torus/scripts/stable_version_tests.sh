#!/bin/ksh

# Run tests which need to be passed by Torus stable versions.
# Author: D. Acreman, June 2008

build_torus()
{
mkdir ${SYSTEM}_${TORUS_RUN_TYPE}
cd ${SYSTEM}_${TORUS_RUN_TYPE}

ln -s ../torus/* .
make depends > compile_log 2>&1

if [[ ${TORUS_RUN_TYPE} == debug ]]; then
    make debug=yes >> compile_log 2>&1 
else
    make >> compile_log 2>&1
fi

cd ..
}

run_benchmark()
{
mkdir run_${SYSTEM}_${TORUS_RUN_TYPE}
cd run_${SYSTEM}_${TORUS_RUN_TYPE}
cp -r ../torus/benchmarks/disc/* .
ln -s ../${SYSTEM}_${TORUS_RUN_TYPE}/torus.${SYSTEM} .

if [[ ${SYSTEM} == intelmac ]]; then
    ./torus.${SYSTEM} > log 2>&1 
elif [[ ${SYSTEM} == ompi ]]; then
    mpirun -np 4 torus.${SYSTEM} > log 2>&1
else
    echo "Unrecognised system ${SYSTEM}"
    exit 1
fi
}

check_benchmark()

{
echo Compiling comparespec code
g95 -o comparespec comparespec.f90
cp test_inc013.dat speca.dat
cp sed100_125.dat specb.dat
echo Comparing the 12.5 degree model...
./comparespec
cp test_inc077.dat speca.dat
cp sed100_775.dat specb.dat
echo Comparing the 77.5 degree model...
./comparespec
cd ..
}

#########################################################################################

# Top level working directory
test_directory=torus_stable_version_tests

# G95 environment variables
export G95_FPU_INVALID=true
export G95_FPU_ZERODIV=true
export G95_FPU_OVERFLOW=true

# Set up directory to use for tests
# Run from current working directory
if [[ -e ${test_directory} ]]; then
    echo "A directory called ${test_directory} already exists"
    echo "Not overwriting existing directory. Aborting ..."
    exit 1
else
    mkdir ${test_directory}
fi

# Get latest code from CVS repository
cd ${test_directory}
cvs -d pinky:/h/th/CVS co torus

for torus_type in fast debug; do
    for sys in ompi intelmac; do

	export SYSTEM=${sys}
	export TORUS_RUN_TYPE=${torus_type}
	build_torus
	run_benchmark
	check_benchmark

    done
done

exit

# End of file ######################################################

