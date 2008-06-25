#!/bin/ksh

# Purpose: run the Torus daily benchmark suite
# Author: D. Acreman

build_and_prepare()

{
this_system=$1

log_file=compile_log_${this_system}.txt
echo "Building executable for ${this_system}"
cd build_${this_system}
ln -s ../torus/* .
export SYSTEM=${this_system}
echo Making dependancies...
/usr/bin/make depends > ${log_file} 2>&1
echo Making torus...
/usr/bin/make debug=yes >> ${log_file} 2>&1
if [[ $? -eq 0 ]]; then
# Count number of warnings. Subtract 2 because there are always warnings
# about include files from the make depends step (run twice).
    num_warn=`grep -i warning ${log_file} | wc -l | awk '{print $1 - 2}`
    echo "Compilation completed with ${num_warn} warnings."
else
    echo "Compilation failed."
fi

# Prepare disc benchmark run directory
cd ../run_${this_system}
ln -s ../build_${this_system}/torus.${this_system}
cp -r ../torus/benchmarks/disc/* .

# Prepare molebench run directory
cd ../run_${this_system}_molebench
mkdir J
ln -s ../build_${this_system}/torus.${this_system}
cp -r ../torus/benchmarks/molebench/* .
cp ../torus/data/hco_benchmark.mol .

cd ..

}


check_benchmark()

{
echo Compiling comparespec code
/sw/bin/g95 -o comparespec comparespec.f90
cp test_inc013.dat speca.dat
cp sed100_125.dat specb.dat
echo Comparing the 12.5 degree model...
./comparespec
cp test_inc077.dat speca.dat
cp sed100_775.dat specb.dat
echo Comparing the 77.5 degree model...
./comparespec
}

check_molebench()
{
echo Compiling compare_molbench code
g95 -o compare_molbench compare_molbench.f90
model_file=`ls results.* | tail -1`
ln -s ${model_file} results.dat
./compare_molbench
}


# Main part of script starts here ------------------------------------------

# Set up 
test_dir=${HOME}/torus_daily_test
sys_to_test="ompi"
export CVSROOT=${USER}@pinky.astro.ex.ac.uk:/h/th/CVS
export CVS_RSH=ssh
export PATH=/sw/bin:/usr/local/bin:${PATH}

# If all the required files are copied in to the working directory
# then this line is not required. 
export TORUS_DATA=${test_dir}/torus/data

export G95_FPU_INVALID=true
export G95_FPU_ZERODIV=true
export G95_FPU_OVERFLOW=true


echo TORUS test harness script started
echo ---------------------------------
echo

if [[ -e ${test_dir} ]]; then
    echo "Found ${test_dir}"
else
    echo "Making ${test_dir}"
    mkdir ${test_dir}
fi

cd ${test_dir}

echo Checking out torus from CVS archive...
rm -rf torus
/usr/bin/cvs -q co torus > cvs_log.txt

for sys in ${sys_to_test}; do
    rm -rf build_${sys} run_${sys} run_${sys}_molebench
    mkdir  build_${sys} run_${sys} run_${sys}_molebench
    build_and_prepare ${sys}
done

echo "g95 environment variables are:"
printenv | grep -i g95

cd run_ompi
echo "Running torus.ompi"
/usr/local/bin/mpirun -np 4 torus.ompi > run_log_ompi.txt 2>&1
check_benchmark
cd ..

cd run_ompi_molebench
echo "Running torus.ompi molebench"
/usr/local/bin/mpirun -np 4 torus.ompi > run_log_ompi.txt 2>&1
check_molebench
cd ..

exit

