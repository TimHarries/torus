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

# Prepare cylindrical polar version of disc benchmark
cd ../run_${sys}_disc_cylindrical
ln -s ../build_${this_system}/torus.${this_system}
cp -r ../torus/benchmarks/disc_cylindrical/parameters.dat .
cp -r ../torus/benchmarks/disc/sed* .
cp -r ../torus/benchmarks/disc/compare* .

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

sphbench()
{
rm -rf sphbench
mkdir sphbench
cd sphbench

echo
echo "Making Torus library"
mkdir build
cd build 
ln -s ../../torus/* .
make depends > compile_log 2>&1
make debug=yes lib >> compile_log 2>&1 

echo "Compiling sphbench"
cp ../../torus/benchmarks/sphbench/*.f90 .
cp ../../torus/benchmarks/sphbench/compile .
./compile >> compile_log 2>&1 
cd ..

mkdir run
cd run
cp ../../torus/benchmarks/sphbench/*.dat . 
cp ../../torus/benchmarks/sphbench/*.txt . 
cp ../../torus/benchmarks/disc/comparespec.f90 . 
cp ../../torus/benchmarks/disc/sed* .
ln -s ../../torus/isochrones/iso* .
ln -s ../build/sphbench .
echo "Running sphbench"
/usr/local/bin/mpirun -np 4 sphbench > run_log 2>&1
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
/usr/bin/cvs -q co torus > cvs_log.txt 2>&1 


for sys in ${sys_to_test}; do

    rm -rf build_${sys} run_${sys} run_${sys}_molebench run_${sys}_disc_cylindrical
    mkdir  build_${sys} run_${sys} run_${sys}_molebench run_${sys}_disc_cylindrical
    build_and_prepare ${sys}

    echo "g95 environment variables are:"
    printenv | grep -i g95

    cd run_${sys}
    echo 
    echo "Running torus.${sys}"
    /usr/local/bin/mpirun -np 4 torus.${sys} > run_log_${sys}.txt 2>&1
    check_benchmark
    cd ..

    cd run_${sys}_molebench
    echo
    echo "Running torus.${sys} molebench"
    /usr/local/bin/mpirun -np 4 torus.${sys} > run_log_${sys}.txt 2>&1
    check_molebench
    cd ..

    sphbench
    echo "" > check_log 
    echo "SPHBENCH RESULTS" >> check_log
    check_benchmark >> check_log 2>&1 
    cd ../..

    cd run_${sys}_disc_cylindrical
    echo
    echo "Running torus.${sys} cylindrical polar disc benchmark"
    /usr/local/bin/mpirun -np 4 torus.${sys} > run_log_${sys}.txt 2>&1
    echo "" > check_log 
    echo "CYLINDIRCAL POLAR DISC BENCHMARK RESULTS" >> check_log
    check_benchmark >> check_log 2>&1 
    cd ..

done

exit

