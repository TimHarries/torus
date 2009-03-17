#!/bin/ksh
# Build the Torus executable for the system defined by the 
# SYSTEM environment variable
make_build()
{
mkdir build
cd    build 

echo "Building Torus for ${SYSTEM}"
log_file=compile_log.txt
ln -s ${TEST_DIR}/torus/* .
/usr/bin/make depends > ${log_file} 2>&1 
/usr/bin/make debug=${USEDEBUGFLAGS} >> ${log_file} 2>&1
if [[ $? -eq 0 ]]; then
# Count number of warnings. Subtract 2 because there are always warnings
# about include files from the make depends step (run twice).
    num_warn=`grep -i warning ${log_file} | wc -l | awk '{print $1 - 2}`
    echo "Compilation completed with ${num_warn} warnings."
else
    echo "Compilation failed."
fi

cd ..
}

# Build the Torus library for the system defined by the 
# SYSTEM environment variable
make_lib()
{
mkdir lib
cd    lib

echo "Building Torus library for ${SYSTEM}"
log_file=compile_log_lib.txt
ln -s ${TEST_DIR}/torus/* .
/usr/bin/make depends > ${log_file} 2>&1 
/usr/bin/make lib debug=${USEDEBUGFLAGS} >> ${log_file} 2>&1
if [[ $? -eq 0 ]]; then
# Count number of warnings. Subtract 2 because there are always warnings
# about include files from the make depends step (run twice).
    num_warn=`grep -i warning ${log_file} | wc -l | awk '{print $1 - 2}'`
    echo "Compilation completed with ${num_warn} warnings."
else
    echo "Compilation failed."
fi
cd ..

}

make_comparespec()
{
echo Compiling comparespec code
mkdir ${WORKING_DIR}/bin
cd    ${WORKING_DIR}/bin
cp ../benchmarks/disc/comparespec.f90 .
${TORUS_FC} -o comparespec comparespec.f90
}

run_bench()
{
cd ${WORKING_DIR}/benchmarks/${THIS_BENCH}
ln -s ${WORKING_DIR}/build/torus.${SYSTEM} . 

case ${SYSTEM} in
    ompi) mpirun -np 4 torus.ompi > run_log_${THIS_BENCH}.txt 2>&1 ;;

    intelmac) ./torus.intelmac  > run_log_${THIS_BENCH}.txt 2>&1 ;;

    zen) mpirun -np 4 torus.zen > run_log_${THIS_BENCH}.txt 2>&1 ;;

    *) echo "Unrecognised SYSTEM type. Aborting"
       exit 1;;
esac

#Tag the tune.dat file 
ln -s tune.dat tune_${THIS_BENCH}.txt 

}

run_hydro()
{
cd ${WORKING_DIR}/benchmarks/hydro
ln -s ${WORKING_DIR}/build/torus.${SYSTEM} . 

case ${SYSTEM} in
    ompi) mpirun -np 3 torus.ompi > run_log_hydro.txt 2>&1 ;;

    zen) mpirun -np 3 torus.zen > run_log_hydro.txt 2>&1 ;;

    *) echo "Unrecognised SYSTEM type. Aborting"
       exit 1;;
esac

}

run_sphbench()
{
cd ${WORKING_DIR}/benchmarks/sphbench
ln -s ${WORKING_DIR}/lib/libtorus.a 
ln -s ${WORKING_DIR}/lib/torus_mod.mod 
ln -s ${TEST_DIR}/torus/isochrones/iso* .
./compile

case ${SYSTEM} in
    ompi) mpirun -np 4 sphbench > run_log_sphbench.txt 2>&1;; 

    intelmac) ./sphbench > run_log_sphbench.txt 2>&1;;

    *) echo "Unrecognised SYSTEM type. Aborting"
       exit 1;;
esac

}

check_benchmark()
{

# Check we have the required benchmark results files
if [[ ! -e sed100_125.dat ]]; then
    cp ../disc/sed100_125.dat .
fi

if [[ ! -e sed100_775.dat ]]; then
    cp ../disc/sed100_775.dat .
fi

echo Comparing the 12.5 degree model...
cp test_inc013.dat speca.dat
cp sed100_125.dat specb.dat
${WORKING_DIR}/bin/comparespec

echo Comparing the 77.5 degree model...
cp test_inc077.dat speca.dat
cp sed100_775.dat specb.dat
${WORKING_DIR}/bin/comparespec
}

check_molebench()
{
echo Compiling compare_molbench code
${TORUS_FC} -o compare_molbench compare_molbench.f90
./compare_molbench
}

check_hydro()
{
echo Compiling compareSod code
${TORUS_FC} -o comparesod compareSod.f90
./comparesod > check_log_hydro.txt
}

run_torus_test_suite()
{
if [[ -e ${TEST_DIR} ]]; then
    echo "Removing old ${TEST_DIR}"
    rm -rf ${TEST_DIR}
fi

echo "Working directory is ${TEST_DIR}"
mkdir -p ${TEST_DIR}
cd ${TEST_DIR}

echo Checking out torus from CVS archive...
/usr/bin/cvs -q co torus > cvs_log.txt 2>&1 

for sys in ${SYS_TO_TEST}; do

    export SYSTEM=${sys}
    export WORKING_DIR=${TEST_DIR}/benchmarks_${SYSTEM}
    mkdir ${WORKING_DIR}
    cd    ${WORKING_DIR} 
    cp -r ${TEST_DIR}/torus/benchmarks . 

# Build code
    make_build
    make_lib
    make_comparespec

# Run hydro benchmark
    echo "Running hydro benchmark"
    run_hydro
    check_hydro

# Run benchmark tests
    echo "Running disc benchmark"
    export THIS_BENCH=disc
    run_bench 
    check_benchmark > check_log_${THIS_BENCH}.txt 2>&1 

#    if [[ ${SYSTEM} == ompi ]]; then
#	echo "Running cylindrical polar disc benchmark"
#	export THIS_BENCH=disc_cylindrical
#	run_bench
#	check_benchmark > check_log_${THIS_BENCH}.txt 2>&1 
#    fi

# SPH-Bench and HII region benchmark have been commented out 
# when running Torus v1.1 stable version tests. 

#    echo "Running SPH-Bench"
#    run_sphbench
#    check_benchmark > check_log_sphbench.txt 2>&1 

#    echo "Running HII region benchmark"
#    export THIS_BENCH=HII_region
#    run_bench

    echo "Running molecular benchmark"
    export THIS_BENCH=molebench 
    mkdir ${WORKING_DIR}/benchmarks/molebench/plots 
    run_bench
    check_molebench > check_log_${THIS_BENCH}.txt 2>&1 

done
}

print_help()
{
echo ""
echo "This script runs the torus test suite. Use the -s option to run the stable version tests."
echo "Use the -d option to run the daily tests (default)."
echo ""
}

########################################################################################################

# Default mode is daily test
export MODE=daily

# Parse command line arguments
while [ $# -gt 0 ]
do
    case "$1" in 
	-s) export MODE=stable;;
	-d) export MODE=daily;;
	-z) export MODE=zen;;
	-h) print_help
	    exit;;
    esac
shift
done

case ${MODE} in 

    daily) export SYS_TO_TEST="ompi"
	   export DEBUG_OPTS="yes"
	   export TORUS_FC="g95"
	   export PATH=/sw/bin:/usr/local/bin:${PATH}
	   echo TORUS daily test suite started on `date`
	   echo -------------------------------------------------------------------
	   echo;;

    stable) export SYS_TO_TEST="ompi intelmac"
            export DEBUG_OPTS="no yes"
	    export TORUS_FC="g95"
	    export PATH=/sw/bin:/usr/local/bin:${PATH}
	    echo TORUS stable version tests started on `date`
	    echo -------------------------------------------------------------------
	    echo;;

    zen) export SYS_TO_TEST="zen"
	 export DEBUG_OPTS="yes"
	 export TORUS_FC="ifort"
	 echo TORUS zen tests started on `date`
	 echo -------------------------------------------------------------------
	 echo;;

    *)  echo "ERROR: unrecognised mode"
	exit 1;;
esac

export CVSROOT=${USER}@pinky.astro.ex.ac.uk:/h/th/CVS
export CVS_RSH=ssh

for opt in ${DEBUG_OPTS}; do
    export USEDEBUGFLAGS=${opt}

# Set name of output directory
    case ${MODE} in 
	daily)  export TEST_DIR=${HOME}/SCRATCH/torus_daily_test;;
	stable) export TEST_DIR=${HOME}/SCRATCH/torus_stable_version_tests/debug=${USEDEBUGFLAGS};;
	zen)    export TEST_DIR=/scratch/acreman/torus_tests/debug=${USEDEBUGFLAGS};;
    esac

    export TORUS_DATA=${TEST_DIR}/torus/data

# Set floating point exception flags for g95
    case ${TORUS_FC} in 
	g95) 
	    case ${USEDEBUGFLAGS} in
		yes) export G95_FPU_INVALID=true
		    export G95_FPU_ZERODIV=true
		    export G95_FPU_OVERFLOW=true;;

		no) export G95_FPU_INVALID=false
		    export G95_FPU_ZERODIV=false
		    export G95_FPU_OVERFLOW=false;;

	    esac

	    echo
	    echo "G95 environment variables are"
	    printenv | grep G95
	    echo

    esac

    run_torus_test_suite

done

exit

