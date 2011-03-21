#!/bin/ksh
# Build the Torus executable for the system defined by the 
# SYSTEM environment variable
make_build()
{
mkdir build
cd    build 

echo "Building Torus for ${SYSTEM}"
log_file=compile_log_${SYSTEM}.txt
ln -s ${TEST_DIR}/torus/* .
/usr/bin/make depends > ${log_file} 2>&1 

case ${SYSTEM} in
    gfortran) /usr/bin/make debug=${USEDEBUGFLAGS} openmp=yes >> ${log_file} 2>&1;;
    ompiosx) /usr/bin/make debug=${USEDEBUGFLAGS} openmp=yes >> ${log_file} 2>&1;;
    *) /usr/bin/make debug=${USEDEBUGFLAGS} >> ${log_file} 2>&1;;
esac

if [[ $? -eq 0 ]]; then
# Count number of warnings. Subtract 2 because there are always warnings
# about include files from the make depends step (run twice).
    num_warn=`grep -i warning ${log_file} | wc -l | awk '{print $1 - 2}`
    echo "Compilation completed with ${num_warn} warnings."
else
    echo "Compilation failed."
    exit 2
fi

cd ..
}

# Build the Torus library for the system defined by the 
# SYSTEM environment variable
#make_lib()
#{
#mkdir lib
#cd    lib
#
#echo "Building Torus library for ${SYSTEM}"
#log_file=compile_log_lib.txt
#ln -s ${TEST_DIR}/torus/* .
#/usr/bin/make depends > ${log_file} 2>&1 
#/usr/bin/make lib debug=${USEDEBUGFLAGS} >> ${log_file} 2>&1
#if [[ $? -eq 0 ]]; then
## Count number of warnings. Subtract 2 because there are always warnings
## about include files from the make depends step (run twice).
#    num_warn=`grep -i warning ${log_file} | wc -l | awk '{print $1 - 2}'`
#    echo "Compilation completed with ${num_warn} warnings."
#else
#    echo "Compilation failed."
#    exit 1
#fi
#cd ..
#
#}

make_comparespec()
{
echo Compiling comparespec code
mkdir ${WORKING_DIR}/bin
cd    ${WORKING_DIR}/bin
cp ../benchmarks/disc/comparespec.f90 .
${TORUS_FC} -o comparespec comparespec.f90
echo
}

run_bench()
{
cd ${WORKING_DIR}/benchmarks/${THIS_BENCH}
ln -s ${WORKING_DIR}/build/torus.${SYSTEM} .
log_file=run_log_${SYSTEM}_${THIS_BENCH}.txt
export TORUS_JOB_DIR=./

case ${SYSTEM} in
    ompi) mpirun -np 4 torus.ompi > ${log_file} 2>&1 ;;
    ompiosx) mpirun -np 2 torus.ompiosx > ${log_file} 2>&1 ;;
    zen) mpirun -np 8 torus.zen > ${log_file} 2>&1 ;;
    nagfor) ./torus.${SYSTEM} > ${log_file} 2>&1 &;;
    *) ./torus.${SYSTEM} > ${log_file} 2>&1 ;;
esac

#Tag the tune.dat file 
ln -s tune.dat tune_${SYSTEM}_${THIS_BENCH}.txt 

}

run_hydro()
{
cd ${WORKING_DIR}/benchmarks/hydro
ln -s ${WORKING_DIR}/build/torus.${SYSTEM} . 
mpirun -np 3 torus.${SYSTEM} > run_log_hydro.txt 2>&1
ln -s tune.dat tune_hydro.txt
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
    g95) ./sphbench > run_log_sphbench.txt 2>&1;;
    zen) mpirun -np 8 sphbench > run_log_sphbench.txt 2>&1 ;;
    *) echo "Unrecognised SYSTEM type. Skipping this test.";;
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
./comparesod
}

check_hII()
{
echo Compiling comparelex code
${TORUS_FC} -o comparelex comparelex.f90
./comparelex
}

prepare_run()
{
if [[ -e ${TEST_DIR} ]]; then
    if [[ ${CLOBBEROK} == yes ]]; then
	echo "Removing old ${TEST_DIR}"
	rm -rf ${TEST_DIR}
    else
	echo "${TEST_DIR} already exisits. Aborting"
	exit 1
    fi
fi

echo "Working directory is ${TEST_DIR}"
mkdir -p ${TEST_DIR}
cd ${TEST_DIR}

echo Checking out torus from SVN archive...
/usr/bin/svn checkout https://repository.astro.ex.ac.uk/torus/trunk/torus/ torus > svn_log.txt 2>&1 
grep "Checked out revision" svn_log.txt
}

run_torus_test_suite()
{

echo 
echo "Running benchmarks"
echo 

for sys in ${SYS_TO_TEST}; do

    export SYSTEM=${sys}
    export WORKING_DIR=${TEST_DIR}/benchmarks_${SYSTEM}
    mkdir ${WORKING_DIR}
    cd    ${WORKING_DIR} 
    cp -r ${TEST_DIR}/torus/benchmarks . 

    if [[ $SYSTEM == ompiosx ]]; then
	export OLD_PATH=${PATH}
	export PATH=/opt/ompi/gfortran/bin:${PATH}
	export OMP_NUM_THREADS=2
    fi

# Build code
    if [[ ${DO_BUILD} == yes ]]; then
	make_build
#        make_lib
    else
	mkdir build
	cd build 
	ln -s ${TORUS_BINARY}
	cd ..
    fi

    make_comparespec

# Run hydro benchmark
    case ${SYSTEM} in
	ompi|zen)  echo "Running hydro benchmark"
	    run_hydro
	    check_hydro  > check_log_${SYSTEM}_hydro.txt 2>&1 
	    cat check_log_${SYSTEM}_hydro.txt
	    echo ;;
	*) echo "Hydro benchmark does not run on this system. Skipping"
	    echo ;;
    esac

# Run 2D disc
    echo "Running disc benchmark"
    export THIS_BENCH=disc
    run_bench 
    if [[ ${SYSTEM} != "nagfor" ]]; then
	check_benchmark > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1 
	cat check_log_${SYSTEM}_${THIS_BENCH}.txt
    fi
    echo

    echo "Running HII region benchmark"
    export THIS_BENCH=HII_region
    run_bench
    if [[ ${SYSTEM} != "nagfor" ]]; then
	check_hII > check_log_${SYSTEM}_hII.txt 2>&1 
	cat check_log_${SYSTEM}_hII.txt
    fi
    echo

    echo "Running molecular benchmark"
    export THIS_BENCH=molebench 
    run_bench
    if [[ ${SYSTEM} != "nagfor" ]]; then
	check_molebench > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1 
	tail check_log_${SYSTEM}_${THIS_BENCH}.txt # Lots of output so tail this file
    fi
    echo

# Only run these tests for MPI systems and not in the daily test
    if [[ ${MODE} != daily ]]; then
	if [[ ${SYSTEM} != "nagfor" && ${SYSTEM} != "g95" ]]; then 
	    echo "Running cylindrical polar disc benchmark"
	    export THIS_BENCH=disc_cylindrical
	    run_bench
	    check_benchmark > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1
	    cat check_log_${SYSTEM}_${THIS_BENCH}.txt
	    echo 

#	echo "Running SPH-Bench"
#	run_sphbench
#	check_benchmark > check_log_sphbench.txt 2>&1 

	fi
    fi

    if [[ $SYSTEM == ompiosx ]]; then
	export PATH=${OLD_PATH}
	unset OMP_NUM_THREADS 
    fi

    if [[ $SYSTEM == nagfor ]]; then
	echo "Waiting for nagfor runs to complete"
	echo
	wait

	cd ${WORKING_DIR}/benchmarks/disc
	check_benchmark > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1 
	cat check_log_${SYSTEM}_${THIS_BENCH}.txt
	echo

	cd ${WORKING_DIR}/benchmarks/HII_region
	check_hII > check_log_${SYSTEM}_hII.txt 2>&1 
	cat check_log_${SYSTEM}_hII.txt
	echo

	cd ${WORKING_DIR}/benchmarks/molebench
	check_molebench > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1 
	tail check_log_${SYSTEM}_${THIS_BENCH}.txt
	echo 

    fi

done
}

build_only_tests()
{

echo 
echo "Running build-only tests"
echo 

for sys in ${BUILD_ONLY}; do

    if [[ $SYSTEM == ompiosx ]]; then
	export OLD_PATH=${PATH}
	export PATH=/opt/ompi/gfortran/bin:${PATH}
    fi

    export SYSTEM=${sys}
    export WORKING_DIR=${TEST_DIR}/build_only_${SYSTEM}
    mkdir ${WORKING_DIR}
    cd    ${WORKING_DIR} 

    make_build
#    make_lib

    if [[ $SYSTEM == ompiosx ]]; then
	export PATH=${OLD_PATH}
    fi

done
}

print_help()
{
echo ""
echo "This script runs the torus test suite."
echo ""
echo "Use the -d option to run the daily tests (default)."
echo "Use the -s option to run the stable version tests."
echo "Use the -z option to run the tests on zen."
echo "Use the -b option to run build tests only"
echo "Use -e followed by full path to a torus executable to use a pre-built binary"
echo ""
}

########################################################################################################

# Default mode is daily test
export MODE=daily
export DO_BUILD=yes
export CLOBBEROK=yes

# If we're running on Zen then set the appropriate mode
this_host=`hostname`
if [[ ${this_host} == service0 ]]; then
    echo "Don't run this script on the log in node!"
    exit 1
elif [[ ${this_host} == service2 ]]; then
    export MODE=zen
fi

# Parse command line arguments
while [ $# -gt 0 ]
do
    case "$1" in 
	-s) export MODE=stable;;
	-d) export MODE=daily;;
	-z) export MODE=zen;;
	-b) export MODE=build;;
	-e) export DO_BUILD=no
	    export CLOBBEROK=no
	    shift 
	    TORUS_BINARY=$1;;
	-h) print_help
	    exit;;
    esac
shift
done

case ${MODE} in 

    daily) export SYS_TO_TEST="ompi gfortran ompiosx"
           export BUILD_ONLY="g95 nagfor"
	   export DEBUG_OPTS="yes"
	   export TORUS_FC="g95"
	   export PATH=~/bin:/usr/local/bin:${PATH}:/usr/bin
	   export NAG_KUSARI_FILE=${HOME}/NAG/nag.licence
	   echo TORUS daily test suite started on `date`
	   echo -------------------------------------------------------------------
	   echo;;

    build) export SYS_TO_TEST=" "
           export BUILD_ONLY="g95 nagfor ompi gfortran ompiosx"
	   export DEBUG_OPTS="yes"
	   export TORUS_FC="g95"
	   export PATH=~/bin:/usr/local/bin:${PATH}:/usr/bin
	   export NAG_KUSARI_FILE=${HOME}/NAG/nag.licence
	   echo TORUS build tests started on `date`
	   echo -------------------------------------------------------------------
	   echo;;

    stable) export SYS_TO_TEST="nagfor ompi ompiosx gfortran"
	    export BUILD_ONLY=""
            export DEBUG_OPTS="yes no"
	    export TORUS_FC="g95"
	    export PATH=~/bin:/usr/local/bin:${PATH}
	    echo TORUS stable version tests started on `date`
	    echo -------------------------------------------------------------------
	    echo;;

    zen) export SYS_TO_TEST="zen"
         export BUILD_ONLY=""
	 export DEBUG_OPTS="no yes"
	 export TORUS_FC="ifort"
	 echo TORUS zen tests started on `date`
	 echo -------------------------------------------------------------------
	 echo;;

    *)  echo "ERROR: unrecognised mode"
	exit 1;;
esac


if [[ $DO_BUILD == no ]]; then
    if [[ ! -x ${TORUS_BINARY} ]]; then
	echo "ERROR: ${TORUS_BINARY} is not an executable file"
	exit 1
    else
	SYS_TO_TEST=`basename ${TORUS_BINARY} | tr '.' ' ' | awk '{print $2}`
	echo "Using pre-built binary ${TORUS_BINARY}"
	echo "System type is ${SYS_TO_TEST}"
	export TORUS_FC="g95"   # Make sure FPE flags are set in case this g95
	export DEBUG_OPTS="yes" # Just sets flags if pre-built binary is used
    fi
fi 

for opt in ${DEBUG_OPTS}; do
    export USEDEBUGFLAGS=${opt}

# Set name of output directory
    case ${MODE} in 
	daily)  export TEST_DIR=${HOME}/SCRATCH/torus_daily_test;;
	stable) export TEST_DIR=${HOME}/SCRATCH/torus_stable_version_tests/debug=${USEDEBUGFLAGS};;
	zen)    export TEST_DIR=/scratch/${USER}/torus_tests/debug=${USEDEBUGFLAGS};;
	build)  export TEST_DIR=${HOME}/SCRATCH/torus_build_tests
    esac

    export TORUS_DATA=${TEST_DIR}/torus/data

# Set floating point exception flags for g95
    case ${TORUS_FC} in 
	g95) 
	    case ${USEDEBUGFLAGS} in
		yes) export G95_FPU_INVALID=true
		    export G95_FPU_ZERODIV=true
		    export G95_FPU_OVERFLOW=true;;
#		    export G95_MEM_INIT=NAN;;

		no) export G95_FPU_INVALID=false
		    export G95_FPU_ZERODIV=false
		    export G95_FPU_OVERFLOW=false;;
#		    export G95_MEM_INIT=NONE;;

	    esac

	    echo
	    echo "G95 environment variables are"
	    printenv | grep G95
	    echo

    esac

# Set up working dir and check out source code
    prepare_run

# Test build but don't run benchmarks
    if [[ ${DO_BUILD} == yes ]]; then
	build_only_tests
    fi

# Run benchmark tests
    run_torus_test_suite

done

exit

