#!/bin/ksh
# Build the Torus executable for the system defined by the 
# SYSTEM environment variable
make_build()
{
mkdir build
cd    build 

echo "Building Torus for ${SYSTEM} OpenMP=${USEOPENMP}"
log_file=compile_log_${SYSTEM}.txt
rsync -a ${TEST_DIR}/torus/* .
rsync -a ${TEST_DIR}/torus/.svn .

/usr/bin/make depends > ${log_file} 2>&1 
/usr/bin/make debug=${USEDEBUGFLAGS} openmp=${USEOPENMP} coverage=${USEGCOV} >> ${log_file} 2>&1

if [[ $? -eq 0 ]]; then
# Count number of warnings. Subtract 2 because there are always warnings
# about include files from the make depends step (run twice).
    num_warn=`grep -i warning ${log_file} | wc -l | awk '{print $1 - 2}`
    echo "Compilation completed with ${num_warn} warnings."
else
    echo "Compilation failed."
    export RETURN_CODE=2
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

echo Compiling check_disc_image
cp ../benchmarks/disc/check_disc_image.f90 . 
${TORUS_FC} -o check_disc_image check_disc_image.f90 -lcfitsio -L${TORUS_FITSLIBS}

echo
}

run_bench()
{
cd ${WORKING_DIR}/benchmarks/${THIS_BENCH}
ln -s ${WORKING_DIR}/build/torus.${SYSTEM} .
log_file=run_log_${SYSTEM}_${THIS_BENCH}.txt
export TORUS_JOB_DIR=./

case ${SYSTEM} in
    ompiosx|zen) mpirun -np ${NUM_MPI_PROC} torus.${SYSTEM} > ${log_file} 2>&1 ;;
    gfortran) ./torus.${SYSTEM} > ${log_file} 2>&1 ;;
    nagfor) ./torus.${SYSTEM} > ${log_file} 2>&1 &;;
    *) echo "Unrecognised SYSTEM type. Skipping this test.";;
esac

#Rename the tune.dat file 
mv tune.dat tune_${SYSTEM}_${THIS_BENCH}.txt 

}

# Run Torus with a domain decomposed grid.  
run_dom_decomp()
{
cd ${WORKING_DIR}/benchmarks/${THIS_BENCH}
ln -s ${WORKING_DIR}/build/torus.${SYSTEM} . 
mpirun -np $1 torus.${SYSTEM} > run_log_${THIS_BENCH}.txt 2>&1
mv tune.dat tune_${THIS_BENCH}.txt
}

setup_sphbench()
{
cd ${WORKING_DIR}/benchmarks/sphbench
./setup_write_sph_file
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

# Check first SED
echo Comparing the 12.5 degree model...
cp test_inc013.dat speca.dat
cp sed100_125.dat specb.dat
${WORKING_DIR}/bin/comparespec

# Check second SED
echo Comparing the 77.5 degree model...
cp test_inc077.dat speca.dat
cp sed100_775.dat specb.dat
${WORKING_DIR}/bin/comparespec

# Check image
${WORKING_DIR}/bin/check_disc_image
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

check_image()
{
echo "Generating analytical solution"
${TORUS_FC} -o cylinder_test cylinder_test.f90
./cylinder_test
echo "Checking Torus result"
${TORUS_FC} -o comparison comparison.f90
./comparison
}

check_it()
{
echo "Compiling check.f90"
${TORUS_FC} -o check check.f90
./check
}

# Check that Torus completed OK. For some tests the results file can be written 
# even if Torus has bugged out. 
check_completion()
{
    grep "Torus completed" run_log_${SYSTEM}_${THIS_BENCH}.txt > /dev/null
    if [[ $? -eq 0 ]]; then
	echo "Torus completed OK"
    else
	echo "WARNING: Torus did not complete"
    fi
}

prepare_run()
{
if [[ -e ${TEST_DIR} ]]; then
    if [[ -e ${TEST_DIR}/lock ]]; then
	echo "Found lock file. Aborting"
	exit 1
    fi
    if [[ ${CLOBBEROK} == yes ]]; then
	echo "Removing old ${TEST_DIR}"
	rm -rf ${TEST_DIR}
    else
	echo "${TEST_DIR} already exists. Aborting"
	exit 1
    fi
fi

echo "Working directory is ${TEST_DIR}"
mkdir -p ${TEST_DIR}
cd ${TEST_DIR}
touch lock

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

# Details of how to run each system are set here
    export USEGCOV=no
    if [[ ${sys} == "ompiosx-openmp" ]]; then
	export SYSTEM=ompiosx
	export USEOPENMP=yes
	export NUM_MPI_PROC=2
	export OMP_NUM_THREADS=2
	echo "Building ompiosx with OpenMP"
    elif [[ ${sys} == "ompiosx" ]]; then
	export SYSTEM=ompiosx
	export USEOPENMP=no
	export NUM_MPI_PROC=4
	export USEGCOV=yes
	echo "Building ompiosx without OpenMP"
    elif [[ ${sys} == "gfortran" ]]; then
	export SYSTEM=gfortran
	export USEOPENMP=yes
	export OMP_NUM_THREADS=4
    elif [[ ${sys} == "nagfor" ]]; then
	export SYSTEM=nagfor
	export USEOPENMP=no
    elif [[ ${sys} == "zen" ]]; then
	export SYSTEM=zen
	export USEOPENMP=no
	export NUM_MPI_PROC=8
    else
	echo "${sys} not recognised. Aborting."
	exit 1
    fi

    export WORKING_DIR=${TEST_DIR}/benchmarks_${sys}
    mkdir ${WORKING_DIR}
    cd    ${WORKING_DIR} 
    cp -r ${TEST_DIR}/torus/benchmarks . 

# Build code
    if [[ ${DO_BUILD} == yes ]]; then
	make_build
    else
	mkdir build
	cd build 
	ln -s ${TORUS_BINARY}
	cd ..
    fi

# Check if we have an executable. If not proceed to the next SYSTEM
    if [[ -x build/torus.${SYSTEM} ]]; then
	echo "Found executable file torus.${SYSTEM}"
    else
	echo "Executable not found. Skipping ${SYSTEM}"
	continue
    fi

    make_comparespec

# Run hydro benchmark
    case ${sys} in
	ompiosx|zen)  echo "Running hydro benchmark"
	    export THIS_BENCH=hydro
	    run_dom_decomp 3
	    check_hydro  > check_log_${SYSTEM}_hydro.txt 2>&1 
	    cat check_log_${SYSTEM}_hydro.txt
	    echo ;;
	*) echo "Hydro benchmark does not run on this system. Skipping"
	    echo ;;
    esac

# Domain decomposed Lexington benchmark
    case ${sys} in
	ompiosx|zen)  echo "Running domain decomposed HII region benchmark"
	    export THIS_BENCH=HII_regionMPI
	    run_dom_decomp 3
	    check_hII > check_log_${SYSTEM}_hII_MPI.txt 2>&1 
	    cat check_log_${SYSTEM}_hII_MPI.txt
	    echo ;;
	*) echo "Domain decomposed HII region does not run on this system. Skipping"
	    echo ;;
    esac

# Imaging test
    case ${sys} in
	ompiosx|zen)  echo "Running imaging benchmark"
	    export THIS_BENCH=cylinder_image_test
	    run_dom_decomp 9
	    check_image > check_log_${SYSTEM}_image.txt 2>&1 
	    cat check_log_${SYSTEM}_image.txt
	    echo ;;
	*) echo "Imaging benchmark does not run on this system. Skipping"
	    echo ;;
    esac

# Gravity solver test
    case ${sys} in
	ompiosx|zen)  echo "Running gravity solver test"
	    export THIS_BENCH=gravtest
	    run_dom_decomp 9
	    check_it > check_log_${SYSTEM}_gravtest.txt 2>&1 
	    cat check_log_${SYSTEM}_gravtest.txt
	    echo ;;
	*) echo "Gravity solver test does not run on this system. Skipping"
	    echo ;;
    esac

# N body test
    echo "Running N body test"
    export THIS_BENCH=nbody
    run_bench
    check_it > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1 
    cat check_log_${SYSTEM}_${THIS_BENCH}.txt
    echo

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
	check_completion
    fi
    echo

    echo "Running SPH to grid test"
    export THIS_BENCH=sphbench
    setup_sphbench
    run_bench
    echo

# Only run this test for MPI systems and not in the daily test
    if [[ ${MODE} != daily ]]; then
	if [[ ${SYSTEM} != "nagfor" ]]; then 
	    echo "Running cylindrical polar disc benchmark"
	    export THIS_BENCH=disc_cylindrical
	    run_bench
	    check_benchmark > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1
	    cat check_log_${SYSTEM}_${THIS_BENCH}.txt
	    echo 
	fi
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
	check_completion
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

# Details of how to build each system are set here
    export USEGCOV=no
    if [[ ${sys} == "ompiosx-openmp" ]]; then
	export SYSTEM=ompiosx
	export USEOPENMP=yes
    elif [[ ${sys} == "ompiosx" ]]; then
	export SYSTEM=ompiosx
	export USEOPENMP=no
    elif [[ ${sys} == "gfortran" ]]; then
	export SYSTEM=gfortran
	export USEOPENMP=yes
    elif [[ ${sys} == "nagfor" ]]; then
	export SYSTEM=nagfor
	export USEOPENMP=no
    elif [[ ${sys} == "zen" ]]; then
	export SYSTEM=zen
	export USEOPENMP=no
    else
	echo "${sys} not recognised. Aborting."
	exit 1
    fi

    export WORKING_DIR=${TEST_DIR}/build_only_${sys}
    mkdir ${WORKING_DIR}
    cd    ${WORKING_DIR} 

    make_build
    echo " " 

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
export RETURN_CODE=0 

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

    daily) export SYS_TO_TEST="gfortran ompiosx ompiosx-openmp"
           export BUILD_ONLY="nagfor"
	   export DEBUG_OPTS="yes"
	   export TORUS_FC="gfortran -g -fcheck=all"
	   export TORUS_FITSLIBS="/Users/acreman/cfitsio"
	   export PATH=~/bin:/usr/local/bin:${PATH}:/usr/bin
	   export NAG_KUSARI_FILE=/Users/acreman/NAG/nag.licence
	   echo TORUS daily test suite started on `date`
	   echo -------------------------------------------------------------------
	   echo;;

    build) export SYS_TO_TEST=" "
           export BUILD_ONLY="nagfor gfortran ompiosx ompiosx-openmp"
	   export DEBUG_OPTS="yes"
	   export TORUS_FC="gfortran -g -fcheck=all"
	   export TORUS_FITSLIBS="/Users/acreman/cfitsio"
	   export PATH=~/bin:/usr/local/bin:${PATH}:/usr/bin
	   export NAG_KUSARI_FILE=/Users/acreman/NAG/nag.licence
	   echo TORUS build tests started on `date`
	   echo -------------------------------------------------------------------
	   echo;;

    stable) export SYS_TO_TEST="nagfor gfortran ompiosx ompiosx-openmp"
	    export BUILD_ONLY=""
            export DEBUG_OPTS="yes no"
	    export TORUS_FC="gfortran -g -fcheck=all"
	    export TORUS_FITSLIBS="/Users/acreman/cfitsio"
	    export PATH=~/bin:/usr/local/bin:${PATH}
	    echo TORUS stable version tests started on `date`
	    echo -------------------------------------------------------------------
	    echo;;

    zen) export SYS_TO_TEST="zen"
         export BUILD_ONLY=""
	 export DEBUG_OPTS="no yes"
	 export TORUS_FC="ifort"
	 export TORUS_FITSLIBS="/home/tjharrie/cfitsio"
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

# Set up working dir and check out source code
    prepare_run

# Test build but don't run benchmarks
    if [[ ${DO_BUILD} == yes ]]; then
	build_only_tests
    fi

# Run benchmark tests
    run_torus_test_suite

done

rm ${TEST_DIR}/lock

exit ${RETURN_CODE}

