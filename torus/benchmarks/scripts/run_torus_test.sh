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
    echo "*** Compile log follows ****************************"
    cat ${log_file}
    echo "*** End compile log ********************************"
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

run_molecularRestart()
{
cd ${WORKING_DIR}/benchmarks/molebench/restart

# Copy/link required files from original run directory
cp ../testDump_HCO+_grid.grid .
cp ../testDump_restart.dat restart.dat
ln -s ../model_1.dat
ln -s ../moltest.dat
ln -s ../compare_molbench.f90
ln -s ../check_cube.f90

log_file=run_log_${SYSTEM}_${THIS_BENCH}.txt
export TORUS_JOB_DIR=./
ln -s ${WORKING_DIR}/build/torus.${SYSTEM} .

case ${SYSTEM} in
    ompiosx) mpirun -np ${NUM_MPI_PROC} torus.${SYSTEM} > ${log_file} 2>&1 ;;
    gfortran) ./torus.${SYSTEM} > ${log_file} 2>&1 ;;
    *) echo "Unrecognised SYSTEM type. Skipping this test.";;
esac

#Rename the tune.dat file 
mv tune.dat tune_${SYSTEM}_${THIS_BENCH}.txt 
}

run_bench()
{
cd ${WORKING_DIR}/benchmarks/${THIS_BENCH}
ln -s ${WORKING_DIR}/build/torus.${SYSTEM} .
log_file=run_log_${SYSTEM}_${THIS_BENCH}.txt
export TORUS_JOB_DIR=./

case ${SYSTEM} in
    ompiosx) mpirun -np ${NUM_MPI_PROC} torus.${SYSTEM} > ${log_file} 2>&1 ;;
    gfortran) ./torus.${SYSTEM} > ${log_file} 2>&1 ;;
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
# The check_completion function expects the run log to be tagged with the SYSTEM
ln -s run_log_${THIS_BENCH}.txt run_log_${SYSTEM}_${THIS_BENCH}.txt
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

echo Compiling check_cube code
${TORUS_FC} -o check_cube check_cube.f90 -lcfitsio -L${TORUS_FITSLIBS}
./check_cube
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

check_angImg()
{
echo "Compiling check.f90 for angular image test"
${TORUS_FC} -o check check.f90 -lcfitsio -L${TORUS_FITSLIBS}
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

echo "Working directory is ${TEST_DIR}"

if [[ -e ${TEST_DIR} ]]; then
    if [[ -e ${TEST_DIR}/lock ]]; then
	echo "Found lock file ${TEST_DIR}/lock. Aborting"
	exit 1
    fi
    echo "Removing old ${TEST_DIR}"
    rm -rf ${TEST_DIR}
fi

if [[ ${TORUS_WORKING_COPY} == none ]]; then 
    mkdir -p ${TEST_DIR}
    cd ${TEST_DIR}
    touch lock
    echo Checking out torus from SVN archive using: ${TORUS_SVN_REVISION} ${TORUS_SVN_PATH}
    /usr/bin/svn checkout ${TORUS_SVN_REVISION} ${TORUS_SVN_PATH} torus > svn_log.txt 2>&1 
    grep "Checked out revision" svn_log.txt
else
    if [[ -e ${TORUS_WORKING_COPY}/torus/torusMainV2.F90 ]]; then
	echo "Taking source code from working copy in ${TORUS_WORKING_COPY}"
	mkdir -p ${TEST_DIR}
	cd ${TEST_DIR}
	touch lock
	ln -s ${TORUS_WORKING_COPY}/torus . 
    else
	echo "Did not find torusMainV2.F90 in ${TORUS_WORKING_COPY}/torus "
	echo "This doesn't look like a working copy. Aborting ..."
	exit 1
    fi
fi
}

run_torus_test_suite()
{

for sys in ${SYS_TO_TEST}; do

# Details of how to run each system are set here
    export USEGCOV=no
# OpenMPI will bind to core by default so we'll use the 2 OpenMP threads available from hyperthreading.
# Binding can be disabled by giving mpirun the '--bind-to none' option if more OpenMP threads are required.
    if [[ ${sys} == "ompiosx-openmp" ]]; then
	export SYSTEM=ompiosx
	export USEOPENMP=yes
	export NUM_MPI_PROC=4
	export OMP_NUM_THREADS=2
    elif [[ ${sys} == "ompiosx" ]]; then
	export SYSTEM=ompiosx
	export USEOPENMP=no
	export NUM_MPI_PROC=8
	export USEGCOV=yes
    elif [[ ${sys} == "gfortran" ]]; then
	export SYSTEM=gfortran
	export USEOPENMP=yes
	export OMP_NUM_THREADS=8
    else
	echo "${sys} not recognised. Aborting."
	exit 1
    fi

    echo
    echo "======================================"
    echo "Running tests for system ${SYSTEM}"
    echo "OpenMP: ${USEOPENMP}"
    echo "======================================"
    echo

    export WORKING_DIR=${TEST_DIR}/benchmarks_${sys}
    mkdir ${WORKING_DIR}
    cd    ${WORKING_DIR} 
    cp -r ${TEST_DIR}/torus/benchmarks . 

# Build code
    make_build

# Check if we have an executable. If not proceed to the next SYSTEM
    if [[ -x build/torus.${SYSTEM} ]]; then
	echo "Found executable file torus.${SYSTEM}"
    else
	echo "Executable not found. Skipping ${SYSTEM}"
	continue
    fi

# Build comparison code for disc benchmark. Be careful moving this as it changes the cwd. 
    make_comparespec

# Run hydro benchmark
    case ${sys} in
	ompiosx)  echo "Running hydro benchmark"
	    export THIS_BENCH=hydro
	    run_dom_decomp 3
	    check_hydro  > check_log_${SYSTEM}_hydro.txt 2>&1 
	    cat check_log_${SYSTEM}_hydro.txt
	    check_completion
	    echo ;;
	*) echo "Hydro benchmark does not run on this system. Skipping"
	    echo ;;
    esac

# Domain decomposed Lexington benchmark
    case ${sys} in
	ompiosx)  echo "Running domain decomposed HII region benchmark"
	    export THIS_BENCH=HII_regionMPI
	    run_dom_decomp 9
	    check_hII > check_log_${SYSTEM}_hII_MPI.txt 2>&1 
	    cat check_log_${SYSTEM}_hII_MPI.txt
	    check_completion
	    echo ;;
	*) echo "Domain decomposed HII region does not run on this system. Skipping"
	    echo ;;
    esac

# Imaging test
    case ${sys} in
	ompiosx)  echo "Running imaging benchmark"
	    export THIS_BENCH=cylinder_image_test
	    run_dom_decomp 9
	    check_image > check_log_${SYSTEM}_image.txt 2>&1 
	    cat check_log_${SYSTEM}_image.txt
            check_completion
	    echo ;;
	*) echo "Imaging benchmark does not run on this system. Skipping"
	    echo ;;
    esac

# Gravity solver test
    case ${sys} in
	ompiosx)  echo "Running 3D gravity solver test"
	    export THIS_BENCH=gravtest
	    run_dom_decomp 9
	    check_it > check_log_${SYSTEM}_gravtest.txt 2>&1 
	    tail check_log_${SYSTEM}_gravtest.txt
	    check_completion
	    echo ;;
	*) echo "Gravity solver test does not run on this system. Skipping"
	    echo ;;
    esac

# Gravity solver test in 2D
    case ${sys} in
	ompiosx)  echo "Running 2D gravity solver test"
	    export THIS_BENCH=gravtest_2d
	    run_dom_decomp 5
	    check_it > check_log_${SYSTEM}_gravtest_2d.txt 2>&1 
	    tail check_log_${SYSTEM}_gravtest_2d.txt
	    check_completion
	    echo ;;
	*) echo "2D Gravity solver test does not run on this system. Skipping"
	    echo ;;
    esac

# N body test
    echo "Running N body test"
    export THIS_BENCH=nbody
    run_bench
    check_it > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1 
    cat check_log_${SYSTEM}_${THIS_BENCH}.txt
    check_completion
    echo

# Run 2D disc
    echo "Running disc benchmark"
    export THIS_BENCH=disc
    run_bench 
    check_benchmark > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1 
    cat check_log_${SYSTEM}_${THIS_BENCH}.txt
    check_completion
    echo

    echo "Running HII region benchmark"
    export THIS_BENCH=HII_region
    run_bench
    check_hII > check_log_${SYSTEM}_hII.txt 2>&1 
    cat check_log_${SYSTEM}_hII.txt
    check_completion
    echo

    echo "Running molecular benchmark"
    export THIS_BENCH=molebench 
    run_bench
    check_molebench > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1 
    tail -20 check_log_${SYSTEM}_${THIS_BENCH}.txt # Lots of output so tail this file
    check_completion
    echo

    echo "Running molecular restart test"
    export THIS_BENCH=moleRestart
    run_molecularRestart
    check_molebench > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1
    ./check_nrays.sh run_log_${SYSTEM}_${THIS_BENCH}.txt >> check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1 
    tail -20 check_log_${SYSTEM}_${THIS_BENCH}.txt # Lots of output so tail this file
    check_completion
    echo

    echo "Running SPH to grid test"
    export THIS_BENCH=sphbench
    setup_sphbench
    run_bench
    ./checkSphToGrid.pl run_log_${SYSTEM}_${THIS_BENCH}.txt > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1 
    cat check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1
    check_completion
    echo

    echo "Running SPH to grid test (binary dump with chemistry)"
    export THIS_BENCH=sphToGridBinary
    cd ${WORKING_DIR}/benchmarks/sphToGridBinary
    ln -s ${HOME}/torus_dev/forTestSuite/SQA0321
    run_bench
    ./checkSphToGridChem.pl run_log_${SYSTEM}_${THIS_BENCH}.txt > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1 
    cat check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1
    check_completion
    echo

    echo "Running restart test"
    export THIS_BENCH=restart
    cd ${WORKING_DIR}/benchmarks/restart
    ln -s ../disc/lucy_grid_tmp.dat
    run_bench
    check_benchmark > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1 
    cat check_log_${SYSTEM}_${THIS_BENCH}.txt
    check_completion
    echo

    echo "Running angular imaging test"
    export THIS_BENCH=angularImageTest
    run_bench
    check_angImg > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1 
    cat check_log_${SYSTEM}_${THIS_BENCH}.txt
    check_completion
    echo

# Only run in stable version tests as this is slow
    if [[ ${MODE} == stable ]]; then
	echo "Running cylindrical polar disc benchmark"
	export THIS_BENCH=disc_cylindrical
	run_bench
	check_benchmark > check_log_${SYSTEM}_${THIS_BENCH}.txt 2>&1
	cat check_log_${SYSTEM}_${THIS_BENCH}.txt
	check_completion
	echo 
    fi

done
}

build_only_tests()
{

for sys in ${BUILD_ONLY}; do

echo
echo "Running build-only test for ${sys}"
echo

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
echo "Use -w followed by full path to a working copy of the code to test a working copy"
echo "Use -r followed by a revision number to test the specified svn version e.g. -r 4300"
echo "Use -p followed by a svn path to check out a branch or tag e.g. -p https://repository.astro.ex.ac.uk/torus/tags/torus3.0.1"
echo "Use the -d option to run the daily tests (default)."
echo "Use the -s option to run the stable version tests."
echo "Use the -b option to run build tests only"
echo ""
}

########################################################################################################

# Default mode is daily test
export MODE=daily
export RETURN_CODE=0 
export TORUS_SVN_PATH=https://repository.astro.ex.ac.uk/torus/trunk/torus/
export TORUS_SVN_REVISION=
export TORUS_WORKING_COPY=none

# Parse command line arguments
while [ $# -gt 0 ]
do
    case "$1" in 
	-s) export MODE=stable;;
	-d) export MODE=daily;;
	-b) export MODE=build;;
	-r) shift
	    export TORUS_SVN_REVISION="-r $1";;
	-p) shift
	    export TORUS_SVN_PATH=$1;;
	-w) shift 
	    export TORUS_WORKING_COPY=$1
	    export MODE=workingcopy;;
	-h) print_help
	    exit;;
    esac
shift
done

# Set platform specific variables. We will assume that we are running on post-zen.
export PATH=/home/torustest/openmpi/bin:/home/torustest/bin:/usr/local/bin:${PATH}:/usr/bin
export TORUS_FITSLIBS="/home/torustest/cfitsio/lib"
export TORUS_FC="gfortran -g -fcheck=all"

case ${MODE} in 

    daily) export SYS_TO_TEST="gfortran ompiosx ompiosx-openmp"
           export BUILD_ONLY=""
	   export DEBUG_OPTS="yes"
	   echo -------------------------------------------------------------------
	   echo TORUS daily test suite started on `date`
	   echo -------------------------------------------------------------------
	   echo;;

    workingcopy) export SYS_TO_TEST="gfortran ompiosx ompiosx-openmp"
           export BUILD_ONLY=""
	   export DEBUG_OPTS="yes"
	   echo -------------------------------------------------------------------
	   echo TORUS working copy tests started on `date`
	   echo -------------------------------------------------------------------
	   echo;;

    build) export SYS_TO_TEST=" "
           export BUILD_ONLY="gfortran ompiosx ompiosx-openmp"
	   export DEBUG_OPTS="yes"
	   echo -------------------------------------------------------------------
	   echo TORUS build tests started on `date`
	   echo -------------------------------------------------------------------
	   echo;;

    stable) export SYS_TO_TEST="gfortran ompiosx ompiosx-openmp"
	    export BUILD_ONLY=""
            export DEBUG_OPTS="yes no"
	    echo -------------------------------------------------------------------
	    echo TORUS stable version tests started on `date`
	    echo -------------------------------------------------------------------
	    echo;;

    *)  echo "ERROR: unrecognised mode"
	exit 1;;
esac

for opt in ${DEBUG_OPTS}; do
    export USEDEBUGFLAGS=${opt}

# Set name of output directory
    case ${MODE} in 
	daily)       export TEST_DIR=/data/torustest/torus_daily_test;;
	workingcopy) export TEST_DIR=${TORUS_WORKING_COPY}/tests;;
	stable)      export TEST_DIR=${HOME}/torus_stable_version_tests/debug=${USEDEBUGFLAGS};;
	build)       export TEST_DIR=${HOME}/torus_build_tests;;
	    *)       echo "Unrecognised MODE"
	             exit 1;;
    esac

    export TORUS_DATA=${TEST_DIR}/torus/data

# Set up working dir and check out source code
    prepare_run

# Test build but don't run benchmarks
    build_only_tests

# Run benchmark tests
    run_torus_test_suite

done

rm ${TEST_DIR}/lock

echo -------------------------------------------------------------------
echo TORUS tests finished at `date`
echo -------------------------------------------------------------------

exit ${RETURN_CODE}
