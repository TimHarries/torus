#!/bin/ksh
# Build the Torus executable
make_build()
{
mkdir build
cd    build 

echo "Building Torus for ${THISCONFIG} configuration"
log_file=compile_log_${THISCONFIG}.txt
rsync --exclude=data -a ${TEST_DIR}/torus/* .
rsync -a ${TEST_DIR}/torus/.svn .

/usr/bin/make depends > ${log_file} 2>&1 
/usr/bin/make debug=${USEDEBUGFLAGS} mpi=${USEMPI} openmp=${USEOPENMP} coverage=${USEGCOV} >> ${log_file} 2>&1

if [[ $? -eq 0 ]]; then
# Count number of warnings. Subtract 2 because there are always warnings
# about include files from the make depends step (run twice).
    num_warn=`grep -i -c warning ${log_file} | awk '{print $1 - 2}`
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

log_file=run_log_${THISCONFIG}_${THIS_BENCH}.txt
export TORUS_JOB_DIR=./
ln -s ${WORKING_DIR}/build/torus.${SYSTEM} .

case ${THISCONFIG} in
    openmp) ./torus.${SYSTEM} > ${log_file} 2>&1 ;;
    mpi)    mpirun -np ${NUM_MPI_PROC} torus.${SYSTEM} > ${log_file} 2>&1 ;;
    hybrid) mpirun -np ${NUM_MPI_PROC} torus.${SYSTEM} > ${log_file} 2>&1 ;;
    *) echo "Unrecognised configuration. Skipping this test.";;
esac

mv tune.dat tune_${THISCONFIG}_${THIS_BENCH}.txt 
}

run_bench()
{
cd ${WORKING_DIR}/benchmarks/${THIS_BENCH}
ln -s ${WORKING_DIR}/build/torus.${SYSTEM} .
log_file=run_log_${THISCONFIG}_${THIS_BENCH}.txt
export TORUS_JOB_DIR=./

case ${THISCONFIG} in
    openmp) ./torus.${SYSTEM} > ${log_file} 2>&1 ;;
    mpi)    mpirun -np ${NUM_MPI_PROC} torus.${SYSTEM} > ${log_file} 2>&1 ;;
    hybrid) mpirun -np ${NUM_MPI_PROC} torus.${SYSTEM} > ${log_file} 2>&1 ;;
    *) echo "Unrecognised configuration. Skipping this test.";;
esac

mv tune.dat tune_${THISCONFIG}_${THIS_BENCH}.txt 
}

# Run Torus with a domain decomposed grid.  
run_dom_decomp()
{
cd ${WORKING_DIR}/benchmarks/${THIS_BENCH}
ln -s ${WORKING_DIR}/build/torus.${SYSTEM} . 
mpirun -np $1 torus.${SYSTEM} > run_log_${THIS_BENCH}.txt 2>&1
mv tune.dat tune_${THIS_BENCH}.txt
# The check_completion function expects the run log to be tagged with the configuration
ln -s run_log_${THIS_BENCH}.txt run_log_${THISCONFIG}_${THIS_BENCH}.txt
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

check_rainbow()
{
echo "Compiling checkrainbow"
${TORUS_FC} -o checkrainbow checkrainbow.f90
./checkrainbow
}

# Check that Torus completed OK. For some tests the results file can be written 
# even if Torus has bugged out. 
check_completion()
{
    grep "Torus completed" run_log_${THISCONFIG}_${THIS_BENCH}.txt > /dev/null
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
    echo "Removing old ${TEST_DIR}"
    rm -rf ${TEST_DIR}
fi

if [[ ${TORUS_WORKING_COPY} == none ]]; then 
    mkdir -p ${TEST_DIR}
    cd ${TEST_DIR}
    touch ${LOCKFILE}
    echo Checking out torus from SVN archive using: ${TORUS_SVN_REVISION} ${TORUS_SVN_PATH}
    /usr/bin/svn checkout ${TORUS_SVN_REVISION} ${TORUS_SVN_PATH} torus > svn_log.txt 2>&1 
    grep "Checked out revision" svn_log.txt
else
    if [[ -e ${TORUS_WORKING_COPY}/torus/torusMainV2.F90 ]]; then
	echo "Taking source code from working copy in ${TORUS_WORKING_COPY}"
	mkdir -p ${TEST_DIR}
	cd ${TEST_DIR}
	touch ${LOCKFILE}
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

for config in ${CONFIG_TO_TEST}; do

# Details of how to run each configuration are set here
    export THISCONFIG=${config}
    export USEGCOV=no
    if [[ ${config} == "hybrid" ]]; then
	export USEOPENMP=yes
	export USEMPI=yes
	# OpenMPI will bind to core by default so we'll use the 2 OpenMP threads available from hyperthreading.
	# Binding can be disabled by giving mpirun the '--bind-to none' option if more OpenMP threads are required.
	export NUM_MPI_PROC=4
	export OMP_NUM_THREADS=2
    elif [[ ${config} == "mpi" ]]; then
	export USEOPENMP=no
	export USEMPI=yes
	export NUM_MPI_PROC=8
	export USEGCOV=yes
    elif [[ ${config} == "openmp" ]]; then
	export USEMPI=no
	export USEOPENMP=yes
	export OMP_NUM_THREADS=8
    else
	echo "${config} not recognised. Aborting."
	exit 1
    fi

    echo
    echo "=============================================="
    echo "Running tests for ${THISCONFIG} configuration"
    echo "=============================================="
    echo

    export WORKING_DIR=${TEST_DIR}/benchmarks_${config}
    mkdir ${WORKING_DIR}
    cd    ${WORKING_DIR} 
    cp -r ${TEST_DIR}/torus/benchmarks . 

# Build code
    make_build

# Check if we have an executable. If not proceed to the next configuration
    if [[ -x build/torus.${SYSTEM} ]]; then
	echo "Found executable."
    else
	echo "Executable not found. Skipping ${THISCONFIG} configuration"
	continue
    fi

# Build comparison code for disc benchmark. Be careful moving this as it changes the cwd. 
    make_comparespec

# Run hydro benchmark
    case ${config} in
	mpi)  echo "Running hydro benchmark"
	    export THIS_BENCH=hydro
	    run_dom_decomp 3
	    check_hydro  > check_log_${config}_hydro.txt 2>&1 
	    cat check_log_${config}_hydro.txt
	    check_completion
	    echo ;;
	*) echo "Hydro benchmark only runs with mpi configuration. Skipping"
	    echo ;;
    esac

# Domain decomposed Lexington benchmark
    case ${config} in
	mpi)  echo "Running domain decomposed HII region benchmark"
	    export THIS_BENCH=HII_regionMPI
	    run_dom_decomp 9
	    check_hII > check_log_${config}_hII_MPI.txt 2>&1 
	    cat check_log_${config}_hII_MPI.txt
	    check_completion
	    echo ;;
	*) echo "Domain decomposed HII region only runs with mpi configuration. Skipping"
	    echo ;;
    esac

# Imaging test
    case ${config} in
	mpi)  echo "Running imaging benchmark"
	    export THIS_BENCH=cylinder_image_test
	    run_dom_decomp 9
	    check_image > check_log_${config}_image.txt 2>&1 
	    cat check_log_${config}_image.txt
            check_completion
	    echo ;;
	*) echo "Imaging benchmark only runs with mpi configuration. Skipping"
	    echo ;;
    esac

# Gravity solver test
    case ${config} in
	mpi)  echo "Running 3D gravity solver test"
	    export THIS_BENCH=gravtest
	    run_dom_decomp 9
	    check_it > check_log_${config}_gravtest.txt 2>&1 
	    tail check_log_${config}_gravtest.txt
	    check_completion
	    echo ;;
	*) echo "Gravity solver test only runs with mpi configuration. Skipping"
	    echo ;;
    esac

# Gravity solver test in 2D
    case ${config} in
	mpi)  echo "Running 2D gravity solver test"
	    export THIS_BENCH=gravtest_2d
	    run_dom_decomp 5
	    check_it > check_log_${config}_gravtest_2d.txt 2>&1 
	    tail check_log_${config}_gravtest_2d.txt
	    check_completion
	    echo ;;
	*) echo "2D Gravity solver test only runs with mpi configuration. Skipping"
	    echo ;;
    esac

# rainbow test
    echo "Running rainbow test"
    export THIS_BENCH=rainbow
    run_bench
    check_rainbow > check_log_${config}_${THIS_BENCH}.txt 2>&1 
    cat check_log_${config}_${THIS_BENCH}.txt
    check_completion
    echo


# N body test
    echo "Running N body test"
    export THIS_BENCH=nbody
    run_bench
    check_it > check_log_${config}_${THIS_BENCH}.txt 2>&1 
    cat check_log_${config}_${THIS_BENCH}.txt
    check_completion
    echo

# Run 2D disc
    echo "Running disc benchmark"
    export THIS_BENCH=disc
    run_bench 
    check_benchmark > check_log_${config}_${THIS_BENCH}.txt 2>&1 
    cat check_log_${config}_${THIS_BENCH}.txt
    check_completion
    echo

    echo "Running HII region benchmark"
    export THIS_BENCH=HII_region
    run_bench
    check_hII > check_log_${config}_hII.txt 2>&1 
    cat check_log_${config}_hII.txt
    check_completion
    echo

    echo "Running molecular benchmark"
    export THIS_BENCH=molebench 
    run_bench
    check_molebench > check_log_${config}_${THIS_BENCH}.txt 2>&1 
    tail -20 check_log_${config}_${THIS_BENCH}.txt # Lots of output so tail this file
    check_completion
    echo

    echo "Running molecular restart test"
    export THIS_BENCH=moleRestart
    run_molecularRestart
    check_molebench > check_log_${config}_${THIS_BENCH}.txt 2>&1
    ./check_nrays.sh run_log_${config}_${THIS_BENCH}.txt >> check_log_${config}_${THIS_BENCH}.txt 2>&1 
    tail -20 check_log_${config}_${THIS_BENCH}.txt # Lots of output so tail this file
    check_completion
    echo

    echo "Running SPH to grid test"
    export THIS_BENCH=sphbench
    setup_sphbench
    run_bench
    ./checkSphToGrid.pl run_log_${config}_${THIS_BENCH}.txt > check_log_${config}_${THIS_BENCH}.txt 2>&1 
    cat check_log_${config}_${THIS_BENCH}.txt 2>&1
    check_completion
    echo

    echo "Running SPH to grid test (binary dump with chemistry)"
    export THIS_BENCH=sphToGridBinary
    cd ${WORKING_DIR}/benchmarks/sphToGridBinary
    ln -s ${HOME}/torus_dev/forTestSuite/SQA0321
    run_bench
    ./checkSphToGridChem.pl run_log_${config}_${THIS_BENCH}.txt > check_log_${config}_${THIS_BENCH}.txt 2>&1 
    cat check_log_${config}_${THIS_BENCH}.txt 2>&1
    check_completion
    echo

    echo "Running restart test"
    export THIS_BENCH=restart
    cd ${WORKING_DIR}/benchmarks/restart
    ln -s ../disc/lucy_grid_tmp.dat
    run_bench
    check_benchmark > check_log_${config}_${THIS_BENCH}.txt 2>&1 
    cat check_log_${config}_${THIS_BENCH}.txt
    check_completion
    echo

    echo "Running angular imaging test"
    export THIS_BENCH=angularImageTest
    run_bench
    check_angImg > check_log_${config}_${THIS_BENCH}.txt 2>&1 
    cat check_log_${config}_${THIS_BENCH}.txt
    check_completion
    echo


# Only run in stable version tests as this is slow
    if [[ ${MODE} == stable ]]; then
	echo "Running cylindrical polar disc benchmark"
	export THIS_BENCH=disc_cylindrical
	run_bench
	check_benchmark > check_log_${config}_${THIS_BENCH}.txt 2>&1
	cat check_log_${config}_${THIS_BENCH}.txt
	check_completion
	echo 
    fi

done
}

build_only_tests()
{
for config in ${BUILD_ONLY}; do

    export THISCONFIG=${config}
    export USEGCOV=no
    if [[ ${config} == "hybrid" ]]; then
	export USEMPI=yes
	export USEOPENMP=yes
    elif [[ ${config} == "mpi" ]]; then
	export USEMPI=yes
	export USEOPENMP=no
    elif [[ ${config} == "openmp" ]]; then
	export USEMPI=no
	export USEOPENMP=yes
    else
	echo "${config} not recognised. Aborting."
	exit 1
    fi

    export WORKING_DIR=${TEST_DIR}/build_only_${config}
    mkdir ${WORKING_DIR}
    cd    ${WORKING_DIR} 

    make_build
    echo " " 

done
}

# Process the output from gcov coverage analysis to work out what fraction of Torus has been 
# exercised by the test suite
process_gcov()
{
cd ${TEST_DIR}/benchmarks_mpi/build

echo "--------------------------"
echo "Processing coverage output"
echo "--------------------------"

rm -f coverage.dat coverage_sorted.dat

for file in *90
do
    gcov ${file} >> coverage.dat 2> /dev/null
done

# Trailing space on grep is to cope with filenames which contain the string File
grep 'File ' -A1 coverage.dat | tr '\n' ' ' | tr '-' '\n' | tr ':' ' ' >> coverage_sorted.dat
mv coverage_sorted.dat coverage_sorted.dat~
sort -r -g --key=5 < coverage_sorted.dat~ | crush > coverage_sorted.dat
export num_prof_lines=`awk '{sum += $7*$5/100} END {print sum}' $1 < coverage_sorted.dat`
export num_lines=`awk '{sum += $7} END {print sum}' $1 < coverage_sorted.dat`
export prof_frac=`echo ${num_prof_lines} ${num_lines} | awk '{print $1/$2}`

echo "Total number of lines: ${num_lines}"
echo "Fraction profiled: ${prof_frac}"
echo

}

# Process the timing information in tune.dat
process_timing()
{
cd ${TEST_DIR}

echo
echo "--------------------------------"
echo "Timing information from tune.dat"
echo "--------------------------------"
echo

for tunefile in `find . -name 'tune*.txt'`; do
    echo ${tunefile}
    tail -1 ${tunefile}
    echo
done
}

check_results()
{
cd ${TEST_DIR}
    
suite_status="PASSED"

echo "Summary of test results: " > header
echo " " >> header

grepper="/bin/grep -c -i"

# Test for success of disc benchmark
num_success=`${grepper} "TORUS: Test successful"  benchmarks_hybrid/benchmarks/disc/check_log_hybrid_disc.txt`
num_success2=`${grepper} "TORUS: Test successful" benchmarks_openmp/benchmarks/disc/check_log_openmp_disc.txt`
num_success3=`${grepper} "TORUS: Test successful" benchmarks_mpi/benchmarks/disc/check_log_mpi_disc.txt`
if [[ ${num_success} -eq 3 && ${num_success2} -eq 3  && ${num_success3} -eq 3 ]]; then
    echo "Disc benchmark successful" >> header 
else
    echo "!! Disc benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

if [[ ${MODE} == "stable" ]]; then
    # Test for success of 3D disc benchmark
    num_success=`${grepper} "TORUS: Test successful"  benchmarks_hybrid/benchmarks/disc_cylindrical/check_log_hybrid_disc_cylindrical.txt`
    num_success2=`${grepper} "TORUS: Test successful" benchmarks_openmp/benchmarks/disc_cylindrical/check_log_openmp_disc_cylindrical.txt`
    num_success3=`${grepper} "TORUS: Test successful" benchmarks_mpi/benchmarks/disc_cylindrical/check_log_mpi_disc_cylindrical.txt`
# SEDs only, no image so look for 2 successful results
    if [[ ${num_success} -eq 2 && ${num_success2} -eq 2  && ${num_success3} -eq 2 ]]; then
	echo "3D Disc benchmark successful" >> header 
    else
	echo "!! 3D Disc benchmark FAILED !!" >> header
	suite_status="FAILED"
    fi
fi

# Test for success of molebench
num_success=`${grepper} "TORUS: Test successful"  benchmarks_hybrid/benchmarks/molebench/check_log_hybrid_molebench.txt`
num_success2=`${grepper} "TORUS: Test successful" benchmarks_openmp/benchmarks/molebench/check_log_openmp_molebench.txt`
num_success3=`${grepper} "TORUS: Test successful" benchmarks_mpi/benchmarks/molebench/check_log_mpi_molebench.txt`
if [[ ${num_success} -eq 2 && ${num_success2} -eq 2 && ${num_success3} -eq 2 ]]; then
    echo "Molecular benchmark successful." >> header
else
    echo "!! Molecular benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of molecular mod restart
num_success=`${grepper} "TORUS: Test successful"  benchmarks_hybrid/benchmarks/molebench/restart/check_log_hybrid_moleRestart.txt`
num_success2=`${grepper} "TORUS: Test successful" benchmarks_openmp/benchmarks/molebench/restart/check_log_openmp_moleRestart.txt`
num_success3=`${grepper} "TORUS: Test successful" benchmarks_mpi/benchmarks/molebench/restart/check_log_mpi_moleRestart.txt`
if [[ ${num_success} -eq 3 && ${num_success2} -eq 3 && ${num_success3} -eq 3 ]]; then
    echo "Molecular restart successful." >> header
else
    echo "!! Molecular restart FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of hydro benchmark
num_success=`${grepper} "TORUS: Test successful" benchmarks_mpi/benchmarks/hydro/check_log_mpi_hydro.txt`
if [[ ${num_success} -eq 1 ]]; then
    echo "Hydro benchmark successful." >> header
else
    echo "!! Hydro benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of hII region benchmark
num_success=`${grepper} "TORUS: Test successful" benchmarks_hybrid/benchmarks/HII_region/check_log_hybrid_hII.txt`
num_success2=`${grepper} "TORUS: Test successful" benchmarks_openmp/benchmarks/HII_region/check_log_openmp_hII.txt`
num_success3=`${grepper} "TORUS: Test successful" benchmarks_mpi/benchmarks/HII_region/check_log_mpi_hII.txt`
if [[ ${num_success} -eq 1 && ${num_success2} -eq 1 && ${num_success3} -eq 1 ]]; then
    echo "HII region benchmark successful." >> header
else
    echo "!! HII region benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of domain decomposed hII region benchmark                                                                  
num_success=`${grepper} "TORUS: Test successful" benchmarks_mpi/benchmarks/HII_regionMPI/check_log_mpi_hII_MPI.txt`
if [[ ${num_success} -eq 1 ]]; then
    echo "MPI HII region benchmark successful." >> header
else
    echo "!! MPI HII region benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of imaging benchmark
num_success=`${grepper} "Test Successful" benchmarks_mpi/benchmarks/cylinder_image_test/check_log_mpi_image.txt`
if [[ ${num_success} -eq 1 ]]; then
    echo "Image benchmark successful. " >> header
else
    echo "!! Image benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of gravity solver test
num_success=`${grepper} "Torus gravity solver test successful" benchmarks_mpi/benchmarks/gravtest/check_log_mpi_gravtest.txt`
if [[ ${num_success} -eq 1 ]]; then
    echo "Gravity test successful. " >> header
else
    echo "!! Gravity test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of 2D gravity solver test
num_success=`${grepper} "Torus gravity solver test successful" benchmarks_mpi/benchmarks/gravtest_2d/check_log_mpi_gravtest_2d.txt`
if [[ ${num_success} -eq 1 ]]; then
    echo "2D gravity test successful. " >> header
else
    echo "!! 2D gravity test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of nbody test
num_success=`${grepper} "Torus nbody test successful" benchmarks_openmp/benchmarks/nbody/check_log_openmp_nbody.txt`
num_success2=`${grepper} "Torus nbody test successful" benchmarks_hybrid/benchmarks/nbody/check_log_hybrid_nbody.txt`
num_success3=`${grepper} "Torus nbody test successful" benchmarks_mpi/benchmarks/nbody/check_log_mpi_nbody.txt`
if [[  ${num_success} -eq 1 && ${num_success2} -eq 1 && ${num_success3} -eq 1 ]]; then
    echo "N body test successful. " >> header
else
    echo "!! N body test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of SPH to grid test
num_success=`${grepper} "TORUS: Test successful" benchmarks_hybrid/benchmarks/sphbench/check_log_hybrid_sphbench.txt`
num_success2=`${grepper} "TORUS: Test successful" benchmarks_openmp/benchmarks/sphbench/check_log_openmp_sphbench.txt`
num_success3=`${grepper} "TORUS: Test successful" benchmarks_mpi/benchmarks/sphbench/check_log_mpi_sphbench.txt`
if [[ ${num_success} -eq 1 && ${num_success2} -eq 1 && ${num_success3} -eq 1 ]]; then
    echo "SPH to grid test successful." >> header
else
    echo "!! SPH to grid test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of SPH to grid test (binary dump with chemistry)
num_success=`${grepper} "TORUS: Test successful" benchmarks_hybrid/benchmarks/sphToGridBinary/check_log_hybrid_sphToGridBinary.txt`
num_success2=`${grepper} "TORUS: Test successful" benchmarks_openmp/benchmarks/sphToGridBinary/check_log_openmp_sphToGridBinary.txt`
num_success3=`${grepper} "TORUS: Test successful" benchmarks_mpi/benchmarks/sphToGridBinary/check_log_mpi_sphToGridBinary.txt`
if [[ ${num_success} -eq 1 && ${num_success2} -eq 1 && ${num_success3} -eq 1 ]]; then
    echo "SPH to grid test with chemistry successful." >> header
else
    echo "!! SPH to grid test with chemistry FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of restart test
num_success=`${grepper} "TORUS: Test successful"  benchmarks_hybrid/benchmarks/restart/check_log_hybrid_restart.txt`
num_success2=`${grepper} "TORUS: Test successful" benchmarks_openmp/benchmarks/restart/check_log_openmp_restart.txt`
num_success3=`${grepper} "TORUS: Test successful" benchmarks_mpi/benchmarks/restart/check_log_mpi_restart.txt`
if [[ ${num_success} -eq 3 && ${num_success2} -eq 3  && ${num_success3} -eq 3 ]]; then
    echo "Restart test successful" >> header 
else
    echo "!! Restart test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of angular image test
num_success=`${grepper} "TORUS: Test successful"  benchmarks_hybrid/benchmarks/angularImageTest/check_log_hybrid_angularImageTest.txt`
num_success2=`${grepper} "TORUS: Test successful" benchmarks_openmp/benchmarks/angularImageTest/check_log_openmp_angularImageTest.txt`
num_success3=`${grepper} "TORUS: Test successful" benchmarks_mpi/benchmarks/angularImageTest/check_log_mpi_angularImageTest.txt`
if [[ ${num_success} -eq 1 && ${num_success2} -eq 1  && ${num_success3} -eq 1 ]]; then
    echo "Angular image test successful" >> header 
else
    echo "!! Angular image test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of rainbow test
num_success=`${grepper} "TORUS: Test successful"  benchmarks_hybrid/benchmarks/rainbow/check_log_hybrid_rainbow.txt`
num_success2=`${grepper} "TORUS: Test successful" benchmarks_openmp/benchmarks/rainbow/check_log_openmp_rainbow.txt`
num_success3=`${grepper} "TORUS: Test successful" benchmarks_mpi/benchmarks/rainbow/check_log_mpi_rainbow.txt`
if [[ ${num_success} -eq 1 && ${num_success2} -eq 1  && ${num_success3} -eq 1 ]]; then
    echo "Rainbow test successful" >> header 
else
    echo "!! Rainbow test FAILED !!" >> header
    suite_status="FAILED"
fi

# We won't attach the output now that it is available on post-zen
echo  >> header
echo "Output from these tests is on post-zen in ${TEST_DIR}" >> header
echo  >> header

# Send mail for daily test or write to terminal for other modes
if [[ ${MODE} == "daily" ]]; then
    mail_to="acreman@astro.ex.ac.uk th@astro.ex.ac.uk aali@astro.ex.ac.uk tdouglas@astro.ex.ac.uk t.haworth@imperial.ac.uk fjmw201@exeter.ac.uk"
# Set up the message body 
    cat header ${TORUS_DAILY_TEST_LOG} > /home/torustest/torus_daily_test_email
    for user in ${mail_to}; do
        /usr/bin/mail -s "Torus test suite: ${suite_status}" ${user} < /home/torustest/torus_daily_test_email
    done
else
    echo "Torus test suite: ${suite_status}"
    cat header
fi
  
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
export SYSTEM=testsuite
export TORUS_RUNNING_LOG=${HOME}/testsuite.log
export TORUS_DAILY_TEST_LOG=/home/torustest/torus_daily_test_log

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

# Check for a lock file from an existing run and bail out if one exists
case ${MODE} in 
    daily)       export LOCKFILE=${HOME}/torus_daily_test/lock;;
    workingcopy) export LOCKFILE=${TORUS_WORKING_COPY}/tests/lock;;
    build)       export LOCKFILE=${HOME}/torus_build_tests/lock;;
    stable)      export LOCKFILE=${HOME}/torus_stable_version_tests/lock;;
    *)           echo "ERROR: unrecognised mode"
	         exit 1;;
esac

if [[ -e ${LOCKFILE} ]]; then
    if [[ ${MODE} == daily ]]; then
	echo `date` "Found lock file. Aborting" >> ${TORUS_RUNNING_LOG}
	exit 1
    else
	echo "Found lock file ${LOCKFILE}. Aborting"
    fi
fi


# Set platform specific variables. We will assume that we are running on post-zen.
export PATH=/home/torustest/openmpi/bin:/home/torustest/bin:/usr/local/bin:${PATH}:/usr/bin
export TORUS_FITSLIBS="/home/torustest/cfitsio/lib"
export TORUS_FC="gfortran -g -fcheck=all"

case ${MODE} in 

    daily) export CONFIG_TO_TEST="openmp mpi hybrid"
           export BUILD_ONLY=""
	   export DEBUG_OPTS="yes"
	   echo `date` "Test suite started" >> ${TORUS_RUNNING_LOG} 
	   echo -------------------------------------------------------------------
	   echo TORUS daily test suite started on `date`
	   echo -------------------------------------------------------------------
	   echo;;

    workingcopy) export CONFIG_TO_TEST="openmp mpi hybrid"
           export BUILD_ONLY=""
	   export DEBUG_OPTS="yes"
	   echo -------------------------------------------------------------------
	   echo TORUS working copy tests started on `date`
	   echo -------------------------------------------------------------------
	   echo;;

    build) export CONFIG_TO_TEST=" "
           export BUILD_ONLY="openmp mpi hybrid"
	   export DEBUG_OPTS="yes"
	   echo -------------------------------------------------------------------
	   echo TORUS build tests started on `date`
	   echo -------------------------------------------------------------------
	   echo;;

    stable) export CONFIG_TO_TEST="openmp mpi hybrid"
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
	daily)       export TEST_DIR=${HOME}/torus_daily_test;;
	workingcopy) export TEST_DIR=${TORUS_WORKING_COPY}/tests;;
	stable)      export TEST_DIR=${HOME}/torus_stable_version_tests/debug=${USEDEBUGFLAGS};;
	build)       export TEST_DIR=${HOME}/torus_build_tests;;
	    *)       echo "Unrecognised MODE"
	             exit 1;;
    esac

    export TORUS_DATA=${TEST_DIR}/torus/data

# Set up working directory and check out source code
    prepare_run

# Test build but don't run benchmarks
    build_only_tests

# Run benchmark tests
    run_torus_test_suite

    if [[ ${MODE} != build ]]; then 
# Process results from coverage analysis
	process_gcov
# Report timing information
	process_timing
    fi

done

echo -------------------------------------------------------------------
echo TORUS test suite finished at `date`
echo -------------------------------------------------------------------

# Check results
if [[ ${MODE} == daily || ${MODE} == workingcopy ]]; then
    check_results
fi

rm ${LOCKFILE}
rm ${TORUS_DAILY_TEST_LOG}

echo `date` "Test suite finished" >> ${TORUS_RUNNING_LOG}

exit ${RETURN_CODE}
