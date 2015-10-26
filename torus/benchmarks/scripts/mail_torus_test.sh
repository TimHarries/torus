#!/bin/ksh

# This script processes the output from a torus test run and reports whether the tests have passed
# D. Acreman

# Default mode is daily test
export MODE=daily
export BASE_DIR=/data/torustest
#mail_to="acreman@astro.ex.ac.uk th@astro.ex.ac.uk haworth@astro.ex.ac.uk claire@astro.ex.ac.uk"

# Parse command line arguments
while [ $# -gt 0 ]
do
    case "$1" in 
	-s) export MODE=stable;;
	-d) export MODE=daily;;
	-test) mail_to="acreman@astro.ex.ac.uk";;
    esac
shift
done

case ${MODE} in 

    daily) export TORUS_TEST_DIR=${BASE_DIR}/torus_daily_test
	export LOG_FILE=${BASE_DIR}/torus_daily_test_log
# This is probably run from cron so set up PATH
	export PATH=/home/torustest/bin:${PATH};;

    stable) export TORUS_TEST_DIR=${BASE_DIR}/torus_stable_version_tests/debug=yes
	export LOG_FILE=${BASE_DIR}/torus_stable_test_log;
# There may not be a log file already so create it so that it can be appended to
	touch ${LOG_FILE};;

    *)  echo "ERROR: unrecognised mode"
	exit 1;;
esac

# Check that we have a test directory
if [[ ! -d ${TORUS_TEST_DIR} ]]; then 
    echo "Test directory ${TORUS_TEST_DIR} not found"
    exit 1
fi

# Check that the test run completed OK 
if [[ ${MODE} == "daily" ]]; then
    if [[ -e ${TORUS_TEST_DIR}/lock ]]; then
# Indicate we are ready to mail out results
	echo "Torus test suite: FAILED (did not complete)" > ${BASE_DIR}/ready
#	for user in ${mail_to}; do
#	    /usr/bin/mail -s "Torus test suite: FAILED (did not complete)" ${user} < ${LOG_FILE}
#	done
# Bail out to avoid messing up the log file
	exit 1
    fi
fi

#
# Process output from coverage analysis. This is from the ompiosx build only.
#
cd ${TORUS_TEST_DIR}/benchmarks_ompiosx/build
echo "--------------------------" >> ${LOG_FILE}
echo "Processing coverage output" >> ${LOG_FILE}
echo "--------------------------" >> ${LOG_FILE}
process_gcov.sh >> ${LOG_FILE}
echo >> ${LOG_FILE}

#
# Make tar file with attachments (daily tests only)
#
cd  ${TORUS_TEST_DIR}
if [[ ${MODE} == "daily" ]]; then
    attach_list=""

    for file in svn_log.txt */*/compile_log* benchmarks_*/benchmarks/*/run_log* benchmarks_*/benchmarks/*/check_log*  benchmarks_*/benchmarks/*/tune_*.txt benchmarks_*/build/coverage_sorted.dat
 
    do
	if [[ ! -e ${file} ]]; then 
	    echo "${file} does not exist"
	else
	    file_size=`du -ks ${file} | awk '{print $1}'`
	    if [[ ${file_size} -gt 1000 ]]; then
		echo "Warning: file ${file} has size ${file_size} k" >> ${LOG_FILE}
		echo "This file will not be mailed." >> ${LOG_FILE}
	    else
		attach_list="${attach_list} ${file}"
	    fi
	fi
    done

    tar cf torus_test_output.tar ${attach_list}
    gzip torus_test_output.tar
fi

#
# Check each benchmark in turn. If any have failed then fail the suite.
#
suite_status="PASSED"

echo "Summary of test results: " > header
echo " " >> header

# Test for success of disc benchmark
num_success=`/bin/grep "TORUS: Test successful"  benchmarks_ompiosx-openmp/benchmarks/disc/check_log_ompiosx_disc.txt | /usr/bin/wc -l`
num_success2=`/bin/grep "TORUS: Test successful" benchmarks_gfortran/benchmarks/disc/check_log_gfortran_disc.txt | /usr/bin/wc -l`
num_success3=`/bin/grep "TORUS: Test successful" benchmarks_ompiosx/benchmarks/disc/check_log_ompiosx_disc.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 3 && ${num_success2} -eq 3  && ${num_success3} -eq 3 ]]; then
    echo "Disc benchmark successful" >> header 
else
    echo "!! Disc benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

if [[ ${MODE} == "stable" ]]; then
    # Test for success of 3D disc benchmark
    num_success=`/bin/grep -c "TORUS: Test successful"  benchmarks_ompiosx-openmp/benchmarks/disc_cylindrical/check_log_ompiosx_disc_cylindrical.txt`
    num_success2=`/bin/grep -c "TORUS: Test successful" benchmarks_gfortran/benchmarks/disc_cylindrical/check_log_gfortran_disc_cylindrical.txt`
    num_success3=`/bin/grep -c "TORUS: Test successful" benchmarks_ompiosx/benchmarks/disc_cylindrical/check_log_ompiosx_disc_cylindrical.txt`
# SEDs only, no image so look for 2 successful results
    if [[ ${num_success} -eq 2 && ${num_success2} -eq 2  && ${num_success3} -eq 2 ]]; then
	echo "3D Disc benchmark successful" >> header 
    else
	echo "!! 3D Disc benchmark FAILED !!" >> header
	suite_status="FAILED"
    fi
fi

# Test for success of molebench
num_success=`/bin/grep -c "TORUS: Test successful"  benchmarks_ompiosx-openmp/benchmarks/molebench/check_log_ompiosx_molebench.txt`
num_success2=`/bin/grep -c "TORUS: Test successful" benchmarks_gfortran/benchmarks/molebench/check_log_gfortran_molebench.txt`
num_success3=`/bin/grep -c "TORUS: Test successful" benchmarks_ompiosx/benchmarks/molebench/check_log_ompiosx_molebench.txt`
if [[ ${num_success} -eq 2 && ${num_success2} -eq 2 && ${num_success3} -eq 2 ]]; then
    echo "Molecular benchmark successful." >> header
else
    echo "!! Molecular benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of molecular mod restart
num_success=`/bin/grep -c "TORUS: Test successful"  benchmarks_ompiosx-openmp/benchmarks/molebench/restart/check_log_ompiosx_moleRestart.txt`
num_success2=`/bin/grep -c "TORUS: Test successful" benchmarks_gfortran/benchmarks/molebench/restart/check_log_gfortran_moleRestart.txt`
num_success3=`/bin/grep -c "TORUS: Test successful" benchmarks_ompiosx/benchmarks/molebench/restart/check_log_ompiosx_moleRestart.txt`
if [[ ${num_success} -eq 2 && ${num_success2} -eq 2 && ${num_success3} -eq 2 ]]; then
    echo "Molecular restart successful." >> header
else
    echo "!! Molecular restart FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of hydro benchmark
num_success=`/bin/grep "TORUS: Test successful" benchmarks_ompiosx/benchmarks/hydro/check_log_ompiosx_hydro.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 ]]; then
    echo "Hydro benchmark successful." >> header
else
    echo "!! Hydro benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of hII region benchmark
num_success=`/bin/grep "TORUS: Test successful" benchmarks_ompiosx-openmp/benchmarks/HII_region/check_log_ompiosx_hII.txt | /usr/bin/wc -l`
num_success2=`/bin/grep "TORUS: Test successful" benchmarks_gfortran/benchmarks/HII_region/check_log_gfortran_hII.txt | /usr/bin/wc -l`
num_success3=`/bin/grep "TORUS: Test successful" benchmarks_ompiosx/benchmarks/HII_region/check_log_ompiosx_hII.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 && ${num_success2} -eq 1 && ${num_success3} -eq 1 ]]; then
    echo "HII region benchmark successful." >> header
else
    echo "!! HII region benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of domain decomposed hII region benchmark                                                                  
num_success=`/bin/grep "TORUS: Test successful" benchmarks_ompiosx/benchmarks/HII_regionMPI/check_log_ompiosx_hII_MPI.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 ]]; then
    echo "MPI HII region benchmark successful." >> header
else
    echo "!! MPI HII region benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of imaging benchmark
num_success=`/bin/grep "Test Successful" benchmarks_ompiosx/benchmarks/cylinder_image_test/check_log_ompiosx_image.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 ]]; then
    echo "Image benchmark successful. " >> header
else
    echo "!! Image benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of gravity solver test
num_success=`/bin/grep "Torus gravity solver test successful" benchmarks_ompiosx/benchmarks/gravtest/check_log_ompiosx_gravtest.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 ]]; then
    echo "Gravity test successful. " >> header
else
    echo "!! Gravity test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of 2D gravity solver test
num_success=`/bin/grep "Torus gravity solver test successful" benchmarks_ompiosx/benchmarks/gravtest_2d/check_log_ompiosx_gravtest_2d.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 ]]; then
    echo "2D gravity test successful. " >> header
else
    echo "!! 2D gravity test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of nbody test
num_success=`/bin/grep "Torus nbody test successful" benchmarks_gfortran/benchmarks/nbody/check_log_gfortran_nbody.txt | /usr/bin/wc -l`
num_success2=`/bin/grep "Torus nbody test successful" benchmarks_ompiosx-openmp/benchmarks/nbody/check_log_ompiosx_nbody.txt | /usr/bin/wc -l`
num_success3=`/bin/grep "Torus nbody test successful" benchmarks_ompiosx/benchmarks/nbody/check_log_ompiosx_nbody.txt | /usr/bin/wc -l`
if [[  ${num_success} -eq 1 && ${num_success2} -eq 1 && ${num_success3} -eq 1 ]]; then
    echo "N body test successful. " >> header
else
    echo "!! N body test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of SPH to grid test
num_success=`/bin/grep "TORUS: Test successful" benchmarks_ompiosx-openmp/benchmarks/sphbench/check_log_ompiosx_sphbench.txt | /usr/bin/wc -l`
num_success2=`/bin/grep "TORUS: Test successful" benchmarks_gfortran/benchmarks/sphbench/check_log_gfortran_sphbench.txt | /usr/bin/wc -l`
num_success3=`/bin/grep "TORUS: Test successful" benchmarks_ompiosx/benchmarks/sphbench/check_log_ompiosx_sphbench.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 && ${num_success2} -eq 1 && ${num_success3} -eq 1 ]]; then
    echo "SPH to grid test successful." >> header
else
    echo "!! SPH to grid test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of SPH to grid test (binary dump with chemistry)
num_success=`/bin/grep -c "TORUS: Test successful" benchmarks_ompiosx-openmp/benchmarks/sphToGridBinary/check_log_ompiosx_sphToGridBinary.txt`
num_success2=`/bin/grep -c "TORUS: Test successful" benchmarks_gfortran/benchmarks/sphToGridBinary/check_log_gfortran_sphToGridBinary.txt`
num_success3=`/bin/grep -c "TORUS: Test successful" benchmarks_ompiosx/benchmarks/sphToGridBinary/check_log_ompiosx_sphToGridBinary.txt`
if [[ ${num_success} -eq 1 && ${num_success2} -eq 1 && ${num_success3} -eq 1 ]]; then
    echo "SPH to grid test with chemistry successful." >> header
else
    echo "!! SPH to grid test with chemistry FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of restart test
num_success=`/bin/grep -c "TORUS: Test successful"  benchmarks_ompiosx-openmp/benchmarks/restart/check_log_ompiosx_restart.txt`
num_success2=`/bin/grep -c "TORUS: Test successful" benchmarks_gfortran/benchmarks/restart/check_log_gfortran_restart.txt`
num_success3=`/bin/grep -c "TORUS: Test successful" benchmarks_ompiosx/benchmarks/restart/check_log_ompiosx_restart.txt`
if [[ ${num_success} -eq 3 && ${num_success2} -eq 3  && ${num_success3} -eq 3 ]]; then
    echo "Restart test successful" >> header 
else
    echo "!! Restart test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of angular image test
num_success=`/bin/grep -c "TORUS: Test successful"  benchmarks_ompiosx-openmp/benchmarks/angularImageTest/check_log_ompiosx_angularImageTest.txt`
num_success2=`/bin/grep -c "TORUS: Test successful" benchmarks_gfortran/benchmarks/angularImageTest/check_log_gfortran_angularImageTest.txt`
num_success3=`/bin/grep -c "TORUS: Test successful" benchmarks_ompiosx/benchmarks/angularImageTest/check_log_ompiosx_angularImageTest.txt`
if [[ ${num_success} -eq 1 && ${num_success2} -eq 1  && ${num_success3} -eq 1 ]]; then
    echo "Angular image test successful" >> header 
else
    echo "!! Angular image test FAILED !!" >> header
    suite_status="FAILED"
fi

# We won't attach the output now that it is available on post-zen
echo  >> header
echo "Output from these tests is on post-zen in ${TORUS_TEST_DIR}" >> header
echo  >> header

# Extract timing information
rm -f timings
echo >> timings
echo "--------------------------------" >> timings
echo "Timing information from tune.dat" >> timings
echo "--------------------------------" >> timings
echo >> timings

for tunefile in `find . -name 'tune*.txt'`; do
    echo ${tunefile} >> timings
    tail -1 ${tunefile} >> timings
    echo >> timings
done

# Send mail for daily test or write to terminal for other modes
if [[ ${MODE} == "daily" ]]; then
# Set up the message body 
    mv ${LOG_FILE} ${TORUS_TEST_DIR}/torus_daily_test_log~
    cat header ${TORUS_TEST_DIR}/torus_daily_test_log~ timings > ${TORUS_TEST_DIR}/torus_daily_test_log 
    export LOG_FILE=${TORUS_TEST_DIR}/torus_daily_test_log
#
    echo "Torus test suite: ${suite_status}" > ${BASE_DIR}/ready
#    for user in ${mail_to}; do
#	/usr/bin/mail -s "Torus test suite: ${suite_status}" -a torus_test_output.tar.gz ${user} < ${LOG_FILE}
#    done
else
    echo "Torus test suite: ${suite_status}"
    cat header
fi

exit

