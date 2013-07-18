#!/bin/ksh

mail_to="acreman@astro.ex.ac.uk th@astro.ex.ac.uk haworth@astro.ex.ac.uk claire@astro.ex.ac.uk"
#mail_to="acreman@astro.ex.ac.uk"

export BASE_DIR=/Users/acreman
export TORUS_TEST_DIR=${BASE_DIR}/SCRATCH/torus_daily_test
export LOG_FILE=${BASE_DIR}/torus_daily_test_log
export PATH=/Users/acreman/bin:${PATH}

if [[ -e ${TORUS_TEST_DIR}/lock ]]; then
  for user in ${mail_to}; do
     /sw/bin/mutt -s "Torus test suite: FAILED (did not complete)" ${user} < ${LOG_FILE}
    exit 1
  done
fi

# Process output from coverage analysis
echo "--------------------------" >> ${LOG_FILE}
echo "Processing coverage output" >> ${LOG_FILE}
echo "--------------------------" >> ${LOG_FILE}
cd ${TORUS_TEST_DIR}/benchmarks_ompiosx/build
process_gcov.sh >> ${LOG_FILE}

cd  ${TORUS_TEST_DIR}

attach_list=""

for file in svn_log.txt build_only_*/*/compile_log* benchmarks_*/*/compile_log* benchmarks_*/benchmarks/*/run_log* benchmarks_*/benchmarks/*/check_log*  benchmarks_*/benchmarks/*/tune_*.txt benchmarks_*/build/coverage_sorted.dat
 
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

suite_status="PASSED"

echo "Summary of test results: " > header
echo " " >> header

# Test for success of disc benchmark
num_success=`/usr/bin/grep "TORUS: Test successful"  benchmarks_ompiosx-openmp/benchmarks/disc/check_log_ompiosx_disc.txt | /usr/bin/wc -l`
num_success2=`/usr/bin/grep "TORUS: Test successful" benchmarks_gfortran/benchmarks/disc/check_log_gfortran_disc.txt | /usr/bin/wc -l`
num_success3=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompiosx/benchmarks/disc/check_log_ompiosx_disc.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 3 && ${num_success2} -eq 3  && ${num_success3} -eq 3 ]]; then
    echo "Disc benchmark successful" >> header 
else
    echo "!! Disc benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of molebench
num_success=`/usr/bin/grep "TORUS: Test successful"  benchmarks_ompiosx-openmp/benchmarks/molebench/check_log_ompiosx_molebench.txt | /usr/bin/wc -l`
num_success2=`/usr/bin/grep "TORUS: Test successful" benchmarks_gfortran/benchmarks/molebench/check_log_gfortran_molebench.txt | /usr/bin/wc -l`
num_success3=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompiosx/benchmarks/molebench/check_log_ompiosx_molebench.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 2 && ${num_success2} -eq 2 && ${num_success3} -eq 2 ]]; then
    echo "Molecular benchmark successful." >> header
else
    echo "!! Molecular benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of hydro benchmark
num_success=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompiosx/benchmarks/hydro/check_log_ompiosx_hydro.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 ]]; then
    echo "Hydro benchmark successful." >> header
else
    echo "!! Hydro benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of hII region benchmark
num_success=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompiosx-openmp/benchmarks/HII_region/check_log_ompiosx_hII.txt | /usr/bin/wc -l`
num_success2=`/usr/bin/grep "TORUS: Test successful" benchmarks_gfortran/benchmarks/HII_region/check_log_gfortran_hII.txt | /usr/bin/wc -l`
num_success3=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompiosx/benchmarks/HII_region/check_log_ompiosx_hII.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 && ${num_success2} -eq 1 && ${num_success3} -eq 1 ]]; then
    echo "HII region benchmark successful." >> header
else
    echo "!! HII region benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of domain decomposed hII region benchmark                                                                  
num_success=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompiosx/benchmarks/HII_regionMPI/check_log_ompiosx_hII_MPI.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 ]]; then
    echo "MPI HII region benchmark successful." >> header
else
    echo "!! MPI HII region benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of imaging benchmark
num_success=`/usr/bin/grep "Test Successful" benchmarks_ompiosx/benchmarks/cylinder_image_test/check_log_ompiosx_image.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 ]]; then
    echo "Image benchmark successful. " >> header
else
    echo "!! Image benchmark FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of gravity solver test
num_success=`/usr/bin/grep "Torus gravity solver test successful" benchmarks_ompiosx/benchmarks/gravtest/check_log_ompiosx_gravtest.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 ]]; then
    echo "Gravity test successful. " >> header
else
    echo "!! Gravity test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of 2D gravity solver test
num_success=`/usr/bin/grep "Torus gravity solver test successful" benchmarks_ompiosx/benchmarks/gravtest_2d/check_log_ompiosx_gravtest_2d.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 ]]; then
    echo "2D gravity test successful. " >> header
else
    echo "!! 2D gravity test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of nbody test
num_success=`/usr/bin/grep "Torus nbody test successful" benchmarks_gfortran/benchmarks/nbody/check_log_gfortran_nbody.txt | /usr/bin/wc -l`
num_success2=`/usr/bin/grep "Torus nbody test successful" benchmarks_ompiosx-openmp/benchmarks/nbody/check_log_ompiosx_nbody.txt | /usr/bin/wc -l`
num_success3=`/usr/bin/grep "Torus nbody test successful" benchmarks_ompiosx/benchmarks/nbody/check_log_ompiosx_nbody.txt | /usr/bin/wc -l`
if [[  ${num_success} -eq 1 && ${num_success2} -eq 1 && ${num_success3} -eq 1 ]]; then
    echo "N body test successful. " >> header
else
    echo "!! N body test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of SPH to grid test
num_success=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompiosx-openmp/benchmarks/sphbench/check_log_ompiosx_sphbench.txt | /usr/bin/wc -l`
num_success2=`/usr/bin/grep "TORUS: Test successful" benchmarks_gfortran/benchmarks/sphbench/check_log_gfortran_sphbench.txt | /usr/bin/wc -l`
num_success3=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompiosx/benchmarks/sphbench/check_log_ompiosx_sphbench.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 && ${num_success2} -eq 1 && ${num_success3} -eq 1 ]]; then
    echo "SPH to grid test successful." >> header
else
    echo "!! SPH to grid test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of restart test
num_success=`/usr/bin/grep "TORUS: Test successful"  benchmarks_ompiosx-openmp/benchmarks/restart/check_log_ompiosx_restart.txt | /usr/bin/wc -l`
num_success2=`/usr/bin/grep "TORUS: Test successful" benchmarks_gfortran/benchmarks/restart/check_log_gfortran_restart.txt | /usr/bin/wc -l`
num_success3=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompiosx/benchmarks/restart/check_log_ompiosx_restart.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 3 && ${num_success2} -eq 3  && ${num_success3} -eq 3 ]]; then
    echo "Restart test successful" >> header 
else
    echo "!! Restart test FAILED !!" >> header
    suite_status="FAILED"
fi

# Test for success of angular image test
num_success=`/usr/bin/grep -c "TORUS: Test successful"  benchmarks_ompiosx-openmp/benchmarks/angularImageTest/check_log_ompiosx_angularImageTest.txt`
num_success2=`/usr/bin/grep -c "TORUS: Test successful" benchmarks_gfortran/benchmarks/angularImageTest/check_log_gfortran_angularImageTest.txt`
num_success3=`/usr/bin/grep -c "TORUS: Test successful" benchmarks_ompiosx/benchmarks/angularImageTest/check_log_ompiosx_angularImageTest.txt`
if [[ ${num_success} -eq 1 && ${num_success2} -eq 1  && ${num_success3} -eq 1 ]]; then
    echo "Angular image test successful" >> header 
else
    echo "!! Angular image test FAILED !!" >> header
    suite_status="FAILED"
fi

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

# Set up the message body 
echo  >> header
mv ${LOG_FILE} ${TORUS_TEST_DIR}/torus_daily_test_log~
cat header ${TORUS_TEST_DIR}/torus_daily_test_log~ timings > ${TORUS_TEST_DIR}/torus_daily_test_log 
export LOG_FILE=${TORUS_TEST_DIR}/torus_daily_test_log

# Send mail 
gzip torus_test_output.tar
for user in ${mail_to}; do
    /sw/bin/mutt -s "Torus test suite: ${suite_status}" -a torus_test_output.tar.gz ${user} < ${LOG_FILE}
done

exit

