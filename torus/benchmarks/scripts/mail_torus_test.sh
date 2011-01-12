#!/bin/ksh

export BASE_DIR=/Users/acreman
export TORUS_TEST_DIR=${BASE_DIR}/SCRATCH/torus_daily_test
export LOG_FILE=${BASE_DIR}/torus_daily_test_log

cd  ${TORUS_TEST_DIR}

attach_list=""

for file in cvs_log.txt build_only_*/*/compile_log* benchmarks_*/*/compile_log* benchmarks_*/benchmarks/*/run_log* benchmarks_*/benchmarks/*/check_log*  benchmarks_*/benchmarks/*/tune_*.txt
 
do
    if [[ ! -e ${file} ]]; then 
	echo "${file} does not exist"
    else
	file_size=`du -ks ${file} | awk '{print $1}'`
	if [[ ${file_size} -gt 1000 ]]; then
	    echo "Warning: file ${file} has size ${file_size} k" >> ${LOG_FILE}
	    echo "This file will not be mailed." >> ${LOG_FILE}
	else
	    attach_list="${attach_list} -a ${file}"
	fi
    fi
done

echo  >> ${LOG_FILE}
echo "Checking if memory has been freed OK ..." >> ${LOG_FILE}
echo  >> ${LOG_FILE}
for file in benchmarks_ompi/benchmarks/*/run_log*
do
    num_not_freed=`/usr/bin/grep "Remaining memory:" ${file} | /usr/bin/wc -l`
    if [[ ${num_not_freed} -gt 0 ]]; then
	this_file=`basename ${file}`
	echo "WARNING: Unfreed memory in ${this_file}" >> ${LOG_FILE}
    fi
done
echo >> ${LOG_FILE}
echo "Finished checking memory" >> ${LOG_FILE}
echo  >> ${LOG_FILE}

subject_line=" " 

# Test for success of disc benchmark
num_success=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompi/benchmarks/disc/check_log_disc.txt | /usr/bin/wc -l`
num_success2=`/usr/bin/grep "TORUS: Test successful" benchmarks_gfortran/benchmarks/disc/check_log_disc.txt | /usr/bin/wc -l`
num_success3=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompiosx/benchmarks/disc/check_log_disc.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 2 && ${num_success2} -eq 2  && ${num_success3} -eq 2 ]]; then
    subject_line="${subject_line} Disc benchmark successful."
else
    subject_line="${subject_line} Disc benchmark failed."
fi

# Test for success of molebench
num_success=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompi/benchmarks/molebench/check_log_molebench.txt | /usr/bin/wc -l`
num_success2=`/usr/bin/grep "TORUS: Test successful" benchmarks_gfortran/benchmarks/molebench/check_log_molebench.txt | /usr/bin/wc -l`
num_success3=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompiosx/benchmarks/molebench/check_log_molebench.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 && ${num_success2} -eq 1 && ${num_success3} -eq 1 ]]; then
    subject_line="${subject_line} Molecular benchmark successful."
else
    subject_line="${subject_line} Molecular benchmark failed."
fi

# Test for success of hydro benchmark
num_success=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompi/benchmarks/hydro/check_log_hydro.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 ]]; then
    subject_line="${subject_line} Hydro benchmark successful. "
else
    subject_line="${subject_line} Hydro benchmark failed. "
fi

# Test for success of hII region benchmark
num_success=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompi/benchmarks/HII_region/check_log_hII.txt | /usr/bin/wc -l`
num_success2=`/usr/bin/grep "TORUS: Test successful" benchmarks_gfortran/benchmarks/HII_region/check_log_hII.txt | /usr/bin/wc -l`
num_success3=`/usr/bin/grep "TORUS: Test successful" benchmarks_ompiosx/benchmarks/HII_region/check_log_hII.txt | /usr/bin/wc -l`
if [[ ${num_success} -eq 1 && ${num_success2} -eq 1 && ${num_success3} -eq 1 ]]; then
    subject_line="${subject_line} HII region benchmark successful. "
else
    subject_line="${subject_line} HII region benchmark failed. "
fi

# Move log file
mv ${LOG_FILE} ${TORUS_TEST_DIR}/torus_daily_test_log 
export LOG_FILE=${TORUS_TEST_DIR}/torus_daily_test_log

# Send mail 
mail_to="acreman@astro.ex.ac.uk th@astro.ex.ac.uk N.J.Mayne@exeter.ac.uk tjh202@exeter.ac.uk"
for user in ${mail_to}; do
    /sw/bin/mutt -s "${subject_line}" ${attach_list} ${user} < ${LOG_FILE}
done

exit

