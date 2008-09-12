#!/bin/ksh

export BASE_DIR=/Users/acreman
export TORUS_TEST_DIR=${BASE_DIR}/torus_daily_test
export LOG_FILE=${BASE_DIR}/torus_daily_test_log
export num_tests=2

cd  ${TORUS_TEST_DIR}

attach_list=""

for file in cvs_log.txt  build_ompi/compile_log_ompi.txt \
            run_ompi/run_log_ompi.txt run_ompi_molebench/run_log_ompi.txt
do
    file_size=`du -ks ${file} | awk '{print $1}'`
    if [[ ${file_size} -gt 1000 ]]; then
	echo "Warning: file ${file} has size ${file_size} k" >> ${LOG_FILE}
	echo "This file will not be mailed." >> ${LOG_FILE}
    else
	attach_list="${attach_list} -a ${file}"
    fi
done

# Test for success of disc benchmark
num_success=`/usr/bin/grep "TORUS: Test successful" ${LOG_FILE} | /usr/bin/wc -l `

if [[ ${num_success} -eq ${num_tests} ]]; then
    subject_line="Disc benchmark successful. "
else
    subject_line="Disc benchmark failed. "
fi

# Test for success of molebench
num_success=`/usr/bin/grep "Test passed" ${LOG_FILE} | /usr/bin/wc -l `

if [[ ${num_success} -eq 1 ]]; then
    subject_line2="Molecular benchmark successful. "
else
    subject_line2="Molecular benchmark failed. "
fi

# Test for success of cylindrical polar disc benchmark
num_success=`/usr/bin/grep "TORUS: Test successful" ${TORUS_TEST_DIR}/run_ompi_disc_cylindrical/check_log | /usr/bin/wc -l `

if [[ ${num_success} -eq ${num_tests} ]]; then
    subject_line3="Cylindrical polar disc benchmark successful. "
else
    subject_line3="Cylindrical polar disc benchmark failed. "
fi

# Join the output files together and get rid of the old log file
cat ${LOG_FILE} ${TORUS_TEST_DIR}/sphbench/run/check_log  ${TORUS_TEST_DIR}/run_ompi_disc_cylindrical/check_log > ${TORUS_TEST_DIR}/torus_daily_test_log 
rm ${LOG_FILE}
export LOG_FILE=${TORUS_TEST_DIR}/torus_daily_test_log

# Send mail 
mail_to="acreman@astro.ex.ac.uk th@astro.ex.ac.uk drundle@astro.ex.ac.uk N.J.Mayne@exeter.ac.uk"
for user in ${mail_to}; do
    /sw/bin/mutt -s "${subject_line} ${subject_line2} ${subject_line3}" ${attach_list} ${user} < ${LOG_FILE}
done

exit

