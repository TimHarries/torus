#!/bin/ksh

export BASE_DIR=/Users/acreman
export TORUS_TEST_DIR=${BASE_DIR}/torus_daily_test
export LOG_FILE=${BASE_DIR}/torus_daily_test_log
export MAIL_DOMAIN="astro.ex.ac.uk"
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

mail_to="acreman th drundle"

# Test for success of disc benchmark
num_success=`/usr/bin/grep "TORUS: Test successful" ${LOG_FILE} | /usr/bin/wc -l `

if [[ ${num_success} -eq ${num_tests} ]]; then
    subject_line="Torus benchmark successful. "
else
    subject_line="Torus benchmark failed. "
fi

# Test for success of molebench
num_success=`/usr/bin/grep "Test passed" ${LOG_FILE} | /usr/bin/wc -l `

if [[ ${num_success} -eq 1 ]]; then
    subject_line2="Molecular benchmark successful. "
else
    subject_line2="Molecular benchmark failed. "
fi


for user in ${mail_to}; do
    /sw/bin/mutt -s "${subject_line} ${subject_line2}" ${attach_list} ${user}@${MAIL_DOMAIN} < ${LOG_FILE} 
done

exit

