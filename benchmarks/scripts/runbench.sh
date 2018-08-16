#!/bin/bash
# 
# Purpose: Script for running Torus benchmarks
# Author:  D. Acreman
#
# To do: testing is only done for the disc benchmark at present
#        but others can be added if required. 

set -u

if [[ $# -eq 0 ]]; then
    bench_list="disc"
else
    bench_list=$@
fi

if [[ ! -d ${TORUS_DATA} ]]; then
    echo "ERROR: TORUS_DATA environment variable does not point to a directory"
    echo "TORUS_DATA=${TORUS_DATA}"
    exit 1
fi

if [[ ! -e ${TORUS_DATA}/bibtex.dat ]]; then
    echo "ERROR: TORUS_DATA environment variable points to a directory but it does not look like a Torus data directory"
    echo "TORUS_DATA=${TORUS_DATA}"
    exit 1
fi

if [[ ! -x bin/torus.openmp ]]; then
    echo "ERROR: did not find an executable torus.openmp in the bin directory"
    echo "Try running buildtorus before running this script"
    exit 1
fi

if [[ -d run ]]; then
    echo "Found run directory"
else
    echo "Making run directory"
    mkdir run
fi
cd run

echo "Benchmarks to run: ${bench_list}"
for benchmark in $bench_list; do
    echo "Running ${benchmark}"

    if [[ -e ${benchmark} ]]; then
	# Keep old runs
	now=`date "+%Y-%m-%d-%H:%M:%S"`
	mv ${benchmark} ${benchmark}.${now}
    fi

    cp -r ../benchmarks/${benchmark} .
    cd ${benchmark}
    ln -s ../../bin/torus.openmp

    ./torus.openmp 2>&1 | tee run_log_${benchmark}.txt

    if [[ -e ./check_${benchmark}.sh ]]; then
	./check_${benchmark}.sh
    else
	echo "Test not configured for this benchmark yet"
    fi

done

echo "All done"
exit 0
