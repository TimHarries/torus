#!/bin/bash

if [[ $# -eq 0 ]]; then
    bench_list="disc"
else
    bench_list=$@
fi

if [[ ! -d ${TORUS_DATA} ]]; then
    echo "TORUS_DATA environment variable does not point to a directory"
    echo "TORUS_DATA=${TORUS_DATA}"
    exit 1
fi

if [[ ! -e ${TORUS_DATA}/bibtex.dat ]]; then
    echo "TORUS_DATA environment variable points to a directory but it does not look like a Torus data directory"
    echo "TORUS_DATA=${TORUS_DATA}"
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
       
done

echo "All done"
exit 0
