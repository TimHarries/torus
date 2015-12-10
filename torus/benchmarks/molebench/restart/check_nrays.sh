#!/bin/bash 

# Check that the number of rays used in the first iteration of the restarted calculation matches the 
# number saved in restart.dat

if [[ $# -ne 1 ]]; then
    echo "Requires one argument (name of log file)"
    exit 3
fi

echo "Checking number of rays against saved value"

# Take the number of saved rays from the restart file in the original run directory
# testDump_restart.dat is written at the end of the first non-fixed ray iteration if
# the  geometry is molebench
if [[ -e ../testDump_restart.dat ]]; then
    nrays_saved=`head -1 ../testDump_restart.dat  | awk '{print $2}'`
    nrays_used=`grep "Maximum fractional change this iteration" $1 | head -1 | awk '{print $14}'`
else
    echo "testDump_restart.dat not found"
    echo "TORUS: test failed"
    exit 2
fi

if [[ ${nrays_saved} -eq ${nrays_used} ]]; then
    echo "Number of rays equal to saved value"
    echo "TORUS: Test successful"
    exit 0
else
    echo "Number of rays NOT equal to saved value: ${nrays_used} vs. ${nrays_saved}"
    echo "TORUS: test failed"
    exit 1
fi
