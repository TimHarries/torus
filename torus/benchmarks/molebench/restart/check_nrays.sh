#!/bin/bash 

# Check that the number of rays used in the first iteration of the restarted calculation matches the 
# number saved in restart.dat

if [[ $# -ne 1 ]]; then
    echo "Requires one argument (name of log file)"
    exit 1
fi

echo "Checking number of rays against saved value"

nrays_saved=`head -1 restart.dat  | awk '{print $2}'`
nrays_used=`grep "Maximum fractional change this iteration" $1 | head -1 | awk '{print $14}'`

if [[ ${nrays_saved} -eq ${nrays_used} ]]; then
    echo "Number of rays equal to saved value"
    echo "TORUS: Test successful"
else
    echo "Number of rays NOT equal to saved value: ${nrays_used} vs. ${nrays_saved}"
    echo "TORUS: test failed"
fi
