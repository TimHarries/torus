#!/usr/bin/ksh

# Convert axis units in a VTK file
# Written by Dave Acreman, October 2009

if [[ $# != 3 ]]; then
    echo "Incorrect number of arguments"
    echo "This script needs three arguments: infile outfile units"
    exit 1
fi

infile=$1
outfile=$2
new_units=$3

num_lines_header=5

# Set up the numerical conversion factor for the new units
case ${new_units} in
    kpc) scalefac=3.086e11
	 echo "Converting to kiloparsecs";;
    pc) scalefac=3.086e8
	 echo "Converting to parsecs";;
    au) scalefac=1.495979e3
	echo "Converting to AU";;
    *) echo "${new_units} is not a unit I understand"
       exit 1
esac

# Copy header to new file
head -${num_lines_header} ${infile} > ${outfile}

# Work out how many points there are
num_pts=`awk '/POINTS/{print $2}' ${outfile}`
echo "Processing ${num_pts} points"

# Copy and scale points
num_lines=$((${num_lines_header}+${num_pts}))
head -${num_lines} ${infile} | tail -${num_pts} | awk -v fac=${scalefac} '{print $1/fac, $2/fac, $3/fac}' >> ${outfile}

# Copy the remaining lines to the new file
total_lines=`wc -l ${infile} | awk '{print $1}'`
num_lines_left=$((${total_lines}-${num_lines}))
tail -${num_lines_left} ${infile} >> ${outfile}

echo "All done"
exit 0 

