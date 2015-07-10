#!/bin/ksh

rundir=/data/torustest/torus_latest

cd ${rundir}

# A lock file prevents multiple tests running at once
if [[ -e lock ]];then
    echo "Found lock file. Aborting ..."
    exit 1
fi
touch lock

rm -f svn_log 

# Make sure there is a local working copy
if [[ -d torus ]]; then
    echo "Found local working copy"
else
    echo "No working copy found. Checking out ... "
    /usr/bin/svn co https://repository.astro.ex.ac.uk/torus/trunk/torus torus > svn_log
fi

# Check previous revison number if available. 
if [[ -e prevRev ]];then
    prevRev=`cat prevRev`
    echo "Previous revision was ${prevRev}"
else
    echo "Previous revision not found. Setting to 1"
    echo 1 > prevRev
fi

# Get the current revision number
/usr/bin/svn --show-updates status torus >> svn_log
thisRev=`grep "Status against revision" svn_log | awk '{print $4}'`
echo "This revision is $thisRev"

# Test stuff here
if [[ ${thisRev} -eq ${prevRev} ]];then
    echo "No commits since last check"
else
    echo "New commit, Checking..."
    /home/torustest/bin/run_torus_test.sh -b > build_log
    if [[ $? -eq 0 ]]; then
	echo "Build successful"
    else
	echo "Build failed"
	/usr/bin/svn update torus
	culprit=`/usr/bin/svn info torus | grep 'Last Changed Author' | awk '{print $4}'`
	echo "Last change was by ${culprit}"
	cd ${HOME}/torus_build_tests
	attach_list=""
	for file in build_only_*/build/compile_log_*.txt; do
	    attach_list="${attach_list} -a ${file}"
	done
	echo ${culprit}@astro.ex.ac.uk > ${rundir}/ready
#	/sw/bin/mutt -s "Torus build failed" ${attach_list} ${culprit}@astro.ex.ac.uk < ~/torus_latest/build_log
	cd ${rundir}
    fi
fi

# Prepare for next time 
rm -f prevRev
echo ${thisRev} > prevRev
rm lock

