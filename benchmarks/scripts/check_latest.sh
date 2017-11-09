#!/bin/bash

rundir=/data/torustest/torus_latest
if [[ -d ${rundir} ]]; then
    cd ${rundir}
else
    mkdir ${rundir}
    cd ${rundir}
fi

# A lock file prevents multiple tests running at once
if [[ -e lock ]];then
    echo "Found lock file. Aborting ..."
    exit 1
fi
touch lock

rm -f git_log 

# Make sure there is an up to date local working copy
if [[ -d torus ]]; then
    echo "Found local working copy"
    echo "Updating working copy"
    cd torus; git pull origin master; cd ..
else
    echo "No working copy found. Checking out ... "
    /usr/bin/git clone git@bitbucket.org:tjharries/torus.git
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
cd torus; /usr/bin/git show >> ../git_log; cd ..
thisRev=`head -1 git_log | awk '{print $2}'`
echo "This revision is $thisRev"

# Test stuff here
if [[ ${thisRev} == ${prevRev} ]];then
    echo "No commits since last check"
else
    echo "New commit, Checking..."
    /home/torustest/bin/run_torus_test.sh -b > build_log
    if [[ $? -eq 0 ]]; then
	echo "Build successful"
    else
	echo "Build failed"
	culprit=`/bin/grep Author: git_log | awk '{print $2}'`
	culprit_email=`/bin/grep Author: git_log | awk '{print $3}'`
	echo "Last change was by ${culprit} ${culprit_email}"
	/usr/bin/mail -s "Torus build failed" "${culprit_email}" < /data/torustest/torus_latest/build_log
	cd ${rundir}
    fi
fi

# Prepare for next time 
rm -f prevRev
echo ${thisRev} > prevRev
rm lock

# End of file ####################################################################################################
