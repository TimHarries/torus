#/bin/csh
thisDir=`pwd`
cd ../torus
echo 'character(len=20), parameter :: svnversion = "-'`svnversion`'"' > ${thisDir}/svn_version.h
cd $thisDir

