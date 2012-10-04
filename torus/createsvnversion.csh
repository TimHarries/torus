#!/bin/csh

# Check whether subversion is available
which svnversion > /dev/null

# Set the version number using svnversion if available or set to unknown otherwise
if ($? == 0) then 
    echo 'character(len=20), parameter :: svnversion = "-'`svnversion`'"' > svn_version.h
else
    echo 'character(len=20), parameter :: svnversion = "-'unknown'"' > svn_version.h
endif

