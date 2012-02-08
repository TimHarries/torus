#!/bin/ksh

# This script builds Torus for MPI/OpenMp/hybrid configurations
# D. Acreman, February 2012

# Known problems
#
# 1. svn version number is not correct
# 2. There's no way to build with -j flag

print_help(){
    echo 
    echo "This script builds Torus for MPI/OpenMp/hybrid configurations."
    echo "The Torus source code needs to be available."
    echo "Exiting builds will not be deleted but will be updated."
    echo "Torus executables will be put in the bin directory"
    echo "Run script with arguments openmp/mpi/hybrid/all."
    echo 
}

####################################
# Start here                       #
####################################

openmp=no
mpi=no
hybrid=no

#
# Parse command line flags
#
while [ $# -gt 0 ]
do
    case "$1" in 
	-h) print_help
	    exit;;
	all) openmp=yes
	    mpi=yes
	    hybrid=yes;;
	openmp) openmp=yes;;
	mpi) mpi=yes;;
	hybrid) hybrid=yes;;
	*) echo "Unrecognised option $1";;
    esac
shift
done

#
# Pre-build checks. 
#

# Is the Torus directory present? 
if [[ -d torus ]]; then
    echo "Found a torus directory"
else
    echo "Did not find a torus directory"
    echo -n "Working directory is "
    pwd
    exit 1
fi

# Do we have a make file? 
if [[ -e torus/Makefile ]]; then
    echo "Found a Makefile in the torus directory"
else
    echo "Did not find torus/Makefile"
fi

# Do we have a bin directory?
if [[ -d bin ]];then 
    echo "Found bin directory"
else
    echo "Making bin directory"
    mkdir bin
fi

#
# Do builds 
#

# OpenMP build 
if [[ $openmp == yes ]]; then
    echo "Building Torus wih OpenMP"
    builddir=build/openmp
    if [[ -d $builddir ]]; then 
	echo "Found existing $builddir"
	cd $builddir
    else
	mkdir -p $builddir
	cd $builddir
	ln -s ../../torus/* . 
    fi
    make depends 
    make SYSTEM=zensingle openmp=yes
    cp torus.zensingle ../../bin/torus.openmp
    cd ../.. 
fi

# MPI build
if [[ $mpi == yes ]]; then
    echo "Building Torus with MPI"
    builddir=build/mpi
    if [[ -d $builddir ]]; then 
	echo "Found existing $builddir"
	cd $builddir
    else
	mkdir -p $builddir
	cd $builddir
	ln -s ../../torus/* . 
    fi
    make depends 
    make SYSTEM=zen
    cp torus.zensingle ../../bin/torus.mpi
    cd ../..
fi

# Hybrid build
if [[ $hybrid == yes ]]; then
    echo "Building hybrid Torus"
    builddir=build/hybrid
    if [[ -d $builddir ]]; then 
	echo "Found existing $builddir"
	cd $builddir
    else
	mkdir -p $builddir
	cd $builddir
	ln -s ../../torus/* . 
    fi
    make depends 
    make SYSTEM=zen openmp=yes
    cp torus.zensingle ../../bin/torus.hybrid
    cd ../..
fi

# Print the help message if nothing was done.
if [[ $mpi == no && $openmp == no && $hybrid == no ]]; then 
    print_help
    exit 1
fi

# End of file
