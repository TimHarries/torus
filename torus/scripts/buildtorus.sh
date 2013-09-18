#!/bin/ksh

# This script builds Torus for MPI/OpenMP/hybrid configurations
# Set up to work on Zen but will make an educated guess for other machines
# The handling of MPI wrappers is basic and assumes that ifort is the back end 
# compiler if it is available.
#
# Original: D. Acreman, February 2012
# Updated: October 2012

print_help(){
    echo 
    echo "This script builds Torus for MPI/OpenMP/hybrid configurations."
    echo "The Torus source code needs to be available."
    echo "Exiting builds will not be deleted but will be updated."
    echo "Torus executables will be put in the bin directory"
    echo "Run script with arguments openmp/mpi/hybrid/all/single."
    echo "Other arguments are passed to the make command"
    echo 
}

####################################
# Start here                       #
####################################

echo
echo "Torus build script"
echo "------------------"
echo 

openmp=no
mpi=no
hybrid=no
single=no

############################
# Parse command line flags #
############################
make_args=
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
	single) single=yes;;
	*) make_args="${make_args} $1";;
    esac
shift
done

#####################
# Pre-build checks. #
#####################

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
    exit 1
fi

# Do we have a bin directory?
if [[ -d bin ]];then 
    echo "Found bin directory"
else
    echo "Making bin directory"
    mkdir bin
fi

#
# Create svn version header
#
cd torus
./createsvnversion.csh
cd ..

####################################
# Choose suitable values of SYSTEM #
####################################

thisHost=`hostname`

if [[ $thisHost == service0 || $thisHost == service2 ]]; then
    echo "This looks like Zen"
    system_mpi=zen
    system_hybrid=zen
    system_openmp=zensingle

elif [[ $thisHost == dirac* ]]; then
    echo "This looks like Complexity"
    system_mpi=complexity
    system_hybrid=complexity
    system_openmp=complexity
    
else

# Look for compiler for OpenMP build 
    torusFortranCompiler=none

    echo "Looking for gfortran"
    which gfortran > /dev/null
    if [[ $? -eq 0 ]]; then
	echo "Found gfortran"
	torusFortranCompiler=gfortran
	system_openmp=gfortran
	system_mpi=ompiosx
	system_hybrid=ompiosx
    else
	echo "gfortran not found"
    fi

    echo "Looking for ifort"
    which ifort > /dev/null
    if [[ $? -eq 0 ]]; then
	echo "Found ifort"
	torusFortranCompiler=ifort
	system_openmp=zensingle
	system_mpi=zen
	system_hybrid=zen
    else
	echo "ifort not found"
    fi

    if [[ ${torusFortranCompiler} == none ]]; then
	echo "No fortran compiler found."
	echo "I'm expecting to find ifort or gfortran"
	exit 1
    fi

# Look for mpif90
    if [[ ${mpi} == yes || ${hybrid} == yes ]]; then
	which mpif90 > /dev/null
	if [[ $? -eq 0 ]]; then
	    echo "Found mpif90"
	else
	    echo "mpif90 not found. Aborting ..."
	    exit 1
	fi 
    fi 

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
    make getsvnver=no SYSTEM=${system_openmp} openmp=yes mpi=no $make_args
    cp torus.${system_openmp} ../../bin/torus.openmp
    cd ../.. 
fi

# Serial build 
if [[ $single == yes ]]; then
    echo "Building Torus wih no parallelisation"
    builddir=build/single
    if [[ -d $builddir ]]; then 
	echo "Found existing $builddir"
	cd $builddir
    else
	mkdir -p $builddir
	cd $builddir
	ln -s ../../torus/* . 
    fi
    make depends 
    make getsvnver=no SYSTEM=${system_openmp} openmp=no mpi=no $make_args
    cp torus.${system_openmp} ../../bin/torus.single
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
    make getsvnver=no SYSTEM=${system_mpi} $make_args
    cp torus.${system_mpi} ../../bin/torus.mpi
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
    make getsvnver=no SYSTEM=${system_hybrid} openmp=yes $make_args
    cp torus.${system_hybrid} ../../bin/torus.hybrid
    cd ../..
fi

# Print the help message if nothing was done.
if [[ $mpi == no && $openmp == no && $hybrid == no && $single == no ]]; then 
    print_help
    exit 1
fi

exit 0

# End of file
