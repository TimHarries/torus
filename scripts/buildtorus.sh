#!/bin/bash
#
# This script builds Torus for MPI/OpenMP/hybrid configurations
#
# Original: D. Acreman, February 2012
# Updated: October 2012
# Updated: Added Isca, November 2016
# Updated: check for and build cfitio library, August 2018

# Error codes: 1 - problem with toolchain
#              2 - build failure

print_help(){
    echo
    echo "This script builds Torus for MPI/OpenMP/hybrid configurations."
    echo "The Torus source code needs to be available."
    echo "Existing builds will not be deleted but will be updated."
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

base_dir=${PWD}

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

# If no build options have been specified then build with OpenMP only
if [[ $mpi == no && $openmp == no && $hybrid == no && $single == no ]]; then
    echo "Building with OpenMP (this is the default option)"
    openmp=yes
fi

#####################
# Pre-build checks. #
#####################

# Is the Torus source directory present?
if [[ -d src ]]; then
    echo "Found src directory"
else
    echo "Did not find src directory"
    echo -n "Working directory is "
    pwd
    exit 1
fi

# Do we have a make file?
if [[ -e src/Makefile ]]; then
    echo "Found a Makefile in the src directory"
else
    echo "Did not find src/Makefile"
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
# Create version header
#
cd src
./creategitversion
cd ..

####################################
# Choose suitable values of SYSTEM #
####################################

thisHost=`hostname -f`

if [[ $thisHost == login*.cluster.local ]]; then
    echo "This looks like Isca"
    export SYSTEM=isca
    torusFortranCompiler=ifort
    
else

# Look for compiler for OpenMP/serial builds
    torusFortranCompiler=none

    echo "Looking for gfortran"
    which gfortran > /dev/null 2>&1
    if [[ $? -eq 0 ]]; then
	gfortranVersion=`gfortran -v 2>&1 | grep 'gcc version' | awk '{print $3}'`
	echo "Found gfortran version ${gfortranVersion}"
	torusFortranCompiler=gfortran
	export SYSTEM=gfortran
    else
	echo "gfortran not found"
    fi

    echo "Looking for ifort"
    which ifort > /dev/null 2>&1
    if [[ $? -eq 0 ]]; then
	ifortVersion=`ifort --version 2>&1 | head -1 | awk '{print $3}'`
	echo "Found ifort version ${ifortVersion}"
	torusFortranCompiler=ifort
	export SYSTEM=ifort
	ifortMajorVersion=`echo ${ifortVersion} | awk 'BEGIN {FS="."}; {print $1}'`
	if [[ ${ifortMajorVersion} -lt 15 ]]; then
	    echo "Intel Fortran is older than v15. Using -openmp for OpenMP builds"
	    make_args="${make_args} OPENMPFLAG=-openmp"
	fi
    else
	echo "ifort not found"
    fi

    if [[ ${torusFortranCompiler} == none ]]; then
	echo "No fortran compiler found."
	echo "I'm expecting to find ifort or gfortran"
	exit 1
    else
	echo "Building using ${torusFortranCompiler}"
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

########################################
# Check that we can link with cfitsio  #
########################################

if [[ -e ${base_dir}/lib/cfitsio/.libs/libcfitsio.a ]]; then
    echo "Found ${base_dir}/lib/cfitsio/libcfitsio.a"
    export LIBRARY_PATH=${LIBRARY_PATH}:/opt/homebrew/lib:${base_dir}/lib/cfitsio/.libs
    make_args="${make_args} curlflag=yes"
else
    echo "Checking that we can link with a cfitsio library"
    fitsTestDir=build/fitsTest
    if [[ -d $fitsTestDir ]]; then
	echo "Removing existing $fitsTestDir"
	rm -r $fitsTestDir
    fi
    mkdir -p $fitsTestDir
    cd $fitsTestDir
    ln -s ../../benchmarks/disc/check_disc_image.f90

    # Set up link flags for specific SYSTEMS
    link_flags=" -L${base_dir}/lib/cfitsio -lcfitsio"
    if [[ $SYSTEM == gfortran || $SYSTEM == ifort ]]; then
	if [[ -d ${HOME}/cfitsio/lib ]]; then
	    link_flags="${link_flags} -L${HOME}/cfitsio/lib"
	fi
    fi
    if [[ $SYSTEM == gfortran ]]; then
	link_flags="${link_flags} -L/usr/local/lib"
    fi

    echo "Testing with link flags: ${link_flags}"
    ${torusFortranCompiler} -o check_disc_image check_disc_image.f90 ${link_flags} > link_test.out 2>&1

    if [[ -x check_disc_image ]]; then
	echo "Linking with cfitsio works"
    else
	${torusFortranCompiler} -o check_disc_image check_disc_image.f90 ${link_flags} -lcurl > /dev/null 2>&1

	if [[ -x check_disc_image ]]; then
	    echo "Linking with cfitsio works (requires -lcurl)"
	    make_args="${make_args} curlflag=yes"
	else
	    echo "Did not find a working cfitsio. Building a cfitsio library ..."
	    cd ../..
	    if [[ ! -d lib ]]; then
		mkdir lib
	    fi
	    cd lib
	    ../scripts/buildcfitsio.sh
	    export LIBRARY_PATH=${LIBRARY_PATH}:${base_dir}/lib/cfitsio
	    make_args="${make_args} curlflag=yes"
	fi
    fi
    cd ${base_dir}
    rm -r build/fitsTest
fi

###############
# Do builds   #
###############

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
	ln -s ../../src/* .
    fi
    make depends
    make getgitver=no openmp=yes mpi=no $make_args
    if [[ -x torus.${SYSTEM} ]]; then
	echo "OpenMP executable built successfully"
	cp torus.${SYSTEM} ../../bin/torus.openmp
    else
	echo "Build failed. Aborting ..."
	exit 2
    fi
    cd ../..
fi

# Serial build
if [[ $single == yes ]]; then
    echo "Building serial executable"
    builddir=build/single
    if [[ -d $builddir ]]; then
	echo "Found existing $builddir"
	cd $builddir
    else
	mkdir -p $builddir
	cd $builddir
	ln -s ../../src/* .
    fi
    make depends
    make getgitver=no openmp=no mpi=no $make_args
    if [[ -x torus.${SYSTEM} ]]; then
	echo "Serial executable built successfully"
	cp torus.${SYSTEM} ../../bin/torus.single
    else
	echo "Build failed. Aborting ..."
	exit 2
    fi
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
	ln -s ../../src/* .
    fi
    make depends
    make getgitver=no mpi=yes openmp=no $make_args
    if [[ -x torus.${SYSTEM} ]]; then
	echo "MPI executable built successfully"
	cp torus.${SYSTEM} ../../bin/torus.mpi
    else
	echo "Build failed. Aborting ..."
	exit 2
    fi
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
	ln -s ../../src/* .
    fi
    make depends
    make getgitver=no mpi=yes openmp=yes $make_args
    if [[ -x torus.${SYSTEM} ]]; then
	echo "Hybrid executable built successfully"
	cp torus.${SYSTEM} ../../bin/torus.hybrid
    else
	echo "Build failed. Aborting ..."
	exit 2
    fi
    cd ../..
fi

exit 0

# End of file ################################################################
