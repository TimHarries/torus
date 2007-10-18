#!/bin/ksh

# Name:    build_sphtorus
# Purpose: To compile a coupled version of sphNG and torus
# Author:  D. Acreman (October 2007)

print_help()
{
echo "This is build_sphtorus.sh"
echo ""
echo "This script builds a coupled version of sphNG and torus as a single executable."
echo ""
echo "Command line switches are:"
echo "-d   Compile with debug flags switched on"
echo "-h   Display help"
echo "-o   Overwrite any existing build"
echo "-mpi Compile for multi-processor running using MPI"
}
# Set up the build information here -----------------------------------------------------------------

sph_cvs=${HOME}/sphNG            # This is the directory where the sphNG source code can be found
torus_cvs=${HOME}/torus          # This is the directory where the torus source code can be found
sphtorus_dir=${HOME}/sphtorus    # Build will be carried out in this directory

#----------------------------------------------------------------------------------------------------
echo "INFO: This is build_sphtorus"
# 0. Handle command line arguments
debug_flag=""
overwrite=0
while [ $# -gt 0 ]
do
    case "$1" in 
	-d) echo "INFO: Compiling with debug flags"
	    debug_flag="debug=yes"  # for torus
	    export DEBUG=yes;;      # for sphNG
	-o) echo "INFO: Will overwrite any exiting build"
            overwrite=1;;
	-h) print_help
	    exit;;
	-mpi) echo "INFO: compiling sphtorus using MPI"
	      export SYSTEM="ompi";;
    esac
shift
done

# 1. Check source code directories can be found and set up directories for build

# 1.1 Check that the directory containing the sphNG code is present
if [[ -e ${sph_cvs} ]]; then
    echo "INFO: ${sph_cvs} exists. Will look for the sphNG code in this directory"
else
    echo "ERROR: ${sph_cvs} does not exist. This script will abort as it can't find the sphNG code"
    exit 1
fi

# 1.2 Check that the directory containing the torus code is present
if [[ -e ${torus_cvs} ]]; then
    echo "INFO: ${torus_cvs} exists. Will look for the torus code in this directory"
else
    echo "ERROR: ${torus_cvs} does not exist. This script will abort as it can't find the torus code"
    exit 1
fi

# 1.3 Set up the build directories
if [[ -e ${sphtorus_dir} ]]; then
    if [ ${overwrite} -eq 1 ]; then
	echo "INFO: Removing old build"
	rm -r ${sphtorus_dir}
    else
	echo "ERROR: ${sphtorus_dir} already exists. Will not overwrite existing build."
	echo "ERROR: Call this script with the -o option if you wish to overwrite the existing build."
	exit 1
    fi
else
    echo "INFO:  ${sphtorus_dir} does not exist. Creating directories required for build"
fi

mkdir -p ${sphtorus_dir}/build/torus
mkdir -p ${sphtorus_dir}/build/sphNG
mkdir -p ${sphtorus_dir}/bin
export TORUS_LIB=${sphtorus_dir}/lib # Used by Makefile
mkdir -p ${TORUS_LIB}

# 2. Build torus as a library

# 2.1 Make symbolic links to torus source code in build directories
#     Using *.f* will include the .raw files
cd ${sphtorus_dir}/build/torus
ln -s ${torus_cvs}/*.f* .
ln -s ${torus_cvs}/*.F* .
ln -s ${torus_cvs}/Makefile .
ln -s ${torus_cvs}/makedepend .

# 2.2 Build torus
echo "INFO: Building torus, SYSTEM=${SYSTEM}"
make depends
make ${debug_flag} sph=yes lib
echo "INFO: Moving libtorus.a to ${sphtorus_dir}/lib"
mv libtorus.a ${sphtorus_dir}/lib

# 3.0 Build sphNG coupled to torus

# 3.1 Make symbolic links to source code in build directory
cd  ${sphtorus_dir}/build/sphNG
ln -s ${sph_cvs}/*.f .
ln -s ${sph_cvs}/*.F .
ln -s ${sph_cvs}/Makefile .
ln -s ${sph_cvs}/idim .
ln -s ${sph_cvs}/igrape .
ln -s ${sph_cvs}/inspho .
ln -s ${sph_cvs}/*.h .
ln -s ${sph_cvs}/COMMONS .

echo "INFO: Building sphtorus, SYSTEM=${SYSTEM}"
make updatedepends
make ${debug_flag} sphtorus
if [[ -e sphtorus ]]; then
    echo "INFO: moving executable sphtorus to ${sphtorus_dir}/bin"
    mv sphtorus ${sphtorus_dir}/bin
else
    echo "ERROR: executable was not created"
    exit 1
fi

echo "INFO: Exiting normally"
exit

