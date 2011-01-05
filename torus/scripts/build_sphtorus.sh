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
echo "-p   Compile with profiling options enabled"
echo "-mpi Compile for multi-processor running using MPI"
echo "-openmp Compile with OpenMP"
}
# Set up the build information here -----------------------------------------------------------------
BASE_DIR=${PWD}

sph_cvs=${BASE_DIR}/sphNG            # This is the directory where the sphNG source code can be found
torus_cvs=${BASE_DIR}/torus          # This is the directory where the torus source code can be found
sphtorus_dir=${BASE_DIR}/sphtorus    # Build will be carried out in this directory

#----------------------------------------------------------------------------------------------------
echo "INFO: This is build_sphtorus"
# 0. Handle command line arguments
debug_flag=""
profile_flag=""
openmp_flag=""
openmp_lib_flag=""
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
	-p) profile_flag="profile=yes"
	    export PROFILE=yes;;
	-mpi) echo "INFO: compiling sphtorus using MPI"
	      export SYSTEM="ompi"
	      export mpi="yes";;   # for sphNG 
        -itac) echo "INFO: Using ITAC profiling"
	       export trace=yes
	       itac_flag="trace=yes";;
        -openmp) echo "INFO: Using OpenMP"
	         openmp_lib_flag=-openmp
                 openmp_flag="openmp=yes";;
    esac
shift
done

case ${SYSTEM} in 
    zen) export mpi="yes"
esac

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
# Specify libraries for linking (used by sphNG Makefile)
export TORUS_LIB="${sphtorus_dir}/lib -ltorus ${openmp_lib_flag}"

if [[ -e ${sphtorus_dir} ]]; then
    if [ ${overwrite} -eq 1 ]; then
	echo "INFO: Removing old build"
	rm -r ${sphtorus_dir}
	mkdir -p ${sphtorus_dir}/build/torus
	mkdir -p ${sphtorus_dir}/build/sphNG
	mkdir -p ${sphtorus_dir}/bin
	mkdir -p ${sphtorus_dir}/lib
	make_links=1
    else
	echo "INFO: ${sphtorus_dir} already exists. Performing incremental build."
	make_links=0
    fi
else
    echo "INFO:  ${sphtorus_dir} does not exist. Creating directories required for build"
    mkdir -p ${sphtorus_dir}/build/torus
    mkdir -p ${sphtorus_dir}/build/sphNG
    mkdir -p ${sphtorus_dir}/bin
    mkdir -p ${sphtorus_dir}/lib
    make_links=1
fi

# 2. Build torus as a library

# 2.1 Make symbolic links to torus source code in build directories
cd ${sphtorus_dir}/build/torus
if [[ ${make_links} -eq 1 ]]; then
    ln -s ${torus_cvs}/* .
fi

# 2.2 Build torus
echo "INFO: Building torus, SYSTEM=${SYSTEM}"
make depends
make ${debug_flag} ${profile_flag} ${itac_flag} ${openmp_flag} cfitsio=no lib
if [[ -e libtorus.a ]]; then
  echo "INFO: Moving libtorus.a to ${sphtorus_dir}/lib"
else
  echo "ERROR: libtorus.a not created"
  exit 1
fi
mv libtorus.a ${sphtorus_dir}/lib

# 3.0 Build sphNG coupled to torus

# 3.1 Make symbolic links to source code in build directory
cd  ${sphtorus_dir}/build/sphNG
if [[ ${make_links} -eq 1 ]]; then
    ln -s ${sph_cvs}/* .
    ln -s ../torus/*.mod .
fi

echo "INFO: Building sphtorus, SYSTEM=${SYSTEM}"
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

