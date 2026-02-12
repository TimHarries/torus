#!/bin/bash

# Purpose: Download and build the cfitsio library
# Author: D. Acreman, August 2018

cfitsio_location="https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c"
cfitsio_file="cfitsio_latest.tar.gz"

# Download source code using ftp or curl
which ftp > /dev/null 2>&1
if [[ $? -eq 0 ]]; then
    echo "Downloading ${cfitsio_file} from ${cfitsio_location} using ftp"
    ftp ${cfitsio_location}/${cfitsio_file}
else
    which curl > /dev/null 2>&1
    if [[ $? -eq 0 ]]; then
	echo "Downloading ${cfitsio_file} from ${cfitsio_location} using curl"
	curl -o ${cfitsio_file} ${cfitsio_location}/${cfitsio_file}
    else
	echo "Could not find ftp or curl. Aborting because I can't download the source code."
	exit 1
    fi
fi

if [ ! -e ${cfitsio_file} ]; then
    echo "${cfitsio_file} was not downloaded"
    exit 1
fi

echo "Unpacking ${cfitsio_file}"
tar zxf ${cfitsio_file}
rm ${cfitsio_file}
if [[ -d cfitsio ]]; then
    echo "Unpacked to cfitsio"
else
    echo "Linking to cfitsio"
    ln -s cfitsio* cfitsio
    cd cfitsio
fi

if [[ -e ./configure ]]; then
    echo "Configuring ..."
    ./configure --enable-static=yes --enable-shared=no
else
    echo "Error: configure not found"
    exit 1
fi
  
# Check return code and run make if the configure command was successful
if [[ $? -eq 0 ]]; then
    echo "Making ..."
    make
else
    echo "Error in configure step."
    exit1
fi

if [[ $? -eq 0 ]]; then
    echo "Make succeeded."
else
    echo "Error in make step."
    exit 1
fi

# End of file ##########################################################################
