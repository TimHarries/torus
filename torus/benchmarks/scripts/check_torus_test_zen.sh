#!/usr/bin/ksh

check_disc()
{
cd ${test_dir}/${config}/benchmarks/disc

${TORUS_FC} -o comparespec comparespec.f90
rm -f check_log.txt

echo "Comparing the 12.5 degree model..." >> check_log.txt
if [[ -e test_inc013.dat ]]; then
    cp test_inc013.dat speca.dat
    cp sed100_125.dat specb.dat
    comparespec >> check_log.txt 2>&1
fi

echo "Comparing the 77.5 degree model..." >> check_log.txt
if [[ -e test_inc077.dat ]]; then
    cp test_inc077.dat speca.dat
    cp sed100_775.dat specb.dat
    comparespec >> check_log.txt 2>&1
fi

${TORUS_FC} -o check_disc_image check_disc_image.f90 -lcfitsio -L${FITSLIB}
echo "Checking image" >> check_log.txt
./check_disc_image >> check_log.txt 2>&1


num_success=`grep -c "TORUS: Test successful" check_log.txt`
if [[ ${num_success} -eq 3 ]]; then
    echo "Disc benchmark successful" 
else
    echo "!! Disc benchmark FAILED !!"
    suite_status="FAILED"
fi
}

###################################################################################

check_disc_cylindrical()
{
cd ${test_dir}/${config}/benchmarks/disc_cylindrical
cp ../disc/sed* .
cp ../disc/comparespec.f90 .

${TORUS_FC} -o comparespec comparespec.f90
rm -f check_log.txt

echo "Comparing the 12.5 degree model..." >> check_log.txt
if [[ -e test_inc013.dat ]]; then
    cp test_inc013.dat speca.dat
    cp sed100_125.dat specb.dat
    comparespec >> check_log.txt 2>&1
fi

echo "Comparing the 77.5 degree model..." >> check_log.txt
if [[ -e test_inc077.dat ]]; then
    cp test_inc077.dat speca.dat
    cp sed100_775.dat specb.dat
    comparespec >> check_log.txt 2>&1
fi

num_success=`grep -c "TORUS: Test successful" check_log.txt`
if [[ ${num_success} -eq 2 ]]; then
    echo "3D Disc benchmark successful" 
else
    echo "!! 3D Disc benchmark FAILED !!"
    suite_status="FAILED"
fi
}

###################################################################################

check_angIm()
{
cd ${test_dir}/${config}/benchmarks/angularImageTest
rm -f check_log.txt
${TORUS_FC} -o check check.f90 -lcfitsio -L${FITSLIB} 
./check >> check_log.txt 2>&1

num_success=`grep -c "TORUS: Test successful" check_log.txt`
if [[ ${num_success} -eq 1 ]]; then
    echo "Angular image test successful" 
else
    echo "!! Angular image test FAILED !!"
    suite_status="FAILED"
fi
}

###################################################################################

check_image()
{
cd ${test_dir}/${config}/benchmarks/cylinder_image_test
rm -f check_log.txt

echo "Generating analytical solution" >> check_log.txt
${TORUS_FC} -o cylinder_test cylinder_test.f90
./cylinder_test >> check_log.txt 2>&1

echo "Checking Torus result" >> check_log.txt
${TORUS_FC} -o comparison comparison.f90
./comparison >> check_log.txt 2>&1

num_success=`grep -c "Test Successful" check_log.txt`
if [[ ${num_success} -eq 1 ]]; then
    echo "Cylinder image test successful" 
else
    echo "!! Cylinder image test FAILED !!"
    suite_status="FAILED"
fi

}

###############################################################################

check_grav()
{

cd ${test_dir}/${config}/benchmarks/${THIS_BENCH}
rm -f check_log.txt

echo "Compiling check.f90" >> check_log.txt
${TORUS_FC} -o check check.f90
./check >> check_log.txt 2>&1

num_success=`grep -c "Torus gravity solver test successful" check_log.txt`
if [[ ${num_success} -eq 1 ]]; then
    echo "${THIS_BENCH} successful" 
else
    echo "!! ${THIS_BENCH} FAILED !!"
    suite_status="FAILED"
fi

}

###############################################################################

check_nbody()
{

cd ${test_dir}/${config}/benchmarks/${THIS_BENCH}
rm -f check_log.txt

echo "Compiling check.f90" >> check_log.txt
${TORUS_FC} -o check check.f90
./check >> check_log.txt 2>&1

num_success=`grep -c "Torus nbody test successful" check_log.txt`
if [[ ${num_success} -eq 1 ]]; then
    echo "nbody test successful" 
else
    echo "!! nbody FAILED !!"
    suite_status="FAILED"
fi

}

###############################################################################

check_hII()
{

cd ${test_dir}/${config}/benchmarks/${THIS_BENCH}
rm -f check_log.txt

echo Compiling comparelex code >> check_log.txt
${TORUS_FC} -o comparelex comparelex.f90
./comparelex >> check_log.txt 2>&1

num_success=`grep -c "TORUS: Test successful" check_log.txt`
if [[ ${num_success} -eq 1 ]]; then
    echo "${THIS_BENCH} benchmark successful."
else
    echo "!! ${THIS_BENCH} benchmark FAILED !!"
    suite_status="FAILED"
fi

}

###############################################################################

check_hydro()
{
cd ${test_dir}/${config}/benchmarks/hydro
rm -f check_log.txt

echo Compiling compareSod code >> check_log.txt
${TORUS_FC} -o comparesod compareSod.f90
./comparesod >> check_log.txt 2>&1

num_success=`grep -c "TORUS: Test successful" check_log.txt`
if [[ ${num_success} -eq 1 ]]; then
    echo "Hydro benchmark successful."
else
    echo "!! Hydro benchmark FAILED !!"
    suite_status="FAILED"
fi
}

###############################################################################

check_molebench()
{
cd ${test_dir}/${config}/benchmarks/molebench
rm -f check_log.txt

echo Compiling compare_molbench code >> check_log.txt
${TORUS_FC} -o compare_molbench compare_molbench.f90
./compare_molbench >> check_log.txt 2>&1

echo Compiling check_cube code >> check_log.txt
${TORUS_FC} -o check_cube check_cube.f90 -lcfitsio -L${FITSLIB}
./check_cube >> check_log.txt 2>&1


num_success=`grep -c "TORUS: Test successful" check_log.txt`
if [[ ${num_success} -eq 2 ]]; then
    echo "Molecular benchmark successful."
else
    echo "!! Molecular benchmark FAILED !!"
    suite_status="FAILED"
fi
}

###############################################################################

check_sphbench()
{
cd ${test_dir}/${config}/benchmarks/sphbench
rm -f check_log.txt
./checkSphToGrid.pl log > check_log.txt 2>&1

num_success=`grep -c "TORUS: Test successful" check_log.txt`
if [[ ${num_success} -eq 1 ]]; then
    echo "SPH to grid test successful."
else
    echo "!! SPH to grid test FAILED !!"
    suite_status="FAILED"
fi
 
}

###############################################################################

# Check the output from Torus test suite on Zen 
# D. Acreman, July 2013

test_dir=/data/acreman/torus_tests
export FITSLIB=/home/tjharrie/cfitsio
export TORUS_FC=ifort

# Any failures should set the suite status to failed
suite_status=PASSED

for config in  mpi  mpi_db  openmp  openmp_db hybrid hybrid_db; do

    echo "${config}:"

# angularImageTest
    check_angIm

# disc (2D disc)
    check_disc

# disc_cylindrical
    check_disc_cylindrical

# HII_region
    export THIS_BENCH=HII_region
    check_hII

# molebench
    check_molebench

# nbody
    export THIS_BENCH=nbody
    check_nbody

# sphbench
    check_sphbench

# Domain decomposed only:

    if [[ $config == "mpi"  || $config == "mpi_db" ]]; then
  
  # cylinder_image_test 
	check_image


# gravtest
	export THIS_BENCH=gravtest
	check_grav

# gravtest_2d
	export THIS_BENCH=gravtest_2d
	check_grav

# HII_regionMPI
	export THIS_BENCH=HII_regionMPI
	check_hII

# hydro
	check_hydro

    fi

    echo

done

echo "Test suite ${suite_status}"
