#!/usr/bin/ksh

print_help(){
echo "Torus parallel scaling test script"
echo "----------------------------------"
echo 
echo "Call this script with openmp, mpi or both to select which type of parallelism to test."
echo
echo "There needs to be a bin directory containing torus.openmp/torus.mpi executables "
echo "in the directory where the script is run from"
echo
echo "Make a template directory called run which contains your parameters file and any"
echo "input files. These will be linked into the actual run directories."
echo 
echo "Results are written to times.dat in the run directory"
echo 
}

# Write a Torque job file for OpenMP jobs
write_openmp_pbs_file(){
cat << EOF > openmp.pbs
#!/usr/bin/ksh
#PBS -l nodes=1
#PBS -l walltime=24:00:00
#PBS -q all
#PBS -o stdout_openmp
#PBS -e stderr_openmp
#PBS -N TorusOpenMP_${TORUS_NUM_THREADS}

export OMP_NUM_THREADS=${TORUS_NUM_THREADS}

cd \${PBS_O_WORKDIR}
time ./torus.openmp

time=\`grep "Torus Main" tune.dat\`
echo "OpenMP_${TORUS_NUM_THREADS}: \${time}" >> ../times.dat

EOF
}

# Write a Torque job file for MPI jobs
write_mpi_pbs_file(){
cat << EOF > mpi.pbs
#!/usr/bin/ksh
#PBS -l nodes=1
#PBS -l walltime=24:00:00
#PBS -q all
#PBS -o stdout_mpi
#PBS -e stderr_mpi
#PBS -N TorusMPI_${TORUS_NUM_PROCS}

export NUMBEROFNODES=1
export NUMPROCS=${TORUS_NUM_PROCS}

cd \${PBS_O_WORKDIR}
mpdboot -n \$NUMBEROFNODES -r ssh -f \$PBS_NODEFILE
time mpiexec -genv I_MPI_DEVICE rdssm:OpenIB-cma -np \$NUMPROCS ./torus.mpi
mpdallexit

time=\`grep "Torus Main" tune.dat\`
echo "MPI_${TORUS_NUM_PROCS}: \${time}" >> ../times.dat

EOF
}

#########################################################################

# Parse command line arguments
export DO_OPENMP=no
export DO_MPI=no
while [ $# -gt 0 ]
do
    case "$1" in 
        openmp) export DO_OPENMP=yes;; 
        mpi) export DO_MPI=yes;;
    esac
shift
done

#
# Sanity checks
#

# Make sure we have something to do
if [[ ${DO_OPENMP} == "no" && ${DO_MPI} == "no" ]]; then
    print_help
    exit 1
fi

# Check there is a bin directory
if [[ -d bin ]];then 
    echo "Found a bin directory"
else
    echo "Did not find a bin directory"
    echo "Aborting"
    exit 1
fi

# Check we have an OpenMP executable if required
if [[ ${DO_OPENMP} == "yes" ]]; then
    if [[ ! -x bin/torus.openmp ]]; then 
	echo "Did not find executable bin/torus.openmp"
	echo "Aborting"
	exit 1
    else
	echo "Found bin/torus.openmp"
    fi
fi

# Check we have an MPI executable if required
if [[ ${DO_MPI} == "yes" ]]; then
    if [[ ! -x bin/torus.mpi ]]; then 
	echo "Did not find executable bin/torus.mpi"
	echo "Aborting"
	exit 1
    else
	echo "Found bin/torus.mpi"
    fi
fi

# Make sure there is a run directory
if [[ -d run ]]; then 
    echo "Found a template run directory"
else
    echo "Did not find a template run directory"
    exit 1
fi

# Make a new directory for running the tests
# An integer is appended to the directory name 
# to avoid overwriting previous results
i=1
while [[ -e torus_scaling${i} ]]; do
    i=$(($i+1))
done 
mkdir torus_scaling${i}
cd torus_scaling${i}

#
# Now run the tests
#

# Run OpenMP test
if [[ ${DO_OPENMP} == "yes" ]]; then
    for threads in 12 24; do 
	export TORUS_NUM_THREADS=$threads
	mkdir openmp${threads}
	cd openmp${threads}
	ln -s ../../run/* . 
	ln -s ../../bin/torus.openmp
	write_openmp_pbs_file
	qsub openmp.pbs
	cd ..
    done
fi

# Run MPI tests
if [[ ${DO_MPI} == "yes" ]]; then
    for procs in 12 24; do 
	export TORUS_NUM_PROCS=$procs
	mkdir mpi${procs}
	cd mpi${procs}
	ln -s ../../run/* . 
	ln -s ../../bin/torus.mpi
	write_mpi_pbs_file
	qsub mpi.pbs
	cd ..
    done
fi

echo "Jobs submitted. Timings will be written to times.dat"
echo "Exiting"

#####################################
