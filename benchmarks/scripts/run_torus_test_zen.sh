#!/usr/bin/ksh

write_qsub_dd()
{
cat <<EOF > torus.pbs
#!/usr/bin/ksh
#PBS -l nodes=1:ppn=${1}
#PBS -l walltime=8:00:00
#PBS -q all
#PBS -o stdout_torus
#PBS -e stderr_torus
#PBS -N Torus_${THIS_CONFIG}

export NUMBEROFNODES=1
export NUMPROCS=${1}

export TORUS_DATA=${TORUS_DATA}
cd \${PBS_O_WORKDIR}

mpdboot -n \$NUMBEROFNODES -r ssh -f \$PBS_NODEFILE
mpiexec  -genv I_MPI_DEVICE rdssm:OpenIB-cma -np \$NUMPROCS torus.zen > log
mpdallexit

EOF
}

#######################################################################################

write_qsub_mpi()
{
cat << EOF > torus.pbs
#!/usr/bin/ksh
#PBS -l nodes=2:ppn=12
#PBS -l walltime=12:00:00
#PBS -q all
#PBS -o stdout_torus
#PBS -e stderr_torus
#PBS -N Torus_${THIS_CONFIG}

export NUMBEROFNODES=2
export NUMPROCS=24

export TORUS_DATA=${TORUS_DATA}
cd \${PBS_O_WORKDIR}

mpdboot -n \$NUMBEROFNODES -r ssh -f \$PBS_NODEFILE
mpiexec  -genv I_MPI_DEVICE rdssm:OpenIB-cma -np \$NUMPROCS torus.zen > log
mpdallexit

EOF
}

#######################################################################################

write_qsub_hybrid()
{
cat << EOF > torus.pbs
#!/usr/bin/ksh
#PBS -l nodes=2:ppn=12
#PBS -l walltime=12:00:00
#PBS -q all
#PBS -o stdout_torus
#PBS -e stderr_torus
#PBS -N Torus_${THIS_CONFIG}

export NUMBEROFNODES=2
export NUMPROCS=2
export OMP_NUM_THREADS=12

export TORUS_DATA=${TORUS_DATA}
cd \${PBS_O_WORKDIR}

mpdboot -n \$NUMBEROFNODES -r ssh -f \$PBS_NODEFILE
mpiexec -perhost 1 -genv I_MPI_DEVICE rdssm:OpenIB-cma  -genv I_MPI_PIN_DOMAIN omp -np \$NUMPROCS torus.zen > log
mpdallexit

EOF
}

#######################################################################################

write_qsub_openmp()
{
cat << EOF > torus.pbs
#!/usr/bin/ksh
#PBS -l nodes=1:ppn=12
#PBS -q all
#PBS -o stdout_torus
#PBS -e stderr_torus
#PBS -N Torus_${THIS_CONFIG}

export OMP_NUM_THREADS=12
export TORUS_DATA=${TORUS_DATA}
cd \${PBS_O_WORKDIR}
./torus.zensingle > log 

EOF
}

#######################################################################################

test_list="mpi mpi_db openmp openmp_db hybrid hybrid_db"
test_dir=/data/${USER}/torus_tests
export TORUS_DATA=${base_dir}/${test_dir}/torus/data

echo
echo "Running Torus test suite on Zen"
echo "-------------------------------"
echo

if [[ -d ${test_dir}/torus ]]; then 
    echo "Found torus directory in ${test_dir}"
else
    echo "Did not find  torus directory in ${test_dir}. Aborting."
    exit 1
fi

# Set up directories for running the tests
echo
cd $test_dir
for test in ${test_list}; do
    echo "Setting up directory for ${test}"
    if [[ -e ${test} ]]; then
	pwd;echo "${test} already exists. Aborting ..."
	exit 1
    fi
    mkdir ${test}
    cp -r torus/benchmarks ${test}
    mkdir ${test}/build
    cd ${test}/build
    ln -s ../../torus/* .
    cd ../..
done

# Build executables for the various different configurations we will be testing
echo 
for config in ${test_list}; do

    case ${config} in
	openmp) export SYSTEM=zensingle
	    dodebug=no
	    doopenmp=yes;;
        openmp_db) export SYSTEM=zensingle
            dodebug=yes
            doopenmp=yes;;
        mpi) export SYSTEM=zen
            dodebug=no
            doopenmp=no;;
        mpi_db) export SYSTEM=zen
            dodebug=yes
            doopenmp=no;;
        hybrid) export SYSTEM=zen
            dodebug=no
            doopenmp=yes;;
        hybrid_db) export SYSTEM=zen
            dodebug=yes
            doopenmp=yes;;
	*) echo "Unrecognised configuration. Aborting ..."
	    exit 1;;
    esac

    echo "Building Torus for ${config} test"
    cd ${test_dir}/${config}/build 
    make depends > compile_log
    make debug=${dodebug} openmp=${doopenmp} >> compile_log

    if [[ -x torus.${SYSTEM} ]]; then
	echo "Found executable torus.${SYSTEM}"
    else
	echo "Did not find  executable torus.${SYSTEM}.  Aborting ..."
	exit 1
    fi
done


# Submit jobs
echo 
echo "Submitting jobs to queue"

# gfortran is required for setting up sphbench
export  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/sw/gcc/lib64:/sw/gcc/lib
export  PATH=${PATH}:/sw/gcc/bin

# Non-domain decomposed benchamrks
for config in mpi mpi_db openmp openmp_db hybrid hybrid_db; do 
    echo "Submitting ${config} runs"
    export THIS_CONFIG=${config}

    case ${config} in
        openmp) export SYSTEM=zensingle;;
        openmp_db) export SYSTEM=zensingle;;
        mpi) export SYSTEM=zen;;
        mpi_db) export SYSTEM=zen;;
        hybrid) export SYSTEM=zen;;
        hybrid_db) export SYSTEM=zen;;
        *) echo "Unrecognised configuration. Aborting ..."
            exit 1;;
    esac

# Leave the restart test as it relies on the disc benchmark finishing first. This could be set up with a qsub dependency.

# Set up sphbench
    echo "Preparing sphbench"
    cd ${test_dir}/${config}/benchmarks/sphbench
    ./setup_write_sph_file

    for bench in disc disc_cylindrical HII_region molebench nbody angularImageTest sphbench; do
	cd ${test_dir}/${config}/benchmarks/${bench}
	ln -s ../../build/torus.${SYSTEM}
	case ${config} in
            openmp) write_qsub_openmp;;
            openmp_db) write_qsub_openmp;;
            mpi) write_qsub_mpi;;
            mpi_db) write_qsub_mpi;;
            hybrid) write_qsub_hybrid;;
            hybrid_db) write_qsub_hybrid;;
            *) echo "Unrecognised configuration. Aborting ..."
		exit 1;;
	esac
	qsub torus.pbs
    done

# Domain decomposed benchmarks (MPI only)
    if [[ $config == "mpi" || $config == "mpi_db" ]]; then 
	echo "Submitting domain decomposed benchmarks"
# 1D
	for bench in HII_regionMPI hydro; do
	    cd ${test_dir}/${config}/benchmarks/${bench}
	    ln -s ../../build/torus.zen
	    write_qsub_dd 3
	    qsub torus.pbs
	done

# 2D
        for bench in gravtest_2d; do
            cd ${test_dir}/${config}/benchmarks/${bench}
            ln -s ../../build/torus.zen
            write_qsub_dd 5
            qsub torus.pbs
        done

# 3D
        for bench in gravtest cylinder_image_test; do
            cd ${test_dir}/${config}/benchmarks/${bench}
            ln -s ../../build/torus.zen
            write_qsub_dd 9
            qsub torus.pbs
        done
    fi
done

echo "All done"
