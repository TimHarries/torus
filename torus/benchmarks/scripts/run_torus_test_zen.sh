#!/usr/bin/ksh

write_qsub_hydro()
{
cat <<EOF > torus.pbs
#!/usr/bin/ksh
#PBS -l nodes=1
#PBS -l walltime=1:00:00
#PBS -q all
#PBS -o stdout_torus
#PBS -e stderr_torus
#PBS -N Torus_MPI

export NUMBEROFNODES=1
export NUMPROCS=3

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
#PBS -l nodes=2
#PBS -l walltime=12:00:00
#PBS -q all
#PBS -o stdout_torus
#PBS -e stderr_torus
#PBS -N Torus_MPI

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
#PBS -l nodes=2
#PBS -l walltime=12:00:00
#PBS -q all
#PBS -o stdout_torus
#PBS -e stderr_torus
#PBS -N Torus_hybrid

export NUMBEROFNODES=2
export NUMPROCS=2
export OMP_NUM_THREADS 12

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
#PBS -l nodes=1
#PBS -q mpiexpress
#PBS -o stdout_torus
#PBS -e stderr_torus
#PBS -N Torus_OpenMP

export OMP_NUM_THREADS=12
export TORUS_DATA=${TORUS_DATA}
cd \${PBS_O_WORKDIR}
./torus.zensingle > log 

EOF
}

#######################################################################################

# To do: 1. Redirect stderr from build to record failures (check build results?)
#        2. Add new benchmarks

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

cd $test_dir

for test in ${test_list}; do
    mkdir ${test}
    cp -r torus/benchmarks ${test}
    mkdir ${test}/build
    cd ${test}/build
    ln -s ../../torus/* .
    cd ../..
done

echo "Building Torus for openmp test"
cd ${test_dir}/openmp/build 
export SYSTEM=zensingle
make depends > compile_log
make debug=no openmp=yes >> compile_log

echo "Building Torus for openmp_db test"
cd ${test_dir}/openmp_db/build
export SYSTEM=zensingle
make depends > compile_log
make debug=yes openmp=yes >> compile_log

echo "Building Torus for mpi test"
cd ${test_dir}/mpi/build
export SYSTEM=zen
make depends > compile_log
make debug=no openmp=no >> compile_log

echo "Building Torus for mpi_db test"
cd /${test_dir}/mpi_db/build
export SYSTEM=zen
make depends > compile_log
make debug=yes openmp=no >> compile_log

echo "Building Torus for hybrid test"
cd ${test_dir}/hybrid/build
export SYSTEM=zen
make depends > compile_log
make debug=no openmp=yes >> compile_log

echo "Building Torus for hybrid_db test"
cd ${test_dir}/hybrid_db/build
export SYSTEM=zen
make depends > compile_log
make debug=yes openmp=yes >> compile_log

# Submit jobs
echo "Submitting jobs to queue"

echo "Submitting mpi runs"
for bench in disc disc_cylindrical HII_region molebench; do
    cd ${test_dir}/mpi/benchmarks/${bench}
    ln -s ../../build/torus.zen 
    write_qsub_mpi
    qsub torus.pbs
done
cd ${test_dir}/mpi/benchmarks/hydro
ln -s ../../build/torus.zen
write_qsub_hydro
qsub torus.pbs

echo "Submitting mpi_db runs"
for bench in disc disc_cylindrical HII_region molebench; do
    cd ${test_dir}/mpi_db/benchmarks/${bench}
    ln -s ../../build/torus.zen
    write_qsub_mpi
    qsub torus.pbs
done
cd ${test_dir}/mpi_db/benchmarks/hydro
ln -s ../../build/torus.zen
write_qsub_hydro
qsub torus.pbs

echo "Submitting openmp runs"
for bench in disc disc_cylindrical HII_region molebench; do
    cd ${test_dir}/openmp/benchmarks/${bench}
    ln -s ../../build/torus.zensingle
    write_qsub_openmp
    qsub torus.pbs
done

echo "Submitting openmp_db runs"
for bench in disc disc_cylindrical HII_region molebench; do
    cd ${test_dir}/openmp_db/benchmarks/${bench}
    ln -s ../../build/torus.zensingle
    write_qsub_openmp
    qsub torus.pbs
done

echo "Submitting hybrid runs"
for bench in disc disc_cylindrical HII_region molebench; do
    cd ${test_dir}/hybrid/benchmarks/${bench}
    ln -s ../../build/torus.zen
    write_qsub_hybrid
    qsub torus.pbs
done

echo "Submitting hybrid_db runs"
for bench in disc disc_cylindrical HII_region molebench; do
    cd ${test_dir}/hybrid_db/benchmarks/${bench}
    ln -s ../../build/torus.zen
    write_qsub_hybrid
    qsub torus.pbs
done

echo "All done"
