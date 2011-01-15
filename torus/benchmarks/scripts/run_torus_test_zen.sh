#!/usr/bin/ksh

write_qsub_hydro()
{
cat <<EOF > torus.pbs
#!/usr/bin/ksh
#PBS -l nodes=1:ppn=8
#PBS -q mpiexpress
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
#PBS -l nodes=8:ppn=8
#PBS -q mpiexpress
#PBS -o stdout_torus
#PBS -e stderr_torus
#PBS -N Torus_MPI

export NUMBEROFNODES=8
export NUMPROCS=64

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
#PBS -l nodes=8:ppn=8
#PBS -q mpiexpress
#PBS -o stdout_torus
#PBS -e stderr_torus
#PBS -N Torus_hybrid

export NUMBEROFNODES=8
export NUMPROCS=64

export TORUS_DATA=${TORUS_DATA}
cd \${PBS_O_WORKDIR}

mpdboot -n \$NUMBEROFNODES -r ssh -f \$PBS_NODEFILE
mpiexec -genv I_MPI_DEVICE rdssm:OpenIB-cma  -genv I_MPI_PIN_DOMAIN node -np \$NUMPROCS torus.zen < input
mpdallexit

EOF
}

#######################################################################################

write_qsub_openmp()
{
cat << EOF > torus.pbs
#!/usr/bin/ksh
#PBS -l nodes=1:ppn=8
#PBS -q mpiexpress
#PBS -o stdout_torus
#PBS -e stderr_torus
#PBS -N Torus_OpenMP

export TORUS_DATA=${TORUS_DATA}
cd \${PBS_O_WORKDIR}
./torus.zensingle > log 

EOF
}

#######################################################################################

export CVSROOT=:ext:${USER}@reduce.astro.ex.ac.uk:/home/cvs/th
export CVS_RSH=ssh
test_list="mpi mpi_db openmp openmp_db hybrid hybrid_db"
test_dir=torus_tests
base_dir=/scratch/acreman
export TORUS_DATA=${base_dir}/${test_dir}/torus/data

echo "Running Torus tests for Zen"

cd ${base_dir}

if [[ -d $test_dir ]]; then
    echo "Test directory $test_dir already exists. Aborting ..."
    exit 1
fi

mkdir $test_dir
cd $test_dir

echo "Checking out Torus from CVS"
cvs -q co torus > cvs_log.txt 2>&1

for test in ${test_list}; do
    mkdir ${test}
    cp -r torus/benchmarks ${test}
    mkdir ${test}/build
    cd ${test}/build
    ln -s ../../torus/* .
    cd ../..
done

echo "Building Torus for openmp test"
cd ${base_dir}/${test_dir}/openmp/build 
export SYSTEM=zensingle
make depends > compile_log
#make debug=no openmp=yes >> compile_log

echo "Building Torus for openmp_db test"
cd ${base_dir}/${test_dir}/openmp_db/build
export SYSTEM=zensingle
make depends > compile_log
#make debug=yes openmp=yes >> compile_log

echo "Building Torus for mpi test"
cd ${base_dir}/${test_dir}/mpi/build
export SYSTEM=zen
make depends > compile_log
make debug=no openmp=no >> compile_log

echo "Building Torus for mpi_db test"
cd ${base_dir}/${test_dir}/mpi_db/build
export SYSTEM=zen
make depends > compile_log
#make debug=yes openmp=no >> compile_log

echo "Building Torus for hybrid test"
cd ${base_dir}/${test_dir}/hybrid/build
export SYSTEM=zen
make depends > compile_log
#make debug=no openmp=yes >> compile_log

echo "Building Torus for hybrid_db test"
cd ${base_dir}/${test_dir}/hybrid_db/build
export SYSTEM=zen
make depends > compile_log
#make debug=yes openmp=yes >> compile_log

# Submit jobs
echo "Submitting jobs to queue"

echo "Submitting mpi runs"
for bench in disc disc_cylindrical HII_region molebench; do
    cd ${base_dir}/${test_dir}/mpi/benchmarks/${bench}
    ln -s ../../build/torus.zen 
    write_qsub_mpi
    qsub torus.pbs
done
cd ${base_dir}/${test_dir}/mpi/benchmarks/hydro
ln -s ../../build/torus.zen
write_qsub_hydro
qsub torus.pbs

exit

echo "Submitting mpi_db runs"
for bench in disc disc_cylindrical HII_region molebench; do
    cd ${base_dir}/${test_dir}/mpi_db/benchmarks/${bench}
    ln -s ../../build/torus.zen
    write_qsub_mpi
    qsub torus.pbs
done
cd ${base_dir}/${test_dir}/mpi_db/benchmarks/hydro
ln -s ../../build/torus.zen
write_qsub_hydro
qsub torus.pbs

echo "Submitting openmp runs"
for bench in disc disc_cylindrical HII_region molebench; do
    cd ${base_dir}/${test_dir}/openmp/benchmarks/${bench}
    ln -s ../../build/torus.zensingle
    write_qsub_openmp
    qsub torus.pbs
done

echo "Submitting openmp_db runs"
for bench in disc disc_cylindrical HII_region molebench; do
    cd ${base_dir}/${test_dir}/openmp_db/benchmarks/${bench}
    ln -s ../../build/torus.zensingle
    write_qsub_openmp
    qsub torus.pbs
done

echo "Submitting hybrid runs"
for bench in disc disc_cylindrical HII_region molebench; do
    cd ${base_dir}/${test_dir}/hybrid/benchmarks/disc
    ln -s ../../build/torus.zen
    write_qsub_hybrid
    qsub torus.pbs
done

echo "Submitting hybrid_db runs"
for bench in disc disc_cylindrical HII_region molebench; do
    cd ${base_dir}/${test_dir}/hybrid_db/benchmarks/disc
    ln -s ../../build/torus.zen
    write_qsub_hybrid
    qsub torus.pbs
done

echo "All done"
