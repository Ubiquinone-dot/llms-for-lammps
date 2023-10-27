#!/bin/bash 
#$ -cwd
#$ -pe smp 8
#$ -l s_rt=100:00:00
#$ -l h_rss=150G
#$ -j y
#$ -o logs/$JOB_ID.log

# -pe smp = number of cores
# -l s_rt = soft run-time limit

###################################
#### create starting structure ####
###################################

if [ -z "$1" ]; then
    echo "No density passed."
    echo "Usage is: qsub submit.sh <density> <anneal_T>"
    echo "Exiting..."
    exit 1
else
    density=$1
fi

if [ -z "$2" ]; then
    echo "No anneal temperature passed."
    echo "Usage is: qsub submit.sh <density> <anneal_T>"
    echo "Exiting..."
    exit 1
else
    anneal_T=$2
fi

struct=$(./generate_structure.py N=200 density=$density t=$anneal_T)
echo Starting Structure: $struct

###################################
######## setup directories ########
###################################

DIR=$(pwd)
rundir=runs/$struct
mkdir -p ${rundir}/dumps
cp $DIR/structures/$struct.data ${rundir}

rsync -aq $rundir/ $TMPDIR/
cd $TMPDIR

###################################
########## load modules ###########
###################################

module load intel/2019
module load mpich3-intel

# some environment variables for parallelisation and memory usage
# LAMMPS mainly uses MPI parallelisation (at least with QUIP), so we 
# 'turn off' the OpenMP parallelisation by setting the number of threads
# to 1 (also for the intel Math Kernel Library (MKL) that handles matrix
# operations etc.)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NMPI=$(expr $NSLOTS / $OMP_NUM_THREADS )
export GFORTRAN_UNBUFFERED_ALL=y
ulimit -s unlimited

lmp_in=$DIR/lammps/infile
lmp_exec=$DIR/lammps/lmp_mpi

###################################
############# lammps ##############
###################################

# timesteps in ps
melt=5
quench=5
warm=50
anneal=100
cool=50

# temperatures in K
melt_T=9000
quench_T=300
# anneal_T=$anneal_T
cool_T=300

rand_seed=$(od -vAn -N4 -td4 < /dev/urandom | sed "s/-//")

pot=$DIR/lammps/carbon.xml

mpirun -np $NMPI $lmp_exec -in ${lmp_in} \
   -var quench $quench \
   -var warm $warm \
   -var anneal $anneal \
   -var cool $cool \
   -var melt $melt \
\
   -var melt_T ${melt_T} \
   -var quench_T ${quench_T} \
   -var anneal_T ${anneal_T} \
   -var cool_T ${cool_T} \
\
   -var structure $struct \
   -var rand_seed ${rand_seed} \
   -var pot $pot \
   -var model gap &

# copy back job data while we wait for it to finish
pid=$! 
while kill -0 $pid 2> /dev/null; do
    sleep 600  # copy every 600s
    rsync -aq $TMPDIR/ $DIR/$rundir/
done
wait $pid
