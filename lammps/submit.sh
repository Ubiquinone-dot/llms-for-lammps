#!/bin/bash -p
#$ -cwd
#$ -pe smp 128
#$ -j y
#$ -P cpu
#$ -l s_rt=24:00:00
#$ -o $JOB_ID.log

# lammps_build="/home/epsilon/vld/univ5120/VLD/amorphous-structure-of-sb2se3/lmp-25Jul2023-PLUMED"
lammps_build=/u/vld/univ5120/VLD/Thesis/llms-for-lammps/lmp
input_file=$1
rundir=$2
NMPI=8

################################################################

DIR=$rundir
module load aocc/3.2.0
module load aocl/3.2.0-aocc
module load mpi/openmpi-x86_64

function Cleanup ()
{
    trap "" SIGUSR1 EXIT SIGTERM SIGKILL # Disable trap now in it
    # Clean up task
    rsync -rlt $TMPDIR/* $DIR/
    exit 0
}
ulimit -s unlimited
trap Cleanup SIGUSR1 EXIT SIGTERM SIGKILL
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NMPI=$(expr $NSLOTS / $OMP_NUM_THREADS )

# Upload files to TMPDIR
INFILE='*' # copy infiles to the cluster (command line arguments here, modify as appropriate)
OUTFILE='*'
rsync -rlt $DIR/* $TMPDIR --exclude='*.log'
cd $TMPDIR # Use temporary directories to avoid i/o wastage from cluster to disk

################################################################


# Run the job
echo "Job started at $(date)\n\n"
echo "args: $1 \n $2"
mpiexec -np $NMPI $lammps_build -in $input_file 


################################################################
PID=$!
while kill -0 $PID 2> /dev/null; do
    rsync -rltq --exclude '*.sh' --exclude '*.in' --exclude '*.out' --exclude 'ompi*' --exclude 'log.lammps' $TMPDIR/ $DIR/
    #sleep 600
done
wait $PID
rsync -rltq --exclude '*.in' --exclude '*.out' --exclude 'ompi*' --exclude 'log.lammps' ./ $DIR
#sleep 200
cd $DIR
pwd
################################################################

echo "Job finished at $(date)\n\n"


