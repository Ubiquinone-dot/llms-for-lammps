#!/bin/bash
#$ -cwd
#$ -j y
#$ -o main_$JOB_ID.log

# Compute resources
#$ -P cpu
#$ -l s_rt=720:00:00
#$ -l h_vmem=2G
#$ -pe smp 4
#$ -N main_agent

# script to submit to the cluster that will manage the deployed jobs
# better than running it locally as your computer may turn off etc.
# Will need to be able to run analysis script so don't give zero resources etc.
# Just don't forget about the job because, yknow, it'll be expensive.



################################################################
DIR=$(pwd)
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
# Run
export OPENAI_API_KEY=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

# conda activate .conda
conda activate /u/vld/univ5120/VLD/Thesis/llms-for-lammps/.conda
python main.py

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


echo "Experiment finished at $(date)"

