#!/bin/sh
# requesting the number of nodes needed
#SBATCH -A lu2025-2-33
#SBATCH -N 1
#SBATCH --tasks-per-node=2
#SBATCH --mem-per-cpu=6200
#SBATCH --qos=test
#
# job time, change for what your job farm requires
#SBATCH -t 00:10:00
#
# job name and output file names
#SBATCH -J jobFarm
#SBATCH -o tmp/master_%j.out
#SBATCH -e tmp/master_%j.err
cat $0

# set the number of jobs - change for your requirements
export NB_of_jobs=2

# create master directory
export MASTER_DIR=JOB_${SLURM_JOB_ID}
mkdir $MASTER_DIR

cp -p input*.mat job.m *.sh $MASTER_DIR
cd $MASTER_DIR

# Loop over the job number

for ((i=0; i<$NB_of_jobs; i++))
do
    srun -Q --exclusive -n 1 -N 1 worker_script.sh $i &> worker_${SLURM_JOB_ID}_${i}.out &
    sleep 1
done

# keep the wait statement, it is important!

wait
