#!/bin/sh
# requesting the number of nodes needed
#SBATCH -A lu2025-2-33
#SBATCH -N 1
#SBATCH --tasks-per-node=2
#SBATCH --mem-per-cpu=6200
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

# Get the absolute path of the current working directory
export WORK_DIR=$PWD

# Create master directory and tmp folder as absolute paths
export MASTER_DIR=$WORK_DIR/JOB_${SLURM_JOB_ID}
mkdir -p $MASTER_DIR

export TMP_OUT=$MASTER_DIR/tmp
mkdir -p $TMP_OUT

module load matlab/2024b

# Run MATLAB to generate input files
matlab -singleCompThread -nodesktop -nodisplay -nosplash -r "genInput; quit"

# Copy necessary files
cp -p input*.mat job.m $MASTER_DIR
cp -p worker_script.sh $TMP_OUT

# Clean up generated input files in the working directory
rm -f input*.mat

# Change to master directory
cd $MASTER_DIR

# Loop over the job number
for ((i=0; i<$NB_of_jobs; i++))
do
    srun -Q --exclusive -n 1 -N 1 $TMP_OUT/worker_script.sh $i &> $TMP_OUT/worker_${SLURM_JOB_ID}_${i}.out &
    sleep 1
done

# Remove job.m after starting the jobs
rm -f job.m

# Keep the wait statement; it is important
wait