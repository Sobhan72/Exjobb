#!/bin/sh
# requesting the number of nodes needed
#SBATCH -A lu2025-2-33
#SBATCH -N 2
#SBATCH --tasks-per-node=36

# job time, change for what your job farm requires
#SBATCH -t 08:00:00
#
# job name and output file names
#SBATCH -J jobFarm
#SBATCH -o tmp/master_%j.out
#SBATCH -e tmp/master_%j.err
cat $0

# set the number of jobs - change for your requirements
export NB_of_jobs=72

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

echo "Done with generating inputs"

# Copy necessary files
cp -p input*.mat job.m $MASTER_DIR
cp -p worker_script.sh $TMP_OUT

# Clean up generated input files in the working directory
rm -f input*.mat

# Change to master directory
cd $MASTER_DIR

# Maximum number of concurrent workers
MAX_PARALLEL=72
running=0

# Loop over the job number with throttling
for ((i=0; i<$NB_of_jobs; i++)); do
    echo "Starting job $i at $(TZ=Europe/Stockholm date +%T)"
    srun -Q --exclusive -n 1 -N 1 $TMP_OUT/worker_script.sh $i \
         &> $TMP_OUT/worker_${SLURM_JOB_ID}_${i}.out &

    ((running++))

    # if we've launched MAX_PARALLEL jobs, wait until one finishes
    if (( running >= MAX_PARALLEL )); then
        wait -n    # wait for *one* job to finish
        ((running--))
    fi
done

# wait for all remaining jobs to finish
wait

echo "Skickat alla jobb"

# Remove job.m after starting the jobs
rm -f job.m
