#!/bin/sh
#SBATCH --partition=slurm_shortgpu
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH -e job_err
#SBATCH -o job_out


cd $SLURM_SUBMIT_DIR

    ./single_thread_out 512 1

mv job_out single_thread.out
mv job_err single_thread.err

