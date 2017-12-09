#!/bin/sh
#SBATCH --partition=slurm_shortgpu
#SBATCH --nodes=1
#SBATCH --tasks=40
#SBATCH --cpus-per-task=1
#SBATCH -e job_err
#SBATCH -o job_out
#SBATCH --gres=gpu:1 # not needed for OpenMP

cd $SLURM_SUBMIT_DIR

    ./problem1a 

mv job_out problem1c.out
mv job_err problem1c.err

