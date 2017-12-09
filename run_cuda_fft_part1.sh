#!/bin/sh
#SBATCH --partition=slurm_shortgpu
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH -e job_err
#SBATCH -o job_out
#SBATCH --gres=gpu:1 # not needed for OpenMP

cd $SLURM_SUBMIT_DIR

for i in {1 2 4 8 16 32};
do 
./cuda_fft_part1 512 $i
done
 

mv job_out cuda_fft_part1.out
mv job_err cuda_fft_part1.err

