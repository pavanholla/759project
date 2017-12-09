#!/bin/sh
#SBATCH --partition=slurm_shortgpu
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH -e job_err
#SBATCH -o job_out


cd $SLURM_SUBMIT_DIR

for i in {1 2 4 8 16 32};
do 
./openmp_out 512 $i
done

mv job_out open_mp.out
mv job_err open_mp.err

