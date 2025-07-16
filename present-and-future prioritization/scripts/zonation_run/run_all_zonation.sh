#!/bin/bash
#SBATCH --account=rrg-gonzalez
#SBATCH --time=02:30:00
#SBATCH --job-name=zonation-scenarios
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=40
#SBATCH --output=zonation_scenarios.out
#SBATCH --mail-user=glaroc@gmail.com
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --array=0-53
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load qt/4.8.7 fftw gdal boost
./run_one_timestep.sh $SLURM_ARRAY_TASK_ID
