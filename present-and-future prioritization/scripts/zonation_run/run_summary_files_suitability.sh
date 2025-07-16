#!/bin/bash
#SBATCH --account=rrg-gonzalez
#SBATCH --time=02:00:00
#SBATCH --job-name=suitability_summaries
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=suitability_summaries
#SBATCH --mail-user=glaroc@gmail.com
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
module load nixpkgs/16.09  gcc/7.3.0 python gdal/3.0.1
parallel -j $SLURM_CPUS_PER_TASK /lustre03/project/6033499/lowlands/ConnectBTSL/Phase-III/zonation/projects/btsl/run/create_summary_files_suitability.sh ::: $(seq 0 4) ::: $(seq 5 8)