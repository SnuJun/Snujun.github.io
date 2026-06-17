#!/bin/bash

#SBATCH --job-name=LPD_AdaHard
#SBATCH --array=1-12
#SBATCH --cpus-per-task=1
##SBATCH --mem=4G
#SBATCH --time=2-00:00:00

#SBATCH --output=logs_revision_v2/LPD_%A_%a.out
#SBATCH --error=logs_revision_v2/LPD_%A_%a.err

cd ~/LPD_section_4_5/LPD_revision_v2

mkdir -p logs_revision_v2

module load R

Rscript main_relative_error_array.R ${SLURM_ARRAY_TASK_ID}
