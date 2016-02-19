#!/bin/bash
#SBATCH --time=16:00:00 --mem-per-cpu=3000
#SBATCH -o /triton/work/jaittoma/dippa-analyysi/log/CV-%a.out
#SBATCH --array=2

## Bash script for running Bayesian regression analysis models on Triton

mkdir -p /local/jaittoma
cd /local/jaittoma
#module load R/3.1.1-openblas
Rscript --vanilla /triton/work/jaittoma/dippa-analyysi/triton/search_cv.R $SLURM_ARRAY_TASK_ID
