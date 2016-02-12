#!/bin/bash
#SBATCH --time=0-01:00:00 --mem-per-cpu=1000
#SBATCH -p short
#SBATCH -o log/regression-%a.out
#SBATCH --array=1

## Bash script for running Bayesian regression analysis models on Triton

#module load R/3.1.1-openblas
Rscript --vanilla regression_multi_triton.R $SLURM_ARRAY_TASK_ID
