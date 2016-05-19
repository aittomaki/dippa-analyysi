#!/bin/bash
#SBATCH --time=4:00:00 --mem-per-cpu=3000
#SBATCH -o /scratch/work/jaittoma/dippa-analyysi/log/CV-%a.out
#SBATCH -p batch
#SBATCH --array=8

## Bash script for running cross-validated variable selection on Triton

#newgrp triton-users #opens new shell, not good
module load R/3.2.4-goolf-triton-2016a
Rscript --default-packages=methods,utils $WRKDIR/dippa-analyysi/triton/search_cv.R $SLURM_ARRAY_TASK_ID

