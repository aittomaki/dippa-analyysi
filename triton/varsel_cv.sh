#!/bin/bash
#SBATCH --time=8:00:00 --mem-per-cpu=3000
#SBATCH -o /triton/work/jaittoma/dippa-analyysi/log/CV-%a.out
#SBATCH --array=8

## Bash script for running cross-validated variable selection on Triton

newgrp triton-users
cd $TMPDIR
module load R
Rscript $WRKDIR/dippa-analyysi/triton/search_cv.R $SLURM_ARRAY_TASK_ID

