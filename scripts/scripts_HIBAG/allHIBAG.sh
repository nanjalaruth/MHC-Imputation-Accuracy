#!/bin/bash

#SBATCH --job-name='HIBAG model for 1000 genome population'
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=32
#SBATCH --output=trial_all-%j-stdout.log
#SBATCH --error=trial_all-%j-stderr.log
#SBATCH --time=14-00:00:00


echo "submitting slurm job"
Rscript trialall_HIBAG.R
