#!/bin/bash

#SBATCH --job-name='HIBAG model for 1kg gambian'
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --output=afr-%j-stdout.log
#SBATCH --error=afr-%j-stderr.log
#SBATCH --time=10-00:00:00


echo "submitting slurm job"
Rscript 1kggwd_HIBAG.R
