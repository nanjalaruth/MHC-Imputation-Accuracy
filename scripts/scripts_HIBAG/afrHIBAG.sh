#!/bin/bash

#SBATCH --job-name='HIBAG model for 1kg africans'
#SBATCH --mem=64GB
#SBATCH --output=afr-%j-stdout.log
#SBATCH --error=afr-%j-stderr.log
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --nodes=20

echo "submitting slurm job"
Rscript afr_HIBAG.R

