#!/bin/bash

#SBATCH --job-name='HIBAG model for 1kg africans'
#SBATCH --mem=64GB
#SBATCH --output=trialafr-%j-stdout.log
#SBATCH --error=trialafr-%j-stderr.log
#SBATCH --time=10-00:00:00


echo "submitting slurm job"
Rscript h3aHIBAG.R
