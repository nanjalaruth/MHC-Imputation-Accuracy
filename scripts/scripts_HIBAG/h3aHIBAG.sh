#!/bin/bash

#SBATCH --job-name='HIBAG model for h3a population'
#SBATCH --mem=90GB
#SBATCH --cpus-per-task=32
#SBATCH --output=h3a-%j-stdout.log
#SBATCH --error=h3a-%j-stderr.log
#SBATCH --time=14-00:00:00


echo "submitting slurm job"
Rscript h3aHIBAG.R
