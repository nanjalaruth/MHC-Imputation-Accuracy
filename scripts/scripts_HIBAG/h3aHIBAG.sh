#!/bin/bash

#SBATCH --job-name='HIBAG model for H3A population'
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=32
#SBATCH --output=h3a-%j-stdout.log
#SBATCH --error=h3a-%j-stderr.log
#SBATCH --time=6-00:00:00
#SBATCH --nodes=20


echo "submitting slurm job"
Rscript h3aHIBAG.R
