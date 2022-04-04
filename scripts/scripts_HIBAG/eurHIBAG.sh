#!/bin/bash

#SBATCH --job-name='HIBAG model for European population'
#SBATCH --mem=64GB
#SBATCH --output=eur-%j-stdout.log
#SBATCH --error=eur-%j-stderr.log
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --nodes=8

echo "submitting slurm job"
Rscript Eur_HIBAG.R 
