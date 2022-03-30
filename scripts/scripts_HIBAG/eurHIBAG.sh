#!/bin/bash

#SBATCH --job-name='HIBAG model for European population'
#SBATCH --mem=64GB
#SBATCH --output=eur-%j-stdout.log
#SBATCH --error=eur-%j-stderr.log
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --nodes=20

echo "submitting slurm job"
Rscript eurHIBAG.R 
