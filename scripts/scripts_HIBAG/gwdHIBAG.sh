#!/bin/bash

#SBATCH --job-name='HIBAG model for 1kg gambian'
#SBATCH --mem=64GB
#SBATCH --output=gwd-%j-stdout.log
#SBATCH --error=gwd-%j-stderr.log
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --nodes=20

echo "submitting slurm job"
Rscript 1kggwd_HIBAG.R 
