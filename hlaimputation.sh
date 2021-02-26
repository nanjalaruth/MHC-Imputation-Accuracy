#!/bin/bash

#SBATCH --job-name='hlaimputation'
#SBATCH --cpus-per-task=2
#SBATCH --mem=5GB
#SBATCH --output=hlaimputation-%j-stdout.log
#SBATCH --error=hlaimputation-%j-stderr.log
#SBATCH --time=5:00:00

cd /scratch3/users/nanje/hlatyping/results

echo "Submitting SLURM job"
nextflow -c /scratch3/users/nanje/MHC-Imputation-Accuracy/nextflow.config run /scratch3/users/nanje/MHC-Imputation-Accuracy/hlaped.nf -profile singularity,slurm -resume
