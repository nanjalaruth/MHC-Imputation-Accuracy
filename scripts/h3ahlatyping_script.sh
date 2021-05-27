#!/bin/bash

#SBATCH --job-name='h3ahlatyping'
#SBATCH --cpus-per-task=2
#SBATCH --mem=5GB
#SBATCH --output=ggvphlatyping-%j-stdout.log
#SBATCH --error=ggvphlatyping-%j-stderr.log
#SBATCH --time=14-00:00:00

cd /scratch3/users/nanje/hlatyping/h3a_hlatypes

echo "Submitting SLURM job"
nextflow -c /users/nanje/Africanpop/h3a.config run nanjalaruth/hlatyping -profile singularity,slurm --input "/cbio/projects/013/custom-bam/*/*.bam" -resume
