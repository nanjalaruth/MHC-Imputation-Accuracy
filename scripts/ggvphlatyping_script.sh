#!/bin/bash

#SBATCH --job-name='ggvphlatyping'
#SBATCH --cpus-per-task=2
#SBATCH --mem=5GB
#SBATCH --output=ggvphlatyping-%j-stdout.log
#SBATCH --error=ggvphlatyping-%j-stderr.log
#SBATCH --time=4-00:00:00

cd /scratch3/users/nanje/hlatyping/ggvp_hlatypes

echo "Submitting SLURM job"
nextflow -c /users/nanje/Africanpop/ggvp.config run nanjalaruth/hlatyping -profile singularity,slurm --input "/cbio/datasets/human/ggvp/b37/bam/*/*.bam" -resume
