#!/bin/bash

#SBATCH --job-name='hlatyping'
#SBATCH --cpus-per-task=2
#SBATCH --mem=5GB
#SBATCH --output=1kghlatyping-%j-stdout.log
#SBATCH --error=1kghlatyping-%j-stderr.log
#SBATCH --time=4-00:00:00

cd /scratch3/users/nanje/hlatyping

echo "Submitting SLURM job"
nextflow -c /users/nanje/Africanpop/1kg-AFR.config run nanjalaruth/hlatyping -profile singularity,slurm --input "/users/nanje/data/link_Afrdata/*/*_R{1,2}.fastq.gz" -resume
