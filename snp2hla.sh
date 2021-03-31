#!/bin/bash

#SBATCH --job-name='phasing'
#SBATCH --cpus-per-task=2
#SBATCH --mem=5GB
#SBATCH --output=phasing-%j-stdout.log
#SBATCH --error=phasing-%j-stderr.log
#SBATCH --time=24:00:00

cd /scratch3/users/nanje/hlatyping/results

echo "Submitting SLURM job"

java -Djava.io.tmpdir=./tempfiles -Xmx60g -jar ./beagle.jar unphased=./my-results/1kg_all_1kg_all_ref.bgl nsamples=4 niterations=10 missing=0 verbose=true maxwindow=1000 log=./my-results/1kg_all_1kg_all_ref.bgl.phasing
