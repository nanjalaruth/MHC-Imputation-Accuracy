#!/bin/bash

#SBATCH --job-name='generate info files'
#SBATCH --cpus-per-task=6
#SBATCH --mem=50GB
#SBATCH --output=info-%j-stdout.log
#SBATCH --error=info-%j-stderr.log
#SBATCH --time=100:00:00

echo "Submitting SLURM job"

for sample in `cat list.txt`
do
	#Convert results from plink to vcf 
	#plink --bfile /scratch3/users/nanje/hlatyping/results/my-results/Imputation/SNP2HLA/Imputed/${sample} \
	#	--recode vcf --out /cbio/projects/001/clients/hlaimpute/datasets/${sample}
	#annotate
	#bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
	#	 /cbio/projects/001/clients/hlaimpute/datasets/${sample}.vcf \
	#	> /cbio/projects/001/clients/hlaimpute/datasets/edited_${sample}.vcf
	#change ID names
	bcftools query -f '%ID\n' /cbio/projects/001/clients/hlaimpute/datasets/edited_${sample}.vcf \
		 > /cbio/projects/001/clients/hlaimpute/datasets/${sample}_info.txt
	#extract allele frequency
	#plink --freq --file /scratch3/users/nanje/hlatyping/results/my-results/Imputation/SNP2HLA/Imputed/${sample} --out /cbio/projects/001/clients/hlaimpute/datasets/${sample}
done
