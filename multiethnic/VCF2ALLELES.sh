echo ID $(bcftools query -l chr6.dose.vcf.gz | tr '\n' ' ') > preprocessed.txt 
bcftools query -f '%ID[ %GT]\n' chr6.dose.vcf.gz  | grep '^HLA' >> preprocessed.txt 
python ../4d_VCF2ALLELES.py preprocessed.txt out
