#extract chrom 6- MHC data only
bcftools view --regions-file /users/nanje/data/mhcfile.tsv /cbio/projects/001/ref-panels/v6/filter-vcf-genotype-refinement/v6.cgp.vf.va.chr6.filter-pass.vcf.gz  -Oz -o /cbio/users/nanje/h3a_vcf/h3a.vcf.gz

#annotate snp ids
#java -jar snpEff/SnpSift.jar annotate -id 00-All.vcf mhc_ggvp.vcf > mhc_ggvp_annotated.vcf

#extract data based on arrays
#bcftools view --regions-file ./arrays/h3Africa_v1_chip.pos ./mhc_ggvp_annotated.vcf.gz -Oz -o h3a_mhc_ggvp.vcf.gz
#bcftools view --regions-file ./arrays/omni.pos ./mhc_ggvp_annotated.vcf.gz -Oz -o omni_mhc_ggvp.vcf.gz

#annotate missing ids
#gunzip mhc_ggvp.vcf.gz 
#bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' mhc_ggvp.vcf > edited_mhc_ggvp.vcf
#bgzip edited_mhc_ggvp.vcf
#tabix edited_mhc_ggvp.vcf.gz

#gunzip h3a_mhc_ggvp.vcf.gz 
#bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' h3a_mhc_ggvp.vcf > edited_h3a_mhc_ggvp.vcf
#bgzip edited_h3a_mhc_ggvp.vcf
#tabix edited_h3a_mhc_ggvp.vcf.gz

#gunzip omni_mhc_ggvp.vcf.gz 
#bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' omni_mhc_ggvp.vcf > edited_omni_mhc_ggvp.vcf
#bgzip edited_omni_ mhc_ggvp.vcf
#tabix edited_omni_mhc_ggvp.vcf.gz

