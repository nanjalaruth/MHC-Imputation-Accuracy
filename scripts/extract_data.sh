#! /bin/bash

#liftover 38 to 37
#CrossMap.py vcf GRCh38_to_GRCh37.chain.gz v6.cgp.vf.va.chr6.filter-pass.vcf \
#/cbio/dbs/gatk/2.8/b37/human_g1k_v37_decoy.fasta h3a_all_hg38.vcf

#replace 6 with chr6 for hg38
#awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' h3ahg38.vcf | bcftools view -Oz -o h3ahg38.vcf.gz
#awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' h3_all_hg38.vcf | bcftools view -Oz -o h3a_all_hg38.vcf.gz
#tabix index the file
#tabix h3a_all_hg38.vcf.gz

#Keep 70% missing genotype
#bcftools view -i 'F_MISSING<0.3' h3a_annotated.vcf.gz >  h3a_clean.vcf.gz

#extract ids from vcf file
#bcftools query -l h3a_annotated.vcf.gz > queried.tsv
#head -n 2000 queried.tsv > samples.txt
#head -n 500 queried.tsv > smp.txt
#head -n 1000 queried.tsv > sp.txt

#extract samples from vcf file 
#bcftools view --samples-file h3a_samples.txt v4.chr6.minawi.ac2_phased_no_ref.vcf.gz --force-samples -Oz -o h3a_africans.vcf.gz
#bcftools view -c 1 -m2 -M2 -v snps --samples-file samples.txt h3a_annotated.vcf.gz -Oz -o h3a.vcf.gz
#bcftools view -c 1 -m2 -M2 -v snps --samples-file smp.txt h3a_annotated.vcf.gz -Oz -o lessh3a.vcf.gz
#bcftools view --samples-file samples.txt h3a_clean.vcf.gz -Oz -o newh3a.vcf.gz
#bcftools view -c 1 -m2 -M2 -v snps --samples-file sp.txt h3a_annotated.vcf.gz -Oz -o medh3a.vcf.gz


#extract chrom 6- MHC data only
#bcftools view --regions-file mhcfile.tsv h3a_africans.vcf.gz  -Oz -o h3a_africans_hla.vcf.gz
#bcftools view --regions-file mhcfile.tsv ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
#-Oz -o 1kg_All_hla.vcf.gz
#bcftools view --regions-file h3afile.tsv \ 
#/cbio/projects/001/ref-panels/v6/filter-vcf-genotype-refinement/v6.cgp.vf.va.chr6.filter-pass.vcf.gz \
#-Oz -o h3a.vcf.gz
#tabix h3a.vcf.gz

#bcftools view --regions-file ./mhcfile.tsv ./ggvp/ggvp.recal-SNP.recal-INDEL.filter-pass.vcf.gz \
# -Oz -o mhc_ggvp.vcf.gz

#cp /cbio/datasets/human/1kg/b37/low-coverage/vcf/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz .
#tabix ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
#bcftools view --regions-file mhcfile.tsv ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -Oz -o 1kg_chr6.vcf.gz

#annotate h3a snp ids
#download ref
#wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/00-All.vcf.gz
#gunzip 00-All.vcf.gz
#annotate using snpsift.jar
#java -jar snpEff/SnpSift.jar annotate -id 00-All.vcf h3a.vcf > h3a_annot.vcf

#bcftools annotate --set-id +'%CHROM\:%POS'  \
#/cbio/projects/001/ref-panels/v6/filter-vcf-genotype-refinement/v6.cgp.vf.va.chr6.filter-pass.vcf.gz \
#> h3a_all_annotated.vcf.gz

#bcftools annotate --set-id +'%CHROM\:%POS' h3a_annot.vcf > h3a_annotated.vcf 
#bgzip h3a_annotated.vcf
#tabix h3a_annotated.vcf.gz

#java -jar snpEff/SnpSift.jar annotate -id 00-All.vcf mhc_ggvp.vcf > mhc_ggvp_annotated.vcf

#extract data based on arrays
#bcftools view --regions-file ./arrays/h3Africa_v1_chip.pos ./mhc_ggvp_annotated.vcf.gz -Oz -o h3a_mhc_ggvp.vcf.gz
#tabix h3a_mhc_ggvp.vcf.gz
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

