#split multiallelic sites
bcftools norm -m- /cbio/projects/001/clients/hlaimpute/datasets/h3a_mhc_ggvp.vcf.gz -Ob -o h3a_splitted_vcf

#phase with Eagle
./eagle --numThreads=16 --vcf="./h3a_splitted_vcf" \
--geneticMapFile="/cbio/dbs/refpanels/eagle/tables/genetic_map_hg19_withX.txt.gz" \
--vcfOutFormat=z --outPrefix="h3a_eagle_phased"

#Convert to m3vcf with minimac3
minimac3 --refHaps Gwd_popref_hla.vcf.gz --processReference --prefix Gwd_popref_phased --rsid

#Imputation
minimac4 --refHaps Gwd_popref_phased.m3vcf.gz --haps omni_eagle_phased.vcf.gz \
 --prefix omni_imputed_minimac4_Gwd --cpus 10 --rsid --passOnly
