#Convert to m3vcf with minimac3
minimac3 --refHaps Gwd_popref_hla.vcf.gz --processReference --prefix Gwd_popref_phased --rsid

#Imputation
minimac4 --refHaps Gwd_popref_phased.m3vcf.gz --haps omni_eagle_phased.vcf.gz \
 --prefix omni_imputed_minimac4_Gwd --cpus 10 --rsid --passOnly
