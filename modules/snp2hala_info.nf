process annotate_vcf {
    tag "Generating snp2hla info file for ${target}_${reference} "
    publishDir "${params.outdir}/SNP2HLA_INFO", mode: 'copy', overwrite: false

     input:
        tuple val(target), val(reference), file(imputed_bed), file(imputed_bim), file(imputed_fam), file(imputed_map), file(imputed_ped), file(imputed_r2)
    output:
        tuple val(target), val(reference), file(imputed_r2), file("${prefix}.txt"), file("${prefix}.vcf.gz")
    script:
        prefix = imputed_bed.baseName
        """
            #Convert results from plink to vcf
            plink --bfile ${prefix} --recode vcf --out ${prefix}
            bcftools +fill-tags ${prefix}.vcf -Oz -o ${prefix}.vcf.gz
            bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%FIRST_ALT\\t%INFO/MAF\\n' -i 'AC > 0' ${prefix}.vcf.gz > ${prefix}.txt
        """
}

process generate_info {
    tag "Generating snp2hla info file for ${target}_${reference}"
    publishDir "${params.outdir}/SNP2HLA_INFO", mode: 'copy', overwrite: false
    
    input:
        tuple val(target), val(reference), file(rsquared), file(ann)
    output:
        tuple val(target), val(reference), file(info_out)
    script:
        info_out = "${target}_${reference}_info.txt"
        template "snp2hla_info.py"
}

// process combine_output {
//     tag "Generating snp2hla info file for ${target}_${reference} "
//     publishDir "${params.outdir}/SNP2HLA_INFO", mode: 'copy', overwrite: false

//      input:
//         tuple val(target), val(reference), file(imputed_bed), file(imputed_bim), file(imputed_fam), file(imputed_map), file(imputed_ped), file(imputed_r2)
//     output:
//         tuple val(target), val(reference), file(imputed_mhc_r2), file(frq), file(id), file("${prefix}.vcf")
//     script:
//         frq = "${target}_${reference}.frq"
//         id = "${target}_${reference}_id.txt"
//         imputed_mhc_r2 = "${target}_${reference}.r2"
//         prefix = imputed_bed.baseName
//         """
//             #Extract imputed SNPs
//             grep "SNP"  ${prefix}.bim >  ${prefix}_SNPs.bim
//             cut -f2  ${prefix}_SNPs.bim > ${prefix}_SNPs.txt
//             #Extract imputed SNPs from SNP file
//             plink --bfile  ${prefix} --extract  ${prefix}_SNPs.txt --make-bed --out  ${prefix}_test
//             #Convert results from plink to vcf
//             plink --bfile ${prefix}_test --recode vcf --out ${prefix}_orig
//             bcftools +fill-tags ${prefix}_orig.vcf -o ${prefix}.vcf -- -t AF,AN,AC,NS
//             #annotate
//             bcftools annotate --set-id '%CHROM\\:%POS\\:%REF\\:%FIRST_ALT' ${prefix}.vcf > edited_${prefix}.vcf
//             #change ID names
//             bcftools query -f '%ID\\n' edited_${prefix}.vcf > ${id}

//             #compute for MAF
//             plink --freq --file ${prefix} --out ${prefix}
//             grep -f ${prefix}_SNPs.txt ${prefix}.frq > ${frq}


//             #Extract MHC r2
//             grep -f ${prefix}_SNPs.txt ${imputed_r2} > ${imputed_mhc_r2}

//         """
// }
