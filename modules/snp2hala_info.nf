process combine_output {
    tag "Generating snp2hla info file for ${target}_${reference} "
    publishDir "${params.outdir}/SNP2HLA_INFO", mode: 'copy', overwrite: false

     input:
        tuple val(target), val(reference), file(imputed_bed), file(imputed_bim), file(imputed_fam), file(imputed_map), file(imputed_ped), file(imputed_r2)
    output:
        tuple val(target), val(reference), file(imputed_r2), file("${frq}.frq"), file(id)
    script:
        frq = "${target}_${reference}"
        id = "${target}_${reference}_id.txt"
        prefix = imputed_bed.baseName
        """
            #Convert results from plink to vcf
            plink --bfile ${prefix} --recode vcf --out ${prefix}
            #annotate
            bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' ${prefix}.vcf > edited_${prefix}.vcf
            #change ID names
            bcftools query -f '%ID\\n' edited_${prefix}.vcf > ${id}

            #compute for MAF
            plink --freq --file ${prefix} --out ${frq}

            #cleanup
            rm -f ${prefix}.vcf edited_${prefix}.vcf
        """
}


process generate_info {
    tag "Generating snp2hla info file for ${target}_${reference} "
    publishDir "${params.outdir}/SNP2HLA_INFO", mode: 'copy', overwrite: false
    
    input:
        tuple val(target), val(reference), file(rsquared), file(frq), file(snpid)
    output:
        tuple val(target), val(reference), file(info_out)
    script:
        info_out = "${target}_${reference}_info.txt"
        template "snp2hla_info.R"
}