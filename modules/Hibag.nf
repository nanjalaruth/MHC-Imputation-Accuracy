process hibag_impute {
    tag "Performing imputation on ${dataset}_${subpop} using HIBAG"
    publishDir "${params.outdir}/Imputation/HIBAG", mode: 'copy', overwrite: false
    label "bigmem"
    
    input:
        tuple val(dataset), val(subpop), file(bed), file(bim), file(fam), val(spop), file(model_a), file(model_b), file(model_c), file(ans_file)
    output:
        tuple val(subpop), val(spop), file("${hibag_out}*"), file("${subpop}*")
    script:
        hibag_out = "${subpop}_${spop}"
        allele_out ="${subpop}"
        template "HIBAG.R"
}

process hibag_impute_eur {
    tag "Performing imputation on ${dataset}_${subpop} using HIBAG"
    publishDir "${params.outdir}/Imputation/HIBAG", mode: 'copy', overwrite: false
    label "medium"
    
    input:
        tuple val(dataset), val(subpop), file(bed), file(bim), file(fam), file(ans_file)
    output:
        tuple val(subpop), file("${hibag_out}*"), file("${subpop}*")
    script:
        hibag_out = "EUR_${subpop}"
        allele_out ="${subpop}"
        template "Eur_HIBAG.R"
}