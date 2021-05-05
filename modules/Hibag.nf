process hibag_impute {
    tag "Performing imputation on ${dataset}_${subpop} using HIBAG"
    publishDir "${params.outdir}"
    label "bigmem"
    
    input:
        tuple val(dataset), val(subpop), file(pbed), file(pbim), file(pfam), val(spop), file(bed), file(bim), file(fam), file(hlatyps), file(hla_a), file(hla_b), file(hla_c)
    output:
        tuple val(subpop), val(spop), file("${hibag_out}*")
    script:
        hibag_out = "${subpop}_${spop}"
        template "HIBAG.R"
}