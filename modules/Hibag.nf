process hibag_impute {
    tag "Performing imputation on ${dataset}_${subpop} using HIBAG"
    publishDir "${params.outdir}"
    label "bigmem"
    
    input:
        tuple val(dataset), val(subpop), file(bed), file(bim), file(fam), file(hlatyps)
    output:
        tuple val(dataset), val(subpop), file("${hibag_out}*")
    script:
        hibag_out = "${dataset}_${subpop}"
        template "HIBAG.R"
        // prefix = "${dataset}_${subpop}_geno"
        // bed = "${prefix}.bed"
        // fam = "${prefix}.fam"
        // bim = "${prefix}.bim"
        // hlatyps = "${subpop}_hibag_HLAType"
}