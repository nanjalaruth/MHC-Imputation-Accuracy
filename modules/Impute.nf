process impute_data {
    tag "Impute_${dataset}_${subpop}"
    publishDir "$params.outdir"
    
    input:
        tuple val(dataset), val(subpop), file(bed), file(bim), file(fam), file(markers)
        file(bgl_phased)
    output:
        tuple val(dataset), val(subpop), file("${imputeoutput}.*")
    script:
        base = bed.baseName
        refoutput = "${dataset}_${subpop}_ref"
        markers = "${refoutput}.markers"
        ref = markers.baseName
        imputeoutput = "${dataset}_${subpop}_imputed"
        template "SNP2HLA.csh"

}