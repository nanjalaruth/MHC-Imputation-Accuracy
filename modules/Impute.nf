process preparedata_imputation {
    tag "get_samplegeno_plink_${dataset}_${subpop}"
    publishDir "${params.outdir}"
    
    input:
        tuple val(dataset), val(subpop), file(sample_geno)
    output:
        tuple val(dataset), val(subpop), file(bed_out), file(bim_out), file(fam_out)
    script:
        prefix = "${dataset}_${subpop}_geno"
        bed_out = "${prefix}.bed"
        bim_out = "${prefix}.bim"
        fam_out = "${prefix}.fam"
        """
            #convert vcf to plink format
            plink2 --vcf ${sample_geno} --make-bed --max-alleles 2 --out ${prefix}

        """
}


process impute_data {
    tag "Impute_${dataset}_${subpop} using ${spop} reference"
    publishDir "$params.outdir"
    
    input:
        tuple val(dataset), val(subpop), file(bed), file(bim), file(fam), val(spop), file(pbed), file(pbim), file(pfam)
    output:
        tuple val(dataset), val(subpop), val(spop), file("${imputeoutput}.*")
    script:
        base = bed.baseName
        ref = pbed.baseName
        imputeoutput = "${dataset}_${subpop}_${spop}"
        template "SNP2HLA.csh"

}

process imputation{
    tag "Impute data from the ${array} array using ${ref} reference"
    publishDir "$params.outdir"
    label "bigmem"
    
    input:
        tuple val(ref), val(array), file(unphased_bgl), file(markers), file(phased_bgl)
    output:
        tuple val(array), val(ref), file(imputeoutput)
    script:
        imputeoutput = "${array}_${ref}.IMPUTED"
        """
         java -Djava.io.tmpdir=${params.bgltemp} -Xmx${task.memory.toGiga()}g -jar /usr/local/bin/beagle.jar markers=${markers} unphased=${unphased_bgl} phased=${phased_bgl} gprobs=true niterations=10 nsamples=4 missing=0 verbose=true maxwindow=1000 lowmem=true seed=994000 out=${imputeoutput} log=${unphased_bgl}.phasing 
        """
}