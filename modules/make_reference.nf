#!/usr/bin/env nextflow
// nextflow.enable.dsl=2

process combine_hlatypes_snp2hla {
    tag "snp2hlatypes_${dataset}"
    publishDir "$params.outdir"

    input:
        tuple val(dataset), val(hlatype_list)
    output:
        tuple val(dataset), file(hlatypes_out)
    script:
        hlatypes_out = "${dataset}.snp2hla.ped"
        hlatype_files = hlatype_list.join(',')
        template "combine_hlatypes.py"
}

process subset_hlatypes_snp2hla {
    tag "subset_hlatypes_snp2hla_${subpop}_${dataset}"
    publishDir "$params.outdir"

    input:
        tuple val(dataset), val(subpop), file(subpop_ids), file(hlatypes)
    output:
        tuple val(dataset), val(subpop), file(hlatype_out)
    script:
        hlatype_out = "${subpop}.snp2hla.ped"
        """
            # extract gambian population
            grep -f ${subpop_ids} ${hlatypes} > ${hlatype_out}
        """
}

process combine_hlatypes_hibag {
    tag "combine_hlatypes_hibag_${dataset}"
    publishDir "$params.outdir"

    input:
        tuple val(dataset), val(hlatype_list)
    output:
        tuple val(dataset), file(hlatypes_out)
    script:
        hlatypes_out = "${dataset}_hibag_hlatypes"
        hlatype_files = hlatype_list.join(',')
        template "combine_hibaghlatypes.py"
}

process subset_hlatypes_hibag {
    tag "subset_hlatypes_hibag_${subpop}_${dataset}"
    publishDir "$params.outdir"

    input:
        tuple val(dataset), val(subpop), file(subpop_ids), file(hlatypes)
    output:
        tuple val(dataset), val(subpop), file(hlatype_out)
    script:
        hlatype_out = "${subpop}_hibag_HLAType"
        """
            # extract the gambian subpopulation
            grep -f ${subpop_ids} ${hlatypes} > subpop_ids_temp
            # remove the distorted first column
            awk '{OFS="\t";for(i=1;i<=NF;i++)if(i!=1)printf\$i OFS;print""}' subpop_ids_temp > t1.txt
            # introduce a new first column
            awk '{OFS="\t";print NR,\$0}' t1.txt > trial
            # introduce the headers
            ( echo -e "\tsample.id\tA.1\tA.2\tB.1\tB.2\tC.1\tC.2\tDQA1.1\tDQA1.2\tDQB1.1\tDQB1.2\tDRB1.1\tDRB1.2"; cat trial ) > ${hlatype_out}
        
        """
}

process get_geno_plink {
    tag "get_geno_plink_${dataset}_${subpop}"
    publishDir "${params.outdir}"
    
    input:
        tuple val(dataset), file(genotypes), val(subpop), file(hlatypes)
    output:
        tuple val(dataset), val(subpop), file(bed_out), file(bim_out), file(fam_out), file(hlatypes)
    script:
        prefix = "${dataset}_${subpop}_geno"
        bed_out = "${prefix}.bed"
        bim_out = "${prefix}.bim"
        fam_out = "${prefix}.fam"
        """
            #convert vcf to plink format
            plink2 --vcf ${genotypes} --make-bed --max-alleles 2 --out MHC
            #remove duplicated snps
            cut -f 2 MHC.bim | sort | uniq -d > 1.dups
            plink2 --bfile MHC --exclude 1.dups --make-bed --out MHC.filt
            #Get the samples that were typed
            cut -f 2 ${hlatypes} > ids.txt
            # extract plink datasets
            plink2 --bfile MHC.filt --keep ids.txt --make-bed --out ${prefix}
            rm MHC.* ids.txt ${prefix}.log
        """
}

process make_snp2hlarefpanel {
    tag "make_snp2hlarefpanel_${dataset}_${subpop}"
    publishDir "$params.outdir"
    
    input:
        tuple val(dataset), val(subpop), file(bed), file(bim), file(fam), file(hla)
    output:
        tuple val(dataset), val(subpop), file("${refoutput}.*")
    script:
        base = bed.baseName
        refoutput = "${dataset}_${subpop}_ref"
        template "MakeReference.csh"

}

process convert_to_beagle {
    tag "convert_to_beagle_${dataset}_${subpop}"
    publishDir "$params.outdir"
    label "medium"

    input:
        tuple val(dataset), val(subpop), file(dat), file(nopheno_ped)
    output:
        tuple val(dataset), val(subpop), file(bgl)
    script:
        refoutput = "${dataset}_${subpop}_ref"
        dat = "${refoutput}.dat"
        nopheno_ped = "${refoutput}.nopheno.ped"
        bgl = "${refoutput}.bgl"
        """
        java -Djava.io.tmpdir=${params.bgltemp} -Xmx${task.memory.toGiga()}g -jar /usr/local/bin/linkage2beagle.jar ${dat} ${nopheno_ped} > ${bgl} 
        
        """
}

process makeref_snp2hla_phasing {
    tag "makeref_snp2hla_phasing_${dataset}_${subpop}"
    publishDir "$params.outdir"
    label "bigmem"
    
    input:
        tuple val(dataset), val(subpop), file(bgl)
    output:
        tuple val(dataset), val(subpop), file("${bgl}.phased")
    script:
        """
        java -Djava.io.tmpdir=${params.scratch} -Xmx${task.memory.toGiga()}g -jar /usr/local/bin/beagle.jar unphased=${bgl} nsamples=4 niterations=10 missing=0 seed=993470 lowmem=true verbose=true maxwindow=1000 log=${bgl}.phasing
        """
}