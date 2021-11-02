process preparedata_imputation {
    tag "get_samplegeno_plink_${dataset}_${subpop}"
    publishDir "${params.outdir}/Genotypes/SampleData", mode: 'copy', overwrite: false
    
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
            plink2 --vcf ${sample_geno} --make-bed --max-alleles 2 --out geno
            #remove duplicated snps
            cut -f 2 geno.bim | sort | uniq -d > 1.dups
            plink2 --bfile geno --exclude 1.dups --make-bed --out ${prefix}

        """
}


process impute_data {
    tag "Impute_${dataset}_${subpop} using ${spop} reference"
    publishDir "${params.outdir}/Imputation/SNP2HLA/PreImpute", mode: 'copy', overwrite: false
    
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
    publishDir "${params.outdir}/Imputation/SNP2HLA/Imputed", mode: 'copy', overwrite: false
    label "biggermem"
    
    input:
        tuple val(ref), val(array), file(unphased_bgl), file(markers), file(phased_bgl)
    output:
        tuple val(array), val(ref), file("${imputeoutput}.*")
    script:
        imputeoutput = "${array}_${ref}_IMPUTED"
        """
         java -Djava.io.tmpdir=${params.bgltemp} -Xmx${task.memory.toGiga()}g -jar /usr/local/bin/beagle.jar markers=${markers} unphased=${unphased_bgl} phased=${phased_bgl} gprobs=true niterations=10 nsamples=4 missing=0 verbose=true maxwindow=1000 lowmem=true seed=994000 out=${imputeoutput} log=${unphased_bgl}.phasing 
        """
}

process posteprob_dosage{
    tag "Convert posterior probability to dosage format"
    publishDir "${params.outdir}/Imputation/SNP2HLA/Imputed", mode: 'copy', overwrite: false

    input:
        tuple val(refenc), val(data), file(gprobs), file(phased), file(r2), file(fam), file(refbim)
    output:
        tuple val(data), val(refenc), file("${output}.*")
    script:
        output = "${data}_${refenc}_IMPUTED"
        """
        #Converting .gprobs to .dosage format
        mv ${phased} ${output}.bgl.phased
        mv ${gprobs} ${output}.bgl.gprobs
        mv ${r2} ${output}.bgl.r2

        /usr/local/bin/ParseDosage.csh ${output}.bgl.gprobs > ${output}.dosage

        #Converting imputation genotypes to PLINK .ped format."; @ i++
        cat ${output}.bgl.phased | java -jar /usr/local/bin/beagle2linkage.jar ${output}.tmp # Buhm
        cut -d ' ' -f6- ${output}.tmp.ped >  ${output}.tmp       # Buhm
        paste -d ' ' ${fam} ${output}.tmp | tr -d "\015" > ${output}.ped
        cut -f1-4 ${refbim} > ${output}.map
        cp ${fam} ${output}.fam

        #Create PLINK bed file
        plink --ped ${output}.ped --map ${output}.map --make-bed --out ${output}
        
        #remove unwanted files
        rm -f ${phased} 
        rm -f ${gprobs} 
        rm -f ${r2}
        """
}

process measureacc {
    tag "Measure accuracy for ${array} data using ${ref} reference"
    publishDir "${params.outdir}/Imputation/SNP2HLA/Accuracy", mode: 'copy', overwrite: false
    
    input:
        tuple file(answer_file), val(array), val(ref), file(bglphased)
    output:
        tuple val(array), val(ref), file("${output}.accuracy"), file("${output}_IMPUTED.hped")
    script:
        output = "${array}_${ref}"
        template "measacc.py"
}