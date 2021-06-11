process makegenetic_map {
    tag "Make genetic map for ${subpop} data using ${spop} reference"
    publishDir "${params.outdir}/Imputation/CookHLA", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), val(subpop), file(dbed), file(dbim), file(dfam), val(spop), file(frqFRQ), file(bed), file(bim), file(fam), file(bglphased), file(markers)
    output:
        tuple val(subpop), val(spop), file(erate), file(mach)  
    script:
        prefix = dbed.baseName
        ref = bed.baseName
        output = "${subpop}_${spop}_GENMAP" 
        erate = "${output}.aver.erate"
        mach = "${output}.mach_step.avg.clpsB" 
        template "makegenmap.py"
}


process cookHLAimpute {
    tag "CookHLA Imputation of ${subpop} datasets using ${spop} reference"
    publishDir "${params.outdir}/Imputation/CookHLA", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), val(subpop), file(dbed), file(dbim), file(dfam), val(spop), file(frqFRQ), file(bed), file(bim), file(fam), file(bglphased), file(markers) 
    output:
        tuple val(subpop), val(spop), file(alleles), file(hped)
    script:
        prefix = dbed.baseName
        ref = bed.baseName
        output = "${subpop}_${spop}_IMPUTED"  
        alleles = "${output}.alleles"
        hped = "${output}.hped"
        template "cook.py"
}

process measureaccuracy {
    tag "Measure accuracy for ${array} data using ${ref} reference"
    publishDir "${params.outdir}/Imputation/CookHLA", mode: 'copy', overwrite: false
    
    input:
        tuple file(answer_file), val(array), val(ref), file(alleles)
    output:
        tuple val(array), val(ref), file("${output}.*") 
    script:
        output = "${array}_${ref}_ACCURACY"
        """
            cd /usr/local/bin
            python -m measureAcc ${answer_file} ${alleles} ${output} 
        """
}