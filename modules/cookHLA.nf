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
        """
            python -m /usr/local/bin/CookHLA/MakeGeneticMap \
            -i ${prefix} \
            -hg 18 \
            -ref ${ref}  \
            -o ${output}
        """
}


process cookHLAimpute {
    tag "CookHLA Imputation of ${subpop} datasets using ${spop} reference"
    publishDir "${params.outdir}/Imputation/CookHLA", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), val(subpop), file(dbed), file(dbim), file(dfam), val(spop), file(frqFRQ), file(bed), file(bim), file(fam), file(bglphased), file(markers), val(dtset), file(mach), file(erate)  
    output:
        tuple val(subpop), val(spop), file(alleles), file(hped)
    script:
        prefix = dbed.baseName
        ref = bed.baseName
        output = "${subpop}_${spop}_IMPUTED"  
        alleles = "${output}.alleles"
        hped = "${output}.hped"
        """
            python /usr/local/bin/CookHLA/CookHLA.py \
            -i ${prefix} \
            -hg 18 \
            -ref ${ref} \
            -o ${output} \
            -gm ${mach} \
            -ae ${erate} \
            -mem 120g \
            -mp 32   # The number of available cores for Multiprocessing.
            # -bgl4   # If you want to use Beagle4 instead of Beagle5.
        """
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
            python -m /usr/local/bin/CookHLA/measureAcc ${answer_file} ${alleles} ${output} 
        """
}