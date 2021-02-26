#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process combine_hibaghlatypes {
    tag "hibaghlatypes_${ref_panel}"
    publishDir "$params.outdir"

    input:
        tuple val(ref_panel), val(hlatype_list)
    output:
        file(hlatypes_out)
    script:
        hlatypes_out = "${ref_panel}.hibag_hlatypes"
        hlatype_files = hlatype_list.join(',')
        template "combine_hibaghlatypes.py"
}
   

// process combine_hlatypes_snp2hla {
//     tag "hlatypes_${ref_panel}"
//     publishDir "$params.outdir"

//     input:
//         tuple val(ref_panel), val(hlatype_list)
//     output:
//         file(hlatypes_out)
//     script:
//         hlatypes_out = "${ref_panel}.snp2hla.ped"
//         hlatype_files = hlatype_list.join(',')
//         template "combine_hlatypes.py"
// }

// process subset_hlatypes_snp2hla {
//     tag "hlatypes_${pop}"
//     publishDir "$params.outdir"

//     input:
//         tuple val(pop), file(gwdid), file(afrhlatypes)
//     output:
//         file(gwdhlatype_out)
//     script:
//         gwdhlatype_out = "${pop}.snp2hla.ped"
//         """
//             # extract gambian population
//             grep -f ${gwdid} ${afrhlatypes} > ${gwdhlatype_out}
//         """
// }


// process combine_gwdhibaghlatypes {
//     tag "combine_afrhibaghlatypes"
//     publishDir "$params.outdir"

//     input:
//         tuple file(h_gwdid), file(h_afrhlatypes)
//     output:
//         file(h_gwdhlatype_out)
//     script:
//         h_gwdhlatype_out = "gwdhibag_HLAType"
//         """
//             # extract the gambian subpopulation
//             grep -f ${h_gwdid} ${h_afrhlatypes} > gwd_hi
//             # remove the distorted first column
//             awk '{OFS="\t";for(i=1;i<=NF;i++)if(i!=1)printf\$i OFS;print""}' gwd_hi > t1.txt
//             # introduce a new first column
//             awk '{OFS="\t";print NR,\$0}' t1.txt > trial
//             # introduce the headers
//             ( echo -e "\tsample.id\tA.1\tA.2\tB.1\tB.2\tC.1\tC.2\tDQA1.1\tDQA1.2\tDQB1.1\tDQB1.2\tDRB1.1\tDRB1.2"; cat trial ) > ${h_gwdhlatype_out}
        
//         """
// }

// process afrsnpgenotypes {
//     tag "afrsnpgenotypes"
//     publishDir "$params.outdir"
    
//     input:
//         tuple file(rds), file(a_hlatypes)
//     output:
//         file(agenotypes_out)
//     script:
//         agenotypes_out = "AFR.*"
//         """
//             #convert vcf to plink format
//             plink2 --vcf ${rds} --make-bed --max-alleles 2 --out MHC
//             #remove duplicated snps
//             cut -f 2 MHC.bim | sort | uniq -d > 1.dups
//             plink2 --bfile MHC --exclude 1.dups --make-bed --out MHC.filt
//             #Get the samples that were typed
//             cut -f 1 ${a_hlatypes} > AFRids.txt
//             #extract AFR plink datasets
//             plink2 --bfile MHC.filt --keep /users/nanje/results/AFRids.txt --make-bed --out AFR
//             rm AFR.log

//         """
// }

// process gwdsnpgenotypes {
//     tag "gwdsnpgenotypes"
//     publishDir "$params.outdir"
    
//     input:
//         tuple file(rd), file(g_hlatypes)
//     output:
//         file(ggenotypes_out)
//     script:
//         ggenotypes_out = "GWD.*"
//         """
//             #convert vcf to plink format
//             plink2 --vcf ${rd} --make-bed --max-alleles 2 --out MHC
//             #remove duplicated snps
//             cut -f 2 MHC.bim | sort | uniq -d > 1.dups
//             plink2 --bfile MHC --exclude 1.dups --make-bed --out MHC.filt
//             #Get the samples that were typed
//             cut -f 1 ${g_hlatypes} > GWDids.txt
//             #extract GWD plink datasets
//             plink2 --bfile MHC.filt --keep /users/nanje/results/GWDids.txt --make-bed --out GWD
//             rm GWD.log

//         """
// }


// process afr_snp2hlarefence {
//     publishDir "$params.outdir"
    
//     input:
//         set file(bed), file(bim), file(fam), file(hla)
//     output:
//         file(afr_snp2hlarefence_out)
//     script:
//         afr_snp2hlarefence_out = "*"
//         base = bed.baseName
//         refoutput = "aref"
//         """
//             MakeReference.csh ${base} ${hla} ${refoutput} plink
//         """
// }
