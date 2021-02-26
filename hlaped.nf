#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// include { combine_hlatypes_snp2hla; subset_hlatypes_snp2hla; combine_hibaghlatypes; combine_gwdhibaghlatypes; afrsnpgenotypes; gwdsnpgenotypes; afr_snp2hlarefence} from './modules/make_reference'
include { combine_hibaghlatypes } from './modules/make_reference'
workflow{
    smp_ch = Channel.fromList(params.hlatype_files)
                    .map{ title, hlatype -> 
                        sub_hlatypes = file(hlatype)
                        return [title, sub_hlatypes] 
                    }
    combine_hibaghlatypes(smp_ch)
}

// HLA TYPES    
// SNP2HLA software
    // Combine afrhlatype files
    // samples_ch = Channel.fromList(params.hlatype_files)
    //                 .map{ name, hlatype -> 
    //                     hlatypes = file(hlatype)
    //                     return [name, hlatypes] 
    //                 }
    // combine_hlatypes_snp2hla(samples_ch)

//     // Get Gambian
    //  ids_ch = Channel.fromList(params.subpop_ids)
    //             .map{ subpop, hlatype ->
    //                 hlatypes = file(hlatype)
    //                 return [subpop, hlatypes]
    //             } 
    //             .combine(combine_hlatypes_snp2hla.out)
    // subset_hlatypes_snp2hla(ids_ch)

// HIBAG software
//     combine afrhlatype files
//     smp_ch = Channel.fromPath(params.files)
//     combine_afrhibaghlatypes(smp_ch.collect()) 

    

//     // Gambian subpopulation

//     id_ch = Channel.fromPath(params.gwdids).combine(combine_afrhibaghlatypes.out)
//     combine_gwdhibaghlatypes(id_ch) 

// // SNP genotypes
//     // African population
//     rds_ch = Channel.fromPath(params.reads).combine(combine_afrhlatypes.out)
//     afrsnpgenotypes(rds_ch)
    
//     //Gambian population
//     read_ch = Channel.fromPath(params.reads).combine(combine_gwdhlatypes.out)
//     gwdsnpgenotypes(read_ch)


// // Make reference
//     //African population
//     aref_ch = (afrsnpgenotypes.out).combine(combine_afrhlatypes.out).view()   
//     afr_snp2hlarefence(aref_ch)
// }