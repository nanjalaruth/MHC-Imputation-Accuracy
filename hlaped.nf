#!/usr/bin/env nextflow

// Enable DSL syntax
nextflow.enable.dsl=2

// Include modules
include {  combine_hlatypes_snp2hla; combine_hlatypes_hibag; 
subset_hlatypes_snp2hla; subset_hlatypes_hibag; get_geno_plink as get_geno_plink_dataset; 
get_geno_plink as get_geno_plink_subpop; make_snp2hlarefpanel as make_snp2hlarefpanel_dataset; 
make_snp2hlarefpanel as make_snp2hlarefpanel_subpop; 
makeref_snp2hla_phasing as makeref_snp2hla_phasing_dataset; 
makeref_snp2hla_phasing as makeref_snp2hla_phasing_subpop; 
convert_to_beagle as convert_to_beagle_dataset; 
convert_to_beagle as convert_to_beagle_subpop} from './modules/make_reference'

include { impute_data as impute_dataset; impute_data as impute_subpop} from './modules/Impute.nf'

include { hibag_impute as hibag_impute_dataset; hibag_impute as hibag_impute_subpop } from './modules/Hibag.nf'

// Main workflow
workflow{

    // Step 1. Getting HLA TYPES 
    // Step 1.1. For SNP2HLA
    // Step 1.1.1. SNP2HLA Africans
    samples_ch = Channel.fromList(params.hlatype_files)
                    .map{ name, hlatype -> 
                        hlatypes = file(hlatype)
                        return [name, hlatypes] 
                    }
    combine_hlatypes_snp2hla(samples_ch)


    // Step 1.1.2. SNP2HLA Subpopulation (Gambian, ....
    subpop_ids_ch = Channel.fromList(params.subpop_ids)
                .map{ dataset, subpop, ids  ->
                    ids_file = file(ids)
                    return [dataset, subpop, ids_file]
                } 
                .combine(combine_hlatypes_snp2hla.out, by:0)
    subset_hlatypes_snp2hla(subpop_ids_ch)


    // Step 1.2. For HIBAG
    // Step 1.2.1. HIBAG African
    smp_ch = Channel.fromList(params.hlatype_files)
                    .map{ title, hlatype -> 
                        hlatypes = file(hlatype)
                        return [title, hlatypes] 
                    }
    combine_hlatypes_hibag(smp_ch)

    // Step 1.2.2. HIBAG_Gambian
    subpop_id_ch = Channel.fromList(params.subpop_ids)
                .map{ dataset, subpop, ids -> 
                    id_file = file(ids)
                    return [dataset, subpop, id_file] 
                }
                .combine(combine_hlatypes_hibag.out, by:0)
    subset_hlatypes_hibag(subpop_id_ch) 


    // Step 2. Getting SNP GENOTYPES
    // Step 2.1. Dataset: 1KG African population, H3A
    dataset_geno_ch = Channel.fromList(params.genotype_files)
                .map{ dataset, genotype -> 
                    genotypes = file(genotype)
                    return [dataset, genotypes] 
                }
                .combine(combine_hlatypes_snp2hla.out, by:0)
                .map{ dataset, genotypes, hlatypes -> 
                    return [dataset, genotypes, dataset, hlatypes] 
                }
    get_geno_plink_dataset(dataset_geno_ch)

    // // //Step 2.2. Subpopulation population: Gambian, Zambiam
    subpop_geno_ch = Channel.fromList(params.genotype_files)
                    .map{ dataset, genotype -> 
                        genotypes = file(genotype)
                        return [dataset, genotypes] 
                    }.combine(subset_hlatypes_snp2hla.out, by:0)
    get_geno_plink_subpop(subpop_geno_ch)

    //SNP2HLA
    // Step 3: Make reference panel
    // Step 3.1: Get first output files
    make_snp2hlarefpanel_dataset(get_geno_plink_dataset.out)
    make_snp2hlarefpanel_subpop(get_geno_plink_subpop.out)

    // Step 3.2 Convert to Beagle
    // Step 3.2.1: Dataset
    beagle_data_in_ch = make_snp2hlarefpanel_dataset.out
    convert_to_beagle_dataset(beagle_data_in_ch)

    // Step 3.2.2: subpop
    beaglein_ch = make_snp2hlarefpanel_subpop.out
    convert_to_beagle_subpop(beaglein_ch)


    // Step 3.3: Phasing
    // Step 3.3.1: Dataset
    makeref_snp2hla_phasing_dataset(convert_to_beagle_dataset.out)

    // Step 3.3.2: Subpop
    makeref_snp2hla_phasing_subpop(convert_to_beagle_subpop.out)

    // Step 4: Imputation
    //Step 4.1: Dataset
    // dataset_testdata_ch = Channel.fromList(params.sample_data)
    //                 .map{ dataset, bed, bim, fam -> [dataset, bed, bim, fam]}
    //                 .combine(make_snp2hlarefpanel_subpop.out, by:0)
    // impute_dataset(dataset_testdata_ch, makeref_snp2hla_phasing_subpop.out)

    //Step 4.2: Subpop
    // subpop_testdata_ch = Channel.fromList(params.sample_data)
    //                 .map{ dataset, genotype -> 
    //                     genotypes = file(genotype)
    //                     return [dataset, genotypes] 
    //                 }.combine(make_snp2hlarefpanel_subpop.out, by:0)
    // impute_subpop(subpop_testdata_geno_ch, makeref_snp2hla_phasing_subpop.out)     


    //HIBAG
    input_ch = get_geno_plink_subpop.out.combine(subset_hlatypes_hibag.out, by:[0,1])
    .map{dataset, subpop, bed, bim, fam, snp2hla_ped, hibag_HLAType -> [dataset, subpop, bed, bim, fam, hibag_HLAType]} 
    hibag_impute_subpop(input_ch)


    input_two_ch = get_geno_plink_dataset.out.combine(combine_hlatypes_hibag.out, by:[0])
    .map{dataset, subpop, bed, bim, fam, snp2hla_ped, hibag_HLAType -> [dataset, subpop, bed, bim, fam, hibag_HLAType]} 
    hibag_impute_dataset(input_two_ch)

} 