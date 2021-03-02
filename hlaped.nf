#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {  combine_hlatypes_snp2hla; combine_hlatypes_hibag; subset_hlatypes_snp2hla; subset_hlatypes_hibag; get_geno_plink as get_geno_plink_dataset; get_geno_plink as get_geno_plink_subpop; make_snp2hlarefpanel as make_snp2hlarefpanel_dataset; make_snp2hlarefpanel as make_snp2hlarefpanel_subpop} from './modules/make_reference'

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


    // Make reference panel
    // African population
    make_snp2hlarefpanel_dataset(get_geno_plink_dataset.out)
    make_snp2hlarefpanel_subpop(get_geno_plink_subpop.out)

}

   
