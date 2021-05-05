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

include { preparedata_imputation; impute_data as impute_dataset; 
impute_data as impute_subpop; imputation as subpop_imputation; 
imputation as dataset_imputation; posteprob_dosage as posteprob_dosage_subpop;
posteprob_dosage as posteprob_dosage_dataset} from './modules/Impute.nf'

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
    beagle_input = make_snp2hlarefpanel_dataset.out.map{dataset, subpop, datas -> [dataset, subpop, datas[11], datas[15]]}
    convert_to_beagle_dataset(beagle_input)

    // Step 3.2.2: subpop
    beagle_subpop_input = make_snp2hlarefpanel_subpop.out.map{dataset, subpop, datas -> [dataset, subpop, datas[11], datas[15]]}
    convert_to_beagle_subpop(beagle_subpop_input)


    // Step 3.3: Phasing
    // Step 3.3.1: Dataset
    makeref_snp2hla_phasing_dataset(convert_to_beagle_dataset.out)

    // Step 3.3.2: Subpop
    makeref_snp2hla_phasing_subpop(convert_to_beagle_subpop.out)

    // Step 4: Imputation
    // Step 4.1: prepare data
    dataset_testdata_ch = Channel.fromList(params.sample_data)
                    .map{ pop, chip, sample_geno -> 
                        genotype = file(sample_geno)
                        return [pop, chip, genotype]
                    }
    // dataset_testdata_ch.view()   
    preparedata_imputation(dataset_testdata_ch)

    // Step 4.2: Pre-imputation
    // Step 4.2.1: Subpop
    bgl_input = make_snp2hlarefpanel_subpop.out.map{dataset, subpop, datas -> [dataset, subpop, datas[9], datas[10], datas[12]]}
    impute_input = preparedata_imputation.out.combine(bgl_input)
    .map {dataset, subpop, bed, bim, fam, pop, spop, pbed, pbim, pfam -> [dataset, subpop, bed, bim, fam, spop, pbed, pbim, pfam]}      
    impute_subpop(impute_input)
      
    // //Step 4.2.2: Dataset
    beagl_input = make_snp2hlarefpanel_dataset.out.map{dataset, subpop, datas -> [dataset, subpop, datas[9], datas[10], datas[12]]}
    imput_input = preparedata_imputation.out.combine(beagl_input)
    .map {dataset, subpop, bed, bim, fam, pop, spop, pbed, pbim, pfam -> [dataset, subpop, bed, bim, fam, spop, pbed, pbim, pfam]}      
    impute_dataset(imput_input)
    

    // Step 4.3: Proceed with imputation
    // Step 4.3.1: Subpop
    inp = make_snp2hlarefpanel_subpop.out.map{dataset, subpop, datas -> [dataset, subpop, datas[14]]}
    .combine(makeref_snp2hla_phasing_subpop.out, by:[0,1])
    .map{dataset, subpop, markers, phased -> [subpop, markers, phased]}
    imp_input = impute_subpop.out.map{datset, spop, pop, data -> [pop, spop, data[4]]}
    .combine(inp, by:0)
    subpop_imputation(imp_input)
    
    // Step 4.3.1: Dataset
    in = make_snp2hlarefpanel_dataset.out.map{dataset, subpop, datas -> [dataset, subpop, datas[14]]}
    .combine(makeref_snp2hla_phasing_dataset.out, by:[0,1])
    .map{dataset, subpop, markers, phased -> [subpop, markers, phased]}
    impt_input = impute_dataset.out.map{datset, spop, pop, data -> [pop, spop, data[4]]}
    .combine(in, by:0)
    dataset_imputation(impt_input)

    // Step 4.4
    // Convert posterior probability to dosage format
    // Step 4.4.1: Subpop
    in = impute_subpop.out
    .map {dtset, subp, ref, result -> [subp, ref, result[12]] }
    ref = make_snp2hlarefpanel_subpop.out.map{dataset, subpop, datas -> [dataset, subpop, datas[10]]}
    dos_inp = subpop_imputation.out
    .map{dataset, reference, results -> [dataset, reference, results[1], results[2], results[3]]}
    .combine(in, by:[0,1])
    .combine(ref, by: 1)
    .map{data, refenc, probs, phased, r2, dfam, set, refbim -> [data, refenc, probs, phased, r2, dfam, refbim]}
    posteprob_dosage_subpop(dos_inp)
    
    // Step 4.4.2: Dataset
    in = impute_dataset.out
    .map {dtset, subp, ref, result -> [subp, ref, result[12]] }
    ref = make_snp2hlarefpanel_dataset.out.map{dataset, subpop, datas -> [dataset, subpop, datas[10]]}
    dos_inp = dataset_imputation.out
    .map{dataset, reference, results -> [dataset, reference, results[1], results[2], results[3]]}
    .combine(in, by:[0,1])
    .combine(ref, by: 1)
    .map{data, refenc, probs, phased, r2, dfam, set, refbim -> [data, refenc, probs, phased, r2, dfam, refbim]}
    posteprob_dosage_dataset(dos_inp)
      
    // MASKED DATA
    // masked_data = Channel.fromList(params.masked_data)
    // .map{dataset, subpop, pbed, pbim, pfam  -> [dataset, subpop, file(pbed), file(pbim), file(pfam)]}
    // bgl_input = make_snp2hlarefpanel_subpop.out.map{dataset, subpop, datas -> [dataset, subpop, datas[9], datas[10], datas[12]]}
    // impute_input = masked_data.combine(bgl_input)
    // .map {dataset, subpop, bed, bim, fam, pop, spop, pbed, pbim, pfam -> [dataset, subpop, bed, bim, fam, spop, pbed, pbim, pfam]}      
    // // impute_input.view() 
    // impute_subpop(impute_input)
    
    // inp = make_snp2hlarefpanel_subpop.out.map{dataset, subpop, datas -> [dataset, subpop, datas[14]]}
    // .combine(makeref_snp2hla_phasing_subpop.out, by:[0,1])
    // .map{dataset, subpop, markers, phased -> [subpop, markers, phased]}
    // imp_input = impute_subpop.out.map{datset, spop, pop, data -> [pop, spop, data[4]]}
    // .combine(inp, by:0)
    // subpop_imputation(imp_input)


    // Train and predict using SOFTWARE 2
    //2. HIBAG
    // 2.1 Subpop
    subpop_modeldata = Channel.fromList(params.subpop_modeldata)
    .map{dataset, subpop, hlaA, hlaB, hlaC  -> [dataset, subpop, file(hlaA), file(hlaB), file(hlaC)]}
    input_ch = get_geno_plink_subpop.out.combine(subset_hlatypes_hibag.out, by:[0,1])
    .map{dataset, subpop, bed, bim, fam, snp2hla_ped, hibag_HLAType -> [dataset, subpop, bed, bim, fam, hibag_HLAType]} 
    .combine(subpop_modeldata, by:[0,1])
    test_data = preparedata_imputation.out.combine(input_ch)
    .map{dataset, subpop, pbed, pbim, pfam, pop, spop, bed, bim, fam, hibag_HLAType, hlaA, hlaB, hlaC -> [dataset, subpop, pbed, pbim, pfam, spop, bed, bim, fam, hibag_HLAType, hlaA, hlaB, hlaC]}
    // test_data.view()
    hibag_impute_subpop(test_data)
    

    // 2.2 Dataset
    dataset_modeldata = Channel.fromList(params.dataset_modeldata)
    .map{dataset, subpop, hlaA, hlaB, hlaC  -> [dataset, subpop, file(hlaA), file(hlaB), file(hlaC)]}
    input_two_ch = get_geno_plink_dataset.out.combine(combine_hlatypes_hibag.out, by:0)
    .map{dataset, subpop, bed, bim, fam, snp2hla_ped, hibag_HLAType -> [dataset, subpop, bed, bim, fam, hibag_HLAType]} 
    .combine(dataset_modeldata, by:[0,1])
    tst_data = preparedata_imputation.out.combine(input_two_ch)
    .map{dataset, subpop, pbed, pbim, pfam, pop, spop, bed, bim, fam, hibag_HLAType, hlaA, hlaB, hlaC -> [dataset, subpop, pbed, pbim, pfam, spop, bed, bim, fam, hibag_HLAType, hlaA, hlaB, hlaC]}
    // tst_data.view()
    hibag_impute_dataset(tst_data)


    //MASKED DATA
    // masked_data = Channel.fromList(params.masked_data)
    // .map{dataset, subpop, pbed, pbim, pfam  -> [dataset, subpop, file(pbed), file(pbim), file(pfam)]}
    // subpop_modeldata = Channel.fromList(params.subpop_modeldata)
    // .map{dataset, subpop, hlaA, hlaB, hlaC  -> [dataset, subpop, file(hlaA), file(hlaB), file(hlaC)]}
    // input_ch = get_geno_plink_subpop.out.combine(subset_hlatypes_hibag.out, by:[0,1])
    // .map{dataset, subpop, bed, bim, fam, snp2hla_ped, hibag_HLAType -> [dataset, subpop, bed, bim, fam, hibag_HLAType]} 
    // .combine(subpop_modeldata, by:[0,1])
    // test_data = masked_data.combine(input_ch)
    // .map{dataset, subpop, pbed, pbim, pfam, pop, spop, bed, bim, fam, hibag_HLAType, hlaA, hlaB, hlaC -> [dataset, subpop, pbed, pbim, pfam, spop, bed, bim, fam, hibag_HLAType, hlaA, hlaB, hlaC]}
    // test_data.view()
    // hibag_impute_subpop(test_data)
}

