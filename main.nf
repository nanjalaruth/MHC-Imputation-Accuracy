#!/usr/bin/env nextflow

// Enable DSL syntax
nextflow.enable.dsl=2

// Include modules
// Make reference
include {  combine_hlatypes_snp2hla as combine_hlatypes_H3A_snp2hla;
combine_hlatypes_snp2hla as combine_hlatypes_1kg_All_snp2hla; combine_hlatypes_hibag; 
combine_hlatypes_hibag as combine_hlatypes_1kg_All_hibag;
combine_hlatypes_hibag as combine_hlatypes_H3A_hibag;;
subset_hlatypes_snp2hla; subset_hlatypes_hibag; get_geno_plink as get_geno_plink_dataset; 
get_geno_plink as get_geno_plink_H3A_dataset; get_geno_plink as get_geno_plink_subpop; 
make_snp2hlarefpanel as make_snp2hlarefpanel_dataset; 
make_snp2hlarefpanel as make_snp2hlarefpanel_H3A_dataset
make_snp2hlarefpanel as make_snp2hlarefpanel_subpop; 
makeref_snp2hla_phasing as makeref_snp2hla_phasing_dataset; 
makeref_snp2hla_phasing as makeref_snp2hla_phasing_H3A_dataset;
makeref_snp2hla_phasing as makeref_snp2hla_phasing_subpop; 
convert_to_beagle as convert_to_beagle_dataset; 
convert_to_beagle as convert_to_beagle_H3A_dataset;
convert_to_beagle as convert_to_beagle_subpop} from './modules/make_reference'

// SNP2HLA Imputation
include { preparedata_imputation as prepare_hg18_data_imputation;
preparedata_imputation as prepare_hg19_data_imputation;
impute_data as impute_dataset; 
impute_data as impute_H3A_dataset; 
impute_data as impute_subpop; imputation as subpop_imputation; 
imputation as dataset_H3A_imputation;
imputation as dataset_imputation; posteprob_dosage as posteprob_dosage_subpop;
posteprob_dosage as posteprob_dosage_dataset;
posteprob_dosage as posteprob_dosage_H3A_dataset;
measureacc as measureacc_H3A_dataset;
measureacc as measureacc_dataset; measureacc as measureacc_subpop} from './modules/Impute.nf'

// HIBAG Imputation
include { hibag_impute as hibag_impute_dataset; hibag_impute as hibag_impute_subpop } from './modules/Hibag.nf'

// CookHLA Imputation
include { makegenetic_map as makegenetic_map_dataset; 
makegenetic_map as makegenetic_map_subpop; makegenetic_map as makegenetic_map_H3A_dataset; 
cookHLAimpute as cookHLAimpute_subpop; cookHLAimpute as cookHLAimpute_H3A_dataset;
cookHLAimpute as cookHLAimpute_dataset; measureaccuracy as measureaccuracy_subpop;
measureaccuracy as measureaccuracy_dataset;
measureaccuracy as measureaccuracy_H3A_dataset} from './modules/cookHLA.nf'

// Generate SNP2HLA Info file
include { annotate_vcf as combine_dataset_H3A_output; annotate_vcf as combine_dataset_output;
annotate_vcf as combine_subpop_output; generate_info as generate_subpop_info;
generate_info as generate_dataset_info; generate_info as generate_dataset_H3A_info } from './modules/snp2hala_info.nf'

// Main workflow
workflow{

    // Step 1. Getting HLA TYPES 
    // Step 1.1. For SNP2HLA
    // Step 1.1.1. 1kg-All
    samples_ch = Channel.fromList(params.hlatype_files)
                    .map{ name, hlatype -> 
                        hlatypes = file(hlatype)
                        return [name, hlatypes] 
                    }
    combine_hlatypes_1kg_All_snp2hla(samples_ch)

    // Step 1.1.2. H3Africa
    samples_ch = Channel.fromList(params.H3Africa_hlatypes)
                    .map{ name, hlatype -> 
                        hlatypes = file(hlatype)
                        return [name, hlatypes] 
                    }
    combine_hlatypes_H3A_snp2hla(samples_ch)

    // // Step 1.1.3. 1kg-All Subpopulation (Gambian, African)
    subpop_ids_ch = Channel.fromList(params.subpop_ids)
                .map{ dataset, subpop, ids  ->
                    ids_file = file(ids)
                    return [dataset, subpop, ids_file]
                } 
                .combine(combine_hlatypes_1kg_All_snp2hla.out, by:0)
    subset_hlatypes_snp2hla(subpop_ids_ch)


    // Step 1.2. For HIBAG
    // Step 1.2.1. 1kg-All
    smp_ch = Channel.fromList(params.hlatype_files)
                    .map{ title, hlatype -> 
                        hlatypes = file(hlatype)
                        return [title, hlatypes] 
                    }
    combine_hlatypes_1kg_All_hibag(smp_ch)

    // Step 1.2.2. H3Africa
    smp_ch = Channel.fromList(params.H3Africa_hlatypes)
                    .map{ title, hlatype -> 
                        hlatypes = file(hlatype)
                        return [title, hlatypes] 
                    }
    combine_hlatypes_H3A_hibag(smp_ch)

    // // Step 1.2.3. 1kg subpopulation (Gambian, African)
    subpop_id_ch = Channel.fromList(params.subpop_ids)
                .map{ dataset, subpop, ids -> 
                    id_file = file(ids)
                    return [dataset, subpop, id_file] 
                }
                .combine(combine_hlatypes_1kg_All_hibag.out, by:0)
    subset_hlatypes_hibag(subpop_id_ch)


    // Step 2. Getting SNP GENOTYPES
    // Step 2.1. 1kg-All
    dataset_geno_ch = Channel.fromList(params.genotype_files)
                .map{ dataset, genotype -> 
                    genotypes = file(genotype)
                    return [dataset, genotypes] 
                }
                .combine(combine_hlatypes_1kg_All_snp2hla.out, by:0)
                .map{ dataset, genotypes, hlatypes -> 
                    return [dataset, genotypes, dataset, hlatypes] 
                }
                .combine(combine_hlatypes_1kg_All_hibag.out, by:0)        
    get_geno_plink_dataset(dataset_geno_ch)


    // Step 2.2. H3Africa
    dataset_geno_ch = Channel.fromList(params.H3Africa_genotypes)
                .map{ dataset, genotype -> 
                    genotypes = file(genotype)
                    return [dataset, genotypes] 
                }
                .combine(combine_hlatypes_H3A_snp2hla.out, by:0)
                .map{ dataset, genotypes, hlatypes -> 
                    return [dataset, genotypes, dataset, hlatypes]
                }
                .combine(combine_hlatypes_H3A_hibag.out, by:0)     
    get_geno_plink_H3A_dataset(dataset_geno_ch)


    // Step 2.2. 1kg-All Subpopulation (Gambian, African)
    subpop_geno_ch = Channel.fromList(params.genotype_files)
                    .map{ dataset, genotype -> 
                        genotypes = file(genotype)
                        return [dataset, genotypes] 
                    }.combine(subset_hlatypes_snp2hla.out, by:0)
                    .map{ dataset, genotype, subpop, snp2hlatype -> [dataset, subpop, genotype, snp2hlatype]}
                    .combine(subset_hlatypes_hibag.out, by:[0,1])
                    .map{ dataset, subpop, genotype, snp2hlatype, hibaghlatypes
                    -> [dataset, genotype, subpop, snp2hlatype, hibaghlatypes]}
    get_geno_plink_subpop(subpop_geno_ch)


    // Step 3: Make reference panel
    // Step 3.1. For SNP2HLA
    // Step 3.1.1 Get first output files
    // Step 3.1.1.1. 1kg-All
    input = get_geno_plink_dataset.out
    .map{dataset, subpop, bed, bim, fam, snphla, hibaghla -> [dataset, subpop, bed, bim, fam, snphla]}
    make_snp2hlarefpanel_dataset(input)

    //  Step 3.1.1.2. H3Africa
    input = get_geno_plink_H3A_dataset.out
    .map{dataset, subpop, bed, bim, fam, snphla, hibaghla -> [dataset, subpop, bed, bim, fam, snphla]}
    make_snp2hlarefpanel_H3A_dataset(input)
    
    //  Step 3.1.1.3. 1kg-All subpop (African&Gambian)
    inp = get_geno_plink_subpop.out
    .map{dataset, subpop, bed, bim, fam, snphla, hibaghla -> [dataset, subpop, bed, bim, fam, snphla]}
    make_snp2hlarefpanel_subpop(inp)


    // Step 3.1.2 Convert to Beagle
    // Step 3.1.2.1: 1kg-All
    beagle_input = make_snp2hlarefpanel_dataset.out.map{dataset, subpop, datas -> [dataset, subpop, datas[11], datas[15]]}
    convert_to_beagle_dataset(beagle_input)

    //Step 3.1.2.2. H3Africa
    beagle_input = make_snp2hlarefpanel_H3A_dataset.out.map{dataset, subpop, datas -> [dataset, subpop, datas[11], datas[15]]}
    convert_to_beagle_H3A_dataset(beagle_input)

    // Step 3.1.2.3. 1kg-All subpop (African&Gambian)
    beagle_subpop_input = make_snp2hlarefpanel_subpop.out.map{dataset, subpop, datas -> [dataset, subpop, datas[11], datas[15]]}
    convert_to_beagle_subpop(beagle_subpop_input)


    // Step 3.1.3: Phasing
    // Step 3.1.3.1: 1kg-All
    makeref_snp2hla_phasing_dataset(convert_to_beagle_dataset.out)

    //Step 3.1.3.2. H3Africa
    makeref_snp2hla_phasing_H3A_dataset(convert_to_beagle_H3A_dataset.out)

    //  Step 3.1.3.3. 1kg-All subpop (African&Gambian)
    makeref_snp2hla_phasing_subpop(convert_to_beagle_subpop.out)


    // Step 4: Imputation
    // Step 4.1: prepare data (hg18 format)
    dataset_testdata_ch = Channel.fromList(params.sample_data_hg18)
                    .map{ pop, chip, sample_geno -> 
                        genotype = file(sample_geno)
                        return [pop, chip, genotype]
                    }  
    prepare_hg18_data_imputation(dataset_testdata_ch)

    // HIBAG test data (hg19 format)
    data_ch = Channel.fromList(params.sample_data_hg19)
                    .map{ pop, chip, sample_geno -> 
                        genotype = file(sample_geno)
                        return [pop, chip, genotype]
                    }  
    prepare_hg19_data_imputation(data_ch)

    // Step 4.2: Pre-imputation
    // Step 4.2.1 1kg-All subpop (African&Gambian)
    bgl_input = make_snp2hlarefpanel_subpop.out.map{dataset, subpop, datas -> [dataset, subpop, datas[9], datas[10], datas[12]]}
    impute_input = prepare_hg18_data_imputation.out.combine(bgl_input)
    .map {dataset, subpop, bed, bim, fam, pop, spop, pbed, pbim, pfam -> [dataset, subpop, bed, bim, fam, spop, pbed, pbim, pfam]}      
    impute_subpop(impute_input)
      
    // Step 4.2.2: 1kg-All
    beagl_input = make_snp2hlarefpanel_dataset.out.map{dataset, subpop, datas -> [dataset, subpop, datas[9], datas[10], datas[12]]}
    imput_input = prepare_hg18_data_imputation.out.combine(beagl_input)
    .map {dataset, subpop, bed, bim, fam, pop, spop, pbed, pbim, pfam -> [dataset, subpop, bed, bim, fam, spop, pbed, pbim, pfam]}      
    impute_dataset(imput_input)
    
    //  Step 4.2.3 H3Africa
    beagl_input = make_snp2hlarefpanel_H3A_dataset.out.map{dataset, subpop, datas -> [dataset, subpop, datas[9], datas[10], datas[12]]}
    imput_input = prepare_hg18_data_imputation.out.combine(beagl_input)
    .map {dataset, subpop, bed, bim, fam, pop, spop, pbed, pbim, pfam -> [dataset, subpop, bed, bim, fam, spop, pbed, pbim, pfam]}      
    impute_H3A_dataset(imput_input)


    // Step 4.3: Proceed with imputation
    // Step 4.3.1: 1kg-All subpop (African&Gambian)
    inp = make_snp2hlarefpanel_subpop.out.map{dataset, subpop, datas -> [dataset, subpop, datas[14]]}
    .combine(makeref_snp2hla_phasing_subpop.out, by:[0,1])
    .map{dataset, subpop, markers, phased -> [subpop, markers, phased]}
    imp_input = impute_subpop.out.map{datset, spop, pop, data -> [pop, spop, data[4]]}
    .combine(inp, by:0)
    subpop_imputation(imp_input)
    
    // Step 4.3.1: 1kg-All
    in = make_snp2hlarefpanel_dataset.out.map{dataset, subpop, datas -> [dataset, subpop, datas[14]]}
    .combine(makeref_snp2hla_phasing_dataset.out, by:[0,1])
    .map{dataset, subpop, markers, phased -> [subpop, markers, phased]}
    impt_input = impute_dataset.out.map{datset, spop, pop, data -> [pop, spop, data[4]]}
    .combine(in, by:0)
    dataset_imputation(impt_input)

    // Step 4.3.2: H3Africa
    in = make_snp2hlarefpanel_H3A_dataset.out.map{dataset, subpop, datas -> [dataset, subpop, datas[14]]}
    .combine(makeref_snp2hla_phasing_H3A_dataset.out, by:[0,1])
    .map{dataset, subpop, markers, phased -> [subpop, markers, phased]}
    impt_input = impute_H3A_dataset.out.map{datset, spop, pop, data -> [pop, spop, data[4]]}
    .combine(in, by:0)
    dataset_H3A_imputation(impt_input)


    // Step 4.4
    // Convert posterior probability to dosage format
    // Step 4.4.1: 1kg-All subpop (African&Gambian)
    in = impute_subpop.out
    .map {dtset, subp, ref, result -> [subp, ref, result[12]] }
    ref = make_snp2hlarefpanel_subpop.out.map{dataset, subpop, datas -> [dataset, subpop, datas[10]]}
    dos_inp = subpop_imputation.out
    .map{dataset, reference, results -> [dataset, reference, results[1], results[2], results[3]]}
    .combine(in, by:[0,1])
    .combine(ref, by: 1)
    .map{data, refenc, probs, phased, r2, dfam, set, refbim -> [data, refenc, probs, phased, r2, dfam, refbim]}
    posteprob_dosage_subpop(dos_inp)
    
    // Step 4.4.2: 1kg-All
    in = impute_dataset.out
    .map {dtset, subp, ref, result -> [subp, ref, result[12]] }
    ref = make_snp2hlarefpanel_dataset.out.map{dataset, subpop, datas -> [dataset, subpop, datas[10]]}
    dos_inp = dataset_imputation.out
    .map{dataset, reference, results -> [dataset, reference, results[1], results[2], results[3]]}
    .combine(in, by:[0,1])
    .combine(ref, by: 1)
    .map{data, refenc, probs, phased, r2, dfam, set, refbim -> [data, refenc, probs, phased, r2, dfam, refbim]}
    posteprob_dosage_dataset(dos_inp)

    // Step 4.4.3: H3Africa
    in = impute_H3A_dataset.out
    .map {dtset, subp, ref, result -> [subp, ref, result[12]] }
    ref = make_snp2hlarefpanel_H3A_dataset.out.map{dataset, subpop, datas -> [dataset, subpop, datas[10]]}
    dos_inp = dataset_H3A_imputation.out
    .map{dataset, reference, results -> [dataset, reference, results[0], results[1], results[2]]}
    .combine(in, by:[0,1])
    .combine(ref, by: 1)
    .map{data, refenc, probs, phased, r2, dfam, set, refbim -> [data, refenc, probs, phased, r2, dfam, refbim]}
    // dos_inp.view()
    posteprob_dosage_H3A_dataset(dos_inp)


    // Step 4.5: MeasureAccuracy
    // Step 4.5.1 1kg-All subpop (African&Gambian)
    ans =  Channel.fromPath(params.answer_file)
    .combine(posteprob_dosage_subpop.out)
    .map{answer, array, ref, rest -> [answer, array, ref, rest[2]]}
    measureacc_subpop(ans)

    // Step 4.5.2 1kg-All
    ans =  Channel.fromPath(params.answer_file)
    .combine(posteprob_dosage_dataset.out)
    .map{answer, array, ref, rest -> [answer, array, ref, rest[2]]}
    measureacc_dataset(ans)

    // Step 4.5.3 H3A
    ans =  Channel.fromPath(params.answer_file)
    .combine(posteprob_dosage_H3A_dataset.out)
    .map{answer, array, ref, rest -> [answer, array, ref, rest[2]]}
    measureacc_H3A_dataset(ans)


    // TRAIN AND PREDICT USING SOFTWARE 2
    //2. HIBAG
    // 2.1 1kg-All subpop (African&Gambian)
    answerfile = Channel.fromPath(params.hibag_answerfile)
    subpop_modeldata = Channel.fromList(params.subpop_modeldata)
    .map{dataset, subpop, hlaA, hlaB, hlaC  -> [dataset, subpop, file(hlaA), file(hlaB), file(hlaC)]}
    test_data = prepare_hg19_data_imputation.out.combine(subpop_modeldata)
    .map{dataset, subpop, bed, bim, fam, dtset, spop, model_a, model_b, model_c -> [dataset, subpop, bed, bim, fam, spop, model_a, model_b, model_c]}
    .combine(answerfile)
    hibag_impute_subpop(test_data)
    

    // 2.2 1kg-All & H3Africa
    answerfile = Channel.fromPath(params.hibag_answerfile)
    dataset_modeldata = Channel.fromList(params.dataset_modeldata)
    .map{dataset, subpop, hlaA, hlaB, hlaC  -> [dataset, subpop, file(hlaA), file(hlaB), file(hlaC)]}
    tst_data = prepare_hg19_data_imputation.out.combine(dataset_modeldata)
    .map{dataset, subpop, bed, bim, fam, dtset, spop, model_a, model_b, model_c -> [dataset, subpop, bed, bim, fam, spop, model_a, model_b, model_c]}
    .combine(answerfile)
    hibag_impute_dataset(tst_data)

    // GENERATE SNP2HLA INFO FILE FOR COMPARISON WITH GENERAL IMPUTATION TOOLS
    // subpop
    file = posteprob_dosage_subpop.out
    .map {target, reference, out -> [target, reference, out[0], out[4], out[6], out[8], out[10], out[3] ]}
    combine_subpop_output(file)

    input = combine_subpop_output.out.map{data, ref, r2, id, vcf -> [data, ref, r2, id]}
    generate_subpop_info(input)
    
    // Dataset
    // 1kg_All
    file = posteprob_dosage_dataset.out
    .map {target, reference, out -> [target, reference, out[0], out[4], out[6], out[8], out[10], out[3]]}
    combine_dataset_output(file)

    input = combine_dataset_output.out.map{data, ref, r2, id, vcf -> [data, ref, r2, id]}
    generate_dataset_info(input)

    // H3A
    file = posteprob_dosage_H3A_dataset.out
    .map {target, reference, out -> [target, reference, out[0], out[4], out[6], out[8], out[10], out[3]]}
    combine_dataset_H3A_output(file)
    
    input = combine_dataset_H3A_output.out.map{data, ref, r2, id, vcf -> [data, ref, r2, id]}
    generate_dataset_H3A_info(input)

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


    // TRAIN AND PREDICT USING SOFTWARE 3: CookHLA
    // Step 1:Make GeneticMap
    // Step 1.1: Subpop
    ref = make_snp2hlarefpanel_subpop.out
    .map{dataset, subpop, inp -> [dataset, subpop, inp[2], inp[9], inp[10], inp[12], inp[14]]}
    .combine(makeref_snp2hla_phasing_subpop.out, by:[0,1])
    data_input = prepare_hg18_data_imputation.out.combine(ref)
    .map{dataset, subpop, bed, bim, fam, dtset, spop, frq, rbed, rbim, rfam, markers, bglphased 
    -> [dataset, subpop, bed, bim, fam, spop, frq, rbed, rbim, rfam, markers, bglphased]}
    // data_input.view()
    makegenetic_map_subpop(data_input)

    // // Step 1.2: Dataset
    ref = make_snp2hlarefpanel_dataset.out
    .map{dataset, subpop, inp -> [dataset, subpop, inp[2], inp[9], inp[10], inp[12], inp[14]]}
    .combine(makeref_snp2hla_phasing_dataset.out, by:[0,1])
    data_input = prepare_hg18_data_imputation.out.combine(ref)
    .map{dataset, subpop, bed, bim, fam, dtset, spop, frq, rbed, rbim, rfam, markers, bglphased 
    -> [dataset, subpop, bed, bim, fam, spop, frq, rbed, rbim, rfam, markers, bglphased]}
    // data_input.view()
    makegenetic_map_dataset(data_input)

    //  Step 1.3: H3Africa
    ref = make_snp2hlarefpanel_H3A_dataset.out
    .map{dataset, subpop, inp -> [dataset, subpop, inp[2], inp[9], inp[10], inp[12], inp[14]]}
    .combine(makeref_snp2hla_phasing_H3A_dataset.out, by:[0,1])
    data_input = prepare_hg18_data_imputation.out.combine(ref)
    .map{dataset, subpop, bed, bim, fam, dtset, spop, frq, rbed, rbim, rfam, markers, bglphased 
    -> [dataset, subpop, bed, bim, fam, spop, frq, rbed, rbim, rfam, markers, bglphased]}
    // data_input.view()
    makegenetic_map_H3A_dataset(data_input)

    // Step 2: Imputation
    // Step 2.1: Dataset
    ref = make_snp2hlarefpanel_dataset.out
    .map{dataset, subpop, inp -> [dataset, subpop, inp[2], inp[9], inp[10], inp[12], inp[14]]}
    .combine(makeref_snp2hla_phasing_dataset.out, by:[0,1])
    data_input = prepare_hg18_data_imputation.out.combine(ref)
    .map{dataset, subpop, bed, bim, fam, dtset, spop, frq, rbed, rbim, rfam, markers, bglphased 
    -> [subpop, spop, bed, bim, fam, frq, rbed, rbim, rfam, markers, bglphased]}
    .combine(makegenetic_map_dataset.out, by:[0,1])
    .map{subpop, spop, bed, bim, fam, frq, rbed, rbim, rfam, markers, bglphased, erate, mach
    -> [subpop, bed, bim, fam, spop, frq, rbed, rbim, rfam, markers, bglphased, erate, mach]}
    // data_input.view()
    // cookHLAimpute_dataset(data_input)

    // H3A dataset
    ref = make_snp2hlarefpanel_H3A_dataset.out
    .map{dataset, subpop, inp -> [dataset, subpop, inp[2], inp[9], inp[10], inp[12], inp[14]]}
    .combine(makeref_snp2hla_phasing_H3A_dataset.out, by:[0,1])
    data_input = prepare_hg18_data_imputation.out.combine(ref)
    .map{dataset, subpop, bed, bim, fam, dtset, spop, frq, rbed, rbim, rfam, markers, bglphased 
    -> [subpop, spop, bed, bim, fam, frq, rbed, rbim, rfam, markers, bglphased]}
    .combine(makegenetic_map_H3A_dataset.out, by:[0,1])
    .map{subpop, spop, bed, bim, fam, frq, rbed, rbim, rfam, markers, bglphased, erate, mach
    -> [subpop, bed, bim, fam, spop, frq, rbed, rbim, rfam, markers, bglphased, erate, mach]}
    // data_input.view()
    cookHLAimpute_H3A_dataset(data_input)


    // Step 2.2: Subpop
    ref = make_snp2hlarefpanel_subpop.out
    .map{dataset, subpop, inp -> [dataset, subpop, inp[2], inp[9], inp[10], inp[12], inp[14]]}
    .combine(makeref_snp2hla_phasing_subpop.out, by:[0,1])
    data_input = prepare_hg18_data_imputation.out.combine(ref)
    .map{dataset, subpop, bed, bim, fam, dtset, spop, frq, rbed, rbim, rfam, markers, bglphased 
    -> [subpop, spop, bed, bim, fam, frq, rbed, rbim, rfam, markers, bglphased]}
    .combine(makegenetic_map_subpop.out, by:[0,1])
    .map{subpop, spop, bed, bim, fam, frq, rbed, rbim, rfam, markers, bglphased, erate, mach
    -> [subpop, bed, bim, fam, spop, frq, rbed, rbim, rfam, markers, bglphased, erate, mach]}
    // data_input.view()
    cookHLAimpute_subpop(data_input)


    // Step 3: Measure accuracy
    // Step 3.1: Subpop
    ans =  Channel.fromPath(params.cookanswer_file)
    .combine(cookHLAimpute_subpop.out)
    .map{true_allele, array, ref, imputed_allele -> [true_allele, array, ref, imputed_allele[4]]}
    // ans.view()
    measureaccuracy_subpop(ans)

    // Step 3.2: Dataset
    // ans =  Channel.fromPath(params.answer_file)
    // .combine(cookHLAimpute_dataset.out)
    // .map{true_allele, array, ref, imputed_allele -> [true_allele, array, ref, imputed_allele[4]]}

    // measureaccuracy_dataset(ans)

    // H3A dataset
     ans =  Channel.fromPath(params.cookanswer_file)
    .combine(cookHLAimpute_H3A_dataset.out)
    .map{true_allele, array, ref, imputed_allele -> [true_allele, array, ref, imputed_allele[4]]}
    // ans.view()
    measureaccuracy_H3A_dataset(ans)

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
}

