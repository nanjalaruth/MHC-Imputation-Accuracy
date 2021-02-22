#!/usr/bin/env nextflow

params.dir = "./datasets"
dir = params.dir
params.pops = ["HAPMAP_CEU"]

Channel.from(params.pops)
    .map { pop ->
        [ file("$dir/${pop}.bed"), file("$dir/${pop}.bim"), file("$dir/${pop}.fam")]
    }
    .set { plink_data }
    
plink_data.subscribe { println "$it" }
