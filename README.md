# HLA Imputation Accuracy Workflow 

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker%20registry-Quay.io-red)](https://quay.io/repository/nanjalaruth/impute-hla?tab=tags)
![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)

## Introduction
Association studies for instance GWAS, traditionally use genotyping arrays to genotype large set of individuals and in this way determine SNPs that are significantly overrepresented in the cases compared to the controls and thus determine association with disease. Genotyping arrays are cheaper than sequencing but can only measure a tag of SNPs from the >300 million SNPs that are available. In order to increase the number of SNPs that can be used for association studies, genotype imputation is performed. Genotype imputation refers to the statistical inference of unobserved genotypes. 

Some regions within the human genome (MHC otherwise known as HLA) are highly variable and maybe difficult to impute. The HLA region has been associated to *autoimmune* diseases such as rheumatoid arthritis and *infectious diseases* such as HIV/AIDS. Accurate imputation of this region is key, as it would help increase the chances of identifying the causal variants of some autoimmune and immune mediated diseases. 

Genotype imputation is a statistical process and thus needs to be assessed to ensure that the predicted genotypes are accurate.

The project focused on assessing the accuracy of imputing HLA Class I alleles in __selected African populations.__ 

Imputation accuracy was based on [SNP2HLA](http://software.broadinstitute.org/mpg/snp2hla/) and [HIBAG](https://github.com/zhengxwen/HIBAG) __imputation tools__, 1kg-All, 1kg-Gwd, 1kg-Afr, H3Africa, prebuilt EUR __reference panels__ and Illumina Omni 2.5 array, H3Africa array __genotyping arrays__  

## Installation 
1. [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
2. [Docker](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04) 
3. [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html)

**N/B** You do not need to install any other tool as `singularity` profile will download the singularity image from https://quay.io/nanjalaruth/impute-hla

## Preparing Input files
### Target Genotype file

The input file must be a VCF file. As the work focuses on the HLA region, you are required to only use SNPs in chr6:29-34Mb. Thus, you can portably prepare only those SNPs in that region as input file.
SNP2HLA uses hg18 data while HIBAG uses hg19 data. You are therefore required to provide both files as input.

### Reference panel

The workflow focuses on running the data on a pre built European reference panel and reference panels that are built in the process of running the pipeline.
To custom make a reference panel, a HLA type file and SNP genotype file are required. The pipeline thus requires a HLA type file and SNP genotype file to make the reference panel.

I focused on 4  custom made reference panels; 1kg-All, 1kg-Gwd, 1kg-Afr, H3Africa. The 1kg-Gwd and 1kg-Afr are a subset of 1kg-All. In case you have a similar dataset, assign the path to the file with the subpopulation rsids to the flag termed `subpop_ids`. This is to bypass the need to get a vcf file that matches the subsetted population.

If you are working with populations that are not linked, provide the paths to the SNP genotype vcf file as demonstrated in the `genotype_files` flag in the test.config file. If you have multiple files, you could assign them to the same flag by creating more lists like what has been done with the sample dataset. Also provide the path to the HLA type file as shown in the `hlatype_files` flag within the test.config file.

**N/B** HLA typing was done using the [Optitype tool](https://github.com/nf-core/hlatyping)

## Running the pipeline
There are 2 ways to run the pipeline:

### Download it from GitHub
#### Steps to follow
- Download the pipeline from GitHub
```
git clone git@github.com:nanjalaruth/MHC-Imputation-Accuracy.git
```
- Move to that folder
```
cd MHC-Imputation-Accuracy
```
- Edit the __*conf/test.config*__ to suit the path to where your datasets are stored.
- Run the command below
```
nextflow run main.nf -c conf/test.config -profile singularity
```

**N/B**
If you are using any job scheduler, be sure to include it in the profile. For example:
```
nextflow run main.nf -c conf/test.config -profile singularity, slurm
```

### Run the pipeline directly from GitHub
`NextFlow` will automatically fetch the pipeline from `GitHub` so you don't need to install it.

#### Steps to follow
- Create your own config file locally to suit the path to your datasets. Use *conf/test.config* as a template.
- Run the command below
```
nextflow run nanjalaruth/MHC-Imputation-Accuracy -c `path to your config file` -profile singularity
```

## In case of any querries, tag me on an issue or send me a message on twitter (@Ruthnanje) or [LinkedIn](https://www.linkedin.com/in/ruth-nanjala-17991117a/)
