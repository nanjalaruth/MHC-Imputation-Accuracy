# HLA Imputation Accuracy Workflow 

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker%20registry-Quay.io-red)](https://quay.io/repository/nanjalaruth/impute-hla?tab=tags)

## Greetings!!
### So excited to have you here!!

## Project description
Genotype imputation refers to the statistical inference of unobserved genotypes.

This project focused on assessing the accuracy of __imputation tools__ ([SNP2HLA](http://software.broadinstitute.org/mpg/snp2hla/) and [HIBAG](https://github.com/zhengxwen/HIBAG)), __reference panels__(1kg-All, 1kg-Gwd, 1kg-Afr, H3Africa) and __genotyping arrays__ (Illumina Omni 2.5 array, H3Africa array) used for imputation of the human leukocyte antigen (HLA) class I alleles in __selected African populations.__ 

The HLA region has been associated to *autoimmune* diseases such as rheumatoid arthritis and *infectious diseases* such as HIV/AIDS. Accurate imputation of this highly polymorphic region would increase the chances of identifying the causal variants of some autoimmune and immune mediated diseases.

## Installation 
1. [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
2. [Docker](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04) 
3. [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html)

## Running the pipeline
The pipeline does not require installation as `NextFlow` will automatically fetch it from `GitHub`.
Please edit the __*nextflow.config*__ to suit the path to where your datasets are stored.

To execute the pipeline run:
```
nextflow run nanjalaruth/MHC-Imputation-Accuracy/main.nf -profile singularity
```
- `singularity` profile will download the singularity image from https://quay.io/nanjalaruth/impute-hla

### Find me on twitter (@Ruthnanje) and [LinkedIn](https://www.linkedin.com/in/ruth-nanjala-17991117a/)
