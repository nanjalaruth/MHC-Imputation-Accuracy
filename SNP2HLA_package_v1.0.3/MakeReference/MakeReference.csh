#!/opt/conda/bin/tcsh

#################################################################################################
#   This script helps prepare a reference dataset for HLA imputation
#
#   Author: Sherman Jia (xiaomingjia@gmail.com)
#   Usage: ./MakeReference.csh SNPS[.bed/.bim/.fam] HLA.ped OUTPUT
#   Where the ped file contains the 4-digit HLA alleles (FID,IID,pID,mID,SEX,PHENO,A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1)
#   and bed/bim/fam contains the SNP data
#
#   HLA PED file should contain HLA alleles in the following (alphabetical) order:
#   HLA-A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1
#
#################################################################################################

# ./MakeReference.csh T1DGC T1DGC.HLA_TYPES.ped T1DGC.REF

if ($#argv != 4) then
    echo "USAGE: ./MakeReference.csh SNPS (.bed/.bim/.fam) HLAfile (FID,IID,pID,mID,SEX,PHENO,A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1) OUTPUT (.bgl.phased/.markers/.bed) plink"; exit 1
endif

set SCRIPTPATH=`dirname $0`

# You may download PLINK from (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml)

# CHECK FOR DEPENDENCIES
if (! -e `which $4`) then
    echo "Please install PLINK (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml) and point to the plink run file."; 
    echo "tcsh: use plink"
    echo "bash: use ./plink"
    exit 1
else if (! -e $SCRIPTPATH/beagle.jar) then
    echo "Please install Beagle 3 (http://faculty.washington.edu/browning/beagle/beagle.html#download) and copy the run file (beagle.3.0.4/beagle.jar) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/linkage2beagle.jar) then
    echo "Please copy linkage2beagle.jar (http://faculty.washington.edu/browning/beagle_utilities/utilities.html) (beagle.3.0.4/utility/linkage2beagle.jar) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/HLAtoSequences.pl) then
    echo "Please copy HLAtoSequences.pl (included with this package) into this directory."; exit 1
else if (! -e $SCRIPTPATH/encodeVariants.pl) then
    echo "Please copy encodeVariants.pl (included with this package) into this directory."; exit 1
else if (! -e $SCRIPTPATH/encodeHLA.pl) then
    echo "Please copy encodeHLA.pl (included with this package) into this directory."; exit 1
else if (! -e $SCRIPTPATH/HLA_DICTIONARY_AA.map) then
    echo "Please copy HLA_DICTIONARY_AA.map (included with this package) into this directory."; exit 1
else if (! -e $SCRIPTPATH/HLA_DICTIONARY_AA.txt) then
    echo "Please copy HLA_DICTIONARY_AA.txt (included with this package) into this directory."; exit 1
else if (! -e $SCRIPTPATH/HLA_DICTIONARY_SNPS.map) then
    echo "Please copy HLA_DICTIONARY_SNPS.map (included with this package) into this directory."; exit 1
else if (! -e $SCRIPTPATH/HLA_DICTIONARY_SNPS.txt) then
    echo "Please copy HLA_DICTIONARY_SNPS.txt (included with this package) into this directory."; exit 1
endif

set SNP_DATA=$1
set HLA_DATA=$2
set OUTPUT=$3
alias plink '$4 --noweb --silent'
alias beagle 'java -Xmx2000m -jar $SCRIPTPATH/beagle.jar'
alias linkage2beagle 'java -Xmx500m -jar $SCRIPTPATH/linkage2beagle.jar'

set ENCODE_AA          = 1
set ENCODE_HLA         = 1
set ENCODE_SNPS        = 1
set EXTRACT_FOUNDERS   = 1
set MERGE              = 1
set QC                 = 1
set PREPARE            = 1
set PHASE              = 1
set CLEANUP            = 1

set i=1

echo "Creating reference panel: $OUTPUT"


# Encode HLA amino acids
if ($ENCODE_AA) then
    echo "[$i] Generating amino acid sequences from HLA types.";  @ i++
    $SCRIPTPATH/HLAtoSequences.pl $HLA_DATA $SCRIPTPATH/HLA_DICTIONARY_AA.txt AA > $OUTPUT.AA.ped
    cp $SCRIPTPATH/HLA_DICTIONARY_AA.map $OUTPUT.AA.map

    echo "[$i] Encoding amino acids positions." ;  @ i++
    $SCRIPTPATH/encodeVariants.pl $OUTPUT.AA.ped $OUTPUT.AA.map $OUTPUT.AA.CODED

    plink --file $OUTPUT.AA.CODED --missing-genotype 0 --make-bed --out $OUTPUT.AA.TMP
    awk '{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}' $OUTPUT.AA.TMP.bim | grep -v INS | cut -f2 > to_remove
    plink --bfile $OUTPUT.AA.TMP --exclude to_remove --make-bed --out $OUTPUT.AA.CODED
    rm $OUTPUT.AA.TMP*; rm to_remove
    rm $OUTPUT.AA.???
endif

# Encode classical HLA alleles into binary format
if ($ENCODE_HLA) then
    echo "[$i] Encoding HLA alleles.";  @ i++
    $SCRIPTPATH/encodeHLA.pl $HLA_DATA $OUTPUT.HLA.map > $OUTPUT.HLA.ped
    plink --file $OUTPUT.HLA --make-bed --out $OUTPUT.HLA
endif

# Encode HLA SNPs
if ($ENCODE_SNPS) then
    echo "[$i] Generating DNA sequences from HLA types.";  @ i++
    $SCRIPTPATH/HLAtoSequences.pl $HLA_DATA $SCRIPTPATH/HLA_DICTIONARY_SNPS.txt SNPS > $OUTPUT.SNPS.ped
    cp $SCRIPTPATH/HLA_DICTIONARY_SNPS.map $OUTPUT.SNPS.map

    echo "[$i] Encoding SNP positions." ;  @ i++
    $SCRIPTPATH/encodeVariants.pl $OUTPUT.SNPS.ped $OUTPUT.SNPS.map $OUTPUT.SNPS.CODED
    plink --file $OUTPUT.SNPS.CODED --missing-genotype 0 --make-bed --out $OUTPUT.SNPS.TMP

    awk '{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}' $OUTPUT.SNPS.TMP.bim | grep -v INS | cut -f2 > to_remove
    plink --bfile $OUTPUT.SNPS.TMP --exclude to_remove --make-bed --out $OUTPUT.SNPS.CODED
    rm $OUTPUT.SNPS.TMP*; rm to_remove
    rm $OUTPUT.SNPS.???
endif

if ($EXTRACT_FOUNDERS) then
    echo "[$i] Extracting founders."; @ i++
    plink --bfile $SNP_DATA --filter-founders --mind 0.3 --alleleACGT --make-bed --out $SNP_DATA.FOUNDERS

    # Initial QC on Reference SNP panel
    plink --bfile $SNP_DATA.FOUNDERS --hardy        --out $SNP_DATA.FOUNDERS.hardy
    plink --bfile $SNP_DATA.FOUNDERS --freq         --out $SNP_DATA.FOUNDERS.freq
    plink --bfile $SNP_DATA.FOUNDERS --missing      --out $SNP_DATA.FOUNDERS.missing
    awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.hardy.hwe      | awk ' $9 < 0.000001 { print $2 }' | sort -u > remove.snps.hardy 
    awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.freq.frq       | awk ' $5 < 0.01 { print $2 } '             > remove.snps.freq
    awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.missing.lmiss  | awk ' $5 > 0.05 { print $2 } '              > remove.snps.missing
    cat remove.snps.*                                            | sort -u                                     > all.remove.snps

    plink --bfile $SNP_DATA.FOUNDERS --allow-no-sex --exclude all.remove.snps --make-bed --out $SNP_DATA.FOUNDERS.QC

    # Founders are identified here as individuals with "0"s in mother and father IDs in .fam file

    plink --bfile $OUTPUT.HLA --filter-founders --maf 0.0001 --make-bed --out $OUTPUT.HLA.FOUNDERS
    plink --bfile $OUTPUT.SNPS.CODED --filter-founders --maf 0.0001 --make-bed --out $OUTPUT.SNPS.FOUNDERS
    plink --bfile $OUTPUT.AA.CODED --filter-founders --maf 0.0001 --make-bed --out $OUTPUT.AA.FOUNDERS

    rm remove.snps.*
endif

# Merging SNP, HLA, and amino acid datasets 
if ($MERGE) then
    echo "[$i] Merging SNP, HLA, and amino acid datasets.";  @ i++
    echo "$OUTPUT.HLA.FOUNDERS.bed $OUTPUT.HLA.FOUNDERS.bim $OUTPUT.HLA.FOUNDERS.fam" > merge_list
    echo "$OUTPUT.AA.FOUNDERS.bed $OUTPUT.AA.FOUNDERS.bim $OUTPUT.AA.FOUNDERS.fam" >> merge_list
    echo "$OUTPUT.SNPS.FOUNDERS.bed $OUTPUT.SNPS.FOUNDERS.bim $OUTPUT.SNPS.FOUNDERS.fam" >> merge_list
    plink --bfile $SNP_DATA.FOUNDERS.QC --merge-list merge_list --make-bed --out $OUTPUT.MERGED.FOUNDERS
    rm $OUTPUT.HLA.???
    rm $OUTPUT.AA.CODED.???
    rm $OUTPUT.SNPS.CODED.???
    rm merge_list
endif

if ($QC) then
    echo "[$i] Performing quality control.";  @ i++
    plink --bfile $OUTPUT.MERGED.FOUNDERS --freq --out $OUTPUT.MERGED.FOUNDERS.FRQ
    awk '{if (NR > 1 && ($5 < 0.0001 || $5 > 0.9999)){print $2}}' $OUTPUT.MERGED.FOUNDERS.FRQ.frq > all.remove.snps
    awk '{if (NR > 1){if (($3 == "A" && $4 == "P") || ($4 == "A" && $3 == "P")){print $2 "\tP"}}}' $OUTPUT.MERGED.FOUNDERS.FRQ.frq > allele.order

    # QC: Maximum per-SNP missing > 0.5, MAF > 0.1%
    plink --bfile $OUTPUT.MERGED.FOUNDERS --reference-allele allele.order --exclude all.remove.snps --geno 0.5 --make-bed --out $OUTPUT

    # Calculate allele frequencies
    plink --bfile $OUTPUT --keep-allele-order --freq --out $OUTPUT.FRQ
    rm $SNP_DATA.FOUNDERS.*
    rm $OUTPUT.MERGED.FOUNDERS.*
    rm $OUTPUT.*.FOUNDERS.???
    rm allele.order
    rm all.remove.snps
endif

if ($PREPARE) then
    # Prepare files for Beagle phasing
    echo "[$i] Preparing files for Beagle. ";  @ i++
    awk '{print $2 " " $4 " " $5 " " $6}' $OUTPUT.bim > $OUTPUT.markers
    plink --bfile $OUTPUT --keep-allele-order --recode --alleleACGT --out $OUTPUT
    awk '{print "M " $2}' $OUTPUT.map > $OUTPUT.dat
    cut -d ' ' -f1-5,7- $OUTPUT.ped > $OUTPUT.nopheno.ped

    echo "[$i] Converting to beagle format.";  @ i++
    java -jar $SCRIPTPATH/linkage2beagle.jar $OUTPUT.dat $OUTPUT.nopheno.ped > $OUTPUT.bgl 
endif 

if ($PHASE) then
    echo "[$i] Phasing reference using Beagle (see progress in $OUTPUT.bgl.log).";  @ i++
    java -jar $SCRIPTPATH/beagle.jar unphased=$OUTPUT.bgl nsamples=4 niterations=10 missing=0 verbose=true maxwindow=1000 log=$OUTPUT.phasing >> $OUTPUT.bgl.log
endif

if ($CLEANUP) then
    rm $OUTPUT.nopheno.ped
    rm $OUTPUT.bgl.gprobs
    rm $OUTPUT.bgl.r2
    #rm $OUTPUT.bgl
    rm $OUTPUT.ped
    rm $OUTPUT.map
    rm $OUTPUT.dat
    #rm $OUTPUT.phasing.log
endif

echo "[$i] Done."
