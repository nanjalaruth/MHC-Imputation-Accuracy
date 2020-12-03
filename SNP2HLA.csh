#!/opt/conda/bin/tcsh

######################################################################################################
#
# SNP2HLA.csh
#
# Author: Sherman Jia (xiaomingjia@gmail.com)
#         + Small modifications by Buhm Han (buhmhan@broadinstitute.org): 8/7/12
# DESCRIPTION: This script runs imputation of HLA amino acids and classical alleles using SNP data.
# 
# INPUTS: 
# 1. Plink dataset
# 2. Reference dataset (beagle format)
# 
# DEPENDENCIES: (download and place in the same folder as this script)
# 1. PLINK (1.07)
# 2. Beagle (3.0.4)
# 3. merge_tables.pl (Perl script to merge files indexed by a specific column)
# 4. linkage2beagle and beagle2linkage (Beagle utilities for PED <-> Beagle format)
# 5. ParseDosage.csh (Converts Beagle posterior probabilities [.gprobs] to dosages in PLINK format [.dosage])
#
# USAGE: ./SNP2HLA.csh DATA (.bed/.bim/.fam) REFERENCE (.bgl.phased/.markers) OUTPUT plink max_memory[mb] window
#
######################################################################################################

if ($#argv < 4) then
    echo "USAGE: ./SNP2HLA.csh DATA (.bed/.bim/.fam) REFERENCE (.bgl.phased/.markers) OUTPUT plink {optional: java_max_memory[mb] marker_window_size}"; exit 1
endif

set SCRIPTPATH=`dirname $0`

set MERGE=$SCRIPTPATH/merge_tables.pl
set PARSEDOSAGE=$SCRIPTPATH/ParseDosage.csh

# CHECK FOR DEPENDENCIES
if (! -e `which $4`) then
    echo "Please install PLINK (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml) and point to the plink run file."; 
    echo "tcsh: use plink"
    echo "bash: use ./plink"
    exit 1
else if (! -e $SCRIPTPATH/beagle.jar) then
    echo "Please install Beagle 3 (http://faculty.washington.edu/browning/beagle/beagle.html#download) and copy the run file (beagle.3.0.4/beagle.jar) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/linkage2beagle.jar) then
    echo "Please copy linkage2beagle.jar included in the beagle 3.0.4 zip file (beagle.3.0.4/utility/linkage2beagle.jar) into $SCRIPTPATH/"; exit 1
else if (! -e $SCRIPTPATH/beagle2linkage.jar) then # We use beagle2linkage (Buhm, 8/13/12)
    echo "Please copy beagle2linkage.jar (http://faculty.washington.edu/browning/beagle_utilities/utilities.html) into $SCRIPTPATH/"; exit 1
else if (! -e $MERGE) then
    echo "Please copy merge_tables.pl (included with this package) into $SCRIPTPATH/"; exit 1
else if (! -e $PARSEDOSAGE) then
    echo "Please copy ParseDosage.csh (included with this package) into $SCRIPTPATH/"; exit 1
endif

# INPUTS
set INPUT=$1
set REFERENCE=$2
set OUTPUT=$3
set PLINK=$4

if ($#argv >= 5) then
    set MEM=$5
else
    set MEM=2000 # Default java memory 2000 Mb (2Gb)
endif

if ($#argv >= 6) then
    set WIN=$6

else
    set WIN=1000 # Default Beagle phasing/imputation window = 1000 markers
endif

#set JAVATMP=$OUTPUT.javatmpdir
#mkdir -p $JAVATMP
alias plink '$PLINK --noweb --silent --allow-no-sex'
#alias beagle 'java -Djava.io.tmpdir=$JAVATMP -Xmx$MEM\m -jar $SCRIPTPATH/beagle.jar'
#alias linkage2beagle 'java -Djava.io.tmpdir=$JAVATMP -Xmx$MEM\m -jar $SCRIPTPATH/linkage2beagle.jar'
#alias beagle2linkage 'java -Djava.io.tmpdir=$JAVATMP -Xmx$MEM\m -jar $SCRIPTPATH/beagle2linkage.jar'

# Functions to run
set EXTRACT_MHC = 1
set FLIP        = 1
set CONVERT_IN  = 1
set IMPUTE      = 1
set CONVERT_OUT = 1
set CLEANUP     = 1

# SET PARAMETERS
set TOLERATED_DIFF = .15
set i = 1

echo ""
echo "SNP2HLA: Performing HLA imputation for dataset $INPUT";
echo "- Java memory = "$MEM"Mb"
echo "- Beagle window size = "$WIN" markers"

set MHC=$OUTPUT.MHC

if ($EXTRACT_MHC) then
    echo "[$i] Extracting SNPs from the MHC."; @ i++
    plink --bfile $INPUT --chr 6 --from-mb 29 --to-mb 34 --maf 0.025 --make-bed --out $OUTPUT.MHC
endif
	
if ($FLIP) then
    echo "[$i] Performing SNP quality control."; @ i++

    # Identifying non-A/T non-C/G SNPs to flip
    echo "SNP 	POS	A1	A2" > $OUTPUT.tmp1
    cut -f2,4- $MHC.bim >> $OUTPUT.tmp1
    echo "SNP 	POSR	A1R	A2R" > $OUTPUT.tmp2
    cut -f2,4- $REFERENCE.bim >> $OUTPUT.tmp2
    $MERGE $OUTPUT.tmp2 $OUTPUT.tmp1 SNP |  grep -v -w NA > $OUTPUT.SNPS.alleles

    awk '{if ($3 != $6 && $3 != $7){print $1}}' $OUTPUT.SNPS.alleles > $OUTPUT.SNPS.toflip1
    plink --bfile $MHC --flip $OUTPUT.SNPS.toflip1 --make-bed --out $MHC.FLP

    # Calculating allele frequencies
    plink --bfile $MHC.FLP --freq --out $MHC.FLP.FRQ
    sed 's/A1/A1I/g' $MHC.FLP.FRQ.frq | sed 's/A2/A2I/g' | sed 's/MAF/MAF_I/g' > $OUTPUT.tmp

    mv $OUTPUT.tmp $MHC.FLP.FRQ
    $MERGE $REFERENCE.FRQ.frq $MHC.FLP.FRQ.frq SNP | grep -v -w NA > $OUTPUT.SNPS.frq
    sed 's/ /\t/g' $OUTPUT.SNPS.frq | awk '{if ($3 != $8){print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $9 "\t" $8 "\t" 1-$10 "\t*"}else{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $9 "\t" $10 "\t."}}' > $OUTPUT.SNPS.frq.parsed
    
    # Finding A/T and C/G SNPs
    awk '{if (($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C")){if ($4 > $7){diff=$4 - $7; if ($4 > 1-$7){corrected=$4-(1-$7)}else{corrected=(1-$7)-$4}}else{diff=$7-$4;if($7 > (1-$4)){corrected=$7-(1-$4)}else{corrected=(1-$4)-$7}};print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" diff "\t" corrected}}' $OUTPUT.SNPS.frq.parsed > $OUTPUT.SNPS.ATCG.frq

    # Identifying A/T and C/G SNPs to flip or remove
    awk '{if ($10 < $9 && $10 < .15){print $1}}' $OUTPUT.SNPS.ATCG.frq > $OUTPUT.SNPS.toflip2
    awk '{if ($4 > 0.4){print $1}}' $OUTPUT.SNPS.ATCG.frq > $OUTPUT.SNPS.toremove

    # Identifying non A/T and non C/G SNPs to remove
    awk '{if (!(($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C"))){if ($4 > $7){diff=$4 - $7;}else{diff=$7-$4}; if (diff > '$TOLERATED_DIFF'){print $1}}}' $OUTPUT.SNPS.frq.parsed >> $OUTPUT.SNPS.toremove
    awk '{if (($2 != "A" && $2 != "C" && $2 != "G" && $2 != "T") || ($3 != "A" && $3 != "C" && $3 != "G" && $3 != "T")){print $1}}' $OUTPUT.SNPS.frq.parsed >> $OUTPUT.SNPS.toremove
    awk '{if (($2 == $5 && $3 != $6) || ($3 == $6 && $2 != $5)){print $1}}' $OUTPUT.SNPS.frq.parsed >> $OUTPUT.SNPS.toremove

    # Making QCd SNP file
    plink --bfile $MHC.FLP --geno 0.2 --exclude $OUTPUT.SNPS.toremove --flip $OUTPUT.SNPS.toflip2 --make-bed --out $MHC.QC
    plink --bfile $MHC.QC --freq --out $MHC.QC.FRQ
    sed 's/A1/A1I/g' $MHC.QC.FRQ.frq | sed 's/A2/A2I/g' | sed 's/MAF/MAF_I/g' > $OUTPUT.tmp
    mv $OUTPUT.tmp $MHC.QC.FRQ.frq
    $MERGE $REFERENCE.FRQ.frq $MHC.QC.FRQ.frq SNP | grep -v -w NA > $OUTPUT.SNPS.QC.frq

    cut -f2 $OUTPUT.SNPS.QC.frq | awk '{if (NR > 1){print $1}}' > $OUTPUT.SNPS.toinclude

    echo "SNP 	POS	A1	A2" > $OUTPUT.tmp1
    cut -f2,4- $MHC.QC.bim >> $OUTPUT.tmp1

    $MERGE $OUTPUT.tmp2 $OUTPUT.tmp1 SNP | awk '{if (NR > 1){if ($5 != "NA"){pos=$5}else{pos=$2}; print "6\t" $1 "\t0\t" pos "\t" $3 "\t" $4}}' > $MHC.QC.bim

    # Recoding QC'd file as ped
    plink --bfile $MHC.QC --extract $OUTPUT.SNPS.toinclude --make-bed --out $MHC.QC.reorder
    plink --bfile $MHC.QC.reorder --recode --out $MHC.QC

    # Making SNP files
    awk '{print "M " $2}' $MHC.QC.map > $MHC.QC.dat
    cut -f2 $MHC.QC.map > $MHC.snps
    cut -d ' ' -f1-5,7- $MHC.QC.ped > $MHC.QC.nopheno.ped

    # Remove temporary files
    rm $OUTPUT.tmp1 $OUTPUT.tmp2
    rm $MHC.FLP*
    rm $MHC.QC.ped $MHC.QC.map
    rm $OUTPUT.SNPS.*
endif

if ($CONVERT_IN) then
    echo "[$i] Convering data to beagle format."; @ i++
    java -jar $SCRIPTPATH/linkage2beagle.jar $MHC.QC.dat $MHC.QC.nopheno.ped > $MHC.QC.bgl
    echo "===============================================================================" >> $OUTPUT.bgl.log

    rm $MHC.QC.reorder*
    rm $MHC.QC.nopheno.ped
endif

if ($IMPUTE) then
    echo "[$i] Performing HLA imputation (see $OUTPUT.bgl.log for progress)."; @ i++

    java -jar $SCRIPTPATH/beagle.jar markers=$REFERENCE.markers unphased=$MHC.QC.bgl phased=$REFERENCE.bgl.phased gprobs=true niterations=10 nsamples=4 missing=0 verbose=true maxwindow=$WIN out=$OUTPUT.IMPUTED log=$OUTPUT.phasing lowmem=true >> $OUTPUT.bgl.log
endif

if ($CONVERT_OUT) then
    echo "[$i] Converting posterior probabilities to PLINK dosage format."; @ i++

    # Converting .gprobs to .dosage format
    mv $OUTPUT.IMPUTED.${MHC:t}.QC.bgl.phased $OUTPUT.bgl.phased
    mv $OUTPUT.IMPUTED.${MHC:t}.QC.bgl.gprobs $OUTPUT.bgl.gprobs
    mv $OUTPUT.IMPUTED.${MHC:t}.QC.bgl.r2 $OUTPUT.bgl.r2

    $PARSEDOSAGE $OUTPUT.bgl.gprobs > $OUTPUT.dosage

    echo "[$i] Converting imputation genotypes to PLINK .ped format."; @ i++
    cat $OUTPUT.bgl.phased | java -jar $SCRIPTPATH/beagle2linkage.jar $OUTPUT.tmp # Buhm
    cut -d ' ' -f6- $OUTPUT.tmp.ped > $OUTPUT.tmp       # Buhm
    paste -d ' ' $MHC.fam $OUTPUT.tmp | tr -d "\015" > $OUTPUT.ped
    cut -f1-4 $REFERENCE.bim > $OUTPUT.map
    cp $MHC.fam $OUTPUT.fam

    # Create PLINK bed file
    plink --ped $OUTPUT.ped --map $OUTPUT.map --make-bed --out $OUTPUT
endif

if ($CLEANUP) then
    rm $OUTPUT.MHC.*
    rm $OUTPUT.tmp*
    rm $OUTPUT.IMPUTED.*.bgl.phased.phased
    rm $OUTPUT.phasing.log
    rm $OUTPUT.ped
    rm $OUTPUT.map
    rm $OUTPUT.log
    #rm -r $JAVATMP
    rm -f plink.log
    echo "DONE!"
    echo ""
endif

