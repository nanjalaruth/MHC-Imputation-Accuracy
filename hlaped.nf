// Merge HLA types
// SNP2HLA version
// All African population

 Channel.fromPath(params.files).set { samples_ch }
 
 process snp2hlatypes {
   publishDir "$params.outdir"
   input:
     file smp from samples_ch.collect()
   output:
     file "snp2hla.PED" into afrsnp2hlatypes
   script:
   """
  # extract tsv file with the first column as the file name
   awk -F, 'NR==1{print "file_name",\$0;next}{print FILENAME , \$0}' $smp | cut -f 1-7 | sed 's/*//g;s/://g' > sample1.txt
   #extract the 2nd line only
   sed -i -n 'n;p' sample1.txt
  #replace the _result.tsv extension in the filename with a tab followed by a 1 (sex)
   sed -i 's/_result.tsv/\t1\t/g' sample1.txt
  # replace any missing value in the columns with a 0
   awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if(\$i ~ /^ *\$/) \$i = 0 }; 1' sample1.txt > sample.txt
  # append 10 columns at the end with 0s
   for i in {1..10}; do sed -i "s/\$/\t0/" sample.txt; done
  # replace the A-Z in the 4th-9th column with nothing
   awk '{OFS="\t";gsub("[A-Z]","",\$4);gsub("[A-Z]","",\$5);gsub("[A-Z]","",\$6);gsub("[A-Z]","",\$7);gsub("[A-Z]","",\$8); gsub("[A-Z]","",\$9)}1' sample.txt >newfile.txt
  # duplicate the 1st column 4 times
   awk -v col=1 -v n=4 'function repeat(v, n,i){for(i=1; i<=n; i++)printf("%s%s",(i==1?"":OFS="\t"),v)}{for(i=1; i<=NF; i++)printf("%s%s",(i==col?repeat(\$i,n):\$i),i==NF?RS:OFS)}' newfile.txt > snp2hla.PED
   
   """
 }
 
// Gambian subpopulation

 Channel.fromPath(params.gwdids).set {ids_ch }
 process gwdsnp2hlatypes {
   publishDir "$params.outdir"
   input:
     file smple from ids_ch
     file g from afrsnp2hlatypes
   output:
     file "gwdsnp2hla.PED" into gwdsnp2hlatypes
   script:
   """
   #extract gambian population
   grep -f $smple $g > gwdsnp2hla.PED
   
   """
 }


// HIBAG version
// All African population

 Channel.fromPath(params.files).set { smp_ch }
 process hibaghlatypes {
   publishDir "$params.outdir"
   input:
     file smp from smp_ch.collect()
   output:
     file "hibag_HLA_Type" into hibaghlatypes
   script:
   """
  # extract tsv file with the first column as the file name
   awk -F, 'NR==1{print "file_name",\$0;next}{print FILENAME , \$0}' $smp | cut -f 1-7 | sed 's/*//g' > sample1.txt
  # extract the 2nd line only
   sed -i -n 'n;p' sample1.txt
  # replace the _result.tsv extension in the filename with nothing
   sed -i 's/_result.tsv//g' sample1.txt
  # replace the A-Z in the 4th-9th column with nothing
   awk '{OFS="\t";gsub("[A-Z]","",\$3);gsub("[A-Z]","",\$4);gsub("[A-Z]","",\$5);gsub("[A-Z]","",\$6);gsub("[A-Z]","",\$7); gsub("[A-Z]","",\$8)}1' sample1.txt >newfile.txt
  # replace any missing value in the columns with a NA
   awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if(\$i ~ /^ *\$/) \$i = "<NA>" }; 1' newfile.txt > sample1.txt
  # append 10 columns at the end with 0s
   for i in {1..6}; do sed -i "s/\$/\t<NA>/" sample1.txt; done
  # remove second column (phenotype)
   awk '{OFS="\t";for(i=1;i<=NF;i++)if(i!=2)printf\$i OFS;print""}' sample1.txt > test1.txt
  # add the rownames
   awk '{OFS="\t";print NR,\$0}' test1.txt > smp
  # add headers
   ( echo -e "\tsample.id\tA.1\tA.2\tB.1\tB.2\tC.1\tC.2\tDQA1.1\tDQA1.2\tDQB1.1\tDQB1.2\tDRB1.1\tDRB1.2"; cat smp ) > hibag_HLA_Type
  
   """
 }

// Gambian subpopulation

 Channel.fromPath(params.gwdids).set {id_ch }
 process gwdhibaghlatypes {
   publishDir "$params.outdir"
   input:
     file alp from id_ch
     file h from hibaghlatypes
   output:
     file "gwd_hibag_HLA_Type"
   script:
   """
   # extract the gambian subpopulation
   grep -f $alp $h > gwd_hi
   # remove the distorted first column
   awk '{OFS="\t";for(i=1;i<=NF;i++)if(i!=1)printf\$i OFS;print""}' gwd_hi > t1.txt
   # introduce a new first column
   awk '{OFS="\t";print NR,\$0}' t1.txt > trial
   # introduce the headers
   ( echo -e "\tsample.id\tA.1\tA.2\tB.1\tB.2\tC.1\tC.2\tDQA1.1\tDQA1.2\tDQB1.1\tDQB1.2\tDRB1.1\tDRB1.2"; cat trial ) > gwd_hibag_HLA_Type
   
   """
 }

// SNP genotypes

 Channel.fromPath(params.reads).set { rds_ch }
 Channel.fromPath(params.gids).set {gids_ch}
 Channel.fromPath(params.aids).set {aids_ch}
 process snpgenotypes {
   publishDir "$params.outdir"
   input:
     file rds from rds_ch
     file c from gids_ch
     file t from aids_ch
   output:
     file "GWD.*" into gwdgenotypes
     file "AFR.*" into afrgenotypes
   script:
   """
   #convert vcf to plink format
   plink2 --vcf $rds --make-bed --max-alleles 2 --out MHC
   #remove duplicated snps
   cut -f 2 MHC.bim | sort | uniq -d > 1.dups
   plink2 --bfile MHC --exclude 1.dups --make-bed --out MHC.filt
   #extract GWD plink datasets
   plink2 --bfile MHC.filt --keep $c --make-bed --out GWD
   rm GWD.log
   #extract AFR plink datasets
   plink2 --bfile MHC.filt --keep $t --make-bed --out AFR
   rm AFR.log
   """
 }



// TRAIN REFERENCE PANEL
// SNP2HLA
// GAMBIANS

// Channel.fromPath(params.makereference).set {csh_ch }
process gwd_snp2hlarefence {
   publishDir "$params.outdir"
   script:
   """
   MakeReference.csh HAPMAP_CEU HAPMAP_CEU_HLA.ped ref plink
   """
 } 

// ALL AFRICANS

// Channel.fromPath(params.makereference).set {sh_ch }
// process afr_snp2hlarefence {
   // publishDir "$params.outdir"
   // input:
     // file d from afrgenotypes
     // file b from afrsnp2hlatypes
  //  output:
     // file "afrreference"
   // script:
  // """
   // MakeReference.csh $d $b afrreference plink
   // """
 // }
