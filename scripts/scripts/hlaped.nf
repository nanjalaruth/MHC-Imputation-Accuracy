params.reads = 'optitype/*/*.tsv'
params.outdir = 'my-results'

// SNP2HLA version
Channel.fromPath(params.reads).set { samples_ch }
process snp2hlatypes {
  publishDir "$params.outdir"

  input:
    file smp from samples_ch.collect()
  output:
    file "snp2hla.PED"
  script:
  """
  #extract tsv file with the first column as the file name
  awk -F, 'NR==1{print "file_name",\$0;next}{print FILENAME , \$0}' $smp | cut -f 1-7 | sed 's/*//g;s/://g' > sample1.txt
  #extract the 2nd line only
  sed -i -n 'n;p' sample1.txt
  #replace the _result.tsv extension in the filename with a tab followed by a 1 (sex)
  sed -i 's/_result.tsv/\t1\t/g' sample1.txt
  #replace any missing value in the columns with a 0
  awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if(\$i ~ /^ *\$/) \$i = 0 }; 1' sample1.txt > sample.txt
  #append 10 columns at the end with 0s
  for i in {1..10}; do sed -i "s/\$/\t0/" sample.txt; done
  #replace the A-Z in the 4th-9th column with nothing
  awk '{OFS="\t";gsub("[A-Z]","",\$4);gsub("[A-Z]","",\$5);gsub("[A-Z]","",\$6);gsub("[A-Z]","",\$7);gsub("[A-Z]","",\$8); gsub("[A-Z]","",\$9)}1' sample.txt >newfile.txt
  #duplicate the 1st column 4 times
  awk -v col=1 -v n=4 'function repeat(v, n,i){for(i=1; i<=n; i++)printf("%s%s",(i==1?"":OFS="\t"),v)}{for(i=1; i<=NF; i++)printf("%s%s",(i==col?repeat(\$i,n):\$i),i==NF?RS:OFS)}' newfile.txt > snp2hla.PED
  
  """
}

// HIBAG version
Channel.fromPath(params.reads).set { smp_ch }
process hibaghlatypes {
  publishDir "$params.outdir"

  input:
    file smp from smp_ch.collect()
  output:
    file "hibag_HLA_Type"
  script:
  """
  #extract tsv file with the first column as the file name
  awk -F, 'NR==1{print "file_name",\$0;next}{print FILENAME , \$0}' $smp | cut -f 1-7 | sed 's/*//g' > sample1.txt
  #extract the 2nd line only
  sed -i -n 'n;p' sample1.txt
  #replace the _result.tsv extension in the filename with nothing
  sed -i 's/_result.tsv//g' sample1.txt
  #replace the A-Z in the 4th-9th column with nothing
  awk '{OFS="\t";gsub("[A-Z]","",\$3);gsub("[A-Z]","",\$4);gsub("[A-Z]","",\$5);gsub("[A-Z]","",\$6);gsub("[A-Z]","",\$7); gsub("[A-Z]","",\$8)}1' sample1.txt >newfile.txt
  #replace any missing value in the columns with a NA
  awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if(\$i ~ /^ *\$/) \$i = "<NA>" }; 1' newfile.txt > sample1.txt
  #append 10 columns at the end with 0s
  for i in {1..6}; do sed -i "s/\$/\t<NA>/" sample1.txt; done
  #remove second column (phenotype)
  awk '{OFS="\t";for(i=1;i<=NF;i++)if(i!=2)printf\$i OFS;print""}' sample1.txt > test1.txt
  #add the rownames
  awk '{OFS="\t";print NR,\$0}' test1.txt > smp
  #add headers
  ( echo -e "\tsample.id\tA.1\tA.2\tB.1\tB.2\tC.1\tC.2\tDQA1.1\tDQA1.2\tDQB1.1\tDQB1.2\tDRB1.1\tDRB1.2"; cat smp ) > hibag_HLA_Type
  
  """
}
