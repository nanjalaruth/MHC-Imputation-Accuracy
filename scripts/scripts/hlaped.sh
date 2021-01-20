#extracting the HLA, removing the colon, letters, asterisk ansd space
tail -n 1  */*.tsv | grep -v "==" | cut -f 1-7|sed 's/[A-Z]//g;s/*//g;s/://g' | awk '$0' > hla.txt

#add 10 columns at the end with 0s
for i in {1..10}; do sed -i "s/$/\t0/" hla.txt; done

#Append the gender 
paste sex.txt hla.txt > hla1.txt

#add 2 columns (paternalid & maternalid) with 0s
for i in {1..2}; do sed -i "s/^/0\t/" hla1.txt; done

#append sample/individualids
paste sampleid.txt hla1.txt > hla2.txt

#add 1st column (familyid) with 0s
sed -i "s/^/0\t/" hla2.txt
