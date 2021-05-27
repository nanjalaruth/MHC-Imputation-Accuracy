#! /bin/bash
#while IFS=$'\t' read -r acc samp sec_samp xxx pop; do \
#sed -i 's|'"${samp}"'|'"${sec_samp}"'|g' GGVP.snp2hla.ped; done < ggvp_metadata.tsv
#sort -k 2 GGVP.snp2hla.ped > sortedGGVP.snp2hla.ped

while IFS=$'\t' read -r acc samp sec_samp xxx pop; do \
sed -i 's|'"${samp}"'|'"${sec_samp}"'|g' GGVP_hibag_hlatypes; done < ggvp_metadata.tsv
sort -k 2 GGVP_hibag_hlatypes > sortedGGVP_hibag_hlatypes
