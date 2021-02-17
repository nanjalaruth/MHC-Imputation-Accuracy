############################################################################################
#
# MakeReference: Creates custom reference panel given HLA types and SNPs for HLA imputation
# Author: Sherman Jia (xiaomingjia@gmail.com)
#
###########################################################################################

Thank you for downloading the MakeReference package. To use this package:

1. Download Plink for your platform (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml). Copy the "plink" run file into the current directory (with MakeReference.csh).

2. Download Beagle (version 3.0.4) .jar files into the current directory.
- "beagle.jar" from http://faculty.washington.edu/browning/beagle/beagle.html#download
- "beagle2linkage.jar" and "linkage2beagle.jar" from  http://faculty.washington.edu/browning/beagle_utilities/utilities.html
We recommend downloading version 3.0.4 for the compatability issue, even if it is not the newest version. Beagle web page described above includes links for all past-version binaries.

3. Run MakeReference with sample data provided (92 CEU founders from HapMap  using the following command:

tcsh: ./MakeReference.csh HAPMAP_CEU HAPMAP_CEU_HLA.ped HM_CEU_REF plink
bash: ./MakeReference.csh HAPMAP_CEU HAPMAP_CEU_HLA.ped HM_CEU_REF ./plink

----------------------------------------------------------------------------------

Files included in this package:

1. MakeReference.csh: Creates a reference panel using SNP and HLA data.
2. HLAtoSequences.pl: Generates amino acid or intragenic SNP variants given HLA alleles
3. encodeHLA.pl: Encodes HLA alleles into binary markers (for
imputation with posterior probabilities using Beagle)
4. encodeVariants.pl: Encodes genetic positions for polymorphic HLA variants
5. Sample HapMap CEU dataset (n=180) with HLA alleles
6. HLA amino acid / protein sequence dictionary: HLA_DICTIONARY_AA.map/txt
7. HLA intragenic SNPs / DNA sequence dictionary: HLA_DICTIONARY_SNP.map/txt

----------------------------------------------------------------------------------

Input files:

1. SNP dataset
2. HLA types for SNP dataset
3. Plink

Output files:

1. New reference panel in Plink format, containing SNPs, HLA alleles, HLA amino acids, and HLA intragenic SNPs (.bed/bim/fam)
2. Allele frequences for all variants in reference panel (.FRQ.frq)
3. Beagle format phased haplotypes (.bgl.phased)
4. A file denoting marker positions and order (.markers) for subsequent imputation

----------------------------------------------------------------------------------

Marker Nomenclature: For binary encodings, P = Present, A = Absent.

1. Classical HLA alleles: HLA_[GENE]_[ALLELE]. 
- HLA_C_0304 = HLA-C:03:04 (four-digit allele)
- HLA_DRB1_07 = HLA-DRB1:07 (two-digit allele)

2. HLA Amino Acids: AA_[GENE]_[AMINO ACID POSITION]_[GENETIC POSITION]_[ALLELE]. 
- AA_A_56_30018678_G = amino acid 56 of HLA-A, genetic position 30018678 (center of codon), allele = G (Gly) of multi-allelic position
- AA_C_291_31345793 = amino acid 291 of HLA-C, genetic position 31345793, bi-allelic (check {OUTPUT}.bim for alleles, P = Present, A = Absent)

3. HLA intragenic SNPS: SNP_[GENE]_[POSITION]_[ALLELE]
- SNP_B_31430319_G = SNP at position 31430319 of HLA-B, allele = G (guanine) of multi-allelic position
- SNP_DRB1_32659974 = SNP at position 32659974 of HLA-DRB1, bi-allelic (check {OUTPUT}.bim for alleles, P = Present, A = Absent)
- SNP_DQB1_32740666_AT = SNP at position 32740666 of HLA-DQB1, alleles = A (adenine) or T (thymine), (check {OUTPUT}.bim for alleles, P = Present, A = Absent)

4. Insertions / deletions: [VARIANT]_[GENE]_[POSITION]_[INSERTION/x=DELETION]
- AA_C_339_31345102_x = deletion at amino acid 339 in HLA-C, genetic position 31345102 (center of codon), (check {OUTPUT}.bim for alleles, P = deletion Present, A = deletion Absent)
- INS_C_295x296_31345779_VLAVLA = insertion between amino acids 295 and 296 of HLA-C, amino acid sequence inserted = VLAVLA, (check {OUTPUT}.bim for alleles, P = insertion Present, A = insertion Absent)
- SNP_DQA1_32717217_x = deletion at genetic position 32717217 of HLA-DQA1, (check {OUTPUT}.bim for alleles, P = deletion Present, A = deletion Absent)
