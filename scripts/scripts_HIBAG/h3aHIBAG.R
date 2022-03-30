#!/opt/conda/bin/Rscript
###############################################################################

#H3A POPULATION

###############################################################################

#BUILD A HIBAG MODEL FOR HLA GENOTYPE IMPUTATION

# Load the HIBAG library
library("HIBAG")

###############################################################################

## HLA-A ALLELE

###############################################################################

#STEP 1
##LOAD SNP GENOTYPES
bed.fn <- "./H3Africa_H3Africa_geno.bed"
fam.fn <- "./H3Africa_H3Africa_geno.fam"
bim.fn <- "./H3Africa_H3Africa_geno.bim"
geno <- hlaBED2Geno(bed.fn, fam.fn, bim.fn, assembly="hg19")

#STEP 2
##LOAD HLA TYPES
HLA_Type_Table <- read.table(file="./H3Africa_hibag_hlatypes.edited")
head(HLA_Type_Table)
dim(HLA_Type_Table)

##Convert the HLA type table fields to a vector
HLA_Type_Table$sample.id <- as.vector(HLA_Type_Table$sample.id)
HLA_Type_Table$A.1 <- as.vector(HLA_Type_Table$A.1)     
HLA_Type_Table$A.2 <- as.vector(HLA_Type_Table$A.2)

#STEP 3
##Train HLA A allele
hlaA.id <- "A"
train.HLA_A <- hlaAllele(HLA_Type_Table$sample.id,
                       H1 = HLA_Type_Table[, paste(hlaA.id, ".1", sep="")],
                       H2 = HLA_Type_Table[, paste(hlaA.id, ".2", sep="")],
                       locus=hlaA.id, assembly="hg19")


#STEP 4
# Select SNPs flanking a region of 500kb on each side 
#or an appropriate flanking size without sacrificing predictive accuracy
region <- 500   # kb
snpid <- hlaFlankingSNP(geno$snp.id,
                        geno$snp.position, hlaA.id, region*1000, assembly="hg19")
length(snpid)


#STEP 5
# TRAIN SNP GENOTYPES
train.geno <- hlaGenoSubset(geno,
                            snp.sel = match(snpid, geno$snp.id))
summary(train.geno)

# #Build model in parallel
library(parallel)
cl <- makeCluster(20)
set.seed(1000)
hlaParallelAttrBagging(cl, train.HLA_A, train.geno, nclassifier=100, 
                      auto.save = "h3a_h3a_HLA_A_Model.RData", stop.cluster=TRUE)

# #Release model
model.obj <- get(load("h3a_h3a_HLA_A_Model.RData"))
model <- hlaModelFromObj(model.obj)
summary(model)


# ###########################################################################

# ## HLA-B ALLELE

# ##########################################################################

# #STEP 1
# ##LOAD SNP GENOTYPES
# # bed.fn <- ("./1kg_all_gambian_geno.bed")
# # fam.fn <- ("./1kg_all_gambian_geno.fam")
# # bim.fn <- ("./1kg_all_gambian_geno.bim")
# # geno <- hlaBED2Geno(bed.fn, fam.fn, bim.fn, assembly="hg19")


# # #STEP 2
# # ##LOAD HLA TYPES
# # HLA_Type_Table <- read.table(file="./gambian_hibag_HLAtype")
# # head(HLA_Type_Table)
# # dim(HLA_Type_Table)

##Convert the HLA type table fields to a vector
HLA_Type_Table$sample.id <- as.vector(HLA_Type_Table$sample.id)
HLA_Type_Table$B.1 <- as.vector(HLA_Type_Table$B.1)
HLA_Type_Table$B.2 <- as.vector(HLA_Type_Table$B.2)


#STEP 3
##Train HLA A allele
hlaB.id <- "B"
train.HLA_B <- hlaAllele(HLA_Type_Table$sample.id,
                         H1 = HLA_Type_Table[, paste(hlaB.id, ".1", sep="")],
                         H2 = HLA_Type_Table[, paste(hlaB.id, ".2", sep="")],
                         locus=hlaB.id, assembly="hg19")

#STEP 4
# Select SNPs flanking a region of 500kb on each side 
#or an appropriate flanking size without sacrificing predictive accuracy
region <- 500   # kb
snpid <- hlaFlankingSNP(geno$snp.id,
                        geno$snp.position, hlaB.id, region*1000, assembly="hg19")
length(snpid)


#STEP 5
# TRAIN SNP GENOTYPES
train.geno <- hlaGenoSubset(geno,
                            snp.sel = match(snpid, geno$snp.id))
summary(train.geno)


#Build HIBAG model
library(parallel)
cl <- makeCluster(20)
set.seed(1000)
hlaParallelAttrBagging(cl, train.HLA_B, train.geno, nclassifier=100, 
                      auto.save = "h3a_h3a_HLA_B_Model.RData", stop.cluster=TRUE)


#Release model
model.obj <- get(load("h3a_h3a_HLA_B_Model.RData"))
model <- hlaModelFromObj(model.obj)
summary(model)


###########################################################################

## HLA-C ALLELE IMPUTATION

############################################################################

#STEP 1
##LOAD SNP GENOTYPES
# bed.fn <- ("./1kg_all_gambian_geno.bed")
# fam.fn <- ("./1kg_all_gambian_geno.fam")
# bim.fn <- ("./1kg_all_gambian_geno.bim")
# geno <- hlaBED2Geno(bed.fn, fam.fn, bim.fn, assembly="hg19")


# #STEP 2
# ##LOAD HLA TYPES
# HLA_Type_Table <- read.table(file="./gambian_hibag_HLAtype")
# head(HLA_Type_Table)
# dim(HLA_Type_Table)

##Convert the HLA type table fields to a vector
HLA_Type_Table$sample.id <- as.vector(HLA_Type_Table$sample.id)
HLA_Type_Table$C.1 <- as.vector(HLA_Type_Table$C.1)
HLA_Type_Table$C.2 <- as.vector(HLA_Type_Table$C.2)

#STEP 3
##Train HLA C allele
hlaC.id <- "C"
train.HLA_C <- hlaAllele(HLA_Type_Table$sample.id,
                         H1 = HLA_Type_Table[, paste(hlaC.id, ".1", sep="")],
                         H2 = HLA_Type_Table[, paste(hlaC.id, ".2", sep="")],
                         locus=hlaC.id, assembly="hg19")

#STEP 4
# Select SNPs flanking a region of 500kb on each side 
#or an appropriate flanking size without sacrificing predictive accuracy
region <- 500   # kb
snpid <- hlaFlankingSNP(geno$snp.id,
                        geno$snp.position, hlaC.id, region*1000, assembly="hg19")
length(snpid)


#STEP 5
# TRAIN SNP GENOTYPES
train.geno <- hlaGenoSubset(geno,
                            snp.sel = match(snpid, geno$snp.id))
summary(train.geno)

#BUILD AND PREDICT IN PARALLEL
library(parallel)
cl <- makeCluster(20)
set.seed(1000)
hlaParallelAttrBagging(cl, train.HLA_C, train.geno, nclassifier=100, 
                      auto.save = "h3a_h3a_HLA_C_Model.RData", stop.cluster=TRUE)

#Release HIBAG model                      
model.obj <- get(load("h3a_h3a_HLA_C_Model.RData"))
model <- hlaModelFromObj(model.obj)
summary(model)