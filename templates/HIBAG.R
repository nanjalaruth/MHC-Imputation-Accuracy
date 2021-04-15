#!/opt/conda/bin/Rscript
###############################################################################

#1000 GENOMES POPULATION

###############################################################################

#BUILD A HIBAG MODEL FOR HLA GENOTYPE IMPUTATION

# Load the HIBAG library
library("HIBAG")

#Set the working directory
# setwd("")

###############################################################################

## HLA-A ALLELE

###############################################################################

#STEP 1
##LOAD SNP GENOTYPES
bed.fn <- "${bed}"
fam.fn <- "${fam}"
bim.fn <- "${bim}"
geno <- hlaBED2Geno(bed.fn, fam.fn, bim.fn, assembly="hg19")

#STEP 2
##LOAD HLA TYPES
HLA_Type_Table <- read.table(file="${hlatyps}")
head(HLA_Type_Table)
dim(HLA_Type_Table)

##Convert the HLA type table fields to a vector
HLA_Type_Table\$sample.id <- as.vector(HLA_Type_Table\$sample.id)
HLA_Type_Table\$A.1 <- as.vector(HLA_Type_Table\$A.1)     
HLA_Type_Table\$A.2 <- as.vector(HLA_Type_Table\$A.2)

#STEP 3
##Train HLA A allele
hlaA.id <- "A"
train.HLA_A <- hlaAllele(HLA_Type_Table\$sample.id,
                       H1 = HLA_Type_Table[, paste(hlaA.id, ".1", sep="")],
                       H2 = HLA_Type_Table[, paste(hlaA.id, ".2", sep="")],
                       locus=hlaA.id, assembly="hg19")


#STEP 4
# Select SNPs flanking a region of 500kb on each side 
#or an appropriate flanking size without sacrificing predictive accuracy
region <- 500   # kb
snpid <- hlaFlankingSNP(geno\$snp.id,
                        geno\$snp.position, hlaA.id, region*1000, assembly="hg19")
length(snpid)


#STEP 5
# TRAIN SNP GENOTYPES
train.geno <- hlaGenoSubset(geno,
                            snp.sel = match(snpid, geno\$snp.id))
summary(train.geno)


#STEP 6
# LOAD THE HIBAG model
model.obj <- get(load("${hla_a}"))
model <- hlaModelFromObj(model.obj)
summary(model)


#STEP 7
#Import a test data file and predict its hla types
# import a PLINK BED file
gbed.fn <- ("${pbed}")
gfam.fn <- ("${pfam}")
gbim.fn <- ("${pbim}")
test.geno <- hlaBED2Geno(gbed.fn, gfam.fn, gbim.fn, assembly="hg19")


#STEP 8
# predict
pred <- hlaPredict(model, test.geno)
head(pred\$value)


#STEP 9
#Compare the predicted HLA allele to the original validation population
comp <- hlaCompareAllele(train.HLA_A, pred, allele.limit=model,
                        call.threshold=0)
comp\$overall


#STEP 10
#REPORT OVERALL ACCURACY, PER ALLELE SENSITIVITY & SPECIFICITY

hlaReport(comp, export.fn="${hibag_out}_HLA_A_accuracy.html", type="html", header=TRUE)


#STEP 11
# Visualize
hlaReportPlot(pred, fig="matching")
hlaReportPlot(model=model, fig="matching")
hlaReportPlot(pred, model=model, fig="matching")
hlaReportPlot(pred, train.HLA_A, fig="call.rate")
hlaReportPlot(pred, train.HLA_A, fig="call.threshold")


###########################################################################

## HLA-B ALLELE

##########################################################################

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
HLA_Type_Table\$sample.id <- as.vector(HLA_Type_Table\$sample.id)
HLA_Type_Table\$B.1 <- as.vector(HLA_Type_Table\$B.1)
HLA_Type_Table\$B.2 <- as.vector(HLA_Type_Table\$B.2)


#STEP 3
##Train HLA A allele
hlaB.id <- "B"
train.HLA_B <- hlaAllele(HLA_Type_Table\$sample.id,
                         H1 = HLA_Type_Table[, paste(hlaB.id, ".1", sep="")],
                         H2 = HLA_Type_Table[, paste(hlaB.id, ".2", sep="")],
                         locus=hlaB.id, assembly="hg19")

#STEP 4
# Select SNPs flanking a region of 500kb on each side 
#or an appropriate flanking size without sacrificing predictive accuracy
region <- 500   # kb
snpid <- hlaFlankingSNP(geno\$snp.id,
                        geno\$snp.position, hlaB.id, region*1000, assembly="hg19")
length(snpid)


#STEP 5
# TRAIN SNP GENOTYPES
train.geno <- hlaGenoSubset(geno,
                            snp.sel = match(snpid, geno\$snp.id))
summary(train.geno)


#STEP 6
# LOAD THE HIBAG model
model.obj <- get(load("${hla_b}"))
model <- hlaModelFromObj(model.obj)
summary(model)


#STEP 7
#Import a test data file and predict its hla types
# import a PLINK BED file
# gbed.fn <- ("${pbed}")
# gfam.fn <- ("${pfam}")
# gbim.fn <- ("${pbim}")
# test.geno <- hlaBED2Geno(gbed.fn, gfam.fn, gbim.fn, assembly="hg19")


#STEP 8
# predict
pred <- hlaPredict(model, test.geno)
head(pred\$value)


#STEP 9
#Compare the predicted HLA allele to the original validation population
comp <- hlaCompareAllele(train.HLA_B, pred, allele.limit=model,
                        call.threshold=0)
comp\$overall


#STEP 10
#REPORT OVERALL ACCURACY, PER ALLELE SENSITIVITY & SPECIFICITY

hlaReport(comp, export.fn="${hibag_out}_HLA_B_accuracy.html", type="html", header=TRUE)


#STEP 11
# visualize
hlaReportPlot(pred, fig="matching")
hlaReportPlot(model=model, fig="matching")
hlaReportPlot(pred, model=model, fig="matching")
hlaReportPlot(pred, train.HLA_B, fig="call.rate")
hlaReportPlot(pred, train.HLA_B, fig="call.threshold")


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
HLA_Type_Table\$sample.id <- as.vector(HLA_Type_Table\$sample.id)
HLA_Type_Table\$C.1 <- as.vector(HLA_Type_Table\$C.1)
HLA_Type_Table\$C.2 <- as.vector(HLA_Type_Table\$C.2)

#STEP 3
##Train HLA C allele
hlaC.id <- "C"
train.HLA_C <- hlaAllele(HLA_Type_Table\$sample.id,
                         H1 = HLA_Type_Table[, paste(hlaC.id, ".1", sep="")],
                         H2 = HLA_Type_Table[, paste(hlaC.id, ".2", sep="")],
                         locus=hlaC.id, assembly="hg19")

#STEP 4
# Select SNPs flanking a region of 500kb on each side 
#or an appropriate flanking size without sacrificing predictive accuracy
region <- 500   # kb
snpid <- hlaFlankingSNP(geno\$snp.id,
                        geno\$snp.position, hlaC.id, region*1000, assembly="hg19")
length(snpid)


#STEP 5
# TRAIN SNP GENOTYPES
train.geno <- hlaGenoSubset(geno,
                            snp.sel = match(snpid, geno\$snp.id))
summary(train.geno)


#STEP 6
# LOAD THE HIBAG model
model.obj <- get(load("${hla_b}"))
model <- hlaModelFromObj(model.obj)
summary(model)


#STEP 7
#Import a test data file and predict its hla types
# import a PLINK BED file
# gbed.fn <- ("${pbed}")
# gfam.fn <- ("${pfam}")
# gbim.fn <- ("${pbim}")
# test.geno <- hlaBED2Geno(gbed.fn, gfam.fn, gbim.fn, assembly="hg19")


#STEP 8
# predict
pred <- hlaPredict(model, test.geno)
head(pred\$value)


#STEP 9
#Compare the predicted HLA allele to the original validation population
comp <- hlaCompareAllele(train.HLA_C, pred, allele.limit=model,
                        call.threshold=0)
comp\$overall


#STEP 10
#REPORT OVERALL ACCURACY, PER ALLELE SENSITIVITY & SPECIFICITY

hlaReport(comp, export.fn="${hibag_out}_HLA_C_accuracy.html", type="html", header=TRUE)


#STEP 11
# visualize
hlaReportPlot(pred, fig="matching")
hlaReportPlot(model=model, fig="matching")
hlaReportPlot(pred, model=model, fig="matching")
hlaReportPlot(pred, train.HLA_C, fig="call.rate")
hlaReportPlot(pred, train.HLA_C, fig="call.threshold")