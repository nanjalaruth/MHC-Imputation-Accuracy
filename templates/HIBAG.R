#!/opt/conda/bin/Rscript
###############################################################################

# Load the HIBAG library
library("HIBAG")

model.objA <- get(load("${model_a}"))
model.objB <- get(load("${model_b}"))
model.objC <- get(load("${model_c}"))

model.a <- hlaModelFromObj(model.objA)
model.b <- hlaModelFromObj(model.objB)
model.c <- hlaModelFromObj(model.objC)


#Import a test data file and predict its hla types
# import a PLINK BED file
bed.fn <- ("${bed}")
fam.fn <- ("${fam}")
bim.fn <- ("${bim}")
test.geno <- hlaBED2Geno(bed.fn, fam.fn, bim.fn, assembly="hg19")


A.guess <- predict(model.a, test.geno, type="response+prob", match.type="Position")
B.guess <- predict(model.b, test.geno, type="response+prob", match.type="Position")
C.guess <- predict(model.c, test.geno, type="response+prob", match.type="Position")

output.table = data.frame(SampleID=A.guess\$value\$sample.id, A_allele1=A.guess\$value\$allele1, A_allele2=A.guess\$value\$allele2, A_prob=A.guess\$value\$prob,
               B_allele1=B.guess\$value\$allele1, B_allele2=B.guess\$value\$allele2, B_prob=B.guess\$value\$prob,
               C_allele1=C.guess\$value\$allele1, C_allele2=C.guess\$value\$allele2, C_prob=C.guess\$value\$prob)

output.file = paste("${hibag_out}_HIBAG_HLA.txt",sep="")
write.table(output.table, output.file, row.names=F, sep="\t", quote=F)
