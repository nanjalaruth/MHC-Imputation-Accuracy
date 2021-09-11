#!/opt/conda/bin/Rscript

# required packages
library(data.table)
library(tidyr)
library(dplyr)

#Read in data
frq_file <- read.csv("${frq}", sep="", header=FALSE)
rsquared_file <- read.delim("${rsquared}", header=FALSE)
id_file <- read.delim("${snpid}", header=FALSE)

#Assign column names
colnames(rsquared_file) <- c("Variants","R2")
colnames(id_file) <- "id"
colnames(frq_file) <- c("CHR", "SNP", "A1", "A2", "MAF", "NCHROBS")

#Combine the info file
info_file <- frq_file %>% 
   bind_cols(id_file, rsquared_file) %>% 
   select(id,A2,A1,MAF,R2) %>%
   dplyr::rename(SNP=id,REF=A2,ALT=A1,Rsq=R2)

# #Duplicate MAF frequency to generate ALT_Frq
info_file\$ALT_Frq<-info_file\$MAF 
# #Rearrange the columns
info_file <- info_file[, c(1, 2, 3, 6, 4, 5)]
# #Add a new column
info_file\$Genotyped <- "Imputed"


write.table(info_file, file= "${info_out}", sep = "\\t", quote = FALSE, row.names = FALSE)