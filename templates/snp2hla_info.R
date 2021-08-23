#!/opt/conda/bin/Rscript

# required packages
library(data.table)
library(tidyr)
library(dplyr)

#Read in data
frq_file <- read.csv(file = "${frq}", sep="")
rsquared_file <- read.delim(file = "${rsquared}", header=FALSE)
id_file <- read.delim(file = "${snpid}", header=FALSE)

#Assign column names to rsquared file
colnames(rsquared_file) <- c("SNP","R2")

#Combine the info file
info_file <- frq_file %>% 
   inner_join(rsquared_file, by = "SNP")  %>% 
   bind_cols(id_file) %>% 
   select(V1,A2,A1,MAF,R2) %>% 
   dplyr::rename(ALT=A1,REF=A2,SNP=V1,Rsq=R2)

write.table(info_file, file= "${info_out}", sep = "\\t", quote = FALSE, row.names = FALSE)