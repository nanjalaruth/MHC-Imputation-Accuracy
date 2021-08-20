#!/users/nanje/miniconda3/bin/Rscript

# required packages
library(optparse)
library(tidyverse)
library(data.table)

# prepares snp2hla info file for comparison with other imputation tools
option_list <- list(
make_option(c("-f", "--frq"), action="store", default= "${frq}", type='list',
              help = "Imputation .frq file"),
make_option(c("-r", "--rsquared"), action="store", default = "${rsquared}", type = 'list',
              help = "Imputation .rsquared file"),
make_option(c("-s", "--snpid"), action = "store", default = "${snpid}", type = 'list',
              help = "Imputation .snpid file"),
make_option(c("-o", "--info_out"), action = "store", default = "${info_out}", type = 'list',
              help = "Combined info file")
)
args <- parse_args(OptionParser(option_list = option_list))

#Read in data
frq_file <- read.csv(args\$f, sep="")
rsquared_file <- read.delim(args\$r, header=FALSE)
id_file <- read.delim(args\$s, header=FALSE)

#Assign column names to rsquared file
colnames(rsquared_file) <- c("SNP","R2")

#Combine the info file
info_file <- frq_file %>% 
   inner_join(rsquared_file, by = "SNP")  %>% 
   bind_cols(id_file) %>% 
   select(V1,A2,A1,MAF,R2) %>% 
   rename(ALT=A1,REF=A2,SNP=V1,Rsq=R2)

write.table(info_file, file= args[4], sep = "\\t", quote = FALSE, row.names = FALSE)