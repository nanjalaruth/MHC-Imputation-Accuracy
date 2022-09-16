library(tidyverse)
library(data.table)

#set working dir
setwd("C:/Users/rnanjala/MHC-Imputation-Accuracy")

##############################################################################
#####Omni array

#Obtain the omni array data subsetted from the ggvp dataset
omni_ggvp <- fread("./ggvp_omni.txt",
                   header = F, data.table = F) %>% 
    setnames(.,colnames(.), "id")

##Obtain the Omni SNP positions
omni <- fread("./omni.pos", 
              header = F, data.table = F) %>% 
  setnames(.,colnames(.), c("chr", "start", "stop")) %>% 
  select(chr,start)
omni$id <- paste(omni$chr, omni$start, sep="_")
new_omni <- omni %>% select(id)


#Obtain the omni array data (MHC region)
#omni <- fread("./omni.pos", 
 #             header = F, data.table = F) %>% 
  #setnames(.,colnames(.), c("chr", "start", "stop")) %>% 
  #filter(chr == 6) %>%
  #filter(start >= 27479609 & start <= 34444100) %>%
  #select(chr,start)
#omni$id <- paste(omni$chr, omni$start, sep="_")
#new_omni <- omni %>% select(id)

#write the results to an output file
write.table(omni_ggvp, file="./omni_ggvp_SNP_pos", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(new_omni, file="./omni_SNP_pos", sep = "\t", quote = FALSE, row.names = FALSE)

library(UpSetR)

movies <- read.csv("./omni_SNP_pos.csv","./omni_ggvp_SNP_pos.csv", header = T, sep = ",")

#plot the venn diagram
#install.packages("VennDiagram")
require(VennDiagram)

#Read in the data
omni_ggvp_SNP_pos <- readLines("./omni_ggvp_SNP_pos")
omni_ggvp_SNP_pos
omni_SNP_pos <- readLines("./omni_SNP_pos")
omni_SNP_pos

venn.plot <- venn.diagram(x = list(omni_SNP_pos, omni_ggvp_SNP_pos),
                          category.names = c("Omni_array", "omni_ggvp_array"),
                          fill = c("darkorchid1", "cyan"),
                          cat.fontface = 2, 
                          lty =0, scaled = TRUE,
                          filename = NULL,
                          main = "Omni array SNPs vs the extracted ggvp omni array snps");
cols <- c("darkorchid1", "cyan")
lg <- legendGrob(labels=c("omni_array","omni_ggvp_array"), pch=rep(19),
                 gp=gpar(col=cols))

library(gridGraphics)
library(gridExtra)
jpeg("./plots/omni.jpg");
grab_grob <- function(){grid.echo();grid.grab()}
grid.draw(venn.plot)
g <- grab_grob()
grid.arrange(g,lg,ncol=2,widths=grid::unit(c(0.7,0.3),"npc"))
dev.off()





###########################################################################
#H3Africa 

####GGVp h3a snp positions
h3a_ggvp <- fread("./ggvp_h3a.txt",
                   header = F, data.table = F) %>%
  setnames(.,colnames(.), "id")


#####H3A SNP positions
h3a <- fread("./h3Africa_v1_chip.pos", 
             header = F, data.table = F) %>% 
  setnames(.,colnames(.), c("chr", "start", "stop")) %>% 
  select(chr,start)
h3a$id <- paste(h3a$chr, h3a$start, sep="_")
new_h3a <- h3a %>% select(id)
new_h3a

####H3A SNP positions (MHC)
#h3a <- fread("./h3Africa_v1_chip.pos", 
 #             header = F, data.table = F) %>% 
  #setnames(.,colnames(.), c("chr", "start", "stop")) %>% 
  #filter(chr == 6) %>%
  #filter(start >= 27479609 & start <= 34444100) %>%
  #select(chr,start)
#h3a$id <- paste(h3a$chr, h3a$start, sep="_")
#new_h3a <- h3a %>% select(id)
#new_h3a


#write the results to an output file
write.table(h3a_ggvp, file="./h3a_ggvp_SNP_pos", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(new_h3a, file="./h3a_SNP_pos", sep = "\t", quote = FALSE, row.names = FALSE)

#Read in the data
h3a_ggvp_SNP_pos <- readLines("./h3a_ggvp_SNP_pos")
h3a_SNP_pos <- readLines("./h3a_SNP_pos")

venn.plot <- venn.diagram(list(h3a_array=h3a_SNP_pos, h3a_ggvp_array=h3a_ggvp_SNP_pos),
                          fill = c("cornflowerblue", "green"),
                          cat.fontface = 2, 
                          lty =0, scaled = TRUE, 
                          filename = NULL,
                          main = "H3Africa array SNPs vs the extracted ggvp H3Africa array snps");

cols <- c("cornflowerblue", "green")
lg <- legendGrob(labels=c("h3a_array","h3a_ggvp_array"), pch=rep(19),
                 gp=gpar(col=cols))

library(gridGraphics)
library(gridExtra)
jpeg("./plots/H3Africa.jpg");
grab_grob <- function(){grid.echo();grid.grab()}
grid.draw(venn.plot)
g <- grab_grob()
grid.arrange(g,lg,ncol=2,widths=grid::unit(c(0.7,0.3),"npc"))
dev.off()

######################################################################################
#### H3A & Omni
venn.plot <- venn.diagram(list(h3a_array=h3a_SNP_pos, h3a_ggvp_array=h3a_ggvp_SNP_pos, omni_array=omni_SNP_pos,omni_ggvp_array=omni_ggvp_SNP_pos),
                          fill = c( "green", "yellow", "darkorchid1", "cyan"),
                          filename = NULL, scaled = TRUE, cat.fontface = 2, 
                          lty =0,
                          main = "H3Africa/Omni array SNPs vs extracted ggvp H3Africa/Omni array snps");
cols <- c("green", "yellow", "darkorchid1", "cyan")
lg <- legendGrob(labels=c("h3a_array", "h3a_ggvp_array", "omni_array","omni_ggvp_array"), pch=rep(19,length(c("omni_SNP_pos","omni_ggvp_SNP_pos"))),
                 gp=gpar(col=cols))

library(gridGraphics)
library(gridExtra)
jpeg("./plots/H3A_omni.jpg");
grab_grob <- function(){grid.echo();grid.grab()}
grid.draw(venn.plot)
g <- grab_grob()
grid.arrange(g,lg,ncol=2,widths=grid::unit(c(0.7,0.3),"npc"))
dev.off()


library(UpSetR)
library(ggplot2)

movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
                   header = T, sep = ";")

movies <- read.csv(readLines("./omni_ggvp_SNP_pos.csv","./omni_SNP_pos"), header = T, sep = " ")
movies <- read.csv("./omni_ggvp_SNP_pos.csv","omni_SNP_pos.csv", header=T, sep=";") 
ggvp <- read.csv("./omni_ggvp_SNP_pos.csv")
all <- read.csv("./omni_SNP_pos.csv")

upset(movies)


Omni_ggvp <- unlist(read.csv("./omni_ggvp_SNP_pos.csv",stringsAsFactors = FALSE,header = TRUE),use.names = FALSE)
Omni <- unlist(read.csv("./omni_SNP_pos.csv",stringsAsFactors = FALSE,header = TRUE),use.names = FALSE)
H3A_ggvp <- unlist(read.csv("./h3a_ggvp_SNP_pos",stringsAsFactors = FALSE,header = TRUE),use.names = FALSE)
H3A <- unlist(read.csv("./h3a_SNP_pos",stringsAsFactors = FALSE,header = TRUE),use.names = FALSE)
sets <- list(Omni_ggvp=Omni_ggvp,Omni=Omni,H3A_ggvp=H3A_ggvp,H3A=H3A)
upset(fromList(sets),sets.bar.color = "gray23", order.by = "freq", nsets = 4,
      point.size = 3.5, line.size = 2, 
      mainbar.y.label = "SNPs Intersections", 
      sets.x.label = "SNPs Per Array")

