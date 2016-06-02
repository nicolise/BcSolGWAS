#reformat bigRR output data
#Nicole E Soltis

#--------------------------------------------------------
rm(list=ls())
#setwd("~/Documents/GitRepos/BcSolGWAS/data/SNP_files")
setwd("~/Projects/BcSolGWAS/data/bigRRout")
#Import data
#reorganize file Sl_LesionSize.HEM.csv
HEMdat <- read.csv("Sl_LesionSize.HEM.csv")

#first remove first 4 rows (threshold data)
HEMthresh <- HEMdat[1:4,]
HEMdat <- HEMdat[-c(1:4),]
HEMdat2 <- HEMdat

#RF-Basically it has two columns containing Chrom and pos separately instead of just one column with eg "III.57894". This is easiest to do with code below:
library(tidyr)
names(HEMdat)

#my problem: some are formatted as Chromosome1.252
#others are formatted as Chromosome1.1.252

HEMpractice <- HEMdat[c(1:200,50000:50300),]
unique(HEMpractice$X.1)

grx <- glob2rx("Chromosome*.*.*")
HEMsub <- HEMdat[grep(grx, HEMpractice$X.1),]
HEMpractice2 <- HEMpractice

HEMdat2$X.1 <- gsub(pattern = "Chromosome1\\.[0-9]\\.", replacement = "Chromosome1.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome2\\.[0-9]\\.", replacement = "Chromosome2.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome3\\.[0-9]\\.", replacement = "Chromosome3.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome4\\.[0-9]\\.", replacement = "Chromosome4.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome5\\.[0-9]\\.", replacement = "Chromosome5.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome6\\.[0-9]\\.", replacement = "Chromosome6.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome7\\.[0-9]\\.", replacement = "Chromosome7.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome8\\.[0-9]\\.", replacement = "Chromosome8.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome9\\.[0-9]\\.", replacement = "Chromosome9.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome10\\.[0-9]\\.", replacement = "Chromosome10.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome11\\.[0-9]\\.", replacement = "Chromosome11.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome12\\.[0-9]\\.", replacement = "Chromosome12.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome13\\.[0-9]\\.", replacement = "Chromosome13.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome14\\.[0-9]\\.", replacement = "Chromosome14.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome15\\.[0-9]\\.", replacement = "Chromosome15.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome16\\.[0-9]\\.", replacement = "Chromosome16.", HEMdat2$X.1)
HEMdat2$X.1 <- gsub(pattern = "Chromosome16\\.[0-9][0-9]\\.", replacement = "Chromosome16.", HEMdat2$X.1)

HEMdat3 <- separate (HEMdat2, X.1, into = c("Chrom", "Pos") )

write.csv(HEMdat3, "Sl_LesionSize.HEM.PlotFormat.csv") 
write.csv(HEMthresh, "Sl_LesionSize.HEM.Thresh.csv")
#read in to RunningBigRRwithRF.R, line 172