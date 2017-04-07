#reformat bigRR output data
#Nicole E Soltis

#--------------------------------------------------------
rm(list=ls())
library(tidyr)
#setwd("~/Documents/GitRepos/BcSolGWAS/data/SNP_files")
setwd("~/Projects/BcSolGWAS/data/GWAS_files/04_bigRRoutput/trueMAF/")
#Import data
#reorganize file Sl_LesionSize.HEM.csv
HEMdat <- read.csv("SlBc_domest_trueMAF20.HEM.csv")

#first remove first 4 rows (threshold data)
HEMdat <- HEMdat[,-c(1)]
HEMthresh <- HEMdat[1:8,]
HEMdat <- HEMdat[-c(1:8),]
HEMdat2 <- HEMdat

#split chromosome and segment
names(HEMdat)
unique(HEMdat$X.1)

HEMdat$X.1 <- gsub(pattern = "Chromosome1\\.", replacement = "Chromosome1.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome2\\.", replacement = "Chromosome2.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome3\\.", replacement = "Chromosome3.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome4\\.", replacement = "Chromosome4.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome5\\.", replacement = "Chromosome5.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome6\\.", replacement = "Chromosome6.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome7\\.", replacement = "Chromosome7.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome8\\.", replacement = "Chromosome8.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome9\\.", replacement = "Chromosome9.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome10\\.", replacement = "Chromosome10.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome11\\.", replacement = "Chromosome11.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome12\\.", replacement = "Chromosome12.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome13\\.", replacement = "Chromosome13.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome14\\.", replacement = "Chromosome14.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome15\\.", replacement = "Chromosome15.0.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.", replacement = "Chromosome16.0.", HEMdat$X.1)
unique(HEMdat$X.1)

#Then everything with 1.0.(1.)number I'll change to 1.(1.)number

HEMdat$X.1 <- gsub(pattern = "Chromosome1\\.0\\.1\\.", replacement = "Chromosome1.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome2\\.0\\.1\\.", replacement = "Chromosome2.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome2\\.0\\.2\\.", replacement = "Chromosome2.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome2\\.0\\.3\\.", replacement = "Chromosome2.3.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome2\\.0\\.4\\.", replacement = "Chromosome2.4.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome2\\.0\\.5\\.", replacement = "Chromosome2.5.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome3\\.0\\.1\\.", replacement = "Chromosome3.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome3\\.0\\.2\\.", replacement = "Chromosome3.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome4\\.0\\.1\\.", replacement = "Chromosome4.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome5\\.0\\.1\\.", replacement = "Chromosome5.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome6\\.0\\.1\\.", replacement = "Chromosome6.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome6\\.0\\.2\\.", replacement = "Chromosome6.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome6\\.0\\.3\\.", replacement = "Chromosome6.3.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome7\\.0\\.1\\.", replacement = "Chromosome7.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome7\\.0\\.2\\.", replacement = "Chromosome7.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome8\\.0\\.1\\.", replacement = "Chromosome8.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome8\\.0\\.2\\.", replacement = "Chromosome8.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome9\\.0\\.1\\.", replacement = "Chromosome9.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome10\\.0\\.1\\.", replacement = "Chromosome10.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome11\\.0\\.1\\.", replacement = "Chromosome11.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome12\\.0\\.1\\.", replacement = "Chromosome12.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome13\\.0\\.1\\.", replacement = "Chromosome13.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome13\\.0\\.2\\.", replacement = "Chromosome13.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome14\\.0\\.1\\.", replacement = "Chromosome14.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome14\\.0\\.2\\.", replacement = "Chromosome14.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome15\\.0\\.1\\.", replacement = "Chromosome15.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome15\\.0\\.2\\.", replacement = "Chromosome15.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome15\\.0\\.3\\.", replacement = "Chromosome15.3.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome15\\.0\\.4\\.", replacement = "Chromosome15.4.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.1\\.", replacement = "Chromosome16.1.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.2\\.", replacement = "Chromosome16.2.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.3\\.", replacement = "Chromosome16.3.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.4\\.", replacement = "Chromosome16.4.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.5\\.", replacement = "Chromosome16.5.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.6\\.", replacement = "Chromosome16.6.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.7\\.", replacement = "Chromosome16.7.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.8\\.", replacement = "Chromosome16.8.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.9\\.", replacement = "Chromosome16.9.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.10\\.", replacement = "Chromosome16.10.", HEMdat$X.1)
HEMdat$X.1 <- gsub(pattern = "Chromosome16\\.0\\.11\\.", replacement = "Chromosome16.11.", HEMdat$X.1)

HEMdat2 <- separate (HEMdat, X.1, into = c("Chrom", "Segment", "Pos") )
#double check
unique(HEMdat2$Chrom)
unique(HEMdat2$Segment)

write.csv(HEMdat2, "SlBc_domest_trueMAF20.HEM.PlotFormat.csv") 
write.csv(HEMthresh, "SlBc_domest_trueMAF20.HEM.Thresh.csv")
#read in to 06_plots