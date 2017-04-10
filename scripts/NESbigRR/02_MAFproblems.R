rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcSolGWAS/data/GWAS_files/")

 SNPsMAF5 <- read.csv("02_csvPrep/snps_maf5.csv")
 SNPsMAF10 <- read.csv("02_csvPrep/snps_maf10.csv")
 SNPsMAF20 <- read.csv("02_csvPrep/snps_maf20.csv")
 
mySNPs <- SNPsMAF5